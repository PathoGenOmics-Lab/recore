use crate::parser::PresenceAbsence;
use crate::search::Hit;
use crate::types::{ClusterIdx, GenomeIdx};
use ahash::{HashMap, HashMapExt};

// ── Rescue types ──

#[derive(Debug, Clone, PartialEq)]
pub enum RescueType {
    /// Gene is present but Panaroo dropped it during graph cleaning.
    Dropped,
    /// Gene is present in the assembly at an unannotated CDS locus.
    MissingCds,
    /// Gene's locus is assigned to a different Panaroo cluster (mis-clustering).
    Split { other_cluster: String },
    /// Gene is broken into multiple fragments (IS, frameshift, assembly break).
    Fragmented {
        n_fragments: usize,
        total_coverage: f64,
    },
    /// Multiple minimap2 hits reconstruct the gene but lack evidence of biological disruption.
    AlignmentSplit {
        n_fragments: usize,
        total_coverage: f64,
    },
    /// Start site variant: gene is there but with a different start annotation.
    StartVariant { length_diff_bp: i64 },
    /// Gene hit overlaps a paralog from a different cluster (identity <95%).
    Paralog { other_cluster: String },
    /// Gene is a pseudogene: present but disrupted (frameshift, premature stop, partial deletion).
    Pseudogene { reason: PseudogeneReason },
}

/// Reason why a gene was classified as a pseudogene.
#[derive(Debug, Clone, PartialEq)]
pub enum PseudogeneReason {
    /// Internal stop codon(s) found in the aligned region.
    InternalStop { n_stops: usize },
    /// Frameshift: alignment has indels not divisible by 3.
    Frameshift,
    /// Partial deletion: high identity but missing 5' or 3' end (coverage 60-80%).
    PartialDeletion { coverage: f64 },
}

impl std::fmt::Display for RescueType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RescueType::Dropped => write!(f, "DROPPED"),
            RescueType::MissingCds => write!(f, "MISSING_CDS"),
            RescueType::Split { other_cluster } => write!(f, "SPLIT:{}", other_cluster),
            RescueType::Fragmented {
                n_fragments,
                total_coverage,
            } => {
                write!(
                    f,
                    "FRAGMENTED:{}frags:{:.0}%cov",
                    n_fragments, total_coverage
                )
            }
            RescueType::AlignmentSplit {
                n_fragments,
                total_coverage,
            } => {
                write!(
                    f,
                    "ALIGNMENT_SPLIT:{}frags:{:.0}%cov",
                    n_fragments, total_coverage
                )
            }
            RescueType::StartVariant { length_diff_bp } => {
                write!(f, "START_VARIANT:{:+}bp", length_diff_bp)
            }
            RescueType::Paralog { other_cluster } => write!(f, "PARALOG:{}", other_cluster),
            RescueType::Pseudogene { ref reason } => match reason {
                PseudogeneReason::InternalStop { n_stops } => {
                    write!(f, "PSEUDOGENE:{}stops", n_stops)
                }
                PseudogeneReason::Frameshift => write!(f, "PSEUDOGENE:frameshift"),
                PseudogeneReason::PartialDeletion { coverage } => {
                    write!(f, "PSEUDOGENE:partial:{:.0}%cov", coverage)
                }
            },
        }
    }
}

// ── Rescue event ──

#[derive(Debug, Clone)]
pub struct RescueEvent {
    pub cluster_idx: ClusterIdx,
    pub genome_idx: GenomeIdx,
    pub rescue_type: RescueType,
    pub identity: f64,
    pub coverage: f64,
    pub contig: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub annotation_match: Option<String>,
    pub confidence: f64,
    pub synteny_ok: Option<bool>,
    pub synteny_evidence: SyntenyEvidence,
    pub fragment_evidence: FragmentEvidence,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct SyntenyEvidence {
    pub score: Option<f64>,
    pub support: usize,
    pub reference_contexts: usize,
    pub n_neighbors: usize,
    pub matched_upstream: usize,
    pub matched_downstream: usize,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct FragmentEvidence {
    pub n_fragments: usize,
    pub total_query_gap_bp: u64,
    pub max_query_gap_bp: u64,
    pub total_target_gap_bp: u64,
    pub max_target_gap_bp: u64,
    pub colinear: Option<bool>,
    pub same_cds: Option<bool>,
    pub note: Option<String>,
}

// ── GFF index ──

#[derive(Debug, Clone)]
pub struct GffGene {
    pub contig: String,
    pub start: u64,
    pub end: u64,
    #[allow(dead_code)]
    pub strand: char,
    pub gene_id: String,
}

pub struct GffIndex {
    pub genes: Vec<GffGene>,
    contig_ranges: HashMap<String, (usize, usize)>,
    gene_id_index: HashMap<String, usize>,
    gene_aliases: HashMap<String, Vec<String>>,
    gene_labels: HashMap<String, Vec<String>>,
}

impl GffIndex {
    #[allow(dead_code)]
    pub fn build(mut genes: Vec<GffGene>) -> Self {
        Self::build_with_aliases(
            genes
                .drain(..)
                .map(|gene| (gene, Vec::new()))
                .collect::<Vec<_>>(),
        )
    }

    pub fn build_with_aliases(mut records: Vec<(GffGene, Vec<String>)>) -> Self {
        let records = records
            .drain(..)
            .map(|(gene, aliases)| (gene, aliases.clone(), aliases))
            .collect();
        Self::build_with_aliases_and_labels(records)
    }

    pub fn build_with_aliases_and_labels(
        mut records: Vec<(GffGene, Vec<String>, Vec<String>)>,
    ) -> Self {
        records.sort_by(|a, b| a.0.contig.cmp(&b.0.contig).then(a.0.start.cmp(&b.0.start)));

        let genes: Vec<GffGene> = records.iter().map(|(gene, _, _)| gene.clone()).collect();

        let mut contig_ranges = HashMap::new();
        if !genes.is_empty() {
            let mut current_contig = &genes[0].contig;
            let mut range_start = 0;
            for (i, g) in genes.iter().enumerate() {
                if &g.contig != current_contig {
                    contig_ranges.insert(current_contig.clone(), (range_start, i));
                    current_contig = &g.contig;
                    range_start = i;
                }
            }
            contig_ranges.insert(current_contig.clone(), (range_start, genes.len()));
        }

        let mut gene_id_index = HashMap::new();
        let mut gene_aliases = HashMap::new();
        let mut gene_labels = HashMap::new();
        for (i, (g, aliases, labels)) in records.iter().enumerate() {
            if !g.gene_id.is_empty() {
                gene_id_index.insert(g.gene_id.clone(), i);
                let mut all_aliases = vec![g.gene_id.clone()];
                for alias in aliases {
                    if !alias.is_empty() && !all_aliases.iter().any(|x| x == alias) {
                        gene_id_index.entry(alias.clone()).or_insert(i);
                        all_aliases.push(alias.clone());
                    }
                }
                let mut all_labels = all_aliases.clone();
                for label in labels {
                    if !label.is_empty() && !all_labels.iter().any(|x| x == label) {
                        all_labels.push(label.clone());
                    }
                }
                gene_aliases.insert(g.gene_id.clone(), all_aliases);
                gene_labels.insert(g.gene_id.clone(), all_labels);
            }
        }

        GffIndex {
            genes,
            contig_ranges,
            gene_id_index,
            gene_aliases,
            gene_labels,
        }
    }

    pub fn get_by_id(&self, gene_id: &str) -> Option<&GffGene> {
        self.gene_id_index.get(gene_id).map(|&i| &self.genes[i])
    }

    pub fn aliases_for_gene<'a>(&'a self, gene: &'a GffGene) -> Vec<&'a str> {
        self.gene_aliases
            .get(&gene.gene_id)
            .map(|aliases| aliases.iter().map(String::as_str).collect())
            .unwrap_or_else(|| vec![gene.gene_id.as_str()])
    }

    pub fn labels_for_gene<'a>(&'a self, gene: &'a GffGene) -> Vec<&'a str> {
        self.gene_labels
            .get(&gene.gene_id)
            .map(|labels| labels.iter().map(String::as_str).collect())
            .unwrap_or_else(|| self.aliases_for_gene(gene))
    }

    pub fn cluster_for_gene(
        &self,
        gene: &GffGene,
        gene_to_cluster: &HashMap<String, ClusterIdx>,
    ) -> Option<ClusterIdx> {
        gene_to_cluster.get(&gene.gene_id).copied().or_else(|| {
            self.gene_aliases
                .get(&gene.gene_id)
                .and_then(|aliases| aliases.iter().find_map(|alias| gene_to_cluster.get(alias)))
                .copied()
        })
    }

    pub fn overlapping(&self, contig: &str, query_start: u64, query_end: u64) -> Vec<&GffGene> {
        let Some(&(range_start, range_end)) = self.contig_ranges.get(contig) else {
            return Vec::new();
        };

        let slice = &self.genes[range_start..range_end];
        let scan_start = slice.partition_point(|g| g.end <= query_start);

        let mut results = Vec::new();
        for g in &slice[scan_start..] {
            if g.start >= query_end {
                break;
            }
            results.push(g);
        }
        results
    }

    pub fn flanking(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> (Option<&GffGene>, Option<&GffGene>) {
        let Some(&(range_start, range_end)) = self.contig_ranges.get(contig) else {
            return (None, None);
        };

        let slice = &self.genes[range_start..range_end];
        let pos = slice.partition_point(|g| g.end <= start);

        let upstream = if pos > 0 { Some(&slice[pos - 1]) } else { None };

        let downstream = slice[pos..].iter().find(|g| g.start >= end);

        (upstream, downstream)
    }

    pub fn context_window(
        &self,
        contig: &str,
        start: u64,
        end: u64,
        center_gene_id: Option<&str>,
        n_neighbors: usize,
    ) -> (Vec<&GffGene>, Vec<&GffGene>) {
        let Some(&(range_start, range_end)) = self.contig_ranges.get(contig) else {
            return (Vec::new(), Vec::new());
        };
        let slice = &self.genes[range_start..range_end];
        if slice.is_empty() {
            return (Vec::new(), Vec::new());
        }

        let center_pos = center_gene_id
            .and_then(|gene_id| self.gene_id_index.get(gene_id).copied())
            .and_then(|idx| {
                if idx >= range_start && idx < range_end {
                    Some(idx - range_start)
                } else {
                    None
                }
            })
            .or_else(|| {
                slice
                    .iter()
                    .enumerate()
                    .filter_map(|(i, gene)| {
                        let overlap = gene.end.min(end).saturating_sub(gene.start.max(start));
                        let query_len = end.saturating_sub(start);
                        let gene_len = gene.end.saturating_sub(gene.start);
                        let substantial = overlap > 0
                            && ((query_len > 0 && overlap as f64 / query_len as f64 >= 0.5)
                                || (gene_len > 0 && overlap as f64 / gene_len as f64 >= 0.5));
                        substantial.then_some((i, overlap))
                    })
                    .max_by_key(|(_, overlap)| *overlap)
                    .map(|(i, _)| i)
            });

        if let Some(center_pos) = center_pos {
            let upstream = slice[..center_pos].iter().rev().take(n_neighbors).collect();
            let downstream = slice
                .iter()
                .skip(center_pos + 1)
                .take(n_neighbors)
                .collect();
            return (upstream, downstream);
        }

        let midpoint = start.saturating_add(end.saturating_sub(start) / 2);
        let pos = slice.partition_point(|g| {
            let gene_mid = g.start.saturating_add(g.end.saturating_sub(g.start) / 2);
            gene_mid < midpoint
        });
        let upstream = slice[..pos].iter().rev().take(n_neighbors).collect();
        let downstream = slice[pos..].iter().take(n_neighbors).collect();
        (upstream, downstream)
    }
}

// ── Gene-to-cluster index ──

/// Build a lookup from gene_id → ClusterIdx for a specific genome.
pub fn build_gene_to_cluster(pa: &PresenceAbsence, gi: GenomeIdx) -> HashMap<String, ClusterIdx> {
    let mut map = HashMap::new();
    for ci_val in 0..pa.n_clusters() {
        let ci = ClusterIdx(ci_val);
        if !pa.is_present(ci, gi) {
            continue;
        }
        let cell = pa.cell(ci, gi);
        for gene_id in cell.split(';') {
            let gene_id = gene_id.trim();
            if !gene_id.is_empty() {
                map.insert(gene_id.to_string(), ci);
                // Also index without Panaroo suffixes for GFF matching.
                for stripped in panaroo_gene_id_suffix_variants(gene_id).into_iter().skip(1) {
                    map.entry(stripped.to_string()).or_insert(ci);
                }
            }
        }
    }
    map
}

pub(crate) fn strip_panaroo_gene_suffix(gene_id: &str) -> Option<&str> {
    if let Some(stripped) = gene_id.strip_suffix("_pseudo") {
        return Some(stripped);
    }
    let (base, suffix) = gene_id.rsplit_once("_dup_")?;
    if !base.is_empty() && !suffix.is_empty() && suffix.chars().all(|ch| ch.is_ascii_digit()) {
        Some(base)
    } else {
        None
    }
}

pub(crate) fn panaroo_gene_id_suffix_variants(gene_id: &str) -> Vec<&str> {
    let mut variants = vec![gene_id];
    let mut current = gene_id;
    while let Some(stripped) = strip_panaroo_gene_suffix(current) {
        if variants.contains(&stripped) {
            break;
        }
        variants.push(stripped);
        current = stripped;
    }
    variants
}

fn gene_id_variants(gene_id: &str) -> impl Iterator<Item = &str> {
    panaroo_gene_id_suffix_variants(gene_id).into_iter()
}

fn lookup_gene_cluster_alias(
    map: &HashMap<String, ClusterIdx>,
    gene_id: &str,
) -> Option<ClusterIdx> {
    gene_id_variants(gene_id).find_map(|id| map.get(id).copied())
}

fn insert_gene_cluster_alias(map: &mut HashMap<String, ClusterIdx>, gene_id: &str, ci: ClusterIdx) {
    for id in gene_id_variants(gene_id) {
        if !id.is_empty() {
            map.entry(id.to_string()).or_insert(ci);
        }
    }
}

/// Build a gene_id → ClusterIdx lookup and conservatively add GFF aliases.
///
/// Panaroo cells can contain a gene name while the GFF uses a locus_tag (or the
/// reverse). If exactly one alias of a GFF CDS already maps to a cluster, all
/// aliases of that CDS are mapped to the same cluster. Ambiguous CDS aliases
/// that point to multiple clusters are left untouched.
pub fn build_gene_to_cluster_with_gff_aliases(
    pa: &PresenceAbsence,
    gi: GenomeIdx,
    gff_index: Option<&GffIndex>,
) -> HashMap<String, ClusterIdx> {
    let mut map = build_gene_to_cluster(pa, gi);
    let Some(gff) = gff_index else {
        return map;
    };

    for gene in &gff.genes {
        let labels = gff.aliases_for_gene(gene);
        let mut candidates: Vec<ClusterIdx> = labels
            .iter()
            .filter_map(|label| lookup_gene_cluster_alias(&map, label))
            .collect();
        candidates.sort_unstable();
        candidates.dedup();
        if candidates.len() != 1 {
            continue;
        }
        let ci = candidates[0];
        for label in labels {
            insert_gene_cluster_alias(&mut map, label, ci);
        }
    }

    map
}

// ── Fragment analysis ──

const BIOLOGICAL_FRAGMENT_MIN_GAP_BP: u64 = 50;
#[cfg(test)]
const SYNTENY_CONTEXT_GENES: usize = 2;
const SYNTENY_MAX_CONTEXT_GENES: usize = 5;
const SYNTENY_REFERENCE_SUPPORT_THRESHOLD: f64 = 0.50;
const MISSING_CDS_MIN_IDENTITY: f64 = 95.0;

#[derive(Debug, Clone)]
struct FragmentChain {
    hits: Vec<Hit>,
    coverage: f64,
    avg_identity: f64,
    n_fragments: usize,
    contig: String,
    strand: char,
    start: u64,
    end: u64,
    evidence: FragmentEvidence,
}

#[derive(Debug, Clone)]
struct ChainSegment {
    query_start: u64,
    query_end: u64,
    target_start: u64,
    target_end: u64,
}

#[cfg(test)]
fn merged_query_coverage(hits: &[&Hit]) -> (f64, usize) {
    if hits.is_empty() {
        return (0.0, 0);
    }

    let mut by_allele: HashMap<(&str, u64), Vec<(u64, u64)>> = HashMap::new();
    for h in hits {
        by_allele
            .entry((&h.query, h.query_length))
            .or_default()
            .push((h.query_start, h.query_end));
    }

    let mut best_coverage = 0.0f64;
    let mut best_fragments = 0usize;

    for ((_qname, qlen), mut intervals) in by_allele {
        if qlen == 0 {
            continue;
        }
        intervals.sort_by_key(|&(s, _)| s);

        let mut merged: Vec<(u64, u64)> = Vec::new();
        for (s, e) in intervals {
            if let Some(last) = merged.last_mut() {
                if s <= last.1 {
                    last.1 = last.1.max(e);
                } else {
                    merged.push((s, e));
                }
            } else {
                merged.push((s, e));
            }
        }

        let covered: u64 = merged.iter().map(|(s, e)| e - s).sum();
        let coverage = (covered as f64 / qlen as f64) * 100.0;

        if coverage > best_coverage {
            best_coverage = coverage;
            best_fragments = merged.len();
        }
    }

    (best_coverage, best_fragments)
}

fn build_fragment_chain(hits: &[&Hit]) -> Option<FragmentChain> {
    let first = *hits.first()?;
    if first.query_length == 0 {
        return None;
    }

    let mut sorted: Vec<&Hit> = hits.to_vec();
    sorted.sort_by(|a, b| {
        a.query_start
            .cmp(&b.query_start)
            .then(a.query_end.cmp(&b.query_end))
            .then(a.target_start.cmp(&b.target_start))
            .then(a.target_end.cmp(&b.target_end))
    });

    let mut segments: Vec<ChainSegment> = Vec::new();
    for hit in &sorted {
        if let Some(last) = segments.last_mut() {
            if hit.query_start <= last.query_end {
                last.query_end = last.query_end.max(hit.query_end);
                last.target_start = last.target_start.min(hit.target_start);
                last.target_end = last.target_end.max(hit.target_end);
                continue;
            }
        }
        segments.push(ChainSegment {
            query_start: hit.query_start,
            query_end: hit.query_end,
            target_start: hit.target_start,
            target_end: hit.target_end,
        });
    }

    let covered: u64 = segments.iter().map(|s| s.query_end - s.query_start).sum();
    let coverage = (covered as f64 / first.query_length as f64) * 100.0;
    let n_fragments = segments.len();

    let mut total_query_gap_bp = 0u64;
    let mut max_query_gap_bp = 0u64;
    let mut total_target_gap_bp = 0u64;
    let mut max_target_gap_bp = 0u64;
    let mut colinear = true;

    for pair in segments.windows(2) {
        let prev = &pair[0];
        let next = &pair[1];
        let query_gap = next.query_start.saturating_sub(prev.query_end);
        total_query_gap_bp += query_gap;
        max_query_gap_bp = max_query_gap_bp.max(query_gap);

        let target_gap = if first.strand == '-' {
            if next.target_end <= prev.target_start {
                prev.target_start - next.target_end
            } else {
                colinear = false;
                0
            }
        } else if next.target_start >= prev.target_end {
            next.target_start - prev.target_end
        } else {
            colinear = false;
            0
        };
        total_target_gap_bp += target_gap;
        max_target_gap_bp = max_target_gap_bp.max(target_gap);
    }

    let total_matches: u64 = sorted.iter().map(|h| h.matches).sum();
    let total_aln: u64 = sorted.iter().map(|h| h.alignment_length).sum();
    let avg_identity = if total_aln > 0 {
        (total_matches as f64 / total_aln as f64) * 100.0
    } else {
        0.0
    };

    let start = sorted.iter().map(|h| h.target_start).min().unwrap_or(0);
    let end = sorted.iter().map(|h| h.target_end).max().unwrap_or(0);

    Some(FragmentChain {
        hits: sorted.into_iter().cloned().collect(),
        coverage,
        avg_identity,
        n_fragments,
        contig: first.target.clone(),
        strand: first.strand,
        start,
        end,
        evidence: FragmentEvidence {
            n_fragments,
            total_query_gap_bp,
            max_query_gap_bp,
            total_target_gap_bp,
            max_target_gap_bp,
            colinear: Some(colinear),
            same_cds: None,
            note: None,
        },
    })
}

fn best_fragment_chain(hits: &[&Hit]) -> Option<FragmentChain> {
    let mut groups: HashMap<(String, u64, String, char), Vec<&Hit>> = HashMap::new();
    for hit in hits {
        groups
            .entry((
                hit.query.clone(),
                hit.query_length,
                hit.target.clone(),
                hit.strand,
            ))
            .or_default()
            .push(*hit);
    }

    groups
        .into_values()
        .filter_map(|group_hits| build_fragment_chain(&group_hits))
        .max_by(|a, b| {
            a.coverage
                .total_cmp(&b.coverage)
                .then_with(|| a.avg_identity.total_cmp(&b.avg_identity))
                .then_with(|| b.n_fragments.cmp(&a.n_fragments))
                .then_with(|| {
                    b.evidence
                        .max_target_gap_bp
                        .cmp(&a.evidence.max_target_gap_bp)
                })
                .then_with(|| b.contig.cmp(&a.contig))
                .then_with(|| b.start.cmp(&a.start))
                .then_with(|| b.end.cmp(&a.end))
                .then_with(|| b.strand.cmp(&a.strand))
                .then_with(|| {
                    let a_query = a.hits.first().map(|hit| hit.query.as_str()).unwrap_or("");
                    let b_query = b.hits.first().map(|hit| hit.query.as_str()).unwrap_or("");
                    b_query.cmp(a_query)
                })
        })
}

fn compatible_single_cds<'a>(
    chain: &FragmentChain,
    gff_index: Option<&'a GffIndex>,
) -> Option<&'a GffGene> {
    let gff = gff_index?;
    let overlapping = gff.overlapping(&chain.contig, chain.start, chain.end);
    let span_len = chain.end.saturating_sub(chain.start) as f64;

    for gene in overlapping {
        let gene_len = gene.end.saturating_sub(gene.start) as i64;
        let span_len_i = chain.end.saturating_sub(chain.start) as i64;
        if (gene_len - span_len_i).unsigned_abs() >= 300 {
            continue;
        }

        let span_overlap = chain
            .end
            .min(gene.end)
            .saturating_sub(chain.start.max(gene.start));
        if span_len > 0.0 && (span_overlap as f64 / span_len) < 0.8 {
            continue;
        }

        let all_hits_in_gene = chain.hits.iter().all(|hit| {
            let hit_len = hit.target_end.saturating_sub(hit.target_start) as f64;
            let overlap = hit
                .target_end
                .min(gene.end)
                .saturating_sub(hit.target_start.max(gene.start)) as f64;
            hit_len > 0.0 && overlap / hit_len >= 0.8
        });
        if all_hits_in_gene {
            return Some(gene);
        }
    }

    None
}

fn classify_matched_gene(
    matched_gene_id: &str,
    identity: f64,
    coverage: f64,
    ci: ClusterIdx,
    pa: &PresenceAbsence,
    gene_to_cluster: &HashMap<String, ClusterIdx>,
) -> (RescueType, f64) {
    let base_confidence = identity / 100.0 * coverage / 100.0;
    if let Some(&other_ci) = gene_to_cluster.get(matched_gene_id) {
        if other_ci != ci {
            if identity >= 95.0 {
                return (
                    RescueType::Split {
                        other_cluster: pa.cluster_name(other_ci).to_string(),
                    },
                    base_confidence * 1.1,
                );
            }
            return (
                RescueType::Paralog {
                    other_cluster: pa.cluster_name(other_ci).to_string(),
                },
                base_confidence * 0.7,
            );
        }
    }
    (RescueType::Dropped, base_confidence)
}

fn detect_start_variant_from_span(
    identity: f64,
    coverage: f64,
    contig: &str,
    start: u64,
    end: u64,
    gff_index: Option<&GffIndex>,
) -> Option<(i64, String)> {
    let gff = gff_index?;

    if identity < 98.0 || coverage < 90.0 {
        return None;
    }

    let overlapping = gff.overlapping(contig, start, end);
    let hit_len = end.saturating_sub(start) as f64;
    let target_aln_len = end.saturating_sub(start) as i64;

    for gene in &overlapping {
        let gene_len = (gene.end - gene.start) as i64;
        let diff = gene_len - target_aln_len;

        if diff.unsigned_abs() < 300 && diff != 0 {
            let overlap_start = start.max(gene.start);
            let overlap_end = end.min(gene.end);
            let overlap = overlap_end.saturating_sub(overlap_start) as f64;

            if hit_len > 0.0 && overlap / hit_len > 0.8 {
                return Some((diff, gene.gene_id.clone()));
            }
        }
    }

    None
}

// ── Start-site variant detection ──

fn detect_start_variant(hit: &Hit, gff_index: Option<&GffIndex>) -> Option<(i64, String)> {
    detect_start_variant_from_span(
        hit.identity(),
        hit.query_coverage(),
        &hit.target,
        hit.target_start,
        hit.target_end,
        gff_index,
    )
}

// ── Main classification ──

pub fn classify_hits(
    hits: &[Hit],
    ci: ClusterIdx,
    gi: GenomeIdx,
    pa: &PresenceAbsence,
    gff_index: Option<&GffIndex>,
    gene_to_cluster: &HashMap<String, ClusterIdx>,
    min_coverage: f64,
) -> Option<RescueEvent> {
    if hits.is_empty() {
        return None;
    }

    let mut sorted_hits: Vec<&Hit> = hits.iter().collect();
    sorted_hits.sort_by(|a, b| {
        b.score()
            .total_cmp(&a.score())
            .then_with(|| b.query_coverage().total_cmp(&a.query_coverage()))
            .then_with(|| b.identity().total_cmp(&a.identity()))
            .then_with(|| b.query_length.cmp(&a.query_length))
            .then(a.target.cmp(&b.target))
            .then(a.target_start.cmp(&b.target_start))
            .then(a.target_end.cmp(&b.target_end))
            .then(a.query.cmp(&b.query))
    });

    // ── Path 1: Single full-coverage hit ──
    if let Some(missing_cds_hit) = sorted_hits
        .iter()
        .copied()
        .filter(|hit| hit.query_coverage() >= min_coverage)
        .find(|hit| is_missing_cds_candidate(hit, gff_index))
    {
        let confidence =
            missing_cds_hit.identity() / 100.0 * missing_cds_hit.query_coverage() / 100.0;
        return Some(RescueEvent {
            cluster_idx: ci,
            genome_idx: gi,
            rescue_type: RescueType::MissingCds,
            identity: missing_cds_hit.identity(),
            coverage: missing_cds_hit.query_coverage(),
            contig: missing_cds_hit.target.clone(),
            start: missing_cds_hit.target_start,
            end: missing_cds_hit.target_end,
            strand: missing_cds_hit.strand,
            annotation_match: None,
            confidence,
            synteny_ok: None,
            synteny_evidence: SyntenyEvidence::default(),
            fragment_evidence: FragmentEvidence {
                note: Some("missing_cds_no_substantial_gff_overlap".to_string()),
                ..FragmentEvidence::default()
            },
        });
    }

    let best = sorted_hits[0];

    if best.query_coverage() >= min_coverage {
        let (rescue_type, annotation_match, confidence) =
            classify_single_hit(best, ci, pa, gff_index, gene_to_cluster);

        let (final_type, final_ann, final_conf) = if matches!(rescue_type, RescueType::Dropped) {
            if let Some((length_diff, matched_gene_id)) = detect_start_variant(best, gff_index) {
                (
                    RescueType::StartVariant {
                        length_diff_bp: length_diff,
                    },
                    Some(matched_gene_id),
                    best.identity() / 100.0 * 0.95,
                )
            } else {
                (rescue_type, annotation_match, confidence)
            }
        } else {
            (rescue_type, annotation_match, confidence)
        };

        return Some(RescueEvent {
            cluster_idx: ci,
            genome_idx: gi,
            rescue_type: final_type,
            identity: best.identity(),
            coverage: best.query_coverage(),
            contig: best.target.clone(),
            start: best.target_start,
            end: best.target_end,
            strand: best.strand,
            annotation_match: final_ann,
            confidence: final_conf,
            synteny_ok: None,
            synteny_evidence: SyntenyEvidence::default(),
            fragment_evidence: FragmentEvidence::default(),
        });
    }

    // ── Path 2: Conservative fragment chain (same allele, contig and strand) ──
    let best_chain = best_fragment_chain(&sorted_hits);
    if let Some(mut chain) = best_chain.clone() {
        if chain.coverage >= min_coverage && chain.n_fragments >= 2 {
            if chain.evidence.colinear != Some(true) {
                return None;
            }

            if let Some(gene) = compatible_single_cds(&chain, gff_index) {
                chain.evidence.same_cds = Some(true);
                chain.evidence.note = Some("same_cds_alignment_split".to_string());
                let (rescue_type, annotation_match, confidence) = {
                    let (rescue_type, confidence) = classify_matched_gene(
                        &gene.gene_id,
                        chain.avg_identity,
                        chain.coverage,
                        ci,
                        pa,
                        gene_to_cluster,
                    );
                    (rescue_type, Some(gene.gene_id.clone()), confidence)
                };

                let (final_type, final_ann, final_conf) = if matches!(
                    rescue_type,
                    RescueType::Dropped
                ) {
                    if let Some((length_diff, matched_gene_id)) = detect_start_variant_from_span(
                        chain.avg_identity,
                        chain.coverage,
                        &chain.contig,
                        chain.start,
                        chain.end,
                        gff_index,
                    ) {
                        (
                            RescueType::StartVariant {
                                length_diff_bp: length_diff,
                            },
                            Some(matched_gene_id),
                            chain.avg_identity / 100.0 * 0.95,
                        )
                    } else {
                        (rescue_type, annotation_match, confidence)
                    }
                } else {
                    (rescue_type, annotation_match, confidence)
                };

                return Some(RescueEvent {
                    cluster_idx: ci,
                    genome_idx: gi,
                    rescue_type: final_type,
                    identity: chain.avg_identity,
                    coverage: chain.coverage,
                    contig: chain.contig,
                    start: chain.start,
                    end: chain.end,
                    strand: chain.strand,
                    annotation_match: final_ann,
                    confidence: final_conf,
                    synteny_ok: None,
                    synteny_evidence: SyntenyEvidence::default(),
                    fragment_evidence: chain.evidence,
                });
            }

            chain.evidence.same_cds = Some(false);

            // If fragments have high identity and small target gaps not divisible by 3,
            // this is stronger evidence for pseudogenization than generic fragmentation.
            let chain_refs: Vec<&Hit> = chain.hits.iter().collect();
            let is_frameshift = chain.avg_identity >= 95.0
                && chain.n_fragments >= 2
                && detect_frameshift_pattern(&chain_refs);

            let has_large_gap = chain.evidence.max_target_gap_bp >= BIOLOGICAL_FRAGMENT_MIN_GAP_BP
                || chain.evidence.max_query_gap_bp >= BIOLOGICAL_FRAGMENT_MIN_GAP_BP;

            let (rescue_type, confidence, note) = if is_frameshift {
                (
                    RescueType::Pseudogene {
                        reason: PseudogeneReason::Frameshift,
                    },
                    (chain.avg_identity / 100.0) * (chain.coverage / 100.0) * 0.5,
                    "frameshift_gap",
                )
            } else if has_large_gap {
                (
                    RescueType::Fragmented {
                        n_fragments: chain.n_fragments,
                        total_coverage: chain.coverage,
                    },
                    (chain.avg_identity / 100.0) * (chain.coverage / 100.0) * 0.8,
                    "biological_gap",
                )
            } else {
                (
                    RescueType::AlignmentSplit {
                        n_fragments: chain.n_fragments,
                        total_coverage: chain.coverage,
                    },
                    (chain.avg_identity / 100.0) * (chain.coverage / 100.0) * 0.9,
                    "alignment_split",
                )
            };
            chain.evidence.note = Some(note.to_string());

            return Some(RescueEvent {
                cluster_idx: ci,
                genome_idx: gi,
                rescue_type,
                identity: chain.avg_identity,
                coverage: chain.coverage,
                contig: chain.contig,
                start: chain.start,
                end: chain.end,
                strand: chain.strand,
                annotation_match: None,
                confidence,
                synteny_ok: None,
                synteny_evidence: SyntenyEvidence::default(),
                fragment_evidence: chain.evidence,
            });
        }
    }

    // ── Path 3: Pseudogene — partial deletion ──
    // High identity hit covering 60-80% of the query → truncated gene (partial deletion).
    // This is below min_coverage so Path 1 didn't catch it, but biologically informative.
    if best.identity() >= 90.0 {
        let cov = best.query_coverage();
        let effective_cov = if let Some(chain) = &best_chain {
            if chain.evidence.colinear == Some(true) {
                chain.coverage
            } else {
                cov
            }
        } else {
            cov
        };
        if effective_cov >= 60.0 && effective_cov < min_coverage {
            let (contig, start, end, strand, fragment_evidence) = if let Some(chain) = best_chain {
                (
                    chain.contig,
                    chain.start,
                    chain.end,
                    chain.strand,
                    chain.evidence,
                )
            } else {
                (
                    best.target.clone(),
                    best.target_start,
                    best.target_end,
                    best.strand,
                    FragmentEvidence::default(),
                )
            };
            return Some(RescueEvent {
                cluster_idx: ci,
                genome_idx: gi,
                rescue_type: RescueType::Pseudogene {
                    reason: PseudogeneReason::PartialDeletion {
                        coverage: effective_cov,
                    },
                },
                identity: best.identity(),
                coverage: effective_cov,
                contig,
                start,
                end,
                strand,
                annotation_match: None,
                confidence: (best.identity() / 100.0) * (effective_cov / 100.0) * 0.4,
                synteny_ok: None,
                synteny_evidence: SyntenyEvidence::default(),
                fragment_evidence,
            });
        }
    }

    None
}

/// Scan a nucleotide sequence for internal stop codons (TAA, TAG, TGA) in all 3 reading frames.
/// Returns the frame with the fewest stops and the count in that frame.
/// Used to detect pseudogenes: a real gene has 0 internal stops in its frame.
pub fn count_internal_stops(seq: &[u8]) -> usize {
    if seq.len() < 6 {
        return 0; // too short to have internal stops
    }

    // Check all 3 frames, return the minimum stop count
    // (the "real" reading frame should have the fewest stops)
    let mut min_stops = usize::MAX;

    for frame in 0..3usize {
        let mut stops = 0usize;
        let mut pos = frame;
        while pos + 3 <= seq.len() {
            let codon = &seq[pos..pos + 3];
            // Last codon is the terminator, not "internal"
            if pos + 3 < seq.len() {
                let upper = [
                    codon[0].to_ascii_uppercase(),
                    codon[1].to_ascii_uppercase(),
                    codon[2].to_ascii_uppercase(),
                ];
                if matches!(
                    upper,
                    [b'T', b'A', b'A'] | [b'T', b'A', b'G'] | [b'T', b'G', b'A']
                ) {
                    stops += 1;
                }
            }
            pos += 3;
        }
        min_stops = min_stops.min(stops);
    }

    min_stops
}

/// Reverse complement a DNA sequence.
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        })
        .collect()
}

/// Check if a rescued hit is a pseudogene by scanning the genome sequence for internal stops.
/// For hits on the minus strand, the sequence is reverse-complemented before scanning.
/// Returns Some(n_stops) if stops found, None if sequence unavailable or clean.
pub fn check_internal_stops_in_hit(
    hit: &Hit,
    genome_seqs: Option<&HashMap<String, Vec<u8>>>,
) -> Option<usize> {
    let seqs = genome_seqs?;
    let contig_seq = seqs.get(&hit.target)?;

    let start = hit.target_start as usize;
    let end = (hit.target_end as usize).min(contig_seq.len());
    if start >= end || end - start < 6 {
        return None;
    }

    let region = &contig_seq[start..end];
    let n_stops = if hit.strand == '-' {
        let rc = reverse_complement(region);
        count_internal_stops(&rc)
    } else {
        count_internal_stops(region)
    };

    if n_stops > 0 {
        Some(n_stops)
    } else {
        None
    }
}

/// Detect frameshift pattern: adjacent fragments with small gaps not divisible by 3.
fn detect_frameshift_pattern(hits: &[&Hit]) -> bool {
    if hits.len() < 2 {
        return false;
    }

    // Sort by target position
    let mut sorted: Vec<&Hit> = hits.to_vec();
    sorted.sort_by_key(|h| h.target_start);

    let mut frameshift_gaps = 0usize;
    let mut total_gaps = 0usize;

    for pair in sorted.windows(2) {
        let gap_start = pair[0].target_end;
        let gap_end = pair[1].target_start;
        if gap_end > gap_start {
            let gap_size = gap_end - gap_start;
            // Small gaps (≤50bp) that aren't multiples of 3 → frameshift signature
            if gap_size <= 50 {
                total_gaps += 1;
                if gap_size % 3 != 0 {
                    frameshift_gaps += 1;
                }
            }
        }
    }

    // At least one gap must be non-divisible by 3
    total_gaps > 0 && frameshift_gaps > 0
}

fn classify_single_hit(
    hit: &Hit,
    ci: ClusterIdx,
    pa: &PresenceAbsence,
    gff_index: Option<&GffIndex>,
    gene_to_cluster: &HashMap<String, ClusterIdx>,
) -> (RescueType, Option<String>, f64) {
    let base_confidence = hit.identity() / 100.0 * hit.query_coverage() / 100.0;

    let Some(gff) = gff_index else {
        return (RescueType::Dropped, None, base_confidence);
    };

    let overlapping = gff.overlapping(&hit.target, hit.target_start, hit.target_end);

    if overlapping.is_empty() {
        return (RescueType::Dropped, None, base_confidence);
    }

    let hit_len = (hit.target_end - hit.target_start) as f64;
    let mut best_gene: Option<(&GffGene, f64)> = None;

    for gene in &overlapping {
        let overlap_start = hit.target_start.max(gene.start);
        let overlap_end = hit.target_end.min(gene.end);
        let overlap = overlap_end.saturating_sub(overlap_start) as f64;
        let gene_len = (gene.end - gene.start) as f64;

        let recip = if hit_len > 0.0 && gene_len > 0.0 {
            (overlap / hit_len).min(overlap / gene_len)
        } else {
            0.0
        };

        if recip > 0.5 && best_gene.is_none_or(|(_, prev)| recip > prev) {
            best_gene = Some((gene, recip));
        }
    }

    let Some((matched_gene, _recip)) = best_gene else {
        return (
            RescueType::Dropped,
            Some(overlapping[0].gene_id.clone()),
            base_confidence,
        );
    };

    if let Some(&other_ci) = gene_to_cluster.get(&matched_gene.gene_id) {
        if other_ci != ci {
            // Gene is annotated but assigned to a different Panaroo cluster.
            // High identity (≥95%) → SPLIT (same gene, mis-clustered by Panaroo).
            // Lower identity → PARALOG (similar gene from a related family).
            if hit.identity() >= 95.0 {
                return (
                    RescueType::Split {
                        other_cluster: pa.cluster_name(other_ci).to_string(),
                    },
                    Some(matched_gene.gene_id.clone()),
                    base_confidence * 1.1,
                );
            } else {
                return (
                    RescueType::Paralog {
                        other_cluster: pa.cluster_name(other_ci).to_string(),
                    },
                    Some(matched_gene.gene_id.clone()),
                    base_confidence * 0.7, // lower confidence: paralog, not real gene
                );
            }
        }
    }

    (
        RescueType::Dropped,
        Some(matched_gene.gene_id.clone()),
        base_confidence,
    )
}

fn is_missing_cds_candidate(hit: &Hit, gff_index: Option<&GffIndex>) -> bool {
    if hit.identity() < MISSING_CDS_MIN_IDENTITY {
        return false;
    }
    let Some(gff) = gff_index else {
        return false;
    };
    let overlapping = gff.overlapping(&hit.target, hit.target_start, hit.target_end);
    if overlapping.is_empty() {
        return true;
    }

    let hit_len = hit.target_end.saturating_sub(hit.target_start) as f64;
    !overlapping.iter().any(|gene| {
        let overlap = hit
            .target_end
            .min(gene.end)
            .saturating_sub(hit.target_start.max(gene.start)) as f64;
        let gene_len = gene.end.saturating_sub(gene.start) as f64;
        hit_len > 0.0 && gene_len > 0.0 && (overlap / hit_len).min(overlap / gene_len) > 0.5
    })
}

// ── Conflict resolution ──

pub fn resolve_conflicts(rescues: &mut Vec<RescueEvent>) {
    if rescues.is_empty() {
        return;
    }

    rescues.sort_by(|a, b| {
        a.genome_idx
            .cmp(&b.genome_idx)
            .then(a.contig.cmp(&b.contig))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
            .then(a.cluster_idx.cmp(&b.cluster_idx))
    });

    let mut to_remove: Vec<bool> = vec![false; rescues.len()];

    for i in 0..rescues.len() {
        if to_remove[i] {
            continue;
        }
        for j in (i + 1)..rescues.len() {
            if to_remove[j] {
                continue;
            }

            let ri = &rescues[i];
            let rj = &rescues[j];

            if ri.genome_idx != rj.genome_idx {
                break;
            }
            if ri.contig != rj.contig {
                break;
            }
            if rj.start >= ri.end {
                break;
            }

            if !is_substantial_overlap(ri, rj) {
                continue;
            }

            if ri.confidence >= rj.confidence {
                to_remove[j] = true;
            } else {
                to_remove[i] = true;
                break;
            }
        }
    }

    let mut idx = 0;
    rescues.retain(|_| {
        let keep = !to_remove[idx];
        idx += 1;
        keep
    });
}

fn is_substantial_overlap(a: &RescueEvent, b: &RescueEvent) -> bool {
    let overlap = a.end.min(b.end).saturating_sub(a.start.max(b.start));
    if overlap == 0 {
        return false;
    }
    let a_len = a.end.saturating_sub(a.start);
    let b_len = b.end.saturating_sub(b.start);
    let shorter = a_len.min(b_len);
    shorter > 0 && (overlap as f64 / shorter as f64) >= 0.50
}

// ── Synteny validation ──

#[derive(Debug, Clone, Default)]
struct SyntenyContext {
    upstream: Vec<ClusterIdx>,
    downstream: Vec<ClusterIdx>,
}

impl SyntenyContext {
    fn has_any_cluster(&self) -> bool {
        !self.upstream.is_empty() || !self.downstream.is_empty()
    }

    fn side_matches(&self, other: &SyntenyContext) -> bool {
        let (_, up_matches) = weighted_side_score(&self.upstream, &other.upstream);
        let (_, down_matches) = weighted_side_score(&self.downstream, &other.downstream);
        let up_comparable = !self.upstream.is_empty() && !other.upstream.is_empty();
        let down_comparable = !self.downstream.is_empty() && !other.downstream.is_empty();

        if up_comparable && down_comparable {
            return up_matches > 0 && down_matches > 0;
        }

        up_matches + down_matches > 0
    }
}

#[derive(Debug, Clone, Default)]
struct SyntenyAssessment {
    passes: bool,
    score: f64,
    support: usize,
    reference_contexts: usize,
    n_neighbors: usize,
    matched_upstream: usize,
    matched_downstream: usize,
}

#[derive(Debug, Clone)]
struct ReferenceSyntenyContext {
    genome_idx: GenomeIdx,
    context: SyntenyContext,
}

#[derive(Debug, Default)]
pub(crate) struct SyntenyReferenceCache {
    contexts: HashMap<(ClusterIdx, usize), Vec<ReferenceSyntenyContext>>,
}

fn weighted_side_score(target: &[ClusterIdx], reference: &[ClusterIdx]) -> (f64, usize) {
    if target.is_empty() || reference.is_empty() {
        return (0.0, 0);
    }

    let mut matched_weight = 0.0;
    let mut denominator = 0.0;
    let mut matches = 0usize;
    let mut used_reference = vec![false; reference.len()];
    for (target_pos, target_cluster) in target.iter().enumerate() {
        let target_weight = 1.0 / (target_pos + 1) as f64;
        denominator += target_weight;
        let best_reference = reference
            .iter()
            .enumerate()
            .filter(|(reference_pos, reference_cluster)| {
                !used_reference[*reference_pos] && *reference_cluster == target_cluster
            })
            .map(|(reference_pos, _)| (reference_pos, 1.0 / (reference_pos + 1) as f64))
            .max_by(|(_, a), (_, b)| a.total_cmp(b));
        if let Some((reference_pos, best_reference_weight)) = best_reference {
            used_reference[reference_pos] = true;
            matches += 1;
            matched_weight += target_weight * best_reference_weight;
        }
    }

    if denominator == 0.0 {
        (0.0, matches)
    } else {
        (matched_weight / denominator, matches)
    }
}

fn context_score(target: &SyntenyContext, reference: &SyntenyContext) -> (f64, usize, usize) {
    let (up_score, up_matches) = weighted_side_score(&target.upstream, &reference.upstream);
    let (down_score, down_matches) = weighted_side_score(&target.downstream, &reference.downstream);

    let comparable_sides =
        usize::from(!target.upstream.is_empty()) + usize::from(!target.downstream.is_empty());
    if comparable_sides == 0 {
        return (0.0, up_matches, down_matches);
    }

    let score = (if target.upstream.is_empty() {
        0.0
    } else {
        up_score
    } + if target.downstream.is_empty() {
        0.0
    } else {
        down_score
    }) / comparable_sides as f64;
    (score, up_matches, down_matches)
}

fn synteny_context(
    gff: &GffIndex,
    gene_to_cluster: &HashMap<String, ClusterIdx>,
    contig: &str,
    start: u64,
    end: u64,
    center_gene_id: Option<&str>,
    n_neighbors: usize,
) -> SyntenyContext {
    let mut context = SyntenyContext::default();
    let (upstream, downstream) =
        gff.context_window(contig, start, end, center_gene_id, n_neighbors);
    for gene in upstream {
        if let Some(cluster) = gff.cluster_for_gene(gene, gene_to_cluster) {
            context.upstream.push(cluster);
        }
    }
    for gene in downstream {
        if let Some(cluster) = gff.cluster_for_gene(gene, gene_to_cluster) {
            context.downstream.push(cluster);
        }
    }
    context
}

#[cfg(test)]
fn compatible_reference_support(
    target_context: &SyntenyContext,
    reference_contexts: &[SyntenyContext],
) -> usize {
    let reference_contexts = labeled_test_contexts(reference_contexts);
    assess_context_support(target_context, &reference_contexts, SYNTENY_CONTEXT_GENES).support
}

#[cfg(test)]
fn has_recurrent_context_support(
    target_context: &SyntenyContext,
    reference_contexts: &[SyntenyContext],
) -> bool {
    let reference_contexts = labeled_test_contexts(reference_contexts);
    assess_context_support(target_context, &reference_contexts, SYNTENY_CONTEXT_GENES).passes
}

#[cfg(test)]
fn labeled_test_contexts(reference_contexts: &[SyntenyContext]) -> Vec<ReferenceSyntenyContext> {
    reference_contexts
        .iter()
        .enumerate()
        .map(|(idx, context)| ReferenceSyntenyContext {
            genome_idx: GenomeIdx(idx),
            context: context.clone(),
        })
        .collect()
}

fn assess_context_support(
    target_context: &SyntenyContext,
    reference_contexts: &[ReferenceSyntenyContext],
    n_neighbors: usize,
) -> SyntenyAssessment {
    let mut support = 0usize;
    let mut support_score_sum = 0.0;
    let mut best_score = 0.0;
    let mut best_up = 0usize;
    let mut best_down = 0usize;
    let reference_genomes: std::collections::BTreeSet<GenomeIdx> = reference_contexts
        .iter()
        .map(|reference_context| reference_context.genome_idx)
        .collect();
    let mut support_genomes = std::collections::BTreeSet::new();

    for reference_context in reference_contexts {
        let (score, up_matches, down_matches) =
            context_score(target_context, &reference_context.context);
        if score > best_score {
            best_score = score;
            best_up = up_matches;
            best_down = down_matches;
        }
        if score >= SYNTENY_REFERENCE_SUPPORT_THRESHOLD
            && target_context.side_matches(&reference_context.context)
        {
            support += 1;
            support_score_sum += score;
            support_genomes.insert(reference_context.genome_idx);
        }
    }

    let score = if support > 0 {
        support_score_sum / support as f64
    } else {
        best_score
    };
    let passes = if reference_genomes.len() <= 3 {
        support > 0
    } else {
        support_genomes.len() >= 2
    };

    SyntenyAssessment {
        passes,
        score,
        support,
        reference_contexts: reference_contexts.len(),
        n_neighbors,
        matched_upstream: best_up,
        matched_downstream: best_down,
    }
}

fn better_synteny_assessment(
    current: Option<SyntenyAssessment>,
    candidate: SyntenyAssessment,
) -> Option<SyntenyAssessment> {
    let Some(current) = current else {
        return Some(candidate);
    };
    if candidate.passes != current.passes {
        return Some(if candidate.passes { candidate } else { current });
    }
    if candidate.score > current.score + f64::EPSILON {
        return Some(candidate);
    }
    if (candidate.score - current.score).abs() <= f64::EPSILON
        && candidate.n_neighbors < current.n_neighbors
    {
        return Some(candidate);
    }
    Some(current)
}

fn reference_contexts_for_cluster<'a>(
    cache: &'a mut SyntenyReferenceCache,
    ci: ClusterIdx,
    n_neighbors: usize,
    pa: &PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    gene_to_cluster_maps: &HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
    verbose: bool,
) -> &'a [ReferenceSyntenyContext] {
    cache
        .contexts
        .entry((ci, n_neighbors))
        .or_insert_with(|| {
            let mut reference_contexts = Vec::new();
            let mut ref_genomes: Vec<GenomeIdx> = pa.present_genomes(ci).iter().copied().collect();
            ref_genomes.sort_by_key(|gi| gi.val());
            for ref_gi in ref_genomes {
                let Some(ref_gff) = gff_indices.get(&ref_gi) else {
                    continue;
                };
                let Some(ref_g2c) = gene_to_cluster_maps.get(&ref_gi) else {
                    continue;
                };

                let ref_cell = pa.cell(ci, ref_gi);
                if ref_cell.is_empty() {
                    continue;
                }
                for ref_gene_id in ref_cell
                    .split(';')
                    .map(str::trim)
                    .filter(|id| !id.is_empty())
                {
                    let ref_gene = ref_gff.get_by_id(ref_gene_id).or_else(|| {
                        panaroo_gene_id_suffix_variants(ref_gene_id)
                            .into_iter()
                            .skip(1)
                            .find_map(|id| ref_gff.get_by_id(id))
                    });
                    let Some(ref_gene) = ref_gene else {
                        if verbose && n_neighbors == 1 {
                            log::debug!(
                                "Synteny: ref gene '{}' from CSV not found in GFF of genome idx={}",
                                ref_gene_id,
                                ref_gi.val()
                            );
                        }
                        continue;
                    };

                    let ref_context = synteny_context(
                        ref_gff,
                        ref_g2c,
                        &ref_gene.contig,
                        ref_gene.start,
                        ref_gene.end,
                        Some(ref_gene.gene_id.as_str()),
                        n_neighbors,
                    );

                    if ref_context.has_any_cluster() {
                        reference_contexts.push(ReferenceSyntenyContext {
                            genome_idx: ref_gi,
                            context: ref_context,
                        });
                    }
                }
            }
            reference_contexts
        })
        .as_slice()
}

#[cfg(test)]
fn validate_synteny(
    rescue: &mut RescueEvent,
    pa: &PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    gene_to_cluster_maps: &HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
    verbose: bool,
) {
    let mut cache = SyntenyReferenceCache::default();
    validate_synteny_with_cache(
        rescue,
        pa,
        gff_indices,
        gene_to_cluster_maps,
        &mut cache,
        verbose,
    );
}

pub(crate) fn validate_synteny_with_cache(
    rescue: &mut RescueEvent,
    pa: &PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    gene_to_cluster_maps: &HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
    cache: &mut SyntenyReferenceCache,
    verbose: bool,
) {
    let ci = rescue.cluster_idx;
    let gi = rescue.genome_idx;

    let Some(target_gff) = gff_indices.get(&gi) else {
        if verbose {
            log::debug!("Synteny skip: no GFF for target genome idx={}", gi.val());
        }
        return;
    };
    let Some(target_g2c) = gene_to_cluster_maps.get(&gi) else {
        return;
    };

    let mut best_assessment: Option<SyntenyAssessment> = None;
    let mut saw_target_context = false;
    let mut saw_reference_context = false;

    for n_neighbors in 1..=SYNTENY_MAX_CONTEXT_GENES {
        let target_context = synteny_context(
            target_gff,
            target_g2c,
            &rescue.contig,
            rescue.start,
            rescue.end,
            rescue.annotation_match.as_deref(),
            n_neighbors,
        );

        if !target_context.has_any_cluster() {
            continue;
        }
        saw_target_context = true;

        let reference_contexts = reference_contexts_for_cluster(
            cache,
            ci,
            n_neighbors,
            pa,
            gff_indices,
            gene_to_cluster_maps,
            verbose,
        );
        if reference_contexts.is_empty() {
            continue;
        }
        saw_reference_context = true;
        let assessment = assess_context_support(&target_context, reference_contexts, n_neighbors);
        best_assessment = better_synteny_assessment(best_assessment, assessment);
    }

    let assessment = if let Some(assessment) = best_assessment {
        assessment
    } else if saw_target_context && !saw_reference_context {
        SyntenyAssessment {
            passes: false,
            score: 0.0,
            support: 0,
            reference_contexts: 0,
            n_neighbors: SYNTENY_MAX_CONTEXT_GENES,
            matched_upstream: 0,
            matched_downstream: 0,
        }
    } else {
        if verbose {
            log::debug!(
                "Synteny skip: no flanking genes mapped to clusters for rescue cluster={} genome={} at {}:{}-{}",
                ci.val(), gi.val(), rescue.contig, rescue.start, rescue.end
            );
            let (target_up, target_down) =
                target_gff.flanking(&rescue.contig, rescue.start, rescue.end);
            if let Some(up) = target_up {
                log::debug!(
                    "  Upstream gene '{}' not found in any Panaroo cluster",
                    up.gene_id
                );
            } else {
                log::debug!("  No upstream gene found in GFF");
            }
            if let Some(down) = target_down {
                log::debug!(
                    "  Downstream gene '{}' not found in any Panaroo cluster",
                    down.gene_id
                );
            } else {
                log::debug!("  No downstream gene found in GFF");
            }
        }
        return;
    };

    rescue.synteny_ok = Some(assessment.passes);
    rescue.synteny_evidence = SyntenyEvidence {
        score: Some(assessment.score),
        support: assessment.support,
        reference_contexts: assessment.reference_contexts,
        n_neighbors: assessment.n_neighbors,
        matched_upstream: assessment.matched_upstream,
        matched_downstream: assessment.matched_downstream,
    };

    if assessment.passes {
        rescue.confidence *= 1.2;
    } else {
        rescue.confidence *= 0.6;
    }
}

#[cfg(test)]
mod tests;
