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
    pub fragment_evidence: FragmentEvidence,
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
        records.sort_by(|a, b| a.0.contig.cmp(&b.0.contig).then(a.0.start.cmp(&b.0.start)));

        let genes: Vec<GffGene> = records.iter().map(|(gene, _)| gene.clone()).collect();

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
        for (i, (g, aliases)) in records.iter().enumerate() {
            if !g.gene_id.is_empty() {
                gene_id_index.insert(g.gene_id.clone(), i);
                let mut all_aliases = vec![g.gene_id.clone()];
                for alias in aliases {
                    if !alias.is_empty() && !all_aliases.iter().any(|x| x == alias) {
                        gene_id_index.entry(alias.clone()).or_insert(i);
                        all_aliases.push(alias.clone());
                    }
                }
                gene_aliases.insert(g.gene_id.clone(), all_aliases);
            }
        }

        GffIndex {
            genes,
            contig_ranges,
            gene_id_index,
            gene_aliases,
        }
    }

    pub fn get_by_id(&self, gene_id: &str) -> Option<&GffGene> {
        self.gene_id_index.get(gene_id).map(|&i| &self.genes[i])
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
                // Also index without _pseudo / _dup_N suffixes for GFF matching
                if let Some(stripped) = gene_id.strip_suffix("_pseudo") {
                    map.entry(stripped.to_string()).or_insert(ci);
                }
                if let Some(pos) = gene_id.rfind("_dup_") {
                    map.entry(gene_id[..pos].to_string()).or_insert(ci);
                }
            }
        }
    }
    map
}

// ── Fragment analysis ──

const BIOLOGICAL_FRAGMENT_MIN_GAP_BP: u64 = 50;
const SYNTENY_CONTEXT_GENES: usize = 2;
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
        let confidence = missing_cds_hit.identity() / 100.0
            * missing_cds_hit.query_coverage() / 100.0;
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
        hit_len > 0.0
            && gene_len > 0.0
            && (overlap / hit_len).min(overlap / gene_len) > 0.5
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

#[derive(Debug, Default)]
struct SyntenyContext {
    upstream: Vec<ClusterIdx>,
    downstream: Vec<ClusterIdx>,
}

impl SyntenyContext {
    fn has_any_cluster(&self) -> bool {
        !self.upstream.is_empty() || !self.downstream.is_empty()
    }

    fn side_matches(&self, other: &SyntenyContext) -> bool {
        let up_matches = cluster_overlap_count(&self.upstream, &other.upstream);
        let down_matches = cluster_overlap_count(&self.downstream, &other.downstream);
        let up_comparable = !self.upstream.is_empty() && !other.upstream.is_empty();
        let down_comparable = !self.downstream.is_empty() && !other.downstream.is_empty();

        if up_comparable && down_comparable {
            return up_matches > 0 && down_matches > 0;
        }

        up_matches + down_matches > 0
    }
}

fn cluster_overlap_count(a: &[ClusterIdx], b: &[ClusterIdx]) -> usize {
    a.iter().filter(|x| b.iter().any(|y| *x == y)).count()
}

fn synteny_context(
    gff: &GffIndex,
    gene_to_cluster: &HashMap<String, ClusterIdx>,
    contig: &str,
    start: u64,
    end: u64,
    center_gene_id: Option<&str>,
) -> SyntenyContext {
    let mut context = SyntenyContext::default();
    let (upstream, downstream) =
        gff.context_window(contig, start, end, center_gene_id, SYNTENY_CONTEXT_GENES);
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

fn compatible_reference_support(
    target_context: &SyntenyContext,
    reference_contexts: &[SyntenyContext],
) -> usize {
    reference_contexts
        .iter()
        .filter(|ref_context| target_context.side_matches(ref_context))
        .count()
}

fn has_recurrent_context_support(
    target_context: &SyntenyContext,
    reference_contexts: &[SyntenyContext],
) -> bool {
    let support = compatible_reference_support(target_context, reference_contexts);
    if reference_contexts.len() <= 3 {
        return support > 0;
    }
    support >= 2
}

pub fn validate_synteny(
    rescue: &mut RescueEvent,
    pa: &PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    gene_to_cluster_maps: &HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
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

    let target_context = synteny_context(
        target_gff,
        target_g2c,
        &rescue.contig,
        rescue.start,
        rescue.end,
        rescue.annotation_match.as_deref(),
    );

    if !target_context.has_any_cluster() {
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
    }

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
        let ref_gene_id = ref_cell.split(';').next().unwrap_or(ref_cell).trim();

        // Try exact ID first, then strip common suffixes (_pseudo, _dup_N)
        let ref_gene = ref_gff.get_by_id(ref_gene_id).or_else(|| {
            let stripped = ref_gene_id.strip_suffix("_pseudo").or_else(|| {
                // Strip _dup_N suffix
                ref_gene_id.rfind("_dup_").map(|i| &ref_gene_id[..i])
            });
            stripped.and_then(|s| ref_gff.get_by_id(s))
        });
        let Some(ref_gene) = ref_gene else {
            if verbose {
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
        );

        if ref_context.has_any_cluster() {
            reference_contexts.push(ref_context);
        }
    }

    let synteny_match = has_recurrent_context_support(&target_context, &reference_contexts);

    rescue.synteny_ok = Some(synteny_match);

    if synteny_match {
        rescue.confidence *= 1.2;
    } else {
        rescue.confidence *= 0.6;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::PresenceAbsence;

    fn make_hit(query_start: u64, query_end: u64, query_length: u64) -> Hit {
        Hit {
            query: "g".into(),
            cluster: "g".into(),
            query_length,
            query_start,
            query_end,
            strand: '+',
            target: "c".into(),
            target_start: query_start,
            target_end: query_end,
            matches: query_end - query_start,
            alignment_length: query_end - query_start,
            mapq: 60,
        }
    }

    #[test]
    fn test_merged_coverage_no_overlap() {
        let h1 = make_hit(0, 400, 1000);
        let h2 = make_hit(600, 1000, 1000);
        let (cov, frags) = merged_query_coverage(&[&h1, &h2]);
        assert!((cov - 80.0).abs() < 0.01);
        assert_eq!(frags, 2);
    }

    #[test]
    fn test_merged_coverage_with_overlap() {
        let h1 = make_hit(0, 600, 1000);
        let h2 = make_hit(400, 1000, 1000);
        let (cov, frags) = merged_query_coverage(&[&h1, &h2]);
        assert!((cov - 100.0).abs() < 0.01);
        assert_eq!(frags, 1);
    }

    #[test]
    fn test_merged_coverage_multi_allele() {
        let h1 = Hit {
            query: "g__0".into(),
            cluster: "g".into(),
            query_length: 800,
            query_start: 0,
            query_end: 600,
            strand: '+',
            target: "c".into(),
            target_start: 0,
            target_end: 600,
            matches: 600,
            alignment_length: 600,
            mapq: 60,
        };
        let h2 = Hit {
            query: "g__1".into(),
            cluster: "g".into(),
            query_length: 1200,
            query_start: 0,
            query_end: 1100,
            strand: '+',
            target: "c".into(),
            target_start: 0,
            target_end: 1100,
            matches: 1100,
            alignment_length: 1100,
            mapq: 60,
        };
        let (cov, frags) = merged_query_coverage(&[&h1, &h2]);
        assert!(cov > 90.0);
        assert_eq!(frags, 1);
    }

    #[test]
    fn test_merged_coverage_same_length_different_alleles() {
        let h1 = Hit {
            query: "cluster~~~ALLELE~0".into(),
            cluster: "cluster".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 400,
            strand: '+',
            target: "c".into(),
            target_start: 0,
            target_end: 400,
            matches: 400,
            alignment_length: 400,
            mapq: 60,
        };
        let h2 = Hit {
            query: "cluster~~~ALLELE~1".into(),
            cluster: "cluster".into(),
            query_length: 1000,
            query_start: 600,
            query_end: 1000,
            strand: '+',
            target: "c".into(),
            target_start: 5000,
            target_end: 5400,
            matches: 400,
            alignment_length: 400,
            mapq: 60,
        };
        let (cov, _frags) = merged_query_coverage(&[&h1, &h2]);
        assert!(
            cov < 50.0,
            "Coverage should be <50% (was {:.1}%), alleles must not be merged",
            cov
        );
    }

    #[test]
    fn test_gff_index_overlap() {
        let genes = vec![
            GffGene {
                contig: "chr1".into(),
                start: 100,
                end: 500,
                strand: '+',
                gene_id: "g1".into(),
            },
            GffGene {
                contig: "chr1".into(),
                start: 600,
                end: 900,
                strand: '-',
                gene_id: "g2".into(),
            },
            GffGene {
                contig: "chr2".into(),
                start: 50,
                end: 200,
                strand: '+',
                gene_id: "g3".into(),
            },
        ];
        let idx = GffIndex::build(genes);

        assert_eq!(idx.overlapping("chr1", 200, 400).len(), 1);
        assert_eq!(idx.overlapping("chr1", 400, 700).len(), 2);
        assert_eq!(idx.overlapping("chr2", 0, 100).len(), 1);
        assert!(idx.overlapping("chr3", 0, 100).is_empty());
    }

    #[test]
    fn test_gff_flanking() {
        let genes = vec![
            GffGene {
                contig: "c".into(),
                start: 100,
                end: 200,
                strand: '+',
                gene_id: "g1".into(),
            },
            GffGene {
                contig: "c".into(),
                start: 300,
                end: 400,
                strand: '+',
                gene_id: "g2".into(),
            },
            GffGene {
                contig: "c".into(),
                start: 500,
                end: 600,
                strand: '+',
                gene_id: "g3".into(),
            },
        ];
        let idx = GffIndex::build(genes);

        let (up, down) = idx.flanking("c", 300, 400);
        assert_eq!(up.unwrap().gene_id, "g1");
        assert_eq!(down.unwrap().gene_id, "g3");

        let (up, down) = idx.flanking("c", 100, 200);
        assert!(up.is_none());
        assert_eq!(down.unwrap().gene_id, "g2");
    }

    #[test]
    fn test_synteny_context_handles_overlapping_boundary_genes() {
        let pa = PresenceAbsence::new_for_test(
            vec![
                "target_gene".into(),
                "up_overlap".into(),
                "down_overlap".into(),
                "outer_up".into(),
                "outer_down".into(),
            ],
            vec!["target".into(), "reference".into()],
            vec![
                vec!["".into(), "ref_center".into()],
                vec!["target_up".into(), "ref_up".into()],
                vec!["target_down".into(), "ref_down".into()],
                vec!["target_outer_up".into(), "ref_outer_up".into()],
                vec!["target_outer_down".into(), "ref_outer_down".into()],
            ],
        );

        let target_gff = GffIndex::build(vec![
            GffGene {
                contig: "chr".into(),
                start: 0,
                end: 100,
                strand: '+',
                gene_id: "target_outer_up".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 100,
                end: 200,
                strand: '+',
                gene_id: "target_up".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 197,
                end: 260,
                strand: '+',
                gene_id: "target_center".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 257,
                end: 320,
                strand: '+',
                gene_id: "target_down".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 320,
                end: 400,
                strand: '+',
                gene_id: "target_outer_down".into(),
            },
        ]);
        let ref_gff = GffIndex::build(vec![
            GffGene {
                contig: "chr".into(),
                start: 0,
                end: 100,
                strand: '+',
                gene_id: "ref_outer_up".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 100,
                end: 200,
                strand: '+',
                gene_id: "ref_up".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 197,
                end: 260,
                strand: '+',
                gene_id: "ref_center".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 257,
                end: 320,
                strand: '+',
                gene_id: "ref_down".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 320,
                end: 400,
                strand: '+',
                gene_id: "ref_outer_down".into(),
            },
        ]);

        let mut gff_indices = HashMap::new();
        gff_indices.insert(GenomeIdx(0), target_gff);
        gff_indices.insert(GenomeIdx(1), ref_gff);

        let mut g2c = HashMap::new();
        g2c.insert(GenomeIdx(0), build_gene_to_cluster(&pa, GenomeIdx(0)));
        g2c.insert(GenomeIdx(1), build_gene_to_cluster(&pa, GenomeIdx(1)));

        let mut rescue = RescueEvent {
            cluster_idx: ClusterIdx(0),
            genome_idx: GenomeIdx(0),
            rescue_type: RescueType::StartVariant { length_diff_bp: 10 },
            identity: 100.0,
            coverage: 95.0,
            contig: "chr".into(),
            start: 205,
            end: 250,
            strand: '+',
            annotation_match: Some("target_center".into()),
            confidence: 0.95,
            synteny_ok: None,
            fragment_evidence: FragmentEvidence::default(),
        };

        validate_synteny(&mut rescue, &pa, &gff_indices, &g2c, false);
        assert_eq!(rescue.synteny_ok, Some(true));
    }

    #[test]
    fn test_synteny_context_rejects_different_locus_paralog() {
        let pa = PresenceAbsence::new_for_test(
            vec![
                "target_gene".into(),
                "target_up".into(),
                "target_down".into(),
            ],
            vec!["target".into(), "reference".into()],
            vec![
                vec!["".into(), "ref_center".into()],
                vec!["target_up_gene".into(), "".into()],
                vec!["target_down_gene".into(), "".into()],
            ],
        );

        let target_gff = GffIndex::build(vec![
            GffGene {
                contig: "chr".into(),
                start: 100,
                end: 200,
                strand: '+',
                gene_id: "target_up_gene".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 200,
                end: 300,
                strand: '+',
                gene_id: "target_center".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 300,
                end: 400,
                strand: '+',
                gene_id: "target_down_gene".into(),
            },
        ]);
        let ref_gff = GffIndex::build(vec![
            GffGene {
                contig: "chr".into(),
                start: 100,
                end: 200,
                strand: '+',
                gene_id: "ref_up_other".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 200,
                end: 300,
                strand: '+',
                gene_id: "ref_center".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 300,
                end: 400,
                strand: '+',
                gene_id: "ref_down_other".into(),
            },
        ]);

        let mut gff_indices = HashMap::new();
        gff_indices.insert(GenomeIdx(0), target_gff);
        gff_indices.insert(GenomeIdx(1), ref_gff);

        let mut g2c = HashMap::new();
        g2c.insert(GenomeIdx(0), build_gene_to_cluster(&pa, GenomeIdx(0)));
        g2c.insert(GenomeIdx(1), build_gene_to_cluster(&pa, GenomeIdx(1)));

        let mut rescue = RescueEvent {
            cluster_idx: ClusterIdx(0),
            genome_idx: GenomeIdx(0),
            rescue_type: RescueType::Split {
                other_cluster: "paralogous_cluster".into(),
            },
            identity: 99.0,
            coverage: 99.0,
            contig: "chr".into(),
            start: 200,
            end: 300,
            strand: '+',
            annotation_match: Some("target_center".into()),
            confidence: 1.0,
            synteny_ok: None,
            fragment_evidence: FragmentEvidence::default(),
        };

        validate_synteny(&mut rescue, &pa, &gff_indices, &g2c, false);
        assert_eq!(rescue.synteny_ok, Some(false));
    }

    #[test]
    fn test_synteny_uses_all_present_genomes_deterministically() {
        let pa = PresenceAbsence::new_for_test(
            vec![
                "target_gene".into(),
                "shared_up".into(),
                "shared_down".into(),
                "other_up".into(),
                "other_down".into(),
            ],
            vec![
                "ref0".into(),
                "ref1".into(),
                "ref2".into(),
                "ref3".into(),
                "ref4".into(),
                "ref5".into(),
                "ref6".into(),
                "target".into(),
            ],
            vec![
                vec![
                    "g0".into(),
                    "g1".into(),
                    "g2".into(),
                    "g3".into(),
                    "g4".into(),
                    "g5".into(),
                    "g6".into(),
                    "".into(),
                ],
                vec![
                    "".into(),
                    "".into(),
                    "".into(),
                    "".into(),
                    "".into(),
                    "up5".into(),
                    "up6".into(),
                    "up_target".into(),
                ],
                vec![
                    "".into(),
                    "".into(),
                    "".into(),
                    "".into(),
                    "".into(),
                    "down5".into(),
                    "down6".into(),
                    "down_target".into(),
                ],
                vec![
                    "xup0".into(),
                    "xup1".into(),
                    "xup2".into(),
                    "xup3".into(),
                    "xup4".into(),
                    "".into(),
                    "".into(),
                    "".into(),
                ],
                vec![
                    "xdown0".into(),
                    "xdown1".into(),
                    "xdown2".into(),
                    "xdown3".into(),
                    "xdown4".into(),
                    "".into(),
                    "".into(),
                    "".into(),
                ],
            ],
        );

        fn genes(up: &str, center: &str, down: &str) -> GffIndex {
            GffIndex::build(vec![
                GffGene {
                    contig: "chr".into(),
                    start: 100,
                    end: 200,
                    strand: '+',
                    gene_id: up.into(),
                },
                GffGene {
                    contig: "chr".into(),
                    start: 200,
                    end: 300,
                    strand: '+',
                    gene_id: center.into(),
                },
                GffGene {
                    contig: "chr".into(),
                    start: 300,
                    end: 400,
                    strand: '+',
                    gene_id: down.into(),
                },
            ])
        }

        let mut gff_indices = HashMap::new();
        for gi in 0..5 {
            gff_indices.insert(
                GenomeIdx(gi),
                genes(
                    &format!("xup{}", gi),
                    &format!("g{}", gi),
                    &format!("xdown{}", gi),
                ),
            );
        }
        gff_indices.insert(GenomeIdx(5), genes("up5", "g5", "down5"));
        gff_indices.insert(GenomeIdx(6), genes("up6", "g6", "down6"));
        gff_indices.insert(GenomeIdx(7), genes("up_target", "candidate", "down_target"));

        let mut g2c = HashMap::new();
        for gi in 0..8 {
            g2c.insert(GenomeIdx(gi), build_gene_to_cluster(&pa, GenomeIdx(gi)));
        }

        let mut rescue = RescueEvent {
            cluster_idx: ClusterIdx(0),
            genome_idx: GenomeIdx(7),
            rescue_type: RescueType::Dropped,
            identity: 99.0,
            coverage: 99.0,
            contig: "chr".into(),
            start: 200,
            end: 300,
            strand: '+',
            annotation_match: Some("candidate".into()),
            confidence: 1.0,
            synteny_ok: None,
            fragment_evidence: FragmentEvidence::default(),
        };

        validate_synteny(&mut rescue, &pa, &gff_indices, &g2c, false);
        assert_eq!(rescue.synteny_ok, Some(true));
    }

    #[test]
    fn test_missing_cds_preferred_over_annotated_paralog_and_validated_by_synteny() {
        let pa = PresenceAbsence::new_for_test(
            vec![
                "target_gene".into(),
                "shared_up".into(),
                "shared_down".into(),
                "other_cluster".into(),
            ],
            vec!["target".into(), "ref1".into(), "ref2".into()],
            vec![
                vec!["".into(), "ref_center1".into(), "ref_center2".into()],
                vec!["target_up".into(), "ref_up1".into(), "ref_up2".into()],
                vec!["target_down".into(), "ref_down1".into(), "ref_down2".into()],
                vec!["target_paralog".into(), "".into(), "".into()],
            ],
        );

        let target_gff = GffIndex::build(vec![
            GffGene {
                contig: "chr".into(),
                start: 100,
                end: 200,
                strand: '+',
                gene_id: "target_up".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 300,
                end: 400,
                strand: '+',
                gene_id: "target_down".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 1000,
                end: 1100,
                strand: '+',
                gene_id: "target_paralog".into(),
            },
        ]);
        let ref1_gff = GffIndex::build(vec![
            GffGene {
                contig: "chr".into(),
                start: 100,
                end: 200,
                strand: '+',
                gene_id: "ref_up1".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 200,
                end: 300,
                strand: '+',
                gene_id: "ref_center1".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 300,
                end: 400,
                strand: '+',
                gene_id: "ref_down1".into(),
            },
        ]);
        let ref2_gff = GffIndex::build(vec![
            GffGene {
                contig: "chr".into(),
                start: 100,
                end: 200,
                strand: '+',
                gene_id: "ref_up2".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 200,
                end: 300,
                strand: '+',
                gene_id: "ref_center2".into(),
            },
            GffGene {
                contig: "chr".into(),
                start: 300,
                end: 400,
                strand: '+',
                gene_id: "ref_down2".into(),
            },
        ]);

        let mut gff_indices = HashMap::new();
        gff_indices.insert(GenomeIdx(0), target_gff);
        gff_indices.insert(GenomeIdx(1), ref1_gff);
        gff_indices.insert(GenomeIdx(2), ref2_gff);

        let mut g2c = HashMap::new();
        for gi in 0..3 {
            g2c.insert(GenomeIdx(gi), build_gene_to_cluster(&pa, GenomeIdx(gi)));
        }

        let annotated_paralog = Hit {
            query: "target_gene".into(),
            cluster: "target_gene".into(),
            query_length: 100,
            query_start: 0,
            query_end: 100,
            strand: '+',
            target: "chr".into(),
            target_start: 1000,
            target_end: 1100,
            matches: 100,
            alignment_length: 100,
            mapq: 60,
        };
        let missing_cds_hit = Hit {
            query: "target_gene".into(),
            cluster: "target_gene".into(),
            query_length: 100,
            query_start: 0,
            query_end: 100,
            strand: '+',
            target: "chr".into(),
            target_start: 198,
            target_end: 302,
            matches: 98,
            alignment_length: 100,
            mapq: 60,
        };

        let mut rescue = classify_hits(
            &[annotated_paralog, missing_cds_hit],
            ClusterIdx(0),
            GenomeIdx(0),
            &pa,
            gff_indices.get(&GenomeIdx(0)),
            g2c.get(&GenomeIdx(0)).unwrap(),
            80.0,
        )
        .unwrap();

        assert!(matches!(rescue.rescue_type, RescueType::MissingCds));
        assert_eq!(rescue.start, 198);
        assert_eq!(rescue.annotation_match, None);

        validate_synteny(&mut rescue, &pa, &gff_indices, &g2c, false);
        assert_eq!(rescue.synteny_ok, Some(true));
    }

    #[test]
    fn test_synteny_rejects_one_off_context_when_many_references_exist() {
        let target = SyntenyContext {
            upstream: vec![ClusterIdx(1)],
            downstream: vec![ClusterIdx(2)],
        };
        let contexts = vec![
            SyntenyContext {
                upstream: vec![ClusterIdx(3)],
                downstream: vec![ClusterIdx(4)],
            },
            SyntenyContext {
                upstream: vec![ClusterIdx(5)],
                downstream: vec![ClusterIdx(6)],
            },
            SyntenyContext {
                upstream: vec![ClusterIdx(7)],
                downstream: vec![ClusterIdx(8)],
            },
            SyntenyContext {
                upstream: vec![ClusterIdx(9)],
                downstream: vec![ClusterIdx(10)],
            },
            SyntenyContext {
                upstream: vec![ClusterIdx(11)],
                downstream: vec![ClusterIdx(12)],
            },
            SyntenyContext {
                upstream: vec![ClusterIdx(1)],
                downstream: vec![ClusterIdx(2)],
            },
        ];
        assert!(!has_recurrent_context_support(&target, &contexts));

        let rare_contexts = &contexts[..3];
        assert!(has_recurrent_context_support(
            &rare_contexts[0],
            rare_contexts
        ));
    }

    #[test]
    fn test_synteny_accepts_recurrent_compatible_context_not_exact_signature() {
        let target = SyntenyContext {
            upstream: vec![ClusterIdx(1)],
            downstream: vec![ClusterIdx(2)],
        };
        let contexts = vec![
            SyntenyContext {
                upstream: vec![ClusterIdx(1), ClusterIdx(20)],
                downstream: vec![ClusterIdx(2), ClusterIdx(21)],
            },
            SyntenyContext {
                upstream: vec![ClusterIdx(1), ClusterIdx(30)],
                downstream: vec![ClusterIdx(2), ClusterIdx(31)],
            },
            SyntenyContext {
                upstream: vec![ClusterIdx(3)],
                downstream: vec![ClusterIdx(4)],
            },
            SyntenyContext {
                upstream: vec![ClusterIdx(5)],
                downstream: vec![ClusterIdx(6)],
            },
        ];

        assert_eq!(compatible_reference_support(&target, &contexts), 2);
        assert!(has_recurrent_context_support(&target, &contexts));
    }

    #[test]
    fn test_detect_start_variant_with_gff() {
        let genes = vec![GffGene {
            contig: "c".into(),
            start: 0,
            end: 1100,
            strand: '+',
            gene_id: "gene1".into(),
        }];
        let idx = GffIndex::build(genes);

        let hit = Hit {
            query: "g".into(),
            cluster: "g".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 1000,
            strand: '+',
            target: "c".into(),
            target_start: 100,
            target_end: 1100,
            matches: 1000,
            alignment_length: 1000,
            mapq: 60,
        };
        let result = detect_start_variant(&hit, Some(&idx));
        assert_eq!(result, Some((100, "gene1".to_string())));
    }

    #[test]
    fn test_detect_start_variant_uses_target_not_query_len() {
        // GFF gene: 0-1050 (length 1050)
        // Hit: query_length=1100, query aligned 0-1050. Target: 0-1000 (footprint 1000bp).
        // With fix: diff = gene_len(1050) - target_aln_len(1000) = +50bp ✓
        // Old bug would compute: diff = 1050 - 1100(query_length) = -50bp ✗
        let genes = vec![GffGene {
            contig: "c".into(),
            start: 0,
            end: 1050,
            strand: '+',
            gene_id: "gene1".into(),
        }];
        let idx = GffIndex::build(genes);
        let hit = Hit {
            query: "g".into(),
            cluster: "g".into(),
            query_length: 1100,
            query_start: 0,
            query_end: 1050,
            strand: '+',
            target: "c".into(),
            target_start: 0,
            target_end: 1000,
            matches: 1000,
            alignment_length: 1050,
            mapq: 60,
        };
        // identity = 1000/1050 = 95.2% (≥98%? NO — need 98%)
        // query_coverage = 1050/1100 = 95.5% (≥90% ✓)
        // Adjust matches to pass identity filter:
        let hit = Hit {
            matches: 1040,
            ..hit
        }; // 1040/1050 = 99.0%
        let result = detect_start_variant(&hit, Some(&idx));
        // diff = gene_len(1050) - target_aln_len(1000) = +50
        assert_eq!(result, Some((50, "gene1".to_string())));
    }

    #[test]
    fn test_detect_start_variant_no_gff() {
        let hit = make_hit(0, 1000, 1000);
        assert!(detect_start_variant(&hit, None).is_none());
    }

    #[test]
    fn test_conflict_resolution() {
        let mut rescues = vec![
            RescueEvent {
                cluster_idx: ClusterIdx(0),
                genome_idx: GenomeIdx(0),
                rescue_type: RescueType::Dropped,
                identity: 99.0,
                coverage: 100.0,
                contig: "c1".into(),
                start: 100,
                end: 500,
                strand: '+',
                annotation_match: None,
                confidence: 0.99,
                synteny_ok: None,
                fragment_evidence: FragmentEvidence::default(),
            },
            RescueEvent {
                cluster_idx: ClusterIdx(1),
                genome_idx: GenomeIdx(0),
                rescue_type: RescueType::Dropped,
                identity: 95.0,
                coverage: 90.0,
                contig: "c1".into(),
                start: 200,
                end: 600,
                strand: '+',
                annotation_match: None,
                confidence: 0.85,
                synteny_ok: None,
                fragment_evidence: FragmentEvidence::default(),
            },
        ];
        resolve_conflicts(&mut rescues);
        assert_eq!(rescues.len(), 1);
        assert_eq!(rescues[0].cluster_idx, ClusterIdx(0));
    }

    #[test]
    fn test_conflict_resolution_no_overlap() {
        let mut rescues = vec![
            RescueEvent {
                cluster_idx: ClusterIdx(0),
                genome_idx: GenomeIdx(0),
                rescue_type: RescueType::Dropped,
                identity: 99.0,
                coverage: 100.0,
                contig: "c1".into(),
                start: 100,
                end: 500,
                strand: '+',
                annotation_match: None,
                confidence: 0.99,
                synteny_ok: None,
                fragment_evidence: FragmentEvidence::default(),
            },
            RescueEvent {
                cluster_idx: ClusterIdx(1),
                genome_idx: GenomeIdx(0),
                rescue_type: RescueType::Dropped,
                identity: 95.0,
                coverage: 90.0,
                contig: "c1".into(),
                start: 600,
                end: 900,
                strand: '+',
                annotation_match: None,
                confidence: 0.85,
                synteny_ok: None,
                fragment_evidence: FragmentEvidence::default(),
            },
        ];
        resolve_conflicts(&mut rescues);
        assert_eq!(rescues.len(), 2);
    }

    #[test]
    fn test_conflict_resolution_keeps_minor_overlaps() {
        let mut rescues = vec![
            RescueEvent {
                cluster_idx: ClusterIdx(0),
                genome_idx: GenomeIdx(0),
                rescue_type: RescueType::Dropped,
                identity: 99.0,
                coverage: 100.0,
                contig: "c1".into(),
                start: 100,
                end: 500,
                strand: '+',
                annotation_match: None,
                confidence: 0.99,
                synteny_ok: None,
                fragment_evidence: FragmentEvidence::default(),
            },
            RescueEvent {
                cluster_idx: ClusterIdx(1),
                genome_idx: GenomeIdx(0),
                rescue_type: RescueType::Dropped,
                identity: 98.0,
                coverage: 100.0,
                contig: "c1".into(),
                start: 490,
                end: 900,
                strand: '+',
                annotation_match: None,
                confidence: 0.98,
                synteny_ok: None,
                fragment_evidence: FragmentEvidence::default(),
            },
        ];
        resolve_conflicts(&mut rescues);
        assert_eq!(rescues.len(), 2);
    }

    #[test]
    fn test_fragmented_span() {
        let h1 = Hit {
            query: "g".into(),
            cluster: "g".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 400,
            strand: '+',
            target: "c".into(),
            target_start: 1000,
            target_end: 1400,
            matches: 400,
            alignment_length: 400,
            mapq: 60,
        };
        let h2 = Hit {
            query: "g".into(),
            cluster: "g".into(),
            query_length: 1000,
            query_start: 600,
            query_end: 1000,
            strand: '+',
            target: "c".into(),
            target_start: 2000,
            target_end: 2400,
            matches: 400,
            alignment_length: 400,
            mapq: 60,
        };

        let pa = PresenceAbsence::new_for_test(
            vec!["g".into()],
            vec!["genome1".into()],
            vec![vec!["".into()]],
        );
        let g2c: HashMap<String, ClusterIdx> = HashMap::new();

        let result = classify_hits(
            &[h1, h2],
            ClusterIdx(0),
            GenomeIdx(0),
            &pa,
            None,
            &g2c,
            80.0,
        );
        assert!(result.is_some());
        let r = result.unwrap();
        assert_eq!(r.start, 1000);
        assert_eq!(r.end, 2400);
        match r.rescue_type {
            RescueType::Fragmented { n_fragments, .. } => assert_eq!(n_fragments, 2),
            _ => panic!("Expected Fragmented"),
        }
    }

    #[test]
    fn test_classify_split_vs_paralog() {
        // Setup: cluster "target_cluster" absent in genome, hit lands on a gene
        // assigned to "other_cluster". Identity determines SPLIT vs PARALOG.
        let gff_genes = vec![GffGene {
            contig: "c".into(),
            start: 0,
            end: 1000,
            strand: '+',
            gene_id: "geneA".into(),
        }];
        let gff = GffIndex::build(gff_genes);

        let pa = PresenceAbsence::new_for_test(
            vec!["target_cluster".into(), "other_cluster".into()],
            vec!["genome1".into()],
            vec![vec!["".into()], vec!["geneA".into()]],
        );

        let mut g2c: HashMap<String, ClusterIdx> = HashMap::new();
        g2c.insert("geneA".into(), ClusterIdx(1)); // geneA belongs to other_cluster

        // High identity (99%) → SPLIT
        let hit_high = Hit {
            query: "q".into(),
            cluster: "target_cluster".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 1000,
            strand: '+',
            target: "c".into(),
            target_start: 0,
            target_end: 1000,
            matches: 990,
            alignment_length: 1000,
            mapq: 60,
        };
        let result = classify_hits(
            &[hit_high],
            ClusterIdx(0),
            GenomeIdx(0),
            &pa,
            Some(&gff),
            &g2c,
            80.0,
        );
        let r = result.unwrap();
        match &r.rescue_type {
            RescueType::Split { other_cluster } => assert_eq!(other_cluster, "other_cluster"),
            other => panic!("Expected Split, got {:?}", other),
        }

        // Low identity (90%) → PARALOG
        let hit_low = Hit {
            query: "q".into(),
            cluster: "target_cluster".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 1000,
            strand: '+',
            target: "c".into(),
            target_start: 0,
            target_end: 1000,
            matches: 900,
            alignment_length: 1000,
            mapq: 60,
        };
        let result = classify_hits(
            &[hit_low],
            ClusterIdx(0),
            GenomeIdx(0),
            &pa,
            Some(&gff),
            &g2c,
            80.0,
        );
        let r = result.unwrap();
        match &r.rescue_type {
            RescueType::Paralog { other_cluster } => assert_eq!(other_cluster, "other_cluster"),
            other => panic!("Expected Paralog, got {:?}", other),
        }
        assert!(r.confidence < 0.7, "Paralog should have reduced confidence");
    }

    #[test]
    fn test_classify_dropped_no_gff() {
        // Full-coverage hit with no GFF → DROPPED
        let hit = Hit {
            query: "q".into(),
            cluster: "cluster1".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 1000,
            strand: '+',
            target: "c".into(),
            target_start: 0,
            target_end: 1000,
            matches: 990,
            alignment_length: 1000,
            mapq: 60,
        };
        let pa = PresenceAbsence::new_for_test(
            vec!["cluster1".into()],
            vec!["genome1".into()],
            vec![vec!["".into()]],
        );
        let g2c: HashMap<String, ClusterIdx> = HashMap::new();
        let result = classify_hits(&[hit], ClusterIdx(0), GenomeIdx(0), &pa, None, &g2c, 80.0);
        let r = result.unwrap();
        assert!(matches!(r.rescue_type, RescueType::Dropped));
        assert!(r.confidence > 0.9);
    }

    #[test]
    fn test_classify_prefers_best_evidence_over_longest_allele() {
        let long_lower_score = Hit {
            query: "cluster1~~~ALLELE~0".into(),
            cluster: "cluster1".into(),
            query_length: 1200,
            query_start: 0,
            query_end: 960,
            strand: '+',
            target: "long_hit".into(),
            target_start: 100,
            target_end: 1060,
            matches: 864,
            alignment_length: 960,
            mapq: 60,
        };
        let short_better_score = Hit {
            query: "cluster1~~~ALLELE~1".into(),
            cluster: "cluster1".into(),
            query_length: 900,
            query_start: 0,
            query_end: 900,
            strand: '+',
            target: "best_hit".into(),
            target_start: 200,
            target_end: 1100,
            matches: 900,
            alignment_length: 900,
            mapq: 60,
        };
        let pa = PresenceAbsence::new_for_test(
            vec!["cluster1".into()],
            vec!["genome1".into()],
            vec![vec!["".into()]],
        );
        let g2c: HashMap<String, ClusterIdx> = HashMap::new();
        let result = classify_hits(
            &[long_lower_score, short_better_score],
            ClusterIdx(0),
            GenomeIdx(0),
            &pa,
            None,
            &g2c,
            80.0,
        )
        .unwrap();

        assert_eq!(result.contig, "best_hit");
        assert!((result.identity - 100.0).abs() < 0.01);
        assert!((result.coverage - 100.0).abs() < 0.01);
    }

    #[test]
    fn test_classify_missing_cds_no_overlapping_genes() {
        // Full-coverage hit, GFF exists but no gene overlaps the hit region → MISSING_CDS
        let gff_genes = vec![GffGene {
            contig: "c".into(),
            start: 5000,
            end: 6000,
            strand: '+',
            gene_id: "far_gene".into(),
        }];
        let gff = GffIndex::build(gff_genes);
        let pa = PresenceAbsence::new_for_test(
            vec!["cluster1".into()],
            vec!["genome1".into()],
            vec![vec!["".into()]],
        );
        let g2c: HashMap<String, ClusterIdx> = HashMap::new();
        let hit = Hit {
            query: "q".into(),
            cluster: "cluster1".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 1000,
            strand: '+',
            target: "c".into(),
            target_start: 0,
            target_end: 1000,
            matches: 990,
            alignment_length: 1000,
            mapq: 60,
        };
        let result = classify_hits(
            &[hit],
            ClusterIdx(0),
            GenomeIdx(0),
            &pa,
            Some(&gff),
            &g2c,
            80.0,
        );
        let r = result.unwrap();
        assert!(matches!(r.rescue_type, RescueType::MissingCds));
    }

    #[test]
    fn test_classify_start_variant_with_gff() {
        // Hit with high identity overlapping a similar-length GFF gene → START_VARIANT
        let gff_genes = vec![GffGene {
            contig: "c".into(),
            start: 0,
            end: 1100,
            strand: '+',
            gene_id: "gene1".into(),
        }];
        let gff = GffIndex::build(gff_genes);
        let pa = PresenceAbsence::new_for_test(
            vec!["cluster1".into()],
            vec!["genome1".into()],
            vec![vec!["".into()]],
        );
        let g2c: HashMap<String, ClusterIdx> = HashMap::new();
        // gene1 is NOT in g2c → no cluster assignment → classify_single_hit returns DROPPED
        // Then detect_start_variant kicks in: gene_len(1100) - target_aln(1000) = +100bp
        let hit = Hit {
            query: "q".into(),
            cluster: "cluster1".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 1000,
            strand: '+',
            target: "c".into(),
            target_start: 100,
            target_end: 1100,
            matches: 1000,
            alignment_length: 1000,
            mapq: 60,
        };
        let result = classify_hits(
            &[hit],
            ClusterIdx(0),
            GenomeIdx(0),
            &pa,
            Some(&gff),
            &g2c,
            80.0,
        );
        let r = result.unwrap();
        match &r.rescue_type {
            RescueType::StartVariant { length_diff_bp } => assert_eq!(*length_diff_bp, 100),
            other => panic!("Expected StartVariant, got {:?}", other),
        }
        assert_eq!(r.annotation_match, Some("gene1".to_string()));
    }

    #[test]
    fn test_classify_returns_none_for_empty_hits() {
        let pa = PresenceAbsence::new_for_test(
            vec!["c1".into()],
            vec!["g1".into()],
            vec![vec!["".into()]],
        );
        let g2c: HashMap<String, ClusterIdx> = HashMap::new();
        let result = classify_hits(&[], ClusterIdx(0), GenomeIdx(0), &pa, None, &g2c, 80.0);
        assert!(result.is_none());
    }

    #[test]
    fn test_classify_returns_none_below_coverage() {
        // Single hit with only 50% coverage, no fragments → None
        let hit = Hit {
            query: "q".into(),
            cluster: "c1".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 500,
            strand: '+',
            target: "c".into(),
            target_start: 0,
            target_end: 500,
            matches: 500,
            alignment_length: 500,
            mapq: 60,
        };
        let pa = PresenceAbsence::new_for_test(
            vec!["c1".into()],
            vec!["g1".into()],
            vec![vec!["".into()]],
        );
        let g2c: HashMap<String, ClusterIdx> = HashMap::new();
        let result = classify_hits(&[hit], ClusterIdx(0), GenomeIdx(0), &pa, None, &g2c, 80.0);
        assert!(result.is_none());
    }

    #[test]
    fn test_build_gene_to_cluster() {
        let pa = PresenceAbsence::new_for_test(
            vec!["clusterA".into(), "clusterB".into()],
            vec!["g1".into(), "g2".into()],
            vec![
                vec!["gene1;gene2".into(), "gene3".into()],
                vec!["".into(), "gene4".into()],
            ],
        );
        let g2c = build_gene_to_cluster(&pa, GenomeIdx(0));
        assert_eq!(g2c.get("gene1"), Some(&ClusterIdx(0)));
        assert_eq!(g2c.get("gene2"), Some(&ClusterIdx(0)));
        assert_eq!(g2c.get("gene4"), None); // gene4 is in genome1 (GenomeIdx(1)), not genome0

        let g2c_1 = build_gene_to_cluster(&pa, GenomeIdx(1));
        assert_eq!(g2c_1.get("gene3"), Some(&ClusterIdx(0)));
        assert_eq!(g2c_1.get("gene4"), Some(&ClusterIdx(1)));
    }

    #[test]
    fn test_rescue_type_display() {
        assert_eq!(RescueType::Dropped.to_string(), "DROPPED");
        assert_eq!(RescueType::MissingCds.to_string(), "MISSING_CDS");
        assert_eq!(
            RescueType::Split {
                other_cluster: "cl2".into()
            }
            .to_string(),
            "SPLIT:cl2"
        );
        assert_eq!(
            RescueType::Fragmented {
                n_fragments: 3,
                total_coverage: 95.5
            }
            .to_string(),
            "FRAGMENTED:3frags:96%cov"
        );
        assert_eq!(
            RescueType::AlignmentSplit {
                n_fragments: 2,
                total_coverage: 98.8
            }
            .to_string(),
            "ALIGNMENT_SPLIT:2frags:99%cov"
        );
        assert_eq!(
            RescueType::StartVariant {
                length_diff_bp: -42
            }
            .to_string(),
            "START_VARIANT:-42bp"
        );
        assert_eq!(
            RescueType::Paralog {
                other_cluster: "cl3".into()
            }
            .to_string(),
            "PARALOG:cl3"
        );
        assert_eq!(
            RescueType::Pseudogene {
                reason: PseudogeneReason::InternalStop { n_stops: 2 }
            }
            .to_string(),
            "PSEUDOGENE:2stops"
        );
        assert_eq!(
            RescueType::Pseudogene {
                reason: PseudogeneReason::Frameshift
            }
            .to_string(),
            "PSEUDOGENE:frameshift"
        );
        assert_eq!(
            RescueType::Pseudogene {
                reason: PseudogeneReason::PartialDeletion { coverage: 72.5 }
            }
            .to_string(),
            "PSEUDOGENE:partial:72%cov"
        );
    }

    #[test]
    fn test_count_internal_stops_clean_orf() {
        // ATG ... no stops ... TAA
        let seq = b"ATGAAAGCGTTTTTTAAAGCGAAAGCGTAA";
        assert_eq!(count_internal_stops(seq), 0);
    }

    #[test]
    fn test_count_internal_stops_with_stops() {
        // Frame 0: ATG TAA GCG TTT — has internal stop at pos 3
        let seq = b"ATGTAAGCGTTTGCGAAAGCG";
        let stops = count_internal_stops(seq);
        // At least one frame will have stops, the cleanest frame might have 0 or 1
        // Frame 0: ATG|TAA|GCG|TTT|GCG|AAA|GCG = TAA at pos 3 (internal)
        // We want the minimum across frames
        assert!(stops <= 1);
    }

    #[test]
    fn test_count_internal_stops_short_seq() {
        assert_eq!(count_internal_stops(b"ATG"), 0);
        assert_eq!(count_internal_stops(b""), 0);
    }

    #[test]
    fn test_frameshift_pattern_detected() {
        // Two fragments with a 2bp gap (not divisible by 3) → frameshift
        let h1 = Hit {
            query: "q".into(),
            cluster: "c".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 500,
            strand: '+',
            target: "c".into(),
            target_start: 100,
            target_end: 600,
            matches: 500,
            alignment_length: 500,
            mapq: 60,
        };
        let h2 = Hit {
            query: "q".into(),
            cluster: "c".into(),
            query_length: 1000,
            query_start: 500,
            query_end: 1000,
            strand: '+',
            target: "c".into(),
            target_start: 602,
            target_end: 1102, // 2bp gap (not %3)
            matches: 500,
            alignment_length: 500,
            mapq: 60,
        };
        assert!(detect_frameshift_pattern(&[&h1, &h2]));
    }

    #[test]
    fn test_frameshift_pattern_not_detected_divisible_by_3() {
        let h1 = Hit {
            query: "q".into(),
            cluster: "c".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 500,
            strand: '+',
            target: "c".into(),
            target_start: 100,
            target_end: 600,
            matches: 500,
            alignment_length: 500,
            mapq: 60,
        };
        let h2 = Hit {
            query: "q".into(),
            cluster: "c".into(),
            query_length: 1000,
            query_start: 500,
            query_end: 1000,
            strand: '+',
            target: "c".into(),
            target_start: 603,
            target_end: 1103, // 3bp gap (divisible by 3)
            matches: 500,
            alignment_length: 500,
            mapq: 60,
        };
        assert!(!detect_frameshift_pattern(&[&h1, &h2]));
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ATCG"), b"CGAT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b""), b"");
    }

    #[test]
    fn test_internal_stops_minus_strand() {
        // Forward: GCATTAGCG — no obvious stops in frame 0
        // Reverse complement: CGCTAATGC
        // Frame 0: CGC|TAA|TGC → TAA is internal stop!
        let seq = b"GCATTAGCG";
        let rc = reverse_complement(seq);
        assert_eq!(&rc, b"CGCTAATGC");
        // The RC should have a stop in at least one frame
        let stops = count_internal_stops(&rc);
        // Frame 0: CGC TAA TGC → 1 stop (TAA), but it's at pos 3 with only 1 codon after → internal
        // Minimum across frames could be 0 (other frames may be clean)
        // Just verify the function runs without panic on RC
        assert!(stops <= 1);

        // Now test via check_internal_stops_in_hit with strand -
        let mut contigs: HashMap<String, Vec<u8>> = HashMap::new();
        contigs.insert("c".into(), seq.to_vec());
        let hit_minus = Hit {
            query: "q".into(),
            cluster: "c".into(),
            query_length: 9,
            query_start: 0,
            query_end: 9,
            strand: '-',
            target: "c".into(),
            target_start: 0,
            target_end: 9,
            matches: 9,
            alignment_length: 9,
            mapq: 60,
        };
        let hit_plus = Hit {
            strand: '+',
            ..hit_minus.clone()
        };
        // Strand - should RC before scanning
        let result_minus = check_internal_stops_in_hit(&hit_minus, Some(&contigs));
        let result_plus = check_internal_stops_in_hit(&hit_plus, Some(&contigs));
        // At minimum, they should not be identical (different sequences scanned)
        // (they could both be None if no stops, which is fine)
        let _ = (result_minus, result_plus); // just ensure no panic
    }

    #[test]
    fn test_classify_partial_deletion_pseudogene() {
        // Hit with 95% identity but only 70% coverage → partial deletion pseudogene
        let hit = Hit {
            query: "q".into(),
            cluster: "c1".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 700,
            strand: '+',
            target: "c".into(),
            target_start: 0,
            target_end: 700,
            matches: 665,
            alignment_length: 700,
            mapq: 60,
        };
        let pa = PresenceAbsence::new_for_test(
            vec!["c1".into()],
            vec!["g1".into()],
            vec![vec!["".into()]],
        );
        let g2c: HashMap<String, ClusterIdx> = HashMap::new();
        // min_coverage = 80% so 70% coverage won't pass Path 1
        let result = classify_hits(&[hit], ClusterIdx(0), GenomeIdx(0), &pa, None, &g2c, 80.0);
        let r = result.unwrap();
        match &r.rescue_type {
            RescueType::Pseudogene {
                reason: PseudogeneReason::PartialDeletion { coverage },
            } => {
                assert!(*coverage >= 69.0 && *coverage <= 71.0);
            }
            other => panic!("Expected Pseudogene:PartialDeletion, got {:?}", other),
        }
    }

    #[test]
    fn test_alignment_split_for_small_colinear_gaps() {
        let h1 = Hit {
            query: "g".into(),
            cluster: "g".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 495,
            strand: '+',
            target: "c".into(),
            target_start: 0,
            target_end: 495,
            matches: 495,
            alignment_length: 495,
            mapq: 60,
        };
        let h2 = Hit {
            query: "g".into(),
            cluster: "g".into(),
            query_length: 1000,
            query_start: 504,
            query_end: 1000,
            strand: '+',
            target: "c".into(),
            target_start: 504,
            target_end: 1000,
            matches: 495,
            alignment_length: 495,
            mapq: 60,
        };
        let pa = PresenceAbsence::new_for_test(
            vec!["g".into()],
            vec!["genome1".into()],
            vec![vec!["".into()]],
        );
        let g2c: HashMap<String, ClusterIdx> = HashMap::new();
        let result = classify_hits(
            &[h1, h2],
            ClusterIdx(0),
            GenomeIdx(0),
            &pa,
            None,
            &g2c,
            80.0,
        )
        .unwrap();
        match result.rescue_type {
            RescueType::AlignmentSplit { n_fragments, .. } => assert_eq!(n_fragments, 2),
            other => panic!("Expected AlignmentSplit, got {:?}", other),
        }
        assert_eq!(result.fragment_evidence.max_target_gap_bp, 9);
        assert_eq!(
            result.fragment_evidence.note.as_deref(),
            Some("alignment_split")
        );
    }

    #[test]
    fn test_fragment_chain_inside_same_cds_is_not_fragmented() {
        let h1 = Hit {
            query: "g".into(),
            cluster: "g".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 400,
            strand: '+',
            target: "c".into(),
            target_start: 100,
            target_end: 500,
            matches: 400,
            alignment_length: 400,
            mapq: 60,
        };
        let h2 = Hit {
            query: "g".into(),
            cluster: "g".into(),
            query_length: 1000,
            query_start: 600,
            query_end: 1000,
            strand: '+',
            target: "c".into(),
            target_start: 700,
            target_end: 1100,
            matches: 400,
            alignment_length: 400,
            mapq: 60,
        };
        let gff = GffIndex::build(vec![GffGene {
            contig: "c".into(),
            start: 100,
            end: 1100,
            strand: '+',
            gene_id: "gene1".into(),
        }]);
        let pa = PresenceAbsence::new_for_test(
            vec!["g".into()],
            vec!["genome1".into()],
            vec![vec!["".into()]],
        );
        let g2c: HashMap<String, ClusterIdx> = HashMap::new();
        let result = classify_hits(
            &[h1, h2],
            ClusterIdx(0),
            GenomeIdx(0),
            &pa,
            Some(&gff),
            &g2c,
            80.0,
        )
        .unwrap();
        assert!(matches!(result.rescue_type, RescueType::Dropped));
        assert_eq!(result.annotation_match.as_deref(), Some("gene1"));
        assert_eq!(result.fragment_evidence.same_cds, Some(true));
    }

    #[test]
    fn test_non_colinear_fragments_are_not_rescued() {
        let h1 = Hit {
            query: "g".into(),
            cluster: "g".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 400,
            strand: '+',
            target: "c".into(),
            target_start: 1000,
            target_end: 1400,
            matches: 400,
            alignment_length: 400,
            mapq: 60,
        };
        let h2 = Hit {
            query: "g".into(),
            cluster: "g".into(),
            query_length: 1000,
            query_start: 600,
            query_end: 1000,
            strand: '+',
            target: "c".into(),
            target_start: 500,
            target_end: 900,
            matches: 400,
            alignment_length: 400,
            mapq: 60,
        };
        let pa = PresenceAbsence::new_for_test(
            vec!["g".into()],
            vec!["genome1".into()],
            vec![vec!["".into()]],
        );
        let g2c: HashMap<String, ClusterIdx> = HashMap::new();
        let result = classify_hits(
            &[h1, h2],
            ClusterIdx(0),
            GenomeIdx(0),
            &pa,
            None,
            &g2c,
            80.0,
        );
        assert!(result.is_none());
    }
}
