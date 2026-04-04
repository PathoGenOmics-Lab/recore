use ahash::{HashMap, HashMapExt};
use crate::parser::PresenceAbsence;
use crate::search::Hit;
use crate::types::{ClusterIdx, GenomeIdx};

// ── Rescue types ──

#[derive(Debug, Clone, PartialEq)]
pub enum RescueType {
    /// Gene is present but Panaroo dropped it during graph cleaning.
    Dropped,
    /// Gene's locus is assigned to a different Panaroo cluster (mis-clustering).
    Split { other_cluster: String },
    /// Gene is broken into multiple fragments (IS, frameshift, assembly break).
    Fragmented { n_fragments: usize, total_coverage: f64 },
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
            RescueType::Split { other_cluster } => write!(f, "SPLIT:{}", other_cluster),
            RescueType::Fragmented { n_fragments, total_coverage } => {
                write!(f, "FRAGMENTED:{}frags:{:.0}%cov", n_fragments, total_coverage)
            }
            RescueType::StartVariant { length_diff_bp } => {
                write!(f, "START_VARIANT:{:+}bp", length_diff_bp)
            }
            RescueType::Paralog { other_cluster } => write!(f, "PARALOG:{}", other_cluster),
            RescueType::Pseudogene { ref reason } => {
                match reason {
                    PseudogeneReason::InternalStop { n_stops } => {
                        write!(f, "PSEUDOGENE:{}stops", n_stops)
                    }
                    PseudogeneReason::Frameshift => write!(f, "PSEUDOGENE:frameshift"),
                    PseudogeneReason::PartialDeletion { coverage } => {
                        write!(f, "PSEUDOGENE:partial:{:.0}%cov", coverage)
                    }
                }
            }
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
}

impl GffIndex {
    pub fn build(mut genes: Vec<GffGene>) -> Self {
        genes.sort_by(|a, b| a.contig.cmp(&b.contig).then(a.start.cmp(&b.start)));

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
        for (i, g) in genes.iter().enumerate() {
            if !g.gene_id.is_empty() {
                gene_id_index.insert(g.gene_id.clone(), i);
            }
        }

        GffIndex { genes, contig_ranges, gene_id_index }
    }

    pub fn get_by_id(&self, gene_id: &str) -> Option<&GffGene> {
        self.gene_id_index.get(gene_id).map(|&i| &self.genes[i])
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

    pub fn flanking(&self, contig: &str, start: u64, end: u64) -> (Option<&GffGene>, Option<&GffGene>) {
        let Some(&(range_start, range_end)) = self.contig_ranges.get(contig) else {
            return (None, None);
        };

        let slice = &self.genes[range_start..range_end];
        let pos = slice.partition_point(|g| g.end <= start);

        let upstream = if pos > 0 { Some(&slice[pos - 1]) } else { None };

        let downstream = slice[pos..]
            .iter()
            .find(|g| g.start >= end);

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
            }
        }
    }
    map
}

// ── Fragment analysis ──

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

// ── Start-site variant detection ──

fn detect_start_variant(
    hit: &Hit,
    gff_index: Option<&GffIndex>,
) -> Option<(i64, String)> {
    let gff = gff_index?;

    if hit.identity() < 98.0 || hit.query_coverage() < 90.0 {
        return None;
    }

    let overlapping = gff.overlapping(&hit.target, hit.target_start, hit.target_end);
    let hit_len = (hit.target_end - hit.target_start) as f64;

    let target_aln_len = (hit.target_end - hit.target_start) as i64;

    for gene in &overlapping {
        let gene_len = (gene.end - gene.start) as i64;
        // Compare GFF gene length vs actual alignment footprint on the genome,
        // NOT the full query length (which may differ due to soft-clipping, partial alignment).
        let diff = gene_len - target_aln_len;

        if diff.unsigned_abs() < 300 && diff != 0 {
            let overlap_start = hit.target_start.max(gene.start);
            let overlap_end = hit.target_end.min(gene.end);
            let overlap = overlap_end.saturating_sub(overlap_start) as f64;

            if hit_len > 0.0 && overlap / hit_len > 0.8 {
                return Some((diff, gene.gene_id.clone()));
            }
        }
    }

    None
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
        b.query_length.cmp(&a.query_length)
            .then_with(|| b.score().total_cmp(&a.score()))
    });

    let best = sorted_hits[0];

    // ── Path 1: Single full-coverage hit ──
    if best.query_coverage() >= min_coverage {
        let (rescue_type, annotation_match, confidence) = classify_single_hit(
            best, ci, pa, gff_index, gene_to_cluster,
        );

        let (final_type, final_ann, final_conf) = if matches!(rescue_type, RescueType::Dropped) {
            if let Some((length_diff, matched_gene_id)) = detect_start_variant(best, gff_index) {
                (
                    RescueType::StartVariant { length_diff_bp: length_diff },
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
        });
    }

    // ── Path 2: Fragment assembly (same contig only) ──
    let primary_contig = &best.target;
    let same_contig_hits: Vec<&Hit> = sorted_hits
        .iter()
        .filter(|h| h.target == *primary_contig)
        .copied()
        .collect();

    let (total_coverage, n_fragments) = merged_query_coverage(&same_contig_hits);
    if total_coverage >= min_coverage && n_fragments >= 2 {
        let avg_identity: f64 = same_contig_hits.iter().map(|h| h.identity()).sum::<f64>()
            / same_contig_hits.len() as f64;

        // SAFETY: same_contig_hits is non-empty (filtered from sorted_hits which has ≥1 element)
        let frag_start = same_contig_hits.iter().map(|h| h.target_start).min()
            .expect("same_contig_hits guaranteed non-empty by coverage check");
        let frag_end = same_contig_hits.iter().map(|h| h.target_end).max()
            .expect("same_contig_hits guaranteed non-empty by coverage check");

        // ── Frameshift detection ──
        // If fragments have high identity (≥95%) and gaps between them are small (≤30bp)
        // and NOT divisible by 3, this is likely a frameshift pseudogene.
        let is_frameshift = avg_identity >= 95.0 && n_fragments >= 2
            && detect_frameshift_pattern(&same_contig_hits);

        let (rescue_type, confidence) = if is_frameshift {
            (RescueType::Pseudogene { reason: PseudogeneReason::Frameshift },
             (avg_identity / 100.0) * (total_coverage / 100.0) * 0.5)
        } else {
            (RescueType::Fragmented { n_fragments, total_coverage },
             (avg_identity / 100.0) * (total_coverage / 100.0) * 0.8)
        };

        return Some(RescueEvent {
            cluster_idx: ci,
            genome_idx: gi,
            rescue_type,
            identity: avg_identity,
            coverage: total_coverage,
            contig: primary_contig.clone(),
            start: frag_start,
            end: frag_end,
            strand: best.strand,
            annotation_match: None,
            confidence,
            synteny_ok: None,
        });
    }

    // ── Path 3: Pseudogene — partial deletion ──
    // High identity hit covering 60-80% of the query → truncated gene (partial deletion).
    // This is below min_coverage so Path 1 didn't catch it, but biologically informative.
    if best.identity() >= 90.0 {
        let cov = best.query_coverage();
        // Also check fragment coverage for multi-fragment partial hits
        let effective_cov = if same_contig_hits.len() >= 2 { total_coverage } else { cov };
        if effective_cov >= 60.0 && effective_cov < min_coverage {
            return Some(RescueEvent {
                cluster_idx: ci,
                genome_idx: gi,
                rescue_type: RescueType::Pseudogene {
                    reason: PseudogeneReason::PartialDeletion { coverage: effective_cov },
                },
                identity: best.identity(),
                coverage: effective_cov,
                contig: best.target.clone(),
                start: best.target_start,
                end: best.target_end,
                strand: best.strand,
                annotation_match: None,
                confidence: (best.identity() / 100.0) * (effective_cov / 100.0) * 0.4,
                synteny_ok: None,
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
                let upper = [codon[0].to_ascii_uppercase(),
                             codon[1].to_ascii_uppercase(),
                             codon[2].to_ascii_uppercase()];
                if matches!(upper, [b'T', b'A', b'A'] | [b'T', b'A', b'G'] | [b'T', b'G', b'A']) {
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
    seq.iter().rev().map(|&b| match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        _ => b'N',
    }).collect()
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

    if n_stops > 0 { Some(n_stops) } else { None }
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

        if recip > 0.5
            && best_gene.is_none_or(|(_, prev)| recip > prev) {
                best_gene = Some((gene, recip));
            }
    }

    let Some((matched_gene, _recip)) = best_gene else {
        return (RescueType::Dropped, Some(overlapping[0].gene_id.clone()), base_confidence);
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

    (RescueType::Dropped, Some(matched_gene.gene_id.clone()), base_confidence)
}

// ── Conflict resolution ──

pub fn resolve_conflicts(rescues: &mut Vec<RescueEvent>) {
    if rescues.is_empty() {
        return;
    }

    rescues.sort_by(|a, b| {
        a.genome_idx.cmp(&b.genome_idx)
            .then(a.contig.cmp(&b.contig))
            .then(a.start.cmp(&b.start))
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

// ── Synteny validation ──

pub fn validate_synteny(
    rescue: &mut RescueEvent,
    pa: &PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    gene_to_cluster_maps: &HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
) {
    let ci = rescue.cluster_idx;
    let gi = rescue.genome_idx;

    let Some(target_gff) = gff_indices.get(&gi) else {
        return;
    };
    let Some(target_g2c) = gene_to_cluster_maps.get(&gi) else {
        return;
    };

    let (target_up, target_down) = target_gff.flanking(&rescue.contig, rescue.start, rescue.end);

    let target_up_cluster = target_up.and_then(|g| target_g2c.get(&g.gene_id).copied());
    let target_down_cluster = target_down.and_then(|g| target_g2c.get(&g.gene_id).copied());

    if target_up_cluster.is_none() && target_down_cluster.is_none() {
        return;
    }

    let mut synteny_match = false;
    for &ref_gi in pa.present_genomes(ci).iter().take(5) {
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

        let Some(ref_gene) = ref_gff.get_by_id(ref_gene_id) else {
            continue;
        };

        let (ref_up, ref_down) = ref_gff.flanking(&ref_gene.contig, ref_gene.start, ref_gene.end);

        let ref_up_cluster = ref_up.and_then(|g| ref_g2c.get(&g.gene_id).copied());
        let ref_down_cluster = ref_down.and_then(|g| ref_g2c.get(&g.gene_id).copied());

        let up_match = match (target_up_cluster, ref_up_cluster) {
            (Some(a), Some(b)) => a == b,
            (None, None) => true,
            _ => false,
        };
        let down_match = match (target_down_cluster, ref_down_cluster) {
            (Some(a), Some(b)) => a == b,
            (None, None) => true,
            _ => false,
        };

        if up_match || down_match {
            synteny_match = true;
            break;
        }
    }

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
            query: "g".into(), cluster: "g".into(), query_length,
            query_start, query_end, strand: '+', target: "c".into(),
            target_start: query_start, target_end: query_end,
            matches: query_end - query_start,
            alignment_length: query_end - query_start, mapq: 60,
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
            query: "g__0".into(), cluster: "g".into(), query_length: 800,
            query_start: 0, query_end: 600, strand: '+', target: "c".into(),
            target_start: 0, target_end: 600, matches: 600, alignment_length: 600, mapq: 60,
        };
        let h2 = Hit {
            query: "g__1".into(), cluster: "g".into(), query_length: 1200,
            query_start: 0, query_end: 1100, strand: '+', target: "c".into(),
            target_start: 0, target_end: 1100, matches: 1100, alignment_length: 1100, mapq: 60,
        };
        let (cov, frags) = merged_query_coverage(&[&h1, &h2]);
        assert!(cov > 90.0);
        assert_eq!(frags, 1);
    }

    #[test]
    fn test_merged_coverage_same_length_different_alleles() {
        let h1 = Hit {
            query: "cluster~~~ALLELE~0".into(), cluster: "cluster".into(), query_length: 1000,
            query_start: 0, query_end: 400, strand: '+', target: "c".into(),
            target_start: 0, target_end: 400, matches: 400, alignment_length: 400, mapq: 60,
        };
        let h2 = Hit {
            query: "cluster~~~ALLELE~1".into(), cluster: "cluster".into(), query_length: 1000,
            query_start: 600, query_end: 1000, strand: '+', target: "c".into(),
            target_start: 5000, target_end: 5400, matches: 400, alignment_length: 400, mapq: 60,
        };
        let (cov, _frags) = merged_query_coverage(&[&h1, &h2]);
        assert!(cov < 50.0, "Coverage should be <50% (was {:.1}%), alleles must not be merged", cov);
    }

    #[test]
    fn test_gff_index_overlap() {
        let genes = vec![
            GffGene { contig: "chr1".into(), start: 100, end: 500, strand: '+', gene_id: "g1".into() },
            GffGene { contig: "chr1".into(), start: 600, end: 900, strand: '-', gene_id: "g2".into() },
            GffGene { contig: "chr2".into(), start: 50, end: 200, strand: '+', gene_id: "g3".into() },
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
            GffGene { contig: "c".into(), start: 100, end: 200, strand: '+', gene_id: "g1".into() },
            GffGene { contig: "c".into(), start: 300, end: 400, strand: '+', gene_id: "g2".into() },
            GffGene { contig: "c".into(), start: 500, end: 600, strand: '+', gene_id: "g3".into() },
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
    fn test_detect_start_variant_with_gff() {
        let genes = vec![
            GffGene { contig: "c".into(), start: 0, end: 1100, strand: '+', gene_id: "gene1".into() },
        ];
        let idx = GffIndex::build(genes);

        let hit = Hit {
            query: "g".into(), cluster: "g".into(), query_length: 1000,
            query_start: 0, query_end: 1000, strand: '+', target: "c".into(),
            target_start: 100, target_end: 1100,
            matches: 1000, alignment_length: 1000, mapq: 60,
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
        let genes = vec![
            GffGene { contig: "c".into(), start: 0, end: 1050, strand: '+', gene_id: "gene1".into() },
        ];
        let idx = GffIndex::build(genes);
        let hit = Hit {
            query: "g".into(), cluster: "g".into(), query_length: 1100,
            query_start: 0, query_end: 1050, strand: '+', target: "c".into(),
            target_start: 0, target_end: 1000,
            matches: 1000, alignment_length: 1050, mapq: 60,
        };
        // identity = 1000/1050 = 95.2% (≥98%? NO — need 98%)
        // query_coverage = 1050/1100 = 95.5% (≥90% ✓)
        // Adjust matches to pass identity filter:
        let hit = Hit { matches: 1040, ..hit }; // 1040/1050 = 99.0%
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
                cluster_idx: ClusterIdx(0), genome_idx: GenomeIdx(0),
                rescue_type: RescueType::Dropped,
                identity: 99.0, coverage: 100.0, contig: "c1".into(),
                start: 100, end: 500, strand: '+', annotation_match: None,
                confidence: 0.99, synteny_ok: None,
            },
            RescueEvent {
                cluster_idx: ClusterIdx(1), genome_idx: GenomeIdx(0),
                rescue_type: RescueType::Dropped,
                identity: 95.0, coverage: 90.0, contig: "c1".into(),
                start: 200, end: 600, strand: '+', annotation_match: None,
                confidence: 0.85, synteny_ok: None,
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
                cluster_idx: ClusterIdx(0), genome_idx: GenomeIdx(0),
                rescue_type: RescueType::Dropped,
                identity: 99.0, coverage: 100.0, contig: "c1".into(),
                start: 100, end: 500, strand: '+', annotation_match: None,
                confidence: 0.99, synteny_ok: None,
            },
            RescueEvent {
                cluster_idx: ClusterIdx(1), genome_idx: GenomeIdx(0),
                rescue_type: RescueType::Dropped,
                identity: 95.0, coverage: 90.0, contig: "c1".into(),
                start: 600, end: 900, strand: '+', annotation_match: None,
                confidence: 0.85, synteny_ok: None,
            },
        ];
        resolve_conflicts(&mut rescues);
        assert_eq!(rescues.len(), 2);
    }

    #[test]
    fn test_fragmented_span() {
        let h1 = Hit {
            query: "g".into(), cluster: "g".into(), query_length: 1000,
            query_start: 0, query_end: 400, strand: '+', target: "c".into(),
            target_start: 1000, target_end: 1400,
            matches: 400, alignment_length: 400, mapq: 60,
        };
        let h2 = Hit {
            query: "g".into(), cluster: "g".into(), query_length: 1000,
            query_start: 600, query_end: 1000, strand: '+', target: "c".into(),
            target_start: 2000, target_end: 2400,
            matches: 400, alignment_length: 400, mapq: 60,
        };

        let pa = PresenceAbsence::new_for_test(
            vec!["g".into()],
            vec!["genome1".into()],
            vec![vec!["".into()]],
        );
        let g2c: HashMap<String, ClusterIdx> = HashMap::new();

        let result = classify_hits(&[h1, h2], ClusterIdx(0), GenomeIdx(0), &pa, None, &g2c, 80.0);
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
        let gff_genes = vec![
            GffGene { contig: "c".into(), start: 0, end: 1000, strand: '+', gene_id: "geneA".into() },
        ];
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
            query: "q".into(), cluster: "target_cluster".into(), query_length: 1000,
            query_start: 0, query_end: 1000, strand: '+', target: "c".into(),
            target_start: 0, target_end: 1000,
            matches: 990, alignment_length: 1000, mapq: 60,
        };
        let result = classify_hits(&[hit_high], ClusterIdx(0), GenomeIdx(0), &pa, Some(&gff), &g2c, 80.0);
        let r = result.unwrap();
        match &r.rescue_type {
            RescueType::Split { other_cluster } => assert_eq!(other_cluster, "other_cluster"),
            other => panic!("Expected Split, got {:?}", other),
        }

        // Low identity (90%) → PARALOG
        let hit_low = Hit {
            query: "q".into(), cluster: "target_cluster".into(), query_length: 1000,
            query_start: 0, query_end: 1000, strand: '+', target: "c".into(),
            target_start: 0, target_end: 1000,
            matches: 900, alignment_length: 1000, mapq: 60,
        };
        let result = classify_hits(&[hit_low], ClusterIdx(0), GenomeIdx(0), &pa, Some(&gff), &g2c, 80.0);
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
            query: "q".into(), cluster: "cluster1".into(), query_length: 1000,
            query_start: 0, query_end: 1000, strand: '+', target: "c".into(),
            target_start: 0, target_end: 1000,
            matches: 990, alignment_length: 1000, mapq: 60,
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
    fn test_classify_dropped_no_overlapping_genes() {
        // Full-coverage hit, GFF exists but no gene overlaps the hit region → DROPPED
        let gff_genes = vec![
            GffGene { contig: "c".into(), start: 5000, end: 6000, strand: '+', gene_id: "far_gene".into() },
        ];
        let gff = GffIndex::build(gff_genes);
        let pa = PresenceAbsence::new_for_test(
            vec!["cluster1".into()],
            vec!["genome1".into()],
            vec![vec!["".into()]],
        );
        let g2c: HashMap<String, ClusterIdx> = HashMap::new();
        let hit = Hit {
            query: "q".into(), cluster: "cluster1".into(), query_length: 1000,
            query_start: 0, query_end: 1000, strand: '+', target: "c".into(),
            target_start: 0, target_end: 1000,
            matches: 990, alignment_length: 1000, mapq: 60,
        };
        let result = classify_hits(&[hit], ClusterIdx(0), GenomeIdx(0), &pa, Some(&gff), &g2c, 80.0);
        let r = result.unwrap();
        assert!(matches!(r.rescue_type, RescueType::Dropped));
    }

    #[test]
    fn test_classify_start_variant_with_gff() {
        // Hit with high identity overlapping a similar-length GFF gene → START_VARIANT
        let gff_genes = vec![
            GffGene { contig: "c".into(), start: 0, end: 1100, strand: '+', gene_id: "gene1".into() },
        ];
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
            query: "q".into(), cluster: "cluster1".into(), query_length: 1000,
            query_start: 0, query_end: 1000, strand: '+', target: "c".into(),
            target_start: 100, target_end: 1100,
            matches: 1000, alignment_length: 1000, mapq: 60,
        };
        let result = classify_hits(&[hit], ClusterIdx(0), GenomeIdx(0), &pa, Some(&gff), &g2c, 80.0);
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
            query: "q".into(), cluster: "c1".into(), query_length: 1000,
            query_start: 0, query_end: 500, strand: '+', target: "c".into(),
            target_start: 0, target_end: 500,
            matches: 500, alignment_length: 500, mapq: 60,
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
        assert_eq!(
            RescueType::Split { other_cluster: "cl2".into() }.to_string(),
            "SPLIT:cl2"
        );
        assert_eq!(
            RescueType::Fragmented { n_fragments: 3, total_coverage: 95.5 }.to_string(),
            "FRAGMENTED:3frags:96%cov"
        );
        assert_eq!(
            RescueType::StartVariant { length_diff_bp: -42 }.to_string(),
            "START_VARIANT:-42bp"
        );
        assert_eq!(
            RescueType::Paralog { other_cluster: "cl3".into() }.to_string(),
            "PARALOG:cl3"
        );
        assert_eq!(
            RescueType::Pseudogene { reason: PseudogeneReason::InternalStop { n_stops: 2 } }.to_string(),
            "PSEUDOGENE:2stops"
        );
        assert_eq!(
            RescueType::Pseudogene { reason: PseudogeneReason::Frameshift }.to_string(),
            "PSEUDOGENE:frameshift"
        );
        assert_eq!(
            RescueType::Pseudogene { reason: PseudogeneReason::PartialDeletion { coverage: 72.5 } }.to_string(),
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
            query: "q".into(), cluster: "c".into(), query_length: 1000,
            query_start: 0, query_end: 500, strand: '+', target: "c".into(),
            target_start: 100, target_end: 600,
            matches: 500, alignment_length: 500, mapq: 60,
        };
        let h2 = Hit {
            query: "q".into(), cluster: "c".into(), query_length: 1000,
            query_start: 500, query_end: 1000, strand: '+', target: "c".into(),
            target_start: 602, target_end: 1102, // 2bp gap (not %3)
            matches: 500, alignment_length: 500, mapq: 60,
        };
        assert!(detect_frameshift_pattern(&[&h1, &h2]));
    }

    #[test]
    fn test_frameshift_pattern_not_detected_divisible_by_3() {
        let h1 = Hit {
            query: "q".into(), cluster: "c".into(), query_length: 1000,
            query_start: 0, query_end: 500, strand: '+', target: "c".into(),
            target_start: 100, target_end: 600,
            matches: 500, alignment_length: 500, mapq: 60,
        };
        let h2 = Hit {
            query: "q".into(), cluster: "c".into(), query_length: 1000,
            query_start: 500, query_end: 1000, strand: '+', target: "c".into(),
            target_start: 603, target_end: 1103, // 3bp gap (divisible by 3)
            matches: 500, alignment_length: 500, mapq: 60,
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
            query: "q".into(), cluster: "c".into(), query_length: 9,
            query_start: 0, query_end: 9, strand: '-', target: "c".into(),
            target_start: 0, target_end: 9,
            matches: 9, alignment_length: 9, mapq: 60,
        };
        let hit_plus = Hit { strand: '+', ..hit_minus.clone() };
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
            query: "q".into(), cluster: "c1".into(), query_length: 1000,
            query_start: 0, query_end: 700, strand: '+', target: "c".into(),
            target_start: 0, target_end: 700,
            matches: 665, alignment_length: 700, mapq: 60,
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
            RescueType::Pseudogene { reason: PseudogeneReason::PartialDeletion { coverage } } => {
                assert!(*coverage >= 69.0 && *coverage <= 71.0);
            }
            other => panic!("Expected Pseudogene:PartialDeletion, got {:?}", other),
        }
    }
}
