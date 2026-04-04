use ahash::{HashMap, HashMapExt};
use crate::parser::PresenceAbsence;
use crate::search::Hit;

// ── Rescue types ──

#[derive(Debug, Clone, PartialEq)]
pub enum RescueType {
    /// Gene is present but Panaroo dropped it during graph cleaning.
    Dropped,
    /// Gene's locus is assigned to a different Panaroo cluster (mis-clustering).
    Split { other_cluster: String },
    /// Gene is broken into multiple fragments (IS, frameshift, assembly break).
    Fragmented { n_fragments: usize, total_coverage: f64 },
    /// Start site variant: gene is there but with a different start codon.
    StartVariant { extension_bp: i64 },
    /// Gene hit overlaps a paralog from a different cluster.
    #[allow(dead_code)]
    Paralog { other_cluster: String },
}

impl std::fmt::Display for RescueType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RescueType::Dropped => write!(f, "DROPPED"),
            RescueType::Split { other_cluster } => write!(f, "SPLIT:{}", other_cluster),
            RescueType::Fragmented { n_fragments, total_coverage } => {
                write!(f, "FRAGMENTED:{}frags:{:.0}%cov", n_fragments, total_coverage)
            }
            RescueType::StartVariant { extension_bp } => {
                write!(f, "START_VARIANT:{:+}bp", extension_bp)
            }
            RescueType::Paralog { other_cluster } => write!(f, "PARALOG:{}", other_cluster),
        }
    }
}

// ── Rescue event ──

#[derive(Debug, Clone)]
pub struct RescueEvent {
    pub cluster_idx: usize,
    pub genome_idx: usize,
    pub rescue_type: RescueType,
    pub identity: f64,
    pub coverage: f64,
    pub contig: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub annotation_match: Option<String>,
    /// Confidence score: 0.0 to 1.0. Used for conflict resolution.
    pub confidence: f64,
    /// Synteny validated (if --synteny enabled and GFF provided).
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

/// Index for fast coordinate overlap queries.
pub struct GffIndex {
    genes: Vec<GffGene>,
    contig_ranges: HashMap<String, (usize, usize)>,
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

        GffIndex { genes, contig_ranges }
    }

    /// Find all genes overlapping [query_start, query_end) on a contig.
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

    /// Get the gene immediately before and after a position on a contig.
    /// Returns (upstream_gene_id, downstream_gene_id).
    pub fn flanking(&self, contig: &str, start: u64, end: u64) -> (Option<String>, Option<String>) {
        let Some(&(range_start, range_end)) = self.contig_ranges.get(contig) else {
            return (None, None);
        };

        let slice = &self.genes[range_start..range_end];

        // Find insertion point
        let pos = slice.partition_point(|g| g.end <= start);

        let upstream = if pos > 0 {
            Some(slice[pos - 1].gene_id.clone())
        } else {
            None
        };

        // Find first gene starting after our region
        let downstream = slice[pos..]
            .iter()
            .find(|g| g.start >= end)
            .map(|g| g.gene_id.clone());

        (upstream, downstream)
    }
}

// ── Gene-to-cluster index ──

/// Build a lookup from gene_id → cluster_idx for a specific genome.
/// Handles Panaroo's semicolon-separated multi-gene cells.
pub fn build_gene_to_cluster(pa: &PresenceAbsence, genome_idx: usize) -> HashMap<String, usize> {
    let mut map = HashMap::new();
    for (ci, present) in pa.present.iter().enumerate() {
        if !present.contains(&genome_idx) {
            continue;
        }
        let cell = &pa.cell_values[ci][genome_idx];
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

/// Compute combined query coverage from multiple hits.
/// Merges overlapping intervals on the query to avoid double-counting.
fn merged_query_coverage(hits: &[&Hit]) -> (f64, usize) {
    if hits.is_empty() {
        return (0.0, 0);
    }

    let query_length = hits[0].query_length;
    if query_length == 0 {
        return (0.0, 0);
    }

    let mut intervals: Vec<(u64, u64)> = hits
        .iter()
        .map(|h| (h.query_start, h.query_end))
        .collect();
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
    let coverage = (covered as f64 / query_length as f64) * 100.0;
    (coverage, merged.len())
}

// ── Start-site variant detection ──

/// Detect if a hit represents a start-site variant.
/// The gene is there, but with a different start → Panaroo puts it in another cluster.
/// Signature: ~100% identity, coverage close to 100% but alignment extends beyond query.
fn detect_start_variant(hit: &Hit) -> Option<i64> {
    if hit.identity() < 98.0 {
        return None;
    }

    let coverage = hit.query_coverage();
    if coverage < 95.0 {
        return None;
    }

    // How much does the target alignment extend beyond the query length?
    let target_span = hit.target_end - hit.target_start;
    let query_span = hit.query_end - hit.query_start;

    if target_span <= query_span {
        return None;
    }

    let extension = target_span as i64 - query_span as i64;
    // Only flag as start variant if extension is small (< 300bp, typical start variation)
    if extension > 0 && extension < 300 {
        Some(extension)
    } else {
        None
    }
}

// ── Main classification ──

/// Classify all hits for a single (cluster, genome) pair.
///
/// Decision tree:
/// 1. Best single hit covers ≥min_coverage → classify as DROPPED/SPLIT/START_VARIANT
/// 2. Multiple partial hits cover ≥min_coverage together → FRAGMENTED
/// 3. Otherwise → no rescue (true absence or below thresholds)
pub fn classify_hits(
    hits: &[Hit],
    cluster_idx: usize,
    genome_idx: usize,
    pa: &PresenceAbsence,
    gff_index: Option<&GffIndex>,
    gene_to_cluster: &HashMap<String, usize>,
    min_coverage: f64,
) -> Option<RescueEvent> {
    if hits.is_empty() {
        return None;
    }

    // Sort hits by score (best first)
    let mut sorted_hits: Vec<&Hit> = hits.iter().collect();
    sorted_hits.sort_by(|a, b| {
        b.score().partial_cmp(&a.score()).unwrap()
    });

    let best = sorted_hits[0];

    // ── Path 1: Single full-coverage hit ──
    if best.query_coverage() >= min_coverage {
        // Check for start-site variant first
        if let Some(extension) = detect_start_variant(best) {
            let confidence = best.identity() / 100.0 * 0.95; // slightly lower than DROPPED
            return Some(RescueEvent {
                cluster_idx,
                genome_idx,
                rescue_type: RescueType::StartVariant { extension_bp: extension },
                identity: best.identity(),
                coverage: best.query_coverage(),
                contig: best.target.clone(),
                start: best.target_start,
                end: best.target_end,
                strand: best.strand,
                annotation_match: None,
                confidence,
                synteny_ok: None,
            });
        }

        let (rescue_type, annotation_match, confidence) = classify_single_hit(
            best,
            cluster_idx,
            pa,
            gff_index,
            gene_to_cluster,
        );

        return Some(RescueEvent {
            cluster_idx,
            genome_idx,
            rescue_type,
            identity: best.identity(),
            coverage: best.query_coverage(),
            contig: best.target.clone(),
            start: best.target_start,
            end: best.target_end,
            strand: best.strand,
            annotation_match,
            confidence,
            synteny_ok: None,
        });
    }

    // ── Path 2: Fragment assembly ──
    let (total_coverage, n_fragments) = merged_query_coverage(&sorted_hits);
    if total_coverage >= min_coverage && n_fragments >= 2 {
        // Confidence is lower for fragmented genes (harder to be sure)
        let avg_identity: f64 = sorted_hits.iter().map(|h| h.identity()).sum::<f64>()
            / sorted_hits.len() as f64;
        let confidence = (avg_identity / 100.0) * (total_coverage / 100.0) * 0.8;

        return Some(RescueEvent {
            cluster_idx,
            genome_idx,
            rescue_type: RescueType::Fragmented {
                n_fragments,
                total_coverage,
            },
            identity: avg_identity,
            coverage: total_coverage,
            contig: best.target.clone(),
            start: best.target_start,
            end: best.target_end,
            strand: best.strand,
            annotation_match: None,
            confidence,
            synteny_ok: None,
        });
    }

    None
}

/// Classify a single full-coverage hit.
fn classify_single_hit(
    hit: &Hit,
    cluster_idx: usize,
    pa: &PresenceAbsence,
    gff_index: Option<&GffIndex>,
    gene_to_cluster: &HashMap<String, usize>,
) -> (RescueType, Option<String>, f64) {
    let base_confidence = hit.identity() / 100.0 * hit.query_coverage() / 100.0;

    let Some(gff) = gff_index else {
        return (RescueType::Dropped, None, base_confidence);
    };

    let overlapping = gff.overlapping(&hit.target, hit.target_start, hit.target_end);

    if overlapping.is_empty() {
        return (RescueType::Dropped, None, base_confidence);
    }

    // Find best reciprocal overlap
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
            && (best_gene.is_none() || recip > best_gene.unwrap().1) {
                best_gene = Some((gene, recip));
            }
    }

    let Some((matched_gene, _recip)) = best_gene else {
        return (RescueType::Dropped, Some(overlapping[0].gene_id.clone()), base_confidence);
    };

    // Is this gene assigned to another cluster? → SPLIT
    if let Some(&other_ci) = gene_to_cluster.get(&matched_gene.gene_id) {
        if other_ci != cluster_idx {
            return (
                RescueType::Split {
                    other_cluster: pa.clusters[other_ci].clone(),
                },
                Some(matched_gene.gene_id.clone()),
                base_confidence * 1.1, // higher confidence: we know exactly where the gene went
            );
        }
    }

    (RescueType::Dropped, Some(matched_gene.gene_id.clone()), base_confidence)
}

// ── Conflict resolution ──

/// Resolve conflicts where multiple clusters claim the same genomic region.
/// The cluster with the highest confidence wins; others are removed.
pub fn resolve_conflicts(rescues: &mut Vec<RescueEvent>) {
    if rescues.is_empty() {
        return;
    }

    // Group rescues by (genome_idx, contig, region).
    // Two rescues conflict if they share the same genome and their target regions overlap.
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

            // Different genome → no conflict
            if ri.genome_idx != rj.genome_idx {
                break; // sorted by genome, so no more conflicts with i
            }

            // Different contig → no conflict
            if ri.contig != rj.contig {
                continue;
            }

            // Past the region → no more overlaps for i
            if rj.start >= ri.end {
                break;
            }

            // Overlapping regions → conflict!
            // The one with higher confidence wins.
            if ri.confidence >= rj.confidence {
                to_remove[j] = true;
            } else {
                to_remove[i] = true;
                break; // i is dead, move on
            }
        }
    }

    // Remove losers
    let mut idx = 0;
    rescues.retain(|_| {
        let keep = !to_remove[idx];
        idx += 1;
        keep
    });
}

// ── Synteny validation ──

/// Validate a rescue by checking flanking gene conservation.
/// For a rescue in genome G at position P, check if the genes flanking P in G
/// are the same genes flanking the cluster in a genome where it IS present.
pub fn validate_synteny(
    rescue: &mut RescueEvent,
    pa: &PresenceAbsence,
    gff_indices: &HashMap<usize, GffIndex>,
) {
    let ci = rescue.cluster_idx;
    let gi = rescue.genome_idx;

    // Get flanking genes of the rescue site in the target genome
    let Some(target_gff) = gff_indices.get(&gi) else {
        return;
    };
    let (target_up, target_down) = target_gff.flanking(&rescue.contig, rescue.start, rescue.end);

    if target_up.is_none() && target_down.is_none() {
        return; // can't validate without flanking info
    }

    // Find a reference genome where this cluster IS present and has GFF
    let mut synteny_match = false;
    for ref_gi in pa.present[ci].iter().take(5) {
        // check up to 5 reference genomes
        let Some(ref_gff) = gff_indices.get(ref_gi) else {
            continue;
        };

        // Get the gene's coordinates in the reference genome
        let ref_cell = &pa.cell_values[ci][*ref_gi];
        if ref_cell.is_empty() {
            continue;
        }

        // Find the gene's position in the reference GFF
        let ref_gene_id = ref_cell.split(';').next().unwrap_or(ref_cell).trim();

        // Look up coordinates in the reference GFF index
        // We need to scan for the gene_id (could optimize with another index)
        let mut ref_start = 0u64;
        let mut ref_end = 0u64;
        let mut ref_contig = String::new();
        let mut found = false;

        for gene in &ref_gff.genes {
            if gene.gene_id == ref_gene_id {
                ref_start = gene.start;
                ref_end = gene.end;
                ref_contig = gene.contig.clone();
                found = true;
                break;
            }
        }

        if !found {
            continue;
        }

        let (ref_up, ref_down) = ref_gff.flanking(&ref_contig, ref_start, ref_end);

        // Compare flanking genes by gene_id
        let up_match = match (&target_up, &ref_up) {
            (Some(a), Some(b)) => a == b,
            (None, None) => true,
            _ => false,
        };
        let down_match = match (&target_down, &ref_down) {
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

    // Boost or penalize confidence based on synteny
    if synteny_match {
        rescue.confidence *= 1.2;
    } else {
        rescue.confidence *= 0.6;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_hit(query_start: u64, query_end: u64, query_length: u64) -> Hit {
        Hit {
            query: "g".into(),
            cluster: "g".into(),
            query_length,
            query_start,
            query_end,
            strand: '+',
            target: "c".into(),
            target_start: query_start, // simplification
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
        assert_eq!(up.as_deref(), Some("g1"));
        assert_eq!(down.as_deref(), Some("g3"));

        let (up, down) = idx.flanking("c", 100, 200);
        assert_eq!(up, None);
        assert_eq!(down.as_deref(), Some("g2"));
    }

    #[test]
    fn test_detect_start_variant() {
        // Normal hit: no start variant
        let normal = make_hit(0, 1000, 1000);
        assert!(detect_start_variant(&normal).is_none());

        // Start variant: target spans more than query
        let sv = Hit {
            query: "g".into(),
            cluster: "g".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 1000,
            strand: '+',
            target: "c".into(),
            target_start: 0,
            target_end: 1100, // 100bp extension
            matches: 1000,
            alignment_length: 1000,
            mapq: 60,
        };
        let ext = detect_start_variant(&sv);
        assert_eq!(ext, Some(100));
    }

    #[test]
    fn test_conflict_resolution() {
        let mut rescues = vec![
            RescueEvent {
                cluster_idx: 0, genome_idx: 0, rescue_type: RescueType::Dropped,
                identity: 99.0, coverage: 100.0, contig: "c1".into(),
                start: 100, end: 500, strand: '+', annotation_match: None,
                confidence: 0.99, synteny_ok: None,
            },
            RescueEvent {
                cluster_idx: 1, genome_idx: 0, rescue_type: RescueType::Dropped,
                identity: 95.0, coverage: 90.0, contig: "c1".into(),
                start: 200, end: 600, strand: '+', annotation_match: None,
                confidence: 0.85, synteny_ok: None,
            },
        ];
        resolve_conflicts(&mut rescues);
        assert_eq!(rescues.len(), 1);
        assert_eq!(rescues[0].cluster_idx, 0); // higher confidence wins
    }
}
