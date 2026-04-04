use ahash::{HashMap, HashMapExt};
use crate::parser::PresenceAbsence;
use crate::search::Hit;

/// Classification of a rescued gene.
#[derive(Debug, Clone, PartialEq)]
pub enum RescueType {
    /// Gene is present in the genome but Panaroo dropped it during graph cleaning.
    Dropped,
    /// Gene's locus is occupied by a different Panaroo cluster (mis-clustering).
    Split { other_cluster: String },
    /// Gene is fragmented: multiple partial hits cover distinct regions of the query.
    /// Likely caused by IS insertion, frameshift, or assembly break.
    Fragmented { n_fragments: usize, total_coverage: f64 },
    /// Gene hit overlaps with a paralog assigned to a different cluster.
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
            RescueType::Paralog { other_cluster } => write!(f, "PARALOG:{}", other_cluster),
        }
    }
}

/// A single gene rescue event.
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
}

/// GFF annotation entry for coordinate-based classification.
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
/// Genes sorted by contig + start for binary-search-based interval queries.
pub struct GffIndex {
    /// Sorted by (contig, start).
    genes: Vec<GffGene>,
    /// contig_name → (start_idx, end_idx) in `genes`.
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

        // Binary search for first gene that could overlap (gene.end > query_start)
        let scan_start = slice.partition_point(|g| g.end <= query_start);

        let mut results = Vec::new();
        for g in &slice[scan_start..] {
            if g.start >= query_end {
                break; // past the query region
            }
            results.push(g);
        }
        results
    }
}

/// Build a lookup from gene_id → cluster_idx for a specific genome.
/// Handles Panaroo's semicolon-separated multi-gene cells.
pub fn build_gene_to_cluster(pa: &PresenceAbsence, genome_idx: usize) -> HashMap<String, usize> {
    let mut map = HashMap::new();
    for (ci, present) in pa.present.iter().enumerate() {
        if !present.contains(&genome_idx) {
            continue;
        }
        let cell = &pa.cell_values[ci][genome_idx];
        // Panaroo can have "gene1;gene2" for paralogs in one cell
        for gene_id in cell.split(';') {
            let gene_id = gene_id.trim();
            if !gene_id.is_empty() {
                map.insert(gene_id.to_string(), ci);
            }
        }
    }
    map
}

/// Compute combined query coverage from multiple partial hits.
/// Merges overlapping intervals on the query to avoid double-counting.
fn merged_query_coverage(hits: &[&Hit]) -> (f64, usize) {
    if hits.is_empty() {
        return (0.0, 0);
    }

    let query_length = hits[0].query_length;
    if query_length == 0 {
        return (0.0, 0);
    }

    // Collect query intervals and merge
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

/// Classify all hits for a single (cluster, genome) pair.
///
/// Decision tree:
/// 1. Single hit with full coverage → check if locus belongs to another cluster (SPLIT) or not (DROPPED)
/// 2. Multiple partial hits → FRAGMENTED (gene broken by IS/frameshift/assembly)
/// 3. Hit overlaps a paralog from a different cluster → PARALOG
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

    // Sort hits by identity × coverage score (best first)
    let mut sorted_hits: Vec<&Hit> = hits.iter().collect();
    sorted_hits.sort_by(|a, b| {
        let score_a = a.identity() * a.query_coverage();
        let score_b = b.identity() * b.query_coverage();
        score_b.partial_cmp(&score_a).unwrap()
    });

    let best = sorted_hits[0];

    // ── Check single best hit ──
    if best.query_coverage() >= min_coverage {
        // Full-length hit. Is it a SPLIT or a DROPPED?
        let (rescue_type, annotation_match) = classify_single_hit(
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
        });
    }

    // ── Check fragment assembly: multiple partial hits ──
    let (total_coverage, n_fragments) = merged_query_coverage(&sorted_hits);
    if total_coverage >= min_coverage && n_fragments >= 2 {
        return Some(RescueEvent {
            cluster_idx,
            genome_idx,
            rescue_type: RescueType::Fragmented {
                n_fragments,
                total_coverage,
            },
            identity: best.identity(),
            coverage: total_coverage,
            contig: best.target.clone(),
            start: best.target_start,
            end: best.target_end,
            strand: best.strand,
            annotation_match: None,
        });
    }

    None // No rescue: true absence or below thresholds
}

/// Classify a single full-coverage hit.
fn classify_single_hit(
    hit: &Hit,
    cluster_idx: usize,
    pa: &PresenceAbsence,
    gff_index: Option<&GffIndex>,
    gene_to_cluster: &HashMap<String, usize>,
) -> (RescueType, Option<String>) {
    // Without GFF, we can't do coordinate-based classification
    let Some(gff) = gff_index else {
        return (RescueType::Dropped, None);
    };

    // Find annotated genes overlapping the hit region
    let overlapping = gff.overlapping(&hit.target, hit.target_start, hit.target_end);

    if overlapping.is_empty() {
        // Hit in an unannotated region → gene was dropped (or novel)
        return (RescueType::Dropped, None);
    }

    // Find the best-overlapping gene (highest reciprocal overlap)
    let hit_len = (hit.target_end - hit.target_start) as f64;
    let mut best_gene: Option<(&GffGene, f64)> = None;

    for gene in &overlapping {
        let overlap_start = hit.target_start.max(gene.start);
        let overlap_end = hit.target_end.min(gene.end);
        let overlap = overlap_end.saturating_sub(overlap_start) as f64;
        let gene_len = (gene.end - gene.start) as f64;

        // Reciprocal overlap: minimum of (overlap/hit_len, overlap/gene_len)
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

    let Some((matched_gene, _)) = best_gene else {
        // Overlapping genes exist but none with >50% reciprocal overlap
        return (RescueType::Dropped, Some(overlapping[0].gene_id.clone()));
    };

    // Check if this gene is assigned to another cluster → SPLIT
    if let Some(&other_ci) = gene_to_cluster.get(&matched_gene.gene_id) {
        if other_ci != cluster_idx {
            return (
                RescueType::Split {
                    other_cluster: pa.clusters[other_ci].clone(),
                },
                Some(matched_gene.gene_id.clone()),
            );
        }
    }

    // Gene is annotated but not in any Panaroo cluster → DROPPED by graph cleaning
    (RescueType::Dropped, Some(matched_gene.gene_id.clone()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merged_coverage_no_overlap() {
        let h1 = Hit {
            query: "g".into(), query_length: 1000, query_start: 0, query_end: 400,
            strand: '+', target: "c".into(), target_start: 0, target_end: 400,
            matches: 400, alignment_length: 400, mapq: 60,
        };
        let h2 = Hit {
            query: "g".into(), query_length: 1000, query_start: 600, query_end: 1000,
            strand: '+', target: "c".into(), target_start: 600, target_end: 1000,
            matches: 400, alignment_length: 400, mapq: 60,
        };
        let (cov, frags) = merged_query_coverage(&[&h1, &h2]);
        assert!((cov - 80.0).abs() < 0.01);
        assert_eq!(frags, 2);
    }

    #[test]
    fn test_merged_coverage_with_overlap() {
        let h1 = Hit {
            query: "g".into(), query_length: 1000, query_start: 0, query_end: 600,
            strand: '+', target: "c".into(), target_start: 0, target_end: 600,
            matches: 600, alignment_length: 600, mapq: 60,
        };
        let h2 = Hit {
            query: "g".into(), query_length: 1000, query_start: 400, query_end: 1000,
            strand: '+', target: "c".into(), target_start: 400, target_end: 1000,
            matches: 600, alignment_length: 600, mapq: 60,
        };
        let (cov, frags) = merged_query_coverage(&[&h1, &h2]);
        assert!((cov - 100.0).abs() < 0.01);
        assert_eq!(frags, 1); // merged into one interval
    }

    #[test]
    fn test_gff_index_overlap() {
        let genes = vec![
            GffGene { contig: "chr1".into(), start: 100, end: 500, strand: '+', gene_id: "g1".into() },
            GffGene { contig: "chr1".into(), start: 600, end: 900, strand: '-', gene_id: "g2".into() },
            GffGene { contig: "chr2".into(), start: 50, end: 200, strand: '+', gene_id: "g3".into() },
        ];
        let idx = GffIndex::build(genes);

        // Query overlapping g1
        let hits = idx.overlapping("chr1", 200, 400);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].gene_id, "g1");

        // Query overlapping both g1 and g2
        let hits = idx.overlapping("chr1", 400, 700);
        assert_eq!(hits.len(), 2);

        // Query on chr2
        let hits = idx.overlapping("chr2", 0, 100);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].gene_id, "g3");

        // Query on non-existent contig
        let hits = idx.overlapping("chr3", 0, 100);
        assert!(hits.is_empty());
    }
}
