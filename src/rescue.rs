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
    /// Start site variant: gene is there but with a different start annotation.
    /// Detected via GFF: hit overlaps a CDS whose length differs from query.
    StartVariant { length_diff_bp: i64 },
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
            RescueType::StartVariant { length_diff_bp } => {
                write!(f, "START_VARIANT:{:+}bp", length_diff_bp)
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
    /// Confidence score: 0.0+. Used for conflict resolution.
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
    pub genes: Vec<GffGene>,
    contig_ranges: HashMap<String, (usize, usize)>,
    /// gene_id → index in `genes` for O(1) lookup by ID.
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

        // Build gene_id → index map
        let mut gene_id_index = HashMap::new();
        for (i, g) in genes.iter().enumerate() {
            if !g.gene_id.is_empty() {
                gene_id_index.insert(g.gene_id.clone(), i);
            }
        }

        GffIndex { genes, contig_ranges, gene_id_index }
    }

    /// Lookup a gene by ID. O(1).
    pub fn get_by_id(&self, gene_id: &str) -> Option<&GffGene> {
        self.gene_id_index.get(gene_id).map(|&i| &self.genes[i])
    }

    /// Find all genes overlapping [query_start, query_end) on a contig.
    pub fn overlapping(&self, contig: &str, query_start: u64, query_end: u64) -> Vec<&GffGene> {
        let Some(&(range_start, range_end)) = self.contig_ranges.get(contig) else {
            return Vec::new();
        };

        // Use bounded slice — partition_point requires monotonic predicate,
        // which only holds within a single contig's sorted genes.
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

/// Compute combined query coverage from multiple hits,
/// considering only the best allele (the one with the highest single-hit coverage).
///
/// Returns (coverage_pct, n_merged_intervals).
fn merged_query_coverage(hits: &[&Hit]) -> (f64, usize) {
    if hits.is_empty() {
        return (0.0, 0);
    }

    // Group hits by query_length (proxy for allele identity).
    // Pick the allele whose hits give the best merged coverage.
    let mut by_length: HashMap<u64, Vec<(u64, u64)>> = HashMap::new();
    for h in hits {
        by_length
            .entry(h.query_length)
            .or_default()
            .push((h.query_start, h.query_end));
    }

    let mut best_coverage = 0.0f64;
    let mut best_fragments = 0usize;

    for (qlen, mut intervals) in by_length {
        if qlen == 0 {
            continue;
        }
        intervals.sort_by_key(|&(s, _)| s);

        // Merge overlapping intervals
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

// ── Start-site variant detection (GFF-based) ──

/// Detect if a hit represents a start-site variant by comparing the hit
/// against the overlapping CDS from GFF. If the CDS has high reciprocal
/// overlap but differs in length by <300bp, it's likely a start variant.
///
/// Returns the length difference (positive = CDS longer than query, negative = shorter).
/// Returns (length_diff, gene_id) if a start variant is detected.
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

    for gene in &overlapping {
        let gene_len = (gene.end - gene.start) as i64;
        let query_len = hit.query_length as i64;
        let diff = gene_len - query_len;

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

    // Sort hits by score (best first). Use total_cmp to handle NaN safely.
    let mut sorted_hits: Vec<&Hit> = hits.iter().collect();
    sorted_hits.sort_by(|a, b| {
        b.score().total_cmp(&a.score())
    });

    let best = sorted_hits[0];

    // ── Path 1: Single full-coverage hit ──
    if best.query_coverage() >= min_coverage {
        // First classify via GFF (SPLIT vs DROPPED)
        let (rescue_type, annotation_match, confidence) = classify_single_hit(
            best,
            cluster_idx,
            pa,
            gff_index,
            gene_to_cluster,
        );

        // Only check START_VARIANT if the gene was classified as DROPPED
        // (not SPLIT — a split gene belongs to another cluster, the start doesn't matter)
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
            cluster_idx,
            genome_idx,
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
    // Group hits by contig and try fragment assembly per contig.
    // Mixing contigs would give nonsensical coordinates.
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
        let confidence = (avg_identity / 100.0) * (total_coverage / 100.0) * 0.8;

        let frag_start = same_contig_hits.iter().map(|h| h.target_start).min().unwrap();
        let frag_end = same_contig_hits.iter().map(|h| h.target_end).max().unwrap();

        return Some(RescueEvent {
            cluster_idx,
            genome_idx,
            rescue_type: RescueType::Fragmented {
                n_fragments,
                total_coverage,
            },
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
                base_confidence * 1.1,
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

            // Different genome → no more possible conflicts with i
            if ri.genome_idx != rj.genome_idx {
                break;
            }

            // Different contig → skip, there may be more on the same contig later
            // (impossible due to sort, but be safe — continue, not break)
            if ri.contig != rj.contig {
                // Sorted by contig, so all remaining same-contig entries for i are done.
                // But there could still be entries on the same contig as j that conflict with j.
                // For conflicts with i specifically: no more on this contig → break inner.
                break;
            }

            // Past the region → no more overlaps for i on this contig
            if rj.start >= ri.end {
                break;
            }

            // Overlapping → conflict. Higher confidence wins.
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

/// Validate a rescue by checking flanking gene conservation.
///
/// Compares Panaroo CLUSTER assignments (not gene IDs) of flanking genes.
/// Gene IDs are genome-specific and never match across genomes.
/// Instead: look up each flanking gene_id → gene_to_cluster → cluster_idx,
/// then compare cluster_idx between rescue genome and reference genome.
pub fn validate_synteny(
    rescue: &mut RescueEvent,
    pa: &PresenceAbsence,
    gff_indices: &HashMap<usize, GffIndex>,
    gene_to_cluster_maps: &HashMap<usize, HashMap<String, usize>>,
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

    // Map flanking gene_ids to cluster indices
    let target_up_cluster = target_up.and_then(|g| target_g2c.get(&g.gene_id).copied());
    let target_down_cluster = target_down.and_then(|g| target_g2c.get(&g.gene_id).copied());

    if target_up_cluster.is_none() && target_down_cluster.is_none() {
        return; // can't validate
    }

    // Check against reference genomes where this cluster IS present
    let mut synteny_match = false;
    for ref_gi in pa.present[ci].iter().take(5) {
        let Some(ref_gff) = gff_indices.get(ref_gi) else {
            continue;
        };
        let Some(ref_g2c) = gene_to_cluster_maps.get(ref_gi) else {
            continue;
        };

        // Find the gene's position in the reference genome
        let ref_cell = &pa.cell_values[ci][*ref_gi];
        if ref_cell.is_empty() {
            continue;
        }
        let ref_gene_id = ref_cell.split(';').next().unwrap_or(ref_cell).trim();

        // Locate gene in reference GFF (O(1) via index)
        let Some(ref_gene) = ref_gff.get_by_id(ref_gene_id) else {
            continue;
        };

        let (ref_up, ref_down) = ref_gff.flanking(&ref_gene.contig, ref_gene.start, ref_gene.end);

        // Map reference flanking genes to cluster indices
        let ref_up_cluster = ref_up.and_then(|g| ref_g2c.get(&g.gene_id).copied());
        let ref_down_cluster = ref_down.and_then(|g| ref_g2c.get(&g.gene_id).copied());

        // Compare CLUSTER indices, not gene IDs
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
        // Two alleles: 800bp and 1200bp. Each has partial hits.
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
        // Should pick allele 1200bp (91.7% coverage, 1 fragment) over 800bp (75%, 1 fragment)
        assert!(cov > 90.0);
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
        assert_eq!(up.unwrap().gene_id, "g1");
        assert_eq!(down.unwrap().gene_id, "g3");

        let (up, down) = idx.flanking("c", 100, 200);
        assert!(up.is_none());
        assert_eq!(down.unwrap().gene_id, "g2");
    }

    #[test]
    fn test_detect_start_variant_with_gff() {
        let genes = vec![
            // CDS is 1100bp, query is 1000bp → 100bp start variant
            GffGene { contig: "c".into(), start: 0, end: 1100, strand: '+', gene_id: "gene1".into() },
        ];
        let idx = GffIndex::build(genes);

        let hit = Hit {
            query: "g".into(), cluster: "g".into(), query_length: 1000,
            query_start: 0, query_end: 1000, strand: '+', target: "c".into(),
            target_start: 100, target_end: 1100, // maps to end of CDS
            matches: 1000, alignment_length: 1000, mapq: 60,
        };
        let result = detect_start_variant(&hit, Some(&idx));
        assert_eq!(result, Some((100, "gene1".to_string()))); // CDS is 100bp longer
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
        assert_eq!(rescues[0].cluster_idx, 0);
    }

    #[test]
    fn test_conflict_resolution_no_overlap() {
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
                start: 600, end: 900, strand: '+', annotation_match: None,
                confidence: 0.85, synteny_ok: None,
            },
        ];
        resolve_conflicts(&mut rescues);
        assert_eq!(rescues.len(), 2); // no overlap → both kept
    }

    #[test]
    fn test_fragmented_span() {
        // Two non-overlapping fragments should report min_start to max_end
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

        use ahash::{HashSet, HashSetExt};
        let pa = PresenceAbsence {
            clusters: vec!["g".into()],
            genomes: vec!["genome1".into()],
            present: vec![HashSet::new()], // g absent in genome1
            cell_values: vec![vec!["".into()]],
            header: vec!["Gene".into(), "genome1".into()],
            metadata: vec![vec!["g".into()]],
            genome_col_indices: vec![1],
        };
        let g2c: HashMap<String, usize> = HashMap::new();

        let result = classify_hits(&[h1, h2], 0, 0, &pa, None, &g2c, 80.0);
        assert!(result.is_some());
        let r = result.unwrap();
        assert_eq!(r.start, 1000); // min of all fragment starts
        assert_eq!(r.end, 2400);   // max of all fragment ends
        match r.rescue_type {
            RescueType::Fragmented { n_fragments, .. } => assert_eq!(n_fragments, 2),
            _ => panic!("Expected Fragmented"),
        }
    }
}
