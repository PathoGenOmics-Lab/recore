use crate::parser::PresenceAbsence;
use crate::search::Hit;

/// Classification of a rescued gene.
#[derive(Debug, Clone, PartialEq)]
pub enum RescueType {
    /// Gene was dropped during Panaroo's graph cleaning.
    Dropped,
    /// Gene was assigned to a different cluster (split).
    Split { original_cluster: String },
    /// Gene overlaps with a paralog in a different cluster.
    #[allow(dead_code)]
    Overlap { overlapping_cluster: String },
}

impl std::fmt::Display for RescueType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RescueType::Dropped => write!(f, "DROPPED"),
            RescueType::Split { original_cluster } => write!(f, "SPLIT:{}", original_cluster),
            RescueType::Overlap { overlapping_cluster } => {
                write!(f, "OVERLAP:{}", overlapping_cluster)
            }
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

/// Classify a hit as DROPPED, SPLIT, or OVERLAP.
///
/// A hit is classified as SPLIT if the genomic region is already occupied
/// by another cluster (i.e., the gene was mis-assigned rather than truly absent).
///
/// For coordinate-based split detection, we need per-gene coordinates from
/// the Panaroo output or GFF. Without them, we can check if any other cluster
/// in this genome shares the same gene ID prefix.
pub fn classify_hit(
    hit: &Hit,
    cluster_idx: usize,
    genome_idx: usize,
    pa: &PresenceAbsence,
    gff_genes: Option<&[(String, u64, u64, String)]>,
) -> RescueEvent {
    let mut rescue_type = RescueType::Dropped;
    let mut annotation_match = None;

    // Check GFF overlap if available
    if let Some(genes) = gff_genes {
        for (contig, gstart, gend, gene_id) in genes {
            if *contig == hit.target && hit.target_start < *gend && hit.target_end > *gstart {
                let overlap_start = hit.target_start.max(*gstart);
                let overlap_end = hit.target_end.min(*gend);
                let overlap_len = overlap_end.saturating_sub(overlap_start);
                let hit_len = hit.target_end - hit.target_start;
                if hit_len > 0 && (overlap_len as f64 / hit_len as f64) > 0.5 {
                    annotation_match = Some(gene_id.clone());
                    break;
                }
            }
        }
    }

    // Check if the gene ID in the hit region is assigned to another cluster
    // by looking at what genes this genome has in other clusters.
    // This is a heuristic: if any other cluster in this genome has a gene
    // and we have a strong hit in the same region, it's likely a SPLIT.
    if let Some(ref ann) = annotation_match {
        for (ci, present_set) in pa.present.iter().enumerate() {
            if ci == cluster_idx {
                continue;
            }
            if !present_set.contains(&genome_idx) {
                continue;
            }
            // Check if the gene_id assigned to this other cluster matches
            let cell = &pa.cell_values[ci][genome_idx];
            if !cell.is_empty() && ann.contains(cell.as_str()) {
                rescue_type = RescueType::Split {
                    original_cluster: pa.clusters[ci].clone(),
                };
                break;
            }
        }
    }

    RescueEvent {
        cluster_idx,
        genome_idx,
        rescue_type,
        identity: hit.identity(),
        coverage: hit.query_coverage(),
        contig: hit.target.clone(),
        start: hit.target_start,
        end: hit.target_end,
        strand: hit.strand,
        annotation_match,
    }
}
