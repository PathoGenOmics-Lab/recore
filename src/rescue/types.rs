use crate::types::{ClusterIdx, GenomeIdx};

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
    /// Where the single alignment this rescue rests on starts on the query,
    /// the cluster's representative CDS. The representative begins at its
    /// start codon, so this fixes the reading frame of the span on the contig
    /// (and the phase of the feature in a corrected GFF). `None` when the
    /// rescue was built from several alignments or the position is unknown.
    pub query_start: Option<u64>,
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
