use super::*;

#[derive(Debug, Clone)]
pub(super) struct OrthologyQuality {
    pub(super) adaptive_identity_floor: f64,
    pub(super) adaptive_coverage_floor: f64,
    pub(super) adaptive_sequence_pass: bool,
    pub(super) integrated_score: f64,
    pub(super) tier: &'static str,
    pub(super) risk_flags: String,
}

pub(super) fn cluster_sequence_lengths(
    pa: &crate::parser::PresenceAbsence,
    pan_ref: &HashMap<String, Vec<Vec<u8>>>,
) -> Vec<usize> {
    (0..pa.n_clusters())
        .map(|ci_val| {
            pan_ref
                .get(pa.cluster_name(ClusterIdx(ci_val)))
                .and_then(|alleles| alleles.iter().map(Vec::len).max())
                .unwrap_or(0)
        })
        .collect()
}

pub(super) fn cluster_prevalences(pa: &crate::parser::PresenceAbsence) -> Vec<f64> {
    let n_genomes = pa.n_genomes();
    (0..pa.n_clusters())
        .map(|ci_val| {
            if n_genomes == 0 {
                0.0
            } else {
                pa.present_count(ClusterIdx(ci_val)) as f64 / n_genomes as f64
            }
        })
        .collect()
}

fn adaptive_sequence_floors(cluster_prevalence: f64, cluster_seq_len: usize) -> (f64, f64) {
    let mut identity_floor: f64 = 90.0;
    let mut coverage_floor: f64 = 80.0;

    if cluster_prevalence < 0.03 {
        identity_floor = 97.0;
        coverage_floor = 90.0;
    } else if cluster_prevalence < 0.10 {
        identity_floor = 95.0;
        coverage_floor = 85.0;
    } else if cluster_prevalence < 0.20 {
        identity_floor = 93.0;
        coverage_floor = 83.0;
    }

    if cluster_seq_len > 0 && cluster_seq_len < 150 {
        identity_floor += 2.0;
        coverage_floor += 5.0;
    } else if cluster_seq_len > 0 && cluster_seq_len < 300 {
        identity_floor += 1.0;
        coverage_floor += 3.0;
    }

    (identity_floor.min(99.0), coverage_floor.min(95.0))
}

fn rescue_orf_penalty(rescue: &RescueEvent) -> f64 {
    match rescue.rescue_type {
        RescueType::Pseudogene { .. } => -0.25,
        RescueType::Fragmented { .. } | RescueType::AlignmentSplit { .. } => -0.15,
        _ => {
            let span = rescue.end.saturating_sub(rescue.start);
            if span > 0 && span % 3 != 0 {
                -0.05
            } else {
                0.05
            }
        }
    }
}

fn state_adjustment(state: &str) -> f64 {
    match state {
        "STRICT_PRESENT" => 0.18,
        "ALIAS_CLUSTER" => 0.14,
        "SEQUENCE_ONLY" => -0.03,
        "MISSING_CDS_UNVALIDATED" => -0.12,
        "PARTIAL_HOMOLOGY" | "PSEUDOGENE" => -0.22,
        "HOMOLOG_PARALOG" => -0.35,
        "REJECTED_CONTEXT" => -0.45,
        _ => -0.10,
    }
}

fn synteny_adjustment(rescue: &RescueEvent) -> f64 {
    match rescue.synteny_ok {
        Some(true) => {
            let score_bonus = rescue
                .synteny_evidence
                .score
                .map(|score| (score - 0.5).max(0.0) * 0.20)
                .unwrap_or(0.0);
            0.20 + score_bonus
        }
        Some(false) => -0.40,
        None => 0.0,
    }
}

fn quality_tier(accepted: bool, state: &str, adaptive_pass: bool, score: f64) -> &'static str {
    if !accepted || state == "REJECTED_CONTEXT" {
        return "rejected";
    }
    if score >= 0.85 && adaptive_pass && matches!(state, "STRICT_PRESENT" | "ALIAS_CLUSTER") {
        "high_confidence"
    } else if score >= 0.70 && adaptive_pass {
        "moderate_confidence"
    } else if score >= 0.50 {
        "review"
    } else {
        "low_confidence"
    }
}

pub(super) fn orthology_quality(
    rescue: &RescueEvent,
    state: &str,
    accepted: bool,
    cluster_prevalence: f64,
    cluster_seq_len: usize,
) -> OrthologyQuality {
    let (identity_floor, coverage_floor) =
        adaptive_sequence_floors(cluster_prevalence, cluster_seq_len);
    let adaptive_sequence_pass =
        rescue.identity >= identity_floor && rescue.coverage >= coverage_floor;
    let same_locus_override = accepted
        && (rescue.synteny_ok == Some(true) || matches!(state, "STRICT_PRESENT" | "ALIAS_CLUSTER"));
    let calibrated_sequence = adaptive_sequence_pass || same_locus_override;

    let sequence_score = ((rescue.identity / 100.0) * (rescue.coverage / 100.0)).clamp(0.0, 1.0);
    let mut score = sequence_score
        + synteny_adjustment(rescue)
        + state_adjustment(state)
        + rescue_orf_penalty(rescue);
    if !calibrated_sequence {
        score -= 0.20;
    }
    score = score.clamp(0.0, 1.0);

    let mut flags = Vec::new();
    if !adaptive_sequence_pass {
        flags.push("adaptive_sequence_below_floor");
    }
    if cluster_seq_len > 0 && cluster_seq_len < 300 {
        flags.push("short_cluster_sequence");
    }
    match rescue.synteny_ok {
        Some(true) => {}
        Some(false) => flags.push("incompatible_synteny"),
        None => flags.push("synteny_not_tested"),
    }
    match rescue.rescue_type {
        RescueType::Paralog { .. } => flags.push("paralog_signal"),
        RescueType::Fragmented { .. } | RescueType::AlignmentSplit { .. } => {
            flags.push("fragmented_or_partial")
        }
        RescueType::Pseudogene { .. } => flags.push("disrupted_orf"),
        RescueType::MissingCds if rescue.synteny_ok != Some(true) => {
            flags.push("missing_cds_without_synteny")
        }
        _ => {}
    }

    OrthologyQuality {
        adaptive_identity_floor: identity_floor,
        adaptive_coverage_floor: coverage_floor,
        adaptive_sequence_pass,
        integrated_score: score,
        tier: quality_tier(accepted, state, calibrated_sequence, score),
        risk_flags: flags.join(";"),
    }
}

pub(super) fn quality_allows_matrix_correction(state: &str, quality: &OrthologyQuality) -> bool {
    match state {
        "STRICT_PRESENT" | "ALIAS_CLUSTER" => quality.integrated_score >= 0.55,
        "SEQUENCE_ONLY" => {
            quality.adaptive_sequence_pass
                && matches!(
                    quality.tier,
                    "high_confidence" | "moderate_confidence" | "review"
                )
                && quality.integrated_score >= 0.60
        }
        "PARTIAL_HOMOLOGY" | "PSEUDOGENE" => {
            quality.adaptive_sequence_pass && quality.integrated_score >= 0.50
        }
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rescue::FragmentEvidence;

    fn rescue(identity: f64, coverage: f64, rescue_type: RescueType) -> RescueEvent {
        RescueEvent {
            cluster_idx: ClusterIdx(0),
            genome_idx: GenomeIdx(0),
            rescue_type,
            identity,
            coverage,
            contig: "ctg".to_string(),
            start: 0,
            end: 300,
            strand: '+',
            annotation_match: None,
            confidence: identity / 100.0 * coverage / 100.0,
            synteny_ok: None,
            synteny_evidence: SyntenyEvidence::default(),
            fragment_evidence: FragmentEvidence::default(),
            query_start: None,
        }
    }

    #[test]
    fn rare_short_sequence_requires_stricter_sequence_support() {
        let r = rescue(94.0, 86.0, RescueType::Dropped);
        let q = orthology_quality(&r, "SEQUENCE_ONLY", true, 0.02, 120);
        assert!(!q.adaptive_sequence_pass);
        assert!(q.risk_flags.contains("adaptive_sequence_below_floor"));
        assert!(!quality_allows_matrix_correction("SEQUENCE_ONLY", &q));
    }

    #[test]
    fn synteny_supported_strict_presence_can_pass_same_locus_override() {
        let mut r = rescue(92.0, 82.0, RescueType::Dropped);
        r.synteny_ok = Some(true);
        r.synteny_evidence.score = Some(0.9);
        let q = orthology_quality(&r, "STRICT_PRESENT", true, 0.02, 120);
        assert!(!q.adaptive_sequence_pass);
        assert!(quality_allows_matrix_correction("STRICT_PRESENT", &q));
    }
}
