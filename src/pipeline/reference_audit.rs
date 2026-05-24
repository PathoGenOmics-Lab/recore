use super::*;

pub(super) fn resolve_reference_genome(
    pa: &crate::parser::PresenceAbsence,
    requested: &str,
) -> Option<GenomeIdx> {
    for gi_val in 0..pa.n_genomes() {
        let gi = GenomeIdx(gi_val);
        if pa.genome_name(gi) == requested {
            return Some(gi);
        }
    }

    let requested_as = extract_as_id(requested);
    for gi_val in 0..pa.n_genomes() {
        let gi = GenomeIdx(gi_val);
        let genome_name = pa.genome_name(gi);
        if let (Some(a), Some(b)) = (requested_as, extract_as_id(genome_name)) {
            if a == b {
                return Some(gi);
            }
        }
    }
    None
}

pub(super) fn reference_status_for(state: &str, accepted: bool) -> &'static str {
    if !accepted {
        return match state {
            "HOMOLOG_PARALOG" => "reference_homolog_paralog_rejected",
            "PARTIAL_PARALOG_EVIDENCE" => "reference_partial_paralog_evidence",
            "EVENT_DERIVED_REFERENCE_CLUSTER" => "reference_event_derived_cluster",
            "REJECTED_CONTEXT" => "reference_incompatible_context",
            _ => "reference_rejected_candidate",
        };
    }
    match state {
        "STRICT_PRESENT" => "reference_rescue_strict_present",
        "ALIAS_CLUSTER" => "reference_alias_cluster",
        "SEQUENCE_ONLY" => "reference_sequence_only",
        "MISSING_CDS_UNVALIDATED" => "reference_missing_cds_unvalidated",
        "PARTIAL_HOMOLOGY" => "reference_partial_like",
        "PARTIAL_PARALOG_EVIDENCE" => "reference_partial_paralog_evidence",
        "EVENT_DERIVED_REFERENCE_CLUSTER" => "reference_event_derived_cluster",
        "PSEUDOGENE" => "reference_pseudogene_like",
        "HOMOLOG_PARALOG" => "reference_homolog_paralog",
        _ => "reference_candidate_present",
    }
}

pub(super) fn reference_publication_status(state: &str, accepted: bool) -> &'static str {
    if !accepted {
        if state == "EVENT_DERIVED_REFERENCE_CLUSTER" {
            return "flag_reference_event_derived_cluster";
        }
        if state == "PARTIAL_PARALOG_EVIDENCE" {
            return "flag_reference_partial_paralog_evidence";
        }
        return "flag_reference_like_ambiguous";
    }
    match state {
        "STRICT_PRESENT" | "ALIAS_CLUSTER" | "SEQUENCE_ONLY" => {
            "remove_from_strict_reference_absent_set"
        }
        "MISSING_CDS_UNVALIDATED" | "PARTIAL_HOMOLOGY" | "PSEUDOGENE" | "HOMOLOG_PARALOG" => {
            "flag_reference_like_ambiguous"
        }
        "PARTIAL_PARALOG_EVIDENCE" => "flag_reference_partial_paralog_evidence",
        "EVENT_DERIVED_REFERENCE_CLUSTER" => "flag_reference_event_derived_cluster",
        _ => "review_reference_status",
    }
}

pub(super) fn reference_recommended_wording(state: &str, accepted: bool) -> &'static str {
    if !accepted {
        if state == "EVENT_DERIVED_REFERENCE_CLUSTER" {
            return "Do not call a clean novel gene; keep the cluster separate and report the inferred reference-gene event.";
        }
        if state == "PARTIAL_PARALOG_EVIDENCE" {
            return "Do not call cleanly novel to the reference; report as partial/paralogous reference-like evidence.";
        }
        return "Do not call cleanly novel to the reference; report as reference-like but rejected or ambiguous evidence.";
    }
    match state {
        "STRICT_PRESENT" | "SEQUENCE_ONLY" => {
            "Remove from strict reference-absent set; the audited reference has sequence evidence."
        }
        "ALIAS_CLUSTER" => {
            "Remove from strict reference-absent set; evidence points to a cluster alias or split."
        }
        "MISSING_CDS_UNVALIDATED" => {
            "Reference-like sequence is present but lacks annotation/synteny-grade support; review before claiming clean absence."
        }
        "PARTIAL_HOMOLOGY" => {
            "Describe as partial reference-like homology rather than fully novel to the reference."
        }
        "PSEUDOGENE" => {
            "Describe as disrupted reference-like sequence rather than fully novel to the reference."
        }
        "HOMOLOG_PARALOG" => {
            "Keep as reference-like paralog evidence; do not treat as clean absence without review."
        }
        "PARTIAL_PARALOG_EVIDENCE" => {
            "Reference has partial/paralogous sequence evidence; keep the exact cluster absent but flag as complex."
        }
        "EVENT_DERIVED_REFERENCE_CLUSTER" => {
            "Exact cluster is absent, but representative sequence supports an event-derived variant of a reference gene."
        }
        _ => "Review before making a reference-absence claim.",
    }
}

pub(super) fn reference_publication_confidence(
    state: &str,
    reference_status: &str,
) -> &'static str {
    match state {
        "STRICT_PRESENT" | "ALIAS_CLUSTER" | "SEQUENCE_ONLY" => "not_reference_absent",
        "PARTIAL_PARALOG_EVIDENCE" | "HOMOLOG_PARALOG" => "complex_paralog",
        "EVENT_DERIVED_REFERENCE_CLUSTER" => "event_derived",
        "PARTIAL_HOMOLOGY" => "fragmented",
        "PSEUDOGENE" => "pseudogene",
        "NO_EVIDENCE" => "clean",
        "REJECTED_CONTEXT" | "MISSING_CDS_UNVALIDATED" => "unresolved",
        _ if reference_status.contains("paralog") => "complex_paralog",
        _ => "unresolved",
    }
}

pub(super) fn reference_shadow_state(
    rescue: &RescueEvent,
    query_cluster: &str,
    annotation_cluster: &str,
) -> &'static str {
    if matches!(rescue.rescue_type, RescueType::Paralog { .. }) {
        return "PARTIAL_PARALOG_EVIDENCE";
    }
    if matches!(
        rescue.rescue_type,
        RescueType::Pseudogene { .. }
            | RescueType::Fragmented { .. }
            | RescueType::AlignmentSplit { .. }
    ) {
        return evidence_state(rescue, true, query_cluster, annotation_cluster, "");
    }
    if !annotation_cluster.is_empty() && annotation_cluster != query_cluster {
        return "PARTIAL_PARALOG_EVIDENCE";
    }
    if rescue.identity >= 90.0 && rescue.coverage < 95.0 {
        return "PARTIAL_PARALOG_EVIDENCE";
    }
    evidence_state(rescue, true, query_cluster, annotation_cluster, "")
}

pub(super) fn best_overlap_gene_id(gff: &GffIndex, hit: &Hit) -> Option<String> {
    let hit_len = hit.target_end.saturating_sub(hit.target_start) as f64;
    if hit_len == 0.0 {
        return None;
    }

    gff.overlapping(&hit.target, hit.target_start, hit.target_end)
        .into_iter()
        .filter_map(|gene| {
            let overlap = hit
                .target_end
                .min(gene.end)
                .saturating_sub(hit.target_start.max(gene.start)) as f64;
            let gene_len = gene.end.saturating_sub(gene.start) as f64;
            if gene_len == 0.0 {
                return None;
            }
            let reciprocal = (overlap / hit_len).min(overlap / gene_len);
            (reciprocal > 0.5).then(|| (gene.gene_id.clone(), reciprocal))
        })
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(gene_id, _)| gene_id)
}

pub(super) fn reference_hit_shadow_candidates(
    hits: &HashMap<(GenomeIdx, ClusterIdx), Vec<Hit>>,
    pa: &crate::parser::PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    gene_to_cluster_maps: &HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
) -> Vec<RescueEvent> {
    let mut candidates = Vec::new();
    for (&(gi, ci), cluster_hits) in hits {
        let gff = gff_indices.get(&gi);
        let g2c = gene_to_cluster_maps.get(&gi);
        for hit in cluster_hits {
            if hit.identity() < 90.0 || hit.query_coverage() < 60.0 {
                continue;
            }
            let annotation_match = gff.and_then(|idx| best_overlap_gene_id(idx, hit));
            let annotation_cluster = annotation_match
                .as_deref()
                .and_then(|gene_id| g2c.and_then(|map| map.get(gene_id)))
                .map(|other_ci| pa.cluster_name(*other_ci).to_string());
            let query_cluster = pa.cluster_name(ci);
            let rescue_type = if let Some(other_cluster) = annotation_cluster {
                if other_cluster != query_cluster {
                    RescueType::Paralog { other_cluster }
                } else if hit.query_coverage() < 95.0 {
                    RescueType::Fragmented {
                        n_fragments: 1,
                        total_coverage: hit.query_coverage(),
                    }
                } else {
                    RescueType::Dropped
                }
            } else if hit.query_coverage() < 95.0 {
                RescueType::Fragmented {
                    n_fragments: 1,
                    total_coverage: hit.query_coverage(),
                }
            } else {
                RescueType::Dropped
            };

            candidates.push(RescueEvent {
                cluster_idx: ci,
                genome_idx: gi,
                rescue_type,
                identity: hit.identity(),
                coverage: hit.query_coverage(),
                contig: hit.target.clone(),
                start: hit.target_start,
                end: hit.target_end,
                strand: hit.strand,
                annotation_match,
                confidence: hit.score(),
                synteny_ok: None,
                synteny_evidence: SyntenyEvidence::default(),
                fragment_evidence: crate::rescue::FragmentEvidence {
                    note: Some("reference_shadow_hit".to_string()),
                    ..crate::rescue::FragmentEvidence::default()
                },
            });
        }
    }
    candidates
}

pub(super) fn reference_candidate_score(rescue: &RescueEvent) -> f64 {
    let synteny_bonus = match rescue.synteny_ok {
        Some(true) => 2.0,
        Some(false) => -1.0,
        None => 0.0,
    };
    rescue.confidence + (rescue.identity / 100.0) + (rescue.coverage / 100.0) + synteny_bonus
}

pub(super) fn best_reference_candidate(
    candidates: &[RescueEvent],
    ci: ClusterIdx,
    reference_gi: GenomeIdx,
) -> Option<&RescueEvent> {
    candidates
        .iter()
        .filter(|r| r.cluster_idx == ci && r.genome_idx == reference_gi)
        .max_by(|a, b| reference_candidate_score(a).total_cmp(&reference_candidate_score(b)))
}

pub(super) const REFERENCE_EVENT_KMER: usize = 15;
pub(super) const REFERENCE_EVENT_TOP_CANDIDATES: usize = 12;
pub(super) const REFERENCE_EVENT_MAX_DP_CELLS: usize = 5_000_000;

pub(super) struct ReferenceClusterEvent {
    pub(super) query_cluster_idx: ClusterIdx,
    pub(super) target_cluster_idx: ClusterIdx,
    pub(super) target_gene_id: String,
    pub(super) event_class: &'static str,
    pub(super) identity: f64,
    pub(super) query_coverage: f64,
    pub(super) reference_coverage: f64,
    pub(super) query_length: usize,
    pub(super) reference_length: usize,
    pub(super) event_size_nt: i64,
    pub(super) event_size_aa: Option<i64>,
    pub(super) matches: usize,
    pub(super) mismatches: usize,
    pub(super) gaps_in_query: usize,
    pub(super) gaps_in_reference: usize,
    pub(super) shared_kmers: usize,
    pub(super) interpretation: &'static str,
}
pub(super) struct PairwiseEventStats {
    identity: f64,
    query_coverage: f64,
    reference_coverage: f64,
    matches: usize,
    mismatches: usize,
    gaps_in_query: usize,
    gaps_in_reference: usize,
}

pub(super) fn longest_cluster_sequence<'a>(
    pan_ref: &'a HashMap<String, Vec<Vec<u8>>>,
    cluster_name: &str,
) -> Option<&'a [u8]> {
    pan_ref
        .get(cluster_name)?
        .iter()
        .max_by_key(|seq| seq.len())
        .map(Vec::as_slice)
}

pub(super) fn base_bits(base: u8) -> Option<u64> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

pub(super) fn canonical_kmer_hash(seq: &[u8]) -> Option<u64> {
    let mut forward = 0u64;
    let mut reverse = 0u64;
    for &base in seq {
        forward = (forward << 2) | base_bits(base)?;
    }
    for &base in seq.iter().rev() {
        let bits = base_bits(base)?;
        reverse = (reverse << 2) | (3 - bits);
    }
    Some(forward.min(reverse))
}

pub(super) fn sequence_kmers(seq: &[u8], k: usize) -> Vec<u64> {
    if seq.len() < k {
        return Vec::new();
    }
    let mut kmers: Vec<u64> = (0..=seq.len() - k)
        .filter_map(|i| canonical_kmer_hash(&seq[i..i + k]))
        .collect();
    kmers.sort_unstable();
    kmers.dedup();
    kmers
}

pub(super) fn global_event_alignment(query: &[u8], reference: &[u8]) -> Option<PairwiseEventStats> {
    let m = query.len();
    let n = reference.len();
    if m == 0 || n == 0 || (m + 1).saturating_mul(n + 1) > REFERENCE_EVENT_MAX_DP_CELLS {
        return None;
    }

    let width = n + 1;
    let idx = |i: usize, j: usize| i * width + j;
    let mut scores = vec![0i32; (m + 1) * (n + 1)];
    let mut trace = vec![0u8; (m + 1) * (n + 1)];

    for i in 1..=m {
        scores[idx(i, 0)] = scores[idx(i - 1, 0)] - 2;
        trace[idx(i, 0)] = 1;
    }
    for j in 1..=n {
        scores[idx(0, j)] = scores[idx(0, j - 1)] - 2;
        trace[idx(0, j)] = 2;
    }

    for i in 1..=m {
        for j in 1..=n {
            let q = query[i - 1].to_ascii_uppercase();
            let r = reference[j - 1].to_ascii_uppercase();
            let diag = scores[idx(i - 1, j - 1)] + if q == r { 2 } else { -1 };
            let up = scores[idx(i - 1, j)] - 2;
            let left = scores[idx(i, j - 1)] - 2;
            if diag >= up && diag >= left {
                scores[idx(i, j)] = diag;
                trace[idx(i, j)] = 0;
            } else if up >= left {
                scores[idx(i, j)] = up;
                trace[idx(i, j)] = 1;
            } else {
                scores[idx(i, j)] = left;
                trace[idx(i, j)] = 2;
            }
        }
    }

    let mut i = m;
    let mut j = n;
    let mut matches = 0usize;
    let mut mismatches = 0usize;
    let mut gaps_in_query = 0usize;
    let mut gaps_in_reference = 0usize;
    while i > 0 || j > 0 {
        match trace[idx(i, j)] {
            0 if i > 0 && j > 0 => {
                if query[i - 1].eq_ignore_ascii_case(&reference[j - 1]) {
                    matches += 1;
                } else {
                    mismatches += 1;
                }
                i -= 1;
                j -= 1;
            }
            1 if i > 0 => {
                gaps_in_reference += 1;
                i -= 1;
            }
            2 if j > 0 => {
                gaps_in_query += 1;
                j -= 1;
            }
            _ if i > 0 => {
                gaps_in_reference += 1;
                i -= 1;
            }
            _ if j > 0 => {
                gaps_in_query += 1;
                j -= 1;
            }
            _ => break,
        }
    }

    let aligned_non_gap = matches + mismatches;
    if aligned_non_gap == 0 {
        return None;
    }
    Some(PairwiseEventStats {
        identity: matches as f64 / aligned_non_gap as f64 * 100.0,
        query_coverage: aligned_non_gap as f64 / m as f64 * 100.0,
        reference_coverage: aligned_non_gap as f64 / n as f64 * 100.0,
        matches,
        mismatches,
        gaps_in_query,
        gaps_in_reference,
    })
}

pub(super) fn classify_reference_event(
    query_len: usize,
    reference_len: usize,
    stats: &PairwiseEventStats,
) -> Option<(&'static str, i64, Option<i64>, &'static str)> {
    if stats.identity < 90.0 || stats.query_coverage < 80.0 || stats.reference_coverage < 40.0 {
        return None;
    }

    let size = query_len as i64 - reference_len as i64;
    let abs_size = size.abs();
    if stats.identity >= 95.0 && stats.query_coverage >= 90.0 && stats.reference_coverage >= 90.0 {
        return Some((
            "same_gene_full_length",
            size,
            (size % 3 == 0).then_some(size / 3),
            "Same-gene sequence evidence; not a clean reference-novel cluster.",
        ));
    }
    if query_len < reference_len && stats.identity >= 95.0 && stats.query_coverage >= 90.0 {
        if abs_size % 3 == 0 {
            return Some((
                "short_in_frame_allele",
                -abs_size,
                Some(-(abs_size / 3)),
                "Short in-frame allele of a reference gene; keep separate in the matrix but do not call a clean novel gene.",
            ));
        }
        return Some((
            "frameshift_or_truncated_reference_homology",
            -abs_size,
            None,
            "Reference-like sequence with an out-of-frame length difference; review as possible frameshift/truncation.",
        ));
    }
    if query_len > reference_len && stats.identity >= 95.0 && stats.reference_coverage >= 90.0 {
        if abs_size % 3 == 0 {
            return Some((
                "long_in_frame_allele",
                abs_size,
                Some(abs_size / 3),
                "Long in-frame allele of a reference gene; keep separate in the matrix but do not call a clean novel gene.",
            ));
        }
        return Some((
            "frameshift_or_extension_reference_homology",
            abs_size,
            None,
            "Reference-like sequence with an out-of-frame extension; review before biological interpretation.",
        ));
    }
    Some((
        "partial_reference_homology",
        size,
        (size % 3 == 0).then_some(size / 3),
        "Partial reference-gene homology; exact cluster is separate, but the sequence is not cleanly novel.",
    ))
}

pub(super) fn is_h37rv_locus_token(token: &str) -> bool {
    let Some(rest) = token.strip_prefix("Rv") else {
        return false;
    };
    if rest.len() < 4 || !rest[..4].chars().all(|ch| ch.is_ascii_digit()) {
        return false;
    }
    rest[4..].chars().all(|ch| ch.is_ascii_alphabetic())
}

pub(super) fn display_reference_gene_id(gene_id: &str) -> String {
    let cleaned = stripped_gene_id(gene_id.trim()).unwrap_or(gene_id.trim());
    cleaned
        .split('_')
        .find(|part| is_h37rv_locus_token(part))
        .unwrap_or(cleaned)
        .to_string()
}

pub(super) fn reference_gene_id(
    pa: &crate::parser::PresenceAbsence,
    ci: ClusterIdx,
    gi: GenomeIdx,
) -> String {
    let gene_id = pa
        .cell(ci, gi)
        .split(';')
        .map(str::trim)
        .find(|gene| !gene.is_empty())
        .unwrap_or("");
    let display = display_reference_gene_id(gene_id);
    if extract_as_id(&display).is_some() && display == pa.genome_name(gi) {
        String::new()
    } else {
        display
    }
}

pub(super) fn build_reference_cluster_events(
    pa: &crate::parser::PresenceAbsence,
    pan_ref: &HashMap<String, Vec<Vec<u8>>>,
    reference_gi: GenomeIdx,
) -> Vec<ReferenceClusterEvent> {
    let mut kmer_index: HashMap<u64, Vec<ClusterIdx>> = HashMap::new();
    let mut reference_clusters = Vec::new();
    for ci_val in 0..pa.n_clusters() {
        let ci = ClusterIdx(ci_val);
        if !pa.is_present(ci, reference_gi) {
            continue;
        }
        let Some(seq) = longest_cluster_sequence(pan_ref, pa.cluster_name(ci)) else {
            continue;
        };
        if seq.len() < REFERENCE_EVENT_KMER {
            continue;
        }
        reference_clusters.push(ci);
        for kmer in sequence_kmers(seq, REFERENCE_EVENT_KMER) {
            kmer_index.entry(kmer).or_default().push(ci);
        }
    }

    let mut events = Vec::new();
    for ci_val in 0..pa.n_clusters() {
        let query_ci = ClusterIdx(ci_val);
        if pa.is_present(query_ci, reference_gi) {
            continue;
        }
        let Some(query_seq) = longest_cluster_sequence(pan_ref, pa.cluster_name(query_ci)) else {
            continue;
        };
        if query_seq.len() < REFERENCE_EVENT_KMER {
            continue;
        }

        let mut candidate_counts: HashMap<ClusterIdx, usize> = HashMap::new();
        for kmer in sequence_kmers(query_seq, REFERENCE_EVENT_KMER) {
            if let Some(candidates) = kmer_index.get(&kmer) {
                for &target_ci in candidates {
                    *candidate_counts.entry(target_ci).or_default() += 1;
                }
            }
        }
        if candidate_counts.is_empty() {
            continue;
        }

        let mut candidates: Vec<(ClusterIdx, usize)> = candidate_counts.into_iter().collect();
        candidates.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        candidates.truncate(REFERENCE_EVENT_TOP_CANDIDATES.min(reference_clusters.len()));

        let mut best_event: Option<ReferenceClusterEvent> = None;
        for (target_ci, shared_kmers) in candidates {
            let Some(reference_seq) = longest_cluster_sequence(pan_ref, pa.cluster_name(target_ci))
            else {
                continue;
            };
            let Some(stats) = global_event_alignment(query_seq, reference_seq) else {
                continue;
            };
            let Some((event_class, event_size_nt, event_size_aa, interpretation)) =
                classify_reference_event(query_seq.len(), reference_seq.len(), &stats)
            else {
                continue;
            };
            let candidate_event = ReferenceClusterEvent {
                query_cluster_idx: query_ci,
                target_cluster_idx: target_ci,
                target_gene_id: reference_gene_id(pa, target_ci, reference_gi),
                event_class,
                identity: stats.identity,
                query_coverage: stats.query_coverage,
                reference_coverage: stats.reference_coverage,
                query_length: query_seq.len(),
                reference_length: reference_seq.len(),
                event_size_nt,
                event_size_aa,
                matches: stats.matches,
                mismatches: stats.mismatches,
                gaps_in_query: stats.gaps_in_query,
                gaps_in_reference: stats.gaps_in_reference,
                shared_kmers,
                interpretation,
            };
            let candidate_score = reference_event_score(&candidate_event);
            if best_event
                .as_ref()
                .map_or(true, |best| candidate_score > reference_event_score(best))
            {
                best_event = Some(candidate_event);
            }
        }
        if let Some(event) = best_event {
            events.push(event);
        }
    }
    events.sort_by(|a, b| {
        a.query_cluster_idx
            .cmp(&b.query_cluster_idx)
            .then(a.target_cluster_idx.cmp(&b.target_cluster_idx))
    });
    events
}

pub(super) fn reference_event_score(event: &ReferenceClusterEvent) -> f64 {
    event.identity
        + event.query_coverage
        + event.reference_coverage
        + (event.shared_kmers as f64 / 10.0)
}

pub(super) struct ReferenceEventRow<'a> {
    reference_name: &'a str,
    rescue: &'a RescueEvent,
    gene_to_cluster_maps: &'a HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
    cluster_prevalences: &'a [f64],
    cluster_seq_lengths: &'a [usize],
    matrix_policy: MatrixPolicy,
    decision: &'a str,
    accepted: bool,
    rejection_reason: &'a str,
}

pub(super) struct ReferenceAuditInputs<'a> {
    pub(super) requested_reference: &'a str,
    pub(super) rescues: &'a [RescueEvent],
    pub(super) synteny_rejected: &'a [RescueEvent],
    pub(super) paralog_rejected: &'a [RescueEvent],
    pub(super) reference_candidates: &'a [RescueEvent],
    pub(super) reference_cluster_events: &'a [ReferenceClusterEvent],
    pub(super) pa: &'a crate::parser::PresenceAbsence,
    pub(super) gene_to_cluster_maps: &'a HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
    pub(super) cluster_prevalences: &'a [f64],
    pub(super) cluster_seq_lengths: &'a [usize],
    pub(super) matrix_policy: MatrixPolicy,
    pub(super) output_dir: &'a std::path::Path,
}

pub(super) fn quality_aware_reference_publication_status(
    state: &str,
    accepted: bool,
    quality: &OrthologyQuality,
) -> &'static str {
    let raw_status = reference_publication_status(state, accepted);
    if accepted
        && raw_status == "remove_from_strict_reference_absent_set"
        && !quality_allows_matrix_correction(state, quality)
    {
        "flag_reference_like_ambiguous"
    } else {
        raw_status
    }
}

pub(super) fn quality_aware_reference_publication_confidence(
    state: &str,
    reference_status: &str,
    accepted: bool,
    quality: &OrthologyQuality,
) -> &'static str {
    let raw_confidence = reference_publication_confidence(state, reference_status);
    if accepted
        && raw_confidence == "not_reference_absent"
        && !quality_allows_matrix_correction(state, quality)
    {
        "unresolved"
    } else {
        raw_confidence
    }
}

pub(super) fn quality_aware_reference_wording(
    state: &str,
    accepted: bool,
    quality: &OrthologyQuality,
) -> &'static str {
    if accepted
        && reference_publication_status(state, accepted)
            == "remove_from_strict_reference_absent_set"
        && !quality_allows_matrix_correction(state, quality)
    {
        "Reference-like sequence is present but below adaptive quality thresholds; flag as ambiguous before removing from a strict reference-absent set."
    } else {
        reference_recommended_wording(state, accepted)
    }
}

pub(super) fn write_reference_event_row(
    w: &mut BufWriter<fs::File>,
    pa: &crate::parser::PresenceAbsence,
    row: ReferenceEventRow<'_>,
) -> Result<()> {
    let rescue = row.rescue;
    let query_cluster = pa.cluster_name(rescue.cluster_idx);
    let annotation_cluster = annotation_cluster_name(rescue, pa, row.gene_to_cluster_maps);
    let state = evidence_state(
        rescue,
        row.accepted,
        query_cluster,
        &annotation_cluster,
        row.rejection_reason,
    );
    let quality = orthology_quality(
        rescue,
        state,
        row.accepted,
        row.cluster_prevalences[rescue.cluster_idx.val()],
        row.cluster_seq_lengths[rescue.cluster_idx.val()],
    );
    let final_status = quality_aware_reference_publication_status(state, row.accepted, &quality);
    let reference_status = reference_status_for(state, row.accepted);
    let publication_confidence = quality_aware_reference_publication_confidence(
        state,
        reference_status,
        row.accepted,
        &quality,
    );
    let dichotomy = if final_status == "remove_from_strict_reference_absent_set" {
        "reference_like_or_present"
    } else {
        "reference_like_ambiguous"
    };

    write_tsv_fields(
        w,
        &[
            query_cluster.to_string(),
            row.reference_name.to_string(),
            dichotomy.to_string(),
            final_status.to_string(),
            reference_status.to_string(),
            matrix_action_with_quality(row.matrix_policy, row.accepted, state, &quality)
                .to_string(),
            state.to_string(),
            evidence_strength(state).to_string(),
            publication_confidence.to_string(),
            quality_aware_reference_wording(state, row.accepted, &quality).to_string(),
            row.decision.to_string(),
            rescue.rescue_type.to_string(),
            row.rejection_reason.to_string(),
            format!("{:.1}", rescue.identity),
            format!("{:.1}", rescue.coverage),
            rescue.contig.clone(),
            (rescue.start + 1).to_string(),
            rescue.end.to_string(),
            rescue.strand.to_string(),
            synteny_label(rescue.synteny_ok).to_string(),
            rescue.annotation_match.as_deref().unwrap_or("").to_string(),
            annotation_cluster,
            rescue
                .fragment_evidence
                .note
                .as_deref()
                .unwrap_or("")
                .to_string(),
        ],
    )?;
    Ok(())
}

pub(super) fn write_reference_shadow_candidate_row(
    w: &mut BufWriter<fs::File>,
    pa: &crate::parser::PresenceAbsence,
    row: ReferenceEventRow<'_>,
) -> Result<()> {
    let rescue = row.rescue;
    let query_cluster = pa.cluster_name(rescue.cluster_idx);
    let annotation_cluster = annotation_cluster_name(rescue, pa, row.gene_to_cluster_maps);
    let state = reference_shadow_state(rescue, query_cluster, &annotation_cluster);
    let final_status = reference_publication_status(state, false);
    let reference_status = reference_status_for(state, false);
    let publication_confidence = reference_publication_confidence(state, reference_status);

    write_tsv_fields(
        w,
        &[
            query_cluster.to_string(),
            row.reference_name.to_string(),
            "reference_like_ambiguous".to_string(),
            final_status.to_string(),
            reference_status.to_string(),
            matrix_action(row.matrix_policy, false, state).to_string(),
            state.to_string(),
            evidence_strength(state).to_string(),
            publication_confidence.to_string(),
            reference_recommended_wording(state, false).to_string(),
            "shadow_candidate".to_string(),
            rescue.rescue_type.to_string(),
            "discarded_or_subthreshold_reference_evidence".to_string(),
            format!("{:.1}", rescue.identity),
            format!("{:.1}", rescue.coverage),
            rescue.contig.clone(),
            (rescue.start + 1).to_string(),
            rescue.end.to_string(),
            rescue.strand.to_string(),
            synteny_label(rescue.synteny_ok).to_string(),
            rescue.annotation_match.as_deref().unwrap_or("").to_string(),
            annotation_cluster,
            rescue
                .fragment_evidence
                .note
                .as_deref()
                .unwrap_or("")
                .to_string(),
        ],
    )?;
    Ok(())
}

pub(super) fn write_reference_cluster_event_row(
    w: &mut BufWriter<fs::File>,
    pa: &crate::parser::PresenceAbsence,
    reference_name: &str,
    event: &ReferenceClusterEvent,
) -> Result<()> {
    let state = "EVENT_DERIVED_REFERENCE_CLUSTER";
    let reference_status = reference_status_for(state, false);
    let publication_confidence = reference_publication_confidence(state, reference_status);
    let event_aa = event
        .event_size_aa
        .map(|value| value.to_string())
        .unwrap_or_default();
    let note = format!(
        "event_class={};target_reference_cluster={};event_size_nt={};event_size_aa={};query_length={};reference_length={};query_coverage={:.1};reference_coverage={:.1};shared_kmers={};{}",
        event.event_class,
        pa.cluster_name(event.target_cluster_idx),
        event.event_size_nt,
        event_aa,
        event.query_length,
        event.reference_length,
        event.query_coverage,
        event.reference_coverage,
        event.shared_kmers,
        event.interpretation
    );

    write_tsv_fields(
        w,
        &[
            pa.cluster_name(event.query_cluster_idx).to_string(),
            reference_name.to_string(),
            "reference_like_ambiguous".to_string(),
            reference_publication_status(state, false).to_string(),
            reference_status.to_string(),
            "audit_only".to_string(),
            state.to_string(),
            evidence_strength(state).to_string(),
            publication_confidence.to_string(),
            reference_recommended_wording(state, false).to_string(),
            "reference_cluster_event".to_string(),
            event.event_class.to_string(),
            "event_derived_reference_cluster".to_string(),
            format!("{:.1}", event.identity),
            format!("{:.1}", event.query_coverage),
            String::new(),
            String::new(),
            String::new(),
            String::new(),
            String::new(),
            event.target_gene_id.clone(),
            pa.cluster_name(event.target_cluster_idx).to_string(),
            note,
        ],
    )?;
    Ok(())
}

pub(super) fn write_reference_genome_audit(inputs: ReferenceAuditInputs<'_>) -> Result<()> {
    let Some(reference_gi) = resolve_reference_genome(inputs.pa, inputs.requested_reference) else {
        warn!(
            "--reference-genome {} did not match any matrix column; skipping reference audit",
            inputs.requested_reference
        );
        return Ok(());
    };
    let reference_name = inputs.pa.genome_name(reference_gi);
    let path = inputs.output_dir.join("reference_genome_audit.tsv");
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "gene_cluster\treference_genome\treference_absence_class\tfinal_publication_status\treference_status\tmatrix_action\tevidence_state\tevidence_strength\tpublication_confidence\trecommended_wording\trescue_decision\trescue_type\trejection_reason\tidentity\tcoverage\tcontig\tstart\tend\tstrand\tsynteny\tannotation_match\tannotation_cluster\tfragment_note"
    )?;

    let rescues_by_key = build_best_rescue_map(inputs.rescues.iter());
    let synteny_rejected_by_key = build_best_rescue_map(inputs.synteny_rejected.iter());
    let paralog_rejected_by_key = build_best_rescue_map(inputs.paralog_rejected.iter());

    for ci_val in 0..inputs.pa.n_clusters() {
        let ci = ClusterIdx(ci_val);
        if inputs.pa.is_present(ci, reference_gi) {
            write_tsv_fields(
                &mut w,
                &[
                    inputs.pa.cluster_name(ci).to_string(),
                    reference_name.to_string(),
                    "reference_like_or_present".to_string(),
                    "remove_from_strict_reference_absent_set".to_string(),
                    "reference_matrix_present".to_string(),
                    "already_present".to_string(),
                    "STRICT_PRESENT".to_string(),
                    "high".to_string(),
                    "not_reference_absent".to_string(),
                    "Not a reference-absent cluster; it is already present in the audited reference matrix cell.".to_string(),
                    "matrix_present".to_string(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                    String::new(),
                ],
            )?;
            continue;
        }

        if let Some(rescue) = rescues_by_key.get(&(ci, reference_gi)).copied() {
            write_reference_event_row(
                &mut w,
                inputs.pa,
                ReferenceEventRow {
                    reference_name,
                    rescue,
                    gene_to_cluster_maps: inputs.gene_to_cluster_maps,
                    cluster_prevalences: inputs.cluster_prevalences,
                    cluster_seq_lengths: inputs.cluster_seq_lengths,
                    matrix_policy: inputs.matrix_policy,
                    decision: "accepted",
                    accepted: true,
                    rejection_reason: "",
                },
            )?;
            continue;
        }

        if let Some(rescue) = synteny_rejected_by_key.get(&(ci, reference_gi)).copied() {
            write_reference_event_row(
                &mut w,
                inputs.pa,
                ReferenceEventRow {
                    reference_name,
                    rescue,
                    gene_to_cluster_maps: inputs.gene_to_cluster_maps,
                    cluster_prevalences: inputs.cluster_prevalences,
                    cluster_seq_lengths: inputs.cluster_seq_lengths,
                    matrix_policy: inputs.matrix_policy,
                    decision: "rejected",
                    accepted: false,
                    rejection_reason: "synteny_fail",
                },
            )?;
            continue;
        }

        if let Some(rescue) = paralog_rejected_by_key.get(&(ci, reference_gi)).copied() {
            write_reference_event_row(
                &mut w,
                inputs.pa,
                ReferenceEventRow {
                    reference_name,
                    rescue,
                    gene_to_cluster_maps: inputs.gene_to_cluster_maps,
                    cluster_prevalences: inputs.cluster_prevalences,
                    cluster_seq_lengths: inputs.cluster_seq_lengths,
                    matrix_policy: inputs.matrix_policy,
                    decision: "rejected",
                    accepted: false,
                    rejection_reason: "paralog_without_synteny",
                },
            )?;
            continue;
        }

        if let Some(rescue) =
            best_reference_candidate(inputs.reference_candidates, ci, reference_gi)
        {
            write_reference_shadow_candidate_row(
                &mut w,
                inputs.pa,
                ReferenceEventRow {
                    reference_name,
                    rescue,
                    gene_to_cluster_maps: inputs.gene_to_cluster_maps,
                    cluster_prevalences: inputs.cluster_prevalences,
                    cluster_seq_lengths: inputs.cluster_seq_lengths,
                    matrix_policy: inputs.matrix_policy,
                    decision: "shadow_candidate",
                    accepted: false,
                    rejection_reason: "discarded_or_subthreshold_reference_evidence",
                },
            )?;
            continue;
        }

        if let Some(event) = inputs
            .reference_cluster_events
            .iter()
            .find(|event| event.query_cluster_idx == ci)
        {
            write_reference_cluster_event_row(&mut w, inputs.pa, reference_name, event)?;
            continue;
        }

        write_tsv_fields(
            &mut w,
            &[
                inputs.pa.cluster_name(ci).to_string(),
                reference_name.to_string(),
                "new_to_reference".to_string(),
                "keep_as_strict_reference_absent".to_string(),
                "clean_absent_from_reference".to_string(),
                "keep_absent".to_string(),
                "NO_EVIDENCE".to_string(),
                evidence_strength("NO_EVIDENCE").to_string(),
                reference_publication_confidence("NO_EVIDENCE", "clean_absent_from_reference")
                    .to_string(),
                "Can be described as absent from the audited reference genome under the current thresholds.".to_string(),
                "none".to_string(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
            ],
        )?;
    }
    w.flush()?;
    info!(
        "Reference genome audit written for {} to {}",
        reference_name,
        path.display()
    );
    Ok(())
}

pub(super) fn write_reference_cluster_events(
    events: &[ReferenceClusterEvent],
    pa: &crate::parser::PresenceAbsence,
    output_dir: &std::path::Path,
) -> Result<()> {
    let path = output_dir.join("reference_cluster_events.tsv");
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "gene_cluster\ttarget_reference_cluster\ttarget_reference_gene\treference_event_class\tidentity\tquery_coverage\treference_coverage\tquery_length\treference_length\tevent_size_nt\tevent_size_aa\tmatches\tmismatches\tgaps_in_query\tgaps_in_reference\tshared_kmers\tmatrix_recommendation\tinterpretation"
    )?;
    for event in events {
        write_tsv_fields(
            &mut w,
            &[
                pa.cluster_name(event.query_cluster_idx).to_string(),
                pa.cluster_name(event.target_cluster_idx).to_string(),
                event.target_gene_id.clone(),
                event.event_class.to_string(),
                format!("{:.1}", event.identity),
                format!("{:.1}", event.query_coverage),
                format!("{:.1}", event.reference_coverage),
                event.query_length.to_string(),
                event.reference_length.to_string(),
                event.event_size_nt.to_string(),
                event
                    .event_size_aa
                    .map(|value| value.to_string())
                    .unwrap_or_default(),
                event.matches.to_string(),
                event.mismatches.to_string(),
                event.gaps_in_query.to_string(),
                event.gaps_in_reference.to_string(),
                event.shared_kmers.to_string(),
                "keep_separate_cluster_annotate_event".to_string(),
                event.interpretation.to_string(),
            ],
        )?;
    }
    w.flush()?;
    info!(
        "Reference cluster event audit: {} rows written to {}",
        events.len(),
        path.display()
    );
    Ok(())
}
