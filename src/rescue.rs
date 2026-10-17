use crate::parser::PresenceAbsence;
use crate::search::Hit;
use crate::types::{ClusterIdx, GenomeIdx};
use ahash::{HashMap, HashMapExt};

mod conflicts;
mod genes;
mod gff;
mod orf;
mod synteny;
mod types;

pub use self::conflicts::resolve_conflicts;
pub use self::genes::build_gene_to_cluster_with_gff_aliases;
pub(crate) use self::genes::panaroo_gene_id_suffix_variants;
pub use self::gff::{GffGene, GffIndex};
pub use self::orf::check_internal_stops_in_hit;
pub(crate) use self::synteny::{validate_synteny_with_cache, SyntenyReferenceCache};
pub use self::types::{
    FragmentEvidence, PseudogeneReason, RescueEvent, RescueType, SyntenyEvidence,
};

#[cfg(test)]
use self::genes::build_gene_to_cluster;
#[cfg(test)]
use self::orf::count_internal_stops;
#[cfg(test)]
use self::orf::reverse_complement;
#[cfg(test)]
use self::synteny::{
    assess_context_support, compatible_reference_support, has_recurrent_context_support,
    validate_synteny, weighted_side_score, ReferenceSyntenyContext, SyntenyContext,
    SYNTENY_CONTEXT_GENES, SYNTENY_REFERENCE_SUPPORT_THRESHOLD,
};

// ── Fragment analysis ──

const BIOLOGICAL_FRAGMENT_MIN_GAP_BP: u64 = 50;
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
            // One segment only when the hits are the same stretch of the gene
            // AND the same stretch of the contig. Touching on the query alone
            // is not enough: the two halves of a gene broken by an insertion
            // abut exactly on the query once alignment is base-level, and sit
            // kilobases apart on the contig. Merging them on the query alone
            // turned a FRAGMENTED gene into one segment and no rescue at all.
            let query_joins = hit.query_start <= last.query_end;
            let target_joins =
                hit.target_start <= last.target_end && last.target_start <= hit.target_end;
            if query_joins && target_joins {
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

    // Union of the query intervals: separate segments may now overlap by a few
    // bases on the query, and summing their lengths would count that twice.
    let mut covered = 0u64;
    let mut reach = 0u64;
    for segment in &segments {
        let from = segment.query_start.max(reach);
        if segment.query_end > from {
            covered += segment.query_end - from;
        }
        reach = reach.max(segment.query_end);
    }
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

        if recip > 0.5 && best_gene.map_or(true, |(_, prev)| recip > prev) {
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

#[cfg(test)]
mod tests;
