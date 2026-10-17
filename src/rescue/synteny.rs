use super::genes::cell_gene_id_variants;
use super::gff::GffIndex;
use super::types::{RescueEvent, SyntenyEvidence};
use crate::parser::PresenceAbsence;
use crate::types::{ClusterIdx, GenomeIdx};
use ahash::HashMap;
use std::collections::BTreeSet;

#[cfg(test)]
pub(super) const SYNTENY_CONTEXT_GENES: usize = 2;
const SYNTENY_MAX_CONTEXT_GENES: usize = 5;
pub(super) const SYNTENY_REFERENCE_SUPPORT_THRESHOLD: f64 = 0.50;

#[derive(Debug, Clone, Default)]
pub(super) struct SyntenyContext {
    pub(super) upstream: Vec<ClusterIdx>,
    pub(super) downstream: Vec<ClusterIdx>,
}

impl SyntenyContext {
    fn has_any_cluster(&self) -> bool {
        !self.upstream.is_empty() || !self.downstream.is_empty()
    }

    fn side_matches(&self, other: &SyntenyContext) -> bool {
        let (_, up_matches) = weighted_side_score(&self.upstream, &other.upstream);
        let (_, down_matches) = weighted_side_score(&self.downstream, &other.downstream);
        let up_comparable = !self.upstream.is_empty() && !other.upstream.is_empty();
        let down_comparable = !self.downstream.is_empty() && !other.downstream.is_empty();

        if up_comparable && down_comparable {
            return up_matches > 0 && down_matches > 0;
        }

        up_matches + down_matches > 0
    }
}

#[derive(Debug, Clone, Default)]
pub(super) struct SyntenyAssessment {
    pub(super) passes: bool,
    pub(super) score: f64,
    pub(super) support: usize,
    pub(super) reference_contexts: usize,
    pub(super) n_neighbors: usize,
    pub(super) matched_upstream: usize,
    pub(super) matched_downstream: usize,
}

#[derive(Debug, Clone)]
pub(super) struct ReferenceSyntenyContext {
    pub(super) genome_idx: GenomeIdx,
    pub(super) context: SyntenyContext,
}

#[derive(Debug, Default)]
pub(crate) struct SyntenyReferenceCache {
    contexts: HashMap<(ClusterIdx, usize), Vec<ReferenceSyntenyContext>>,
    /// Whether any genome carrying the cluster has one of its genes in a GFF.
    located: HashMap<ClusterIdx, bool>,
}

/// True when at least one genome carrying the cluster has a GFF in which one of
/// the cluster's genes can be found, so its neighbourhood could be read.
fn cluster_reference_located(
    cache: &mut SyntenyReferenceCache,
    ci: ClusterIdx,
    pa: &PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
) -> bool {
    *cache.located.entry(ci).or_insert_with(|| {
        pa.present_genomes(ci).iter().any(|&ref_gi| {
            let Some(ref_gff) = gff_indices.get(&ref_gi) else {
                return false;
            };
            pa.cell(ci, ref_gi)
                .split(';')
                .map(str::trim)
                .filter(|id| !id.is_empty())
                .any(|id| {
                    ref_gff.get_by_id(id).is_some()
                        || cell_gene_id_variants(id)
                            .into_iter()
                            .skip(1)
                            .any(|variant| ref_gff.get_by_id(variant).is_some())
                })
        })
    })
}

pub(super) fn weighted_side_score(target: &[ClusterIdx], reference: &[ClusterIdx]) -> (f64, usize) {
    if target.is_empty() || reference.is_empty() {
        return (0.0, 0);
    }

    let mut matched_weight = 0.0;
    let mut denominator = 0.0;
    let mut matches = 0usize;
    let mut used_reference = vec![false; reference.len()];
    for (target_pos, target_cluster) in target.iter().enumerate() {
        let target_weight = 1.0 / (target_pos + 1) as f64;
        denominator += target_weight;
        let best_reference = reference
            .iter()
            .enumerate()
            .filter(|(reference_pos, reference_cluster)| {
                !used_reference[*reference_pos] && *reference_cluster == target_cluster
            })
            .map(|(reference_pos, _)| (reference_pos, 1.0 / (reference_pos + 1) as f64))
            .max_by(|(_, a), (_, b)| a.total_cmp(b));
        if let Some((reference_pos, best_reference_weight)) = best_reference {
            used_reference[reference_pos] = true;
            matches += 1;
            matched_weight += target_weight * best_reference_weight;
        }
    }

    if denominator == 0.0 {
        (0.0, matches)
    } else {
        (matched_weight / denominator, matches)
    }
}

fn context_score(target: &SyntenyContext, reference: &SyntenyContext) -> (f64, usize, usize) {
    let (up_score, up_matches) = weighted_side_score(&target.upstream, &reference.upstream);
    let (down_score, down_matches) = weighted_side_score(&target.downstream, &reference.downstream);

    let comparable_sides =
        usize::from(!target.upstream.is_empty()) + usize::from(!target.downstream.is_empty());
    if comparable_sides == 0 {
        return (0.0, up_matches, down_matches);
    }

    let score = (if target.upstream.is_empty() {
        0.0
    } else {
        up_score
    } + if target.downstream.is_empty() {
        0.0
    } else {
        down_score
    }) / comparable_sides as f64;
    (score, up_matches, down_matches)
}

fn synteny_context(
    gff: &GffIndex,
    gene_to_cluster: &HashMap<String, ClusterIdx>,
    contig: &str,
    start: u64,
    end: u64,
    center_gene_id: Option<&str>,
    n_neighbors: usize,
) -> SyntenyContext {
    let mut context = SyntenyContext::default();
    let (upstream, downstream) =
        gff.context_window(contig, start, end, center_gene_id, n_neighbors);
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

#[cfg(test)]
pub(super) fn compatible_reference_support(
    target_context: &SyntenyContext,
    reference_contexts: &[SyntenyContext],
) -> usize {
    let reference_contexts = labeled_test_contexts(reference_contexts);
    assess_context_support(target_context, &reference_contexts, SYNTENY_CONTEXT_GENES).support
}

#[cfg(test)]
pub(super) fn has_recurrent_context_support(
    target_context: &SyntenyContext,
    reference_contexts: &[SyntenyContext],
) -> bool {
    let reference_contexts = labeled_test_contexts(reference_contexts);
    assess_context_support(target_context, &reference_contexts, SYNTENY_CONTEXT_GENES).passes
}

#[cfg(test)]
fn labeled_test_contexts(reference_contexts: &[SyntenyContext]) -> Vec<ReferenceSyntenyContext> {
    reference_contexts
        .iter()
        .enumerate()
        .map(|(idx, context)| ReferenceSyntenyContext {
            genome_idx: GenomeIdx(idx),
            context: context.clone(),
        })
        .collect()
}

pub(super) fn assess_context_support(
    target_context: &SyntenyContext,
    reference_contexts: &[ReferenceSyntenyContext],
    n_neighbors: usize,
) -> SyntenyAssessment {
    let mut support = 0usize;
    let mut support_score_sum = 0.0;
    let mut best_score = 0.0;
    let mut best_up = 0usize;
    let mut best_down = 0usize;
    let reference_genomes: BTreeSet<GenomeIdx> = reference_contexts
        .iter()
        .map(|reference_context| reference_context.genome_idx)
        .collect();
    let mut support_genomes = BTreeSet::new();

    for reference_context in reference_contexts {
        let (score, up_matches, down_matches) =
            context_score(target_context, &reference_context.context);
        if score > best_score {
            best_score = score;
            best_up = up_matches;
            best_down = down_matches;
        }
        if score >= SYNTENY_REFERENCE_SUPPORT_THRESHOLD
            && target_context.side_matches(&reference_context.context)
        {
            support += 1;
            support_score_sum += score;
            support_genomes.insert(reference_context.genome_idx);
        }
    }

    let score = if support > 0 {
        support_score_sum / support as f64
    } else {
        best_score
    };
    let passes = if reference_genomes.len() <= 3 {
        support > 0
    } else {
        support_genomes.len() >= 2
    };

    SyntenyAssessment {
        passes,
        score,
        support,
        reference_contexts: reference_contexts.len(),
        n_neighbors,
        matched_upstream: best_up,
        matched_downstream: best_down,
    }
}

fn better_synteny_assessment(
    current: Option<SyntenyAssessment>,
    candidate: SyntenyAssessment,
) -> Option<SyntenyAssessment> {
    let Some(current) = current else {
        return Some(candidate);
    };
    if candidate.passes != current.passes {
        return Some(if candidate.passes { candidate } else { current });
    }
    if candidate.score > current.score + f64::EPSILON {
        return Some(candidate);
    }
    if (candidate.score - current.score).abs() <= f64::EPSILON
        && candidate.n_neighbors < current.n_neighbors
    {
        return Some(candidate);
    }
    Some(current)
}

fn reference_contexts_for_cluster<'a>(
    cache: &'a mut SyntenyReferenceCache,
    ci: ClusterIdx,
    n_neighbors: usize,
    pa: &PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    gene_to_cluster_maps: &HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
    verbose: bool,
) -> &'a [ReferenceSyntenyContext] {
    cache
        .contexts
        .entry((ci, n_neighbors))
        .or_insert_with(|| {
            let mut reference_contexts = Vec::new();
            let mut ref_genomes = pa.present_genomes(ci).to_vec();
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
                for ref_gene_id in ref_cell
                    .split(';')
                    .map(str::trim)
                    .filter(|id| !id.is_empty())
                {
                    let ref_gene = ref_gff.get_by_id(ref_gene_id).or_else(|| {
                        cell_gene_id_variants(ref_gene_id)
                            .into_iter()
                            .skip(1)
                            .find_map(|id| ref_gff.get_by_id(id))
                    });
                    let Some(ref_gene) = ref_gene else {
                        if verbose && n_neighbors == 1 {
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
                        n_neighbors,
                    );

                    if ref_context.has_any_cluster() {
                        reference_contexts.push(ReferenceSyntenyContext {
                            genome_idx: ref_gi,
                            context: ref_context,
                        });
                    }
                }
            }
            reference_contexts
        })
        .as_slice()
}

#[cfg(test)]
pub(super) fn validate_synteny(
    rescue: &mut RescueEvent,
    pa: &PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    gene_to_cluster_maps: &HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
    verbose: bool,
) {
    let mut cache = SyntenyReferenceCache::default();
    validate_synteny_with_cache(
        rescue,
        pa,
        gff_indices,
        gene_to_cluster_maps,
        &mut cache,
        verbose,
    );
}

pub(crate) fn validate_synteny_with_cache(
    rescue: &mut RescueEvent,
    pa: &PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    gene_to_cluster_maps: &HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
    cache: &mut SyntenyReferenceCache,
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

    let mut best_assessment: Option<SyntenyAssessment> = None;
    let mut saw_target_context = false;
    let mut saw_reference_context = false;

    for n_neighbors in 1..=SYNTENY_MAX_CONTEXT_GENES {
        let target_context = synteny_context(
            target_gff,
            target_g2c,
            &rescue.contig,
            rescue.start,
            rescue.end,
            rescue.annotation_match.as_deref(),
            n_neighbors,
        );

        if !target_context.has_any_cluster() {
            continue;
        }
        saw_target_context = true;

        let reference_contexts = reference_contexts_for_cluster(
            cache,
            ci,
            n_neighbors,
            pa,
            gff_indices,
            gene_to_cluster_maps,
            verbose,
        );
        if reference_contexts.is_empty() {
            continue;
        }
        saw_reference_context = true;
        let assessment = assess_context_support(&target_context, reference_contexts, n_neighbors);
        best_assessment = better_synteny_assessment(best_assessment, assessment);
    }

    let assessment = if let Some(assessment) = best_assessment {
        assessment
    } else if saw_target_context
        && !saw_reference_context
        && cluster_reference_located(cache, ci, pa, gff_indices)
    {
        // The reference copies were found, but none sits next to a clustered
        // gene: their neighbourhood shares nothing with the target's.
        SyntenyAssessment {
            passes: false,
            score: 0.0,
            support: 0,
            reference_contexts: 0,
            n_neighbors: SYNTENY_MAX_CONTEXT_GENES,
            matched_upstream: 0,
            matched_downstream: 0,
        }
    } else if saw_target_context && !saw_reference_context {
        // No genome carrying the cluster has it in a GFF, so there is nothing
        // to compare the target's context with. That is not an incompatible
        // context: it used to be a FAIL with score 0.000 whenever the other
        // genomes simply had no annotation.
        if verbose {
            log::debug!(
                "Synteny not assessed: no reference context for cluster={} genome={}",
                ci.val(),
                gi.val()
            );
        }
        return;
    } else {
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
    };

    rescue.synteny_ok = Some(assessment.passes);
    rescue.synteny_evidence = SyntenyEvidence {
        score: Some(assessment.score),
        support: assessment.support,
        reference_contexts: assessment.reference_contexts,
        n_neighbors: assessment.n_neighbors,
        matched_upstream: assessment.matched_upstream,
        matched_downstream: assessment.matched_downstream,
    };

    if assessment.passes {
        rescue.confidence *= 1.2;
    } else {
        rescue.confidence *= 0.6;
    }
}
