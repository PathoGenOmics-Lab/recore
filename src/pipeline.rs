use ahash::{HashMap, HashMapExt};
use anyhow::{Context, Result};
use log::{debug, info, warn};
use std::collections::BTreeMap;
use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::cli::{Args, MatrixPolicy};
use crate::nscluster;
use crate::parser::{
    extract_as_id, find_annotation_file, find_genome_files, load_genome_contig_names,
    normalize_gff_records_to_fasta, parse_gff_genes, parse_pan_reference, parse_presence_absence,
    resolve_genome_files_for_pa, validate_inputs,
};
use crate::relationships::write_cluster_relationships;
use crate::rescue::{
    build_gene_to_cluster_with_gff_aliases, check_internal_stops_in_hit, classify_hits,
    panaroo_gene_id_suffix_variants, resolve_conflicts, validate_synteny_with_cache, GffGene,
    GffIndex, PseudogeneReason, RescueEvent, RescueType, SyntenyEvidence, SyntenyReferenceCache,
};
use crate::search::{batch_search, Hit};
use crate::types::{ClusterIdx, GenomeIdx};

mod evolution;
mod output;
mod reference_audit;

use self::evolution::{
    load_phylo_tree, nsc_presence_updates, search_mobile_element_reference,
    write_cluster_locus_contexts, write_evolutionary_events, EvolutionaryEventInputs,
    MobileElementSearchInput,
};
use self::output::{
    count_core_after, count_soft_core_after, write_corrected_csv, write_corrected_gffs,
    write_multistate_matrix, write_nsc_report, write_run_manifest, write_summary,
    MultistateMatrixInputs, SummaryStats,
};
use self::reference_audit::{
    build_reference_cluster_events, reference_hit_shadow_candidates, resolve_reference_genome,
    write_reference_cluster_events, write_reference_genome_audit, ReferenceAuditInputs,
};

#[cfg(test)]
use self::evolution::{
    nearest_mobile_element_for_rescue, nsc_role_for_cluster, resolve_phylo_tips,
    summarize_phylogenetic_pattern, MobileElementHit, NewickParser, PhyloTree,
};
#[cfg(test)]
use self::output::copy_text_maybe_gzip_with_header;
#[cfg(test)]
use self::reference_audit::{
    best_reference_candidate, display_reference_gene_id, reference_publication_confidence,
    reference_publication_status, reference_shadow_state, reference_status_for,
};

fn rescue_type_priority(rescue_type: &RescueType) -> u8 {
    match rescue_type {
        RescueType::Split { .. } => 7,
        RescueType::MissingCds => 6,
        RescueType::StartVariant { .. } => 5,
        RescueType::Dropped => 4,
        RescueType::AlignmentSplit { .. } => 3,
        RescueType::Fragmented { .. } => 2,
        RescueType::Pseudogene { .. } => 1,
        RescueType::Paralog { .. } => 0,
    }
}

fn rescue_selection_score(rescue: &RescueEvent) -> f64 {
    let synteny_bonus = match rescue.synteny_ok {
        Some(true) => 2.0,
        Some(false) => -2.0,
        None => 0.0,
    };
    rescue.confidence + rescue.identity / 100.0 + rescue.coverage / 100.0 + synteny_bonus
}

fn rescue_preferred_over(candidate: &RescueEvent, current: &RescueEvent) -> bool {
    let score_order = rescue_selection_score(candidate).total_cmp(&rescue_selection_score(current));
    if score_order != std::cmp::Ordering::Equal {
        return score_order == std::cmp::Ordering::Greater;
    }

    rescue_type_priority(&candidate.rescue_type)
        .cmp(&rescue_type_priority(&current.rescue_type))
        .then_with(|| candidate.contig.cmp(&current.contig).reverse())
        .then_with(|| candidate.start.cmp(&current.start).reverse())
        .then_with(|| candidate.end.cmp(&current.end).reverse())
        == std::cmp::Ordering::Greater
}

fn build_best_rescue_map<'a, I>(rescues: I) -> HashMap<(ClusterIdx, GenomeIdx), &'a RescueEvent>
where
    I: IntoIterator<Item = &'a RescueEvent>,
{
    let mut map: HashMap<(ClusterIdx, GenomeIdx), &RescueEvent> = HashMap::new();
    for rescue in rescues {
        let key = (rescue.cluster_idx, rescue.genome_idx);
        match map.get(&key).copied() {
            Some(current) if rescue_preferred_over(rescue, current) => {
                map.insert(key, rescue);
            }
            None => {
                map.insert(key, rescue);
            }
            _ => {}
        }
    }
    map
}

pub fn run(args: &Args) -> Result<()> {
    validate_args(args)?;

    // ── 1. Parse inputs ──
    let pa = parse_presence_absence(&args.presence_absence)?;
    let pan_ref = parse_pan_reference(&args.pan_reference, args.multi_allele)?;
    let raw_genome_files = find_genome_files(&args.genomes)?;
    let genome_files = resolve_genome_files_for_pa(&pa, &raw_genome_files);
    let alias_resolved = genome_files
        .iter()
        .filter(|(name, _)| !raw_genome_files.contains_key(name.as_str()))
        .count();
    if alias_resolved > 0 {
        info!(
            "Resolved {} genome FASTA files by embedded AS identifier",
            alias_resolved
        );
    }
    // ── 1b. Cross-validate inputs (includes path traversal checks) ──
    validate_inputs(&pa, &pan_ref, &genome_files)?;

    let phylo_tree = if let Some(tree_path) = args.tree.as_deref() {
        Some(load_phylo_tree(tree_path, &pa)?)
    } else {
        None
    };

    // ── 2. Filter gaps by prevalence (streaming, no intermediate Vec) ──
    let min_prevalence = args.min_prevalence;
    let mut gaps_by_genome: HashMap<GenomeIdx, Vec<ClusterIdx>> = HashMap::new();
    let mut total_gaps = 0usize;
    let mut filtered_gaps = 0usize;
    let mut skipped_by_prevalence = 0usize;

    let n_genomes = pa.n_genomes();
    let n_clusters = pa.n_clusters();

    for ci_val in 0..n_clusters {
        let ci = ClusterIdx(ci_val);
        let prevalence = pa.present_count(ci) as f64 / n_genomes as f64;
        for gi_val in 0..n_genomes {
            let gi = GenomeIdx(gi_val);
            if !pa.is_present(ci, gi) {
                total_gaps += 1;
                if prevalence < min_prevalence {
                    skipped_by_prevalence += 1;
                } else {
                    gaps_by_genome.entry(gi).or_default().push(ci);
                    filtered_gaps += 1;
                }
            }
        }
    }

    let core_before = pa.core_size();
    info!(
        "Gaps: {} total, {} after prevalence filter (≥{:.0}%), {} skipped",
        total_gaps,
        filtered_gaps,
        min_prevalence * 100.0,
        skipped_by_prevalence
    );
    info!(
        "Core genome (before): {} / {} clusters",
        core_before, n_clusters
    );

    // ── 3. Filter to genomes with FASTA files ──
    let mut searchable_genomes: Vec<(GenomeIdx, Vec<ClusterIdx>)> = Vec::new();
    let mut missing_genomes = 0usize;
    let mut missing_fasta_gaps = 0usize;
    let mut missing_names: Vec<String> = Vec::new();
    for (gi, cluster_indices) in gaps_by_genome {
        let genome_name = pa.genome_name(gi);
        if genome_files.contains_key(genome_name) {
            searchable_genomes.push((gi, cluster_indices));
        } else {
            missing_genomes += 1;
            missing_fasta_gaps += cluster_indices.len();
            if args.verbose {
                missing_names.push(genome_name.to_string());
            }
        }
    }
    if missing_genomes > 0 {
        warn!(
            "{} genomes with gaps have no FASTA file (skipped)",
            missing_genomes
        );
        for name in &missing_names {
            debug!("  Missing FASTA: {}", name);
        }
    }
    let searchable_gaps: usize = searchable_genomes
        .iter()
        .map(|(_, clusters)| clusters.len())
        .sum();

    // ── 4. Load GFF annotations ──
    let mut gff_paths: HashMap<GenomeIdx, PathBuf> = HashMap::new();
    let mut gff_contig_remaps: Vec<(String, String, String, usize)> = Vec::new();
    let gff_indices: HashMap<GenomeIdx, GffIndex> = if let Some(ref ann_dir) = args.annotations {
        let mut data = HashMap::new();
        for gi_val in 0..n_genomes {
            let gi = GenomeIdx(gi_val);
            let genome_name = pa.genome_name(gi);
            let Some(path) = find_annotation_file(ann_dir, genome_name) else {
                continue;
            };
            match parse_gff_genes(&path) {
                Ok(mut raw_genes) => {
                    if let Some(fasta_path) = genome_files.get(genome_name) {
                        match load_genome_contig_names(fasta_path) {
                            Ok(contig_names) => {
                                if let Some((from, to, n_records)) =
                                    normalize_gff_records_to_fasta(&mut raw_genes, &contig_names)
                                {
                                    gff_contig_remaps.push((
                                        genome_name.to_string(),
                                        from,
                                        to,
                                        n_records,
                                    ));
                                }
                            }
                            Err(e) => {
                                warn!("Cannot read FASTA contig names for {}: {}", genome_name, e)
                            }
                        }
                    }
                    let genes: Vec<(GffGene, Vec<String>, Vec<String>)> = raw_genes
                        .into_iter()
                        .map(|(contig, start, end, strand, gene_id, aliases, labels)| {
                            (
                                GffGene {
                                    contig,
                                    start,
                                    end,
                                    strand,
                                    gene_id,
                                },
                                aliases,
                                labels,
                            )
                        })
                        .collect();
                    data.insert(gi, GffIndex::build_with_aliases_and_labels(genes));
                    gff_paths.insert(gi, path.clone());
                }
                Err(e) => warn!("Cannot parse GFF for {}: {}", genome_name, e),
            }
        }
        let loaded = data.len();
        let missing = n_genomes - loaded;
        if loaded > 0 {
            info!("Loaded GFF annotations for {} genomes", loaded);
        }
        if !gff_contig_remaps.is_empty() {
            let examples: Vec<String> = gff_contig_remaps
                .iter()
                .take(3)
                .map(|(genome, from, to, _)| format!("{}:{}→{}", genome, from, to))
                .collect();
            info!(
                "Normalized single-contig GFF seqids for {} genomes ({})",
                gff_contig_remaps.len(),
                examples.join(", ")
            );
        }
        if missing > 0 && loaded > 0 {
            warn!(
                "{} genomes have no matching GFF file in {}",
                missing,
                ann_dir.display()
            );
            if args.verbose {
                let mut shown = 0;
                for gi_val in 0..n_genomes {
                    let gi = GenomeIdx(gi_val);
                    if !data.contains_key(&gi) {
                        debug!("  Missing GFF: {}", pa.genome_name(gi));
                        shown += 1;
                        if shown >= 10 {
                            debug!("  ... and {} more", missing - shown);
                            break;
                        }
                    }
                }
            }
        } else if loaded == 0 {
            warn!("No GFF files found in {} — synteny and annotation-based classification will be skipped", ann_dir.display());
            warn!(
                "Expected files like: {}.gff or {}.gff3",
                pa.genome_name(GenomeIdx(0)),
                pa.genome_name(GenomeIdx(0))
            );
        }
        data
    } else {
        HashMap::new()
    };

    // ── 5. Build gene_to_cluster maps ONCE ──
    let gene_to_cluster_maps: HashMap<GenomeIdx, HashMap<String, ClusterIdx>> = {
        let mut maps = HashMap::new();
        for (gi, _) in &searchable_genomes {
            maps.entry(*gi).or_insert_with(|| {
                build_gene_to_cluster_with_gff_aliases(&pa, *gi, gff_indices.get(gi))
            });
        }
        if args.synteny || !gff_indices.is_empty() {
            for &gi in gff_indices.keys() {
                maps.entry(gi).or_insert_with(|| {
                    build_gene_to_cluster_with_gff_aliases(&pa, gi, gff_indices.get(&gi))
                });
            }
        }
        maps
    };

    // ── 6. Batch search: ONE minimap2 invocation ──
    let hits_map = batch_search(
        &pan_ref,
        &genome_files,
        &searchable_genomes,
        pa.genome_names(),
        pa.cluster_names(),
        args.threads,
        args.min_identity,
        &args.preset,
        false,
        50,
    )?;

    let reference_hit_candidates = if let Some(reference_genome) = args.reference_genome.as_deref()
    {
        if let Some(reference_gi) = resolve_reference_genome(&pa, reference_genome) {
            let reference_name = pa.genome_name(reference_gi);
            if genome_files.contains_key(reference_name) {
                let reference_gaps: Vec<ClusterIdx> = (0..n_clusters)
                    .map(ClusterIdx)
                    .filter(|&ci| !pa.is_present(ci, reference_gi))
                    .collect();
                if reference_gaps.is_empty() {
                    Vec::new()
                } else {
                    info!(
                        "Running reference shadow audit for {} absent clusters in {}",
                        reference_gaps.len(),
                        reference_name
                    );
                    let reference_search = vec![(reference_gi, reference_gaps)];
                    let reference_hits = batch_search(
                        &pan_ref,
                        &genome_files,
                        &reference_search,
                        pa.genome_names(),
                        pa.cluster_names(),
                        args.threads,
                        args.min_identity,
                        &args.preset,
                        true,
                        5000,
                    )?;
                    reference_hit_shadow_candidates(
                        &reference_hits,
                        &pa,
                        &gff_indices,
                        &gene_to_cluster_maps,
                    )
                }
            } else {
                warn!(
                    "Reference shadow audit skipped for {}: no resolved FASTA file",
                    reference_name
                );
                Vec::new()
            }
        } else {
            Vec::new()
        }
    } else {
        Vec::new()
    };

    let reference_cluster_events = if let Some(reference_genome) = args.reference_genome.as_deref()
    {
        if let Some(reference_gi) = resolve_reference_genome(&pa, reference_genome) {
            build_reference_cluster_events(&pa, &pan_ref, reference_gi)
        } else {
            Vec::new()
        }
    } else {
        Vec::new()
    };

    // ── 7. Classify hits (parallel) ──
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .context("Failed to build rayon thread pool")?;

    let classify_inputs: Vec<(ClusterIdx, GenomeIdx, &[crate::search::Hit])> = hits_map
        .iter()
        .map(|((gi, ci), hits)| (*ci, *gi, hits.as_slice()))
        .collect();

    let empty_map: HashMap<String, ClusterIdx> = HashMap::new();

    let rescues: Vec<RescueEvent> = pool.install(|| {
        use rayon::prelude::*;
        classify_inputs
            .par_iter()
            .filter_map(|&(ci, gi, hits)| {
                let gff_index = gff_indices.get(&gi);
                let g2c = gene_to_cluster_maps.get(&gi).unwrap_or(&empty_map);
                classify_hits(hits, ci, gi, &pa, gff_index, g2c, args.min_coverage)
            })
            .collect()
    });

    let mut rescues = rescues;
    info!("Raw rescues before filtering: {}", rescues.len());

    // ── 7a. Exact missing-CDS pass ──
    //
    // Repetitive families can hide true CDS omissions when the best annotated
    // paralog consumes minimap2's global secondary-hit budget. Run an exact,
    // high-secondary pass, but only import MISSING_CDS candidates; all other
    // extra alignments remain ignored by the binary correction path.
    if args.synteny && !gff_indices.is_empty() && !searchable_genomes.is_empty() {
        let missing_cds_chunk_size = 64usize;
        info!(
            "Running exact missing-CDS pass (-c) in {}-genome chunks...",
            missing_cds_chunk_size
        );
        let before = rescues.len();
        let mut seen_missing: ahash::HashSet<(ClusterIdx, GenomeIdx, String, u64, u64)> = rescues
            .iter()
            .map(|r| {
                (
                    r.cluster_idx,
                    r.genome_idx,
                    r.contig.clone(),
                    r.start,
                    r.end,
                )
            })
            .collect();

        for chunk in searchable_genomes.chunks(missing_cds_chunk_size) {
            let max_secondary_missing_cds = chunk.len().saturating_mul(24).clamp(50, 5000);
            let missing_hits_map = batch_search(
                &pan_ref,
                &genome_files,
                chunk,
                pa.genome_names(),
                pa.cluster_names(),
                args.threads,
                args.min_identity,
                &args.preset,
                true,
                max_secondary_missing_cds,
            )?;
            let missing_inputs: Vec<(ClusterIdx, GenomeIdx, &[crate::search::Hit])> =
                missing_hits_map
                    .iter()
                    .map(|((gi, ci), hits)| (*ci, *gi, hits.as_slice()))
                    .collect();
            let missing_rescues: Vec<RescueEvent> = pool.install(|| {
                use rayon::prelude::*;
                missing_inputs
                    .par_iter()
                    .filter_map(|&(ci, gi, hits)| {
                        let gff_index = gff_indices.get(&gi);
                        let g2c = gene_to_cluster_maps.get(&gi).unwrap_or(&empty_map);
                        classify_hits(hits, ci, gi, &pa, gff_index, g2c, args.min_coverage)
                    })
                    .filter(|rescue| matches!(rescue.rescue_type, RescueType::MissingCds))
                    .collect()
            });

            for rescue in missing_rescues {
                let key = (
                    rescue.cluster_idx,
                    rescue.genome_idx,
                    rescue.contig.clone(),
                    rescue.start,
                    rescue.end,
                );
                if seen_missing.insert(key) {
                    rescues.push(rescue);
                }
            }
        }
        info!(
            "Exact missing-CDS pass added {} candidates",
            rescues.len() - before
        );
    }

    // ── 7b. Internal stop codon scanning for sequence-only rescues ──
    // Load genome sequences on-demand for rescued genomes only
    {
        let mut rescues_by_genome: HashMap<GenomeIdx, Vec<usize>> = HashMap::default();
        for (idx, rescue) in rescues.iter().enumerate() {
            if matches!(
                rescue.rescue_type,
                RescueType::Dropped | RescueType::MissingCds
            ) {
                rescues_by_genome
                    .entry(rescue.genome_idx)
                    .or_default()
                    .push(idx);
            }
        }

        if !rescues_by_genome.is_empty() {
            info!(
                "Scanning {} genomes for internal stop codons in sequence-only rescues...",
                rescues_by_genome.len()
            );
            let mut upgraded = 0usize;

            for (gi, rescue_indices) in rescues_by_genome {
                let genome_name = pa.genome_name(gi);
                let Some(fasta_path) = genome_files.get(genome_name) else {
                    continue;
                };
                let Ok(contigs) = crate::parser::load_genome_contigs(fasta_path) else {
                    continue;
                };

                for rescue_idx in rescue_indices {
                    let rescue = &mut rescues[rescue_idx];
                    if let Some(n_stops) = check_internal_stops_in_hit(
                        // Build a minimal hit from rescue coords for the check
                        &crate::search::Hit {
                            query: String::new(),
                            cluster: String::new(),
                            query_length: 0,
                            query_start: 0,
                            query_end: 0,
                            strand: rescue.strand,
                            target: rescue.contig.clone(),
                            target_start: rescue.start,
                            target_end: rescue.end,
                            matches: 0,
                            alignment_length: 0,
                            mapq: 0,
                        },
                        Some(&contigs),
                    ) {
                        rescue.rescue_type = RescueType::Pseudogene {
                            reason: PseudogeneReason::InternalStop { n_stops },
                        };
                        rescue.confidence *= 0.5; // Pseudogenes get lower confidence
                        upgraded += 1;
                    }
                }
            }

            if upgraded > 0 {
                info!(
                    "Upgraded {} sequence-only rescues to PSEUDOGENE:internal_stop",
                    upgraded
                );
            }
        }
    }

    // ── 8. Synteny validation ──
    if args.synteny && !gff_indices.is_empty() {
        info!("Validating synteny...");
        let mut synteny_cache = SyntenyReferenceCache::default();
        for rescue in &mut rescues {
            validate_synteny_with_cache(
                rescue,
                &pa,
                &gff_indices,
                &gene_to_cluster_maps,
                &mut synteny_cache,
                args.verbose,
            );
        }
        let validated = rescues
            .iter()
            .filter(|r| r.synteny_ok == Some(true))
            .count();
        let failed = rescues
            .iter()
            .filter(|r| r.synteny_ok == Some(false))
            .count();
        let skipped = rescues.iter().filter(|r| r.synteny_ok.is_none()).count();
        info!(
            "Synteny: {} validated, {} failed, {} skipped (no flanking genes found)",
            validated, failed, skipped
        );
        if validated == 0 && failed == 0 && skipped > 0 {
            warn!("Synteny could not validate any rescue. Possible causes:");
            warn!("  - GFF gene IDs don't match the IDs in Panaroo's CSV cells");
            warn!("  - GFF has no genes flanking the rescued regions");
            warn!("  Rerun with -v to see detailed synteny debug output");
        }
    } else if args.synteny && gff_indices.is_empty() {
        warn!("Synteny requested (--synteny) but no GFF annotations loaded — skipping");
    }

    let reference_audit_candidates = if args.reference_genome.is_some() {
        let mut candidates = reference_hit_candidates;
        candidates.extend(rescues.clone());
        candidates
    } else {
        Vec::new()
    };

    // Explicit synteny failures are contradictory orthology evidence. Keep
    // them in an audit report, but do not use them to correct the matrix.
    let mut synteny_rejected = Vec::new();
    if args.synteny {
        let before_rejection = rescues.len();
        let old_rescues = std::mem::take(&mut rescues);
        let (rejected, kept): (Vec<_>, Vec<_>) = old_rescues
            .into_iter()
            .partition(|r| r.synteny_ok == Some(false));
        rescues = kept;
        synteny_rejected = rejected;
        if !synteny_rejected.is_empty() {
            info!(
                "Synteny rejection: {} → {} rescues ({} explicit FAIL candidates written to synteny_rejected_report.tsv)",
                before_rejection,
                rescues.len(),
                synteny_rejected.len()
            );
        }
    }

    // Paralog-only sequence hits are useful audit evidence, but they are not
    // orthology evidence unless the local gene context explicitly supports
    // the same locus.
    let before_paralog_rejection = rescues.len();
    let old_rescues = std::mem::take(&mut rescues);
    let (paralog_rejected, kept): (Vec<_>, Vec<_>) = old_rescues
        .into_iter()
        .partition(is_paralog_without_synteny_support);
    rescues = kept;
    if !paralog_rejected.is_empty() {
        info!(
            "Paralog rejection: {} → {} rescues ({} paralog-only candidates written to paralog_rejected_report.tsv)",
            before_paralog_rejection,
            rescues.len(),
            paralog_rejected.len()
        );
    }

    // ── 9. Conflict resolution ──
    let before_conflicts = rescues.len();
    resolve_conflicts(&mut rescues);
    if before_conflicts != rescues.len() {
        info!(
            "Conflict resolution: {} → {} rescues ({} conflicts resolved)",
            before_conflicts,
            rescues.len(),
            before_conflicts - rescues.len()
        );
    }

    // ── 10. Count rescue types ──
    let mut dropped = 0usize;
    let mut missing_cds = 0usize;
    let mut split = 0usize;
    let mut fragmented = 0usize;
    let mut alignment_split = 0usize;
    let mut start_variant = 0usize;
    let mut paralog = 0usize;
    let mut pseudogene = 0usize;
    for r in &rescues {
        match &r.rescue_type {
            RescueType::Dropped => dropped += 1,
            RescueType::MissingCds => missing_cds += 1,
            RescueType::Split { .. } => split += 1,
            RescueType::Fragmented { .. } => fragmented += 1,
            RescueType::AlignmentSplit { .. } => alignment_split += 1,
            RescueType::StartVariant { .. } => start_variant += 1,
            RescueType::Paralog { .. } => paralog += 1,
            RescueType::Pseudogene { .. } => pseudogene += 1,
        }
    }
    info!(
        "Final rescues: {} (DROPPED: {}, MISSING_CDS: {}, SPLIT: {}, FRAGMENTED: {}, ALIGNMENT_SPLIT: {}, START_VARIANT: {}, PARALOG: {}, PSEUDOGENE: {})",
        rescues.len(), dropped, missing_cds, split, fragmented, alignment_split, start_variant, paralog, pseudogene
    );

    if args.verbose && !rescues.is_empty() {
        debug!("Top rescues (first 10):");
        for r in rescues.iter().take(10) {
            let cluster_name = pa.cluster_name(r.cluster_idx);
            let genome_name = pa.genome_name(r.genome_idx);
            debug!(
                "  {} in {} → {:?} (id={:.1}%, cov={:.1}%, conf={:.2}, synteny={:?})",
                cluster_name,
                genome_name,
                r.rescue_type,
                r.identity,
                r.coverage,
                r.confidence,
                r.synteny_ok
            );
        }
    }

    // ── 11. Nucleotide Similarity Clustering (optional) ──
    let nsc_result = if args.nscluster {
        let prevalences: Vec<f64> = (0..n_clusters)
            .map(|ci_val| pa.present_count(ClusterIdx(ci_val)) as f64 / n_genomes as f64)
            .collect();

        let result = match args.nsc_engine.as_str() {
            "mmseqs2" => {
                info!(
                    "Running NSC with mmseqs2 linclust (min_id={:.0}%, cov={:.0}%)...",
                    args.nsc_min_seq_id * 100.0,
                    args.nsc_coverage * 100.0
                );
                nscluster::run_nscluster_mmseqs2(
                    &pan_ref,
                    pa.cluster_names(),
                    &prevalences,
                    args.nsc_min_seq_id,
                    args.nsc_coverage,
                    args.threads,
                )?
            }
            _ => {
                info!(
                    "Running NSC with MinHash (k={}, threshold={:.2})...",
                    args.nsc_kmer_size, args.nsc_threshold
                );
                nscluster::run_nscluster(
                    &pan_ref,
                    pa.cluster_names(),
                    &prevalences,
                    args.nsc_kmer_size,
                    args.nsc_threshold,
                    args.threads,
                )
            }
        };

        if result.n_clusters > 0 {
            info!(
                "NSC: {} clusters merged into {} NS groups ({} genes total)",
                result.total_merged_genes, result.n_clusters, result.total_merged_genes
            );
        } else {
            info!("NSC: no redundant clusters found");
        }
        Some(result)
    } else {
        None
    };

    // ── 12. Generate outputs ──
    fs::create_dir_all(&args.output)?;

    // Sort rescues for deterministic output order regardless of thread scheduling
    rescues.sort_by(|a, b| {
        a.cluster_idx
            .cmp(&b.cluster_idx)
            .then(a.genome_idx.cmp(&b.genome_idx))
    });
    synteny_rejected.sort_by(|a, b| {
        a.cluster_idx
            .cmp(&b.cluster_idx)
            .then(a.genome_idx.cmp(&b.genome_idx))
    });
    let mut paralog_rejected = paralog_rejected;
    paralog_rejected.sort_by(|a, b| {
        a.cluster_idx
            .cmp(&b.cluster_idx)
            .then(a.genome_idx.cmp(&b.genome_idx))
    });
    let mobile_element_hits = if let Some(path) = args.mobile_element_reference.as_deref() {
        search_mobile_element_reference(MobileElementSearchInput {
            path,
            genome_files: &genome_files,
            pa: &pa,
            rescues: &rescues,
            synteny_rejected: &synteny_rejected,
            paralog_rejected: &paralog_rejected,
            threads: args.threads,
            min_identity: args.mobile_element_min_identity,
            min_coverage: args.mobile_element_min_coverage,
            preset: &args.preset,
        })?
    } else {
        Vec::new()
    };

    let all_rescue_map = build_best_rescue_map(rescues.iter());
    let matrix_rescue_map = build_best_rescue_map(rescues.iter().filter(|r| {
        let state = accepted_rescue_evidence_state(r, &pa, &gene_to_cluster_maps);
        matrix_policy_corrects(args.matrix_policy, state)
    }));

    write_report(&rescues, &pa, &args.output)?;
    if !synteny_rejected.is_empty() {
        write_named_report(
            &synteny_rejected,
            &pa,
            &args.output,
            "synteny_rejected_report.tsv",
        )?;
    }
    if !paralog_rejected.is_empty() {
        write_named_report(
            &paralog_rejected,
            &pa,
            &args.output,
            "paralog_rejected_report.tsv",
        )?;
    }
    write_orthology_evidence(
        &rescues,
        &synteny_rejected,
        &paralog_rejected,
        &pa,
        &gene_to_cluster_maps,
        args.matrix_policy,
        &args.output,
    )?;
    if let Some(reference_genome) = args.reference_genome.as_deref() {
        write_reference_genome_audit(ReferenceAuditInputs {
            requested_reference: reference_genome,
            rescues: &rescues,
            synteny_rejected: &synteny_rejected,
            paralog_rejected: &paralog_rejected,
            reference_candidates: &reference_audit_candidates,
            reference_cluster_events: &reference_cluster_events,
            pa: &pa,
            gene_to_cluster_maps: &gene_to_cluster_maps,
            matrix_policy: args.matrix_policy,
            output_dir: &args.output,
        })?;
        if !reference_cluster_events.is_empty() {
            write_reference_cluster_events(&reference_cluster_events, &pa, &args.output)?;
        }
    }
    let context_summaries = if !gff_indices.is_empty() {
        let (context_rows, summaries) =
            write_cluster_locus_contexts(&pa, &gff_indices, &gene_to_cluster_maps, &args.output)?;
        info!(
            "Cluster locus context audit: {} context rows across {} clusters",
            context_rows,
            summaries.len()
        );
        summaries
    } else {
        Vec::new()
    };
    write_corrected_csv(&pa, &matrix_rescue_map, &args.output, &nsc_result)?;
    write_multistate_matrix(MultistateMatrixInputs {
        pa: &pa,
        accepted_rescue_map: &all_rescue_map,
        nsc_rescue_map: &matrix_rescue_map,
        synteny_rejected: &synteny_rejected,
        paralog_rejected: &paralog_rejected,
        gene_to_cluster_maps: &gene_to_cluster_maps,
        output_dir: &args.output,
        nsc_result: &nsc_result,
    })?;
    let n_relationships = write_cluster_relationships(&pa, &rescues, &nsc_result, &args.output)?;
    if n_relationships > 0 {
        info!("Cluster relationships: {} rows", n_relationships);
    }
    write_evolutionary_events(EvolutionaryEventInputs {
        pa: &pa,
        matrix_rescue_map: &matrix_rescue_map,
        rescues: &rescues,
        synteny_rejected: &synteny_rejected,
        paralog_rejected: &paralog_rejected,
        reference_cluster_events: &reference_cluster_events,
        context_summaries: &context_summaries,
        gff_indices: &gff_indices,
        mobile_element_hits: &mobile_element_hits,
        mobile_element_window: args.mobile_element_window,
        phylo_tree: phylo_tree.as_ref(),
        reference_genome: args.reference_genome.as_deref(),
        nsc_result: &nsc_result,
        output_dir: &args.output,
    })?;
    if args.write_corrected_gff {
        let matrix_rescues: Vec<RescueEvent> =
            matrix_rescue_map.values().map(|r| (*r).clone()).collect();
        write_corrected_gffs(&matrix_rescues, &pa, &gff_paths, &args.output)?;
    }

    let stats = SummaryStats {
        core_before,
        total_gaps,
        filtered_gaps,
        skipped_by_prevalence,
        searchable_gaps,
        missing_fasta_gaps,
        dropped,
        missing_cds,
        split,
        fragmented,
        alignment_split,
        start_variant,
        paralog,
        pseudogene,
        synteny_rejected: synteny_rejected.len(),
        paralog_rejected: paralog_rejected.len(),
        matrix_policy: args.matrix_policy,
        matrix_corrected: matrix_rescue_map.len(),
        nsc_matrix_corrected: nsc_presence_updates(&pa, &nsc_result, &matrix_rescue_map).len(),
        synteny_enabled: args.synteny,
        nsc_clusters: nsc_result.as_ref().map(|r| r.n_clusters).unwrap_or(0),
        nsc_merged_genes: nsc_result
            .as_ref()
            .map(|r| r.total_merged_genes)
            .unwrap_or(0),
    };
    write_summary(
        &pa,
        &matrix_rescue_map,
        &rescues,
        &stats,
        &nsc_result,
        &args.output,
    )?;
    write_run_manifest(args, &pa, &stats, &args.output)?;

    // Write NSC report if clusters were merged
    if let Some(ref nsc) = nsc_result {
        if !nsc.merged_clusters.is_empty() {
            write_nsc_report(nsc, pa.cluster_names(), &args.output)?;
        }
    }

    let core_after = count_core_after(&pa, &matrix_rescue_map, &nsc_result);
    let soft_core_after = count_soft_core_after(&pa, &matrix_rescue_map, &nsc_result);
    let core_diff = core_after as i64 - core_before as i64;
    let soft_core_before = pa.soft_core_size();
    let soft_diff = soft_core_after as i64 - soft_core_before as i64;
    info!(
        "Core genome: {} → {} ({:+})",
        core_before, core_after, core_diff
    );
    info!(
        "Soft-core (≥95%): {} → {} ({:+})",
        soft_core_before, soft_core_after, soft_diff
    );
    info!("Output: {}", args.output.display());

    Ok(())
}

fn validate_percent_arg(name: &str, value: f64) -> Result<()> {
    if !value.is_finite() || !(0.0..=100.0).contains(&value) {
        anyhow::bail!("{} must be a finite percentage in the range 0..=100", name);
    }
    Ok(())
}

fn validate_fraction_arg(name: &str, value: f64) -> Result<()> {
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        anyhow::bail!("{} must be a finite fraction in the range 0..=1", name);
    }
    Ok(())
}

fn validate_args(args: &Args) -> Result<()> {
    if args.threads == 0 {
        anyhow::bail!("--threads must be at least 1");
    }
    validate_percent_arg("--min-identity", args.min_identity)?;
    validate_percent_arg("--min-coverage", args.min_coverage)?;
    validate_percent_arg(
        "--mobile-element-min-identity",
        args.mobile_element_min_identity,
    )?;
    validate_percent_arg(
        "--mobile-element-min-coverage",
        args.mobile_element_min_coverage,
    )?;
    validate_fraction_arg("--min-prevalence", args.min_prevalence)?;
    validate_fraction_arg("--nsc-threshold", args.nsc_threshold)?;
    validate_fraction_arg("--nsc-min-seq-id", args.nsc_min_seq_id)?;
    validate_fraction_arg("--nsc-coverage", args.nsc_coverage)?;
    if args.nsc_kmer_size == 0 {
        anyhow::bail!("--nsc-kmer-size must be at least 1");
    }
    if !matches!(args.nsc_engine.as_str(), "minhash" | "mmseqs2") {
        anyhow::bail!("--nsc-engine must be either 'minhash' or 'mmseqs2'");
    }
    Ok(())
}

// ── Output helpers ──

fn write_report(
    rescues: &[RescueEvent],
    pa: &crate::parser::PresenceAbsence,
    output_dir: &std::path::Path,
) -> Result<()> {
    write_named_report(rescues, pa, output_dir, "rescue_report.tsv")
}

fn write_named_report(
    rescues: &[RescueEvent],
    pa: &crate::parser::PresenceAbsence,
    output_dir: &std::path::Path,
    filename: &str,
) -> Result<()> {
    let path = output_dir.join(filename);
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "gene_cluster\tgenome\trescue_type\tidentity\tcoverage\tcontig\tstart\tend\tstrand\tconfidence\tsynteny\tsynteny_score\tsynteny_support_n\tsynteny_reference_contexts\tsynteny_neighbors\tsynteny_upstream_matches\tsynteny_downstream_matches\torf_integrity_signal\tannotation_match\tfragment_n\tfragment_query_gap_bp\tfragment_target_gap_bp\tfragment_colinear\tfragment_same_cds\tfragment_note"
    )?;
    for r in rescues {
        let frag = &r.fragment_evidence;
        let fragment_n = if frag.n_fragments > 0 {
            frag.n_fragments.to_string()
        } else {
            String::new()
        };
        let query_gap = if frag.n_fragments > 0 {
            frag.total_query_gap_bp.to_string()
        } else {
            String::new()
        };
        let target_gap = if frag.n_fragments > 0 {
            frag.total_target_gap_bp.to_string()
        } else {
            String::new()
        };
        let colinear = match frag.colinear {
            Some(true) => "TRUE",
            Some(false) => "FALSE",
            None => "",
        };
        let same_cds = match frag.same_cds {
            Some(true) => "TRUE",
            Some(false) => "FALSE",
            None => "",
        };
        write_tsv_fields(
            &mut w,
            &[
                pa.cluster_name(r.cluster_idx).to_string(),
                pa.genome_name(r.genome_idx).to_string(),
                r.rescue_type.to_string(),
                format!("{:.1}", r.identity),
                format!("{:.1}", r.coverage),
                r.contig.clone(),
                (r.start + 1).to_string(),
                r.end.to_string(),
                r.strand.to_string(),
                format!("{:.3}", r.confidence),
                synteny_label(r.synteny_ok).to_string(),
                synteny_score_value(&r.synteny_evidence),
                synteny_count_value(r.synteny_evidence.support, &r.synteny_evidence),
                synteny_count_value(r.synteny_evidence.reference_contexts, &r.synteny_evidence),
                synteny_count_value(r.synteny_evidence.n_neighbors, &r.synteny_evidence),
                synteny_count_value(r.synteny_evidence.matched_upstream, &r.synteny_evidence),
                synteny_count_value(r.synteny_evidence.matched_downstream, &r.synteny_evidence),
                rescue_orf_signal(r).to_string(),
                r.annotation_match.as_deref().unwrap_or("").to_string(),
                fragment_n,
                query_gap,
                target_gap,
                colinear.to_string(),
                same_cds.to_string(),
                r.fragment_evidence
                    .note
                    .as_deref()
                    .unwrap_or("")
                    .to_string(),
            ],
        )?;
    }
    w.flush()?;
    Ok(())
}

fn clean_tsv_field(value: &str) -> String {
    value
        .chars()
        .map(|ch| match ch {
            '\t' | '\n' | '\r' | '\0' => ' ',
            _ => ch,
        })
        .collect()
}

fn write_tsv_fields(w: &mut BufWriter<fs::File>, fields: &[String]) -> Result<()> {
    let cleaned: Vec<String> = fields.iter().map(|field| clean_tsv_field(field)).collect();
    writeln!(w, "{}", cleaned.join("\t"))?;
    Ok(())
}

fn synteny_label(synteny_ok: Option<bool>) -> &'static str {
    match synteny_ok {
        Some(true) => "PASS",
        Some(false) => "FAIL",
        None => "",
    }
}

fn synteny_score_value(evidence: &SyntenyEvidence) -> String {
    evidence
        .score
        .map(|score| format!("{:.3}", score))
        .unwrap_or_default()
}

fn synteny_count_value(value: usize, evidence: &SyntenyEvidence) -> String {
    if evidence.score.is_some() {
        value.to_string()
    } else {
        String::new()
    }
}

// Keep the documented Rust MSRV; u64::is_multiple_of is newer.
#[allow(clippy::manual_is_multiple_of)]
fn rescue_orf_signal(rescue: &RescueEvent) -> &'static str {
    match rescue.rescue_type {
        RescueType::Pseudogene { .. } => "disrupted_orf",
        RescueType::Fragmented { .. } | RescueType::AlignmentSplit { .. } => "fragmented_orf",
        _ => {
            let span = rescue.end.saturating_sub(rescue.start);
            if span == 0 {
                "unknown"
            } else if span % 3 == 0 {
                "in_frame_span"
            } else {
                "out_of_frame_span_review"
            }
        }
    }
}

fn stripped_gene_id(gene_id: &str) -> Option<&str> {
    panaroo_gene_id_suffix_variants(gene_id).into_iter().nth(1)
}

fn gff_gene_for_id<'a>(gff: &'a GffIndex, gene_id: &str) -> Option<&'a GffGene> {
    gff.get_by_id(gene_id)
        .or_else(|| stripped_gene_id(gene_id).and_then(|id| gff.get_by_id(id)))
}

fn gff_gene_for_cell<'a>(gff: &'a GffIndex, cell: &str) -> Option<&'a GffGene> {
    cell.split(';')
        .map(str::trim)
        .filter(|id| !id.is_empty())
        .find_map(|id| gff_gene_for_id(gff, id))
}

fn annotation_cluster_name(
    rescue: &RescueEvent,
    pa: &crate::parser::PresenceAbsence,
    gene_to_cluster_maps: &HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
) -> String {
    let Some(annotation_match) = rescue.annotation_match.as_deref() else {
        return String::new();
    };
    let Some(g2c) = gene_to_cluster_maps.get(&rescue.genome_idx) else {
        return String::new();
    };
    g2c.get(annotation_match)
        .map(|ci| pa.cluster_name(*ci).to_string())
        .unwrap_or_default()
}

fn accepted_orthology_status(
    rescue: &RescueEvent,
    query_cluster: &str,
    annotation_cluster: &str,
) -> &'static str {
    if rescue.synteny_ok == Some(true) {
        return "accepted_synteny";
    }
    if !annotation_cluster.is_empty() && annotation_cluster == query_cluster {
        return "accepted_annotation_same_cluster";
    }
    match rescue.rescue_type {
        RescueType::Split { .. } => "accepted_split_alias",
        RescueType::StartVariant { .. } => "accepted_boundary_variant",
        RescueType::AlignmentSplit { .. } => "accepted_alignment_split",
        RescueType::Fragmented { .. } => "accepted_fragment_evidence",
        RescueType::Pseudogene { .. } => "accepted_pseudogene_evidence",
        RescueType::Paralog { .. } => "accepted_paralog_sequence_evidence",
        RescueType::MissingCds => {
            if rescue.synteny_ok == Some(true) {
                "accepted_missing_cds_synteny"
            } else {
                "accepted_missing_cds_sequence"
            }
        }
        RescueType::Dropped => "accepted_sequence_evidence",
    }
}

fn is_paralog_without_synteny_support(rescue: &RescueEvent) -> bool {
    matches!(rescue.rescue_type, RescueType::Paralog { .. }) && rescue.synteny_ok != Some(true)
}

fn evidence_state(
    rescue: &RescueEvent,
    accepted: bool,
    query_cluster: &str,
    annotation_cluster: &str,
    rejection_reason: &str,
) -> &'static str {
    if !accepted {
        if rejection_reason == "synteny_fail" {
            return "REJECTED_CONTEXT";
        }
        if matches!(rescue.rescue_type, RescueType::Paralog { .. }) {
            return "HOMOLOG_PARALOG";
        }
        return "REJECTED_CONTEXT";
    }

    match rescue.rescue_type {
        RescueType::Split { .. } => "ALIAS_CLUSTER",
        RescueType::Paralog { .. } => {
            if rescue.synteny_ok == Some(true) {
                "STRICT_PRESENT"
            } else {
                "HOMOLOG_PARALOG"
            }
        }
        RescueType::Pseudogene { .. } => "PSEUDOGENE",
        RescueType::Fragmented { .. } | RescueType::AlignmentSplit { .. } => "PARTIAL_HOMOLOGY",
        RescueType::MissingCds => {
            if rescue.synteny_ok == Some(true) {
                "STRICT_PRESENT"
            } else {
                "MISSING_CDS_UNVALIDATED"
            }
        }
        RescueType::Dropped | RescueType::StartVariant { .. } => {
            if rescue.synteny_ok == Some(true)
                || (!annotation_cluster.is_empty() && annotation_cluster == query_cluster)
            {
                "STRICT_PRESENT"
            } else {
                "SEQUENCE_ONLY"
            }
        }
    }
}

fn matrix_policy_corrects(policy: MatrixPolicy, state: &str) -> bool {
    match policy {
        MatrixPolicy::Strict => matches!(state, "STRICT_PRESENT" | "ALIAS_CLUSTER"),
        MatrixPolicy::Balanced => {
            matches!(state, "STRICT_PRESENT" | "ALIAS_CLUSTER" | "SEQUENCE_ONLY")
        }
        MatrixPolicy::Sensitive => matches!(
            state,
            "STRICT_PRESENT"
                | "ALIAS_CLUSTER"
                | "SEQUENCE_ONLY"
                | "PARTIAL_HOMOLOGY"
                | "PSEUDOGENE"
        ),
    }
}

fn matrix_action(policy: MatrixPolicy, accepted: bool, state: &str) -> &'static str {
    if accepted && matrix_policy_corrects(policy, state) {
        "correct_presence"
    } else if accepted {
        "audit_only_policy"
    } else {
        "audit_only"
    }
}

fn accepted_rescue_evidence_state(
    rescue: &RescueEvent,
    pa: &crate::parser::PresenceAbsence,
    gene_to_cluster_maps: &HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
) -> &'static str {
    let query_cluster = pa.cluster_name(rescue.cluster_idx);
    let annotation_cluster = annotation_cluster_name(rescue, pa, gene_to_cluster_maps);
    evidence_state(rescue, true, query_cluster, &annotation_cluster, "")
}

fn publication_status(accepted: bool, state: &str) -> &'static str {
    if !accepted {
        return "audit_only_not_corrected";
    }
    match state {
        "STRICT_PRESENT" | "ALIAS_CLUSTER" => "publication_ready_correction",
        "SEQUENCE_ONLY" | "PARTIAL_HOMOLOGY" | "PSEUDOGENE" => {
            "corrected_but_review_for_publication"
        }
        "HOMOLOG_PARALOG" => "corrected_paralog_review",
        _ => "corrected_review",
    }
}

fn evidence_strength(state: &str) -> &'static str {
    match state {
        "STRICT_PRESENT" | "ALIAS_CLUSTER" => "high",
        "SEQUENCE_ONLY" => "moderate",
        "PARTIAL_HOMOLOGY" | "PSEUDOGENE" => "review",
        "PARTIAL_PARALOG_EVIDENCE" => "review",
        "EVENT_DERIVED_REFERENCE_CLUSTER" => "review",
        "MISSING_CDS_UNVALIDATED" => "review",
        "HOMOLOG_PARALOG" | "REJECTED_CONTEXT" => "low",
        "NO_EVIDENCE" => "none",
        _ => "review",
    }
}

fn recommended_wording(state: &str, accepted: bool) -> &'static str {
    if !accepted {
        return "Do not correct the matrix; keep as audit evidence only.";
    }
    match state {
        "STRICT_PRESENT" => "Treat as present at the same locus.",
        "ALIAS_CLUSTER" => {
            "Treat as a cluster alias or split; avoid counting as an independent novel gene."
        }
        "SEQUENCE_ONLY" => "Treat as sequence-supported presence; review if used as a key claim.",
        "PARTIAL_HOMOLOGY" => {
            "Describe as partial or fragmented homology, not clean gene presence."
        }
        "PSEUDOGENE" => "Describe as disrupted or pseudogenized sequence evidence.",
        "MISSING_CDS_UNVALIDATED" => {
            "Sequence is present but the CDS is missing from the annotation; require synteny support before matrix correction."
        }
        "HOMOLOG_PARALOG" => "Keep distinct unless synteny confirms the same locus.",
        "PARTIAL_PARALOG_EVIDENCE" => {
            "Treat as partial/paralogous reference-like evidence; do not convert to binary presence."
        }
        "EVENT_DERIVED_REFERENCE_CLUSTER" => {
            "Keep as a separate cluster, but interpret as an event-derived reference-gene variant."
        }
        _ => "Review before publication.",
    }
}

fn write_orthology_evidence(
    rescues: &[RescueEvent],
    synteny_rejected: &[RescueEvent],
    paralog_rejected: &[RescueEvent],
    pa: &crate::parser::PresenceAbsence,
    gene_to_cluster_maps: &HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
    matrix_policy: MatrixPolicy,
    output_dir: &std::path::Path,
) -> Result<()> {
    let path = output_dir.join("orthology_evidence.tsv");
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "decision\torthology_status\tevidence_state\tmatrix_action\tpublication_status\tevidence_strength\trecommended_wording\tgene_cluster\tgenome\trescue_type\tidentity\tcoverage\tcontig\tstart\tend\tstrand\tconfidence\tsynteny\tsynteny_score\tsynteny_support_n\tsynteny_reference_contexts\tsynteny_neighbors\tsynteny_upstream_matches\tsynteny_downstream_matches\torf_integrity_signal\tannotation_match\tannotation_cluster\trejection_reason\tfragment_note"
    )?;

    for rescue in rescues {
        let query_cluster = pa.cluster_name(rescue.cluster_idx);
        let annotation_cluster = annotation_cluster_name(rescue, pa, gene_to_cluster_maps);
        let status = accepted_orthology_status(rescue, query_cluster, &annotation_cluster);
        let state = evidence_state(rescue, true, query_cluster, &annotation_cluster, "");
        write_tsv_fields(
            &mut w,
            &[
                "accepted".to_string(),
                status.to_string(),
                state.to_string(),
                matrix_action(matrix_policy, true, state).to_string(),
                publication_status(true, state).to_string(),
                evidence_strength(state).to_string(),
                recommended_wording(state, true).to_string(),
                query_cluster.to_string(),
                pa.genome_name(rescue.genome_idx).to_string(),
                rescue.rescue_type.to_string(),
                format!("{:.1}", rescue.identity),
                format!("{:.1}", rescue.coverage),
                rescue.contig.clone(),
                (rescue.start + 1).to_string(),
                rescue.end.to_string(),
                rescue.strand.to_string(),
                format!("{:.3}", rescue.confidence),
                synteny_label(rescue.synteny_ok).to_string(),
                synteny_score_value(&rescue.synteny_evidence),
                synteny_count_value(rescue.synteny_evidence.support, &rescue.synteny_evidence),
                synteny_count_value(
                    rescue.synteny_evidence.reference_contexts,
                    &rescue.synteny_evidence,
                ),
                synteny_count_value(
                    rescue.synteny_evidence.n_neighbors,
                    &rescue.synteny_evidence,
                ),
                synteny_count_value(
                    rescue.synteny_evidence.matched_upstream,
                    &rescue.synteny_evidence,
                ),
                synteny_count_value(
                    rescue.synteny_evidence.matched_downstream,
                    &rescue.synteny_evidence,
                ),
                rescue_orf_signal(rescue).to_string(),
                rescue.annotation_match.as_deref().unwrap_or("").to_string(),
                annotation_cluster,
                String::new(),
                rescue
                    .fragment_evidence
                    .note
                    .as_deref()
                    .unwrap_or("")
                    .to_string(),
            ],
        )?;
    }

    for (rejected, status, reason) in [
        (
            synteny_rejected,
            "rejected_incompatible_context",
            "synteny_fail",
        ),
        (
            paralog_rejected,
            "rejected_paralog_without_synteny",
            "paralog_without_synteny",
        ),
    ] {
        for rescue in rejected {
            let query_cluster = pa.cluster_name(rescue.cluster_idx);
            let annotation_cluster = annotation_cluster_name(rescue, pa, gene_to_cluster_maps);
            let state = evidence_state(rescue, false, query_cluster, &annotation_cluster, reason);
            write_tsv_fields(
                &mut w,
                &[
                    "rejected".to_string(),
                    status.to_string(),
                    state.to_string(),
                    matrix_action(matrix_policy, false, state).to_string(),
                    publication_status(false, state).to_string(),
                    evidence_strength(state).to_string(),
                    recommended_wording(state, false).to_string(),
                    query_cluster.to_string(),
                    pa.genome_name(rescue.genome_idx).to_string(),
                    rescue.rescue_type.to_string(),
                    format!("{:.1}", rescue.identity),
                    format!("{:.1}", rescue.coverage),
                    rescue.contig.clone(),
                    (rescue.start + 1).to_string(),
                    rescue.end.to_string(),
                    rescue.strand.to_string(),
                    format!("{:.3}", rescue.confidence),
                    synteny_label(rescue.synteny_ok).to_string(),
                    synteny_score_value(&rescue.synteny_evidence),
                    synteny_count_value(rescue.synteny_evidence.support, &rescue.synteny_evidence),
                    synteny_count_value(
                        rescue.synteny_evidence.reference_contexts,
                        &rescue.synteny_evidence,
                    ),
                    synteny_count_value(
                        rescue.synteny_evidence.n_neighbors,
                        &rescue.synteny_evidence,
                    ),
                    synteny_count_value(
                        rescue.synteny_evidence.matched_upstream,
                        &rescue.synteny_evidence,
                    ),
                    synteny_count_value(
                        rescue.synteny_evidence.matched_downstream,
                        &rescue.synteny_evidence,
                    ),
                    rescue_orf_signal(rescue).to_string(),
                    rescue.annotation_match.as_deref().unwrap_or("").to_string(),
                    annotation_cluster,
                    reason.to_string(),
                    rescue
                        .fragment_evidence
                        .note
                        .as_deref()
                        .unwrap_or("")
                        .to_string(),
                ],
            )?;
        }
    }
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests;
