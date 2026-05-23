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
    resolve_conflicts, validate_synteny_with_cache, GffGene, GffIndex, PseudogeneReason,
    RescueEvent, RescueType, SyntenyEvidence, SyntenyReferenceCache,
};
use crate::search::{batch_search, Hit};
use crate::types::{ClusterIdx, GenomeIdx};

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

    let all_rescue_map: HashMap<(ClusterIdx, GenomeIdx), &RescueEvent> = rescues
        .iter()
        .map(|r| ((r.cluster_idx, r.genome_idx), r))
        .collect();
    let matrix_rescue_map: HashMap<(ClusterIdx, GenomeIdx), &RescueEvent> = rescues
        .iter()
        .filter(|r| {
            let state = accepted_rescue_evidence_state(r, &pa, &gene_to_cluster_maps);
            matrix_policy_corrects(args.matrix_policy, state)
        })
        .map(|r| ((r.cluster_idx, r.genome_idx), r))
        .collect();

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
    gene_id
        .strip_suffix("_pseudo")
        .or_else(|| gene_id.rfind("_dup_").map(|pos| &gene_id[..pos]))
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

fn resolve_reference_genome(
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

fn reference_status_for(state: &str, accepted: bool) -> &'static str {
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

fn reference_publication_status(state: &str, accepted: bool) -> &'static str {
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

fn reference_recommended_wording(state: &str, accepted: bool) -> &'static str {
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

fn reference_publication_confidence(state: &str, reference_status: &str) -> &'static str {
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

fn reference_shadow_state(
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

fn best_overlap_gene_id(gff: &GffIndex, hit: &Hit) -> Option<String> {
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

fn reference_hit_shadow_candidates(
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

fn reference_candidate_score(rescue: &RescueEvent) -> f64 {
    let synteny_bonus = match rescue.synteny_ok {
        Some(true) => 2.0,
        Some(false) => -1.0,
        None => 0.0,
    };
    rescue.confidence + (rescue.identity / 100.0) + (rescue.coverage / 100.0) + synteny_bonus
}

fn best_reference_candidate(
    candidates: &[RescueEvent],
    ci: ClusterIdx,
    reference_gi: GenomeIdx,
) -> Option<&RescueEvent> {
    candidates
        .iter()
        .filter(|r| r.cluster_idx == ci && r.genome_idx == reference_gi)
        .max_by(|a, b| reference_candidate_score(a).total_cmp(&reference_candidate_score(b)))
}

const REFERENCE_EVENT_KMER: usize = 15;
const REFERENCE_EVENT_TOP_CANDIDATES: usize = 12;
const REFERENCE_EVENT_MAX_DP_CELLS: usize = 5_000_000;

#[derive(Debug, Clone)]
struct ReferenceClusterEvent {
    query_cluster_idx: ClusterIdx,
    target_cluster_idx: ClusterIdx,
    target_gene_id: String,
    event_class: &'static str,
    identity: f64,
    query_coverage: f64,
    reference_coverage: f64,
    query_length: usize,
    reference_length: usize,
    event_size_nt: i64,
    event_size_aa: Option<i64>,
    matches: usize,
    mismatches: usize,
    gaps_in_query: usize,
    gaps_in_reference: usize,
    shared_kmers: usize,
    interpretation: &'static str,
}

#[derive(Debug, Clone)]
struct PairwiseEventStats {
    identity: f64,
    query_coverage: f64,
    reference_coverage: f64,
    matches: usize,
    mismatches: usize,
    gaps_in_query: usize,
    gaps_in_reference: usize,
}

fn longest_cluster_sequence<'a>(
    pan_ref: &'a HashMap<String, Vec<Vec<u8>>>,
    cluster_name: &str,
) -> Option<&'a [u8]> {
    pan_ref
        .get(cluster_name)?
        .iter()
        .max_by_key(|seq| seq.len())
        .map(Vec::as_slice)
}

fn base_bits(base: u8) -> Option<u64> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn canonical_kmer_hash(seq: &[u8]) -> Option<u64> {
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

fn sequence_kmers(seq: &[u8], k: usize) -> Vec<u64> {
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

fn global_event_alignment(query: &[u8], reference: &[u8]) -> Option<PairwiseEventStats> {
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

fn classify_reference_event(
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

fn is_h37rv_locus_token(token: &str) -> bool {
    let Some(rest) = token.strip_prefix("Rv") else {
        return false;
    };
    if rest.len() < 4 || !rest[..4].chars().all(|ch| ch.is_ascii_digit()) {
        return false;
    }
    rest[4..].chars().all(|ch| ch.is_ascii_alphabetic())
}

fn display_reference_gene_id(gene_id: &str) -> String {
    let cleaned = stripped_gene_id(gene_id.trim()).unwrap_or(gene_id.trim());
    cleaned
        .split('_')
        .find(|part| is_h37rv_locus_token(part))
        .unwrap_or(cleaned)
        .to_string()
}

fn reference_gene_id(pa: &crate::parser::PresenceAbsence, ci: ClusterIdx, gi: GenomeIdx) -> String {
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

fn build_reference_cluster_events(
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
                .is_none_or(|best| candidate_score > reference_event_score(best))
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

fn reference_event_score(event: &ReferenceClusterEvent) -> f64 {
    event.identity
        + event.query_coverage
        + event.reference_coverage
        + (event.shared_kmers as f64 / 10.0)
}

struct ReferenceEventRow<'a> {
    reference_name: &'a str,
    rescue: &'a RescueEvent,
    gene_to_cluster_maps: &'a HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
    matrix_policy: MatrixPolicy,
    decision: &'a str,
    accepted: bool,
    rejection_reason: &'a str,
}

struct ReferenceAuditInputs<'a> {
    requested_reference: &'a str,
    rescues: &'a [RescueEvent],
    synteny_rejected: &'a [RescueEvent],
    paralog_rejected: &'a [RescueEvent],
    reference_candidates: &'a [RescueEvent],
    reference_cluster_events: &'a [ReferenceClusterEvent],
    pa: &'a crate::parser::PresenceAbsence,
    gene_to_cluster_maps: &'a HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
    matrix_policy: MatrixPolicy,
    output_dir: &'a std::path::Path,
}

fn write_reference_event_row(
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
    let final_status = reference_publication_status(state, row.accepted);
    let reference_status = reference_status_for(state, row.accepted);
    let publication_confidence = reference_publication_confidence(state, reference_status);
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
            matrix_action(row.matrix_policy, row.accepted, state).to_string(),
            state.to_string(),
            evidence_strength(state).to_string(),
            publication_confidence.to_string(),
            reference_recommended_wording(state, row.accepted).to_string(),
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

fn write_reference_shadow_candidate_row(
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

fn write_reference_cluster_event_row(
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

fn write_reference_genome_audit(inputs: ReferenceAuditInputs<'_>) -> Result<()> {
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

        if let Some(rescue) = inputs
            .rescues
            .iter()
            .find(|r| r.cluster_idx == ci && r.genome_idx == reference_gi)
        {
            write_reference_event_row(
                &mut w,
                inputs.pa,
                ReferenceEventRow {
                    reference_name,
                    rescue,
                    gene_to_cluster_maps: inputs.gene_to_cluster_maps,
                    matrix_policy: inputs.matrix_policy,
                    decision: "accepted",
                    accepted: true,
                    rejection_reason: "",
                },
            )?;
            continue;
        }

        if let Some(rescue) = inputs
            .synteny_rejected
            .iter()
            .find(|r| r.cluster_idx == ci && r.genome_idx == reference_gi)
        {
            write_reference_event_row(
                &mut w,
                inputs.pa,
                ReferenceEventRow {
                    reference_name,
                    rescue,
                    gene_to_cluster_maps: inputs.gene_to_cluster_maps,
                    matrix_policy: inputs.matrix_policy,
                    decision: "rejected",
                    accepted: false,
                    rejection_reason: "synteny_fail",
                },
            )?;
            continue;
        }

        if let Some(rescue) = inputs
            .paralog_rejected
            .iter()
            .find(|r| r.cluster_idx == ci && r.genome_idx == reference_gi)
        {
            write_reference_event_row(
                &mut w,
                inputs.pa,
                ReferenceEventRow {
                    reference_name,
                    rescue,
                    gene_to_cluster_maps: inputs.gene_to_cluster_maps,
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

fn write_reference_cluster_events(
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

#[derive(Debug, Clone)]
struct PhyloNode {
    label: Option<String>,
    genome_idx: Option<GenomeIdx>,
    children: Vec<PhyloNode>,
}

#[derive(Debug, Clone)]
struct PhyloTree {
    root: PhyloNode,
}

#[derive(Debug, Clone, Default)]
struct PhyloPatternSummary {
    assessed: bool,
    matched_tips: usize,
    present_tips: usize,
    absent_tips: usize,
    parsimony_changes: usize,
    gains: usize,
    losses: usize,
    root_state: &'static str,
    carrier_monophyly: bool,
    max_carrier_clade_fraction: f64,
    pattern: &'static str,
}

#[derive(Debug, Clone, Default)]
struct FitchNodeSummary {
    state_set: u8,
    matched_tips: usize,
    present_tips: usize,
    changes: usize,
    children: Vec<FitchNodeSummary>,
}

struct NewickParser<'a> {
    input: &'a str,
    pos: usize,
}

impl<'a> NewickParser<'a> {
    fn new(input: &'a str) -> Self {
        Self { input, pos: 0 }
    }

    fn parse(mut self) -> Result<PhyloNode> {
        self.skip_ignored();
        let root = self.parse_subtree()?;
        self.skip_ignored();
        if self.peek() == Some(';') {
            self.pos += 1;
        }
        self.skip_ignored();
        if self.pos < self.input.len() {
            anyhow::bail!("Unexpected Newick content at byte {}", self.pos);
        }
        Ok(root)
    }

    fn parse_subtree(&mut self) -> Result<PhyloNode> {
        self.skip_ignored();
        if self.peek() == Some('(') {
            self.pos += 1;
            let mut children = Vec::new();
            loop {
                children.push(self.parse_subtree()?);
                self.skip_ignored();
                match self.peek() {
                    Some(',') => {
                        self.pos += 1;
                    }
                    Some(')') => {
                        self.pos += 1;
                        break;
                    }
                    other => anyhow::bail!(
                        "Expected ',' or ')' in Newick at byte {}, got {:?}",
                        self.pos,
                        other
                    ),
                }
            }
            let label = self.parse_optional_label();
            self.skip_branch_length();
            Ok(PhyloNode {
                label,
                genome_idx: None,
                children,
            })
        } else {
            let label = self.parse_optional_label();
            if label.as_deref().unwrap_or("").is_empty() {
                anyhow::bail!("Expected Newick leaf label at byte {}", self.pos);
            }
            self.skip_branch_length();
            Ok(PhyloNode {
                label,
                genome_idx: None,
                children: Vec::new(),
            })
        }
    }

    fn parse_optional_label(&mut self) -> Option<String> {
        self.skip_ignored();
        match self.peek() {
            Some('\'') => {
                self.pos += 1;
                let start = self.pos;
                while let Some(ch) = self.peek() {
                    if ch == '\'' {
                        let label = self.input[start..self.pos].to_string();
                        self.pos += 1;
                        return Some(label);
                    }
                    self.pos += ch.len_utf8();
                }
                Some(self.input[start..].to_string())
            }
            Some(':') | Some(',') | Some(')') | Some(';') | None => None,
            _ => {
                let start = self.pos;
                while let Some(ch) = self.peek() {
                    if matches!(ch, ':' | ',' | ')' | ';' | '[') || ch.is_whitespace() {
                        break;
                    }
                    self.pos += ch.len_utf8();
                }
                (self.pos > start).then(|| self.input[start..self.pos].to_string())
            }
        }
    }

    fn skip_branch_length(&mut self) {
        self.skip_ignored();
        if self.peek() != Some(':') {
            return;
        }
        self.pos += 1;
        while let Some(ch) = self.peek() {
            if matches!(ch, ',' | ')' | ';' | '[') || ch.is_whitespace() {
                break;
            }
            self.pos += ch.len_utf8();
        }
    }

    fn skip_ignored(&mut self) {
        loop {
            while self.peek().is_some_and(char::is_whitespace) {
                self.pos += self.peek().unwrap().len_utf8();
            }
            if self.peek() != Some('[') {
                break;
            }
            self.pos += 1;
            while let Some(ch) = self.peek() {
                self.pos += ch.len_utf8();
                if ch == ']' {
                    break;
                }
            }
        }
    }

    fn peek(&self) -> Option<char> {
        self.input[self.pos..].chars().next()
    }
}

fn load_phylo_tree(path: &Path, pa: &crate::parser::PresenceAbsence) -> Result<PhyloTree> {
    let text =
        fs::read_to_string(path).with_context(|| format!("Cannot read tree {}", path.display()))?;
    let mut root = NewickParser::new(&text)
        .parse()
        .with_context(|| format!("Cannot parse Newick tree {}", path.display()))?;
    let (matched_tips, unmatched_tips) = resolve_phylo_tips(&mut root, pa);
    info!(
        "Loaded tree {} with {} matched tips and {} unmatched tips",
        path.display(),
        matched_tips,
        unmatched_tips
    );
    if matched_tips == 0 {
        warn!("No tree tips matched matrix genome names; phylogenetic audit will be uninformative");
    } else if matched_tips < pa.n_genomes() {
        warn!(
            "Tree tip matching is incomplete: {}/{} matrix genomes matched; unmatched genomes are ignored in evolutionary event calls",
            matched_tips,
            pa.n_genomes()
        );
    }
    Ok(PhyloTree { root })
}

fn resolve_phylo_tips(node: &mut PhyloNode, pa: &crate::parser::PresenceAbsence) -> (usize, usize) {
    if node.children.is_empty() {
        node.genome_idx = node
            .label
            .as_deref()
            .and_then(|label| resolve_tree_tip(label, pa));
        return if node.genome_idx.is_some() {
            (1, 0)
        } else {
            (0, 1)
        };
    }
    let mut matched = 0usize;
    let mut unmatched = 0usize;
    for child in &mut node.children {
        let (m, u) = resolve_phylo_tips(child, pa);
        matched += m;
        unmatched += u;
    }
    (matched, unmatched)
}

fn resolve_tree_tip(label: &str, pa: &crate::parser::PresenceAbsence) -> Option<GenomeIdx> {
    for gi_val in 0..pa.n_genomes() {
        let gi = GenomeIdx(gi_val);
        if pa.genome_name(gi) == label {
            return Some(gi);
        }
    }
    let tip_as = extract_as_id(label);
    for gi_val in 0..pa.n_genomes() {
        let gi = GenomeIdx(gi_val);
        if let (Some(a), Some(b)) = (tip_as, extract_as_id(pa.genome_name(gi))) {
            if a == b {
                return Some(gi);
            }
        }
    }
    None
}

fn fitch_postorder(node: &PhyloNode, present_by_genome: &[bool]) -> FitchNodeSummary {
    const ABSENT: u8 = 0b01;
    const PRESENT: u8 = 0b10;

    if node.children.is_empty() {
        let Some(gi) = node.genome_idx else {
            return FitchNodeSummary::default();
        };
        let present = present_by_genome.get(gi.val()).copied().unwrap_or(false);
        return FitchNodeSummary {
            state_set: if present { PRESENT } else { ABSENT },
            matched_tips: 1,
            present_tips: usize::from(present),
            changes: 0,
            children: Vec::new(),
        };
    }

    let mut children = Vec::new();
    let mut state_set = 0u8;
    let mut matched_tips = 0usize;
    let mut present_tips = 0usize;
    let mut changes = 0usize;
    for child in &node.children {
        let summary = fitch_postorder(child, present_by_genome);
        if summary.matched_tips == 0 {
            children.push(summary);
            continue;
        }
        matched_tips += summary.matched_tips;
        present_tips += summary.present_tips;
        changes += summary.changes;
        if state_set == 0 {
            state_set = summary.state_set;
        } else {
            let intersection = state_set & summary.state_set;
            if intersection == 0 {
                state_set |= summary.state_set;
                changes += 1;
            } else {
                state_set = intersection;
            }
        }
        children.push(summary);
    }

    FitchNodeSummary {
        state_set,
        matched_tips,
        present_tips,
        changes,
        children,
    }
}

fn state_set_contains(state_set: u8, present: bool) -> bool {
    let bit = if present { 0b10 } else { 0b01 };
    state_set & bit != 0
}

fn choose_root_state(summary: &FitchNodeSummary) -> bool {
    if summary.state_set == 0b10 {
        true
    } else if summary.state_set == 0b01 {
        false
    } else {
        summary.present_tips * 2 >= summary.matched_tips
    }
}

fn assign_fitch_states(
    summary: &FitchNodeSummary,
    parent_state: Option<bool>,
    gains: &mut usize,
    losses: &mut usize,
) -> Option<bool> {
    if summary.matched_tips == 0 || summary.state_set == 0 {
        return None;
    }
    let state = if let Some(parent) = parent_state {
        if state_set_contains(summary.state_set, parent) {
            parent
        } else {
            !parent
        }
    } else {
        choose_root_state(summary)
    };

    if let Some(parent) = parent_state {
        if !parent && state {
            *gains += 1;
        } else if parent && !state {
            *losses += 1;
        }
    }
    for child in &summary.children {
        assign_fitch_states(child, Some(state), gains, losses);
    }
    Some(state)
}

fn carrier_clade_stats(summary: &FitchNodeSummary, total_present: usize) -> (bool, f64) {
    if total_present == 0 {
        return (false, 0.0);
    }
    if total_present == 1 {
        return (true, 1.0);
    }

    fn visit(summary: &FitchNodeSummary, total_present: usize, best: &mut usize, exact: &mut bool) {
        if summary.matched_tips > 0
            && summary.present_tips > 0
            && summary.present_tips == summary.matched_tips
        {
            *best = (*best).max(summary.present_tips);
            if summary.present_tips == total_present {
                *exact = true;
            }
        }
        for child in &summary.children {
            visit(child, total_present, best, exact);
        }
    }

    let mut best = 0usize;
    let mut exact = false;
    visit(summary, total_present, &mut best, &mut exact);
    (exact, best as f64 / total_present as f64)
}

fn classify_phylogenetic_pattern(summary: &PhyloPatternSummary) -> &'static str {
    if !summary.assessed {
        return "not_assessed_no_tree";
    }
    if summary.present_tips == 0 {
        return "absent_from_tree";
    }
    if summary.absent_tips == 0 {
        return "ancestral_retained";
    }
    if summary.present_tips == 1 {
        return "singleton_or_noise";
    }
    if summary.parsimony_changes >= 6 {
        return "homoplastic_gain_loss";
    }
    if summary.root_state == "present" && summary.gains == 0 && summary.losses > 0 {
        return "ancestral_loss";
    }
    if summary.root_state == "absent" && summary.losses == 0 && summary.gains == 1 {
        if summary.carrier_monophyly {
            return "lineage_specific_gain";
        }
        return "single_gain_with_incomplete_clade_support";
    }
    if summary.gains > 1 && summary.losses == 0 {
        return "recurrent_gain";
    }
    if summary.losses > 1 && summary.gains == 0 {
        return "recurrent_loss";
    }
    if summary.carrier_monophyly {
        return "clade_associated_presence";
    }
    "dispersed_gain_loss"
}

fn summarize_phylogenetic_pattern(
    tree: Option<&PhyloTree>,
    present_by_genome: &[bool],
) -> PhyloPatternSummary {
    let Some(tree) = tree else {
        return PhyloPatternSummary {
            pattern: "not_assessed_no_tree",
            root_state: "unknown",
            ..PhyloPatternSummary::default()
        };
    };
    let fitch = fitch_postorder(&tree.root, present_by_genome);
    if fitch.matched_tips == 0 {
        return PhyloPatternSummary {
            pattern: "not_assessed_no_matched_tips",
            root_state: "unknown",
            ..PhyloPatternSummary::default()
        };
    }
    let mut gains = 0usize;
    let mut losses = 0usize;
    let root_present = assign_fitch_states(&fitch, None, &mut gains, &mut losses).unwrap_or(false);
    let (carrier_monophyly, max_carrier_clade_fraction) =
        carrier_clade_stats(&fitch, fitch.present_tips);
    let mut summary = PhyloPatternSummary {
        assessed: true,
        matched_tips: fitch.matched_tips,
        present_tips: fitch.present_tips,
        absent_tips: fitch.matched_tips.saturating_sub(fitch.present_tips),
        parsimony_changes: fitch.changes,
        gains,
        losses,
        root_state: if root_present { "present" } else { "absent" },
        carrier_monophyly,
        max_carrier_clade_fraction,
        pattern: "unclassified",
    };
    summary.pattern = classify_phylogenetic_pattern(&summary);
    summary
}

#[derive(Debug, Clone)]
struct MobileElementHit {
    genome_idx: GenomeIdx,
    query_name: String,
    contig: String,
    start: u64,
    end: u64,
    strand: char,
    identity: f64,
    coverage: f64,
}

#[derive(Debug, Clone)]
struct MobileElementEvidence {
    genome_name: String,
    source: &'static str,
    kind: &'static str,
    label: String,
    distance_bp: u64,
    identity: Option<f64>,
    coverage: Option<f64>,
}

#[derive(Debug, Clone)]
struct MechanismCall {
    mechanism: &'static str,
    detail: String,
    support: &'static str,
    wording: &'static str,
    mobile: Option<MobileElementEvidence>,
}

struct EvolutionaryEventInputs<'a> {
    pa: &'a crate::parser::PresenceAbsence,
    matrix_rescue_map: &'a HashMap<(ClusterIdx, GenomeIdx), &'a RescueEvent>,
    rescues: &'a [RescueEvent],
    synteny_rejected: &'a [RescueEvent],
    paralog_rejected: &'a [RescueEvent],
    reference_cluster_events: &'a [ReferenceClusterEvent],
    context_summaries: &'a [ClusterContextSummary],
    gff_indices: &'a HashMap<GenomeIdx, GffIndex>,
    mobile_element_hits: &'a [MobileElementHit],
    mobile_element_window: u64,
    phylo_tree: Option<&'a PhyloTree>,
    reference_genome: Option<&'a str>,
    nsc_result: &'a Option<nscluster::NsClusterResult>,
    output_dir: &'a Path,
}

fn nsc_presence_updates(
    pa: &crate::parser::PresenceAbsence,
    nsc_result: &Option<nscluster::NsClusterResult>,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
) -> ahash::HashSet<(usize, usize)> {
    nsc_result
        .as_ref()
        .map(|r| {
            let (_, updates) = nscluster::apply_ns_merges(
                &r.merge_map,
                |ci, gi| pa.is_present(ci, gi) || rescue_map.contains_key(&(ci, gi)),
                pa.n_genomes(),
            );
            updates
                .into_iter()
                .map(|(ci, gi)| (ci.val(), gi.val()))
                .collect()
        })
        .unwrap_or_default()
}

fn corrected_cluster_present(
    pa: &crate::parser::PresenceAbsence,
    ci: ClusterIdx,
    gi: GenomeIdx,
    matrix_rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    nsc_presence: &ahash::HashSet<(usize, usize)>,
) -> bool {
    pa.is_present(ci, gi)
        || matrix_rescue_map.contains_key(&(ci, gi))
        || nsc_presence.contains(&(ci.val(), gi.val()))
}

fn mobile_kind_from_labels(labels: &[&str]) -> Option<(&'static str, String)> {
    let joined = labels.join("|");
    let normalized = joined
        .to_ascii_lowercase()
        .replace("%20", " ")
        .replace('_', " ");
    let kind = if normalized.contains("is6110") {
        "IS6110"
    } else if normalized.contains("is1081") {
        "IS1081"
    } else if normalized.contains("insertion sequence")
        || normalized.contains("transposase")
        || normalized.contains("is element")
    {
        "insertion_sequence"
    } else if normalized.contains("integrase") || normalized.contains("phage") {
        "phage_or_integrase"
    } else if normalized.contains("recombinase") {
        "recombinase"
    } else {
        return None;
    };
    Some((kind, joined))
}

struct MobileElementSearchInput<'a> {
    path: &'a Path,
    genome_files: &'a HashMap<String, PathBuf>,
    pa: &'a crate::parser::PresenceAbsence,
    rescues: &'a [RescueEvent],
    synteny_rejected: &'a [RescueEvent],
    paralog_rejected: &'a [RescueEvent],
    threads: usize,
    min_identity: f64,
    min_coverage: f64,
    preset: &'a str,
}

fn search_mobile_element_reference(
    input: MobileElementSearchInput<'_>,
) -> Result<Vec<MobileElementHit>> {
    let mobile_ref = parse_pan_reference(input.path, true)
        .with_context(|| format!("Cannot parse mobile-element FASTA {}", input.path.display()))?;
    if mobile_ref.is_empty() {
        return Ok(Vec::new());
    }

    let mut mobile_names: Vec<String> = mobile_ref.keys().cloned().collect();
    mobile_names.sort();
    let mobile_indices: Vec<ClusterIdx> = (0..mobile_names.len()).map(ClusterIdx).collect();

    let mut genomes: ahash::HashSet<GenomeIdx> = ahash::HashSet::default();
    for rescue in input
        .rescues
        .iter()
        .chain(input.synteny_rejected.iter())
        .chain(input.paralog_rejected.iter())
    {
        genomes.insert(rescue.genome_idx);
    }
    if genomes.is_empty() {
        return Ok(Vec::new());
    }
    let mut searchable_genomes: Vec<(GenomeIdx, Vec<ClusterIdx>)> = genomes
        .into_iter()
        .map(|gi| (gi, mobile_indices.clone()))
        .collect();
    searchable_genomes.sort_by_key(|(gi, _)| gi.val());

    info!(
        "Running mobile-element sequence audit: {} queries × {} genomes",
        mobile_names.len(),
        searchable_genomes.len()
    );
    let hits = batch_search(
        &mobile_ref,
        input.genome_files,
        &searchable_genomes,
        input.pa.genome_names(),
        &mobile_names,
        input.threads,
        input.min_identity,
        input.preset,
        true,
        5000,
    )?;

    let mut mobile_hits = Vec::new();
    for ((gi, ci), cluster_hits) in hits {
        for hit in cluster_hits {
            if hit.query_coverage() < input.min_coverage {
                continue;
            }
            let identity = hit.identity();
            let coverage = hit.query_coverage();
            mobile_hits.push(MobileElementHit {
                genome_idx: gi,
                query_name: mobile_names[ci.val()].clone(),
                contig: hit.target,
                start: hit.target_start,
                end: hit.target_end,
                strand: hit.strand,
                identity,
                coverage,
            });
        }
    }
    mobile_hits.sort_by(|a, b| {
        a.genome_idx
            .cmp(&b.genome_idx)
            .then(a.contig.cmp(&b.contig))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
            .then(a.query_name.cmp(&b.query_name))
    });
    info!(
        "Mobile-element sequence audit retained {} hits",
        mobile_hits.len()
    );
    Ok(mobile_hits)
}

fn interval_distance(start_a: u64, end_a: u64, start_b: u64, end_b: u64) -> u64 {
    if end_a < start_b {
        start_b - end_a
    } else {
        start_a.saturating_sub(end_b)
    }
}

fn nearest_mobile_element_for_rescue(
    rescue: &RescueEvent,
    pa: &crate::parser::PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    mobile_hits: &[MobileElementHit],
    mobile_window: u64,
) -> Option<MobileElementEvidence> {
    let mut best: Option<MobileElementEvidence> = None;

    for hit in mobile_hits {
        if hit.genome_idx != rescue.genome_idx || hit.contig != rescue.contig {
            continue;
        }
        let distance = interval_distance(rescue.start, rescue.end, hit.start, hit.end);
        if distance > mobile_window {
            continue;
        }
        let candidate = MobileElementEvidence {
            genome_name: pa.genome_name(rescue.genome_idx).to_string(),
            source: "sequence",
            kind: "mobile_element_sequence",
            label: format!(
                "{}:{}-{}({};id={:.1};cov={:.1})",
                hit.query_name,
                hit.start + 1,
                hit.end,
                hit.strand,
                hit.identity,
                hit.coverage
            ),
            distance_bp: distance,
            identity: Some(hit.identity),
            coverage: Some(hit.coverage),
        };
        if best
            .as_ref()
            .is_none_or(|current| mobile_evidence_is_better(&candidate, current))
        {
            best = Some(candidate);
        }
    }

    let Some(gff) = gff_indices.get(&rescue.genome_idx) else {
        return best;
    };
    for gene in &gff.genes {
        if gene.contig != rescue.contig {
            continue;
        }
        let distance = interval_distance(rescue.start, rescue.end, gene.start, gene.end);
        if distance > mobile_window {
            continue;
        }
        let labels = gff.labels_for_gene(gene);
        let Some((kind, label)) = mobile_kind_from_labels(&labels) else {
            continue;
        };
        let candidate = MobileElementEvidence {
            genome_name: pa.genome_name(rescue.genome_idx).to_string(),
            source: "annotation",
            kind,
            label,
            distance_bp: distance,
            identity: None,
            coverage: None,
        };
        if best
            .as_ref()
            .is_none_or(|current| mobile_evidence_is_better(&candidate, current))
        {
            best = Some(candidate);
        }
    }
    best
}

fn nearest_mobile_element_for_cluster(
    rescues: &[&RescueEvent],
    pa: &crate::parser::PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    mobile_hits: &[MobileElementHit],
    mobile_window: u64,
) -> Option<MobileElementEvidence> {
    rescues
        .iter()
        .filter_map(|rescue| {
            nearest_mobile_element_for_rescue(rescue, pa, gff_indices, mobile_hits, mobile_window)
        })
        .min_by(|a, b| {
            a.distance_bp
                .cmp(&b.distance_bp)
                .then(mobile_source_priority(a.source).cmp(&mobile_source_priority(b.source)))
                .then(a.genome_name.cmp(&b.genome_name))
                .then(a.label.cmp(&b.label))
        })
}

fn mobile_source_priority(source: &str) -> u8 {
    if source == "sequence" {
        0
    } else {
        1
    }
}

fn mobile_evidence_is_better(
    candidate: &MobileElementEvidence,
    current: &MobileElementEvidence,
) -> bool {
    candidate
        .distance_bp
        .cmp(&current.distance_bp)
        .then(mobile_source_priority(candidate.source).cmp(&mobile_source_priority(current.source)))
        .then(candidate.label.cmp(&current.label))
        .is_lt()
}

fn nsc_role_for_cluster(ci: ClusterIdx, nsc_result: &Option<nscluster::NsClusterResult>) -> String {
    let Some(nsc) = nsc_result else {
        return String::new();
    };
    if nsc.merge_map.contains_key(&ci) {
        return "nsc_member_merged_into_representative".to_string();
    }
    for cluster in &nsc.merged_clusters {
        if cluster.members.contains(&ci) {
            return "nsc_representative".to_string();
        }
    }
    String::new()
}

fn rescue_type_counts(rescues: &[&RescueEvent]) -> String {
    if rescues.is_empty() {
        return String::new();
    }
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for rescue in rescues {
        *counts.entry(rescue.rescue_type.to_string()).or_default() += 1;
    }
    counts
        .into_iter()
        .map(|(kind, count)| format!("{}:{}", kind, count))
        .collect::<Vec<_>>()
        .join(";")
}

fn orf_signal_counts(rescues: &[&RescueEvent]) -> String {
    if rescues.is_empty() {
        return String::new();
    }
    let mut counts: BTreeMap<&'static str, usize> = BTreeMap::new();
    for rescue in rescues {
        *counts.entry(rescue_orf_signal(rescue)).or_default() += 1;
    }
    counts
        .into_iter()
        .map(|(signal, count)| format!("{}:{}", signal, count))
        .collect::<Vec<_>>()
        .join(";")
}

fn synteny_mean_score(rescues: &[&RescueEvent]) -> String {
    let scores: Vec<f64> = rescues
        .iter()
        .filter_map(|rescue| rescue.synteny_evidence.score)
        .collect();
    if scores.is_empty() {
        String::new()
    } else {
        format!("{:.3}", scores.iter().sum::<f64>() / scores.len() as f64)
    }
}

fn synteny_pass_fraction(rescues: &[&RescueEvent]) -> String {
    let assessed = rescues
        .iter()
        .filter(|rescue| rescue.synteny_ok.is_some())
        .count();
    if assessed == 0 {
        return String::new();
    }
    let passed = rescues
        .iter()
        .filter(|rescue| rescue.synteny_ok == Some(true))
        .count();
    format!("{:.3}", passed as f64 / assessed as f64)
}

fn paralog_locus_risk(
    context: Option<&ClusterContextSummary>,
    accepted: &[&RescueEvent],
    synteny_rejected: &[&RescueEvent],
    paralog_rejected: &[&RescueEvent],
) -> &'static str {
    let has_paralog = accepted
        .iter()
        .any(|rescue| matches!(rescue.rescue_type, RescueType::Paralog { .. }))
        || !paralog_rejected.is_empty();
    if has_paralog || context.is_some_and(|c| c.interpretation == "possible_multilocus_or_paralog")
    {
        "high"
    } else if !synteny_rejected.is_empty()
        || context.is_some_and(|c| c.interpretation == "mostly_stable_review_minor_contexts")
    {
        "moderate"
    } else if context.is_some_and(|c| {
        matches!(
            c.interpretation,
            "stable_single_locus" | "stable_single_locus_minor_noise"
        )
    }) {
        "low"
    } else {
        "unknown"
    }
}

fn phylogenetic_reliability(phylo: &PhyloPatternSummary) -> &'static str {
    if !phylo.assessed {
        return "not_assessed";
    }
    if matches!(
        phylo.pattern,
        "singleton_or_noise" | "homoplastic_gain_loss" | "dispersed_gain_loss"
    ) {
        return "low";
    }
    if phylo.parsimony_changes >= 6 {
        return "low";
    }
    if matches!(
        phylo.pattern,
        "lineage_specific_gain"
            | "clade_associated_presence"
            | "ancestral_retained"
            | "ancestral_loss"
    ) && phylo.parsimony_changes <= 2
    {
        return "high";
    }
    "moderate"
}

struct MechanismInput<'a> {
    pa: &'a crate::parser::PresenceAbsence,
    phylo: &'a PhyloPatternSummary,
    reference_state: &'a str,
    reference_event: Option<&'a ReferenceClusterEvent>,
    accepted: &'a [&'a RescueEvent],
    synteny_rejected: &'a [&'a RescueEvent],
    paralog_rejected: &'a [&'a RescueEvent],
    context: Option<&'a ClusterContextSummary>,
    mobile: Option<MobileElementEvidence>,
    nsc_role: &'a str,
}

fn classify_mechanism(input: MechanismInput<'_>) -> MechanismCall {
    let MechanismInput {
        pa,
        phylo,
        reference_state,
        reference_event,
        accepted,
        synteny_rejected,
        paralog_rejected,
        context,
        mobile,
        nsc_role,
    } = input;
    let has_pseudogene = accepted
        .iter()
        .any(|r| matches!(r.rescue_type, RescueType::Pseudogene { .. }));
    let has_fragment = accepted.iter().any(|r| {
        matches!(
            r.rescue_type,
            RescueType::Fragmented { .. } | RescueType::AlignmentSplit { .. }
        )
    });
    let has_split_or_boundary = accepted.iter().any(|r| {
        matches!(
            r.rescue_type,
            RescueType::Split { .. } | RescueType::StartVariant { .. }
        )
    });
    let has_paralog = accepted
        .iter()
        .any(|r| matches!(r.rescue_type, RescueType::Paralog { .. }))
        || !paralog_rejected.is_empty();
    let unstable_context =
        context.is_some_and(|c| c.interpretation == "possible_multilocus_or_paralog");
    let mobile_is = mobile.as_ref().is_some_and(|m| {
        (m.source == "sequence" || matches!(m.kind, "IS6110" | "IS1081" | "insertion_sequence"))
            && m.distance_bp <= 250
    });

    if let Some(event) = reference_event {
        if (has_paralog || unstable_context)
            && matches!(
                event.event_class,
                "partial_reference_homology" | "same_gene_full_length"
            )
        {
            return MechanismCall {
                mechanism: "duplication_or_paralog_expansion",
                detail: format!(
                    "reference_event={};target_reference_cluster={};paralog_or_multilocus_context",
                    event.event_class,
                    pa.cluster_name(event.target_cluster_idx)
                ),
                support: "low",
                wording: "Reference-like sequence is present in a paralogous or unstable context; interpret as copy-number/paralog or cluster-definition complexity before claiming gain/loss.",
                mobile,
            };
        }
        let mechanism = match event.event_class {
            "short_in_frame_allele" | "long_in_frame_allele" | "same_gene_full_length" => {
                "reference_gene_allelic_derivative"
            }
            "frameshift_or_truncated_reference_homology"
            | "frameshift_or_extension_reference_homology" => {
                "reference_gene_disruption_or_truncation"
            }
            _ => "partial_reference_homology",
        };
        return MechanismCall {
            mechanism,
            detail: format!(
                "{};target_reference_cluster={};target_reference_gene={}",
                event.event_class,
                pa.cluster_name(event.target_cluster_idx),
                event.target_gene_id
            ),
            support: if event.identity >= 95.0 { "high" } else { "moderate" },
            wording: "Keep as a separate matrix cluster, but describe as an event-derived reference-gene variant rather than a clean novel gene.",
            mobile,
        };
    }

    if mobile_is && (has_pseudogene || has_fragment) {
        return MechanismCall {
            mechanism: "IS_interruption_supported",
            detail: "fragmented_or_disrupted_rescue_with_nearby_insertion_sequence".to_string(),
            support: "high",
            wording: "Interpret as an IS-supported gene interruption/disruption candidate; keep matrix action governed by the orthology evidence.",
            mobile,
        };
    }

    if has_pseudogene {
        return MechanismCall {
            mechanism: "pseudogenization",
            detail: "pseudogene_rescue_evidence".to_string(),
            support: "moderate",
            wording: "Describe as disrupted homologous sequence, not as intact gene presence.",
            mobile,
        };
    }

    if has_fragment {
        return MechanismCall {
            mechanism: "gene_fragmentation_or_structural_break",
            detail: "fragmented_or_alignment_split_rescue_evidence".to_string(),
            support: "moderate",
            wording: "Describe as partial/fragmented homology and review the local structure before making a gene gain/loss claim.",
            mobile,
        };
    }

    if has_paralog || unstable_context {
        return MechanismCall {
            mechanism: "duplication_or_paralog_expansion",
            detail: "paralog_evidence_or_multilocus_context".to_string(),
            support: if has_paralog { "moderate" } else { "low" },
            wording: "Treat as a copy-number/paralog or multi-locus signal; avoid interpreting as a single clean gain/loss without locus review.",
            mobile,
        };
    }

    if !synteny_rejected.is_empty() {
        return MechanismCall {
            mechanism: "translocation_or_context_shift_candidate",
            detail: "sequence_evidence_rejected_by_synteny".to_string(),
            support: "low",
            wording: "Sequence similarity exists in an incompatible context; review as possible translocation, paralog, or annotation artefact.",
            mobile,
        };
    }

    if has_split_or_boundary || !nsc_role.is_empty() {
        return MechanismCall {
            mechanism: "annotation_or_clustering_artifact",
            detail: if nsc_role.is_empty() {
                "split_or_start_boundary_evidence".to_string()
            } else {
                nsc_role.to_string()
            },
            support: "moderate",
            wording: "Treat as a cluster split, boundary shift, or redundant annotation signal rather than an independent evolutionary gene event.",
            mobile,
        };
    }

    if reference_state == "absent"
        && matches!(
            phylo.pattern,
            "lineage_specific_gain"
                | "clade_associated_presence"
                | "single_gain_with_incomplete_clade_support"
        )
    {
        return MechanismCall {
            mechanism: "clean_insertion_or_acquisition_candidate",
            detail: "reference_absent_clade_associated_presence".to_string(),
            support: if phylo.carrier_monophyly { "moderate" } else { "low" },
            wording: "Candidate accessory acquisition/insertion relative to the audited reference; validate flanking sequence before claiming a novel locus.",
            mobile,
        };
    }

    if reference_state == "absent" && phylo.root_state == "present" {
        return MechanismCall {
            mechanism: "reference_or_lineage_loss_candidate",
            detail: "inferred_ancestral_presence_with_reference_absence".to_string(),
            support: "moderate",
            wording: "Candidate loss in the reference/background lineage rather than a de novo gain in all other genomes.",
            mobile,
        };
    }

    if matches!(phylo.pattern, "ancestral_loss" | "recurrent_loss") {
        return MechanismCall {
            mechanism: "deletion_loss_candidate",
            detail: phylo.pattern.to_string(),
            support: "moderate",
            wording: "Presence is inferred ancestral with one or more losses; inspect flanks for deletion scars or mobile-element mediation.",
            mobile,
        };
    }

    if phylo.pattern == "homoplastic_gain_loss" {
        return MechanismCall {
            mechanism: "unstable_or_homoplastic_accessory_signal",
            detail: "many_parsimony_changes".to_string(),
            support: "review",
            wording: "Highly recurrent presence/absence can reflect real mobility, recombination, or cluster/annotation instability; review before highlighting.",
            mobile,
        };
    }

    if phylo.pattern == "ancestral_retained" {
        return MechanismCall {
            mechanism: "ancestral_retained_no_disruption",
            detail: "present_across_matched_tree".to_string(),
            support: "high",
            wording: "Consistent with conserved ancestral presence in the analysed tree.",
            mobile,
        };
    }

    MechanismCall {
        mechanism: "unresolved_accessory_variation",
        detail: phylo.pattern.to_string(),
        support: "review",
        wording: "Evolutionary process is not resolved from matrix, synteny, reference-event, and optional tree evidence alone.",
        mobile,
    }
}

fn write_evolutionary_events(inputs: EvolutionaryEventInputs<'_>) -> Result<()> {
    let nsc_presence = nsc_presence_updates(inputs.pa, inputs.nsc_result, inputs.matrix_rescue_map);
    let reference_gi = inputs
        .reference_genome
        .and_then(|name| resolve_reference_genome(inputs.pa, name));

    let mut accepted_by_cluster: Vec<Vec<&RescueEvent>> = vec![Vec::new(); inputs.pa.n_clusters()];
    for rescue in inputs.rescues {
        accepted_by_cluster[rescue.cluster_idx.val()].push(rescue);
    }
    let mut synteny_by_cluster: Vec<Vec<&RescueEvent>> = vec![Vec::new(); inputs.pa.n_clusters()];
    for rescue in inputs.synteny_rejected {
        synteny_by_cluster[rescue.cluster_idx.val()].push(rescue);
    }
    let mut paralog_by_cluster: Vec<Vec<&RescueEvent>> = vec![Vec::new(); inputs.pa.n_clusters()];
    for rescue in inputs.paralog_rejected {
        paralog_by_cluster[rescue.cluster_idx.val()].push(rescue);
    }
    let mut reference_events_by_cluster: Vec<Option<&ReferenceClusterEvent>> =
        vec![None; inputs.pa.n_clusters()];
    for event in inputs.reference_cluster_events {
        reference_events_by_cluster[event.query_cluster_idx.val()] = Some(event);
    }
    let mut contexts_by_cluster: Vec<Option<&ClusterContextSummary>> =
        vec![None; inputs.pa.n_clusters()];
    for context in inputs.context_summaries {
        contexts_by_cluster[context.cluster_idx.val()] = Some(context);
    }

    let path = inputs.output_dir.join("evolutionary_events.tsv");
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "gene_cluster\tn_present_corrected\tprevalence_corrected\ttree_matched_tips\treference_state\tphylogenetic_pattern\tphylogenetic_reliability\tancestral_state\tphylogenetic_changes\tphylogenetic_gains\tphylogenetic_losses\tcarrier_monophyly\tmax_carrier_clade_fraction\tmolecular_mechanism\tmechanism_detail\tsupport_level\tis_interruption_supported\tmobile_element_source\tmobile_element_kind\tmobile_element_label\tmobile_element_genome\tmobile_element_distance_bp\tmobile_element_identity\tmobile_element_coverage\treference_event_class\ttarget_reference_cluster\ttarget_reference_gene\tcontext_interpretation\tparalog_locus_risk\tsynteny_mean_score\tsynteny_pass_fraction\torf_integrity_summary\trescue_evidence_counts\trecommended_wording"
    )?;

    for ci_val in 0..inputs.pa.n_clusters() {
        let ci = ClusterIdx(ci_val);
        let present_by_genome: Vec<bool> = (0..inputs.pa.n_genomes())
            .map(|gi_val| {
                corrected_cluster_present(
                    inputs.pa,
                    ci,
                    GenomeIdx(gi_val),
                    inputs.matrix_rescue_map,
                    &nsc_presence,
                )
            })
            .collect();
        let n_present = present_by_genome.iter().filter(|present| **present).count();
        let prevalence = if inputs.pa.n_genomes() > 0 {
            n_present as f64 / inputs.pa.n_genomes() as f64
        } else {
            0.0
        };
        let reference_state = reference_gi
            .map(|gi| {
                if corrected_cluster_present(
                    inputs.pa,
                    ci,
                    gi,
                    inputs.matrix_rescue_map,
                    &nsc_presence,
                ) {
                    "present"
                } else {
                    "absent"
                }
            })
            .unwrap_or("not_assessed");
        let phylo = summarize_phylogenetic_pattern(inputs.phylo_tree, &present_by_genome);
        let accepted = &accepted_by_cluster[ci_val];
        let synteny_rejected = &synteny_by_cluster[ci_val];
        let paralog_rejected = &paralog_by_cluster[ci_val];
        let all_rescue_evidence: Vec<&RescueEvent> = accepted
            .iter()
            .copied()
            .chain(synteny_rejected.iter().copied())
            .chain(paralog_rejected.iter().copied())
            .collect();
        let mobile = nearest_mobile_element_for_cluster(
            &all_rescue_evidence,
            inputs.pa,
            inputs.gff_indices,
            inputs.mobile_element_hits,
            inputs.mobile_element_window,
        );
        let nsc_role = nsc_role_for_cluster(ci, inputs.nsc_result);
        let mechanism = classify_mechanism(MechanismInput {
            pa: inputs.pa,
            phylo: &phylo,
            reference_state,
            reference_event: reference_events_by_cluster[ci_val],
            accepted,
            synteny_rejected,
            paralog_rejected,
            context: contexts_by_cluster[ci_val],
            mobile,
            nsc_role: &nsc_role,
        });
        let is_interruption = mechanism.mechanism == "IS_interruption_supported";
        let (
            mobile_source,
            mobile_kind,
            mobile_label,
            mobile_genome,
            mobile_distance,
            mobile_identity,
            mobile_coverage,
        ) = if let Some(mobile) = mechanism.mobile.as_ref() {
            (
                mobile.source.to_string(),
                mobile.kind.to_string(),
                mobile.label.clone(),
                mobile.genome_name.clone(),
                mobile.distance_bp.to_string(),
                mobile
                    .identity
                    .map(|value| format!("{:.1}", value))
                    .unwrap_or_default(),
                mobile
                    .coverage
                    .map(|value| format!("{:.1}", value))
                    .unwrap_or_default(),
            )
        } else {
            (
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
            )
        };
        let reference_event = reference_events_by_cluster[ci_val];
        let context_interpretation = contexts_by_cluster[ci_val]
            .map(|context| context.interpretation)
            .unwrap_or("");
        let paralog_risk = paralog_locus_risk(
            contexts_by_cluster[ci_val],
            accepted,
            synteny_rejected,
            paralog_rejected,
        );
        write_tsv_fields(
            &mut w,
            &[
                inputs.pa.cluster_name(ci).to_string(),
                n_present.to_string(),
                format!("{:.4}", prevalence),
                if phylo.assessed {
                    phylo.matched_tips.to_string()
                } else {
                    String::new()
                },
                reference_state.to_string(),
                phylo.pattern.to_string(),
                phylogenetic_reliability(&phylo).to_string(),
                phylo.root_state.to_string(),
                if phylo.assessed {
                    phylo.parsimony_changes.to_string()
                } else {
                    String::new()
                },
                if phylo.assessed {
                    phylo.gains.to_string()
                } else {
                    String::new()
                },
                if phylo.assessed {
                    phylo.losses.to_string()
                } else {
                    String::new()
                },
                if phylo.assessed {
                    phylo.carrier_monophyly.to_string()
                } else {
                    String::new()
                },
                if phylo.assessed {
                    format!("{:.3}", phylo.max_carrier_clade_fraction)
                } else {
                    String::new()
                },
                mechanism.mechanism.to_string(),
                mechanism.detail,
                mechanism.support.to_string(),
                is_interruption.to_string(),
                mobile_source,
                mobile_kind,
                mobile_label,
                mobile_genome,
                mobile_distance,
                mobile_identity,
                mobile_coverage,
                reference_event
                    .map(|event| event.event_class.to_string())
                    .unwrap_or_default(),
                reference_event
                    .map(|event| inputs.pa.cluster_name(event.target_cluster_idx).to_string())
                    .unwrap_or_default(),
                reference_event
                    .map(|event| event.target_gene_id.clone())
                    .unwrap_or_default(),
                context_interpretation.to_string(),
                paralog_risk.to_string(),
                synteny_mean_score(&all_rescue_evidence),
                synteny_pass_fraction(&all_rescue_evidence),
                orf_signal_counts(&all_rescue_evidence),
                rescue_type_counts(accepted),
                mechanism.wording.to_string(),
            ],
        )?;
    }
    w.flush()?;
    info!("Evolutionary event audit written to {}", path.display());
    Ok(())
}

#[derive(Debug, Clone, Default)]
struct ContextAccumulator {
    count: usize,
    examples: Vec<String>,
}

#[derive(Debug, Clone)]
struct ClusterContextSummary {
    cluster_idx: ClusterIdx,
    n_present: usize,
    n_mapped: usize,
    n_unresolved: usize,
    n_unique_contexts: usize,
    dominant_count: usize,
    dominant_fraction: f64,
    interpretation: &'static str,
}

fn context_signature(
    gff: &GffIndex,
    genes: &[&GffGene],
    gene_to_cluster: &HashMap<String, ClusterIdx>,
    pa: &crate::parser::PresenceAbsence,
) -> String {
    let names: Vec<&str> = genes
        .iter()
        .filter_map(|gene| gff.cluster_for_gene(gene, gene_to_cluster))
        .map(|ci| pa.cluster_name(ci))
        .collect();
    if names.is_empty() {
        "NA".to_string()
    } else {
        names.join("|")
    }
}

fn context_interpretation(
    n_mapped: usize,
    n_unique_contexts: usize,
    dominant_fraction: f64,
) -> &'static str {
    if n_mapped < 3 {
        "insufficient_context"
    } else if n_unique_contexts <= 1 {
        "stable_single_locus"
    } else if dominant_fraction >= 0.95 {
        "stable_single_locus_minor_noise"
    } else if dominant_fraction >= 0.80 {
        "mostly_stable_review_minor_contexts"
    } else {
        "possible_multilocus_or_paralog"
    }
}

fn write_cluster_locus_contexts(
    pa: &crate::parser::PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    gene_to_cluster_maps: &HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
    output_dir: &std::path::Path,
) -> Result<(usize, Vec<ClusterContextSummary>)> {
    let mut context_rows: BTreeMap<(usize, String, String), ContextAccumulator> = BTreeMap::new();
    let mut summaries = Vec::new();

    for ci_val in 0..pa.n_clusters() {
        let ci = ClusterIdx(ci_val);
        let mut n_mapped = 0usize;
        let mut n_unresolved = 0usize;
        let mut context_counts: HashMap<(String, String), usize> = HashMap::new();

        let mut genomes: Vec<GenomeIdx> = pa.present_genomes(ci).iter().copied().collect();
        genomes.sort_by_key(|gi| gi.val());

        for gi in genomes {
            let Some(gff) = gff_indices.get(&gi) else {
                continue;
            };
            let Some(g2c) = gene_to_cluster_maps.get(&gi) else {
                continue;
            };
            let cell = pa.cell(ci, gi);
            let Some(gene) = gff_gene_for_cell(gff, cell) else {
                n_unresolved += 1;
                continue;
            };

            let (upstream, downstream) = gff.context_window(
                &gene.contig,
                gene.start,
                gene.end,
                Some(gene.gene_id.as_str()),
                2,
            );
            let upstream_sig = context_signature(gff, &upstream, g2c, pa);
            let downstream_sig = context_signature(gff, &downstream, g2c, pa);
            let key = (upstream_sig.clone(), downstream_sig.clone());
            *context_counts.entry(key).or_insert(0) += 1;
            n_mapped += 1;

            let row_key = (ci.val(), upstream_sig, downstream_sig);
            let row = context_rows.entry(row_key).or_default();
            row.count += 1;
            if row.examples.len() < 5 {
                row.examples
                    .push(format!("{}:{}", pa.genome_name(gi), gene.gene_id));
            }
        }

        let n_unique_contexts = context_counts.len();
        let dominant_count = context_counts.values().copied().max().unwrap_or(0);
        let dominant_fraction = if n_mapped > 0 {
            dominant_count as f64 / n_mapped as f64
        } else {
            0.0
        };
        let interpretation = context_interpretation(n_mapped, n_unique_contexts, dominant_fraction);
        summaries.push(ClusterContextSummary {
            cluster_idx: ci,
            n_present: pa.present_count(ci),
            n_mapped,
            n_unresolved,
            n_unique_contexts,
            dominant_count,
            dominant_fraction,
            interpretation,
        });
    }

    let context_path = output_dir.join("cluster_locus_context.tsv");
    let context_file = fs::File::create(&context_path)?;
    let mut context_w = BufWriter::new(context_file);
    writeln!(
        context_w,
        "gene_cluster\tcontext_rank\tn_genomes\tfraction_of_mapped_contexts\tupstream_context\tdownstream_context\texample_genomes"
    )?;
    let mut written_context_rows = 0usize;
    let mut rows_by_cluster: BTreeMap<usize, Vec<(String, String, ContextAccumulator)>> =
        BTreeMap::new();
    for ((ci, upstream, downstream), acc) in context_rows {
        rows_by_cluster
            .entry(ci)
            .or_default()
            .push((upstream, downstream, acc));
    }
    for (ci_val, mut rows) in rows_by_cluster {
        rows.sort_by(|a, b| {
            b.2.count
                .cmp(&a.2.count)
                .then(a.0.cmp(&b.0))
                .then(a.1.cmp(&b.1))
        });
        let mapped = summaries[ci_val].n_mapped.max(1);
        for (rank, (upstream, downstream, acc)) in rows.into_iter().enumerate() {
            write_tsv_fields(
                &mut context_w,
                &[
                    pa.cluster_name(ClusterIdx(ci_val)).to_string(),
                    (rank + 1).to_string(),
                    acc.count.to_string(),
                    format!("{:.3}", acc.count as f64 / mapped as f64),
                    upstream,
                    downstream,
                    acc.examples.join(";"),
                ],
            )?;
            written_context_rows += 1;
        }
    }
    context_w.flush()?;

    let summary_path = output_dir.join("cluster_locus_summary.tsv");
    let summary_file = fs::File::create(&summary_path)?;
    let mut summary_w = BufWriter::new(summary_file);
    writeln!(
        summary_w,
        "gene_cluster\tn_present\tn_context_mapped\tn_unresolved_gene_ids\tn_unique_contexts\tdominant_context_count\tdominant_context_fraction\tinterpretation"
    )?;
    for summary in &summaries {
        write_tsv_fields(
            &mut summary_w,
            &[
                pa.cluster_name(summary.cluster_idx).to_string(),
                summary.n_present.to_string(),
                summary.n_mapped.to_string(),
                summary.n_unresolved.to_string(),
                summary.n_unique_contexts.to_string(),
                summary.dominant_count.to_string(),
                format!("{:.3}", summary.dominant_fraction),
                summary.interpretation.to_string(),
            ],
        )?;
    }
    summary_w.flush()?;

    Ok((written_context_rows, summaries))
}

fn write_corrected_csv(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    output_dir: &std::path::Path,
    nsc_result: &Option<nscluster::NsClusterResult>,
) -> Result<()> {
    // Build set of clusters to skip (merged into another by NSC)
    let clusters_to_skip: ahash::HashSet<usize> = nsc_result
        .as_ref()
        .map(|r| r.merge_map.keys().map(|ci| ci.val()).collect())
        .unwrap_or_default();

    // Build NSC presence updates: (rep_ci, gi) pairs to mark present
    let nsc_presence = nsc_presence_updates(pa, nsc_result, rescue_map);
    let path = output_dir.join("gene_presence_absence_recore.csv");
    let file = fs::File::create(&path)?;
    let writer = BufWriter::with_capacity(256 * 1024, file);
    let mut csv_writer = csv::Writer::from_writer(writer);

    csv_writer.write_record(pa.header())?;

    let genome_idx_set: ahash::HashSet<usize> = pa.genome_col_indices().iter().copied().collect();
    let meta_col_indices: Vec<usize> = (0..pa.header().len())
        .filter(|i| !genome_idx_set.contains(i))
        .collect();

    let isolates_col = pa.header().iter().position(|h| h == "No. isolates");
    let sequences_col = pa.header().iter().position(|h| h == "No. sequences");
    let avg_seq_col = pa
        .header()
        .iter()
        .position(|h| h == "Avg sequences per isolate");

    let n_cols = pa.header().len();
    let mut row: Vec<String> = vec![String::new(); n_cols];
    let mut rescue_buf = String::new();

    for ci_val in 0..pa.n_clusters() {
        // Skip clusters merged by NSC (their presence is merged into representative)
        if clusters_to_skip.contains(&ci_val) {
            continue;
        }

        let ci = ClusterIdx(ci_val);

        for (mi, &col_idx) in meta_col_indices.iter().enumerate() {
            let meta = pa.metadata(ci);
            if mi < meta.len() {
                row[col_idx].clone_from(&meta[mi]);
            } else {
                row[col_idx].clear();
            }
        }

        let mut rescue_count = 0usize;
        let mut nsc_count = 0usize;
        for (gi_val, &col_idx) in pa.genome_col_indices().iter().enumerate() {
            let gi = GenomeIdx(gi_val);
            if let Some(r) = rescue_map.get(&(ci, gi)) {
                rescue_buf.clear();
                use std::fmt::Write as FmtWrite;
                // 1-based coords for the corrected CSV (biologist-friendly)
                write!(rescue_buf, "RESCUED:{}:{}-{}", r.contig, r.start + 1, r.end)
                    .expect("write to String is infallible");
                row[col_idx].clone_from(&rescue_buf);
                rescue_count += 1;
            } else if nsc_presence.contains(&(ci_val, gi_val)) {
                // Gene absent in this cluster but present in a merged NSC member
                row[col_idx].clear();
                row[col_idx].push_str("NSC_MERGED");
                nsc_count += 1;
            } else {
                row[col_idx].clear();
                row[col_idx].push_str(pa.cell(ci, gi));
            }
        }

        let total_new = rescue_count + nsc_count;
        if total_new > 0 {
            let new_count = pa.present_count(ci) + total_new;
            use std::fmt::Write as FmtWrite;
            if let Some(col) = isolates_col {
                row[col].clear();
                write!(row[col], "{}", new_count).expect("write to String is infallible");
            }
            if let Some(col) = sequences_col {
                row[col].clear();
                write!(row[col], "{}", new_count).expect("write to String is infallible");
            }
            if let Some(col) = avg_seq_col {
                row[col].clear();
                write!(row[col], "1.00").expect("write to String is infallible");
            }
        }

        csv_writer.write_record(row.iter().map(|s| s.as_str()))?;
    }
    csv_writer.flush()?;
    Ok(())
}

struct MultistateMatrixInputs<'a> {
    pa: &'a crate::parser::PresenceAbsence,
    accepted_rescue_map: &'a HashMap<(ClusterIdx, GenomeIdx), &'a RescueEvent>,
    nsc_rescue_map: &'a HashMap<(ClusterIdx, GenomeIdx), &'a RescueEvent>,
    synteny_rejected: &'a [RescueEvent],
    paralog_rejected: &'a [RescueEvent],
    gene_to_cluster_maps: &'a HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
    output_dir: &'a std::path::Path,
    nsc_result: &'a Option<nscluster::NsClusterResult>,
}

fn write_multistate_matrix(inputs: MultistateMatrixInputs<'_>) -> Result<()> {
    let pa = inputs.pa;
    let clusters_to_skip: ahash::HashSet<usize> = inputs
        .nsc_result
        .as_ref()
        .map(|r| r.merge_map.keys().map(|ci| ci.val()).collect())
        .unwrap_or_default();

    let nsc_presence = nsc_presence_updates(pa, inputs.nsc_result, inputs.nsc_rescue_map);

    let synteny_rejected_map: HashMap<(ClusterIdx, GenomeIdx), &RescueEvent> = inputs
        .synteny_rejected
        .iter()
        .map(|r| ((r.cluster_idx, r.genome_idx), r))
        .collect();
    let paralog_rejected_map: HashMap<(ClusterIdx, GenomeIdx), &RescueEvent> = inputs
        .paralog_rejected
        .iter()
        .map(|r| ((r.cluster_idx, r.genome_idx), r))
        .collect();

    let path = inputs
        .output_dir
        .join("gene_presence_absence_recore_states.csv");
    let file = fs::File::create(&path)?;
    let writer = BufWriter::with_capacity(256 * 1024, file);
    let mut csv_writer = csv::Writer::from_writer(writer);
    csv_writer.write_record(pa.header())?;

    let genome_idx_set: ahash::HashSet<usize> = pa.genome_col_indices().iter().copied().collect();
    let meta_col_indices: Vec<usize> = (0..pa.header().len())
        .filter(|i| !genome_idx_set.contains(i))
        .collect();

    let n_cols = pa.header().len();
    let mut row: Vec<String> = vec![String::new(); n_cols];

    for ci_val in 0..pa.n_clusters() {
        if clusters_to_skip.contains(&ci_val) {
            continue;
        }

        let ci = ClusterIdx(ci_val);
        for (mi, &col_idx) in meta_col_indices.iter().enumerate() {
            let meta = pa.metadata(ci);
            if mi < meta.len() {
                row[col_idx].clone_from(&meta[mi]);
            } else {
                row[col_idx].clear();
            }
        }

        for (gi_val, &col_idx) in pa.genome_col_indices().iter().enumerate() {
            let gi = GenomeIdx(gi_val);
            row[col_idx].clear();
            if pa.is_present(ci, gi) {
                row[col_idx].push_str("PRESENT_INPUT");
            } else if let Some(rescue) = inputs.accepted_rescue_map.get(&(ci, gi)) {
                row[col_idx].push_str(accepted_rescue_evidence_state(
                    rescue,
                    pa,
                    inputs.gene_to_cluster_maps,
                ));
            } else if let Some(rescue) = synteny_rejected_map.get(&(ci, gi)) {
                let query_cluster = pa.cluster_name(rescue.cluster_idx);
                let annotation_cluster =
                    annotation_cluster_name(rescue, pa, inputs.gene_to_cluster_maps);
                row[col_idx].push_str(evidence_state(
                    rescue,
                    false,
                    query_cluster,
                    &annotation_cluster,
                    "synteny_fail",
                ));
            } else if let Some(rescue) = paralog_rejected_map.get(&(ci, gi)) {
                let query_cluster = pa.cluster_name(rescue.cluster_idx);
                let annotation_cluster =
                    annotation_cluster_name(rescue, pa, inputs.gene_to_cluster_maps);
                row[col_idx].push_str(evidence_state(
                    rescue,
                    false,
                    query_cluster,
                    &annotation_cluster,
                    "paralog_without_synteny",
                ));
            } else if nsc_presence.contains(&(ci_val, gi_val)) {
                row[col_idx].push_str("ALIAS_CLUSTER");
            } else {
                row[col_idx].push_str("ABSENT");
            }
        }

        csv_writer.write_record(row.iter().map(|s| s.as_str()))?;
    }
    csv_writer.flush()?;
    Ok(())
}

fn copy_text_maybe_gzip_with_header(path: &std::path::Path, writer: &mut impl Write) -> Result<()> {
    let file = fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    let mut reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(flate2::read::GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut prefix = Vec::with_capacity(256);
    let mut line = Vec::with_capacity(4096);
    let mut has_gff_header = false;
    let mut ended_newline = true;

    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        ended_newline = line.ends_with(b"\n");

        let text = String::from_utf8_lossy(&line);
        let trimmed = text.trim_start();
        if trimmed.is_empty() {
            prefix.extend_from_slice(&line);
            continue;
        }
        has_gff_header = trimmed.starts_with("##gff-version");

        prefix.extend_from_slice(&line);
        break;
    }

    if !has_gff_header {
        writeln!(writer, "##gff-version 3")?;
    }
    writer.write_all(&prefix)?;

    let mut buf = Vec::with_capacity(64 * 1024);
    loop {
        buf.clear();
        let n = reader.read_until(b'\n', &mut buf)?;
        if n == 0 {
            break;
        }
        ended_newline = buf.ends_with(b"\n");
        writer.write_all(&buf)?;
    }

    if !ended_newline {
        writeln!(writer)?;
    }
    Ok(())
}

fn gff_escape(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '%' => out.push_str("%25"),
            ';' => out.push_str("%3B"),
            '=' => out.push_str("%3D"),
            ',' => out.push_str("%2C"),
            '\t' => out.push_str("%09"),
            '\n' | '\r' => out.push(' '),
            ' ' => out.push_str("%20"),
            _ => out.push(ch),
        }
    }
    out
}

fn safe_id_component(value: &str) -> String {
    value
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.' | ':') {
                ch
            } else {
                '_'
            }
        })
        .collect()
}

fn write_corrected_gffs(
    rescues: &[RescueEvent],
    pa: &crate::parser::PresenceAbsence,
    gff_paths: &HashMap<GenomeIdx, PathBuf>,
    output_dir: &std::path::Path,
) -> Result<()> {
    let gff_dir = output_dir.join("corrected_gff");
    fs::create_dir_all(&gff_dir)?;

    let mut by_genome: HashMap<GenomeIdx, Vec<&RescueEvent>> = HashMap::new();
    for rescue in rescues {
        by_genome.entry(rescue.genome_idx).or_default().push(rescue);
    }

    let mut genome_indices: Vec<GenomeIdx> = by_genome.keys().copied().collect();
    genome_indices.sort();

    for gi in genome_indices {
        let genome_name = pa.genome_name(gi);
        let path = gff_dir.join(format!("{}.recore.gff3", genome_name));
        let file = fs::File::create(&path)?;
        let mut w = BufWriter::new(file);

        if let Some(source_gff) = gff_paths.get(&gi) {
            copy_text_maybe_gzip_with_header(source_gff, &mut w)?;
        } else {
            writeln!(w, "##gff-version 3")?;
        }

        writeln!(w, "###")?;
        writeln!(w, "# Recore rescue features appended below")?;

        let mut rows = by_genome.remove(&gi).unwrap_or_default();
        rows.sort_by(|a, b| {
            a.contig
                .cmp(&b.contig)
                .then(a.start.cmp(&b.start))
                .then(a.end.cmp(&b.end))
                .then(a.cluster_idx.cmp(&b.cluster_idx))
        });

        for rescue in rows {
            let cluster = pa.cluster_name(rescue.cluster_idx);
            let feature_id = format!(
                "recore_{}_{}_{}_{}",
                safe_id_component(cluster),
                safe_id_component(&rescue.contig),
                rescue.start + 1,
                rescue.end
            );
            let synteny = match rescue.synteny_ok {
                Some(true) => "PASS",
                Some(false) => "FAIL",
                None => "NA",
            };
            let mut attrs = vec![
                format!("ID={}", gff_escape(&feature_id)),
                format!("gene_cluster={}", gff_escape(cluster)),
                format!(
                    "rescue_type={}",
                    gff_escape(&rescue.rescue_type.to_string())
                ),
                format!("identity={:.1}", rescue.identity),
                format!("coverage={:.1}", rescue.coverage),
                format!("confidence={:.3}", rescue.confidence),
                format!("synteny={}", synteny),
            ];
            if let Some(annotation_match) = rescue.annotation_match.as_deref() {
                attrs.push(format!("annotation_match={}", gff_escape(annotation_match)));
            }
            if rescue.fragment_evidence.n_fragments > 0 {
                attrs.push(format!(
                    "fragment_n={}",
                    rescue.fragment_evidence.n_fragments
                ));
                attrs.push(format!(
                    "fragment_query_gap_bp={}",
                    rescue.fragment_evidence.total_query_gap_bp
                ));
                attrs.push(format!(
                    "fragment_target_gap_bp={}",
                    rescue.fragment_evidence.total_target_gap_bp
                ));
                if let Some(note) = rescue.fragment_evidence.note.as_deref() {
                    attrs.push(format!("fragment_note={}", gff_escape(note)));
                }
            }
            attrs.push("Note=Recore_rescue_feature".to_string());

            writeln!(
                w,
                "{}\tRecore\tCDS\t{}\t{}\t.\t{}\t0\t{}",
                rescue.contig,
                rescue.start + 1,
                rescue.end,
                rescue.strand,
                attrs.join(";")
            )?;
        }
        w.flush()?;
    }

    info!("Corrected GFF3 files written to {}", gff_dir.display());
    Ok(())
}

fn count_core_after(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    nsc_result: &Option<nscluster::NsClusterResult>,
) -> usize {
    let n = pa.n_genomes();
    let nsc_presence = nsc_presence_updates(pa, nsc_result, rescue_map);
    let clusters_to_skip: ahash::HashSet<usize> = nsc_result
        .as_ref()
        .map(|r| r.merge_map.keys().map(|ci| ci.val()).collect())
        .unwrap_or_default();

    (0..pa.n_clusters())
        .filter(|&ci_val| {
            if clusters_to_skip.contains(&ci_val) {
                return false;
            }
            let ci = ClusterIdx(ci_val);
            let present = pa.present_count(ci);
            if present == n {
                return true;
            }
            let extra = (0..n)
                .filter(|&gi_val| {
                    let gi = GenomeIdx(gi_val);
                    !pa.is_present(ci, gi)
                        && (rescue_map.contains_key(&(ci, gi))
                            || nsc_presence.contains(&(ci_val, gi_val)))
                })
                .count();
            present + extra == n
        })
        .count()
}

fn count_soft_core_after(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    nsc_result: &Option<nscluster::NsClusterResult>,
) -> usize {
    let n = pa.n_genomes();
    let threshold = (n as f64 * 0.95).ceil() as usize;
    let nsc_presence = nsc_presence_updates(pa, nsc_result, rescue_map);
    let clusters_to_skip: ahash::HashSet<usize> = nsc_result
        .as_ref()
        .map(|r| r.merge_map.keys().map(|ci| ci.val()).collect())
        .unwrap_or_default();

    (0..pa.n_clusters())
        .filter(|&ci_val| {
            if clusters_to_skip.contains(&ci_val) {
                return false;
            }
            let ci = ClusterIdx(ci_val);
            let present = pa.present_count(ci);
            if present >= threshold {
                return true;
            }
            let extra = (0..n)
                .filter(|&gi_val| {
                    let gi = GenomeIdx(gi_val);
                    !pa.is_present(ci, gi)
                        && (rescue_map.contains_key(&(ci, gi))
                            || nsc_presence.contains(&(ci_val, gi_val)))
                })
                .count();
            present + extra >= threshold
        })
        .count()
}

struct SummaryStats {
    core_before: usize,
    total_gaps: usize,
    filtered_gaps: usize,
    skipped_by_prevalence: usize,
    searchable_gaps: usize,
    missing_fasta_gaps: usize,
    dropped: usize,
    missing_cds: usize,
    split: usize,
    fragmented: usize,
    alignment_split: usize,
    start_variant: usize,
    paralog: usize,
    pseudogene: usize,
    synteny_rejected: usize,
    paralog_rejected: usize,
    matrix_policy: MatrixPolicy,
    matrix_corrected: usize,
    nsc_matrix_corrected: usize,
    synteny_enabled: bool,
    nsc_clusters: usize,
    nsc_merged_genes: usize,
}

fn write_nsc_report(
    nsc: &nscluster::NsClusterResult,
    cluster_names: &[String],
    output_dir: &std::path::Path,
) -> Result<()> {
    let path = output_dir.join("nscluster_report.tsv");
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);
    writeln!(w, "ns_cluster\trepresentative\tn_members\tmembers")?;
    for (i, nsc_cl) in nsc.merged_clusters.iter().enumerate() {
        let member_names: Vec<&str> = nsc_cl
            .members
            .iter()
            .map(|ci| cluster_names[ci.val()].as_str())
            .collect();
        write_tsv_fields(
            &mut w,
            &[
                format!("NS_{}", i),
                nsc_cl.representative.clone(),
                nsc_cl.members.len().to_string(),
                member_names.join(";"),
            ],
        )?;
    }
    w.flush()?;
    Ok(())
}

fn write_summary(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    rescues: &[RescueEvent],
    stats: &SummaryStats,
    nsc_result: &Option<nscluster::NsClusterResult>,
    output_dir: &std::path::Path,
) -> Result<()> {
    let core_after = count_core_after(pa, rescue_map, nsc_result);
    let soft_core_before = pa.soft_core_size();
    let soft_core_after = count_soft_core_after(pa, rescue_map, nsc_result);

    let path = output_dir.join("summary.txt");
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);
    writeln!(w, "Recore Summary")?;
    writeln!(w, "{}", "=".repeat(50))?;
    writeln!(w, "Clusters:              {}", pa.n_clusters())?;
    writeln!(w, "Genomes:               {}", pa.n_genomes())?;
    writeln!(w, "Total gaps:            {}", stats.total_gaps)?;
    writeln!(w, "Gaps passing filter:   {}", stats.filtered_gaps)?;
    writeln!(w, "Gaps searched:         {}", stats.searchable_gaps)?;
    if stats.missing_fasta_gaps > 0 {
        writeln!(w, "Gaps skipped (FASTA):  {}", stats.missing_fasta_gaps)?;
    }
    writeln!(w, "Skipped (prevalence):  {}", stats.skipped_by_prevalence)?;
    writeln!(w, "Rescued:               {}", rescues.len())?;
    writeln!(w, "Matrix policy:         {:?}", stats.matrix_policy)?;
    let total_matrix_corrected = stats.matrix_corrected + stats.nsc_matrix_corrected;
    writeln!(w, "Matrix-corrected:      {}", total_matrix_corrected)?;
    writeln!(w, "  rescue cells:        {}", stats.matrix_corrected)?;
    if stats.nsc_matrix_corrected > 0 {
        writeln!(w, "  NSC merged cells:    {}", stats.nsc_matrix_corrected)?;
    }
    if stats.synteny_rejected > 0 {
        writeln!(w, "Rejected (synteny):    {}", stats.synteny_rejected)?;
    }
    if stats.paralog_rejected > 0 {
        writeln!(w, "Rejected (paralog):    {}", stats.paralog_rejected)?;
    }
    writeln!(
        w,
        "Remaining gaps:        {}",
        stats.searchable_gaps.saturating_sub(total_matrix_corrected)
    )?;
    writeln!(w)?;
    let core_diff = core_after as i64 - stats.core_before as i64;
    let soft_diff = soft_core_after as i64 - soft_core_before as i64;
    writeln!(
        w,
        "Core genome (strict):  {} → {} ({:+})",
        stats.core_before, core_after, core_diff
    )?;
    writeln!(
        w,
        "Soft-core (≥95%):      {} → {} ({:+})",
        soft_core_before, soft_core_after, soft_diff
    )?;
    writeln!(w)?;
    writeln!(w, "Rescue types:")?;
    writeln!(w, "  DROPPED:       {}", stats.dropped)?;
    writeln!(w, "  MISSING_CDS:   {}", stats.missing_cds)?;
    writeln!(w, "  SPLIT:         {}", stats.split)?;
    writeln!(w, "  FRAGMENTED:    {}", stats.fragmented)?;
    writeln!(w, "  ALIGNMENT_SPLIT: {}", stats.alignment_split)?;
    writeln!(w, "  START_VARIANT: {}", stats.start_variant)?;
    writeln!(w, "  PARALOG:       {}", stats.paralog)?;
    writeln!(w, "  PSEUDOGENE:    {}", stats.pseudogene)?;
    if stats.synteny_enabled {
        let syn_pass = rescues
            .iter()
            .filter(|r| r.synteny_ok == Some(true))
            .count();
        let syn_fail = rescues
            .iter()
            .filter(|r| r.synteny_ok == Some(false))
            .count();
        writeln!(w)?;
        writeln!(w, "Synteny validation:")?;
        writeln!(w, "  PASS:  {}", syn_pass)?;
        writeln!(w, "  FAIL:  {}", syn_fail)?;
        if stats.synteny_rejected > 0 {
            writeln!(w, "  REJECTED_FAIL: {}", stats.synteny_rejected)?;
        }
    }
    if stats.paralog_rejected > 0 {
        writeln!(w)?;
        writeln!(w, "Orthology rejection:")?;
        writeln!(w, "  PARALOG_WITHOUT_SYNTENY: {}", stats.paralog_rejected)?;
    }
    if stats.nsc_clusters > 0 {
        writeln!(w)?;
        writeln!(w, "Nucleotide similarity clustering:")?;
        writeln!(w, "  NS clusters:     {}", stats.nsc_clusters)?;
        writeln!(w, "  Genes merged:    {}", stats.nsc_merged_genes)?;
        writeln!(
            w,
            "  Clusters removed: {}",
            stats.nsc_merged_genes - stats.nsc_clusters
        )?;
    }
    w.flush()?;
    Ok(())
}

fn json_escape(value: &str) -> String {
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\0' => out.push_str("\\u0000"),
            c if c.is_control() => {
                use std::fmt::Write as _;
                write!(out, "\\u{:04x}", c as u32).expect("write to String is infallible");
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn json_path(path: &std::path::Path) -> String {
    json_escape(&path.display().to_string())
}

fn write_run_manifest(
    args: &Args,
    pa: &crate::parser::PresenceAbsence,
    stats: &SummaryStats,
    output_dir: &std::path::Path,
) -> Result<()> {
    let path = output_dir.join("run_manifest.json");
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);
    let total_matrix_corrected = stats.matrix_corrected + stats.nsc_matrix_corrected;
    let rescued = stats.dropped
        + stats.missing_cds
        + stats.split
        + stats.fragmented
        + stats.alignment_split
        + stats.start_variant
        + stats.paralog
        + stats.pseudogene;

    writeln!(w, "{{")?;
    writeln!(w, "  \"tool\": \"recore\",")?;
    writeln!(
        w,
        "  \"tool_version\": {},",
        json_escape(env!("CARGO_PKG_VERSION"))
    )?;
    writeln!(w, "  \"schema_version\": \"recore.outputs.v2\",")?;
    writeln!(w, "  \"inputs\": {{")?;
    writeln!(
        w,
        "    \"presence_absence\": {},",
        json_path(&args.presence_absence)
    )?;
    writeln!(
        w,
        "    \"pan_reference\": {},",
        json_path(&args.pan_reference)
    )?;
    writeln!(w, "    \"genomes\": {},", json_path(&args.genomes))?;
    if let Some(path) = args.annotations.as_deref() {
        writeln!(w, "    \"annotations\": {},", json_path(path))?;
    } else {
        writeln!(w, "    \"annotations\": null,")?;
    }
    if let Some(tree) = args.tree.as_deref() {
        writeln!(w, "    \"tree\": {},", json_path(tree))?;
    } else {
        writeln!(w, "    \"tree\": null,")?;
    }
    if let Some(mobile_ref) = args.mobile_element_reference.as_deref() {
        writeln!(
            w,
            "    \"mobile_element_reference\": {}",
            json_path(mobile_ref)
        )?;
    } else {
        writeln!(w, "    \"mobile_element_reference\": null")?;
    }
    writeln!(w, "  }},")?;

    writeln!(w, "  \"parameters\": {{")?;
    writeln!(w, "    \"min_identity\": {:.3},", args.min_identity)?;
    writeln!(w, "    \"min_coverage\": {:.3},", args.min_coverage)?;
    writeln!(w, "    \"min_prevalence\": {:.6},", args.min_prevalence)?;
    writeln!(w, "    \"preset\": {},", json_escape(&args.preset))?;
    writeln!(w, "    \"threads\": {},", args.threads)?;
    writeln!(
        w,
        "    \"matrix_policy\": {},",
        json_escape(&format!("{:?}", args.matrix_policy))
    )?;
    writeln!(w, "    \"multi_allele\": {},", args.multi_allele)?;
    writeln!(w, "    \"synteny\": {},", args.synteny)?;
    writeln!(
        w,
        "    \"write_corrected_gff\": {},",
        args.write_corrected_gff
    )?;
    if let Some(reference_genome) = args.reference_genome.as_deref() {
        writeln!(
            w,
            "    \"reference_genome\": {},",
            json_escape(reference_genome)
        )?;
    } else {
        writeln!(w, "    \"reference_genome\": null,")?;
    }
    writeln!(w, "    \"nscluster\": {},", args.nscluster)?;
    writeln!(w, "    \"nsc_engine\": {},", json_escape(&args.nsc_engine))?;
    writeln!(w, "    \"nsc_kmer_size\": {},", args.nsc_kmer_size)?;
    writeln!(w, "    \"nsc_threshold\": {:.6},", args.nsc_threshold)?;
    writeln!(w, "    \"nsc_min_seq_id\": {:.6},", args.nsc_min_seq_id)?;
    writeln!(w, "    \"nsc_coverage\": {:.6},", args.nsc_coverage)?;
    writeln!(
        w,
        "    \"mobile_element_min_identity\": {:.3},",
        args.mobile_element_min_identity
    )?;
    writeln!(
        w,
        "    \"mobile_element_min_coverage\": {:.3},",
        args.mobile_element_min_coverage
    )?;
    writeln!(
        w,
        "    \"mobile_element_window\": {}",
        args.mobile_element_window
    )?;
    writeln!(w, "  }},")?;

    writeln!(w, "  \"counts\": {{")?;
    writeln!(w, "    \"clusters\": {},", pa.n_clusters())?;
    writeln!(w, "    \"genomes\": {},", pa.n_genomes())?;
    writeln!(w, "    \"total_gaps\": {},", stats.total_gaps)?;
    writeln!(w, "    \"gaps_passing_filter\": {},", stats.filtered_gaps)?;
    writeln!(w, "    \"gaps_searched\": {},", stats.searchable_gaps)?;
    writeln!(
        w,
        "    \"skipped_by_prevalence\": {},",
        stats.skipped_by_prevalence
    )?;
    writeln!(
        w,
        "    \"missing_fasta_gaps\": {},",
        stats.missing_fasta_gaps
    )?;
    writeln!(w, "    \"rescued\": {},", rescued)?;
    writeln!(
        w,
        "    \"matrix_corrected_total\": {},",
        total_matrix_corrected
    )?;
    writeln!(
        w,
        "    \"matrix_corrected_rescue_cells\": {},",
        stats.matrix_corrected
    )?;
    writeln!(
        w,
        "    \"matrix_corrected_nsc_cells\": {},",
        stats.nsc_matrix_corrected
    )?;
    writeln!(w, "    \"synteny_rejected\": {},", stats.synteny_rejected)?;
    writeln!(w, "    \"paralog_rejected\": {},", stats.paralog_rejected)?;
    writeln!(w, "    \"nsc_clusters\": {},", stats.nsc_clusters)?;
    writeln!(w, "    \"nsc_merged_genes\": {}", stats.nsc_merged_genes)?;
    writeln!(w, "  }}")?;
    writeln!(w, "}}")?;
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_args_for_validation() -> Args {
        Args {
            presence_absence: "pa.csv".into(),
            pan_reference: "pan_ref.fa".into(),
            genomes: "genomes".into(),
            annotations: None,
            output: "recore_output".into(),
            min_identity: 90.0,
            min_coverage: 80.0,
            threads: 4,
            preset: "asm5".into(),
            min_prevalence: 0.1,
            multi_allele: false,
            synteny: false,
            write_corrected_gff: false,
            matrix_policy: MatrixPolicy::Balanced,
            reference_genome: None,
            tree: None,
            mobile_element_reference: None,
            mobile_element_min_identity: 85.0,
            mobile_element_min_coverage: 50.0,
            mobile_element_window: 1000,
            nscluster: false,
            nsc_kmer_size: 21,
            nsc_threshold: 0.8,
            nsc_engine: "minhash".into(),
            verbose: false,
            nsc_min_seq_id: 0.9,
            nsc_coverage: 0.8,
        }
    }

    fn rescue_with_type(rescue_type: RescueType, synteny_ok: Option<bool>) -> RescueEvent {
        RescueEvent {
            cluster_idx: ClusterIdx(0),
            genome_idx: GenomeIdx(0),
            rescue_type,
            identity: 90.0,
            coverage: 90.0,
            contig: "chr".into(),
            start: 0,
            end: 100,
            strand: '+',
            annotation_match: None,
            confidence: 0.8,
            synteny_ok,
            synteny_evidence: SyntenyEvidence::default(),
            fragment_evidence: crate::rescue::FragmentEvidence::default(),
        }
    }

    #[test]
    fn invalid_args_are_rejected_before_running_pipeline() {
        let mut args = default_args_for_validation();
        args.threads = 0;
        assert!(validate_args(&args).is_err());

        let mut args = default_args_for_validation();
        args.min_identity = 101.0;
        assert!(validate_args(&args).is_err());

        let mut args = default_args_for_validation();
        args.min_prevalence = -0.1;
        assert!(validate_args(&args).is_err());

        let mut args = default_args_for_validation();
        args.nsc_kmer_size = 0;
        assert!(validate_args(&args).is_err());

        let mut args = default_args_for_validation();
        args.nsc_engine = "minhsah".into();
        assert!(validate_args(&args).is_err());
    }

    #[test]
    fn default_args_pass_validation() {
        let args = default_args_for_validation();
        validate_args(&args).unwrap();
    }

    #[test]
    fn rescue_orf_signal_flags_disrupted_and_frame_review_spans() {
        let mut rescue = rescue_with_type(RescueType::Dropped, Some(true));
        rescue.start = 0;
        rescue.end = 99;
        assert_eq!(rescue_orf_signal(&rescue), "in_frame_span");

        rescue.end = 100;
        assert_eq!(rescue_orf_signal(&rescue), "out_of_frame_span_review");

        rescue.rescue_type = RescueType::Pseudogene {
            reason: PseudogeneReason::Frameshift,
        };
        assert_eq!(rescue_orf_signal(&rescue), "disrupted_orf");
    }

    #[test]
    fn paralog_requires_explicit_synteny_support() {
        let no_synteny = rescue_with_type(
            RescueType::Paralog {
                other_cluster: "other".into(),
            },
            None,
        );
        assert!(is_paralog_without_synteny_support(&no_synteny));

        let failed_synteny = rescue_with_type(
            RescueType::Paralog {
                other_cluster: "other".into(),
            },
            Some(false),
        );
        assert!(is_paralog_without_synteny_support(&failed_synteny));

        let passed_synteny = rescue_with_type(
            RescueType::Paralog {
                other_cluster: "other".into(),
            },
            Some(true),
        );
        assert!(!is_paralog_without_synteny_support(&passed_synteny));

        let dropped = rescue_with_type(RescueType::Dropped, None);
        assert!(!is_paralog_without_synteny_support(&dropped));
    }

    #[test]
    fn tsv_fields_are_single_line_safe() {
        assert_eq!(clean_tsv_field("a\tb\nc\rd\0e"), "a b c d e");
    }

    #[test]
    fn sensitive_policy_does_not_correct_paralog_state() {
        assert!(!matrix_policy_corrects(
            MatrixPolicy::Sensitive,
            "HOMOLOG_PARALOG"
        ));
    }

    #[test]
    fn reference_genome_resolution_accepts_embedded_as_id() {
        let pa = crate::parser::PresenceAbsence::new_for_test(
            vec!["group_1".into()],
            vec!["L4_AS00146".into(), "AS00999".into()],
            vec![vec!["gene1".into(), "".into()]],
        );

        assert_eq!(resolve_reference_genome(&pa, "AS00146"), Some(GenomeIdx(0)));
        assert_eq!(
            resolve_reference_genome(&pa, "L4_AS00146"),
            Some(GenomeIdx(0))
        );
        assert_eq!(resolve_reference_genome(&pa, "AS00000"), None);
    }

    #[test]
    fn evidence_states_separate_strict_alias_and_paralog() {
        let strict = rescue_with_type(RescueType::Dropped, Some(true));
        assert_eq!(
            evidence_state(&strict, true, "group_1", "", ""),
            "STRICT_PRESENT"
        );

        let alias = rescue_with_type(
            RescueType::Split {
                other_cluster: "group_2".into(),
            },
            Some(true),
        );
        assert_eq!(
            evidence_state(&alias, true, "group_1", "group_2", ""),
            "ALIAS_CLUSTER"
        );

        let paralog = rescue_with_type(
            RescueType::Paralog {
                other_cluster: "group_2".into(),
            },
            None,
        );
        assert_eq!(
            evidence_state(
                &paralog,
                false,
                "group_1",
                "group_2",
                "paralog_without_synteny"
            ),
            "HOMOLOG_PARALOG"
        );
    }

    #[test]
    fn reference_shadow_candidate_flags_partial_paralog_evidence() {
        let partial_hit = RescueEvent {
            coverage: 91.0,
            identity: 100.0,
            annotation_match: Some("rv_duf222".into()),
            ..rescue_with_type(
                RescueType::Split {
                    other_cluster: "group_537".into(),
                },
                None,
            )
        };

        let state = reference_shadow_state(&partial_hit, "group_1044", "group_537");
        assert_eq!(state, "PARTIAL_PARALOG_EVIDENCE");
        assert_eq!(
            reference_publication_status(state, false),
            "flag_reference_partial_paralog_evidence"
        );
        assert_eq!(
            reference_publication_confidence(state, reference_status_for(state, false)),
            "complex_paralog"
        );
    }

    #[test]
    fn reference_hit_shadow_candidates_capture_partial_annotated_paralog() {
        let pa = crate::parser::PresenceAbsence::new_for_test(
            vec!["group_1044".into(), "group_537".into()],
            vec!["L4_AS00146".into()],
            vec![vec!["".into()], vec!["rv_duf222".into()]],
        );
        let gff = GffIndex::build(vec![GffGene {
            contig: "chr".into(),
            start: 0,
            end: 954,
            strand: '+',
            gene_id: "rv_duf222".into(),
        }]);
        let mut gff_indices = HashMap::new();
        gff_indices.insert(GenomeIdx(0), gff);
        let mut gene_map = HashMap::new();
        gene_map.insert("rv_duf222".to_string(), ClusterIdx(1));
        let mut gene_to_cluster_maps = HashMap::new();
        gene_to_cluster_maps.insert(GenomeIdx(0), gene_map);
        let mut hits = HashMap::new();
        hits.insert(
            (GenomeIdx(0), ClusterIdx(0)),
            vec![Hit {
                query: "group_1044".into(),
                cluster: "group_1044".into(),
                query_length: 972,
                query_start: 0,
                query_end: 885,
                strand: '+',
                target: "chr".into(),
                target_start: 0,
                target_end: 885,
                matches: 885,
                alignment_length: 885,
                mapq: 0,
            }],
        );

        let candidates =
            reference_hit_shadow_candidates(&hits, &pa, &gff_indices, &gene_to_cluster_maps);

        assert_eq!(candidates.len(), 1);
        assert!(matches!(
            candidates[0].rescue_type,
            RescueType::Paralog { .. }
        ));
        assert_eq!(
            reference_shadow_state(&candidates[0], "group_1044", "group_537"),
            "PARTIAL_PARALOG_EVIDENCE"
        );
    }

    #[test]
    fn reference_cluster_events_detect_short_in_frame_allele() {
        let left = b"ATGAAACCCGGGTTTAAACCCGGGTTTAAACCC";
        let middle = b"GCTGCTGCTGCTGCTGCTGCTGCTGCTGCT";
        let right = b"GGGAAATTTCCCGGGAAATTTCCCGGGAAATAA";
        let mut reference_seq = Vec::new();
        reference_seq.extend_from_slice(left);
        reference_seq.extend_from_slice(middle);
        reference_seq.extend_from_slice(right);
        let mut short_seq = Vec::new();
        short_seq.extend_from_slice(left);
        short_seq.extend_from_slice(right);

        assert_eq!((reference_seq.len() - short_seq.len()) % 3, 0);

        let pa = crate::parser::PresenceAbsence::new_for_test(
            vec!["group_short".into(), "group_ref".into()],
            vec!["L4_AS00146".into()],
            vec![vec!["".into()], vec!["Rv0601c".into()]],
        );
        let mut pan_ref = HashMap::new();
        pan_ref.insert("group_short".to_string(), vec![short_seq]);
        pan_ref.insert("group_ref".to_string(), vec![reference_seq]);

        let events = build_reference_cluster_events(&pa, &pan_ref, GenomeIdx(0));

        assert_eq!(events.len(), 1);
        let event = &events[0];
        assert_eq!(event.query_cluster_idx, ClusterIdx(0));
        assert_eq!(event.target_cluster_idx, ClusterIdx(1));
        assert_eq!(event.target_gene_id, "Rv0601c");
        assert_eq!(event.event_class, "short_in_frame_allele");
        assert_eq!(event.event_size_nt, -30);
        assert_eq!(event.event_size_aa, Some(-10));
        assert!(event.identity >= 99.0);
    }

    #[test]
    fn reference_gene_display_prefers_h37rv_locus() {
        assert_eq!(display_reference_gene_id("L4_AS00146_Rv0601c"), "Rv0601c");
        assert_eq!(
            display_reference_gene_id("L4_AS00146_Rv0601c_pseudo"),
            "Rv0601c"
        );
        assert_eq!(display_reference_gene_id("custom_gene"), "custom_gene");
    }

    #[test]
    fn newick_fitch_detects_clade_gain() {
        let pa = crate::parser::PresenceAbsence::new_for_test(
            vec!["g".into()],
            vec!["s1".into(), "s2".into(), "s3".into(), "s4".into()],
            vec![vec!["".into(), "g_s2".into(), "g_s3".into(), "".into()]],
        );
        let mut root = NewickParser::new("(s1,(s2:0.1,s3:0.1),s4);")
            .parse()
            .unwrap();
        let (matched, unmatched) = resolve_phylo_tips(&mut root, &pa);
        assert_eq!(matched, 4);
        assert_eq!(unmatched, 0);
        let tree = PhyloTree { root };
        let present = vec![false, true, true, false];

        let summary = summarize_phylogenetic_pattern(Some(&tree), &present);

        assert_eq!(summary.present_tips, 2);
        assert_eq!(summary.gains, 1);
        assert_eq!(summary.losses, 0);
        assert!(summary.carrier_monophyly);
        assert_eq!(summary.pattern, "lineage_specific_gain");
    }

    #[test]
    fn mobile_element_detection_uses_gff_aliases() {
        let pa = crate::parser::PresenceAbsence::new_for_test(
            vec!["g".into()],
            vec!["s1".into()],
            vec![vec!["".into()]],
        );
        let gff = GffIndex::build_with_aliases(vec![(
            GffGene {
                contig: "chr".into(),
                start: 150,
                end: 300,
                strand: '+',
                gene_id: "mobile1".into(),
            },
            vec!["IS6110 transposase".into()],
        )]);
        let mut gff_indices = HashMap::new();
        gff_indices.insert(GenomeIdx(0), gff);
        let rescue = RescueEvent {
            cluster_idx: ClusterIdx(0),
            genome_idx: GenomeIdx(0),
            rescue_type: RescueType::Fragmented {
                n_fragments: 2,
                total_coverage: 91.0,
            },
            identity: 99.0,
            coverage: 91.0,
            contig: "chr".into(),
            start: 100,
            end: 220,
            strand: '+',
            annotation_match: None,
            confidence: 0.9,
            synteny_ok: None,
            synteny_evidence: SyntenyEvidence::default(),
            fragment_evidence: crate::rescue::FragmentEvidence::default(),
        };

        let mobile =
            nearest_mobile_element_for_rescue(&rescue, &pa, &gff_indices, &[], 1_000).unwrap();

        assert_eq!(mobile.source, "annotation");
        assert_eq!(mobile.kind, "IS6110");
        assert_eq!(mobile.distance_bp, 0);
    }

    #[test]
    fn mobile_element_sequence_hits_are_name_agnostic_and_prioritized() {
        let pa = crate::parser::PresenceAbsence::new_for_test(
            vec!["g".into()],
            vec!["s1".into()],
            vec![vec!["".into()]],
        );
        let mut gff_indices = HashMap::new();
        gff_indices.insert(
            GenomeIdx(0),
            GffIndex::build_with_aliases(vec![(
                GffGene {
                    contig: "chr".into(),
                    start: 90,
                    end: 210,
                    strand: '+',
                    gene_id: "anonymous_cds".into(),
                },
                vec!["hypothetical protein".into()],
            )]),
        );
        let rescue = RescueEvent {
            cluster_idx: ClusterIdx(0),
            genome_idx: GenomeIdx(0),
            rescue_type: RescueType::Fragmented {
                n_fragments: 2,
                total_coverage: 91.0,
            },
            identity: 99.0,
            coverage: 91.0,
            contig: "chr".into(),
            start: 100,
            end: 220,
            strand: '+',
            annotation_match: None,
            confidence: 0.9,
            synteny_ok: None,
            synteny_evidence: SyntenyEvidence::default(),
            fragment_evidence: crate::rescue::FragmentEvidence::default(),
        };
        let mobile_hits = vec![MobileElementHit {
            genome_idx: GenomeIdx(0),
            query_name: "ISfinder_anon_001".into(),
            contig: "chr".into(),
            start: 120,
            end: 260,
            strand: '+',
            identity: 97.5,
            coverage: 82.0,
        }];

        let mobile =
            nearest_mobile_element_for_rescue(&rescue, &pa, &gff_indices, &mobile_hits, 1_000)
                .unwrap();

        assert_eq!(mobile.source, "sequence");
        assert_eq!(mobile.kind, "mobile_element_sequence");
        assert_eq!(mobile.distance_bp, 0);
        assert!(mobile.label.contains("ISfinder_anon_001"));
    }

    #[test]
    fn nsc_role_labels_member_and_representative() {
        let mut merge_map = HashMap::new();
        merge_map.insert(ClusterIdx(2), ClusterIdx(1));
        let nsc = Some(nscluster::NsClusterResult {
            merged_clusters: vec![nscluster::NsCluster {
                members: vec![ClusterIdx(1), ClusterIdx(2)],
                representative: "group_1".into(),
            }],
            total_merged_genes: 2,
            n_clusters: 1,
            merge_map,
        });

        assert_eq!(
            nsc_role_for_cluster(ClusterIdx(2), &nsc),
            "nsc_member_merged_into_representative"
        );
        assert_eq!(
            nsc_role_for_cluster(ClusterIdx(1), &nsc),
            "nsc_representative"
        );
        assert_eq!(nsc_role_for_cluster(ClusterIdx(3), &nsc), "");
    }

    #[test]
    fn nsc_presence_union_includes_matrix_rescued_members() {
        let pa = crate::parser::PresenceAbsence::new_for_test(
            vec!["rep".into(), "member".into()],
            vec!["g1".into(), "g2".into()],
            vec![vec!["rep_g1".into(), "".into()], vec!["".into(), "".into()]],
        );
        let mut merge_map = HashMap::new();
        merge_map.insert(ClusterIdx(1), ClusterIdx(0));
        let nsc = Some(nscluster::NsClusterResult {
            merged_clusters: vec![nscluster::NsCluster {
                members: vec![ClusterIdx(0), ClusterIdx(1)],
                representative: "rep".into(),
            }],
            total_merged_genes: 2,
            n_clusters: 1,
            merge_map,
        });
        let member_rescue = RescueEvent {
            cluster_idx: ClusterIdx(1),
            genome_idx: GenomeIdx(1),
            ..rescue_with_type(RescueType::Dropped, Some(true))
        };
        let mut rescue_map = HashMap::new();
        rescue_map.insert(
            (member_rescue.cluster_idx, member_rescue.genome_idx),
            &member_rescue,
        );

        let updates = nsc_presence_updates(&pa, &nsc, &rescue_map);

        assert!(updates.contains(&(0, 1)));
    }

    #[test]
    fn reference_shadow_candidate_prefers_synteny_supported_evidence() {
        let weak = RescueEvent {
            cluster_idx: ClusterIdx(1),
            genome_idx: GenomeIdx(2),
            confidence: 1.5,
            synteny_ok: Some(false),
            ..rescue_with_type(RescueType::Dropped, Some(false))
        };
        let supported = RescueEvent {
            cluster_idx: ClusterIdx(1),
            genome_idx: GenomeIdx(2),
            confidence: 0.9,
            synteny_ok: Some(true),
            ..rescue_with_type(RescueType::Dropped, Some(true))
        };
        let candidates = vec![weak, supported];

        let best = best_reference_candidate(&candidates, ClusterIdx(1), GenomeIdx(2)).unwrap();
        assert_eq!(best.synteny_ok, Some(true));
    }

    #[test]
    fn matrix_policy_controls_binary_correction() {
        assert!(matrix_policy_corrects(
            MatrixPolicy::Strict,
            "STRICT_PRESENT"
        ));
        assert!(matrix_policy_corrects(
            MatrixPolicy::Strict,
            "ALIAS_CLUSTER"
        ));
        assert!(!matrix_policy_corrects(
            MatrixPolicy::Strict,
            "SEQUENCE_ONLY"
        ));

        assert!(matrix_policy_corrects(
            MatrixPolicy::Balanced,
            "SEQUENCE_ONLY"
        ));
        assert!(!matrix_policy_corrects(
            MatrixPolicy::Balanced,
            "PSEUDOGENE"
        ));
        assert!(!matrix_policy_corrects(
            MatrixPolicy::Balanced,
            "PARTIAL_HOMOLOGY"
        ));
        assert!(!matrix_policy_corrects(
            MatrixPolicy::Balanced,
            "MISSING_CDS_UNVALIDATED"
        ));

        assert!(matrix_policy_corrects(
            MatrixPolicy::Sensitive,
            "PSEUDOGENE"
        ));
        assert!(matrix_policy_corrects(
            MatrixPolicy::Sensitive,
            "PARTIAL_HOMOLOGY"
        ));
        assert!(!matrix_policy_corrects(
            MatrixPolicy::Sensitive,
            "MISSING_CDS_UNVALIDATED"
        ));
    }

    #[test]
    fn gff_stream_copy_adds_missing_header() {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            input,
            "chr1\tAMAP\tCDS\t1\t90\t.\t+\t0\tID=gene1;Name=gene1"
        )
        .unwrap();

        let mut output = Vec::new();
        copy_text_maybe_gzip_with_header(input.path(), &mut output).unwrap();
        let text = String::from_utf8(output).unwrap();

        assert!(text.starts_with("##gff-version 3\n"));
        assert!(text.contains("ID=gene1"));
    }

    #[test]
    fn gff_stream_copy_preserves_existing_header_after_blank_prefix() {
        let mut input = tempfile::NamedTempFile::new().unwrap();
        writeln!(input).unwrap();
        writeln!(input, "##gff-version 3").unwrap();
        writeln!(
            input,
            "chr1\tAMAP\tCDS\t1\t90\t.\t+\t0\tID=gene1;Name=gene1"
        )
        .unwrap();

        let mut output = Vec::new();
        copy_text_maybe_gzip_with_header(input.path(), &mut output).unwrap();
        let text = String::from_utf8(output).unwrap();

        assert_eq!(text.matches("##gff-version 3").count(), 1);
        assert!(text.starts_with('\n'));
        assert!(text.contains("ID=gene1"));
    }
}
