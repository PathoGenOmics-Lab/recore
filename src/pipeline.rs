use ahash::{HashMap, HashMapExt};
use anyhow::{Context, Result};
use log::{debug, info, warn};
use std::collections::BTreeMap;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use crate::cli::Args;
use crate::nscluster;
use crate::parser::{
    find_annotation_file, find_genome_files, load_genome_contig_names,
    normalize_gff_records_to_fasta, parse_gff_genes, parse_pan_reference, parse_presence_absence,
    resolve_genome_files_for_pa, validate_inputs,
};
use crate::relationships::write_cluster_relationships;
use crate::rescue::{
    build_gene_to_cluster, check_internal_stops_in_hit, classify_hits, resolve_conflicts,
    validate_synteny, GffGene, GffIndex, PseudogeneReason, RescueEvent, RescueType,
};
use crate::search::batch_search;
use crate::types::{ClusterIdx, GenomeIdx};

pub fn run(args: &Args) -> Result<()> {
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

    let cluster_to_idx = pa.cluster_index();

    // ── 1b. Cross-validate inputs (includes path traversal checks) ──
    validate_inputs(&pa, &pan_ref, &genome_files)?;

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
    let mut missing_names: Vec<String> = Vec::new();
    for (gi, cluster_indices) in gaps_by_genome {
        let genome_name = pa.genome_name(gi);
        if genome_files.contains_key(genome_name) {
            searchable_genomes.push((gi, cluster_indices));
        } else {
            missing_genomes += 1;
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
                    let genes: Vec<GffGene> = raw_genes
                        .into_iter()
                        .map(|(contig, start, end, strand, gene_id)| GffGene {
                            contig,
                            start,
                            end,
                            strand,
                            gene_id,
                        })
                        .collect();
                    data.insert(gi, GffIndex::build(genes));
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
            maps.entry(*gi)
                .or_insert_with(|| build_gene_to_cluster(&pa, *gi));
        }
        if args.synteny || !gff_indices.is_empty() {
            for &gi in gff_indices.keys() {
                maps.entry(gi)
                    .or_insert_with(|| build_gene_to_cluster(&pa, gi));
            }
        }
        maps
    };

    // ── 6. Batch search: ONE minimap2 invocation ──
    // Adapt types for batch_search (uses usize internally for genome files lookup)
    let searchable_usize: Vec<(usize, Vec<usize>)> = searchable_genomes
        .iter()
        .map(|(gi, cis)| (gi.val(), cis.iter().map(|ci| ci.val()).collect()))
        .collect();

    let hits_map = batch_search(
        &pan_ref,
        &genome_files,
        &searchable_usize,
        pa.genome_names(),
        pa.cluster_names(),
        args.threads,
        args.min_identity,
        &args.preset,
    )?;

    // ── 7. Classify hits (parallel) ──
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .context("Failed to build rayon thread pool")?;

    let classify_inputs: Vec<(ClusterIdx, GenomeIdx, &[crate::search::Hit])> = hits_map
        .iter()
        .filter_map(|((gi_val, cluster_name), hits)| {
            let ci = *cluster_to_idx.get(cluster_name.as_str())?;
            let gi = GenomeIdx(*gi_val);
            Some((ci, gi, hits.as_slice()))
        })
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

    // ── 7b. Internal stop codon scanning for DROPPED rescues ──
    // Load genome sequences on-demand for rescued genomes only
    {
        let rescued_genomes: ahash::HashSet<GenomeIdx> = rescues
            .iter()
            .filter(|r| matches!(r.rescue_type, RescueType::Dropped))
            .map(|r| r.genome_idx)
            .collect();

        if !rescued_genomes.is_empty() {
            info!(
                "Scanning {} genomes for internal stop codons in DROPPED rescues...",
                rescued_genomes.len()
            );
            let mut upgraded = 0usize;

            for gi in &rescued_genomes {
                let genome_name = pa.genome_name(*gi);
                let Some(fasta_path) = genome_files.get(genome_name) else {
                    continue;
                };
                let Ok(contigs) = crate::parser::load_genome_contigs(fasta_path) else {
                    continue;
                };

                for rescue in rescues.iter_mut() {
                    if rescue.genome_idx != *gi
                        || !matches!(rescue.rescue_type, RescueType::Dropped)
                    {
                        continue;
                    }
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
                    "Upgraded {} DROPPED rescues to PSEUDOGENE:internal_stop",
                    upgraded
                );
            }
        }
    }

    // ── 8. Synteny validation ──
    if args.synteny && !gff_indices.is_empty() {
        info!("Validating synteny...");
        for rescue in &mut rescues {
            validate_synteny(
                rescue,
                &pa,
                &gff_indices,
                &gene_to_cluster_maps,
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
    let mut split = 0usize;
    let mut fragmented = 0usize;
    let mut alignment_split = 0usize;
    let mut start_variant = 0usize;
    let mut paralog = 0usize;
    let mut pseudogene = 0usize;
    for r in &rescues {
        match &r.rescue_type {
            RescueType::Dropped => dropped += 1,
            RescueType::Split { .. } => split += 1,
            RescueType::Fragmented { .. } => fragmented += 1,
            RescueType::AlignmentSplit { .. } => alignment_split += 1,
            RescueType::StartVariant { .. } => start_variant += 1,
            RescueType::Paralog { .. } => paralog += 1,
            RescueType::Pseudogene { .. } => pseudogene += 1,
        }
    }
    info!(
        "Final rescues: {} (DROPPED: {}, SPLIT: {}, FRAGMENTED: {}, ALIGNMENT_SPLIT: {}, START_VARIANT: {}, PARALOG: {}, PSEUDOGENE: {})",
        rescues.len(), dropped, split, fragmented, alignment_split, start_variant, paralog, pseudogene
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

    let rescue_map: HashMap<(ClusterIdx, GenomeIdx), &RescueEvent> = rescues
        .iter()
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
    write_orthology_evidence(
        &rescues,
        &synteny_rejected,
        &pa,
        &gene_to_cluster_maps,
        &args.output,
    )?;
    if !gff_indices.is_empty() {
        let (context_rows, summary_rows) =
            write_cluster_locus_contexts(&pa, &gff_indices, &gene_to_cluster_maps, &args.output)?;
        info!(
            "Cluster locus context audit: {} context rows across {} clusters",
            context_rows, summary_rows
        );
    }
    write_corrected_csv(&pa, &rescue_map, &args.output, &nsc_result)?;
    let n_relationships = write_cluster_relationships(&pa, &rescues, &nsc_result, &args.output)?;
    if n_relationships > 0 {
        info!("Cluster relationships: {} rows", n_relationships);
    }
    if args.write_corrected_gff {
        write_corrected_gffs(&rescues, &pa, &gff_paths, &args.output)?;
    }

    let stats = SummaryStats {
        core_before,
        total_gaps,
        filtered_gaps,
        skipped_by_prevalence,
        dropped,
        split,
        fragmented,
        alignment_split,
        start_variant,
        paralog,
        pseudogene,
        synteny_rejected: synteny_rejected.len(),
        synteny_enabled: args.synteny,
        nsc_clusters: nsc_result.as_ref().map(|r| r.n_clusters).unwrap_or(0),
        nsc_merged_genes: nsc_result
            .as_ref()
            .map(|r| r.total_merged_genes)
            .unwrap_or(0),
    };
    write_summary(
        &pa,
        &rescue_map,
        &rescues,
        &stats,
        &nsc_result,
        &args.output,
    )?;

    // Write NSC report if clusters were merged
    if let Some(ref nsc) = nsc_result {
        if !nsc.merged_clusters.is_empty() {
            write_nsc_report(nsc, pa.cluster_names(), &args.output)?;
        }
    }

    let core_after = count_core_after(&pa, &rescue_map, &nsc_result);
    let soft_core_after = count_soft_core_after(&pa, &rescue_map, &nsc_result);
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
        "gene_cluster\tgenome\trescue_type\tidentity\tcoverage\tcontig\tstart\tend\tstrand\tconfidence\tsynteny\tannotation_match\tfragment_n\tfragment_query_gap_bp\tfragment_target_gap_bp\tfragment_colinear\tfragment_same_cds\tfragment_note"
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
        writeln!(
            w,
            "{}\t{}\t{}\t{:.1}\t{:.1}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            pa.cluster_name(r.cluster_idx),
            pa.genome_name(r.genome_idx),
            r.rescue_type,
            r.identity,
            r.coverage,
            r.contig,
            r.start + 1, // convert 0-based to 1-based for TSV output
            r.end,
            r.strand,
            r.confidence,
            match r.synteny_ok {
                Some(true) => "PASS",
                Some(false) => "FAIL",
                None => "",
            },
            r.annotation_match.as_deref().unwrap_or(""),
            fragment_n,
            query_gap,
            target_gap,
            colinear,
            same_cds,
            r.fragment_evidence.note.as_deref().unwrap_or(""),
        )?;
    }
    w.flush()?;
    Ok(())
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
        RescueType::Dropped => "accepted_sequence_evidence",
    }
}

fn write_orthology_evidence(
    rescues: &[RescueEvent],
    synteny_rejected: &[RescueEvent],
    pa: &crate::parser::PresenceAbsence,
    gene_to_cluster_maps: &HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
    output_dir: &std::path::Path,
) -> Result<()> {
    let path = output_dir.join("orthology_evidence.tsv");
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "decision\torthology_status\tgene_cluster\tgenome\trescue_type\tidentity\tcoverage\tcontig\tstart\tend\tstrand\tconfidence\tsynteny\tannotation_match\tannotation_cluster\trejection_reason\tfragment_note"
    )?;

    for (accepted, rescue) in rescues
        .iter()
        .map(|r| (true, r))
        .chain(synteny_rejected.iter().map(|r| (false, r)))
    {
        let query_cluster = pa.cluster_name(rescue.cluster_idx);
        let annotation_cluster = annotation_cluster_name(rescue, pa, gene_to_cluster_maps);
        let status = if accepted {
            accepted_orthology_status(rescue, query_cluster, &annotation_cluster)
        } else {
            "rejected_incompatible_context"
        };
        let rejection_reason = if accepted { "" } else { "synteny_fail" };
        let decision = if accepted { "accepted" } else { "rejected" };
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{:.1}\t{:.1}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}",
            decision,
            status,
            query_cluster,
            pa.genome_name(rescue.genome_idx),
            rescue.rescue_type,
            rescue.identity,
            rescue.coverage,
            rescue.contig,
            rescue.start + 1,
            rescue.end,
            rescue.strand,
            rescue.confidence,
            match rescue.synteny_ok {
                Some(true) => "PASS",
                Some(false) => "FAIL",
                None => "",
            },
            rescue.annotation_match.as_deref().unwrap_or(""),
            annotation_cluster,
            rejection_reason,
            rescue.fragment_evidence.note.as_deref().unwrap_or(""),
        )?;
    }
    w.flush()?;
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
    genes: &[&GffGene],
    gene_to_cluster: &HashMap<String, ClusterIdx>,
    pa: &crate::parser::PresenceAbsence,
) -> String {
    let names: Vec<&str> = genes
        .iter()
        .filter_map(|gene| gene_to_cluster.get(&gene.gene_id).copied())
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
) -> Result<(usize, usize)> {
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
            let upstream_sig = context_signature(&upstream, g2c, pa);
            let downstream_sig = context_signature(&downstream, g2c, pa);
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
            writeln!(
                context_w,
                "{}\t{}\t{}\t{:.3}\t{}\t{}\t{}",
                pa.cluster_name(ClusterIdx(ci_val)),
                rank + 1,
                acc.count,
                acc.count as f64 / mapped as f64,
                upstream,
                downstream,
                acc.examples.join(";"),
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
        writeln!(
            summary_w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}",
            pa.cluster_name(summary.cluster_idx),
            summary.n_present,
            summary.n_mapped,
            summary.n_unresolved,
            summary.n_unique_contexts,
            summary.dominant_count,
            summary.dominant_fraction,
            summary.interpretation,
        )?;
    }
    summary_w.flush()?;

    Ok((written_context_rows, summaries.len()))
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
    let nsc_presence: ahash::HashSet<(usize, usize)> = nsc_result
        .as_ref()
        .map(|r| {
            let (_, updates) = nscluster::apply_ns_merges(
                &r.merge_map,
                |ci, gi| pa.is_present(ci, gi),
                pa.n_genomes(),
            );
            updates
                .into_iter()
                .map(|(ci, gi)| (ci.val(), gi.val()))
                .collect()
        })
        .unwrap_or_default();
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
                row[col_idx].clone_from(&pa.cell(ci, gi).to_owned());
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

fn read_text_maybe_gzip(path: &std::path::Path) -> Result<String> {
    use std::io::Read;

    let file = fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    let mut text = String::new();
    if path.extension().is_some_and(|e| e == "gz") {
        let mut decoder = flate2::read::GzDecoder::new(file);
        decoder.read_to_string(&mut text)?;
    } else {
        let mut reader = std::io::BufReader::new(file);
        reader.read_to_string(&mut text)?;
    }
    Ok(text)
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
            let original = read_text_maybe_gzip(source_gff)?;
            if original.trim_start().starts_with("##gff-version") {
                write!(w, "{}", original)?;
            } else {
                writeln!(w, "##gff-version 3")?;
                write!(w, "{}", original)?;
            }
            if !original.ends_with('\n') {
                writeln!(w)?;
            }
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
    let nsc_presence: ahash::HashSet<(usize, usize)> = nsc_result
        .as_ref()
        .map(|r| {
            let (_, updates) =
                nscluster::apply_ns_merges(&r.merge_map, |ci, gi| pa.is_present(ci, gi), n);
            updates
                .into_iter()
                .map(|(ci, gi)| (ci.val(), gi.val()))
                .collect()
        })
        .unwrap_or_default();
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
    let nsc_presence: ahash::HashSet<(usize, usize)> = nsc_result
        .as_ref()
        .map(|r| {
            let (_, updates) =
                nscluster::apply_ns_merges(&r.merge_map, |ci, gi| pa.is_present(ci, gi), n);
            updates
                .into_iter()
                .map(|(ci, gi)| (ci.val(), gi.val()))
                .collect()
        })
        .unwrap_or_default();
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
    dropped: usize,
    split: usize,
    fragmented: usize,
    alignment_split: usize,
    start_variant: usize,
    paralog: usize,
    pseudogene: usize,
    synteny_rejected: usize,
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
        writeln!(
            w,
            "NS_{}\t{}\t{}\t{}",
            i,
            nsc_cl.representative,
            nsc_cl.members.len(),
            member_names.join(";")
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
    writeln!(w, "Gaps searched:         {}", stats.filtered_gaps)?;
    writeln!(w, "Skipped (prevalence):  {}", stats.skipped_by_prevalence)?;
    writeln!(w, "Rescued:               {}", rescues.len())?;
    if stats.synteny_rejected > 0 {
        writeln!(w, "Rejected (synteny):    {}", stats.synteny_rejected)?;
    }
    writeln!(
        w,
        "Remaining gaps:        {}",
        stats.filtered_gaps - rescues.len()
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
