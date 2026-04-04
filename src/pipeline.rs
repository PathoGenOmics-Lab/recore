use ahash::{HashMap, HashMapExt};
use anyhow::Result;
use log::{info, warn};
use rayon::prelude::*;
use std::fs;
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use crate::cli::Args;
use crate::parser::{
    find_genome_files, parse_gff_genes, parse_pan_reference, parse_presence_absence,
};
use crate::rescue::{
    build_gene_to_cluster, classify_hits, resolve_conflicts, validate_synteny,
    GffGene, GffIndex, RescueEvent, RescueType,
};
use crate::search::search_genome;

pub fn run(args: &Args) -> Result<()> {
    // ── 1. Parse inputs ──
    let pa = parse_presence_absence(&args.presence_absence)?;
    let pan_ref = parse_pan_reference(&args.pan_reference, args.multi_allele)?;
    let genome_files = find_genome_files(&args.genomes)?;

    // Cluster name → index for O(1) lookup
    let cluster_to_idx: HashMap<String, usize> = pa
        .clusters
        .iter()
        .enumerate()
        .map(|(i, c)| (c.clone(), i))
        .collect();

    // ── 2. Filter gaps by prevalence ──
    let all_gaps = pa.gaps();

    // Apply prevalence filter: skip clusters present in too few genomes
    let min_prevalence = args.min_prevalence;
    let mut gaps_by_genome: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut filtered_gaps = 0usize;
    let mut skipped_by_prevalence = 0usize;

    for &(ci, gi) in &all_gaps {
        if pa.prevalence(ci) < min_prevalence {
            skipped_by_prevalence += 1;
            continue;
        }
        gaps_by_genome.entry(gi).or_default().push(ci);
        filtered_gaps += 1;
    }

    let core_before = pa.core_size();
    info!(
        "Gaps: {} total, {} after prevalence filter (≥{:.0}%), {} skipped",
        all_gaps.len(),
        filtered_gaps,
        min_prevalence * 100.0,
        skipped_by_prevalence
    );
    info!(
        "Core genome (before): {} / {} clusters",
        core_before,
        pa.clusters.len()
    );

    // Check genome file availability
    let mut searchable_genomes: Vec<(usize, Vec<usize>)> = Vec::new();
    let mut missing_genomes = 0usize;
    for (gi, cluster_indices) in &gaps_by_genome {
        let genome_name = &pa.genomes[*gi];
        if genome_files.contains_key(genome_name) {
            searchable_genomes.push((*gi, cluster_indices.clone()));
        } else {
            missing_genomes += 1;
        }
    }
    if missing_genomes > 0 {
        warn!(
            "{} genomes with gaps have no FASTA file (skipped)",
            missing_genomes
        );
    }

    // ── 3. Load GFF annotations ──
    let gff_indices: HashMap<usize, GffIndex> = if let Some(ref ann_dir) = args.annotations {
        let mut data = HashMap::new();
        for &(gi, _) in &searchable_genomes {
            let genome_name = &pa.genomes[gi];
            let gff_path = ann_dir.join(format!("{}.gff", genome_name));
            let alt_path = ann_dir.join(format!("{}.gff3", genome_name));
            let path = if gff_path.exists() {
                gff_path
            } else if alt_path.exists() {
                alt_path
            } else {
                continue;
            };
            match parse_gff_genes(&path) {
                Ok(raw_genes) => {
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
                }
                Err(e) => {
                    warn!("Cannot parse GFF for {}: {}", genome_name, e);
                }
            }
        }
        if !data.is_empty() {
            info!("Loaded GFF annotations for {} genomes", data.len());
        }
        data
    } else {
        HashMap::new()
    };

    // ── 4. Search and rescue (parallel per genome) ──
    let all_rescues: Mutex<Vec<RescueEvent>> = Mutex::new(Vec::new());
    let progress = AtomicUsize::new(0);
    let total = searchable_genomes.len();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .ok();

    searchable_genomes
        .par_iter()
        .for_each(|(gi, cluster_indices)| {
            let genome_name = &pa.genomes[*gi];
            let genome_fasta = &genome_files[genome_name];

            // Build query set: cluster → allele sequences (borrowed)
            let mut queries: HashMap<String, Vec<&[u8]>> = HashMap::new();
            for &ci in cluster_indices {
                let cluster_name = &pa.clusters[ci];
                if let Some(alleles) = pan_ref.get(cluster_name) {
                    let borrowed: Vec<&[u8]> = alleles.iter().map(|a| a.as_slice()).collect();
                    queries.insert(cluster_name.clone(), borrowed);
                }
            }

            if queries.is_empty() {
                return;
            }

            // Search
            let hits_per_cluster = match search_genome(
                &queries,
                genome_fasta,
                1,
                args.min_identity,
                &args.preset,
            ) {
                Ok(h) => h,
                Err(e) => {
                    warn!("minimap2 failed for {}: {}", genome_name, e);
                    return;
                }
            };

            // Build gene_id → cluster index for this genome
            let gene_to_cluster = build_gene_to_cluster(&pa, *gi);
            let gff_index = gff_indices.get(gi);

            // Classify
            let mut rescues = Vec::new();
            for (cluster_name, hits) in &hits_per_cluster {
                let Some(&ci) = cluster_to_idx.get(cluster_name.as_str()) else {
                    continue;
                };
                if let Some(rescue) = classify_hits(
                    hits,
                    ci,
                    *gi,
                    &pa,
                    gff_index,
                    &gene_to_cluster,
                    args.min_coverage,
                ) {
                    rescues.push(rescue);
                }
            }

            if !rescues.is_empty() {
                all_rescues.lock().unwrap().extend(rescues);
            }

            let count = progress.fetch_add(1, Ordering::Relaxed) + 1;
            if count.is_multiple_of(50) || count == total {
                info!("Processed {}/{} genomes", count, total);
            }
        });

    let mut rescues = all_rescues.into_inner().unwrap();
    info!("Raw rescues before filtering: {}", rescues.len());

    // ── 5. Synteny validation (if requested) ──
    if args.synteny && !gff_indices.is_empty() {
        info!("Validating synteny...");
        // Build gene_to_cluster maps for all genomes with GFF
        let mut gene_to_cluster_maps: HashMap<usize, HashMap<String, usize>> = HashMap::new();
        for &gi in gff_indices.keys() {
            gene_to_cluster_maps.insert(gi, build_gene_to_cluster(&pa, gi));
        }
        for rescue in &mut rescues {
            validate_synteny(rescue, &pa, &gff_indices, &gene_to_cluster_maps);
        }
        let validated = rescues.iter().filter(|r| r.synteny_ok == Some(true)).count();
        let failed = rescues.iter().filter(|r| r.synteny_ok == Some(false)).count();
        info!("Synteny: {} validated, {} failed", validated, failed);
    }

    // ── 6. Conflict resolution ──
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

    // ── 7. Count rescue types ──
    let mut dropped = 0usize;
    let mut split = 0usize;
    let mut fragmented = 0usize;
    let mut start_variant = 0usize;
    let mut paralog = 0usize;
    for r in &rescues {
        match &r.rescue_type {
            RescueType::Dropped => dropped += 1,
            RescueType::Split { .. } => split += 1,
            RescueType::Fragmented { .. } => fragmented += 1,
            RescueType::StartVariant { .. } => start_variant += 1,
            RescueType::Paralog { .. } => paralog += 1,
        }
    }
    info!(
        "Final rescues: {} (DROPPED: {}, SPLIT: {}, FRAGMENTED: {}, START_VARIANT: {}, PARALOG: {})",
        rescues.len(), dropped, split, fragmented, start_variant, paralog
    );

    // ── 8. Generate outputs ──
    fs::create_dir_all(&args.output)?;

    // Build rescue lookup
    let mut rescue_map: HashMap<(usize, usize), &RescueEvent> = HashMap::new();
    for r in &rescues {
        rescue_map.insert((r.cluster_idx, r.genome_idx), r);
    }

    // 8a. Rescue report TSV
    let report_path = args.output.join("rescue_report.tsv");
    let report_file = fs::File::create(&report_path)?;
    let mut report = BufWriter::new(report_file);
    writeln!(
        report,
        "gene_cluster\tgenome\trescue_type\tidentity\tcoverage\tcontig\tstart\tend\tstrand\tconfidence\tsynteny\tannotation_match"
    )?;
    for r in &rescues {
        writeln!(
            report,
            "{}\t{}\t{}\t{:.1}\t{:.1}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}",
            pa.clusters[r.cluster_idx],
            pa.genomes[r.genome_idx],
            r.rescue_type,
            r.identity,
            r.coverage,
            r.contig,
            r.start,
            r.end,
            r.strand,
            r.confidence,
            match r.synteny_ok {
                Some(true) => "PASS",
                Some(false) => "FAIL",
                None => "",
            },
            r.annotation_match.as_deref().unwrap_or(""),
        )?;
    }
    report.flush()?;

    // 8b. Corrected presence/absence CSV
    let corrected_path = args.output.join("gene_presence_absence_recore.csv");
    let corrected_file = fs::File::create(&corrected_path)?;
    let mut corrected = BufWriter::new(corrected_file);

    write!(corrected, "Gene")?;
    for g in &pa.genomes {
        write!(corrected, ",{}", g)?;
    }
    writeln!(corrected)?;

    for (ci, cluster_name) in pa.clusters.iter().enumerate() {
        write!(corrected, "{}", cluster_name)?;
        for gi in 0..pa.genomes.len() {
            if let Some(r) = rescue_map.get(&(ci, gi)) {
                write!(
                    corrected,
                    ",RESCUED:{}:{}-{}",
                    r.contig, r.start, r.end
                )?;
            } else {
                write!(corrected, ",{}", pa.cell_values[ci][gi])?;
            }
        }
        writeln!(corrected)?;
    }
    corrected.flush()?;

    // 8c. Summary
    let core_after = {
        let n = pa.genomes.len();
        let mut count = 0;
        for ci in 0..pa.clusters.len() {
            let mut all_present = true;
            for gi in 0..n {
                if !pa.present[ci].contains(&gi) && !rescue_map.contains_key(&(ci, gi)) {
                    all_present = false;
                    break;
                }
            }
            if all_present {
                count += 1;
            }
        }
        count
    };

    let soft_core_after = {
        let n = pa.genomes.len();
        let threshold = (n as f64 * 0.95).ceil() as usize;
        let mut count = 0;
        for ci in 0..pa.clusters.len() {
            let present_count = pa.present[ci].len()
                + (0..n)
                    .filter(|gi| !pa.present[ci].contains(gi) && rescue_map.contains_key(&(ci, *gi)))
                    .count();
            if present_count >= threshold {
                count += 1;
            }
        }
        count
    };

    let summary_path = args.output.join("summary.txt");
    let summary_file = fs::File::create(&summary_path)?;
    let mut summary = BufWriter::new(summary_file);
    writeln!(summary, "Recore Summary")?;
    writeln!(summary, "{}", "=".repeat(50))?;
    writeln!(summary, "Clusters:              {}", pa.clusters.len())?;
    writeln!(summary, "Genomes:               {}", pa.genomes.len())?;
    writeln!(summary, "Total gaps:            {}", all_gaps.len())?;
    writeln!(summary, "Gaps searched:         {}", filtered_gaps)?;
    writeln!(summary, "Skipped (prevalence):  {}", skipped_by_prevalence)?;
    writeln!(summary, "Rescued:               {}", rescues.len())?;
    writeln!(summary, "Remaining gaps:        {}", filtered_gaps - rescues.len())?;
    writeln!(summary)?;
    writeln!(summary, "Core genome (strict):  {} → {} (+{})", core_before, core_after, core_after - core_before)?;
    writeln!(summary, "Soft-core (≥95%):      {} → {} (+{})", pa.soft_core_size(), soft_core_after, soft_core_after - pa.soft_core_size())?;
    writeln!(summary)?;
    writeln!(summary, "Rescue types:")?;
    writeln!(summary, "  DROPPED:       {}", dropped)?;
    writeln!(summary, "  SPLIT:         {}", split)?;
    writeln!(summary, "  FRAGMENTED:    {}", fragmented)?;
    writeln!(summary, "  START_VARIANT: {}", start_variant)?;
    writeln!(summary, "  PARALOG:       {}", paralog)?;
    if args.synteny {
        let syn_pass = rescues.iter().filter(|r| r.synteny_ok == Some(true)).count();
        let syn_fail = rescues.iter().filter(|r| r.synteny_ok == Some(false)).count();
        writeln!(summary)?;
        writeln!(summary, "Synteny validation:")?;
        writeln!(summary, "  PASS:  {}", syn_pass)?;
        writeln!(summary, "  FAIL:  {}", syn_fail)?;
    }
    summary.flush()?;

    info!(
        "Core genome: {} → {} (+{})",
        core_before,
        core_after,
        core_after - core_before
    );
    info!(
        "Soft-core (≥95%): {} → {} (+{})",
        pa.soft_core_size(),
        soft_core_after,
        soft_core_after - pa.soft_core_size()
    );
    info!("Output: {}", args.output.display());

    Ok(())
}
