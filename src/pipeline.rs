use ahash::{HashMap, HashMapExt};
use anyhow::Result;
use log::{info, warn};
use std::fs;
use std::io::{BufWriter, Write};

use crate::cli::Args;
use crate::parser::{
    find_genome_files, parse_gff_genes, parse_pan_reference, parse_presence_absence,
};
use crate::rescue::{
    build_gene_to_cluster, classify_hits, resolve_conflicts, validate_synteny,
    GffGene, GffIndex, RescueEvent, RescueType,
};
use crate::search::batch_search;

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

    // ── 2. Filter gaps by prevalence (no intermediate Vec) ──
    let min_prevalence = args.min_prevalence;
    let mut gaps_by_genome: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut total_gaps = 0usize;
    let mut filtered_gaps = 0usize;
    let mut skipped_by_prevalence = 0usize;

    let n_genomes = pa.genomes.len();
    for (ci, present_set) in pa.present.iter().enumerate() {
        let prevalence = present_set.len() as f64 / n_genomes as f64;
        for gi in 0..n_genomes {
            if !present_set.contains(&gi) {
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
        core_before,
        pa.clusters.len()
    );

    // Check genome file availability
    let mut searchable_genomes: Vec<(usize, Vec<usize>)> = Vec::new();
    let mut missing_genomes = 0usize;
    for (gi, cluster_indices) in gaps_by_genome {
        let genome_name = &pa.genomes[gi];
        if genome_files.contains_key(genome_name) {
            searchable_genomes.push((gi, cluster_indices));
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
        for (gi, genome_name) in pa.genomes.iter().enumerate() {
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

    // ── 4. Build gene_to_cluster maps ONCE (reused in rescue + synteny) ──
    // Only build for genomes that have GFF or have gaps (need for classify_hits)
    let gene_to_cluster_maps: HashMap<usize, HashMap<String, usize>> = {
        let mut maps = HashMap::new();
        // Genomes with gaps (for classify_hits)
        for (gi, _) in &searchable_genomes {
            maps.entry(*gi)
                .or_insert_with(|| build_gene_to_cluster(&pa, *gi));
        }
        // Genomes with GFF (for synteny validation)
        if args.synteny {
            for &gi in gff_indices.keys() {
                maps.entry(gi)
                    .or_insert_with(|| build_gene_to_cluster(&pa, gi));
            }
        }
        maps
    };

    // ── 5. Batch search: ONE minimap2 invocation for all genomes ──
    let hits_map = batch_search(
        &pan_ref,
        &genome_files,
        &searchable_genomes,
        &pa.genomes,
        &pa.clusters,
        args.threads,
        args.min_identity,
        &args.preset,
    )?;

    // ── 6. Classify hits into rescue events ──
    let mut rescues: Vec<RescueEvent> = Vec::new();

    for ((gi, cluster_name), hits) in &hits_map {
        let Some(&ci) = cluster_to_idx.get(cluster_name.as_str()) else {
            continue;
        };
        let gff_index = gff_indices.get(gi);
        let gene_to_cluster = gene_to_cluster_maps
            .get(gi)
            .map(|m| m as &HashMap<String, usize>);

        // Empty map fallback (genome has no GFF)
        let empty_map = HashMap::new();
        let g2c = gene_to_cluster.unwrap_or(&empty_map);

        if let Some(rescue) = classify_hits(
            hits,
            ci,
            *gi,
            &pa,
            gff_index,
            g2c,
            args.min_coverage,
        ) {
            rescues.push(rescue);
        }
    }

    info!("Raw rescues before filtering: {}", rescues.len());

    // ── 7. Synteny validation (if requested) ──
    if args.synteny && !gff_indices.is_empty() {
        info!("Validating synteny...");
        for rescue in &mut rescues {
            validate_synteny(rescue, &pa, &gff_indices, &gene_to_cluster_maps);
        }
        let validated = rescues.iter().filter(|r| r.synteny_ok == Some(true)).count();
        let failed = rescues.iter().filter(|r| r.synteny_ok == Some(false)).count();
        info!("Synteny: {} validated, {} failed", validated, failed);
    }

    // ── 8. Conflict resolution ──
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

    // ── 9. Count rescue types ──
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

    // ── 10. Generate outputs ──
    fs::create_dir_all(&args.output)?;

    // Build rescue lookup
    let rescue_map: HashMap<(usize, usize), &RescueEvent> = rescues
        .iter()
        .map(|r| ((r.cluster_idx, r.genome_idx), r))
        .collect();

    write_report(&rescues, &pa, &args.output)?;
    write_corrected_csv(&pa, &rescue_map, &args.output)?;

    let stats = SummaryStats {
        core_before,
        total_gaps,
        filtered_gaps,
        skipped_by_prevalence,
        dropped,
        split,
        fragmented,
        start_variant,
        paralog,
        synteny_enabled: args.synteny,
    };
    write_summary(&pa, &rescue_map, &rescues, &stats, &args.output)?;

    // Log final stats
    let core_after = count_core_after(&pa, &rescue_map);
    let soft_core_after = count_soft_core_after(&pa, &rescue_map);
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

// ── Output helpers ──

fn write_report(
    rescues: &[RescueEvent],
    pa: &crate::parser::PresenceAbsence,
    output_dir: &std::path::Path,
) -> Result<()> {
    let path = output_dir.join("rescue_report.tsv");
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "gene_cluster\tgenome\trescue_type\tidentity\tcoverage\tcontig\tstart\tend\tstrand\tconfidence\tsynteny\tannotation_match"
    )?;
    for r in rescues {
        writeln!(
            w,
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
    w.flush()?;
    Ok(())
}

fn write_corrected_csv(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(usize, usize), &RescueEvent>,
    output_dir: &std::path::Path,
) -> Result<()> {
    let path = output_dir.join("gene_presence_absence_recore.csv");
    let file = fs::File::create(&path)?;
    let writer = BufWriter::new(file);
    let mut csv_writer = csv::Writer::from_writer(writer);

    csv_writer.write_record(&pa.header)?;

    let genome_idx_set: ahash::HashSet<usize> = pa.genome_col_indices.iter().copied().collect();
    let meta_col_indices: Vec<usize> = (0..pa.header.len())
        .filter(|i| !genome_idx_set.contains(i))
        .collect();

    let isolates_col = pa.header.iter().position(|h| h == "No. isolates");
    let sequences_col = pa.header.iter().position(|h| h == "No. sequences");
    let avg_seq_col = pa.header.iter().position(|h| h == "Avg sequences per isolate");

    // Pre-allocate row buffer (reuse across iterations to reduce allocations)
    let n_cols = pa.header.len();
    let mut row: Vec<String> = vec![String::new(); n_cols];
    // Pre-allocate a rescue format buffer
    let mut rescue_buf = String::new();

    for ci in 0..pa.clusters.len() {
        // Fill metadata columns (by reference where possible, clone only when needed)
        for (mi, &col_idx) in meta_col_indices.iter().enumerate() {
            if mi < pa.metadata[ci].len() {
                row[col_idx].clone_from(&pa.metadata[ci][mi]);
            } else {
                row[col_idx].clear();
            }
        }

        // Fill genome columns
        let mut rescue_count = 0usize;
        for (gi, &col_idx) in pa.genome_col_indices.iter().enumerate() {
            if let Some(r) = rescue_map.get(&(ci, gi)) {
                rescue_buf.clear();
                use std::fmt::Write as FmtWrite;
                write!(rescue_buf, "RESCUED:{}:{}-{}", r.contig, r.start, r.end).unwrap();
                row[col_idx].clone_from(&rescue_buf);
                rescue_count += 1;
            } else {
                row[col_idx].clone_from(&pa.cell_values[ci][gi]);
            }
        }

        // Recalculate counts if rescues happened
        if rescue_count > 0 {
            let new_count = pa.present[ci].len() + rescue_count;
            if let Some(col) = isolates_col {
                row[col].clear();
                use std::fmt::Write as FmtWrite;
                write!(row[col], "{}", new_count).unwrap();
            }
            if let Some(col) = sequences_col {
                row[col].clear();
                use std::fmt::Write as FmtWrite;
                write!(row[col], "{}", new_count).unwrap();
            }
            if let Some(col) = avg_seq_col {
                row[col].clear();
                use std::fmt::Write as FmtWrite;
                // After rescue, each isolate has exactly 1 sequence
                write!(row[col], "1.00").unwrap();
            }
        }

        // Write row as references (no clone into csv writer)
        csv_writer.write_record(row.iter().map(|s| s.as_str()))?;
    }
    csv_writer.flush()?;
    Ok(())
}

fn count_core_after(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(usize, usize), &RescueEvent>,
) -> usize {
    let n = pa.genomes.len();
    let mut count = 0;
    for ci in 0..pa.clusters.len() {
        let present = pa.present[ci].len();
        if present == n {
            count += 1;
            continue;
        }
        // Only check rescue map if close to core
        let rescued = (0..n)
            .filter(|gi| !pa.present[ci].contains(gi) && rescue_map.contains_key(&(ci, *gi)))
            .count();
        if present + rescued == n {
            count += 1;
        }
    }
    count
}

fn count_soft_core_after(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(usize, usize), &RescueEvent>,
) -> usize {
    let n = pa.genomes.len();
    let threshold = (n as f64 * 0.95).ceil() as usize;
    let mut count = 0;
    for ci in 0..pa.clusters.len() {
        let present = pa.present[ci].len();
        if present >= threshold {
            count += 1;
            continue;
        }
        let rescued = (0..n)
            .filter(|gi| !pa.present[ci].contains(gi) && rescue_map.contains_key(&(ci, *gi)))
            .count();
        if present + rescued >= threshold {
            count += 1;
        }
    }
    count
}

struct SummaryStats {
    core_before: usize,
    total_gaps: usize,
    filtered_gaps: usize,
    skipped_by_prevalence: usize,
    dropped: usize,
    split: usize,
    fragmented: usize,
    start_variant: usize,
    paralog: usize,
    synteny_enabled: bool,
}

fn write_summary(
    pa: &crate::parser::PresenceAbsence,
    rescue_map: &HashMap<(usize, usize), &RescueEvent>,
    rescues: &[RescueEvent],
    stats: &SummaryStats,
    output_dir: &std::path::Path,
) -> Result<()> {
    let core_after = count_core_after(pa, rescue_map);
    let soft_core_before = pa.soft_core_size();
    let soft_core_after = count_soft_core_after(pa, rescue_map);

    let path = output_dir.join("summary.txt");
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);
    writeln!(w, "Recore Summary")?;
    writeln!(w, "{}", "=".repeat(50))?;
    writeln!(w, "Clusters:              {}", pa.clusters.len())?;
    writeln!(w, "Genomes:               {}", pa.genomes.len())?;
    writeln!(w, "Total gaps:            {}", stats.total_gaps)?;
    writeln!(w, "Gaps searched:         {}", stats.filtered_gaps)?;
    writeln!(w, "Skipped (prevalence):  {}", stats.skipped_by_prevalence)?;
    writeln!(w, "Rescued:               {}", rescues.len())?;
    writeln!(w, "Remaining gaps:        {}", stats.filtered_gaps - rescues.len())?;
    writeln!(w)?;
    writeln!(w, "Core genome (strict):  {} → {} (+{})", stats.core_before, core_after, core_after - stats.core_before)?;
    writeln!(w, "Soft-core (≥95%):      {} → {} (+{})", soft_core_before, soft_core_after, soft_core_after - soft_core_before)?;
    writeln!(w)?;
    writeln!(w, "Rescue types:")?;
    writeln!(w, "  DROPPED:       {}", stats.dropped)?;
    writeln!(w, "  SPLIT:         {}", stats.split)?;
    writeln!(w, "  FRAGMENTED:    {}", stats.fragmented)?;
    writeln!(w, "  START_VARIANT: {}", stats.start_variant)?;
    writeln!(w, "  PARALOG:       {}", stats.paralog)?;
    if stats.synteny_enabled {
        let syn_pass = rescues.iter().filter(|r| r.synteny_ok == Some(true)).count();
        let syn_fail = rescues.iter().filter(|r| r.synteny_ok == Some(false)).count();
        writeln!(w)?;
        writeln!(w, "Synteny validation:")?;
        writeln!(w, "  PASS:  {}", syn_pass)?;
        writeln!(w, "  FAIL:  {}", syn_fail)?;
    }
    w.flush()?;
    Ok(())
}
