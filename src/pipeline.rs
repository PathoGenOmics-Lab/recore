use ahash::{HashMap, HashMapExt};
use anyhow::Result;
use log::{info, warn};
use rayon::prelude::*;
use std::fs;
use std::io::Write;
use std::sync::Mutex;

use crate::cli::Args;
use crate::parser::{find_genome_files, parse_gff_genes, parse_pan_reference, parse_presence_absence};
use crate::rescue::{classify_hit, RescueEvent, RescueType};
use crate::search::search_genome;

pub fn run(args: &Args) -> Result<()> {
    // ── Parse inputs ──
    let pa = parse_presence_absence(&args.presence_absence)?;
    let pan_ref = parse_pan_reference(&args.pan_reference)?;
    let genome_files = find_genome_files(&args.genomes)?;

    // Find gaps grouped by genome
    let gaps = pa.gaps();
    let mut gaps_by_genome: HashMap<usize, Vec<usize>> = HashMap::new();
    for &(ci, gi) in &gaps {
        gaps_by_genome.entry(gi).or_default().push(ci);
    }

    info!(
        "Total gaps: {} across {} genomes (core before: {})",
        gaps.len(),
        gaps_by_genome.len(),
        pa.core_size()
    );

    // Check genome file availability
    let mut searchable_genomes: Vec<(usize, Vec<usize>)> = Vec::new();
    for (gi, cluster_indices) in &gaps_by_genome {
        let genome_name = &pa.genomes[*gi];
        if genome_files.contains_key(genome_name) {
            searchable_genomes.push((*gi, cluster_indices.clone()));
        } else {
            warn!("No FASTA file for genome '{}', skipping", genome_name);
        }
    }

    // Load GFF annotations if provided
    let gff_data: Option<HashMap<usize, Vec<(String, u64, u64, String)>>> =
        if let Some(ref ann_dir) = args.annotations {
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
                    Ok(genes) => {
                        data.insert(gi, genes);
                    }
                    Err(e) => {
                        warn!("Cannot parse GFF for {}: {}", genome_name, e);
                    }
                }
            }
            info!("Loaded GFF annotations for {} genomes", data.len());
            Some(data)
        } else {
            None
        };

    // ── Search and rescue (parallel per genome) ──
    let all_rescues: Mutex<Vec<RescueEvent>> = Mutex::new(Vec::new());
    let progress: Mutex<usize> = Mutex::new(0);
    let total = searchable_genomes.len();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .ok(); // ignore if already set

    searchable_genomes
        .par_iter()
        .for_each(|(gi, cluster_indices)| {
            let genome_name = &pa.genomes[*gi];
            let genome_fasta = &genome_files[genome_name];

            // Build query set for this genome
            let mut queries: HashMap<String, Vec<u8>> = HashMap::new();
            for &ci in cluster_indices {
                let cluster_name = &pa.clusters[ci];
                if let Some(seq) = pan_ref.get(cluster_name) {
                    queries.insert(cluster_name.clone(), seq.clone());
                }
            }

            if queries.is_empty() {
                return;
            }

            // Search with minimap2 (1 thread per genome, parallelism at genome level)
            let hits = match search_genome(&queries, genome_fasta, 1, args.min_identity, args.min_coverage) {
                Ok(h) => h,
                Err(e) => {
                    warn!("minimap2 failed for {}: {}", genome_name, e);
                    return;
                }
            };

            // Classify hits
            let gff_genes = gff_data.as_ref().and_then(|d| d.get(gi));
            let mut rescues = Vec::new();
            for hit in &hits {
                let ci = pa
                    .clusters
                    .iter()
                    .position(|c| c == &hit.query)
                    .unwrap();
                let rescue = classify_hit(hit, ci, *gi, &pa, gff_genes.map(|v| v.as_slice()));
                rescues.push(rescue);
            }

            if !rescues.is_empty() {
                all_rescues.lock().unwrap().extend(rescues);
            }

            // Progress
            let mut count = progress.lock().unwrap();
            *count += 1;
            if *count % 50 == 0 || *count == total {
                info!("Processed {}/{} genomes", *count, total);
            }
        });

    let rescues = all_rescues.into_inner().unwrap();
    info!("Total rescues: {}", rescues.len());

    // ── Count rescue types ──
    let mut dropped = 0usize;
    let mut split = 0usize;
    let mut overlap = 0usize;
    for r in &rescues {
        match &r.rescue_type {
            RescueType::Dropped => dropped += 1,
            RescueType::Split { .. } => split += 1,
            RescueType::Overlap { .. } => overlap += 1,
        }
    }
    info!("  DROPPED: {}, SPLIT: {}, OVERLAP: {}", dropped, split, overlap);

    // ── Generate outputs ──
    fs::create_dir_all(&args.output)?;

    // 1. Rescue report TSV
    let report_path = args.output.join("rescue_report.tsv");
    let mut report = fs::File::create(&report_path)?;
    writeln!(
        report,
        "gene_cluster\tgenome\trescue_type\tidentity\tcoverage\tcontig\tstart\tend\tstrand\tannotation_match"
    )?;
    for r in &rescues {
        writeln!(
            report,
            "{}\t{}\t{}\t{:.1}\t{:.1}\t{}\t{}\t{}\t{}\t{}",
            pa.clusters[r.cluster_idx],
            pa.genomes[r.genome_idx],
            r.rescue_type,
            r.identity,
            r.coverage,
            r.contig,
            r.start,
            r.end,
            r.strand,
            r.annotation_match.as_deref().unwrap_or(""),
        )?;
    }

    // 2. Corrected presence/absence CSV
    let corrected_path = args.output.join("gene_presence_absence_recore.csv");
    let mut corrected = fs::File::create(&corrected_path)?;

    // Header
    write!(corrected, "Gene")?;
    for g in &pa.genomes {
        write!(corrected, ",{}", g)?;
    }
    writeln!(corrected)?;

    // Build rescue lookup: (cluster_idx, genome_idx) -> rescue info
    let mut rescue_map: HashMap<(usize, usize), String> = HashMap::new();
    for r in &rescues {
        rescue_map.insert(
            (r.cluster_idx, r.genome_idx),
            format!("RESCUED:{}:{}-{}", r.contig, r.start, r.end),
        );
    }

    for (ci, cluster_name) in pa.clusters.iter().enumerate() {
        write!(corrected, "{}", cluster_name)?;
        for gi in 0..pa.genomes.len() {
            let val = if let Some(rescued) = rescue_map.get(&(ci, gi)) {
                rescued.as_str()
            } else {
                pa.cell_values[ci][gi].as_str()
            };
            write!(corrected, ",{}", val)?;
        }
        writeln!(corrected)?;
    }

    // 3. Summary
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

    let summary_path = args.output.join("summary.txt");
    let mut summary = fs::File::create(&summary_path)?;
    writeln!(summary, "Recore Summary")?;
    writeln!(summary, "{}", "=".repeat(40))?;
    writeln!(summary, "Clusters:           {}", pa.clusters.len())?;
    writeln!(summary, "Genomes:            {}", pa.genomes.len())?;
    writeln!(summary, "Total gaps:         {}", gaps.len())?;
    writeln!(summary, "Rescued:            {}", rescues.len())?;
    writeln!(summary, "Remaining gaps:     {}", gaps.len() - rescues.len())?;
    writeln!(summary, "Core size (before): {}", pa.core_size())?;
    writeln!(summary, "Core size (after):  {} (+{})", core_after, core_after - pa.core_size())?;
    writeln!(summary)?;
    writeln!(summary, "Rescue types:")?;
    writeln!(summary, "  DROPPED: {}", dropped)?;
    writeln!(summary, "  SPLIT:   {}", split)?;
    writeln!(summary, "  OVERLAP: {}", overlap)?;

    info!("Core genome: {} → {} (+{})", pa.core_size(), core_after, core_after - pa.core_size());
    info!("Output: {}", args.output.display());

    Ok(())
}
