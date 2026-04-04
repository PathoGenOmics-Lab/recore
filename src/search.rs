use ahash::HashMap;
use anyhow::{Context, Result};
use log::{debug, info};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::Command;
use tempfile::NamedTempFile;

/// Separator between genome name and original contig name in concatenated targets.
/// Chosen to be extremely unlikely in real contig names.
const CONTIG_SEP: &str = "|||";

/// A single minimap2 hit (PAF format).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct Hit {
    pub query: String,       // allele query name
    pub cluster: String,     // original cluster name (without allele suffix)
    pub query_length: u64,
    pub query_start: u64,
    pub query_end: u64,
    pub strand: char,
    pub target: String,      // original contig name (genome prefix stripped)
    pub target_start: u64,
    pub target_end: u64,
    pub matches: u64,
    pub alignment_length: u64,
    pub mapq: u32,
}

impl Hit {
    pub fn identity(&self) -> f64 {
        if self.alignment_length == 0 {
            return 0.0;
        }
        (self.matches as f64 / self.alignment_length as f64) * 100.0
    }

    pub fn query_coverage(&self) -> f64 {
        if self.query_length == 0 {
            return 0.0;
        }
        ((self.query_end - self.query_start) as f64 / self.query_length as f64) * 100.0
    }

    /// Combined score for ranking: identity × coverage / 10000.
    pub fn score(&self) -> f64 {
        self.identity() * self.query_coverage() / 10000.0
    }
}

#[allow(clippy::too_many_arguments)]
/// Batch search: run ONE minimap2 invocation against ALL genomes.
///
/// Strategy:
/// 1. Write all unique query sequences to a single temp FASTA
/// 2. Concatenate all genome FASTAs into one target, prefixing contig names
///    with `genome_name|||original_contig` to avoid collisions
/// 3. Run minimap2 once (multi-threaded)
/// 4. Demux hits back to per-(genome, cluster) groups
/// 5. Filter: only return hits for (cluster, genome) pairs that are actual gaps
///
/// Returns: HashMap<(genome_idx, cluster_name), Vec<Hit>>
pub fn batch_search(
    queries: &HashMap<String, Vec<Vec<u8>>>,    // cluster → allele sequences (owned)
    genome_files: &HashMap<String, std::path::PathBuf>,  // genome_name → fasta path
    searchable_genomes: &[(usize, Vec<usize>)], // (genome_idx, cluster_indices_to_search)
    genome_names: &[String],                     // pa.genomes
    cluster_names: &[String],                    // pa.clusters
    threads: usize,
    min_identity: f64,
    preset: &str,
) -> Result<HashMap<(usize, String), Vec<Hit>>> {
    if searchable_genomes.is_empty() {
        return Ok(HashMap::default());
    }

    // ── 1. Collect unique query clusters needed across all genomes ──
    let mut needed_clusters: ahash::HashSet<usize> = ahash::HashSet::default();
    for (_gi, cluster_indices) in searchable_genomes {
        for &ci in cluster_indices {
            needed_clusters.insert(ci);
        }
    }

    // ── 2. Write query FASTA ──
    let mut query_file = NamedTempFile::new().context("Cannot create temp query file")?;
    let mut name_to_cluster: HashMap<String, String> = HashMap::default();
    let mut total_queries = 0usize;

    for &ci in &needed_clusters {
        let cluster_name = &cluster_names[ci];
        let Some(alleles) = queries.get(cluster_name) else {
            continue;
        };
        for (ai, seq) in alleles.iter().enumerate() {
            let query_name = if alleles.len() == 1 {
                cluster_name.clone()
            } else {
                format!("{}~~~ALLELE~{}", cluster_name, ai)
            };
            name_to_cluster.insert(query_name.clone(), cluster_name.clone());
            writeln!(query_file, ">{}", query_name)?;
            query_file.write_all(seq)?;
            writeln!(query_file)?;
            total_queries += 1;
        }
    }
    query_file.flush()?;

    if total_queries == 0 {
        return Ok(HashMap::default());
    }

    // ── 3. Concatenate genome FASTAs with prefixed contig names ──
    let mut target_file = NamedTempFile::new().context("Cannot create temp target file")?;
    {
        let mut writer = BufWriter::new(&mut target_file);
        let mut genomes_written = 0usize;

        // Only include genomes that actually have gaps to search
        let genome_indices_needed: ahash::HashSet<usize> = searchable_genomes
            .iter()
            .map(|(gi, _)| *gi)
            .collect();

        for &gi in &genome_indices_needed {
            let genome_name = &genome_names[gi];
            let Some(fasta_path) = genome_files.get(genome_name) else {
                continue;
            };

            let file = std::fs::File::open(fasta_path)
                .with_context(|| format!("Cannot open genome {}", fasta_path.display()))?;
            let reader = BufReader::new(file);

            for line in reader.lines() {
                let line = line?;
                if let Some(header) = line.strip_prefix('>') {
                    // Prefix contig name: genome_name|||original_header
                    writeln!(writer, ">{}{}{}", genome_name, CONTIG_SEP, header)?;
                } else {
                    writer.write_all(line.as_bytes())?;
                    writeln!(writer)?;
                }
            }
            genomes_written += 1;
        }
        writer.flush()?;
        info!(
            "Prepared batch search: {} queries × {} genomes",
            total_queries, genomes_written
        );
    }

    // ── 4. Run minimap2 ONCE ──
    let output = Command::new("minimap2")
        .args([
            "-x", preset,
            "-t", &threads.to_string(),
            "--secondary=yes",
            "-N", "50",
            "-p", "0.5",
            &target_file.path().to_string_lossy(),
            &query_file.path().to_string_lossy(),
        ])
        .output()
        .context("Failed to run minimap2 — is it installed?")?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        anyhow::bail!("minimap2 failed: {}", stderr);
    }

    if !output.stderr.is_empty() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        for line in stderr.lines() {
            if line.contains("[WARNING]") || line.contains("[W::") {
                log::warn!("minimap2: {}", line);
            }
        }
    }

    // ── 5. Build gap lookup: which (genome_idx, cluster_idx) pairs are gaps? ──
    let mut gap_set: ahash::HashSet<(usize, usize)> = ahash::HashSet::default();
    for (gi, cluster_indices) in searchable_genomes {
        for &ci in cluster_indices {
            gap_set.insert((*gi, ci));
        }
    }

    // genome_name → genome_idx for O(1) lookup during PAF parsing
    let genome_name_to_idx: HashMap<&str, usize> = genome_names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();

    // cluster_name → cluster_idx for gap filtering
    let cluster_name_to_idx: HashMap<&str, usize> = cluster_names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();

    // ── 6. Parse PAF and demux hits ──
    let stdout = String::from_utf8_lossy(&output.stdout);
    let mut hits: HashMap<(usize, String), Vec<Hit>> = HashMap::default();
    let mut total_hits = 0usize;
    let mut filtered_hits = 0usize;

    for line in stdout.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            continue;
        }
        total_hits += 1;

        // Parse target: "genome_name|||original_contig rest_of_header"
        let target_field = fields[5];
        let Some(sep_pos) = target_field.find(CONTIG_SEP) else {
            continue;
        };
        let genome_name = &target_field[..sep_pos];
        // Original contig: strip prefix, take first whitespace-delimited token
        let original_header = &target_field[sep_pos + CONTIG_SEP.len()..];
        let original_contig = original_header
            .split_whitespace()
            .next()
            .unwrap_or(original_header);

        // Resolve genome index
        let Some(&gi) = genome_name_to_idx.get(genome_name) else {
            continue;
        };

        // Resolve cluster
        let query_name = fields[0];
        let cluster_name = name_to_cluster
            .get(query_name)
            .map(|s| s.as_str())
            .unwrap_or(query_name);

        // Filter: only keep hits for actual gaps
        if let Some(&ci) = cluster_name_to_idx.get(cluster_name) {
            if !gap_set.contains(&(gi, ci)) {
                continue; // gene already present in this genome, skip
            }
        }

        let hit = Hit {
            query: query_name.to_string(),
            cluster: cluster_name.to_string(),
            query_length: fields[1].parse().unwrap_or(0),
            query_start: fields[2].parse().unwrap_or(0),
            query_end: fields[3].parse().unwrap_or(0),
            strand: fields[4].chars().next().unwrap_or('+'),
            target: original_contig.to_string(),
            target_start: fields[7].parse().unwrap_or(0),
            target_end: fields[8].parse().unwrap_or(0),
            matches: fields[9].parse().unwrap_or(0),
            alignment_length: fields[10].parse().unwrap_or(0),
            mapq: fields[11].parse().unwrap_or(0),
        };

        if hit.identity() < min_identity {
            continue;
        }

        filtered_hits += 1;
        hits.entry((gi, cluster_name.to_string()))
            .or_default()
            .push(hit);
    }

    debug!(
        "minimap2 batch: {} total hits → {} passing filters → {} (genome, cluster) pairs",
        total_hits,
        filtered_hits,
        hits.len()
    );

    Ok(hits)
}
