use ahash::HashMap;
use anyhow::{Context, Result};
use log::debug;
use std::io::Write;
use std::path::Path;
use std::process::Command;
use tempfile::NamedTempFile;

/// A single minimap2 hit (PAF format).
#[derive(Debug, Clone)]
pub struct Hit {
    pub query: String,
    pub query_length: u64,
    pub query_start: u64,
    pub query_end: u64,
    pub strand: char,
    pub target: String,
    pub target_start: u64,
    pub target_end: u64,
    pub matches: u64,
    pub alignment_length: u64,
    #[allow(dead_code)]
    pub mapq: u32,
}

impl Hit {
    /// Alignment identity as percentage.
    pub fn identity(&self) -> f64 {
        if self.alignment_length == 0 {
            return 0.0;
        }
        (self.matches as f64 / self.alignment_length as f64) * 100.0
    }

    /// Query coverage as percentage.
    pub fn query_coverage(&self) -> f64 {
        if self.query_length == 0 {
            return 0.0;
        }
        ((self.query_end - self.query_start) as f64 / self.query_length as f64) * 100.0
    }
}

/// Run minimap2 on a set of query sequences against a genome.
/// Returns the best hit per query that passes thresholds.
pub fn search_genome(
    queries: &HashMap<String, Vec<u8>>,
    genome_fasta: &Path,
    threads: usize,
    min_identity: f64,
    min_coverage: f64,
) -> Result<Vec<Hit>> {
    if queries.is_empty() {
        return Ok(Vec::new());
    }

    // Write queries to temp FASTA
    let mut query_file = NamedTempFile::new().context("Cannot create temp file for queries")?;
    for (name, seq) in queries {
        writeln!(query_file, ">{}", name)?;
        query_file.write_all(seq)?;
        writeln!(query_file)?;
    }
    query_file.flush()?;

    // Run minimap2
    let output = Command::new("minimap2")
        .args([
            "-x", "asm20",
            "-t", &threads.to_string(),
            "--secondary=yes",
            "-N", "5",
            genome_fasta.to_str().unwrap(),
            query_file.path().to_str().unwrap(),
        ])
        .output()
        .context("Failed to run minimap2 — is it installed?")?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        anyhow::bail!("minimap2 failed: {}", stderr);
    }

    // Parse PAF output
    let stdout = String::from_utf8_lossy(&output.stdout);
    let mut best_per_query: HashMap<String, (Hit, f64)> = HashMap::default();

    for line in stdout.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            continue;
        }

        let hit = Hit {
            query: fields[0].to_string(),
            query_length: fields[1].parse().unwrap_or(0),
            query_start: fields[2].parse().unwrap_or(0),
            query_end: fields[3].parse().unwrap_or(0),
            strand: fields[4].chars().next().unwrap_or('+'),
            target: fields[5].to_string(),
            target_start: fields[7].parse().unwrap_or(0),
            target_end: fields[8].parse().unwrap_or(0),
            matches: fields[9].parse().unwrap_or(0),
            alignment_length: fields[10].parse().unwrap_or(0),
            mapq: fields[11].parse().unwrap_or(0),
        };

        let ident = hit.identity();
        let cov = hit.query_coverage();

        if ident < min_identity || cov < min_coverage {
            continue;
        }

        let score = ident * cov;
        let entry = best_per_query.get(&hit.query);
        if entry.is_none() || score > entry.unwrap().1 {
            best_per_query.insert(hit.query.clone(), (hit, score));
        }
    }

    let hits: Vec<Hit> = best_per_query.into_values().map(|(h, _)| h).collect();
    debug!(
        "minimap2: {} queries → {} hits passing filters",
        queries.len(),
        hits.len()
    );

    Ok(hits)
}
