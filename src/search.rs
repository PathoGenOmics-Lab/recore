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

/// Run minimap2 and return ALL hits per query that pass identity threshold.
/// Coverage filtering happens downstream (we need partial hits for fragment detection).
pub fn search_genome(
    queries: &HashMap<String, &[u8]>,
    genome_fasta: &Path,
    threads: usize,
    min_identity: f64,
    preset: &str,
) -> Result<HashMap<String, Vec<Hit>>> {
    if queries.is_empty() {
        return Ok(HashMap::default());
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
            "-x", preset,
            "-t", &threads.to_string(),
            "--secondary=yes",
            "-N", "50",  // generous: catch paralogs and fragments
            "-p", "0.5", // report secondary if score >= 50% of primary
            genome_fasta.to_str().unwrap(),
            query_file.path().to_str().unwrap(),
        ])
        .output()
        .context("Failed to run minimap2 — is it installed?")?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        anyhow::bail!("minimap2 failed: {}", stderr);
    }

    // Parse PAF output — collect ALL hits per query
    let stdout = String::from_utf8_lossy(&output.stdout);
    let mut hits_per_query: HashMap<String, Vec<Hit>> = HashMap::default();

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

        // Only filter by identity; coverage filtering is context-dependent
        if hit.identity() < min_identity {
            continue;
        }

        hits_per_query
            .entry(hit.query.clone())
            .or_default()
            .push(hit);
    }

    debug!(
        "minimap2: {} queries → {} with hits",
        queries.len(),
        hits_per_query.len()
    );

    Ok(hits_per_query)
}
