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
    #[allow(dead_code)]
    pub query: String,       // allele query name
    #[allow(dead_code)]
    pub cluster: String,     // original cluster name (without allele suffix)
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

/// Run minimap2 and return ALL hits per cluster that pass identity threshold.
/// Multi-allele queries are named "cluster__0", "cluster__1", etc.
/// Hits are grouped back by original cluster name.
pub fn search_genome(
    queries: &HashMap<String, Vec<&[u8]>>,
    genome_fasta: &Path,
    threads: usize,
    min_identity: f64,
    preset: &str,
) -> Result<HashMap<String, Vec<Hit>>> {
    if queries.is_empty() {
        return Ok(HashMap::default());
    }

    // Write queries to temp FASTA with allele-indexed names
    let mut query_file = NamedTempFile::new().context("Cannot create temp file for queries")?;
    let mut name_to_cluster: HashMap<String, String> = HashMap::default();

    for (cluster_name, alleles) in queries {
        for (ai, seq) in alleles.iter().enumerate() {
            let query_name = if alleles.len() == 1 {
                cluster_name.clone()
            } else {
                format!("{}__{}", cluster_name, ai)
            };
            name_to_cluster.insert(query_name.clone(), cluster_name.clone());
            writeln!(query_file, ">{}", query_name)?;
            query_file.write_all(seq)?;
            writeln!(query_file)?;
        }
    }
    query_file.flush()?;

    let output = Command::new("minimap2")
        .args([
            "-x", preset,
            "-t", &threads.to_string(),
            "--secondary=yes",
            "-N", "50",
            "-p", "0.5",
            &genome_fasta.to_string_lossy(),
            &query_file.path().to_string_lossy(),
        ])
        .output()
        .context("Failed to run minimap2 — is it installed?")?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        anyhow::bail!("minimap2 failed: {}", stderr);
    }

    // Log minimap2 warnings even on success
    if !output.stderr.is_empty() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        for line in stderr.lines() {
            if line.contains("[WARNING]") || line.contains("[W::") {
                log::warn!("minimap2: {}", line);
            }
        }
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    let mut hits_per_cluster: HashMap<String, Vec<Hit>> = HashMap::default();

    for line in stdout.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            continue;
        }

        let query_name = fields[0].to_string();
        let cluster_name = name_to_cluster
            .get(&query_name)
            .cloned()
            .unwrap_or_else(|| query_name.clone());

        let hit = Hit {
            query: query_name,
            cluster: cluster_name.clone(),
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

        if hit.identity() < min_identity {
            continue;
        }

        hits_per_cluster
            .entry(cluster_name)
            .or_default()
            .push(hit);
    }

    debug!(
        "minimap2: {} queries → {} clusters with hits",
        queries.len(),
        hits_per_cluster.len()
    );

    Ok(hits_per_cluster)
}
