use ahash::HashMap;
use anyhow::{Context, Result};
use log::{debug, info};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::{Command, Stdio};
use tempfile::NamedTempFile;

/// Separator between genome name and original contig name in concatenated targets.
const CONTIG_SEP: &str = "|||";

/// A single minimap2 hit (PAF format).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct Hit {
    pub query: String,   // allele query name
    pub cluster: String, // original cluster name (without allele suffix)
    pub query_length: u64,
    pub query_start: u64,
    pub query_end: u64,
    pub strand: char,
    pub target: String, // original contig name (genome prefix stripped)
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
/// 3. Pre-build minimap2 index (.mmi) — reusable, avoids re-indexing
/// 4. Run minimap2 alignment with piped stdout (streaming PAF, no RAM spike)
/// 5. Demux hits back to per-(genome, cluster) groups
/// 6. Filter: only return hits for (cluster, genome) pairs that are actual gaps
///
/// Returns: HashMap<(genome_idx, cluster_idx), Vec<Hit>>
pub fn batch_search(
    queries: &HashMap<String, Vec<Vec<u8>>>,
    genome_files: &HashMap<String, std::path::PathBuf>,
    searchable_genomes: &[(usize, Vec<usize>)],
    genome_names: &[String],
    cluster_names: &[String],
    threads: usize,
    min_identity: f64,
    preset: &str,
) -> Result<HashMap<(usize, usize), Vec<Hit>>> {
    if searchable_genomes.is_empty() {
        return Ok(HashMap::default());
    }

    // ── 1. Collect unique query clusters needed across all genomes ──
    let mut needed_cluster_set: ahash::HashSet<usize> = ahash::HashSet::default();
    for (_gi, cluster_indices) in searchable_genomes {
        for &ci in cluster_indices {
            needed_cluster_set.insert(ci);
        }
    }
    let mut needed_clusters: Vec<usize> = needed_cluster_set.into_iter().collect();
    needed_clusters.sort_unstable();

    // ── 2. Write query FASTA ──
    let mut query_file = NamedTempFile::new().context("Cannot create temp query file")?;
    let mut query_to_cluster_idx: HashMap<String, usize> = HashMap::default();
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
            query_to_cluster_idx.insert(query_name.clone(), ci);
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
        let mut writer = BufWriter::with_capacity(256 * 1024, &mut target_file);
        let mut genomes_written = 0usize;

        let mut genome_indices_needed: Vec<usize> =
            searchable_genomes.iter().map(|(gi, _)| *gi).collect();
        genome_indices_needed.sort_unstable();
        genome_indices_needed.dedup();

        for &gi in &genome_indices_needed {
            let genome_name = &genome_names[gi];
            let Some(fasta_path) = genome_files.get(genome_name) else {
                continue;
            };

            let file = std::fs::File::open(fasta_path)
                .with_context(|| format!("Cannot open genome {}", fasta_path.display()))?;
            let reader = BufReader::with_capacity(128 * 1024, file);

            for line in reader.lines() {
                let line = line?;
                if let Some(header) = line.strip_prefix('>') {
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

    // ── 4. Pre-build minimap2 index (.mmi) ──
    let index_file = NamedTempFile::new().context("Cannot create temp index file")?;
    let index_path = index_file.path().to_string_lossy().to_string();

    let idx_output = Command::new("minimap2")
        .args([
            "-x",
            preset,
            "-t",
            &threads.to_string(),
            "-d",
            &index_path,
            &target_file.path().to_string_lossy(),
        ])
        .stderr(Stdio::piped())
        .stdout(Stdio::null())
        .output()
        .context("Failed to run minimap2 indexing — is it installed?")?;

    if !idx_output.status.success() {
        let stderr = String::from_utf8_lossy(&idx_output.stderr);
        anyhow::bail!("minimap2 indexing failed: {}", stderr);
    }
    info!("Built minimap2 index");

    // ── 5. Run minimap2 alignment with streamed stdout ──
    let mut child = Command::new("minimap2")
        .args([
            "-x",
            preset,
            "-t",
            &threads.to_string(),
            "--secondary=yes",
            "-N",
            "50",
            "-p",
            "0.5",
            &index_path,
            &query_file.path().to_string_lossy(),
        ])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .context("Failed to spawn minimap2 — is it installed?")?;

    let stdout = child.stdout.take().context("No stdout from minimap2")?;
    let paf_reader = BufReader::with_capacity(256 * 1024, stdout);

    // Capture stderr in a background thread so it doesn't block/deadlock
    let stderr = child.stderr.take();
    let stderr_handle = std::thread::spawn(move || -> String {
        let Some(stderr) = stderr else {
            return String::new();
        };
        let reader = BufReader::new(stderr);
        let mut buf = String::new();
        for line in reader.lines().map_while(Result::ok) {
            buf.push_str(&line);
            buf.push('\n');
        }
        buf
    });

    // ── 6. Build gap lookup ──
    let mut gap_set: ahash::HashSet<(usize, usize)> = ahash::HashSet::default();
    for (gi, cluster_indices) in searchable_genomes {
        for &ci in cluster_indices {
            gap_set.insert((*gi, ci));
        }
    }

    let genome_name_to_idx: HashMap<&str, usize> = genome_names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();

    // ── 7. Stream-parse PAF (no full stdout in RAM) ──
    let mut hits: HashMap<(usize, usize), Vec<Hit>> = HashMap::default();
    let mut total_hits = 0usize;
    let mut filtered_hits = 0usize;

    for line_result in paf_reader.lines() {
        let line = line_result?;
        let Some(fields) = paf_required_fields(&line) else {
            continue;
        };
        total_hits += 1;

        // Parse target: "genome_name|||original_contig rest_of_header"
        let target_field = fields[5];
        let Some(sep_pos) = target_field.find(CONTIG_SEP) else {
            continue;
        };
        let genome_name = &target_field[..sep_pos];
        let original_header = &target_field[sep_pos + CONTIG_SEP.len()..];
        let original_contig = original_header
            .split_whitespace()
            .next()
            .unwrap_or(original_header);

        let Some(&gi) = genome_name_to_idx.get(genome_name) else {
            continue;
        };

        let query_name = fields[0];
        let Some(&ci) = query_to_cluster_idx.get(query_name) else {
            continue;
        };
        let cluster_name = &cluster_names[ci];

        // Filter: only keep hits for actual gaps
        if !gap_set.contains(&(gi, ci)) {
            continue;
        }

        let hit = Hit {
            query: query_name.to_string(),
            cluster: cluster_name.clone(),
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
        hits.entry((gi, ci)).or_default().push(hit);
    }

    // Wait for minimap2 to finish and collect stderr
    let status = child.wait().context("minimap2 process failed")?;
    let stderr_output = stderr_handle.join().unwrap_or_default();

    // Log warnings from stderr regardless of exit code
    for line in stderr_output.lines() {
        if line.contains("[WARNING]") || line.contains("[W::") {
            log::warn!("minimap2: {}", line);
        }
    }

    if !status.success() {
        // Include first 5 stderr lines in error message for diagnostics
        let context: String = stderr_output.lines().take(5).collect::<Vec<_>>().join("\n");
        anyhow::bail!("minimap2 alignment exited with {}\n{}", status, context);
    }

    debug!(
        "minimap2 batch: {} total hits → {} passing filters → {} (genome, cluster) pairs",
        total_hits,
        filtered_hits,
        hits.len()
    );

    Ok(hits)
}

fn paf_required_fields(line: &str) -> Option<[&str; 12]> {
    let mut fields = line.split('\t');
    Some([
        fields.next()?,
        fields.next()?,
        fields.next()?,
        fields.next()?,
        fields.next()?,
        fields.next()?,
        fields.next()?,
        fields.next()?,
        fields.next()?,
        fields.next()?,
        fields.next()?,
        fields.next()?,
    ])
}

/// Parse a single PAF line into a Hit, resolving genome prefix and cluster name.
/// Returns None if the line is malformed, has no contig separator, or the genome is unknown.
/// Exposed for unit testing.
#[cfg(test)]
fn parse_paf_line(
    line: &str,
    genome_name_to_idx: &HashMap<&str, usize>,
    name_to_cluster: &HashMap<String, String>,
) -> Option<(usize, Hit)> {
    let fields = paf_required_fields(line)?;

    let target_field = fields[5];
    let sep_pos = target_field.find(CONTIG_SEP)?;
    let genome_name = &target_field[..sep_pos];
    let original_header = &target_field[sep_pos + CONTIG_SEP.len()..];
    let original_contig = original_header
        .split_whitespace()
        .next()
        .unwrap_or(original_header);

    let &gi = genome_name_to_idx.get(genome_name)?;

    let query_name = fields[0];
    let cluster_name = name_to_cluster
        .get(query_name)
        .map(|s| s.as_str())
        .unwrap_or(query_name);

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

    Some((gi, hit))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ahash::HashMapExt;

    // ── Hit method tests ──

    #[test]
    fn test_hit_identity() {
        let h = Hit {
            query: "q".into(),
            cluster: "c".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 1000,
            strand: '+',
            target: "t".into(),
            target_start: 0,
            target_end: 1000,
            matches: 950,
            alignment_length: 1000,
            mapq: 60,
        };
        assert!((h.identity() - 95.0).abs() < 0.01);
    }

    #[test]
    fn test_hit_identity_zero_alignment() {
        let h = Hit {
            query: "q".into(),
            cluster: "c".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 0,
            strand: '+',
            target: "t".into(),
            target_start: 0,
            target_end: 0,
            matches: 0,
            alignment_length: 0,
            mapq: 0,
        };
        assert_eq!(h.identity(), 0.0);
    }

    #[test]
    fn test_hit_query_coverage() {
        let h = Hit {
            query: "q".into(),
            cluster: "c".into(),
            query_length: 2000,
            query_start: 0,
            query_end: 1000,
            strand: '+',
            target: "t".into(),
            target_start: 0,
            target_end: 1000,
            matches: 1000,
            alignment_length: 1000,
            mapq: 60,
        };
        assert!((h.query_coverage() - 50.0).abs() < 0.01);
    }

    #[test]
    fn test_hit_query_coverage_zero_length() {
        let h = Hit {
            query: "q".into(),
            cluster: "c".into(),
            query_length: 0,
            query_start: 0,
            query_end: 0,
            strand: '+',
            target: "t".into(),
            target_start: 0,
            target_end: 0,
            matches: 0,
            alignment_length: 0,
            mapq: 0,
        };
        assert_eq!(h.query_coverage(), 0.0);
    }

    #[test]
    fn test_hit_score() {
        let h = Hit {
            query: "q".into(),
            cluster: "c".into(),
            query_length: 1000,
            query_start: 0,
            query_end: 1000,
            strand: '+',
            target: "t".into(),
            target_start: 0,
            target_end: 1000,
            matches: 1000,
            alignment_length: 1000,
            mapq: 60,
        };
        assert!((h.score() - 1.0).abs() < 0.001); // 100% × 100% = 1.0
    }

    // ── PAF parsing tests ──

    fn make_genome_idx() -> HashMap<&'static str, usize> {
        let mut m = HashMap::new();
        m.insert("genome1", 0);
        m.insert("genome2", 1);
        m
    }

    fn make_name_to_cluster() -> HashMap<String, String> {
        let mut m = HashMap::new();
        m.insert("clA~~~allele1".into(), "clA".into());
        m
    }

    #[test]
    fn test_parse_paf_line_valid() {
        let g2i = make_genome_idx();
        let n2c = make_name_to_cluster();
        // Standard PAF: query, qlen, qstart, qend, strand, target, tlen, tstart, tend, matches, alnlen, mapq
        let line = "clA\t1000\t0\t1000\t+\tgenome1|||contig_1\t5000000\t100\t1100\t990\t1000\t60";
        let result = parse_paf_line(line, &g2i, &n2c);
        assert!(result.is_some());
        let (gi, hit) = result.unwrap();
        assert_eq!(gi, 0);
        assert_eq!(hit.target, "contig_1");
        assert_eq!(hit.cluster, "clA");
        assert_eq!(hit.target_start, 100);
        assert_eq!(hit.target_end, 1100);
    }

    #[test]
    fn test_parse_paf_line_with_allele_name() {
        let g2i = make_genome_idx();
        let n2c = make_name_to_cluster();
        let line = "clA~~~allele1\t800\t0\t800\t-\tgenome2|||chr\t4000000\t500\t1300\t780\t800\t30";
        let result = parse_paf_line(line, &g2i, &n2c);
        let (gi, hit) = result.unwrap();
        assert_eq!(gi, 1);
        assert_eq!(hit.cluster, "clA"); // resolved from allele to cluster
        assert_eq!(hit.query, "clA~~~allele1");
        assert_eq!(hit.strand, '-');
    }

    #[test]
    fn test_parse_paf_line_no_separator() {
        let g2i = make_genome_idx();
        let n2c = make_name_to_cluster();
        // Target without ||| separator → should return None
        let line = "clA\t1000\t0\t1000\t+\tplain_contig\t5000000\t100\t1100\t990\t1000\t60";
        assert!(parse_paf_line(line, &g2i, &n2c).is_none());
    }

    #[test]
    fn test_parse_paf_line_unknown_genome() {
        let g2i = make_genome_idx();
        let n2c = make_name_to_cluster();
        let line =
            "clA\t1000\t0\t1000\t+\tunknown_genome|||contig_1\t5000000\t100\t1100\t990\t1000\t60";
        assert!(parse_paf_line(line, &g2i, &n2c).is_none());
    }

    #[test]
    fn test_parse_paf_line_too_few_fields() {
        let g2i = make_genome_idx();
        let n2c = make_name_to_cluster();
        let line = "clA\t1000\t0\t1000\t+\tgenome1|||contig\t5000000";
        assert!(parse_paf_line(line, &g2i, &n2c).is_none());
    }

    #[test]
    fn test_parse_paf_line_contig_with_spaces() {
        let g2i = make_genome_idx();
        let n2c = make_name_to_cluster();
        // Target header with extra description after space → only first token as contig name
        let line =
            "clA\t1000\t0\t1000\t+\tgenome1|||ctg1 length=5000\t5000000\t100\t1100\t990\t1000\t60";
        let (_, hit) = parse_paf_line(line, &g2i, &n2c).unwrap();
        assert_eq!(hit.target, "ctg1");
    }

    // ── CONTIG_SEP constant test ──

    #[test]
    fn test_contig_sep_is_unlikely_in_real_names() {
        assert_eq!(CONTIG_SEP, "|||");
        assert!(!"normal_contig_name".contains(CONTIG_SEP));
        assert!(!"NC_000962.3".contains(CONTIG_SEP));
    }
}
