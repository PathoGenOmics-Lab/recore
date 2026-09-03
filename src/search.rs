use ahash::HashMap;
use anyhow::{Context, Result};
use log::{debug, info};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::{Command, Stdio};
use tempfile::NamedTempFile;

use crate::parser::open_maybe_gzip;
use crate::types::{ClusterIdx, GenomeIdx};

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
    searchable_genomes: &[(GenomeIdx, Vec<ClusterIdx>)],
    genome_names: &[String],
    cluster_names: &[String],
    threads: usize,
    min_identity: f64,
    preset: &str,
    exact_base_alignment: bool,
    max_secondary: usize,
) -> Result<HashMap<(GenomeIdx, ClusterIdx), Vec<Hit>>> {
    if searchable_genomes.is_empty() {
        return Ok(HashMap::default());
    }

    // ── 1. Collect unique query clusters needed across all genomes ──
    let mut needed_cluster_bits = vec![false; cluster_names.len()];
    for (_gi, cluster_indices) in searchable_genomes {
        for &ci in cluster_indices {
            needed_cluster_bits[ci.val()] = true;
        }
    }
    let needed_clusters: Vec<ClusterIdx> = needed_cluster_bits
        .iter()
        .enumerate()
        .filter_map(|(ci, needed)| needed.then_some(ClusterIdx(ci)))
        .collect();

    // ── 2. Write query FASTA ──
    let mut query_file = NamedTempFile::new().context("Cannot create temp query file")?;
    let mut query_to_cluster_idx: HashMap<String, ClusterIdx> = HashMap::default();
    let mut total_queries = 0usize;

    for &ci in &needed_clusters {
        let cluster_name = &cluster_names[ci.val()];
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

    let mut genome_indices_needed: Vec<GenomeIdx> =
        searchable_genomes.iter().map(|(gi, _)| *gi).collect();
    genome_indices_needed.sort_unstable();
    genome_indices_needed.dedup();
    let max_secondary = max_secondary.max(1).to_string();

    // ── 3. Concatenate genome FASTAs with prefixed contig names ──
    let mut target_file = NamedTempFile::new().context("Cannot create temp target file")?;
    {
        let mut writer = BufWriter::with_capacity(256 * 1024, &mut target_file);
        let mut genomes_written = 0usize;

        for &gi in &genome_indices_needed {
            let genome_name = &genome_names[gi.val()];
            let Some(fasta_path) = genome_files.get(genome_name) else {
                continue;
            };

            let mut reader = open_maybe_gzip(fasta_path)
                .with_context(|| format!("Cannot open genome {}", fasta_path.display()))?;
            // Context on the read as well as on the open. A compressed file
            // fails while being decoded, long after it opened cleanly, and
            // "unexpected end of file" with no path is close to useless to
            // somebody holding five thousand genomes.
            write_prefixed_fasta(&mut reader, &mut writer, genome_name)
                .with_context(|| format!("Cannot read genome {}", fasta_path.display()))?;
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
    let mut aln_args = vec![
        "-x".to_string(),
        preset.to_string(),
        "-t".to_string(),
        threads.to_string(),
    ];
    if exact_base_alignment {
        aln_args.push("-c".to_string());
    }
    aln_args.extend([
        "--secondary=yes".to_string(),
        "-N".to_string(),
        max_secondary,
        "-p".to_string(),
        "0.5".to_string(),
        index_path,
        query_file.path().to_string_lossy().to_string(),
    ]);

    let mut child = Command::new("minimap2")
        .args(&aln_args)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .context("Failed to spawn minimap2 — is it installed?")?;

    let stdout = child.stdout.take().context("No stdout from minimap2")?;
    let mut paf_reader = BufReader::with_capacity(256 * 1024, stdout);

    // Capture stderr in a background thread so it doesn't block/deadlock
    let stderr = child.stderr.take();
    let stderr_handle = std::thread::spawn(move || -> String {
        let Some(stderr) = stderr else {
            return String::new();
        };
        let mut reader = BufReader::new(stderr);
        let mut line = String::with_capacity(512);
        let mut context = String::new();
        let mut n_context = 0usize;
        loop {
            line.clear();
            let Ok(n) = reader.read_line(&mut line) else {
                break;
            };
            if n == 0 {
                break;
            }
            let trimmed = line.trim_end_matches(['\r', '\n']);
            if trimmed.contains("[WARNING]") || trimmed.contains("[W::") {
                log::warn!("minimap2: {}", trimmed);
            }
            if n_context < 5 {
                context.push_str(trimmed);
                context.push('\n');
                n_context += 1;
            }
        }
        context
    });

    // ── 6. Build gap lookup ──
    let mut gap_lookup = vec![vec![false; cluster_names.len()]; genome_names.len()];
    for (gi, cluster_indices) in searchable_genomes {
        for &ci in cluster_indices {
            gap_lookup[gi.val()][ci.val()] = true;
        }
    }

    let genome_name_to_idx: HashMap<&str, GenomeIdx> = genome_names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), GenomeIdx(i)))
        .collect();

    // ── 7. Stream-parse PAF (no full stdout in RAM) ──
    let mut hits: HashMap<(GenomeIdx, ClusterIdx), Vec<Hit>> = HashMap::default();
    let mut total_hits = 0usize;
    let mut filtered_hits = 0usize;
    let mut line = String::with_capacity(1024);

    loop {
        line.clear();
        let n = paf_reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        let line_trimmed = line.trim_end_matches(['\r', '\n']);
        let Some(fields) = paf_required_fields(line_trimmed) else {
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
        let cluster_name = &cluster_names[ci.val()];

        // Filter: only keep hits for actual gaps
        if !gap_lookup[gi.val()][ci.val()] {
            continue;
        }

        let Some(hit) = hit_from_paf_fields(&fields, query_name, cluster_name, original_contig)
        else {
            continue;
        };

        if hit.identity() < min_identity {
            continue;
        }

        filtered_hits += 1;
        hits.entry((gi, ci)).or_default().push(hit);
    }

    // Wait for minimap2 to finish and collect stderr
    let status = child.wait().context("minimap2 process failed")?;
    let stderr_context = stderr_handle.join().unwrap_or_default();

    if !status.success() {
        anyhow::bail!(
            "minimap2 alignment exited with {}\n{}",
            status,
            stderr_context
        );
    }

    debug!(
        "minimap2 batch: {} total hits → {} passing filters → {} (genome, cluster) pairs",
        total_hits,
        filtered_hits,
        hits.len()
    );

    Ok(hits)
}

fn write_prefixed_fasta<R: BufRead + ?Sized, W: Write>(
    reader: &mut R,
    writer: &mut W,
    genome_name: &str,
) -> Result<()> {
    let mut line = Vec::with_capacity(8192);
    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }

        if line.starts_with(b">") {
            let mut header_end = line.len();
            while header_end > 1 && matches!(line[header_end - 1], b'\n' | b'\r') {
                header_end -= 1;
            }
            writer.write_all(b">")?;
            writer.write_all(genome_name.as_bytes())?;
            writer.write_all(CONTIG_SEP.as_bytes())?;
            writer.write_all(&line[1..header_end])?;
            writer.write_all(b"\n")?;
        } else {
            writer.write_all(&line)?;
            if !line.ends_with(b"\n") {
                writer.write_all(b"\n")?;
            }
        }
    }
    Ok(())
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

fn parse_u64_field(value: &str) -> Option<u64> {
    value.parse().ok()
}

fn parse_u32_field(value: &str) -> Option<u32> {
    value.parse().ok()
}

fn hit_from_paf_fields(
    fields: &[&str; 12],
    query_name: &str,
    cluster_name: &str,
    original_contig: &str,
) -> Option<Hit> {
    let query_length = parse_u64_field(fields[1])?;
    let query_start = parse_u64_field(fields[2])?;
    let query_end = parse_u64_field(fields[3])?;
    let strand = match fields[4] {
        "+" => '+',
        "-" => '-',
        _ => return None,
    };
    let target_length = parse_u64_field(fields[6])?;
    let target_start = parse_u64_field(fields[7])?;
    let target_end = parse_u64_field(fields[8])?;
    let matches = parse_u64_field(fields[9])?;
    let alignment_length = parse_u64_field(fields[10])?;
    let mapq = parse_u32_field(fields[11])?;

    if query_start > query_end
        || query_end > query_length
        || target_start > target_end
        || target_end > target_length
        || matches > alignment_length
    {
        return None;
    }

    Some(Hit {
        query: query_name.to_string(),
        cluster: cluster_name.to_string(),
        query_length,
        query_start,
        query_end,
        strand,
        target: original_contig.to_string(),
        target_start,
        target_end,
        matches,
        alignment_length,
        mapq,
    })
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

    let hit = hit_from_paf_fields(&fields, query_name, cluster_name, original_contig)?;

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
    fn test_parse_paf_line_rejects_invalid_numeric_fields() {
        let g2i = make_genome_idx();
        let n2c = make_name_to_cluster();
        let line = "clA\t1000\tbad\t1000\t+\tgenome1|||contig\t5000000\t100\t1100\t990\t1000\t60";
        assert!(parse_paf_line(line, &g2i, &n2c).is_none());
    }

    #[test]
    fn test_parse_paf_line_rejects_invalid_coordinate_order() {
        let g2i = make_genome_idx();
        let n2c = make_name_to_cluster();
        let line = "clA\t1000\t900\t100\t+\tgenome1|||contig\t5000000\t100\t1100\t990\t1000\t60";
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

    #[test]
    fn test_write_prefixed_fasta_reuses_input_lines() {
        let input = b">ctg1 description\nACGT\n>ctg2\nTTAA";
        let mut reader = std::io::Cursor::new(input);
        let mut output = Vec::new();

        write_prefixed_fasta(&mut reader, &mut output, "genome1").unwrap();
        let text = String::from_utf8(output).unwrap();

        assert_eq!(
            text,
            ">genome1|||ctg1 description\nACGT\n>genome1|||ctg2\nTTAA\n"
        );
    }

    #[test]
    fn test_write_prefixed_fasta_handles_crlf_headers() {
        let input = b">ctg1 description\r\nACGT\r\n";
        let mut reader = std::io::Cursor::new(input);
        let mut output = Vec::new();

        write_prefixed_fasta(&mut reader, &mut output, "genome1").unwrap();
        let text = String::from_utf8(output).unwrap();

        assert_eq!(text, ">genome1|||ctg1 description\nACGT\r\n");
    }

    // ── CONTIG_SEP constant test ──

    #[test]
    fn test_contig_sep_is_unlikely_in_real_names() {
        assert_eq!(CONTIG_SEP, "|||");
        assert!(!"normal_contig_name".contains(CONTIG_SEP));
        assert!(!"NC_000962.3".contains(CONTIG_SEP));
    }
}
