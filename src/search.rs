use ahash::HashMap;
use anyhow::{Context, Result};
use log::{debug, info};
use std::io::{BufRead, BufWriter, Write};
use std::process::{Command, Stdio};
use tempfile::NamedTempFile;

use crate::parser::open_maybe_gzip;
use crate::types::{ClusterIdx, GenomeIdx};

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

/// Minimum peak alignment score once base-level alignment (`-c`) is on. The asm
/// presets set `-s200`, and with a match score of 1 that silently drops every
/// gene shorter than about 200 bp: small proteins were unrescuable in every
/// exact pass.
const MIN_PEAK_SCORE: &str = "40";

/// Bonus for an alignment that reaches either end of the query. It has to beat
/// the asm5 mismatch penalty of 19, otherwise one SNP in the last twenty bases
/// of a gene clips the alignment short of the gene's end, and the clipped span
/// then reads as a truncated gene or a different start codon.
const END_BONUS: &str = "40";

#[allow(clippy::too_many_arguments)]
/// Search every gap genome for its own absent clusters, one minimap2 run per
/// genome, with base-level alignment so identities are exact.
///
/// One genome at a time, not all genomes concatenated into a single target.
/// The single target let genomes that carry a cluster compete with the genome
/// where it is missing, under one shared `-N` secondary budget and one `-p`
/// ratio. Once more than about fifty genomes carried a gene its full-length
/// carrier copies used up the budget and the gap genome's hit was dropped, so
/// on a few hundred genomes almost nothing was rescued, with exit 0; and a gap
/// genome's weaker copy, or the pieces of a gene split across contigs, fell
/// under half the carriers' chain score and were cut by `-p 0.5`. Searched
/// alone, a genome's hits only compete with its own other copies.
///
/// The target each run sees is written by recore's own FASTA reader rather
/// than handed to minimap2 as the user's file, so the aligner and the rest of
/// the pipeline agree on what the assembly contains, and a malformed FASTA
/// fails here with its path instead of producing hits on the wrong contig.
///
/// `max_secondary` bounds the secondary alignments per query within one
/// genome. Returns hits keyed by (genome, cluster), only for the requested gaps.
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
    // One job per genome. A genome listed twice has its gap lists merged, and
    // the BTreeMap keeps the job order, and so the log, deterministic.
    let mut gaps_by_genome: std::collections::BTreeMap<GenomeIdx, Vec<ClusterIdx>> =
        std::collections::BTreeMap::new();
    for (gi, cluster_indices) in searchable_genomes {
        gaps_by_genome
            .entry(*gi)
            .or_default()
            .extend(cluster_indices.iter().copied());
    }
    let jobs: Vec<(GenomeIdx, &std::path::Path, Vec<ClusterIdx>)> = gaps_by_genome
        .into_iter()
        .filter_map(|(gi, mut clusters)| {
            clusters.sort_unstable();
            clusters.dedup();
            let path = genome_files.get(&genome_names[gi.val()])?;
            Some((gi, path.as_path(), clusters))
        })
        .collect();
    if jobs.is_empty() {
        return Ok(HashMap::default());
    }

    let total_queries: usize = jobs
        .iter()
        .flat_map(|(_, _, clusters)| clusters)
        .filter_map(|ci| queries.get(&cluster_names[ci.val()]))
        .map(Vec::len)
        .sum();
    info!(
        "Prepared search: {} query sequences across {} genomes, one minimap2 run per genome",
        total_queries,
        jobs.len()
    );

    let max_secondary = max_secondary.max(1).to_string();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .context("Failed to build the search thread pool")?;
    let results: Vec<Result<Vec<(ClusterIdx, Hit)>>> = pool.install(|| {
        use rayon::prelude::*;
        jobs.par_iter()
            .map(|(_, path, clusters)| {
                search_one_genome(
                    queries,
                    path,
                    clusters,
                    cluster_names,
                    min_identity,
                    preset,
                    exact_base_alignment,
                    &max_secondary,
                )
            })
            .collect()
    });

    let mut hits: HashMap<(GenomeIdx, ClusterIdx), Vec<Hit>> = HashMap::default();
    let mut kept = 0usize;
    for ((gi, path, _), result) in jobs.iter().zip(results) {
        let found = result.with_context(|| {
            format!(
                "Search failed for genome {} ({})",
                genome_names[gi.val()],
                path.display()
            )
        })?;
        kept += found.len();
        for (ci, hit) in found {
            hits.entry((*gi, ci)).or_default().push(hit);
        }
    }
    debug!(
        "minimap2: {} hits kept across {} (genome, cluster) pairs",
        kept,
        hits.len()
    );
    Ok(hits)
}

/// Align one genome's absent clusters against that genome alone.
#[allow(clippy::too_many_arguments)]
fn search_one_genome(
    queries: &HashMap<String, Vec<Vec<u8>>>,
    fasta_path: &std::path::Path,
    clusters: &[ClusterIdx],
    cluster_names: &[String],
    min_identity: f64,
    preset: &str,
    exact_base_alignment: bool,
    max_secondary: &str,
) -> Result<Vec<(ClusterIdx, Hit)>> {
    let mut query_file = NamedTempFile::new().context("Cannot create temp query file")?;
    let mut query_to_cluster_idx: HashMap<String, ClusterIdx> = HashMap::default();
    {
        let mut writer = BufWriter::new(&mut query_file);
        for &ci in clusters {
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
                writeln!(writer, ">{}", query_name)?;
                writer.write_all(seq)?;
                writeln!(writer)?;
                query_to_cluster_idx.insert(query_name, ci);
            }
        }
        writer.flush()?;
    }
    if query_to_cluster_idx.is_empty() {
        return Ok(Vec::new());
    }

    let mut target_file = NamedTempFile::new().context("Cannot create temp target file")?;
    {
        let mut writer = BufWriter::with_capacity(256 * 1024, &mut target_file);
        let mut reader = open_maybe_gzip(fasta_path)?;
        write_search_target(&mut reader, &mut writer)
            .with_context(|| format!("Cannot read genome {}", fasta_path.display()))?;
        writer
            .flush()
            .context("Cannot write the temporary search target")?;
    }

    let mut args: Vec<&str> = vec!["-x", preset, "-t", "1"];
    if exact_base_alignment {
        args.extend(["-c", "-s", MIN_PEAK_SCORE, "--end-bonus", END_BONUS]);
    }
    let target_path = target_file.path().to_string_lossy().to_string();
    let query_path = query_file.path().to_string_lossy().to_string();
    args.extend([
        "--secondary=yes",
        "-N",
        max_secondary,
        "-p",
        "0.5",
        &target_path,
        &query_path,
    ]);

    let output = Command::new("minimap2")
        .args(&args)
        .stdin(Stdio::null())
        .output()
        .context("Failed to run minimap2: is it installed and on PATH?")?;
    let stderr = String::from_utf8_lossy(&output.stderr);
    for line in stderr.lines() {
        if line.contains("[WARNING]") || line.contains("[W::") {
            log::warn!("minimap2 ({}): {}", fasta_path.display(), line);
        }
    }
    if !output.status.success() {
        let context: Vec<&str> = stderr.lines().take(5).collect();
        anyhow::bail!(
            "minimap2 exited with {}\n{}",
            output.status,
            context.join("\n")
        );
    }

    // Target names are the contig tokens written above, so they are UTF-8
    // and contain no whitespace. Query names are recore's own.
    let stdout =
        std::str::from_utf8(&output.stdout).context("minimap2 printed output that is not UTF-8")?;
    let mut hits = Vec::new();
    for line in stdout.lines() {
        let Some(fields) = paf_required_fields(line) else {
            continue;
        };
        let Some(&ci) = query_to_cluster_idx.get(fields[0]) else {
            continue;
        };
        let Some(hit) =
            hit_from_paf_fields(&fields, fields[0], &cluster_names[ci.val()], fields[5])
        else {
            continue;
        };
        if hit.identity() < min_identity {
            continue;
        }
        hits.push((ci, hit));
    }
    Ok(hits)
}

/// Copy a genome FASTA into the plain, normalised form minimap2 is given.
///
/// Each header is cut to its first whitespace-delimited token, which is the
/// contig name everything else in recore matches on, and line endings are
/// normalised. It refuses, naming the line, the inputs that used to produce
/// confident wrong answers: sequence before the first header (it was spliced
/// onto the previous genome's last contig), an empty or non-UTF-8 contig name,
/// and a contig name used twice in one assembly, where every later lookup by
/// name would silently pick one of the two.
fn write_search_target<R: BufRead + ?Sized, W: Write>(
    reader: &mut R,
    writer: &mut W,
) -> Result<()> {
    let mut line = Vec::with_capacity(8192);
    let mut seen: ahash::HashSet<String> = ahash::HashSet::default();
    let mut in_record = false;
    let mut line_no = 0usize;
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        line_no += 1;
        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        if let Some(header) = line.strip_prefix(b">") {
            let token = header
                .split(|b| b.is_ascii_whitespace())
                .next()
                .unwrap_or_default();
            if token.is_empty() {
                anyhow::bail!("line {}: FASTA header with no contig name", line_no);
            }
            let name = std::str::from_utf8(token)
                .map_err(|_| anyhow::anyhow!("line {}: contig name is not valid UTF-8", line_no))?;
            if !seen.insert(name.to_string()) {
                anyhow::bail!(
                    "line {}: contig name '{}' appears twice in this assembly",
                    line_no,
                    name
                );
            }
            writer.write_all(b">")?;
            writer.write_all(token)?;
            writer.write_all(b"\n")?;
            in_record = true;
        } else if line.iter().all(u8::is_ascii_whitespace) {
            continue;
        } else if !in_record {
            anyhow::bail!(
                "line {}: sequence before the first '>' header; this is not a FASTA file \
                 or its first header is missing",
                line_no
            );
        } else {
            writer.write_all(&line)?;
            writer.write_all(b"\n")?;
        }
    }
    if !in_record {
        anyhow::bail!("no FASTA records: the file is empty or has no '>' header");
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

/// Parse a single PAF line into a Hit, resolving allele query names to their
/// cluster. Each search targets one genome, so the target field is the contig.
/// Exposed for unit testing.
#[cfg(test)]
fn parse_paf_line(line: &str, name_to_cluster: &HashMap<String, String>) -> Option<Hit> {
    let fields = paf_required_fields(line)?;
    let query_name = fields[0];
    let cluster_name = name_to_cluster
        .get(query_name)
        .map(|s| s.as_str())
        .unwrap_or(query_name);
    hit_from_paf_fields(&fields, query_name, cluster_name, fields[5])
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

    fn make_name_to_cluster() -> HashMap<String, String> {
        let mut m = HashMap::new();
        m.insert("clA~~~allele1".into(), "clA".into());
        m
    }

    #[test]
    fn test_parse_paf_line_valid() {
        let n2c = make_name_to_cluster();
        // Standard PAF: query, qlen, qstart, qend, strand, target, tlen, tstart, tend, matches, alnlen, mapq
        let line = "clA\t1000\t0\t1000\t+\tcontig_1\t5000000\t100\t1100\t990\t1000\t60";
        let hit = parse_paf_line(line, &n2c).unwrap();
        assert_eq!(hit.target, "contig_1");
        assert_eq!(hit.cluster, "clA");
        assert_eq!(hit.target_start, 100);
        assert_eq!(hit.target_end, 1100);
    }

    #[test]
    fn test_parse_paf_line_with_allele_name() {
        let n2c = make_name_to_cluster();
        let line = "clA~~~allele1\t800\t0\t800\t-\tchr\t4000000\t500\t1300\t780\t800\t30";
        let hit = parse_paf_line(line, &n2c).unwrap();
        assert_eq!(hit.cluster, "clA"); // resolved from allele to cluster
        assert_eq!(hit.query, "clA~~~allele1");
        assert_eq!(hit.strand, '-');
    }

    #[test]
    fn test_parse_paf_line_too_few_fields() {
        let n2c = make_name_to_cluster();
        let line = "clA\t1000\t0\t1000\t+\tcontig\t5000000";
        assert!(parse_paf_line(line, &n2c).is_none());
    }

    #[test]
    fn test_parse_paf_line_rejects_invalid_numeric_fields() {
        let n2c = make_name_to_cluster();
        let line = "clA\t1000\tbad\t1000\t+\tcontig\t5000000\t100\t1100\t990\t1000\t60";
        assert!(parse_paf_line(line, &n2c).is_none());
    }

    #[test]
    fn test_parse_paf_line_rejects_invalid_coordinate_order() {
        let n2c = make_name_to_cluster();
        let line = "clA\t1000\t900\t100\t+\tcontig\t5000000\t100\t1100\t990\t1000\t60";
        assert!(parse_paf_line(line, &n2c).is_none());
    }

    // ── Search target writer ──

    fn target_of(input: &[u8]) -> Result<String> {
        let mut reader = std::io::Cursor::new(input);
        let mut output = Vec::new();
        write_search_target(&mut reader, &mut output)?;
        Ok(String::from_utf8(output).unwrap())
    }

    #[test]
    fn test_search_target_keeps_the_contig_token_only() {
        let text = target_of(b">ctg1 description words\nACGT\nAC\n>ctg2\nTTAA").unwrap();
        assert_eq!(text, ">ctg1\nACGT\nAC\n>ctg2\nTTAA\n");
    }

    #[test]
    fn test_search_target_normalises_crlf_and_blank_lines() {
        let text = target_of(b">ctg1 description\r\nACGT\r\n\r\n>ctg2\r\nTT\r\n").unwrap();
        assert_eq!(text, ">ctg1\nACGT\n>ctg2\nTT\n");
    }

    /// Sequence before the first header used to be spliced onto the previous
    /// genome's last contig in the shared target, so it produced hits there.
    #[test]
    fn test_search_target_refuses_sequence_before_a_header() {
        let err = target_of(b"ACGTACGT\n>ctg1\nACGT\n")
            .unwrap_err()
            .to_string();
        assert!(err.contains("before the first '>' header"), "{err}");
    }

    #[test]
    fn test_search_target_refuses_an_empty_file() {
        let err = target_of(b"").unwrap_err().to_string();
        assert!(err.contains("no FASTA records"), "{err}");
    }

    /// Two contigs with one name make every later lookup by name ambiguous.
    #[test]
    fn test_search_target_refuses_duplicate_contig_names() {
        let err = target_of(b">ctg1\nACGT\n>ctg1 again\nTTTT\n")
            .unwrap_err()
            .to_string();
        assert!(err.contains("'ctg1' appears twice"), "{err}");
    }

    /// Only the contig token has to be UTF-8: a Latin-1 description is fine.
    #[test]
    fn test_search_target_tolerates_non_utf8_descriptions_only() {
        let text = target_of(b">ctg1 caf\xe9\nACGT\n").unwrap();
        assert_eq!(text, ">ctg1\nACGT\n");
        let err = target_of(b">ctg\xe91\nACGT\n").unwrap_err().to_string();
        assert!(err.contains("not valid UTF-8"), "{err}");
    }
}
