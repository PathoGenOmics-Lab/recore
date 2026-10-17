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
    /// Length of the target contig, from PAF column 7; 0 when unknown. Needed
    /// to see that a piece of a gene reaches the end of its contig.
    pub target_length: u64,
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
/// genome. With `cache_dir`, results are kept and reused between runs (see
/// [`SearchCache`]). Returns hits keyed by (genome, cluster), only for the
/// requested gaps.
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
    cache_dir: Option<&std::path::Path>,
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

    let settings = SearchSettings {
        min_identity,
        preset,
        exact_base_alignment,
        max_secondary: max_secondary.max(1).to_string(),
    };
    let cache = cache_dir
        .map(|dir| SearchCache::open(dir, &settings))
        .transpose()?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .context("Failed to build the search thread pool")?;
    let results: Vec<Result<GenomeSearch>> = pool.install(|| {
        use rayon::prelude::*;
        jobs.par_iter()
            .map(|(_, path, clusters)| {
                search_one_genome(
                    queries,
                    path,
                    clusters,
                    cluster_names,
                    &settings,
                    cache.as_ref(),
                )
            })
            .collect()
    });

    let mut hits: HashMap<(GenomeIdx, ClusterIdx), Vec<Hit>> = HashMap::default();
    let mut kept = 0usize;
    let mut cache_use = CacheUse::default();
    let mut genomes_without_minimap2 = 0usize;
    for ((gi, path, _), result) in jobs.iter().zip(results) {
        let (found, used) = result.with_context(|| {
            format!(
                "Search failed for genome {} ({})",
                genome_names[gi.val()],
                path.display()
            )
        })?;
        if used.searched == 0 && used.reused > 0 {
            genomes_without_minimap2 += 1;
        }
        cache_use.reused += used.reused;
        cache_use.searched += used.searched;
        kept += found.len();
        for (ci, hit) in found {
            hits.entry((*gi, ci)).or_default().push(hit);
        }
    }
    if let Some(cache) = cache.as_ref() {
        info!(
            "Search cache {}: {} of {} query sequences reused, {} searched; {} of {} genomes needed no minimap2 run",
            cache.dir.display(),
            cache_use.reused,
            cache_use.reused + cache_use.searched,
            cache_use.searched,
            genomes_without_minimap2,
            jobs.len()
        );
    }
    debug!(
        "minimap2: {} hits kept across {} (genome, cluster) pairs",
        kept,
        hits.len()
    );
    Ok(hits)
}

/// The settings that reach minimap2 or filter what it returns.
struct SearchSettings<'a> {
    min_identity: f64,
    preset: &'a str,
    exact_base_alignment: bool,
    max_secondary: String,
}

impl SearchSettings<'_> {
    fn minimap2_args(&self) -> Vec<&str> {
        let mut args: Vec<&str> = vec!["-x", self.preset, "-t", "1"];
        if self.exact_base_alignment {
            args.extend(["-c", "-s", MIN_PEAK_SCORE, "--end-bonus", END_BONUS]);
        }
        args.extend(["--secondary=yes", "-N", &self.max_secondary, "-p", "0.5"]);
        args
    }
}

/// One sequence to align against a genome: a cluster's representative, or one
/// of its alleles.
struct QueryItem<'a> {
    name: String,
    ci: ClusterIdx,
    seq: &'a [u8],
}

/// A query's hits from its alignment, and from realigning what that alignment
/// left uncovered, each in the order minimap2 reported them.
#[derive(Debug, Clone, Default)]
struct ItemHits {
    first: Vec<Hit>,
    remainder: Vec<Hit>,
}

impl ItemHits {
    /// The same hits under another query's name, for a sequence shared by two
    /// queries or read back from the cache.
    fn named(mut self, query: &str, cluster: &str) -> Self {
        for hit in self.first.iter_mut().chain(self.remainder.iter_mut()) {
            hit.query = query.to_string();
            hit.cluster = cluster.to_string();
        }
        self
    }
}

/// One genome's hits, by cluster, and how much of its search the cache spared.
type GenomeSearch = (Vec<(ClusterIdx, Hit)>, CacheUse);

#[derive(Debug, Clone, Copy, Default)]
struct CacheUse {
    reused: usize,
    searched: usize,
}

/// Align one genome's absent clusters against that genome alone.
fn search_one_genome(
    queries: &HashMap<String, Vec<Vec<u8>>>,
    fasta_path: &std::path::Path,
    clusters: &[ClusterIdx],
    cluster_names: &[String],
    settings: &SearchSettings<'_>,
    cache: Option<&SearchCache>,
) -> Result<GenomeSearch> {
    let mut items: Vec<QueryItem> = Vec::new();
    for &ci in clusters {
        let cluster_name = &cluster_names[ci.val()];
        let Some(alleles) = queries.get(cluster_name) else {
            continue;
        };
        for (ai, seq) in alleles.iter().enumerate() {
            let name = if alleles.len() == 1 {
                cluster_name.clone()
            } else {
                format!("{}~~~ALLELE~{}", cluster_name, ai)
            };
            items.push(QueryItem {
                name,
                ci,
                seq: seq.as_slice(),
            });
        }
    }
    if items.is_empty() {
        return Ok((Vec::new(), CacheUse::default()));
    }

    let Some(cache) = cache else {
        let all: Vec<&QueryItem> = items.iter().collect();
        let per_item = search_items(&all, fasta_path, cluster_names, settings)?;
        return Ok((in_search_order(&items, &per_item), CacheUse::default()));
    };

    let path = cache.file_for(fasta_path)?;
    let mut stored = cache.load(&path);
    let hashes: Vec<u128> = items.iter().map(|item| stable_hash(&[item.seq])).collect();
    // Each sequence not in the cache is searched once, however many queries
    // share it: identical sequences align identically.
    let mut pending: Vec<usize> = Vec::new();
    let mut pending_hashes: ahash::HashSet<u128> = ahash::HashSet::default();
    for (k, hash) in hashes.iter().enumerate() {
        if !stored.contains_key(hash) && pending_hashes.insert(*hash) {
            pending.push(k);
        }
    }
    if !pending.is_empty() {
        let to_search: Vec<&QueryItem> = pending.iter().map(|&k| &items[k]).collect();
        let found = search_items(&to_search, fasta_path, cluster_names, settings)?;
        for (&k, hits) in pending.iter().zip(found) {
            stored.insert(hashes[k], hits);
        }
        cache.store(&path, &stored)?;
    }
    let per_item: Vec<ItemHits> = items
        .iter()
        .zip(&hashes)
        .map(|(item, hash)| {
            stored
                .get(hash)
                .cloned()
                .unwrap_or_default()
                .named(&item.name, &cluster_names[item.ci.val()])
        })
        .collect();
    let used = CacheUse {
        reused: items.len() - pending.len(),
        searched: pending.len(),
    };
    Ok((in_search_order(&items, &per_item), used))
}

/// Hits in the order a search of every query at once reports them: first
/// alignments in query order, then the realigned remainders by query name. A
/// cached run hands its hits on in this order too, so it decides exactly what
/// an uncached run decides.
fn in_search_order(items: &[QueryItem], per_item: &[ItemHits]) -> Vec<(ClusterIdx, Hit)> {
    let mut out = Vec::new();
    for (item, hits) in items.iter().zip(per_item) {
        out.extend(hits.first.iter().cloned().map(|hit| (item.ci, hit)));
    }
    let mut by_name: Vec<usize> = (0..items.len()).collect();
    by_name.sort_by(|&a, &b| items[a].name.cmp(&items[b].name));
    for k in by_name {
        out.extend(
            per_item[k]
                .remainder
                .iter()
                .cloned()
                .map(|hit| (items[k].ci, hit)),
        );
    }
    out
}

/// Run minimap2 for these queries against one genome. Returns each query's
/// hits, aligned with `items`.
fn search_items(
    items: &[&QueryItem],
    fasta_path: &std::path::Path,
    cluster_names: &[String],
    settings: &SearchSettings<'_>,
) -> Result<Vec<ItemHits>> {
    let mut per_item = vec![ItemHits::default(); items.len()];
    let mut query_file = NamedTempFile::new().context("Cannot create temp query file")?;
    let mut index_of: HashMap<&str, usize> = HashMap::default();
    {
        let mut writer = BufWriter::new(&mut query_file);
        for (k, item) in items.iter().enumerate() {
            writeln!(writer, ">{}", item.name)?;
            writer.write_all(item.seq)?;
            writeln!(writer)?;
            index_of.insert(item.name.as_str(), k);
        }
        writer.flush()?;
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

    let base_args = settings.minimap2_args();
    let target_path = target_file.path().to_string_lossy().to_string();

    let paf = run_minimap2(&base_args, &target_path, query_file.path(), fasta_path)?;
    for line in paf.lines() {
        let Some(fields) = paf_required_fields(line) else {
            continue;
        };
        let Some(&k) = index_of.get(fields[0]) else {
            continue;
        };
        let cluster = &cluster_names[items[k].ci.val()];
        let Some(hit) = hit_from_paf_fields(&fields, fields[0], cluster, fields[5]) else {
            continue;
        };
        if hit.identity() >= settings.min_identity {
            per_item[k].first.push(hit);
        }
    }

    // Realign what a partial hit left uncovered. minimap2 chains the two
    // halves of a gene broken by a short insertion into one chain; base-level
    // alignment then stops at the insertion and drops the rest of the chain,
    // so only one half was reported and the gene could not be seen as
    // FRAGMENTED. Aligned on its own, the missing stretch finds its half.
    if settings.exact_base_alignment {
        let mut remainder_file =
            NamedTempFile::new().context("Cannot create temp remainder query file")?;
        let mut remainders: HashMap<String, (usize, u64, u64)> = HashMap::default();
        {
            let mut writer = BufWriter::new(&mut remainder_file);
            // Sorted, so the remainder queries and minimap2's output come in
            // the same order on every run.
            let mut ordered: Vec<usize> = (0..items.len()).collect();
            ordered.sort_unstable_by(|&a, &b| items[a].name.cmp(&items[b].name));
            for k in ordered {
                let intervals: Vec<(u64, u64)> = per_item[k]
                    .first
                    .iter()
                    .map(|h| (h.query_start, h.query_end))
                    .collect();
                if intervals.is_empty() {
                    continue;
                }
                let seq = items[k].seq;
                let full_len = seq.len() as u64;
                for (from, to) in uncovered_stretches(intervals, full_len) {
                    if to - from < REMAINDER_MIN_BP {
                        continue;
                    }
                    let name = format!("{}~~~REM~{}", items[k].name, from);
                    writeln!(writer, ">{}", name)?;
                    writer.write_all(&seq[from as usize..to as usize])?;
                    writeln!(writer)?;
                    remainders.insert(name, (k, from, full_len));
                }
            }
            writer.flush()?;
        }
        if !remainders.is_empty() {
            let paf = run_minimap2(&base_args, &target_path, remainder_file.path(), fasta_path)?;
            for line in paf.lines() {
                let Some(fields) = paf_required_fields(line) else {
                    continue;
                };
                let Some(&(k, offset, full_len)) = remainders.get(fields[0]) else {
                    continue;
                };
                let cluster = &cluster_names[items[k].ci.val()];
                let Some(mut hit) =
                    hit_from_paf_fields(&fields, &items[k].name, cluster, fields[5])
                else {
                    continue;
                };
                hit.query_start += offset;
                hit.query_end += offset;
                hit.query_length = full_len;
                if hit.identity() >= settings.min_identity {
                    per_item[k].remainder.push(hit);
                }
            }
        }
    }
    Ok(per_item)
}

/// minimap2 results kept between runs: one file per genome and search setting,
/// holding the final hits of every query sequence searched in that genome.
///
/// A search is a function of the genome, the query sequence, the settings that
/// reach minimap2 or filter its output (preset, base-level alignment, the
/// secondary budget, and --min-identity, which also decides what the second
/// pass realigns), and of minimap2 and recore themselves. A file is named after
/// all of those but the query; inside, each query is keyed by its sequence's
/// own hash, and a sequence that found nothing is recorded as searched. A rerun
/// that changes only what happens after the search (the matrix policy, the
/// coverage threshold, synteny) or which gaps are searched (the prevalence
/// floor) aligns only the sequences it has not seen; anything that could change
/// a hit names another file.
pub struct SearchCache {
    dir: std::path::PathBuf,
    settings: String,
}

const SEARCH_CACHE_FORMAT: &str = "recore-search-cache-v1";
const CACHE_FILE_END: &str = "# end";

impl SearchCache {
    fn open(dir: &std::path::Path, settings: &SearchSettings<'_>) -> Result<Self> {
        std::fs::create_dir_all(dir)
            .with_context(|| format!("Cannot create search cache directory {}", dir.display()))?;
        let settings = format!(
            "{} recore={} minimap2={} preset={} exact={} secondary={} min_identity_bits={:016x} peak={} end_bonus={} remainder_min={}",
            SEARCH_CACHE_FORMAT,
            env!("CARGO_PKG_VERSION"),
            minimap2_version()?,
            settings.preset,
            settings.exact_base_alignment,
            settings.max_secondary,
            settings.min_identity.to_bits(),
            MIN_PEAK_SCORE,
            END_BONUS,
            REMAINDER_MIN_BP
        );
        Ok(Self {
            dir: dir.to_path_buf(),
            settings,
        })
    }

    /// The cache file for a genome, named after its bytes and the settings.
    fn file_for(&self, fasta_path: &std::path::Path) -> Result<std::path::PathBuf> {
        use std::io::Read;
        let mut file = std::fs::File::open(fasta_path)
            .with_context(|| format!("Cannot open genome {}", fasta_path.display()))?;
        let mut hasher = Fnv128::new();
        let mut buf = vec![0u8; 1 << 20];
        loop {
            let n = file
                .read(&mut buf)
                .with_context(|| format!("Cannot read genome {}", fasta_path.display()))?;
            if n == 0 {
                break;
            }
            hasher.write(&buf[..n]);
        }
        let genome = hasher.finish();
        let name = stable_hash(&[&genome.to_le_bytes(), self.settings.as_bytes()]);
        Ok(self.dir.join(format!("{name:032x}.tsv")))
    }

    /// The stored hits by sequence hash. A file that cannot be read back is
    /// treated as empty, so its genome is searched again and the file rewritten.
    fn load(&self, path: &std::path::Path) -> HashMap<u128, ItemHits> {
        let Ok(text) = std::fs::read_to_string(path) else {
            return HashMap::default();
        };
        match parse_cache_file(&text, &self.settings) {
            Some(stored) => stored,
            None => {
                log::warn!(
                    "Search cache file {} is unreadable or from other settings; searching that genome again",
                    path.display()
                );
                HashMap::default()
            }
        }
    }

    /// Replace the file in one step, so a run that stops halfway, or another
    /// run on the same cache, never sees half a file.
    fn store(&self, path: &std::path::Path, stored: &HashMap<u128, ItemHits>) -> Result<()> {
        let mut file = NamedTempFile::new_in(&self.dir)
            .with_context(|| format!("Cannot write to search cache {}", self.dir.display()))?;
        {
            let mut w = BufWriter::new(&mut file);
            writeln!(w, "# {}", self.settings)?;
            let mut hashes: Vec<&u128> = stored.keys().collect();
            hashes.sort_unstable();
            for hash in hashes {
                writeln!(w, "S\t{hash:032x}")?;
                let hits = &stored[hash];
                for (pass, list) in [("1", &hits.first), ("2", &hits.remainder)] {
                    for h in list {
                        writeln!(
                            w,
                            "{pass}\t{hash:032x}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            h.query_length,
                            h.query_start,
                            h.query_end,
                            h.strand,
                            h.target,
                            h.target_start,
                            h.target_end,
                            h.target_length,
                            h.matches,
                            h.alignment_length,
                            h.mapq
                        )?;
                    }
                }
            }
            // A file without this last line was cut short.
            writeln!(w, "{CACHE_FILE_END}")?;
            w.flush()?;
        }
        file.persist(path)
            .with_context(|| format!("Cannot write search cache file {}", path.display()))?;
        Ok(())
    }
}

/// A cache file's hits by sequence hash, or None when any line is malformed,
/// the file was cut short, or it was written under other settings. Cut short
/// mid-line, a number could otherwise read back as a smaller one.
fn parse_cache_file(text: &str, settings: &str) -> Option<HashMap<u128, ItemHits>> {
    let body = text.strip_suffix('\n')?.strip_suffix(CACHE_FILE_END)?;
    let mut lines = body.lines();
    if lines.next()?.strip_prefix("# ")? != settings {
        return None;
    }
    let mut stored: HashMap<u128, ItemHits> = HashMap::default();
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        let hash = u128::from_str_radix(fields.get(1)?, 16).ok()?;
        match fields[0] {
            "S" if fields.len() == 2 => {
                stored.entry(hash).or_default();
            }
            pass @ ("1" | "2") if fields.len() == 13 => {
                let mut chars = fields[5].chars();
                let strand = chars.next()?;
                if chars.next().is_some() {
                    return None;
                }
                let hit = Hit {
                    query: String::new(),
                    cluster: String::new(),
                    query_length: fields[2].parse().ok()?,
                    query_start: fields[3].parse().ok()?,
                    query_end: fields[4].parse().ok()?,
                    strand,
                    target: fields[6].to_string(),
                    target_start: fields[7].parse().ok()?,
                    target_end: fields[8].parse().ok()?,
                    target_length: fields[9].parse().ok()?,
                    matches: fields[10].parse().ok()?,
                    alignment_length: fields[11].parse().ok()?,
                    mapq: fields[12].parse().ok()?,
                };
                let entry = stored.get_mut(&hash)?;
                if pass == "1" {
                    entry.first.push(hit);
                } else {
                    entry.remainder.push(hit);
                }
            }
            _ => return None,
        }
    }
    Some(stored)
}

/// minimap2's version, asked once per run: its hits can change between versions.
fn minimap2_version() -> Result<String> {
    static VERSION: std::sync::OnceLock<String> = std::sync::OnceLock::new();
    if let Some(version) = VERSION.get() {
        return Ok(version.clone());
    }
    let output = Command::new("minimap2")
        .arg("--version")
        .stdin(Stdio::null())
        .output()
        .context("Failed to run minimap2: is it installed and on PATH?")?;
    let version = String::from_utf8_lossy(&output.stdout).trim().to_string();
    if !output.status.success() || version.is_empty() {
        anyhow::bail!("minimap2 --version did not report a version");
    }
    Ok(VERSION.get_or_init(|| version).clone())
}

/// 128-bit FNV-1a: stable across platforms and Rust versions, which the
/// standard library's hasher does not promise, and that is what a file name
/// read back by a later run needs.
struct Fnv128(u128);

impl Fnv128 {
    fn new() -> Self {
        Self(0x6c62_272e_07bb_0142_62b8_2175_6295_c58d)
    }

    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.0 ^= u128::from(byte);
            self.0 = self
                .0
                .wrapping_mul(0x0000_0000_0100_0000_0000_0000_0000_013b);
        }
    }

    fn finish(&self) -> u128 {
        self.0
    }
}

/// Hash of the parts in order, each followed by its length so that a boundary
/// cannot move between parts.
fn stable_hash(parts: &[&[u8]]) -> u128 {
    let mut hasher = Fnv128::new();
    for part in parts {
        hasher.write(part);
        hasher.write(&(part.len() as u64).to_le_bytes());
    }
    hasher.finish()
}

/// A stretch of a query this short is not worth a second alignment: it is an
/// end trimmed by alignment, not a missing half.
const REMAINDER_MIN_BP: u64 = 100;

/// The parts of `[0, len)` that no interval covers.
fn uncovered_stretches(mut intervals: Vec<(u64, u64)>, len: u64) -> Vec<(u64, u64)> {
    intervals.sort_unstable();
    let mut gaps = Vec::new();
    let mut reach = 0u64;
    for (start, end) in intervals {
        if start > reach {
            gaps.push((reach, start.min(len)));
        }
        reach = reach.max(end);
    }
    if reach < len {
        gaps.push((reach, len));
    }
    gaps
}

fn run_minimap2(
    base_args: &[&str],
    target_path: &str,
    query_path: &std::path::Path,
    fasta_path: &std::path::Path,
) -> Result<String> {
    let query_path = query_path.to_string_lossy().to_string();
    let output = Command::new("minimap2")
        .args(base_args)
        .arg(target_path)
        .arg(&query_path)
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
    // Target names are the contig tokens recore wrote, so they are UTF-8 and
    // contain no whitespace. Query names are recore's own.
    String::from_utf8(output.stdout).context("minimap2 printed output that is not UTF-8")
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
        target_length,
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

    fn cached_hit(query_start: u64, target: &str) -> Hit {
        Hit {
            query: String::new(),
            cluster: String::new(),
            query_length: 900,
            query_start,
            query_end: 900,
            strand: '-',
            target: target.into(),
            target_start: 10,
            target_end: 910,
            target_length: 5000,
            matches: 899,
            alignment_length: 900,
            mapq: 60,
        }
    }

    /// FNV-1a 128 against its published test vectors: a cache file name has to
    /// be the same on every platform and in every later build.
    #[test]
    fn fnv128_matches_the_reference_vectors() {
        let fnv = |bytes: &[u8]| {
            let mut h = Fnv128::new();
            h.write(bytes);
            h.finish()
        };
        assert_eq!(fnv(b""), 0x6c62272e07bb014262b821756295c58d);
        assert_eq!(fnv(b"a"), 0xd228cb696f1a8caf78912b704e4a8964);
        // A boundary between parts cannot move without changing the hash.
        assert_ne!(stable_hash(&[b"ab", b"c"]), stable_hash(&[b"a", b"bc"]));
    }

    /// A cache file reads back exactly what was stored, sequences with no hit
    /// included, and anything cut short or from other settings reads as nothing.
    #[test]
    fn cache_file_round_trips_and_refuses_damage() {
        let dir = tempfile::tempdir().unwrap();
        let cache = SearchCache {
            dir: dir.path().to_path_buf(),
            settings: "recore-search-cache-v1 test".to_string(),
        };
        let mut stored: HashMap<u128, ItemHits> = HashMap::new();
        stored.insert(
            7,
            ItemHits {
                first: vec![cached_hit(0, "ctg1")],
                remainder: vec![cached_hit(450, "ctg2")],
            },
        );
        stored.insert(9, ItemHits::default());
        let path = dir.path().join("genome.tsv");
        cache.store(&path, &stored).unwrap();

        let back = cache.load(&path);
        assert_eq!(back.len(), 2);
        assert!(back[&9].first.is_empty() && back[&9].remainder.is_empty());
        let hits = &back[&7];
        assert_eq!(hits.first.len(), 1);
        assert_eq!(hits.remainder.len(), 1);
        let (a, b) = (&hits.first[0], &stored[&7].first[0]);
        assert_eq!(
            (
                a.query_start,
                a.strand,
                a.target.as_str(),
                a.target_length,
                a.matches,
                a.mapq
            ),
            (
                b.query_start,
                b.strand,
                b.target.as_str(),
                b.target_length,
                b.matches,
                b.mapq
            )
        );
        assert_eq!(hits.remainder[0].query_start, 450);

        let text = std::fs::read_to_string(&path).unwrap();
        // Cut short anywhere, even on a line boundary, it is refused.
        assert!(parse_cache_file(&text[..text.len() / 2], &cache.settings).is_none());
        let without_end = text.trim_end_matches("# end\n");
        assert!(parse_cache_file(without_end, &cache.settings).is_none());
        // Written under other settings, it is not this cache's.
        assert!(parse_cache_file(&text, "recore-search-cache-v1 other").is_none());
        // A hit for a sequence never marked as searched is damage too.
        let orphan = text.replace("S\t00000000000000000000000000000007\n", "");
        assert_ne!(orphan, text);
        assert!(parse_cache_file(&orphan, &cache.settings).is_none());
    }

    /// Cached or not, hits come out in the order one search of every query
    /// reports them: first alignments in query order, remainders by name.
    #[test]
    fn hits_come_out_in_search_order() {
        let seq = b"ACGT".to_vec();
        let items = vec![
            QueryItem {
                name: "b".into(),
                ci: ClusterIdx(0),
                seq: &seq,
            },
            QueryItem {
                name: "a".into(),
                ci: ClusterIdx(1),
                seq: &seq,
            },
        ];
        let named = |start| Hit {
            query: String::new(),
            ..cached_hit(start, "t")
        };
        let per_item = vec![
            ItemHits {
                first: vec![named(1)],
                remainder: vec![named(2)],
            },
            ItemHits {
                first: vec![named(3)],
                remainder: vec![named(4)],
            },
        ];
        let order: Vec<(usize, u64)> = in_search_order(&items, &per_item)
            .into_iter()
            .map(|(ci, hit)| (ci.val(), hit.query_start))
            .collect();
        assert_eq!(order, vec![(0, 1), (1, 3), (1, 4), (0, 2)]);
    }

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
            target_length: 0,
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
            target_length: 0,
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
            target_length: 0,
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
            target_length: 0,
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
            target_length: 0,
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

    #[test]
    fn test_uncovered_stretches() {
        assert_eq!(uncovered_stretches(vec![(0, 516)], 1035), vec![(516, 1035)]);
        assert_eq!(uncovered_stretches(vec![(605, 1215)], 1215), vec![(0, 605)]);
        assert_eq!(
            uncovered_stretches(vec![(100, 400), (300, 600)], 1000),
            vec![(0, 100), (600, 1000)]
        );
        assert!(uncovered_stretches(vec![(0, 1000)], 1000).is_empty());
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
