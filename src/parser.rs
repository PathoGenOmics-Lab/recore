use ahash::{HashMap, HashMapExt, HashSet, HashSetExt};
use anyhow::{Context, Result};
use log::info;
use std::fs;
use std::path::{Path, PathBuf};

use crate::types::{ClusterIdx, GenomeIdx};

// ── Presence/absence matrix (encapsulated) ──

/// Presence/absence matrix from Panaroo/Roary output.
///
/// Fields are private to enforce consistent access and prevent
/// accidental mutation of internal indices.
pub struct PresenceAbsence {
    clusters: Vec<String>,
    genomes: Vec<String>,
    present_genomes: Vec<Vec<GenomeIdx>>,
    present_bits: Vec<Vec<bool>>,
    cell_values: Vec<Vec<String>>,
    header: Vec<String>,
    metadata: Vec<Vec<String>>,
    genome_col_indices: Vec<usize>,
}

impl PresenceAbsence {
    // ── Dimensions ──

    #[inline]
    pub fn n_clusters(&self) -> usize {
        self.clusters.len()
    }

    #[inline]
    pub fn n_genomes(&self) -> usize {
        self.genomes.len()
    }

    // ── Name lookups ──

    #[inline]
    pub fn cluster_name(&self, ci: ClusterIdx) -> &str {
        &self.clusters[ci.val()]
    }

    #[inline]
    pub fn genome_name(&self, gi: GenomeIdx) -> &str {
        &self.genomes[gi.val()]
    }

    pub fn cluster_names(&self) -> &[String] {
        &self.clusters
    }

    pub fn genome_names(&self) -> &[String] {
        &self.genomes
    }

    // ── Presence queries ──

    #[inline]
    pub fn is_present(&self, ci: ClusterIdx, gi: GenomeIdx) -> bool {
        self.present_bits[ci.val()][gi.val()]
    }

    #[inline]
    pub fn present_count(&self, ci: ClusterIdx) -> usize {
        self.present_genomes[ci.val()].len()
    }

    pub fn present_genomes(&self, ci: ClusterIdx) -> &[GenomeIdx] {
        &self.present_genomes[ci.val()]
    }

    pub fn presence_row(&self, ci: ClusterIdx) -> &[bool] {
        &self.present_bits[ci.val()]
    }

    // ── Cell values (for gene IDs and CSV reconstruction) ──

    #[inline]
    pub fn cell(&self, ci: ClusterIdx, gi: GenomeIdx) -> &str {
        &self.cell_values[ci.val()][gi.val()]
    }

    // ── CSV structure (for output reconstruction) ──

    pub fn header(&self) -> &[String] {
        &self.header
    }

    pub fn metadata(&self, ci: ClusterIdx) -> &[String] {
        &self.metadata[ci.val()]
    }

    pub fn genome_col_indices(&self) -> &[usize] {
        &self.genome_col_indices
    }

    // ── Derived statistics ──

    #[allow(dead_code)]
    pub fn gaps(&self) -> Vec<(ClusterIdx, GenomeIdx)> {
        let mut gaps = Vec::new();
        for ci in 0..self.clusters.len() {
            let ci = ClusterIdx(ci);
            for gi in 0..self.genomes.len() {
                let gi = GenomeIdx(gi);
                if !self.is_present(ci, gi) {
                    gaps.push((ci, gi));
                }
            }
        }
        gaps
    }

    #[allow(dead_code)]
    pub fn prevalence(&self, ci: ClusterIdx) -> f64 {
        if self.genomes.is_empty() {
            return 0.0;
        }
        self.present_count(ci) as f64 / self.genomes.len() as f64
    }

    pub fn core_size(&self) -> usize {
        let n = self.genomes.len();
        self.present_genomes.iter().filter(|p| p.len() == n).count()
    }

    pub fn soft_core_size(&self) -> usize {
        let threshold = (self.genomes.len() as f64 * 0.95).ceil() as usize;
        self.present_genomes
            .iter()
            .filter(|p| p.len() >= threshold)
            .count()
    }

    /// Constructor for tests only — minimal PA with all genes absent.
    #[cfg(test)]
    pub fn new_for_test(
        clusters: Vec<String>,
        genomes: Vec<String>,
        cell_values: Vec<Vec<String>>,
    ) -> Self {
        let n_cl = clusters.len();
        let n_ge = genomes.len();
        let mut present_genomes = Vec::with_capacity(n_cl);
        let mut present_bits = Vec::with_capacity(n_cl);
        for row in &cell_values {
            let mut present_row = Vec::new();
            let mut bit_row = vec![false; n_ge];
            for (gi, v) in row.iter().enumerate() {
                if !v.is_empty() {
                    present_row.push(GenomeIdx(gi));
                    bit_row[gi] = true;
                }
            }
            present_genomes.push(present_row);
            present_bits.push(bit_row);
        }
        let metadata = vec![vec![]; n_cl];
        let genome_col_indices: Vec<usize> = (0..n_ge).collect();
        let mut header = vec!["Gene".to_string()];
        header.extend(genomes.iter().cloned());
        PresenceAbsence {
            clusters,
            genomes,
            present_genomes,
            present_bits,
            cell_values,
            header,
            metadata,
            genome_col_indices,
        }
    }

    /// Test constructor that also carries Panaroo's Gene and Annotation
    /// metadata columns, so the report sees annotations the way a real run does.
    #[cfg(test)]
    pub fn new_for_test_with_annotations(
        clusters: Vec<String>,
        annotations: Vec<String>,
        genomes: Vec<String>,
        cell_values: Vec<Vec<String>>,
    ) -> Self {
        let mut pa = Self::new_for_test(clusters, genomes, cell_values);
        pa.metadata = pa
            .clusters
            .iter()
            .zip(annotations)
            .map(|(name, annotation)| vec![name.clone(), annotation])
            .collect();
        let mut header = vec!["Gene".to_string(), "Annotation".to_string()];
        header.extend(pa.genomes.iter().cloned());
        pa.header = header;
        pa.genome_col_indices = (2..pa.genomes.len() + 2).collect();
        pa
    }

    /// Build cluster_name → ClusterIdx index.
    pub fn cluster_index(&self) -> HashMap<String, ClusterIdx> {
        self.clusters
            .iter()
            .enumerate()
            .map(|(i, c)| (c.clone(), ClusterIdx(i)))
            .collect()
    }
}

// ── Panaroo CSV parsing ──

const METADATA_COLS: &[&str] = &[
    "Gene",
    "Non-unique Gene name",
    "Annotation",
    "No. isolates",
    "No. sequences",
    "Avg sequences per isolate",
    "Genome Fragment",
    "Order within Fragment",
    "Accessory Fragment",
    "Accessory Order with Fragment",
    "Accessory Order",
    "QC",
    "Min group size nuc",
    "Max group size nuc",
    "Avg group size nuc",
];

pub fn parse_presence_absence(path: &Path) -> Result<PresenceAbsence> {
    let mut reader =
        csv::Reader::from_path(path).with_context(|| format!("Cannot open {}", path.display()))?;

    let headers: Vec<String> = reader.headers()?.iter().map(|s| s.to_string()).collect();

    let meta_set: HashSet<&str> = METADATA_COLS.iter().copied().collect();
    let gene_col = headers
        .iter()
        .position(|h| h == "Gene")
        .context("Missing 'Gene' column in presence/absence CSV")?;

    let genome_indices: Vec<usize> = headers
        .iter()
        .enumerate()
        .filter(|(_, h)| !meta_set.contains(h.as_str()))
        .map(|(i, _)| i)
        .collect();

    let genomes: Vec<String> = genome_indices.iter().map(|&i| headers[i].clone()).collect();

    let genome_idx_set: HashSet<usize> = genome_indices.iter().copied().collect();
    let meta_col_indices: Vec<usize> = (0..headers.len())
        .filter(|i| !genome_idx_set.contains(i))
        .collect();

    let mut clusters = Vec::new();
    let mut present_genomes = Vec::new();
    let mut present_bits = Vec::new();
    let mut cell_values = Vec::new();
    let mut metadata = Vec::new();

    for result in reader.records() {
        let record = result?;
        let cluster_name = record.get(gene_col).unwrap_or("").to_string();

        let mut present_row = Vec::new();
        let mut bit_row = vec![false; genomes.len()];
        let mut row_values = Vec::with_capacity(genomes.len());

        for (gi, &col_idx) in genome_indices.iter().enumerate() {
            let val = record.get(col_idx).unwrap_or("").trim().to_string();
            if !val.is_empty() {
                present_row.push(GenomeIdx(gi));
                bit_row[gi] = true;
            }
            row_values.push(val);
        }

        let meta_values: Vec<String> = meta_col_indices
            .iter()
            .map(|&i| record.get(i).unwrap_or("").to_string())
            .collect();

        clusters.push(cluster_name);
        present_genomes.push(present_row);
        present_bits.push(bit_row);
        cell_values.push(row_values);
        metadata.push(meta_values);
    }

    info!(
        "Parsed presence/absence: {} clusters × {} genomes (core: {}, soft-core 95%: {})",
        clusters.len(),
        genomes.len(),
        {
            let n = genomes.len();
            present_genomes.iter().filter(|p| p.len() == n).count()
        },
        {
            let t = (genomes.len() as f64 * 0.95).ceil() as usize;
            present_genomes.iter().filter(|p| p.len() >= t).count()
        }
    );

    Ok(PresenceAbsence {
        clusters,
        genomes,
        present_genomes,
        present_bits,
        cell_values,
        header: headers,
        metadata,
        genome_col_indices: genome_indices,
    })
}

// ── FASTA / GFF parsing ──

/// Parse a reference FASTA keyed by the first token of each header.
///
/// For the pan-genome reference use [`parse_pan_reference_for`], which resolves
/// headers against the matrix's cluster names; this keyed-as-written form is
/// for references with free names, such as `--mobile-element-reference`.
pub fn parse_pan_reference(
    path: &Path,
    multi_allele: bool,
) -> Result<HashMap<String, Vec<Vec<u8>>>> {
    parse_reference_records(path, multi_allele, None)
}

/// Parse the pan-genome reference, naming each record after the matrix cluster
/// it belongs to.
///
/// Three header layouts are in use, and each record is resolved against the
/// matrix's own cluster names in this order:
///
/// 1. `>cluster ...`, Panaroo's layout: the first token is the cluster.
/// 2. `>seqid cluster`, Roary's layout: the second token is the cluster. This
///    used to match nothing, and a Roary run exited 0 with no rescues at all.
/// 3. `>cluster~~~ALLELE~n`: allele naming, where the cluster is the part before
///    `~~~`. Only as a fallback: Panaroo writes `~~~` inside real cluster names
///    when it merges genes (`pstB~~~phoT`), and cutting every name at `~~~`
///    either lost those clusters or, worse, filed them under another cluster
///    (`pstB~~~phoT` replaced `pstB`'s own sequence).
///
/// A header that matches none of these keeps its first token, and is reported
/// as not being in the matrix.
pub fn parse_pan_reference_for(
    path: &Path,
    multi_allele: bool,
    clusters: &HashSet<&str>,
) -> Result<HashMap<String, Vec<Vec<u8>>>> {
    parse_reference_records(path, multi_allele, Some(clusters))
}

fn resolve_reference_name(header: &str, clusters: Option<&HashSet<&str>>) -> String {
    let mut tokens = header.split_whitespace();
    let first = tokens.next().unwrap_or(header);
    let Some(clusters) = clusters else {
        return first.to_string();
    };
    if clusters.contains(first) {
        return first.to_string();
    }
    if let Some(second) = tokens.next() {
        if clusters.contains(second) {
            return second.to_string();
        }
    }
    if let Some(pos) = first.find("~~~") {
        let prefix = &first[..pos];
        if clusters.contains(prefix) {
            return prefix.to_string();
        }
    }
    first.to_string()
}

fn parse_reference_records(
    path: &Path,
    multi_allele: bool,
    clusters: Option<&HashSet<&str>>,
) -> Result<HashMap<String, Vec<Vec<u8>>>> {
    use std::io::{BufRead, BufReader};

    let file = fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    let mut reader = BufReader::new(file);

    let mut refs: HashMap<String, Vec<Vec<u8>>> = HashMap::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();
    let mut records = 0usize;

    fn flush_record(
        name: String,
        seq: Vec<u8>,
        refs: &mut HashMap<String, Vec<Vec<u8>>>,
        multi_allele: bool,
    ) {
        if seq.is_empty() {
            return;
        }
        let entry = refs.entry(name).or_default();
        if multi_allele {
            entry.push(seq);
        } else if entry.is_empty() || entry[0].len() < seq.len() {
            if entry.is_empty() {
                entry.push(seq);
            } else {
                entry[0] = seq;
            }
        }
    }

    // Bytes, not lines: only the name tokens have to be UTF-8, and a stray
    // Latin-1 byte in a description used to abort the run with a bare
    // "stream did not contain valid UTF-8" that named no file.
    let mut line = Vec::with_capacity(4096);
    let mut line_no = 0usize;
    loop {
        line.clear();
        if reader
            .read_until(b'\n', &mut line)
            .with_context(|| format!("Cannot read {}", path.display()))?
            == 0
        {
            break;
        }
        line_no += 1;
        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        if let Some(header) = line.strip_prefix(b">") {
            if let Some(name) = current_name.take() {
                flush_record(
                    name,
                    std::mem::take(&mut current_seq),
                    &mut refs,
                    multi_allele,
                );
            }
            current_seq.clear();
            records += 1;
            // The name is in the first two tokens; anything after them is free
            // text and may be in any encoding.
            let names: Vec<&[u8]> = header
                .split(|b| b.is_ascii_whitespace())
                .filter(|t| !t.is_empty())
                .take(2)
                .collect();
            let mut decoded = Vec::with_capacity(names.len());
            for token in names {
                match std::str::from_utf8(token) {
                    Ok(t) => decoded.push(t),
                    Err(_) if !decoded.is_empty() => break,
                    Err(_) => anyhow::bail!(
                        "{}: line {}: sequence name is not valid UTF-8",
                        path.display(),
                        line_no
                    ),
                }
            }
            if decoded.is_empty() {
                anyhow::bail!(
                    "{}: line {}: FASTA header with no name",
                    path.display(),
                    line_no
                );
            }
            current_name = Some(resolve_reference_name(&decoded.join(" "), clusters));
        } else {
            current_seq.extend(line.iter().filter(|b| !b.is_ascii_whitespace()));
        }
    }

    if let Some(name) = current_name {
        flush_record(name, current_seq, &mut refs, multi_allele);
    }

    let total_alleles: usize = refs.values().map(|v| v.len()).sum();
    info!(
        "Parsed {} clusters ({} total alleles, {} records) from {}",
        refs.len(),
        total_alleles,
        records,
        path.display()
    );
    Ok(refs)
}

pub fn find_genome_files(dir: &Path) -> Result<HashMap<String, PathBuf>> {
    let extensions = [".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz"];
    let mut files: HashMap<String, PathBuf> = HashMap::new();
    let mut as_id_candidates: HashMap<String, PathBuf> = HashMap::new();
    let mut ambiguous_as_ids: HashSet<String> = HashSet::new();

    for entry in
        fs::read_dir(dir).with_context(|| format!("Cannot read directory {}", dir.display()))?
    {
        let entry = entry?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }

        let Some(fname) = path.file_name() else {
            continue;
        };
        let filename = fname.to_string_lossy().to_string();

        let mut name = filename.clone();
        for ext in &extensions {
            if name.ends_with(ext) {
                name = name[..name.len() - ext.len()].to_string();
                break;
            }
        }

        if name != filename {
            files.entry(name.clone()).or_insert_with(|| path.clone());
            if let Some(as_id) = extract_as_id(&name) {
                if as_id != name {
                    match as_id_candidates.get(as_id) {
                        Some(existing) if existing != &path => {
                            ambiguous_as_ids.insert(as_id.to_string());
                        }
                        Some(_) => {}
                        None => {
                            as_id_candidates.insert(as_id.to_string(), path.clone());
                        }
                    }
                }
            }
        }
    }

    for (as_id, path) in as_id_candidates {
        if !ambiguous_as_ids.contains(&as_id) {
            files.entry(as_id).or_insert(path);
        }
    }

    // Count files, not keys: an AS identifier is a second key for the same file,
    // and counting keys reported twice the real number.
    let distinct: HashSet<&PathBuf> = files.values().collect();
    info!("Found {} genome files in {}", distinct.len(), dir.display());
    Ok(files)
}

/// Extract the canonical AS identifier embedded in a genome or file name.
pub fn extract_as_id(name: &str) -> Option<&str> {
    let bytes = name.as_bytes();
    if bytes.len() < 7 {
        return None;
    }
    for i in 0..=bytes.len() - 7 {
        if bytes[i] == b'A'
            && bytes[i + 1] == b'S'
            && bytes[i + 2..i + 7].iter().all(u8::is_ascii_digit)
            && bytes.get(i + 7).map_or(true, |ch| !ch.is_ascii_digit())
        {
            return Some(&name[i..i + 7]);
        }
    }
    None
}

/// Resolve a genome FASTA by exact column name, then by embedded AS identifier.
pub fn resolve_genome_file<'a>(
    genome_files: &'a HashMap<String, PathBuf>,
    genome_name: &str,
) -> Option<&'a PathBuf> {
    genome_files
        .get(genome_name)
        .or_else(|| extract_as_id(genome_name).and_then(|as_id| genome_files.get(as_id)))
}

/// AS identifiers carried by more than one matrix column.
fn shared_as_ids(names: &[String]) -> HashSet<&str> {
    let mut seen: HashSet<&str> = HashSet::new();
    let mut shared: HashSet<&str> = HashSet::new();
    for as_id in names.iter().filter_map(|n| extract_as_id(n)) {
        if !seen.insert(as_id) {
            shared.insert(as_id);
        }
    }
    shared
}

/// Whether a column may take a file through its AS identifier rather than its
/// own name, and if not, why. The fallback used to take any file carrying the
/// identifier, so a column with no assembly of its own was given another
/// column's, and recore reported rescues from the wrong genome with no warning.
fn as_fallback_refusal(
    name: &str,
    as_id: &str,
    path: &Path,
    shared: &HashSet<&str>,
    exact_owner: &HashMap<PathBuf, &str>,
) -> Option<String> {
    if shared.contains(as_id) {
        return Some(format!(
            "the identifier {} is carried by more than one matrix column",
            as_id
        ));
    }
    if let Some(owner) = exact_owner.get(path) {
        if *owner != name {
            return Some(format!("it is column {}'s own file", owner));
        }
    }
    None
}

/// Build a map keyed by presence/absence genome names, with AS-ID fallback.
pub fn resolve_genome_files_for_pa(
    pa: &PresenceAbsence,
    genome_files: &HashMap<String, PathBuf>,
) -> HashMap<String, PathBuf> {
    use log::warn;
    let names = pa.genome_names();
    let shared = shared_as_ids(names);
    let exact_owner: HashMap<PathBuf, &str> = names
        .iter()
        .filter_map(|n| genome_files.get(n).map(|p| (p.clone(), n.as_str())))
        .collect();
    let mut resolved = HashMap::new();
    for genome_name in names {
        if let Some(path) = genome_files.get(genome_name) {
            resolved.insert(genome_name.clone(), path.clone());
            continue;
        }
        let Some(as_id) = extract_as_id(genome_name) else {
            continue;
        };
        let Some(path) = genome_files.get(as_id) else {
            continue;
        };
        if let Some(reason) = as_fallback_refusal(genome_name, as_id, path, &shared, &exact_owner) {
            warn!(
                "Genome {} has no FASTA of its own; not using {} through its AS identifier, because {}",
                genome_name,
                path.display(),
                reason
            );
            continue;
        }
        resolved.insert(genome_name.clone(), path.clone());
    }
    resolved
}

/// The GFF for each matrix column, by exact name or, under the same rules as
/// for FASTA files, by embedded AS identifier.
pub fn resolve_annotation_files_for_pa(
    pa: &PresenceAbsence,
    ann_dir: &Path,
) -> HashMap<String, PathBuf> {
    use log::warn;
    let names = pa.genome_names();
    let shared = shared_as_ids(names);
    let exact = |name: &str| {
        [".gff", ".gff3", ".gff.gz", ".gff3.gz"]
            .iter()
            .map(|suffix| ann_dir.join(format!("{}{}", name, suffix)))
            .find(|path| path.exists())
    };
    let mut resolved: HashMap<String, PathBuf> = HashMap::new();
    let mut exact_owner: HashMap<PathBuf, &str> = HashMap::new();
    for name in names {
        if let Some(path) = exact(name) {
            exact_owner.insert(path.clone(), name.as_str());
            resolved.insert(name.clone(), path);
        }
    }
    for name in names {
        if resolved.contains_key(name) {
            continue;
        }
        let (Some(as_id), Some(path)) = (extract_as_id(name), find_annotation_file(ann_dir, name))
        else {
            continue;
        };
        if let Some(reason) = as_fallback_refusal(name, as_id, &path, &shared, &exact_owner) {
            warn!(
                "Genome {} has no GFF of its own; not using {} through its AS identifier, because {}",
                name,
                path.display(),
                reason
            );
            continue;
        }
        resolved.insert(name.clone(), path);
    }
    resolved
}

/// Resolve a GFF file by exact genome name, then by embedded AS identifier.
pub fn find_annotation_file(ann_dir: &Path, genome_name: &str) -> Option<PathBuf> {
    let mut names = vec![genome_name.to_string()];
    if let Some(as_id) = extract_as_id(genome_name) {
        if as_id != genome_name {
            names.push(as_id.to_string());
        }
    }

    for name in names {
        for suffix in [".gff", ".gff3", ".gff.gz", ".gff3.gz"] {
            let path = ann_dir.join(format!("{}{}", name, suffix));
            if path.exists() {
                return Some(path);
            }
        }
    }

    let as_id = extract_as_id(genome_name)?;
    let mut matches = Vec::new();
    let entries = fs::read_dir(ann_dir).ok()?;
    for entry in entries.flatten() {
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let filename = path.file_name()?.to_string_lossy();
        let valid_suffix = [".gff", ".gff3", ".gff.gz", ".gff3.gz"]
            .iter()
            .any(|suffix| filename.ends_with(suffix));
        if valid_suffix && extract_as_id(&filename).is_some_and(|candidate| candidate == as_id) {
            matches.push(path);
        }
    }
    if matches.len() == 1 {
        return matches.pop();
    }
    None
}

/// The contig name of a FASTA header line (without its '>'): the first
/// whitespace-delimited token, which has to be UTF-8. The rest of the line is
/// free text in any encoding. Matches what the search target writer accepts,
/// so every reader of a genome agrees on its contig names.
fn fasta_contig_name(header: &[u8], path: &Path, line_no: usize) -> Result<String> {
    let token = header
        .split(|b| b.is_ascii_whitespace())
        .next()
        .unwrap_or_default();
    if token.is_empty() {
        anyhow::bail!(
            "{}: line {}: FASTA header with no contig name",
            path.display(),
            line_no
        );
    }
    std::str::from_utf8(token).map(str::to_string).map_err(|_| {
        anyhow::anyhow!(
            "{}: line {}: contig name is not valid UTF-8",
            path.display(),
            line_no
        )
    })
}

/// Visit each line of a FASTA as bytes, with line endings removed.
///
/// Bytes rather than `lines()`: one Latin-1 byte in a header's description used
/// to fail the whole read, and the internal-stop scan, which skipped a genome it
/// could not read, then let pseudogenes into the matrix as intact genes.
fn for_each_fasta_line(
    path: &Path,
    mut visit: impl FnMut(&[u8], usize) -> Result<()>,
) -> Result<()> {
    use std::io::BufRead;
    let mut reader = open_maybe_gzip(path)?;
    let mut line = Vec::with_capacity(8192);
    let mut line_no = 0usize;
    loop {
        line.clear();
        if reader
            .read_until(b'\n', &mut line)
            .with_context(|| format!("Cannot read {}", path.display()))?
            == 0
        {
            return Ok(());
        }
        line_no += 1;
        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        visit(&line, line_no)?;
    }
}

/// Load genome FASTA as contig_name → sequence map.
/// Used for internal stop codon scanning in pseudogene detection.
pub fn load_genome_contigs(path: &Path) -> Result<HashMap<String, Vec<u8>>> {
    let mut contigs: HashMap<String, Vec<u8>> = HashMap::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = Vec::new();
    for_each_fasta_line(path, |line, line_no| {
        if let Some(header) = line.strip_prefix(b">") {
            if let Some(name) = current_name.take() {
                contigs.insert(name, std::mem::take(&mut current_seq));
            }
            current_name = Some(fasta_contig_name(header, path, line_no)?);
        } else {
            current_seq.extend(line.iter().filter(|b| !b.is_ascii_whitespace()));
        }
        Ok(())
    })?;
    if let Some(name) = current_name {
        contigs.insert(name, current_seq);
    }
    Ok(contigs)
}

/// Load only FASTA contig names, preserving file order.
pub fn load_genome_contig_names(path: &Path) -> Result<Vec<String>> {
    let mut names = Vec::new();
    for_each_fasta_line(path, |line, line_no| {
        if let Some(header) = line.strip_prefix(b">") {
            names.push(fasta_contig_name(header, path, line_no)?);
        }
        Ok(())
    })?;
    Ok(names)
}

pub type GffRecord = (String, u64, u64, char, String, Vec<String>, Vec<String>);

/// If a single-contig GFF and single-contig FASTA use different contig names,
/// rewrite GFF seqids to the FASTA contig name used by minimap2 hits.
pub fn normalize_gff_records_to_fasta(
    records: &mut [GffRecord],
    fasta_contigs: &[String],
) -> Option<(String, String, usize)> {
    if records.is_empty() || fasta_contigs.len() != 1 {
        return None;
    }

    let mut gff_contigs = ahash::HashSet::new();
    for record in records.iter() {
        gff_contigs.insert(record.0.clone());
        if gff_contigs.len() > 1 {
            return None;
        }
    }

    let from = gff_contigs.into_iter().next()?;
    let to = fasta_contigs[0].clone();
    if from == to {
        return None;
    }

    let mut changed = 0usize;
    for record in records.iter_mut() {
        record.0.clone_from(&to);
        changed += 1;
    }
    Some((from, to, changed))
}

/// Open a file with transparent gzip decompression.
/// If the path ends with `.gz`, decompresses on-the-fly via flate2.
///
/// This is the only place in the crate that decides whether a file is
/// compressed. It used to be one of three, each with its own copy of the
/// extension test and its own decoder, and two of the three were wrong in the
/// same way for the same reason.
///
/// `MultiGzDecoder`, not `GzDecoder`. A gzip file is allowed to be a run of
/// concatenated members, and BGZF, what bgzip, samtools, bcftools and tabix
/// write, is exactly that: one member per 64 KB block, so a real bgzipped
/// genome holds hundreds of them. `GzDecoder` decodes the first member and then
/// reports a clean EOF, which handed this tool a fraction of the assembly and
/// no error at all. Measured on a bgzipped FASTA: 65,259 bytes of 848,890, and
/// a run that exited 0 with a summary asserting a core genome it had never
/// seen. A gene absent because its contig was never read is precisely the
/// artefact Recore exists to remove.
///
/// The cost is real and worth stating rather than discovering. This is stricter
/// than what came before: `MultiGzDecoder` refuses trailing bytes after the
/// final member, where `GzDecoder` ignored them. NUL padding, a stray newline,
/// appended junk, and one byte is enough. Such a file was read yesterday and
/// now fails. A user who wants to know in advance can ask `gzip -t`, which
/// prints "trailing garbage ignored" and exits 2, so the check does catch it.
/// Python's `gzip` is the one that will not warn them: it reads NUL padding in
/// silence and only raises on the junk-suffixed shape.
///
/// That is the intended direction even so. A hard error beats a silently
/// partial genome in a tool whose output is a scientific claim, and the failure
/// now names the file that caused it. If tolerating trailing bytes ever matters
/// more than refusing them, the shape that buys both is a loop over
/// `flate2::bufread::GzDecoder` that peeks for the gzip magic at each member
/// boundary, at the price of owning a decoder loop here.
pub(crate) fn open_maybe_gzip(path: &Path) -> Result<Box<dyn std::io::BufRead>> {
    use std::io::BufReader;
    let file = fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    if path.extension().is_some_and(|e| e == "gz") {
        let decoder = flate2::read::MultiGzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

pub fn parse_gff_genes(path: &Path) -> Result<Vec<GffRecord>> {
    use std::io::BufRead;

    let mut reader = open_maybe_gzip(path)?;

    let mut genes = Vec::new();
    // Bytes, decoded lossily: one non-UTF-8 byte in a free-text attribute such
    // as product= used to fail the whole file, and the genome silently lost its
    // annotation, its synteny and its corrected GFF records.
    let mut buf = Vec::with_capacity(4096);
    loop {
        buf.clear();
        if reader
            .read_until(b'\n', &mut buf)
            .with_context(|| format!("Cannot read {}", path.display()))?
            == 0
        {
            break;
        }
        while matches!(buf.last(), Some(b'\n' | b'\r')) {
            buf.pop();
        }
        let line = String::from_utf8_lossy(&buf);
        if line.starts_with('#') || !line.contains('\t') {
            continue;
        }

        let mut fields = line.split('\t');
        let Some(contig_field) = fields.next() else {
            continue;
        };
        let _source = fields.next();
        let Some(feature_type) = fields.next() else {
            continue;
        };
        if feature_type != "CDS" {
            continue;
        }
        let Some(start_field) = fields.next() else {
            continue;
        };
        let Some(end_field) = fields.next() else {
            continue;
        };
        let _score = fields.next();
        let Some(strand_field) = fields.next() else {
            continue;
        };
        let _phase = fields.next();
        let Some(attributes) = fields.next() else {
            continue;
        };

        let contig = contig_field.to_string();
        // GFF3 is 1-based closed; convert to 0-based half-open to match minimap2 PAF coords
        let Ok(start_1based) = start_field.parse::<u64>() else {
            continue;
        };
        let Ok(end) = end_field.parse::<u64>() else {
            continue;
        };
        if start_1based == 0 || end < start_1based {
            continue;
        }
        let start = start_1based - 1;
        let strand = strand_field.chars().next().unwrap_or('+');

        let mut id = String::new();
        let mut locus_tag = String::new();
        let mut gene_name = String::new();
        let mut name = String::new();
        let mut protein_id = String::new();
        let mut old_locus_tag = String::new();
        let mut dbxref = String::new();
        let mut product = String::new();
        let mut note = String::new();
        for attr in attributes.split(';') {
            let attr = attr.trim();
            if let Some(val) = attr.strip_prefix("ID=") {
                id = decode_gff3_value(val);
            } else if let Some(val) = attr.strip_prefix("locus_tag=") {
                locus_tag = decode_gff3_value(val);
            } else if let Some(val) = attr.strip_prefix("gene=") {
                gene_name = decode_gff3_value(val);
            } else if let Some(val) = attr.strip_prefix("Name=") {
                name = decode_gff3_value(val);
            } else if let Some(val) = attr.strip_prefix("protein_id=") {
                protein_id = decode_gff3_value(val);
            } else if let Some(val) = attr.strip_prefix("old_locus_tag=") {
                old_locus_tag = decode_gff3_value(val);
            } else if let Some(val) = attr.strip_prefix("Dbxref=") {
                dbxref = decode_gff3_value(val);
            } else if let Some(val) = attr.strip_prefix("product=") {
                product = decode_gff3_value(val);
            } else if let Some(val) = attr.strip_prefix("Note=") {
                note = decode_gff3_value(val);
            }
        }

        let gene_id = [
            locus_tag.as_str(),
            gene_name.as_str(),
            name.as_str(),
            id.as_str(),
            protein_id.as_str(),
            old_locus_tag.as_str(),
        ]
        .into_iter()
        .find(|v| !v.is_empty())
        .unwrap_or("")
        .to_string();

        let mut aliases = Vec::new();
        let dbxref_aliases: Vec<String> = dbxref
            .split(',')
            .filter_map(|entry| {
                entry
                    .trim()
                    .strip_prefix("AMAP_ids:")
                    .or_else(|| entry.trim().strip_prefix("GeneID:"))
                    .or_else(|| entry.trim().strip_prefix("locus_tag:"))
                    .map(str::to_string)
            })
            .collect();
        for value in [
            id.as_str(),
            locus_tag.as_str(),
            gene_name.as_str(),
            name.as_str(),
            protein_id.as_str(),
            old_locus_tag.as_str(),
        ] {
            if !value.is_empty() && !aliases.iter().any(|x: &String| x == value) {
                aliases.push(value.to_string());
            }
        }
        for value in dbxref_aliases {
            if !value.is_empty() && !aliases.iter().any(|x: &String| x == &value) {
                aliases.push(value);
            }
        }

        let mut labels = aliases.clone();
        for value in [product.as_str(), note.as_str()] {
            if !value.is_empty() && !labels.iter().any(|x: &String| x == value) {
                labels.push(value.to_string());
            }
        }

        genes.push((contig, start, end, strand, gene_id, aliases, labels));
    }

    Ok(genes)
}

fn decode_gff3_value(value: &str) -> String {
    let bytes = value.as_bytes();
    let mut out = Vec::with_capacity(value.len());
    let mut i = 0usize;
    while i < bytes.len() {
        if bytes[i] == b'%' && i + 2 < bytes.len() {
            let hi = bytes[i + 1];
            let lo = bytes[i + 2];
            if let (Some(hi), Some(lo)) = (hex_value(hi), hex_value(lo)) {
                out.push(hi << 4 | lo);
                i += 3;
                continue;
            }
        }
        out.push(bytes[i]);
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

fn hex_value(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        b'A'..=b'F' => Some(byte - b'A' + 10),
        _ => None,
    }
}

/// Check if a name is safe to use in file path construction.
/// Rejects names containing path separators, parent directory references, or null bytes.
fn is_safe_name(name: &str) -> bool {
    !name.is_empty()
        && !name.contains('/')
        && !name.contains('\\')
        && !name.contains('\0')
        && name != "."
        && name != ".."
        && !name.starts_with("..") // e.g. "../etc/passwd"
        && !name.contains("/../")
        && !name.contains("\\..\\")
}

/// Cross-validate inputs: warn about mismatches between CSV clusters and pan-reference.
/// Also validates that genome/cluster names are safe for path construction.
pub fn validate_inputs(
    pa: &PresenceAbsence,
    pan_ref: &HashMap<String, Vec<Vec<u8>>>,
    genome_files: &HashMap<String, PathBuf>,
) -> Result<()> {
    use log::warn;

    // ── Security: validate genome names for path traversal ──
    if pa.n_genomes() == 0 {
        anyhow::bail!("Presence/absence CSV has no genome columns");
    }
    if pa.n_clusters() == 0 {
        anyhow::bail!("Presence/absence CSV has no cluster rows");
    }
    if pan_ref.is_empty() {
        anyhow::bail!("pan_genome_reference FASTA has no usable sequences");
    }

    let mut seen_genomes: HashSet<&str> = HashSet::new();
    for gi_val in 0..pa.n_genomes() {
        let gi = GenomeIdx(gi_val);
        let name = pa.genome_name(gi);
        if !is_safe_name(name) {
            anyhow::bail!(
                "Unsafe genome name in CSV column header: {:?} (contains path separators or traversal patterns)",
                name
            );
        }
        if !seen_genomes.insert(name) {
            anyhow::bail!(
                "Duplicate genome column in presence/absence CSV: {:?}",
                name
            );
        }
    }

    // ── Security: validate cluster names ──
    let mut seen_clusters: HashSet<&str> = HashSet::new();
    for ci_val in 0..pa.n_clusters() {
        let ci = ClusterIdx(ci_val);
        let name = pa.cluster_name(ci);
        if !is_safe_name(name) {
            anyhow::bail!(
                "Unsafe cluster name in CSV: {:?} (contains path separators or traversal patterns)",
                name
            );
        }
        if !seen_clusters.insert(name) {
            anyhow::bail!("Duplicate cluster name in presence/absence CSV: {:?}", name);
        }
    }

    // Clusters in CSV but not in pan-reference
    let mut csv_only = 0usize;
    for ci_val in 0..pa.n_clusters() {
        let ci = ClusterIdx(ci_val);
        if !pan_ref.contains_key(pa.cluster_name(ci)) {
            csv_only += 1;
        }
    }
    if csv_only == pa.n_clusters() {
        // Nothing would be searched, and the run used to finish with exit 0 and
        // no rescues, which reads exactly like a clean dataset.
        anyhow::bail!(
            "No cluster in the presence/absence CSV has a sequence in the pan-genome \
             reference. Recore matches a FASTA header's first token (Panaroo: \
             '>cluster'), or its second token (Roary: '>seqid cluster'), against the \
             CSV's Gene column."
        );
    }
    if csv_only > 0 {
        warn!(
            "{} clusters in presence/absence CSV have no sequence in pan_genome_reference.fa (will be skipped)",
            csv_only
        );
    }

    // Clusters in pan-reference but not in CSV
    let cluster_idx = pa.cluster_index();
    let ref_only: usize = pan_ref
        .keys()
        .filter(|k| !cluster_idx.contains_key(k.as_str()))
        .count();
    if ref_only > 0 {
        warn!(
            "{} clusters in pan_genome_reference.fa are not in presence/absence CSV (will be ignored)",
            ref_only
        );
    }

    // Genomes in CSV with no FASTA
    let mut no_fasta = 0usize;
    for gi_val in 0..pa.n_genomes() {
        let gi = GenomeIdx(gi_val);
        if resolve_genome_file(genome_files, pa.genome_name(gi)).is_none() {
            no_fasta += 1;
        }
    }
    if no_fasta > 0 {
        warn!(
            "{} genomes in CSV have no matching FASTA file in the genomes directory",
            no_fasta
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_temp(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f
    }

    /// Write a gzip file made of one member per element, which is the shape
    /// bgzip writes: one member per block, hundreds of them in a real genome.
    /// Two members is the smallest fixture that can tell the difference, and it
    /// has to be built by hand because a single `GzEncoder` never produces one.
    fn write_multi_member_gzip(path: &Path, members: &[&str]) {
        use flate2::write::GzEncoder;
        use flate2::Compression;

        let mut file = std::fs::File::create(path).unwrap();
        for member in members {
            let mut gz = GzEncoder::new(Vec::new(), Compression::fast());
            gz.write_all(member.as_bytes()).unwrap();
            file.write_all(&gz.finish().unwrap()).unwrap();
        }
    }

    #[test]
    fn test_parse_pa() {
        let csv = write_temp(
            "Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,\
            Avg sequences per isolate,Genome Fragment,Order within Fragment,\
            Accessory Order,QC,Min group size nuc,Max group size nuc,\
            Avg group size nuc,g1,g2,g3\n\
            gA,gA,prod,3,3,1.0,0,0,0,,100,100,100,id1,id2,id3\n\
            gB,gB,prod,2,2,1.0,0,0,0,,200,200,200,id4,,id5\n",
        );
        let pa = parse_presence_absence(csv.path()).unwrap();
        assert_eq!(pa.n_clusters(), 2);
        assert_eq!(pa.n_genomes(), 3);
        assert_eq!(pa.core_size(), 1);
        assert_eq!(pa.gaps().len(), 1);
    }

    #[test]
    fn test_parse_pan_ref_single() {
        let fa =
            write_temp(">clA~~~g1 prod\nATGCATGC\n>clA~~~g2 prod\nATGCATGCATGC\n>clB\nGGGGAAAA\n");
        let clusters: HashSet<&str> = ["clA", "clB"].into_iter().collect();
        let refs = parse_pan_reference_for(fa.path(), false, &clusters).unwrap();
        assert_eq!(refs.len(), 2);
        assert_eq!(refs["clA"].len(), 1);
        assert_eq!(refs["clA"][0].len(), 12);
        assert_eq!(refs["clB"][0], b"GGGGAAAA");
    }

    #[test]
    fn test_parse_pan_ref_multi() {
        let fa =
            write_temp(">clA~~~g1 prod\nATGCATGC\n>clA~~~g2 prod\nATGCATGCATGC\n>clB\nGGGGAAAA\n");
        let clusters: HashSet<&str> = ["clA", "clB"].into_iter().collect();
        let refs = parse_pan_reference_for(fa.path(), true, &clusters).unwrap();
        assert_eq!(refs["clA"].len(), 2);
    }

    /// Roary writes '>seqid cluster'. Keyed by the first token, a Roary run
    /// matched no cluster and exited 0 with no rescues.
    #[test]
    fn test_parse_pan_ref_roary_layout() {
        let fa = write_temp(">L1_00001 clA\nATGCATGC\n>L1_00002 clB\nGGGG\n");
        let clusters: HashSet<&str> = ["clA", "clB"].into_iter().collect();
        let refs = parse_pan_reference_for(fa.path(), false, &clusters).unwrap();
        assert_eq!(refs["clA"][0], b"ATGCATGC");
        assert_eq!(refs["clB"][0], b"GGGG");
    }

    /// Panaroo writes '~~~' inside real cluster names. Cutting every name at
    /// '~~~' lost 'pstB~~~phoT', and filed its sequence as an allele of 'pstB',
    /// replacing pstB's own as the longest.
    #[test]
    fn test_parse_pan_ref_keeps_tilde_cluster_names() {
        let fa = write_temp(">pstB\nAAAA\n>pstB~~~phoT\nCCCCCCCC\n");
        let clusters: HashSet<&str> = ["pstB", "pstB~~~phoT"].into_iter().collect();
        let refs = parse_pan_reference_for(fa.path(), false, &clusters).unwrap();
        assert_eq!(refs["pstB"][0], b"AAAA");
        assert_eq!(refs["pstB~~~phoT"][0], b"CCCCCCCC");
    }

    #[test]
    fn test_parse_pan_ref_tolerates_non_utf8_descriptions() {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(b">clA caf\xe9 product\nACGT\n").unwrap();
        let clusters: HashSet<&str> = ["clA"].into_iter().collect();
        let refs = parse_pan_reference_for(f.path(), false, &clusters).unwrap();
        assert_eq!(refs["clA"][0], b"ACGT");
    }

    /// A column with no files of its own used to borrow another column's
    /// assembly through a shared AS identifier, and rescues came from the
    /// wrong genome.
    #[test]
    fn test_resolve_genome_files_does_not_borrow_another_columns_assembly() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("L4_AS00401.fasta"), ">c\nACGT\n").unwrap();
        std::fs::write(dir.path().join("L1_AS00777.fasta"), ">c\nACGT\n").unwrap();
        let pa = PresenceAbsence::new_for_test(
            vec!["c".into()],
            vec![
                "L4_AS00401".into(),
                "L4_AS00401_v2".into(),
                "AS00777".into(),
            ],
            vec![vec!["x".into(), "".into(), "".into()]],
        );
        let files = find_genome_files(dir.path()).unwrap();
        let resolved = resolve_genome_files_for_pa(&pa, &files);
        assert!(resolved.contains_key("L4_AS00401"));
        assert!(
            !resolved.contains_key("L4_AS00401_v2"),
            "AS00401 is another column's genome"
        );
        // A unique identifier with no exact owner still resolves.
        assert!(resolved.contains_key("AS00777"));
    }

    #[test]
    fn test_gaps() {
        let csv = write_temp(
            "Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,\
            Avg sequences per isolate,Genome Fragment,Order within Fragment,\
            Accessory Order,QC,Min group size nuc,Max group size nuc,\
            Avg group size nuc,g1,g2\n\
            gA,gA,p,2,2,1,0,0,0,,1,1,1,x,y\n\
            gB,gB,p,1,1,1,0,0,0,,1,1,1,x,\n\
            gC,gC,p,0,0,0,0,0,0,,1,1,1,,\n",
        );
        let pa = parse_presence_absence(csv.path()).unwrap();
        assert_eq!(pa.gaps().len(), 3);
    }

    #[test]
    fn test_prevalence() {
        let csv = write_temp(
            "Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,\
            Avg sequences per isolate,Genome Fragment,Order within Fragment,\
            Accessory Order,QC,Min group size nuc,Max group size nuc,\
            Avg group size nuc,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10\n\
            gA,gA,p,10,10,1,0,0,0,,1,1,1,a,a,a,a,a,a,a,a,a,a\n\
            gB,gB,p,1,1,1,0,0,0,,1,1,1,b,,,,,,,,,\n",
        );
        let pa = parse_presence_absence(csv.path()).unwrap();
        assert!((pa.prevalence(ClusterIdx(0)) - 1.0).abs() < 0.001);
        assert!((pa.prevalence(ClusterIdx(1)) - 0.1).abs() < 0.001);
    }

    #[test]
    fn test_safe_name() {
        assert!(is_safe_name("AS00001"));
        assert!(is_safe_name("genome_123"));
        assert!(is_safe_name("my-genome.v2"));
        assert!(!is_safe_name(""));
        assert!(!is_safe_name("."));
        assert!(!is_safe_name(".."));
        assert!(!is_safe_name("../etc/passwd"));
        assert!(!is_safe_name("foo/bar"));
        assert!(!is_safe_name("foo\\bar"));
        assert!(!is_safe_name("foo\0bar"));
        assert!(!is_safe_name("..\\Windows\\System32"));
    }

    #[test]
    fn test_extract_as_id() {
        assert_eq!(extract_as_id("AS00146"), Some("AS00146"));
        assert_eq!(extract_as_id("L4_AS00146"), Some("AS00146"));
        assert_eq!(extract_as_id("prefix_AS00146_suffix"), Some("AS00146"));
        assert_eq!(extract_as_id("AS001461"), None);
        assert_eq!(extract_as_id("L4_00146"), None);
    }

    #[test]
    fn test_resolve_genome_file_by_as_suffix() {
        let mut files = HashMap::new();
        files.insert("AS00146".to_string(), PathBuf::from("AS00146.fna"));
        let resolved = resolve_genome_file(&files, "L4_AS00146").unwrap();
        assert_eq!(resolved, &PathBuf::from("AS00146.fna"));
    }

    #[test]
    fn test_find_genome_files_indexes_embedded_as_id() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("L4_AS00146.fna");
        std::fs::write(&path, b">ctg\nATGC\n").unwrap();

        let files = find_genome_files(dir.path()).unwrap();
        assert_eq!(files.get("L4_AS00146"), Some(&path));
        assert_eq!(files.get("AS00146"), Some(&path));
    }

    #[test]
    fn test_find_genome_files_keeps_ambiguous_as_id_unaliased() {
        let dir = tempfile::tempdir().unwrap();
        let path_a = dir.path().join("L4_AS00146.fna");
        let path_b = dir.path().join("other_AS00146.fna");
        std::fs::write(&path_a, b">ctg\nATGC\n").unwrap();
        std::fs::write(&path_b, b">ctg\nATGC\n").unwrap();

        let files = find_genome_files(dir.path()).unwrap();
        assert_eq!(files.get("L4_AS00146"), Some(&path_a));
        assert_eq!(files.get("other_AS00146"), Some(&path_b));
        assert_eq!(files.get("AS00146"), None);
    }

    #[test]
    fn test_find_annotation_file_by_embedded_as_id() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("L4_AS00146.gff3");
        std::fs::write(&path, b"##gff-version 3\n").unwrap();

        assert_eq!(find_annotation_file(dir.path(), "AS00146"), Some(path));
    }

    #[test]
    fn test_find_annotation_file_ignores_ambiguous_embedded_as_id() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("L4_AS00146.gff3"), b"##gff-version 3\n").unwrap();
        std::fs::write(dir.path().join("other_AS00146.gff3"), b"##gff-version 3\n").unwrap();

        assert_eq!(find_annotation_file(dir.path(), "AS00146"), None);
    }

    #[test]
    fn test_load_genome_contig_names() {
        let fa = write_temp(">AS00146 extra words\nATGC\n>plasmid\nGGGG\n");
        let contigs = load_genome_contig_names(fa.path()).unwrap();
        assert_eq!(contigs, vec!["AS00146", "plasmid"]);
    }

    #[test]
    fn test_normalize_single_contig_gff_to_fasta() {
        let mut records = vec![
            (
                "contig_1".to_string(),
                0,
                10,
                '+',
                "gene1".to_string(),
                vec!["gene1".to_string()],
                vec!["gene1".to_string()],
            ),
            (
                "contig_1".to_string(),
                20,
                30,
                '+',
                "gene2".to_string(),
                vec!["gene2".to_string()],
                vec!["gene2".to_string()],
            ),
        ];
        let remap = normalize_gff_records_to_fasta(&mut records, &[String::from("AS00146")]);
        assert_eq!(
            remap,
            Some(("contig_1".to_string(), "AS00146".to_string(), 2))
        );
        assert!(records.iter().all(|record| record.0 == "AS00146"));
    }

    #[test]
    fn test_parse_gff_genes_prefers_locus_tag_and_keeps_aliases() {
        let gff = write_temp(
            "##gff-version 3\n\
             chr\tRefSeq\tCDS\t10\t99\t.\t+\t0\tID=cds-Rv0001;locus_tag=Rv0001;gene=dnaA;Name=chromosomal_replication_initiator;product=hypothetical%20protein\n",
        );
        let genes = parse_gff_genes(gff.path()).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].4, "Rv0001");
        assert!(genes[0].5.iter().any(|alias| alias == "cds-Rv0001"));
        assert!(genes[0].5.iter().any(|alias| alias == "dnaA"));
        assert!(!genes[0]
            .5
            .iter()
            .any(|alias| alias == "hypothetical%20protein"));
        assert!(genes[0]
            .6
            .iter()
            .any(|label| label == "hypothetical protein"));
    }

    #[test]
    fn test_parse_gff_decodes_percent_encoded_identifiers() {
        let gff = write_temp(
            "##gff-version 3\n\
             chr\tRefSeq\tCDS\t10\t99\t.\t+\t0\tID=cds-Rv%3A0001;locus_tag=Rv%3A0001;product=foo%20bar\n",
        );
        let genes = parse_gff_genes(gff.path()).unwrap();
        assert_eq!(genes[0].4, "Rv:0001");
        assert!(genes[0].5.iter().any(|alias| alias == "Rv:0001"));
        assert!(genes[0].6.iter().any(|label| label == "foo bar"));
    }

    #[test]
    fn test_parse_gff_genes_skips_invalid_coordinates() {
        let gff = write_temp(
            "##gff-version 3\n\
             chr\tRefSeq\tCDS\tbad\t99\t.\t+\t0\tID=bad_start;locus_tag=bad_start\n\
             chr\tRefSeq\tCDS\t50\t10\t.\t+\t0\tID=inverted;locus_tag=inverted\n\
             chr\tRefSeq\tCDS\t10\t99\t.\t+\t0\tID=good;locus_tag=good\n",
        );
        let genes = parse_gff_genes(gff.path()).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].4, "good");
        assert_eq!(genes[0].1, 9);
        assert_eq!(genes[0].2, 99);
    }

    #[test]
    fn test_validate_inputs_rejects_traversal() {
        let pa = PresenceAbsence::new_for_test(
            vec!["cluster1".into()],
            vec!["../etc/passwd".into()],
            vec![vec!["gene1".into()]],
        );
        let pan_ref: HashMap<String, Vec<Vec<u8>>> = {
            let mut m = HashMap::new();
            m.insert("cluster1".into(), vec![b"ATGC".to_vec()]);
            m
        };
        let genome_files: HashMap<String, PathBuf> = HashMap::new();
        let result = validate_inputs(&pa, &pan_ref, &genome_files);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Unsafe genome name"));
    }

    #[test]
    fn test_validate_inputs_rejects_duplicate_names() {
        let pan_ref: HashMap<String, Vec<Vec<u8>>> = {
            let mut m = HashMap::new();
            m.insert("cluster1".into(), vec![b"ATGC".to_vec()]);
            m.insert("cluster2".into(), vec![b"ATGC".to_vec()]);
            m
        };
        let genome_files: HashMap<String, PathBuf> = HashMap::new();

        let pa = PresenceAbsence::new_for_test(
            vec!["cluster1".into()],
            vec!["g1".into(), "g1".into()],
            vec![vec!["gene1".into(), "".into()]],
        );
        let err = validate_inputs(&pa, &pan_ref, &genome_files).unwrap_err();
        assert!(err.to_string().contains("Duplicate genome column"));

        let pa = PresenceAbsence::new_for_test(
            vec!["cluster1".into(), "cluster1".into()],
            vec!["g1".into()],
            vec![vec!["gene1".into()], vec!["gene2".into()]],
        );
        let err = validate_inputs(&pa, &pan_ref, &genome_files).unwrap_err();
        assert!(err.to_string().contains("Duplicate cluster name"));
    }

    /// One Latin-1 byte in a header's description used to fail this read, and
    /// the stop scan skipped the genome without a word: its pseudogenes went
    /// into the matrix as intact genes.
    #[test]
    fn test_load_genome_contigs_tolerates_non_utf8_descriptions() {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(b">ctg1 isolated in S\xe3o Paulo\nATCG\nATCG\n>ctg2\nGG\n")
            .unwrap();
        let contigs = load_genome_contigs(f.path()).unwrap();
        assert_eq!(contigs["ctg1"], b"ATCGATCG");
        assert_eq!(contigs["ctg2"], b"GG");
        assert_eq!(
            load_genome_contig_names(f.path()).unwrap(),
            vec!["ctg1", "ctg2"]
        );
    }

    /// A Latin-1 byte in a free-text attribute used to lose the genome's whole
    /// annotation.
    #[test]
    fn test_parse_gff_genes_tolerates_non_utf8_attributes() {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(
            b"##gff-version 3\nchr\tsrc\tCDS\t10\t99\t.\t+\t0\tID=a;locus_tag=a;product=subunit \xdf\n\
              chr\tsrc\tCDS\t200\t299\t.\t+\t0\tID=b;locus_tag=b\n",
        )
        .unwrap();
        let genes = parse_gff_genes(f.path()).unwrap();
        let ids: Vec<&str> = genes.iter().map(|gene| gene.4.as_str()).collect();
        assert_eq!(ids, vec!["a", "b"]);
    }

    #[test]
    fn test_load_genome_contigs_plain() {
        let fa = write_temp(">contig1\nATCGATCG\n>contig2\nGGGGAAAA\n");
        let contigs = load_genome_contigs(fa.path()).unwrap();
        assert_eq!(contigs.len(), 2);
        assert_eq!(contigs["contig1"], b"ATCGATCG");
        assert_eq!(contigs["contig2"], b"GGGGAAAA");
    }

    /// The regression test for the multi-member case. Its single-member twin
    /// below passed throughout, which is the whole reason the defect survived:
    /// the only compressed fixture in the project was the one shape that
    /// `GzDecoder` reads correctly.
    #[test]
    fn test_load_genome_contigs_multi_member_gz() {
        let dir = tempfile::tempdir().unwrap();
        let gz_path = dir.path().join("genome.fasta.gz");
        write_multi_member_gzip(&gz_path, &[">ctg1\nATCGATCG\n", ">ctg2\nTTTTCCCC\n"]);

        let contigs = load_genome_contigs(&gz_path).unwrap();
        assert_eq!(contigs.len(), 2);
        assert_eq!(contigs["ctg1"], b"ATCGATCG");
        assert_eq!(contigs["ctg2"], b"TTTTCCCC");
    }

    /// Pinned at the source as well as at the callers, so that rerouting a
    /// caller cannot quietly undo it.
    #[test]
    fn test_open_maybe_gzip_reads_every_member() {
        use std::io::BufRead;

        let dir = tempfile::tempdir().unwrap();
        let gz_path = dir.path().join("genome.fasta.gz");
        write_multi_member_gzip(&gz_path, &[">ctg1\nATCGATCG\n", ">ctg2\nTTTTCCCC\n"]);

        let reader = open_maybe_gzip(&gz_path).unwrap();
        let lines: Vec<String> = reader.lines().map(|line| line.unwrap()).collect();

        assert_eq!(lines, vec![">ctg1", "ATCGATCG", ">ctg2", "TTTTCCCC"]);
    }

    /// Annotations take the same path, and a truncated GFF is worse than a
    /// short one: the run kept its rescues and silently reclassified them.
    #[test]
    fn test_parse_gff_genes_multi_member_gz() {
        let dir = tempfile::tempdir().unwrap();
        let gz_path = dir.path().join("annotation.gff3.gz");
        write_multi_member_gzip(
            &gz_path,
            &[
                "##gff-version 3\nchr\tRefSeq\tCDS\t10\t99\t.\t+\t0\tID=first;locus_tag=first\n",
                "chr\tRefSeq\tCDS\t200\t299\t.\t+\t0\tID=second;locus_tag=second\n",
            ],
        );

        let genes = parse_gff_genes(&gz_path).unwrap();
        let ids: Vec<&str> = genes.iter().map(|gene| gene.4.as_str()).collect();
        assert_eq!(ids, vec!["first", "second"]);
    }

    #[test]
    fn test_load_genome_contigs_gz() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;

        let dir = tempfile::tempdir().unwrap();
        let gz_path = dir.path().join("genome.fasta.gz");
        {
            let file = std::fs::File::create(&gz_path).unwrap();
            let mut gz = GzEncoder::new(file, Compression::fast());
            gz.write_all(b">ctg1\nATCGATCG\n>ctg2\nTTTTCCCC\n").unwrap();
            gz.finish().unwrap();
        }
        let contigs = load_genome_contigs(&gz_path).unwrap();
        assert_eq!(contigs.len(), 2);
        assert_eq!(contigs["ctg1"], b"ATCGATCG");
        assert_eq!(contigs["ctg2"], b"TTTTCCCC");
    }
}
