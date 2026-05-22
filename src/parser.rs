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
    present: Vec<HashSet<GenomeIdx>>,
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
        self.present[ci.val()].contains(&gi)
    }

    #[inline]
    pub fn present_count(&self, ci: ClusterIdx) -> usize {
        self.present[ci.val()].len()
    }

    pub fn present_genomes(&self, ci: ClusterIdx) -> &HashSet<GenomeIdx> {
        &self.present[ci.val()]
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
        self.present.iter().filter(|p| p.len() == n).count()
    }

    pub fn soft_core_size(&self) -> usize {
        let threshold = (self.genomes.len() as f64 * 0.95).ceil() as usize;
        self.present.iter().filter(|p| p.len() >= threshold).count()
    }

    /// Constructor for tests only — minimal PA with all genes absent.
    #[cfg(test)]
    pub fn new_for_test(
        clusters: Vec<String>,
        genomes: Vec<String>,
        cell_values: Vec<Vec<String>>,
    ) -> Self {
        use ahash::HashSetExt;
        let n_cl = clusters.len();
        let n_ge = genomes.len();
        let present = (0..n_cl)
            .map(|ci| {
                let mut s = HashSet::new();
                for (gi, v) in cell_values[ci].iter().enumerate() {
                    if !v.is_empty() {
                        s.insert(GenomeIdx(gi));
                    }
                }
                s
            })
            .collect();
        let metadata = vec![vec![]; n_cl];
        let genome_col_indices: Vec<usize> = (0..n_ge).collect();
        let mut header = vec!["Gene".to_string()];
        header.extend(genomes.iter().cloned());
        PresenceAbsence {
            clusters,
            genomes,
            present,
            cell_values,
            header,
            metadata,
            genome_col_indices,
        }
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
    let mut present = Vec::new();
    let mut cell_values = Vec::new();
    let mut metadata = Vec::new();

    for result in reader.records() {
        let record = result?;
        let cluster_name = record.get(gene_col).unwrap_or("").to_string();

        let mut present_set = HashSet::new();
        let mut row_values = Vec::with_capacity(genomes.len());

        for (gi, &col_idx) in genome_indices.iter().enumerate() {
            let val = record.get(col_idx).unwrap_or("").trim().to_string();
            if !val.is_empty() {
                present_set.insert(GenomeIdx(gi));
            }
            row_values.push(val);
        }

        let meta_values: Vec<String> = meta_col_indices
            .iter()
            .map(|&i| record.get(i).unwrap_or("").to_string())
            .collect();

        clusters.push(cluster_name);
        present.push(present_set);
        cell_values.push(row_values);
        metadata.push(meta_values);
    }

    info!(
        "Parsed presence/absence: {} clusters × {} genomes (core: {}, soft-core 95%: {})",
        clusters.len(),
        genomes.len(),
        {
            let n = genomes.len();
            present.iter().filter(|p| p.len() == n).count()
        },
        {
            let t = (genomes.len() as f64 * 0.95).ceil() as usize;
            present.iter().filter(|p| p.len() >= t).count()
        }
    );

    Ok(PresenceAbsence {
        clusters,
        genomes,
        present,
        cell_values,
        header: headers,
        metadata,
        genome_col_indices: genome_indices,
    })
}

// ── FASTA / GFF parsing ──

pub fn parse_pan_reference(
    path: &Path,
    multi_allele: bool,
) -> Result<HashMap<String, Vec<Vec<u8>>>> {
    use std::io::{BufRead, BufReader};

    let file = fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut refs: HashMap<String, Vec<Vec<u8>>> = HashMap::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();

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

    for line_result in reader.lines() {
        let line = line_result?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(name) = current_name.take() {
                flush_record(
                    name,
                    std::mem::take(&mut current_seq),
                    &mut refs,
                    multi_allele,
                );
            }
            current_seq.clear();

            let id = header.split_whitespace().next().unwrap_or(header);
            let name = if let Some(pos) = id.find("~~~") {
                &id[..pos]
            } else {
                id
            };
            current_name = Some(name.to_string());
        } else {
            current_seq.extend(line.trim().as_bytes());
        }
    }

    if let Some(name) = current_name {
        flush_record(name, current_seq, &mut refs, multi_allele);
    }

    let total_alleles: usize = refs.values().map(|v| v.len()).sum();
    info!(
        "Parsed {} clusters ({} total alleles) from {}",
        refs.len(),
        total_alleles,
        path.display()
    );
    Ok(refs)
}

pub fn find_genome_files(dir: &Path) -> Result<HashMap<String, PathBuf>> {
    let extensions = [".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz"];
    let mut files: HashMap<String, PathBuf> = HashMap::new();

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
            files.entry(name).or_insert(path);
        }
    }

    info!("Found {} genome files in {}", files.len(), dir.display());
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

/// Build a map keyed by presence/absence genome names, with AS-ID fallback.
pub fn resolve_genome_files_for_pa(
    pa: &PresenceAbsence,
    genome_files: &HashMap<String, PathBuf>,
) -> HashMap<String, PathBuf> {
    let mut resolved = HashMap::new();
    for genome_name in pa.genome_names() {
        if let Some(path) = resolve_genome_file(genome_files, genome_name) {
            resolved.insert(genome_name.clone(), path.clone());
        }
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
    None
}

/// Load genome FASTA as contig_name → sequence map.
/// Used for internal stop codon scanning in pseudogene detection.
pub fn load_genome_contigs(path: &Path) -> Result<HashMap<String, Vec<u8>>> {
    use std::io::BufRead;
    let reader = open_maybe_gzip(path)?;
    let mut contigs: HashMap<String, Vec<u8>> = HashMap::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = Vec::new();

    for line_result in reader.lines() {
        let line = line_result?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(name) = current_name.take() {
                contigs.insert(name, std::mem::take(&mut current_seq));
            }
            let id = header.split_whitespace().next().unwrap_or(header);
            current_name = Some(id.to_string());
            current_seq.clear();
        } else {
            current_seq.extend(line.trim().as_bytes());
        }
    }
    if let Some(name) = current_name {
        contigs.insert(name, current_seq);
    }
    Ok(contigs)
}

/// Load only FASTA contig names, preserving file order.
pub fn load_genome_contig_names(path: &Path) -> Result<Vec<String>> {
    use std::io::BufRead;
    let reader = open_maybe_gzip(path)?;
    let mut names = Vec::new();
    for line_result in reader.lines() {
        let line = line_result?;
        if let Some(header) = line.strip_prefix('>') {
            let id = header.split_whitespace().next().unwrap_or(header);
            names.push(id.to_string());
        }
    }
    Ok(names)
}

pub type GffRecord = (String, u64, u64, char, String, Vec<String>);

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
pub(crate) fn open_maybe_gzip(path: &Path) -> Result<Box<dyn std::io::BufRead>> {
    use std::io::BufReader;
    let file = fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;
    if path.extension().is_some_and(|e| e == "gz") {
        let decoder = flate2::read::GzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

pub fn parse_gff_genes(path: &Path) -> Result<Vec<GffRecord>> {
    use std::io::BufRead;

    let reader = open_maybe_gzip(path)?;

    let mut genes = Vec::new();
    for line_result in reader.lines() {
        let line = line_result?;
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
        let start: u64 = start_field.parse::<u64>().unwrap_or(1).saturating_sub(1);
        let end: u64 = end_field.parse().unwrap_or(0); // end is already exclusive after -1+1
        let strand = strand_field.chars().next().unwrap_or('+');

        let mut id = String::new();
        let mut locus_tag = String::new();
        let mut gene_name = String::new();
        let mut name = String::new();
        let mut protein_id = String::new();
        let mut old_locus_tag = String::new();
        for attr in attributes.split(';') {
            let attr = attr.trim();
            if let Some(val) = attr.strip_prefix("ID=") {
                id = val.to_string();
            } else if let Some(val) = attr.strip_prefix("locus_tag=") {
                locus_tag = val.to_string();
            } else if let Some(val) = attr.strip_prefix("gene=") {
                gene_name = val.to_string();
            } else if let Some(val) = attr.strip_prefix("Name=") {
                name = val.to_string();
            } else if let Some(val) = attr.strip_prefix("protein_id=") {
                protein_id = val.to_string();
            } else if let Some(val) = attr.strip_prefix("old_locus_tag=") {
                old_locus_tag = val.to_string();
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

        genes.push((contig, start, end, strand, gene_id, aliases));
    }

    Ok(genes)
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
    for gi_val in 0..pa.n_genomes() {
        let gi = GenomeIdx(gi_val);
        let name = pa.genome_name(gi);
        if !is_safe_name(name) {
            anyhow::bail!(
                "Unsafe genome name in CSV column header: {:?} (contains path separators or traversal patterns)",
                name
            );
        }
    }

    // ── Security: validate cluster names ──
    for ci_val in 0..pa.n_clusters() {
        let ci = ClusterIdx(ci_val);
        let name = pa.cluster_name(ci);
        if !is_safe_name(name) {
            anyhow::bail!(
                "Unsafe cluster name in CSV: {:?} (contains path separators or traversal patterns)",
                name
            );
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
        let refs = parse_pan_reference(fa.path(), false).unwrap();
        assert_eq!(refs.len(), 2);
        assert_eq!(refs["clA"].len(), 1);
        assert_eq!(refs["clA"][0].len(), 12);
        assert_eq!(refs["clB"][0], b"GGGGAAAA");
    }

    #[test]
    fn test_parse_pan_ref_multi() {
        let fa =
            write_temp(">clA~~~g1 prod\nATGCATGC\n>clA~~~g2 prod\nATGCATGCATGC\n>clB\nGGGGAAAA\n");
        let refs = parse_pan_reference(fa.path(), true).unwrap();
        assert_eq!(refs["clA"].len(), 2);
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
            ),
            (
                "contig_1".to_string(),
                20,
                30,
                '+',
                "gene2".to_string(),
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
             chr\tRefSeq\tCDS\t10\t99\t.\t+\t0\tID=cds-Rv0001;locus_tag=Rv0001;gene=dnaA;Name=chromosomal_replication_initiator\n",
        );
        let genes = parse_gff_genes(gff.path()).unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].4, "Rv0001");
        assert!(genes[0].5.iter().any(|alias| alias == "cds-Rv0001"));
        assert!(genes[0].5.iter().any(|alias| alias == "dnaA"));
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
    fn test_load_genome_contigs_plain() {
        let fa = write_temp(">contig1\nATCGATCG\n>contig2\nGGGGAAAA\n");
        let contigs = load_genome_contigs(fa.path()).unwrap();
        assert_eq!(contigs.len(), 2);
        assert_eq!(contigs["contig1"], b"ATCGATCG");
        assert_eq!(contigs["contig2"], b"GGGGAAAA");
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
