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
            clusters, genomes, present, cell_values, header, metadata, genome_col_indices,
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
    let mut reader = csv::Reader::from_path(path)
        .with_context(|| format!("Cannot open {}", path.display()))?;

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

    let genomes: Vec<String> = genome_indices
        .iter()
        .map(|&i| headers[i].clone())
        .collect();

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

    let file = fs::File::open(path)
        .with_context(|| format!("Cannot open {}", path.display()))?;
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
                flush_record(name, std::mem::take(&mut current_seq), &mut refs, multi_allele);
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

    for entry in fs::read_dir(dir)
        .with_context(|| format!("Cannot read directory {}", dir.display()))?
    {
        let entry = entry?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }

        let filename = path.file_name().unwrap().to_string_lossy().to_string();

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

pub type GffRecord = (String, u64, u64, char, String);

/// Open a file with transparent gzip decompression.
/// If the path ends with `.gz`, decompresses on-the-fly via flate2.
fn open_maybe_gzip(path: &Path) -> Result<Box<dyn std::io::BufRead>> {
    use std::io::BufReader;
    let file = fs::File::open(path)
        .with_context(|| format!("Cannot open {}", path.display()))?;
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
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 || fields[2] != "CDS" {
            continue;
        }

        let contig = fields[0].to_string();
        let start: u64 = fields[3].parse().unwrap_or(0);
        let end: u64 = fields[4].parse().unwrap_or(0);
        let strand = fields[6].chars().next().unwrap_or('+');

        let mut gene_id = String::new();
        for attr in fields[8].split(';') {
            let attr = attr.trim();
            if let Some(val) = attr.strip_prefix("ID=") {
                gene_id = val.to_string();
                break;
            } else if let Some(val) = attr.strip_prefix("locus_tag=") {
                gene_id = val.to_string();
            }
        }

        genes.push((contig, start, end, strand, gene_id));
    }

    Ok(genes)
}

/// Cross-validate inputs: warn about mismatches between CSV clusters and pan-reference.
pub fn validate_inputs(
    pa: &PresenceAbsence,
    pan_ref: &HashMap<String, Vec<Vec<u8>>>,
    genome_files: &HashMap<String, PathBuf>,
) {
    use log::warn;

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
    let ref_only: usize = pan_ref.keys()
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
        if !genome_files.contains_key(pa.genome_name(gi)) {
            no_fasta += 1;
        }
    }
    if no_fasta > 0 {
        warn!(
            "{} genomes in CSV have no matching FASTA file in the genomes directory",
            no_fasta
        );
    }
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
        let fa = write_temp(
            ">clA~~~g1 prod\nATGCATGC\n>clA~~~g2 prod\nATGCATGCATGC\n>clB\nGGGGAAAA\n",
        );
        let refs = parse_pan_reference(fa.path(), false).unwrap();
        assert_eq!(refs.len(), 2);
        assert_eq!(refs["clA"].len(), 1);
        assert_eq!(refs["clA"][0].len(), 12);
        assert_eq!(refs["clB"][0], b"GGGGAAAA");
    }

    #[test]
    fn test_parse_pan_ref_multi() {
        let fa = write_temp(
            ">clA~~~g1 prod\nATGCATGC\n>clA~~~g2 prod\nATGCATGCATGC\n>clB\nGGGGAAAA\n",
        );
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
}
