use ahash::{HashMap, HashMapExt, HashSet, HashSetExt};
use anyhow::{Context, Result};
use log::info;
use std::fs;
use std::path::{Path, PathBuf};

/// Presence/absence matrix.
pub struct PresenceAbsence {
    /// All gene cluster names, in order.
    pub clusters: Vec<String>,
    /// All genome names (column headers), in order.
    pub genomes: Vec<String>,
    /// For each cluster index, set of genome indices where the gene is present.
    pub present: Vec<HashSet<usize>>,
    /// Original cell values: cell_values[cluster_idx][genome_idx] = gene_id string.
    pub cell_values: Vec<Vec<String>>,
    /// Full CSV header (all columns, including metadata).
    pub header: Vec<String>,
    /// Metadata columns per cluster: metadata[cluster_idx] = vec of metadata values.
    /// These are the non-genome columns (Annotation, No. isolates, etc.) preserved
    /// for faithful output.
    pub metadata: Vec<Vec<String>>,
    /// Indices of genome columns in the original header.
    pub genome_col_indices: Vec<usize>,
}

impl PresenceAbsence {
    /// Find all (cluster_idx, genome_idx) pairs where the gene is absent.
    pub fn gaps(&self) -> Vec<(usize, usize)> {
        let mut gaps = Vec::new();
        for (ci, present_set) in self.present.iter().enumerate() {
            for gi in 0..self.genomes.len() {
                if !present_set.contains(&gi) {
                    gaps.push((ci, gi));
                }
            }
        }
        gaps
    }

    /// Cluster prevalence: fraction of genomes where cluster is present.
    pub fn prevalence(&self, ci: usize) -> f64 {
        if self.genomes.is_empty() {
            return 0.0;
        }
        self.present[ci].len() as f64 / self.genomes.len() as f64
    }

    /// Core size: clusters present in ALL genomes.
    pub fn core_size(&self) -> usize {
        let n = self.genomes.len();
        self.present.iter().filter(|p| p.len() == n).count()
    }

    /// Soft-core size: clusters present in ≥95% of genomes.
    pub fn soft_core_size(&self) -> usize {
        let threshold = (self.genomes.len() as f64 * 0.95).ceil() as usize;
        self.present.iter().filter(|p| p.len() >= threshold).count()
    }
}

/// Panaroo metadata columns to skip when identifying genome columns.
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

/// Parse Panaroo/Roary gene_presence_absence.csv.
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

    // Metadata column indices (everything that's not a genome column)
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
                present_set.insert(gi);
            }
            row_values.push(val);
        }

        // Store metadata columns
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

/// Parse pan_genome_reference.fa → HashMap<cluster_name, Vec<sequence>>.
/// In multi-allele mode, keeps ALL alleles per cluster.
/// In single mode, keeps only the longest.
pub fn parse_pan_reference(
    path: &Path,
    multi_allele: bool,
) -> Result<HashMap<String, Vec<Vec<u8>>>> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("Cannot read {}", path.display()))?;

    let mut refs: HashMap<String, Vec<Vec<u8>>> = HashMap::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = Vec::new();

    for line in content.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some(name) = current_name.take() {
                if !current_seq.is_empty() {
                    let entry = refs.entry(name).or_default();
                    if multi_allele {
                        entry.push(current_seq.clone());
                    } else if entry.is_empty() || entry[0].len() < current_seq.len() {
                        // Single mode: keep only the longest allele
                        if entry.is_empty() {
                            entry.push(current_seq.clone());
                        } else {
                            entry[0] = current_seq.clone();
                        }
                    }
                }
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

    // Last record
    if let Some(name) = current_name {
        if !current_seq.is_empty() {
            let entry = refs.entry(name).or_default();
            if multi_allele {
                entry.push(current_seq);
            } else if entry.is_empty() || entry[0].len() < current_seq.len() {
                if entry.is_empty() {
                    entry.push(current_seq);
                } else {
                    entry[0] = current_seq;
                }
            }
        }
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

/// Find genome FASTA files in a directory.
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

/// (contig, start, end, strand, gene_id)
pub type GffRecord = (String, u64, u64, char, String);

/// Parse CDS features from a GFF file.
pub fn parse_gff_genes(path: &Path) -> Result<Vec<GffRecord>> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("Cannot read GFF {}", path.display()))?;

    let mut genes = Vec::new();
    for line in content.lines() {
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
        assert_eq!(pa.clusters.len(), 2);
        assert_eq!(pa.genomes.len(), 3);
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
        assert_eq!(refs["clA"][0].len(), 12); // kept the longest
        assert_eq!(refs["clB"][0], b"GGGGAAAA");
    }

    #[test]
    fn test_parse_pan_ref_multi() {
        let fa = write_temp(
            ">clA~~~g1 prod\nATGCATGC\n>clA~~~g2 prod\nATGCATGCATGC\n>clB\nGGGGAAAA\n",
        );
        let refs = parse_pan_reference(fa.path(), true).unwrap();
        assert_eq!(refs["clA"].len(), 2); // both alleles
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
        assert!((pa.prevalence(0) - 1.0).abs() < 0.001); // gA: 100%
        assert!((pa.prevalence(1) - 0.1).abs() < 0.001); // gB: 10%
    }
}
