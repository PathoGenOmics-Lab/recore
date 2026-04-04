use ahash::{HashMap, HashMapExt, HashSet, HashSetExt};
use anyhow::{Context, Result};
use log::info;
use std::fs;
use std::path::{Path, PathBuf};

/// Presence/absence matrix: gene_cluster -> {genome -> gene_id or empty}.
/// We store which genomes have gaps (gene absent).
pub struct PresenceAbsence {
    /// All gene cluster names, in order.
    pub clusters: Vec<String>,
    /// All genome names (column headers), in order.
    pub genomes: Vec<String>,
    /// For each cluster index, set of genome indices where the gene is present.
    pub present: Vec<HashSet<usize>>,
    /// Original cell values: present[cluster_idx][genome_idx] = gene_id string.
    pub cell_values: Vec<Vec<String>>,
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

    /// Number of genomes where a cluster is present.
    #[allow(dead_code)]
    pub fn cluster_count(&self, ci: usize) -> usize {
        self.present[ci].len()
    }

    /// Core size: clusters present in ALL genomes.
    pub fn core_size(&self) -> usize {
        let n = self.genomes.len();
        self.present.iter().filter(|p| p.len() == n).count()
    }
}

/// Panaroo metadata column names to skip when identifying genome columns.
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

    // Genome columns: indices and names
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

    let mut clusters = Vec::new();
    let mut present = Vec::new();
    let mut cell_values = Vec::new();

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

        clusters.push(cluster_name);
        present.push(present_set);
        cell_values.push(row_values);
    }

    info!(
        "Parsed presence/absence matrix: {} clusters × {} genomes",
        clusters.len(),
        genomes.len()
    );

    Ok(PresenceAbsence {
        clusters,
        genomes,
        present,
        cell_values,
    })
}

/// Parse pan_genome_reference.fa → HashMap<cluster_name, sequence>.
/// Panaroo format: >cluster_name~~~gene_id description
pub fn parse_pan_reference(path: &Path) -> Result<HashMap<String, Vec<u8>>> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("Cannot read {}", path.display()))?;

    let mut refs: HashMap<String, Vec<u8>> = HashMap::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = Vec::new();

    for line in content.lines() {
        if let Some(header) = line.strip_prefix('>') {
            // Save previous record
            if let Some(name) = current_name.take() {
                let existing = refs.get(&name);
                if existing.is_none() || existing.unwrap().len() < current_seq.len() {
                    refs.insert(name, current_seq.clone());
                }
            }
            current_seq.clear();

            // Parse name: cluster_name~~~gene_id or just cluster_name
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
        let existing = refs.get(&name);
        if existing.is_none() || existing.unwrap().len() < current_seq.len() {
            refs.insert(name, current_seq);
        }
    }

    info!(
        "Parsed {} representative sequences from {}",
        refs.len(),
        path.display()
    );
    Ok(refs)
}

/// Find genome FASTA files in a directory.
/// Returns HashMap<genome_name, file_path>.
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

        // Strip extension to get genome name
        let mut name = filename.clone();
        for ext in &extensions {
            if name.ends_with(ext) {
                name = name[..name.len() - ext.len()].to_string();
                break;
            }
        }

        // Only accept known extensions
        if name != filename {
            files.entry(name).or_insert(path);
        }
    }

    info!("Found {} genome files in {}", files.len(), dir.display());
    Ok(files)
}

/// Parse CDS features from a GFF file.
/// Returns Vec<(contig, start, end, gene_id)>.
pub fn parse_gff_genes(path: &Path) -> Result<Vec<(String, u64, u64, String)>> {
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

        // Extract ID or locus_tag from attributes
        let mut gene_id = String::new();
        for attr in fields[8].split(';') {
            if let Some(val) = attr.strip_prefix("ID=") {
                gene_id = val.to_string();
                break;
            } else if let Some(val) = attr.strip_prefix("locus_tag=") {
                gene_id = val.to_string();
            }
        }

        genes.push((contig, start, end, gene_id));
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
        assert_eq!(pa.core_size(), 1); // gA is core, gB is not
        assert_eq!(pa.gaps().len(), 1); // gB missing in g2
    }

    #[test]
    fn test_parse_pan_ref() {
        let fa = write_temp(">clusterA~~~gene1 product\nATGCATGC\n>clusterB\nGGGGAAAA\n");
        let refs = parse_pan_reference(fa.path()).unwrap();
        assert_eq!(refs.len(), 2);
        assert_eq!(refs["clusterA"], b"ATGCATGC");
        assert_eq!(refs["clusterB"], b"GGGGAAAA");
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
        let gaps = pa.gaps();
        // gB missing in g2 (1 gap), gC missing in g1 and g2 (2 gaps)
        assert_eq!(gaps.len(), 3);
    }
}
