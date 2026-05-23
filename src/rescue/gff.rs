use crate::types::ClusterIdx;
use ahash::{HashMap, HashMapExt};

#[derive(Debug, Clone)]
pub struct GffGene {
    pub contig: String,
    pub start: u64,
    pub end: u64,
    #[allow(dead_code)]
    pub strand: char,
    pub gene_id: String,
}

pub struct GffIndex {
    pub genes: Vec<GffGene>,
    contig_ranges: HashMap<String, (usize, usize)>,
    gene_id_index: HashMap<String, usize>,
    gene_aliases: HashMap<String, Vec<String>>,
    gene_labels: HashMap<String, Vec<String>>,
}

impl GffIndex {
    #[allow(dead_code)]
    pub fn build(mut genes: Vec<GffGene>) -> Self {
        Self::build_with_aliases(
            genes
                .drain(..)
                .map(|gene| (gene, Vec::new()))
                .collect::<Vec<_>>(),
        )
    }

    pub fn build_with_aliases(mut records: Vec<(GffGene, Vec<String>)>) -> Self {
        let records = records
            .drain(..)
            .map(|(gene, aliases)| (gene, aliases.clone(), aliases))
            .collect();
        Self::build_with_aliases_and_labels(records)
    }

    pub fn build_with_aliases_and_labels(
        mut records: Vec<(GffGene, Vec<String>, Vec<String>)>,
    ) -> Self {
        records.sort_by(|a, b| a.0.contig.cmp(&b.0.contig).then(a.0.start.cmp(&b.0.start)));

        let genes: Vec<GffGene> = records.iter().map(|(gene, _, _)| gene.clone()).collect();

        let mut contig_ranges = HashMap::new();
        if !genes.is_empty() {
            let mut current_contig = &genes[0].contig;
            let mut range_start = 0;
            for (i, g) in genes.iter().enumerate() {
                if &g.contig != current_contig {
                    contig_ranges.insert(current_contig.clone(), (range_start, i));
                    current_contig = &g.contig;
                    range_start = i;
                }
            }
            contig_ranges.insert(current_contig.clone(), (range_start, genes.len()));
        }

        let mut gene_id_index = HashMap::new();
        let mut gene_aliases = HashMap::new();
        let mut gene_labels = HashMap::new();
        for (i, (g, aliases, labels)) in records.iter().enumerate() {
            if !g.gene_id.is_empty() {
                gene_id_index.insert(g.gene_id.clone(), i);
                let mut all_aliases = vec![g.gene_id.clone()];
                for alias in aliases {
                    if !alias.is_empty() && !all_aliases.iter().any(|x| x == alias) {
                        gene_id_index.entry(alias.clone()).or_insert(i);
                        all_aliases.push(alias.clone());
                    }
                }
                let mut all_labels = all_aliases.clone();
                for label in labels {
                    if !label.is_empty() && !all_labels.iter().any(|x| x == label) {
                        all_labels.push(label.clone());
                    }
                }
                gene_aliases.insert(g.gene_id.clone(), all_aliases);
                gene_labels.insert(g.gene_id.clone(), all_labels);
            }
        }

        GffIndex {
            genes,
            contig_ranges,
            gene_id_index,
            gene_aliases,
            gene_labels,
        }
    }

    pub fn get_by_id(&self, gene_id: &str) -> Option<&GffGene> {
        self.gene_id_index.get(gene_id).map(|&i| &self.genes[i])
    }

    pub fn aliases_for_gene<'a>(&'a self, gene: &'a GffGene) -> Vec<&'a str> {
        self.gene_aliases
            .get(&gene.gene_id)
            .map(|aliases| aliases.iter().map(String::as_str).collect())
            .unwrap_or_else(|| vec![gene.gene_id.as_str()])
    }

    pub fn labels_for_gene<'a>(&'a self, gene: &'a GffGene) -> Vec<&'a str> {
        self.gene_labels
            .get(&gene.gene_id)
            .map(|labels| labels.iter().map(String::as_str).collect())
            .unwrap_or_else(|| self.aliases_for_gene(gene))
    }

    pub fn cluster_for_gene(
        &self,
        gene: &GffGene,
        gene_to_cluster: &HashMap<String, ClusterIdx>,
    ) -> Option<ClusterIdx> {
        gene_to_cluster.get(&gene.gene_id).copied().or_else(|| {
            self.gene_aliases
                .get(&gene.gene_id)
                .and_then(|aliases| aliases.iter().find_map(|alias| gene_to_cluster.get(alias)))
                .copied()
        })
    }

    pub fn overlapping(&self, contig: &str, query_start: u64, query_end: u64) -> Vec<&GffGene> {
        let Some(&(range_start, range_end)) = self.contig_ranges.get(contig) else {
            return Vec::new();
        };

        let slice = &self.genes[range_start..range_end];
        let scan_start = slice.partition_point(|g| g.end <= query_start);

        let mut results = Vec::new();
        for g in &slice[scan_start..] {
            if g.start >= query_end {
                break;
            }
            results.push(g);
        }
        results
    }

    pub fn flanking(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> (Option<&GffGene>, Option<&GffGene>) {
        let Some(&(range_start, range_end)) = self.contig_ranges.get(contig) else {
            return (None, None);
        };

        let slice = &self.genes[range_start..range_end];
        let pos = slice.partition_point(|g| g.end <= start);

        let upstream = if pos > 0 { Some(&slice[pos - 1]) } else { None };
        let downstream = slice[pos..].iter().find(|g| g.start >= end);

        (upstream, downstream)
    }

    pub fn context_window(
        &self,
        contig: &str,
        start: u64,
        end: u64,
        center_gene_id: Option<&str>,
        n_neighbors: usize,
    ) -> (Vec<&GffGene>, Vec<&GffGene>) {
        let Some(&(range_start, range_end)) = self.contig_ranges.get(contig) else {
            return (Vec::new(), Vec::new());
        };
        let slice = &self.genes[range_start..range_end];
        if slice.is_empty() {
            return (Vec::new(), Vec::new());
        }

        let center_pos = center_gene_id
            .and_then(|gene_id| self.gene_id_index.get(gene_id).copied())
            .and_then(|idx| {
                if idx >= range_start && idx < range_end {
                    Some(idx - range_start)
                } else {
                    None
                }
            })
            .or_else(|| {
                slice
                    .iter()
                    .enumerate()
                    .filter_map(|(i, gene)| {
                        let overlap = gene.end.min(end).saturating_sub(gene.start.max(start));
                        let query_len = end.saturating_sub(start);
                        let gene_len = gene.end.saturating_sub(gene.start);
                        let substantial = overlap > 0
                            && ((query_len > 0 && overlap as f64 / query_len as f64 >= 0.5)
                                || (gene_len > 0 && overlap as f64 / gene_len as f64 >= 0.5));
                        substantial.then_some((i, overlap))
                    })
                    .max_by_key(|(_, overlap)| *overlap)
                    .map(|(i, _)| i)
            });

        if let Some(center_pos) = center_pos {
            let upstream = slice[..center_pos].iter().rev().take(n_neighbors).collect();
            let downstream = slice
                .iter()
                .skip(center_pos + 1)
                .take(n_neighbors)
                .collect();
            return (upstream, downstream);
        }

        let midpoint = start.saturating_add(end.saturating_sub(start) / 2);
        let pos = slice.partition_point(|g| {
            let gene_mid = g.start.saturating_add(g.end.saturating_sub(g.start) / 2);
            gene_mid < midpoint
        });
        let upstream = slice[..pos].iter().rev().take(n_neighbors).collect();
        let downstream = slice[pos..].iter().take(n_neighbors).collect();
        (upstream, downstream)
    }
}
