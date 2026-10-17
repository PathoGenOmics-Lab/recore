use super::*;

/// The index as it was before it was made compact: every name a String, one
/// map from names to genes and one map each for aliases and labels keyed by
/// gene ID. The compact index has to answer every question exactly as this did.
mod model {
    use crate::types::ClusterIdx;
    use ahash::{HashMap, HashMapExt};

    #[derive(Debug, Clone)]
    pub struct Gene {
        pub contig: String,
        pub start: u64,
        pub end: u64,
        pub strand: char,
        pub gene_id: String,
    }

    pub struct Index {
        pub genes: Vec<Gene>,
        contig_ranges: HashMap<String, (usize, usize)>,
        prefix_max_end: Vec<u64>,
        max_gene_len: HashMap<String, u64>,
        pub gene_id_index: HashMap<String, usize>,
        gene_aliases: HashMap<String, Vec<String>>,
        gene_labels: HashMap<String, Vec<String>>,
    }

    impl Index {
        pub fn build(mut records: Vec<(Gene, Vec<String>, Vec<String>)>) -> Self {
            records.sort_by(|a, b| a.0.contig.cmp(&b.0.contig).then(a.0.start.cmp(&b.0.start)));
            let genes: Vec<Gene> = records.iter().map(|(gene, _, _)| gene.clone()).collect();
            let mut contig_ranges = HashMap::new();
            if !genes.is_empty() {
                let mut current = &genes[0].contig;
                let mut range_start = 0;
                for (i, g) in genes.iter().enumerate() {
                    if &g.contig != current {
                        contig_ranges.insert(current.clone(), (range_start, i));
                        current = &g.contig;
                        range_start = i;
                    }
                }
                contig_ranges.insert(current.clone(), (range_start, genes.len()));
            }
            let mut prefix_max_end: Vec<u64> = Vec::new();
            let mut max_gene_len: HashMap<String, u64> = HashMap::new();
            for (i, g) in genes.iter().enumerate() {
                let new_contig = i == 0 || genes[i - 1].contig != g.contig;
                let running = if new_contig {
                    g.end
                } else {
                    prefix_max_end[i - 1].max(g.end)
                };
                prefix_max_end.push(running);
                let longest = max_gene_len.entry(g.contig.clone()).or_insert(0);
                *longest = (*longest).max(g.end.saturating_sub(g.start));
            }
            let mut gene_id_index = HashMap::new();
            let mut gene_aliases = HashMap::new();
            let mut gene_labels = HashMap::new();
            for (i, (g, aliases, labels)) in records.iter().enumerate() {
                if !g.gene_id.is_empty() {
                    gene_id_index.insert(g.gene_id.clone(), i);
                    let mut all_aliases = vec![g.gene_id.clone()];
                    for alias in aliases {
                        if !alias.is_empty() && !all_aliases.contains(alias) {
                            gene_id_index.entry(alias.clone()).or_insert(i);
                            all_aliases.push(alias.clone());
                        }
                    }
                    let mut all_labels = all_aliases.clone();
                    for label in labels {
                        if !label.is_empty() && !all_labels.contains(label) {
                            all_labels.push(label.clone());
                        }
                    }
                    gene_aliases.insert(g.gene_id.clone(), all_aliases);
                    gene_labels.insert(g.gene_id.clone(), all_labels);
                }
            }
            Index {
                genes,
                contig_ranges,
                prefix_max_end,
                max_gene_len,
                gene_id_index,
                gene_aliases,
                gene_labels,
            }
        }

        pub fn get_by_id(&self, gene_id: &str) -> Option<&Gene> {
            self.gene_id_index.get(gene_id).map(|&i| &self.genes[i])
        }

        pub fn aliases_for_gene<'a>(&'a self, gene: &'a Gene) -> Vec<&'a str> {
            self.gene_aliases
                .get(&gene.gene_id)
                .map(|aliases| aliases.iter().map(String::as_str).collect())
                .unwrap_or_else(|| vec![gene.gene_id.as_str()])
        }

        pub fn labels_for_gene<'a>(&'a self, gene: &'a Gene) -> Vec<&'a str> {
            self.gene_labels
                .get(&gene.gene_id)
                .map(|labels| labels.iter().map(String::as_str).collect())
                .unwrap_or_else(|| self.aliases_for_gene(gene))
        }

        pub fn cluster_for_gene(
            &self,
            gene: &Gene,
            gene_to_cluster: &HashMap<String, ClusterIdx>,
        ) -> Option<ClusterIdx> {
            gene_to_cluster.get(&gene.gene_id).copied().or_else(|| {
                self.gene_aliases
                    .get(&gene.gene_id)
                    .and_then(|aliases| aliases.iter().find_map(|alias| gene_to_cluster.get(alias)))
                    .copied()
            })
        }

        pub fn overlapping(&self, contig: &str, query_start: u64, query_end: u64) -> Vec<&Gene> {
            let Some(&(range_start, range_end)) = self.contig_ranges.get(contig) else {
                return Vec::new();
            };
            let slice = &self.genes[range_start..range_end];
            let running_end = &self.prefix_max_end[range_start..range_end];
            let scan_start = running_end.partition_point(|&end| end <= query_start);
            let mut results = Vec::new();
            for g in &slice[scan_start..] {
                if g.start >= query_end {
                    break;
                }
                if g.end > query_start {
                    results.push(g);
                }
            }
            results
        }

        pub fn flanking(
            &self,
            contig: &str,
            start: u64,
            end: u64,
        ) -> (Option<&Gene>, Option<&Gene>) {
            let Some(&(range_start, range_end)) = self.contig_ranges.get(contig) else {
                return (None, None);
            };
            let slice = &self.genes[range_start..range_end];
            let down = slice.partition_point(|g| g.start < end);
            let downstream = slice.get(down);
            let longest = self.max_gene_len.get(contig).copied().unwrap_or(0);
            let mut upstream: Option<&Gene> = None;
            for g in slice[..slice.partition_point(|g| g.start < start)]
                .iter()
                .rev()
            {
                if let Some(best) = upstream {
                    if g.start.saturating_add(longest) < best.end {
                        break;
                    }
                }
                if g.end <= start && upstream.is_none_or(|best| g.end > best.end) {
                    upstream = Some(g);
                }
            }
            (upstream, downstream)
        }

        pub fn context_window(
            &self,
            contig: &str,
            start: u64,
            end: u64,
            center_gene_id: Option<&str>,
            n_neighbors: usize,
        ) -> (Vec<&Gene>, Vec<&Gene>) {
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
}

/// xorshift64*: reproducible draws without a dependency.
struct Draws(u64);

impl Draws {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }

    fn pick<'a>(&mut self, from: &[&'a str]) -> &'a str {
        from[self.below(from.len() as u64) as usize]
    }
}

type Key = (String, u64, u64, char, String);

fn key_of(gene: &GffGene<'_>) -> Key {
    (
        gene.contig.to_string(),
        gene.start,
        gene.end,
        gene.strand,
        gene.gene_id.to_string(),
    )
}

fn model_key(gene: &model::Gene) -> Key {
    (
        gene.contig.clone(),
        gene.start,
        gene.end,
        gene.strand,
        gene.gene_id.clone(),
    )
}

const GENE_IDS: &[&str] = &["", "g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8"];
const ALIASES: &[&str] = &[
    "", "g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "a1", "a2", "a3", "a4", "a5", "a6",
];
const LABELS: &[&str] = &[
    "",
    "transposase IS6110",
    "integrase",
    "hypothetical protein",
    "phage tail protein",
    "site-specific recombinase",
    "insertion_sequence element",
];
const CONTIGS: &[&str] = &["c1", "c2", "c3"];

/// The compact index and its per-name clusters answer every question the old
/// index and the old map of names answered, the same way: names claimed by
/// two genes, aliases that are another gene's ID, genes with no ID, nested
/// genes, corrections that name one gene ID twice.
#[test]
fn compact_index_answers_like_the_maps_it_replaced() {
    let mut draws = Draws(0x9e37_79b9_7f4a_7c15);
    for trial in 0..400 {
        let n_records = draws.below(40) as usize;
        let mut records = Vec::new();
        let mut model_records = Vec::new();
        for _ in 0..n_records {
            let contig = draws.pick(CONTIGS).to_string();
            let start = draws.below(2000);
            let end = start + draws.below(600);
            let strand = if draws.below(2) == 0 { '+' } else { '-' };
            let gene_id = draws.pick(GENE_IDS).to_string();
            let aliases: Vec<String> = (0..draws.below(5))
                .map(|_| draws.pick(ALIASES).to_string())
                .collect();
            let labels: Vec<String> = (0..draws.below(3))
                .map(|_| draws.pick(LABELS).to_string())
                .collect();
            records.push((
                GffRecord {
                    contig: contig.clone(),
                    start,
                    end,
                    strand,
                    gene_id: gene_id.clone(),
                },
                aliases.clone(),
                labels.clone(),
            ));
            model_records.push((
                model::Gene {
                    contig,
                    start,
                    end,
                    strand,
                    gene_id,
                },
                aliases,
                labels,
            ));
        }
        let index = GffIndex::build_with_aliases_and_labels(records);
        let old = model::Index::build(model_records);

        let genes: Vec<GffGene> = index.genes().collect();
        let keys: Vec<Key> = genes.iter().map(key_of).collect();
        let model_keys: Vec<Key> = old.genes.iter().map(model_key).collect();
        assert_eq!(keys, model_keys, "trial {trial}: gene order");

        for name in ALIASES.iter().chain(["zz", "c1"].iter()) {
            assert_eq!(
                index.get_by_id(name).as_ref().map(key_of),
                old.get_by_id(name).map(model_key),
                "trial {trial}: get_by_id({name:?})"
            );
        }
        for (gene, model_gene) in genes.iter().zip(&old.genes) {
            assert_eq!(
                index.aliases_for_gene(gene),
                old.aliases_for_gene(model_gene),
                "trial {trial}: aliases of {model_gene:?}"
            );
            let expected = mobile_kind_from_labels(&old.labels_for_gene(model_gene));
            assert_eq!(
                index.mobile_annotation(gene),
                expected
                    .as_ref()
                    .map(|(kind, label)| (*kind, label.as_str())),
                "trial {trial}: mobile annotation of {model_gene:?}"
            );
        }

        // A map from names to clusters: some of the annotation's names, some
        // it does not use, and sometimes the empty name.
        let mut map: HashMap<String, ClusterIdx> = HashMap::new();
        for name in ALIASES.iter().chain(["x1", "x2"].iter()) {
            if draws.below(3) == 0 {
                map.insert(name.to_string(), ClusterIdx(draws.below(5) as usize));
            }
        }
        let mut clusters = GeneClusters::from_map(&index, &map);
        let check = |clusters: &GeneClusters, map: &HashMap<String, ClusterIdx>, stage: &str| {
            for name in old.gene_id_index.keys().map(String::as_str).chain([""]) {
                assert_eq!(
                    clusters.get(name),
                    map.get(name).copied(),
                    "trial {trial} {stage}: cluster of {name:?}"
                );
            }
            for (gene, model_gene) in genes.iter().zip(&old.genes) {
                assert_eq!(
                    clusters.cluster_for_gene(gene),
                    old.cluster_for_gene(model_gene, map),
                    "trial {trial} {stage}: cluster of {model_gene:?}"
                );
            }
        };
        check(&clusters, &map, "input");

        // Corrections go in as entries added only when the name is absent.
        let corrections: Vec<(&str, ClusterIdx)> = (0..draws.below(4))
            .map(|_| (draws.pick(GENE_IDS), ClusterIdx(draws.below(5) as usize)))
            .collect();
        let mut corrected = map.clone();
        for &(name, ci) in &corrections {
            corrected.entry(name.to_string()).or_insert(ci);
        }
        clusters.fill(&corrections);
        check(&clusters, &corrected, "corrected");

        for _ in 0..20 {
            let contig = draws.pick(&["c1", "c2", "c3", "c4"]);
            let start = draws.below(2400);
            let end = start + draws.below(700);
            let as_keys = |genes: Vec<GffGene>| genes.iter().map(key_of).collect::<Vec<_>>();
            let as_model_keys =
                |genes: Vec<&model::Gene>| genes.into_iter().map(model_key).collect::<Vec<_>>();
            assert_eq!(
                as_keys(index.overlapping(contig, start, end)),
                as_model_keys(old.overlapping(contig, start, end)),
                "trial {trial}: overlapping {contig}:{start}-{end}"
            );
            let (up, down) = index.flanking(contig, start, end);
            let (model_up, model_down) = old.flanking(contig, start, end);
            assert_eq!(up.as_ref().map(key_of), model_up.map(model_key));
            assert_eq!(down.as_ref().map(key_of), model_down.map(model_key));
            let center = if draws.below(2) == 0 {
                None
            } else {
                Some(draws.pick(ALIASES))
            };
            let n = draws.below(4) as usize;
            let (up, down) = index.context_window(contig, start, end, center, n);
            let (model_up, model_down) = old.context_window(contig, start, end, center, n);
            assert_eq!(
                (as_keys(up), as_keys(down)),
                (as_model_keys(model_up), as_model_keys(model_down)),
                "trial {trial}: context of {contig}:{start}-{end} around {center:?}"
            );
        }
    }
}

/// What the index costs per gene on a PGAP-style annotation, where every CDS
/// has six names: under 250 bytes, against about 1,200 when each name was a
/// String in several maps. A guard against storing names per gene again.
#[test]
fn index_stays_compact_on_pgap_style_annotation() {
    let n_genes = 4000u64;
    let records: Vec<(GffRecord, Vec<String>, Vec<String>)> = (0..n_genes)
        .map(|i| {
            let locus = format!("G00000_{i:05}");
            (
                GffRecord {
                    contig: "G00000_1".into(),
                    start: i * 1200,
                    end: i * 1200 + 1000,
                    strand: '+',
                    gene_id: locus.clone(),
                },
                vec![
                    format!("cds-{locus}"),
                    locus,
                    format!("gen{i:05}"),
                    format!("WP_{i:09}.1"),
                    format!("Rv{i:04}"),
                    format!("{}", 1_000_000 + i),
                ],
                vec!["DNA-directed RNA polymerase subunit beta".into()],
            )
        })
        .collect();
    let index = GffIndex::build_with_aliases_and_labels(records);
    let per_gene = index.heap_bytes() as u64 / n_genes;
    assert!(per_gene < 250, "{per_gene} bytes per gene");
    let map: HashMap<String, ClusterIdx> = HashMap::new();
    let clusters = GeneClusters::from_map(&index, &map);
    assert_eq!(clusters.by_name.len() as u64, 6 * n_genes);
}
