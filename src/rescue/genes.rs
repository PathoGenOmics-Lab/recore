use super::gff::GffIndex;
use crate::parser::PresenceAbsence;
use crate::types::{ClusterIdx, GenomeIdx};
use ahash::{HashMap, HashMapExt};

/// Build a lookup from gene_id to ClusterIdx for a specific genome.
pub fn build_gene_to_cluster(pa: &PresenceAbsence, gi: GenomeIdx) -> HashMap<String, ClusterIdx> {
    let mut map = HashMap::new();
    for ci_val in 0..pa.n_clusters() {
        let ci = ClusterIdx(ci_val);
        if !pa.is_present(ci, gi) {
            continue;
        }
        let cell = pa.cell(ci, gi);
        for gene_id in cell.split(';') {
            let gene_id = gene_id.trim();
            if !gene_id.is_empty() {
                map.insert(gene_id.to_string(), ci);
                for variant in cell_gene_id_variants(gene_id).into_iter().skip(1) {
                    map.entry(variant.to_string()).or_insert(ci);
                }
            }
        }
    }
    map
}

pub(crate) fn strip_panaroo_gene_suffix(gene_id: &str) -> Option<&str> {
    if let Some(stripped) = gene_id.strip_suffix("_pseudo") {
        return Some(stripped);
    }
    let (base, suffix) = gene_id.rsplit_once("_dup_")?;
    if !base.is_empty() && !suffix.is_empty() && suffix.chars().all(|ch| ch.is_ascii_digit()) {
        Some(base)
    } else {
        None
    }
}

pub(crate) fn panaroo_gene_id_suffix_variants(gene_id: &str) -> Vec<&str> {
    let mut variants = vec![gene_id];
    let mut current = gene_id;
    while let Some(stripped) = strip_panaroo_gene_suffix(current) {
        if variants.contains(&stripped) {
            break;
        }
        variants.push(stripped);
        current = stripped;
    }
    variants
}

/// Every name a gene ID from a matrix cell may go by in the GFF, the ID itself
/// first. Then Panaroo's `_pseudo` and `_dup_N` suffixes stripped; the loci of
/// a PIRATE fission, which it writes `(locus1:locus2)`; and the parts after
/// each `-` of a PanTA ID, which is `sample-contig-gene`. The cells keep the
/// tool's own IDs, and a later variant only answers when the ones before it do
/// not, so an exact match always wins. Without the PanTA parts no cell gene
/// matched a GFF, and synteny validated nothing on PanTA input.
pub(crate) fn cell_gene_id_variants(gene_id: &str) -> Vec<&str> {
    let mut variants = panaroo_gene_id_suffix_variants(gene_id);
    if let Some(inner) = gene_id
        .strip_prefix('(')
        .and_then(|rest| rest.strip_suffix(')'))
    {
        for locus in inner.split(':').map(str::trim) {
            if !locus.is_empty() && !variants.contains(&locus) {
                variants.push(locus);
            }
        }
    }
    let mut rest = gene_id;
    while let Some((_, tail)) = rest.split_once('-') {
        if !tail.is_empty() && !variants.contains(&tail) {
            variants.push(tail);
        }
        rest = tail;
    }
    variants
}

fn gene_id_variants(gene_id: &str) -> impl Iterator<Item = &str> {
    panaroo_gene_id_suffix_variants(gene_id).into_iter()
}

fn lookup_gene_cluster_alias(
    map: &HashMap<String, ClusterIdx>,
    gene_id: &str,
) -> Option<ClusterIdx> {
    gene_id_variants(gene_id).find_map(|id| map.get(id).copied())
}

fn insert_gene_cluster_alias(map: &mut HashMap<String, ClusterIdx>, gene_id: &str, ci: ClusterIdx) {
    for id in gene_id_variants(gene_id) {
        if !id.is_empty() {
            map.entry(id.to_string()).or_insert(ci);
        }
    }
}

/// Build a gene_id to ClusterIdx lookup and conservatively add GFF aliases.
///
/// Panaroo cells can contain a gene name while the GFF uses a locus_tag (or the
/// reverse). If exactly one alias of a GFF CDS already maps to a cluster, all
/// aliases of that CDS are mapped to the same cluster. Ambiguous CDS aliases
/// that point to multiple clusters are left untouched.
pub fn build_gene_to_cluster_with_gff_aliases(
    pa: &PresenceAbsence,
    gi: GenomeIdx,
    gff_index: Option<&GffIndex>,
) -> HashMap<String, ClusterIdx> {
    let mut map = build_gene_to_cluster(pa, gi);
    let Some(gff) = gff_index else {
        return map;
    };

    for gene in &gff.genes {
        let labels = gff.aliases_for_gene(gene);
        let mut candidates: Vec<ClusterIdx> = labels
            .iter()
            .filter_map(|label| lookup_gene_cluster_alias(&map, label))
            .collect();
        candidates.sort_unstable();
        candidates.dedup();
        if candidates.len() != 1 {
            continue;
        }
        let ci = candidates[0];
        for label in labels {
            insert_gene_cluster_alias(&mut map, label, ci);
        }
    }

    map
}
