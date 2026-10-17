//! Cluster relationship reporting.
//!
//! This module does not change the corrected matrix. It summarizes evidence
//! that two Panaroo clusters may represent the same locus, a paralogous locus,
//! or a nucleotide-similarity merge candidate. The goal is reproducible review
//! rather than aggressive automatic collapsing.

use ahash::{HashMap, HashMapExt};
use anyhow::Result;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::nscluster::NsClusterResult;
use crate::parser::PresenceAbsence;
use crate::rescue::{panaroo_gene_id_suffix_variants, RescueEvent, RescueType};
use crate::types::ClusterIdx;

#[derive(Debug, Clone)]
struct RelationshipAccumulator {
    query_cluster: String,
    related_cluster: String,
    relationship_type: String,
    evidence_type: String,
    identities: Vec<f64>,
    coverages: Vec<f64>,
    synteny_pass: usize,
    synteny_fail: usize,
    synteny_not_tested: usize,
    examples: Vec<String>,
}

#[derive(Debug, Clone, Default)]
pub struct ClusterComponentInfo {
    pub component_id: String,
    pub component_size: usize,
    pub role: String,
    pub relation_types: String,
    pub recommendation: String,
}

#[derive(Debug, Clone)]
struct ComponentEdge {
    a: usize,
    b: usize,
    relation_type: String,
}

#[derive(Debug, Clone)]
struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            let root = self.find(self.parent[x]);
            self.parent[x] = root;
        }
        self.parent[x]
    }

    fn union(&mut self, a: usize, b: usize) {
        let mut ra = self.find(a);
        let mut rb = self.find(b);
        if ra == rb {
            return;
        }
        if self.rank[ra] < self.rank[rb] {
            std::mem::swap(&mut ra, &mut rb);
        }
        self.parent[rb] = ra;
        if self.rank[ra] == self.rank[rb] {
            self.rank[ra] += 1;
        }
    }
}

impl RelationshipAccumulator {
    fn new(
        query_cluster: String,
        related_cluster: String,
        relationship_type: String,
        evidence_type: String,
    ) -> Self {
        Self {
            query_cluster,
            related_cluster,
            relationship_type,
            evidence_type,
            identities: Vec::new(),
            coverages: Vec::new(),
            synteny_pass: 0,
            synteny_fail: 0,
            synteny_not_tested: 0,
            examples: Vec::new(),
        }
    }

    fn add_rescue(&mut self, rescue: &RescueEvent, genome_name: &str) {
        self.identities.push(rescue.identity);
        self.coverages.push(rescue.coverage);
        match rescue.synteny_ok {
            Some(true) => self.synteny_pass += 1,
            Some(false) => self.synteny_fail += 1,
            None => self.synteny_not_tested += 1,
        }
        if self.examples.len() < 5 {
            self.examples.push(format!(
                "{}:{}:{}-{}",
                genome_name,
                rescue.contig,
                rescue.start + 1,
                rescue.end
            ));
        }
    }
}

// Keep the documented Rust MSRV; usize::is_multiple_of is newer.
#[allow(clippy::manual_is_multiple_of)]
fn median(values: &mut [f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    values.sort_by(|a, b| a.total_cmp(b));
    let mid = values.len() / 2;
    if values.len() % 2 == 0 {
        (values[mid - 1] + values[mid]) / 2.0
    } else {
        values[mid]
    }
}

fn fmt_float(value: f64) -> String {
    if value.is_nan() {
        String::new()
    } else {
        format!("{:.1}", value)
    }
}

fn clean_tsv_field(value: &str) -> String {
    value
        .chars()
        .map(|ch| match ch {
            '\t' | '\n' | '\r' | '\0' => ' ',
            _ => ch,
        })
        .collect()
}

fn write_tsv_fields(w: &mut BufWriter<fs::File>, fields: &[String]) -> Result<()> {
    let cleaned: Vec<String> = fields.iter().map(|field| clean_tsv_field(field)).collect();
    writeln!(w, "{}", cleaned.join("\t"))?;
    Ok(())
}

fn synteny_support(pass: usize, fail: usize, not_tested: usize) -> &'static str {
    match (pass, fail, not_tested) {
        (0, 0, _) => "not_tested",
        (p, 0, _) if p > 0 => "pass",
        (0, f, _) if f > 0 => "fail",
        _ => "mixed",
    }
}

fn recommendation(
    relationship_type: &str,
    identity_median: f64,
    coverage_median: f64,
    synteny: &str,
) -> &'static str {
    match relationship_type {
        "nucleotide_similarity_cluster" => "collapse_alias",
        "split" if identity_median >= 95.0 && coverage_median >= 85.0 && synteny == "pass" => {
            "collapse_alias"
        }
        "split" if identity_median >= 95.0 && coverage_median >= 85.0 => {
            "review_alias_without_full_synteny"
        }
        "split" => "review_split",
        "annotation_match_cluster"
            if identity_median >= 95.0 && coverage_median >= 85.0 && synteny == "pass" =>
        {
            "collapse_alias"
        }
        "annotation_match_cluster" if identity_median >= 95.0 && coverage_median >= 85.0 => {
            "review_alias_without_full_synteny"
        }
        "annotation_match_cluster" => "review_annotation_match",
        "paralog" => "review_paralog_keep_distinct_by_default",
        _ => "review",
    }
}

/// Write cluster-level relationships inferred from rescue classifications and
/// optional nucleotide similarity clustering.
///
/// Returns the number of relationship rows written. The file is always written
/// so downstream workflows can depend on its presence.
pub fn write_cluster_relationships(
    pa: &PresenceAbsence,
    rescues: &[RescueEvent],
    nsc_result: &Option<NsClusterResult>,
    output_dir: &Path,
) -> Result<usize> {
    let mut accs: HashMap<(String, String, String), RelationshipAccumulator> = HashMap::new();

    for rescue in rescues {
        let (related, rel_type, evidence_type) = match &rescue.rescue_type {
            RescueType::Split { other_cluster } => (
                other_cluster.clone(),
                "split".to_string(),
                rescue.rescue_type.to_string(),
            ),
            RescueType::Paralog { other_cluster } => (
                other_cluster.clone(),
                "paralog".to_string(),
                rescue.rescue_type.to_string(),
            ),
            _ => continue,
        };
        let query = pa.cluster_name(rescue.cluster_idx).to_string();
        let key = (query.clone(), related.clone(), rel_type.clone());
        let genome_name = pa.genome_name(rescue.genome_idx);
        accs.entry(key)
            .or_insert_with(|| {
                RelationshipAccumulator::new(query, related, rel_type, evidence_type)
            })
            .add_rescue(rescue, genome_name);
    }

    let mut gene_to_cluster_by_genome: HashMap<usize, HashMap<String, usize>> = HashMap::new();
    for rescue in rescues {
        if matches!(
            rescue.rescue_type,
            RescueType::Split { .. } | RescueType::Paralog { .. }
        ) {
            continue;
        }
        let Some(annotation_match) = rescue.annotation_match.as_deref() else {
            continue;
        };
        let genome_map = gene_to_cluster_by_genome
            .entry(rescue.genome_idx.val())
            .or_insert_with(|| build_gene_cluster_map(pa, rescue.genome_idx));
        let Some(&related_idx) = genome_map.get(annotation_match) else {
            continue;
        };
        if related_idx == rescue.cluster_idx.val() {
            continue;
        }

        let query = pa.cluster_name(rescue.cluster_idx).to_string();
        let related = pa
            .cluster_name(crate::types::ClusterIdx(related_idx))
            .to_string();
        let rel_type = "annotation_match_cluster".to_string();
        let key = (query.clone(), related.clone(), rel_type.clone());
        let genome_name = pa.genome_name(rescue.genome_idx);
        accs.entry(key)
            .or_insert_with(|| {
                RelationshipAccumulator::new(
                    query,
                    related,
                    rel_type,
                    "annotation_match".to_string(),
                )
            })
            .add_rescue(rescue, genome_name);
    }

    if let Some(nsc) = nsc_result {
        for (&member, &representative) in &nsc.merge_map {
            let query = pa.cluster_name(member).to_string();
            let related = pa.cluster_name(representative).to_string();
            let rel_type = "nucleotide_similarity_cluster".to_string();
            let key = (query.clone(), related.clone(), rel_type.clone());
            accs.entry(key).or_insert_with(|| {
                RelationshipAccumulator::new(query, related, rel_type, "nscluster".to_string())
            });
        }
    }

    let mut rows: Vec<_> = accs.into_values().collect();
    rows.sort_by(|a, b| {
        a.query_cluster
            .cmp(&b.query_cluster)
            .then(a.related_cluster.cmp(&b.related_cluster))
            .then(a.relationship_type.cmp(&b.relationship_type))
    });
    let row_count = rows.len();

    let path = output_dir.join("cluster_relationships.tsv");
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "query_cluster\trelated_cluster\trelationship_type\tevidence_type\tn_genomes_supporting\tsynteny_support\tsynteny_pass_n\tsynteny_fail_n\tsynteny_not_tested_n\tidentity_min\tidentity_median\tidentity_max\tcoverage_min\tcoverage_median\tcoverage_max\trecommendation\texample_regions"
    )?;

    for mut row in rows {
        let n = row.identities.len();
        let identity_min = row
            .identities
            .iter()
            .copied()
            .reduce(f64::min)
            .unwrap_or(f64::NAN);
        let identity_max = row
            .identities
            .iter()
            .copied()
            .reduce(f64::max)
            .unwrap_or(f64::NAN);
        let coverage_min = row
            .coverages
            .iter()
            .copied()
            .reduce(f64::min)
            .unwrap_or(f64::NAN);
        let coverage_max = row
            .coverages
            .iter()
            .copied()
            .reduce(f64::max)
            .unwrap_or(f64::NAN);
        let identity_median = median(&mut row.identities);
        let coverage_median = median(&mut row.coverages);
        let synteny = synteny_support(row.synteny_pass, row.synteny_fail, row.synteny_not_tested);
        let rec = recommendation(
            &row.relationship_type,
            identity_median,
            coverage_median,
            synteny,
        );

        write_tsv_fields(
            &mut w,
            &[
                row.query_cluster,
                row.related_cluster,
                row.relationship_type,
                row.evidence_type,
                n.to_string(),
                synteny.to_string(),
                row.synteny_pass.to_string(),
                row.synteny_fail.to_string(),
                row.synteny_not_tested.to_string(),
                fmt_float(identity_min),
                fmt_float(identity_median),
                fmt_float(identity_max),
                fmt_float(coverage_min),
                fmt_float(coverage_median),
                fmt_float(coverage_max),
                rec.to_string(),
                row.examples.join(";"),
            ],
        )?;
    }
    w.flush()?;
    Ok(row_count)
}

pub fn write_cluster_components(
    pa: &PresenceAbsence,
    rescues: &[RescueEvent],
    nsc_result: &Option<NsClusterResult>,
    output_dir: &Path,
) -> Result<Vec<Option<ClusterComponentInfo>>> {
    let name_to_idx = pa.cluster_index();
    let mut edges = Vec::new();

    for rescue in rescues {
        let related = match &rescue.rescue_type {
            RescueType::Split { other_cluster } => Some((other_cluster.as_str(), "split")),
            RescueType::Paralog { other_cluster } => Some((other_cluster.as_str(), "paralog")),
            _ => None,
        };
        if let Some((related_name, relation_type)) = related {
            if let Some(&other_ci) = name_to_idx.get(related_name) {
                edges.push(ComponentEdge {
                    a: rescue.cluster_idx.val(),
                    b: other_ci.val(),
                    relation_type: relation_type.to_string(),
                });
            }
        }
    }

    let mut gene_to_cluster_by_genome: HashMap<usize, HashMap<String, usize>> = HashMap::new();
    for rescue in rescues {
        if matches!(
            rescue.rescue_type,
            RescueType::Split { .. } | RescueType::Paralog { .. }
        ) {
            continue;
        }
        let Some(annotation_match) = rescue.annotation_match.as_deref() else {
            continue;
        };
        let genome_map = gene_to_cluster_by_genome
            .entry(rescue.genome_idx.val())
            .or_insert_with(|| build_gene_cluster_map(pa, rescue.genome_idx));
        let Some(&related_idx) = genome_map.get(annotation_match) else {
            continue;
        };
        if related_idx != rescue.cluster_idx.val() {
            edges.push(ComponentEdge {
                a: rescue.cluster_idx.val(),
                b: related_idx,
                relation_type: "annotation_match_cluster".to_string(),
            });
        }
    }

    let mut nsc_members = BTreeSet::new();
    let mut nsc_representatives = BTreeSet::new();
    if let Some(nsc) = nsc_result {
        for (&member, &representative) in &nsc.merge_map {
            edges.push(ComponentEdge {
                a: member.val(),
                b: representative.val(),
                relation_type: "nucleotide_similarity_cluster".to_string(),
            });
            nsc_members.insert(member.val());
            nsc_representatives.insert(representative.val());
        }
    }

    let mut uf = UnionFind::new(pa.n_clusters());
    for edge in &edges {
        uf.union(edge.a, edge.b);
    }

    let mut members_by_root: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for ci_val in 0..pa.n_clusters() {
        let root = uf.find(ci_val);
        members_by_root.entry(root).or_default().push(ci_val);
    }
    let mut relation_types_by_root: BTreeMap<usize, BTreeSet<String>> = BTreeMap::new();
    for edge in &edges {
        let root = uf.find(edge.a);
        relation_types_by_root
            .entry(root)
            .or_default()
            .insert(edge.relation_type.clone());
    }

    let component_path = output_dir.join("cluster_components.tsv");
    let component_file = fs::File::create(&component_path)?;
    let mut component_w = BufWriter::new(component_file);
    writeln!(
        component_w,
        "component_id\tn_clusters\trelation_types\trecommendation\tmember_clusters"
    )?;

    let member_path = output_dir.join("cluster_component_members.tsv");
    let member_file = fs::File::create(&member_path)?;
    let mut member_w = BufWriter::new(member_file);
    writeln!(
        member_w,
        "gene_cluster\tcomponent_id\tcomponent_size\tcomponent_role\trelation_types\trecommendation"
    )?;

    let mut component_info = vec![None; pa.n_clusters()];
    let mut component_idx = 1usize;
    for (root, mut members) in members_by_root {
        if members.len() < 2 {
            continue;
        }
        members.sort_unstable();
        let relation_types = relation_types_by_root
            .get(&root)
            .map(|types| types.iter().cloned().collect::<Vec<_>>().join(";"))
            .unwrap_or_default();
        let recommendation = component_recommendation(&relation_types);
        let component_id = format!("component_{component_idx:04}");
        component_idx += 1;

        write_tsv_fields(
            &mut component_w,
            &[
                component_id.clone(),
                members.len().to_string(),
                relation_types.clone(),
                recommendation.to_string(),
                members
                    .iter()
                    .map(|&ci| pa.cluster_name(ClusterIdx(ci)))
                    .collect::<Vec<_>>()
                    .join(";"),
            ],
        )?;

        for &ci in &members {
            let role = component_role(ci, &relation_types, &nsc_members, &nsc_representatives);
            let info = ClusterComponentInfo {
                component_id: component_id.clone(),
                component_size: members.len(),
                role: role.to_string(),
                relation_types: relation_types.clone(),
                recommendation: recommendation.to_string(),
            };
            write_tsv_fields(
                &mut member_w,
                &[
                    pa.cluster_name(ClusterIdx(ci)).to_string(),
                    info.component_id.clone(),
                    info.component_size.to_string(),
                    info.role.clone(),
                    info.relation_types.clone(),
                    info.recommendation.clone(),
                ],
            )?;
            component_info[ci] = Some(info);
        }
    }

    component_w.flush()?;
    member_w.flush()?;
    Ok(component_info)
}

fn component_recommendation(relation_types: &str) -> &'static str {
    if relation_types.contains("paralog") {
        "review_component_keep_paralogs_distinct"
    } else if relation_types.contains("nucleotide_similarity_cluster")
        || relation_types.contains("split")
        || relation_types.contains("annotation_match_cluster")
    {
        "review_or_collapse_aliases"
    } else {
        "review"
    }
}

fn component_role(
    ci: usize,
    relation_types: &str,
    nsc_members: &BTreeSet<usize>,
    nsc_representatives: &BTreeSet<usize>,
) -> &'static str {
    if relation_types.contains("paralog") {
        "paralog_component_member"
    } else if nsc_representatives.contains(&ci) {
        "nsc_representative"
    } else if nsc_members.contains(&ci) {
        "nsc_member"
    } else {
        "related_cluster_member"
    }
}

fn build_gene_cluster_map(
    pa: &PresenceAbsence,
    gi: crate::types::GenomeIdx,
) -> HashMap<String, usize> {
    let mut map = HashMap::new();
    for ci_val in 0..pa.n_clusters() {
        let ci = crate::types::ClusterIdx(ci_val);
        if !pa.is_present(ci, gi) {
            continue;
        }
        for gene_id in pa.cell(ci, gi).split(';') {
            let gene_id = gene_id.trim();
            if gene_id.is_empty() {
                continue;
            }
            map.insert(gene_id.to_string(), ci_val);
            for stripped in panaroo_gene_id_suffix_variants(gene_id).into_iter().skip(1) {
                map.entry(stripped.to_string()).or_insert(ci_val);
            }
        }
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::PresenceAbsence;
    use crate::rescue::{FragmentEvidence, RescueType, SyntenyEvidence};
    use crate::types::{ClusterIdx, GenomeIdx};
    use tempfile::tempdir;

    #[test]
    fn split_relationship_is_summarized() {
        let pa = PresenceAbsence::new_for_test(
            vec!["group_a".into(), "group_b".into()],
            vec!["genome_1".into()],
            vec![vec!["".into()], vec!["gene_b".into()]],
        );
        let rescues = vec![RescueEvent {
            cluster_idx: ClusterIdx(0),
            genome_idx: GenomeIdx(0),
            rescue_type: RescueType::Split {
                other_cluster: "group_b".into(),
            },
            identity: 99.2,
            coverage: 98.5,
            contig: "ctg".into(),
            start: 10,
            end: 110,
            strand: '+',
            annotation_match: Some("gene_b".into()),
            confidence: 1.1,
            synteny_ok: Some(true),
            synteny_evidence: SyntenyEvidence::default(),
            fragment_evidence: FragmentEvidence::default(),
            query_start: None,
        }];
        let dir = tempdir().unwrap();
        let n = write_cluster_relationships(&pa, &rescues, &None, dir.path()).unwrap();
        assert_eq!(n, 1);
        let text = std::fs::read_to_string(dir.path().join("cluster_relationships.tsv")).unwrap();
        assert!(text.contains("group_a\tgroup_b\tsplit"));
        assert!(text.contains("collapse_alias"));
    }

    #[test]
    fn annotation_match_relationship_is_summarized() {
        let pa = PresenceAbsence::new_for_test(
            vec!["group_a".into(), "group_b".into()],
            vec!["genome_1".into()],
            vec![vec!["".into()], vec!["gene_b".into()]],
        );
        let rescues = vec![RescueEvent {
            cluster_idx: ClusterIdx(0),
            genome_idx: GenomeIdx(0),
            rescue_type: RescueType::Dropped,
            identity: 98.0,
            coverage: 96.0,
            contig: "ctg".into(),
            start: 10,
            end: 110,
            strand: '+',
            annotation_match: Some("gene_b".into()),
            confidence: 1.0,
            synteny_ok: None,
            synteny_evidence: SyntenyEvidence::default(),
            fragment_evidence: FragmentEvidence::default(),
            query_start: None,
        }];
        let dir = tempdir().unwrap();
        let n = write_cluster_relationships(&pa, &rescues, &None, dir.path()).unwrap();
        assert_eq!(n, 1);
        let text = std::fs::read_to_string(dir.path().join("cluster_relationships.tsv")).unwrap();
        assert!(text.contains("group_a\tgroup_b\tannotation_match_cluster"));
        assert!(text.contains("review_alias_without_full_synteny"));
    }

    #[test]
    fn cluster_components_group_related_clusters() {
        let pa = PresenceAbsence::new_for_test(
            vec!["group_a".into(), "group_b".into(), "group_c".into()],
            vec!["genome_1".into()],
            vec![
                vec!["".into()],
                vec!["gene_b".into()],
                vec!["gene_c".into()],
            ],
        );
        let rescues = vec![
            RescueEvent {
                cluster_idx: ClusterIdx(0),
                genome_idx: GenomeIdx(0),
                rescue_type: RescueType::Split {
                    other_cluster: "group_b".into(),
                },
                identity: 99.2,
                coverage: 98.5,
                contig: "ctg".into(),
                start: 10,
                end: 110,
                strand: '+',
                annotation_match: Some("gene_b".into()),
                confidence: 1.1,
                synteny_ok: Some(true),
                synteny_evidence: SyntenyEvidence::default(),
                fragment_evidence: FragmentEvidence::default(),
                query_start: None,
            },
            RescueEvent {
                cluster_idx: ClusterIdx(2),
                genome_idx: GenomeIdx(0),
                rescue_type: RescueType::Dropped,
                identity: 99.2,
                coverage: 98.5,
                contig: "ctg".into(),
                start: 210,
                end: 310,
                strand: '+',
                annotation_match: None,
                confidence: 1.1,
                synteny_ok: Some(true),
                synteny_evidence: SyntenyEvidence::default(),
                fragment_evidence: FragmentEvidence::default(),
                query_start: None,
            },
        ];
        let dir = tempdir().unwrap();
        let components = write_cluster_components(&pa, &rescues, &None, dir.path()).unwrap();
        assert_eq!(components[0].as_ref().unwrap().component_size, 2);
        assert_eq!(components[1].as_ref().unwrap().component_size, 2);
        assert!(components[2].is_none());
        let text = std::fs::read_to_string(dir.path().join("cluster_components.tsv")).unwrap();
        assert!(text.contains("group_a;group_b"));
        assert!(text.contains("review_or_collapse_aliases"));
    }
}
