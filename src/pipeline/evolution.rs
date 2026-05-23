use super::reference_audit::ReferenceClusterEvent;
use super::*;

#[derive(Debug, Clone)]
pub(super) struct PhyloNode {
    pub(super) label: Option<String>,
    pub(super) genome_idx: Option<GenomeIdx>,
    pub(super) children: Vec<PhyloNode>,
}

#[derive(Debug, Clone)]
pub(super) struct PhyloTree {
    pub(super) root: PhyloNode,
}

#[derive(Debug, Clone, Default)]
pub(super) struct PhyloPatternSummary {
    pub(super) assessed: bool,
    pub(super) matched_tips: usize,
    pub(super) present_tips: usize,
    pub(super) absent_tips: usize,
    pub(super) parsimony_changes: usize,
    pub(super) gains: usize,
    pub(super) losses: usize,
    pub(super) root_state: &'static str,
    pub(super) carrier_monophyly: bool,
    pub(super) max_carrier_clade_fraction: f64,
    pub(super) pattern: &'static str,
}

#[derive(Debug, Clone, Default)]
pub(super) struct FitchNodeSummary {
    pub(super) state_set: u8,
    pub(super) matched_tips: usize,
    pub(super) present_tips: usize,
    pub(super) changes: usize,
    pub(super) children: Vec<FitchNodeSummary>,
}

pub(super) struct NewickParser<'a> {
    pub(super) input: &'a str,
    pub(super) pos: usize,
}

impl<'a> NewickParser<'a> {
    pub(super) fn new(input: &'a str) -> Self {
        Self { input, pos: 0 }
    }

    pub(super) fn parse(mut self) -> Result<PhyloNode> {
        self.skip_ignored();
        let root = self.parse_subtree()?;
        self.skip_ignored();
        if self.peek() == Some(';') {
            self.pos += 1;
        }
        self.skip_ignored();
        if self.pos < self.input.len() {
            anyhow::bail!("Unexpected Newick content at byte {}", self.pos);
        }
        Ok(root)
    }

    fn parse_subtree(&mut self) -> Result<PhyloNode> {
        self.skip_ignored();
        if self.peek() == Some('(') {
            self.pos += 1;
            let mut children = Vec::new();
            loop {
                children.push(self.parse_subtree()?);
                self.skip_ignored();
                match self.peek() {
                    Some(',') => {
                        self.pos += 1;
                    }
                    Some(')') => {
                        self.pos += 1;
                        break;
                    }
                    other => anyhow::bail!(
                        "Expected ',' or ')' in Newick at byte {}, got {:?}",
                        self.pos,
                        other
                    ),
                }
            }
            let label = self.parse_optional_label();
            self.skip_branch_length();
            Ok(PhyloNode {
                label,
                genome_idx: None,
                children,
            })
        } else {
            let label = self.parse_optional_label();
            if label.as_deref().unwrap_or("").is_empty() {
                anyhow::bail!("Expected Newick leaf label at byte {}", self.pos);
            }
            self.skip_branch_length();
            Ok(PhyloNode {
                label,
                genome_idx: None,
                children: Vec::new(),
            })
        }
    }

    fn parse_optional_label(&mut self) -> Option<String> {
        self.skip_ignored();
        match self.peek() {
            Some('\'') => {
                self.pos += 1;
                let start = self.pos;
                while let Some(ch) = self.peek() {
                    if ch == '\'' {
                        let label = self.input[start..self.pos].to_string();
                        self.pos += 1;
                        return Some(label);
                    }
                    self.pos += ch.len_utf8();
                }
                Some(self.input[start..].to_string())
            }
            Some(':') | Some(',') | Some(')') | Some(';') | None => None,
            _ => {
                let start = self.pos;
                while let Some(ch) = self.peek() {
                    if matches!(ch, ':' | ',' | ')' | ';' | '[') || ch.is_whitespace() {
                        break;
                    }
                    self.pos += ch.len_utf8();
                }
                (self.pos > start).then(|| self.input[start..self.pos].to_string())
            }
        }
    }

    fn skip_branch_length(&mut self) {
        self.skip_ignored();
        if self.peek() != Some(':') {
            return;
        }
        self.pos += 1;
        while let Some(ch) = self.peek() {
            if matches!(ch, ',' | ')' | ';' | '[') || ch.is_whitespace() {
                break;
            }
            self.pos += ch.len_utf8();
        }
    }

    fn skip_ignored(&mut self) {
        loop {
            while self.peek().is_some_and(char::is_whitespace) {
                self.pos += self.peek().unwrap().len_utf8();
            }
            if self.peek() != Some('[') {
                break;
            }
            self.pos += 1;
            while let Some(ch) = self.peek() {
                self.pos += ch.len_utf8();
                if ch == ']' {
                    break;
                }
            }
        }
    }

    fn peek(&self) -> Option<char> {
        self.input[self.pos..].chars().next()
    }
}

pub(super) fn load_phylo_tree(
    path: &Path,
    pa: &crate::parser::PresenceAbsence,
) -> Result<PhyloTree> {
    let text =
        fs::read_to_string(path).with_context(|| format!("Cannot read tree {}", path.display()))?;
    let mut root = NewickParser::new(&text)
        .parse()
        .with_context(|| format!("Cannot parse Newick tree {}", path.display()))?;
    let (matched_tips, unmatched_tips) = resolve_phylo_tips(&mut root, pa);
    info!(
        "Loaded tree {} with {} matched tips and {} unmatched tips",
        path.display(),
        matched_tips,
        unmatched_tips
    );
    if matched_tips == 0 {
        warn!("No tree tips matched matrix genome names; phylogenetic audit will be uninformative");
    } else if matched_tips < pa.n_genomes() {
        warn!(
            "Tree tip matching is incomplete: {}/{} matrix genomes matched; unmatched genomes are ignored in evolutionary event calls",
            matched_tips,
            pa.n_genomes()
        );
    }
    Ok(PhyloTree { root })
}

pub(super) fn resolve_phylo_tips(
    node: &mut PhyloNode,
    pa: &crate::parser::PresenceAbsence,
) -> (usize, usize) {
    if node.children.is_empty() {
        node.genome_idx = node
            .label
            .as_deref()
            .and_then(|label| resolve_tree_tip(label, pa));
        return if node.genome_idx.is_some() {
            (1, 0)
        } else {
            (0, 1)
        };
    }
    let mut matched = 0usize;
    let mut unmatched = 0usize;
    for child in &mut node.children {
        let (m, u) = resolve_phylo_tips(child, pa);
        matched += m;
        unmatched += u;
    }
    (matched, unmatched)
}

pub(super) fn resolve_tree_tip(
    label: &str,
    pa: &crate::parser::PresenceAbsence,
) -> Option<GenomeIdx> {
    for gi_val in 0..pa.n_genomes() {
        let gi = GenomeIdx(gi_val);
        if pa.genome_name(gi) == label {
            return Some(gi);
        }
    }
    let tip_as = extract_as_id(label);
    for gi_val in 0..pa.n_genomes() {
        let gi = GenomeIdx(gi_val);
        if let (Some(a), Some(b)) = (tip_as, extract_as_id(pa.genome_name(gi))) {
            if a == b {
                return Some(gi);
            }
        }
    }
    None
}

pub(super) fn fitch_postorder(node: &PhyloNode, present_by_genome: &[bool]) -> FitchNodeSummary {
    const ABSENT: u8 = 0b01;
    const PRESENT: u8 = 0b10;

    if node.children.is_empty() {
        let Some(gi) = node.genome_idx else {
            return FitchNodeSummary::default();
        };
        let present = present_by_genome.get(gi.val()).copied().unwrap_or(false);
        return FitchNodeSummary {
            state_set: if present { PRESENT } else { ABSENT },
            matched_tips: 1,
            present_tips: usize::from(present),
            changes: 0,
            children: Vec::new(),
        };
    }

    let mut children = Vec::new();
    let mut state_set = 0u8;
    let mut matched_tips = 0usize;
    let mut present_tips = 0usize;
    let mut changes = 0usize;
    for child in &node.children {
        let summary = fitch_postorder(child, present_by_genome);
        if summary.matched_tips == 0 {
            children.push(summary);
            continue;
        }
        matched_tips += summary.matched_tips;
        present_tips += summary.present_tips;
        changes += summary.changes;
        if state_set == 0 {
            state_set = summary.state_set;
        } else {
            let intersection = state_set & summary.state_set;
            if intersection == 0 {
                state_set |= summary.state_set;
                changes += 1;
            } else {
                state_set = intersection;
            }
        }
        children.push(summary);
    }

    FitchNodeSummary {
        state_set,
        matched_tips,
        present_tips,
        changes,
        children,
    }
}

pub(super) fn state_set_contains(state_set: u8, present: bool) -> bool {
    let bit = if present { 0b10 } else { 0b01 };
    state_set & bit != 0
}

pub(super) fn choose_root_state(summary: &FitchNodeSummary) -> bool {
    if summary.state_set == 0b10 {
        true
    } else if summary.state_set == 0b01 {
        false
    } else {
        summary.present_tips * 2 >= summary.matched_tips
    }
}

pub(super) fn assign_fitch_states(
    summary: &FitchNodeSummary,
    parent_state: Option<bool>,
    gains: &mut usize,
    losses: &mut usize,
) -> Option<bool> {
    if summary.matched_tips == 0 || summary.state_set == 0 {
        return None;
    }
    let state = if let Some(parent) = parent_state {
        if state_set_contains(summary.state_set, parent) {
            parent
        } else {
            !parent
        }
    } else {
        choose_root_state(summary)
    };

    if let Some(parent) = parent_state {
        if !parent && state {
            *gains += 1;
        } else if parent && !state {
            *losses += 1;
        }
    }
    for child in &summary.children {
        assign_fitch_states(child, Some(state), gains, losses);
    }
    Some(state)
}

pub(super) fn carrier_clade_stats(summary: &FitchNodeSummary, total_present: usize) -> (bool, f64) {
    if total_present == 0 {
        return (false, 0.0);
    }
    if total_present == 1 {
        return (true, 1.0);
    }

    fn visit(summary: &FitchNodeSummary, total_present: usize, best: &mut usize, exact: &mut bool) {
        if summary.matched_tips > 0
            && summary.present_tips > 0
            && summary.present_tips == summary.matched_tips
        {
            *best = (*best).max(summary.present_tips);
            if summary.present_tips == total_present {
                *exact = true;
            }
        }
        for child in &summary.children {
            visit(child, total_present, best, exact);
        }
    }

    let mut best = 0usize;
    let mut exact = false;
    visit(summary, total_present, &mut best, &mut exact);
    (exact, best as f64 / total_present as f64)
}

pub(super) fn classify_phylogenetic_pattern(summary: &PhyloPatternSummary) -> &'static str {
    if !summary.assessed {
        return "not_assessed_no_tree";
    }
    if summary.present_tips == 0 {
        return "absent_from_tree";
    }
    if summary.absent_tips == 0 {
        return "ancestral_retained";
    }
    if summary.present_tips == 1 {
        return "singleton_or_noise";
    }
    if summary.parsimony_changes >= 6 {
        return "homoplastic_gain_loss";
    }
    if summary.root_state == "present" && summary.gains == 0 && summary.losses > 0 {
        return "ancestral_loss";
    }
    if summary.root_state == "absent" && summary.losses == 0 && summary.gains == 1 {
        if summary.carrier_monophyly {
            return "lineage_specific_gain";
        }
        return "single_gain_with_incomplete_clade_support";
    }
    if summary.gains > 1 && summary.losses == 0 {
        return "recurrent_gain";
    }
    if summary.losses > 1 && summary.gains == 0 {
        return "recurrent_loss";
    }
    if summary.carrier_monophyly {
        return "clade_associated_presence";
    }
    "dispersed_gain_loss"
}

pub(super) fn summarize_phylogenetic_pattern(
    tree: Option<&PhyloTree>,
    present_by_genome: &[bool],
) -> PhyloPatternSummary {
    let Some(tree) = tree else {
        return PhyloPatternSummary {
            pattern: "not_assessed_no_tree",
            root_state: "unknown",
            ..PhyloPatternSummary::default()
        };
    };
    let fitch = fitch_postorder(&tree.root, present_by_genome);
    if fitch.matched_tips == 0 {
        return PhyloPatternSummary {
            pattern: "not_assessed_no_matched_tips",
            root_state: "unknown",
            ..PhyloPatternSummary::default()
        };
    }
    let mut gains = 0usize;
    let mut losses = 0usize;
    let root_present = assign_fitch_states(&fitch, None, &mut gains, &mut losses).unwrap_or(false);
    let (carrier_monophyly, max_carrier_clade_fraction) =
        carrier_clade_stats(&fitch, fitch.present_tips);
    let mut summary = PhyloPatternSummary {
        assessed: true,
        matched_tips: fitch.matched_tips,
        present_tips: fitch.present_tips,
        absent_tips: fitch.matched_tips.saturating_sub(fitch.present_tips),
        parsimony_changes: fitch.changes,
        gains,
        losses,
        root_state: if root_present { "present" } else { "absent" },
        carrier_monophyly,
        max_carrier_clade_fraction,
        pattern: "unclassified",
    };
    summary.pattern = classify_phylogenetic_pattern(&summary);
    summary
}

#[derive(Debug, Clone)]
pub(super) struct MobileElementHit {
    pub(super) genome_idx: GenomeIdx,
    pub(super) query_name: String,
    pub(super) contig: String,
    pub(super) start: u64,
    pub(super) end: u64,
    pub(super) strand: char,
    pub(super) identity: f64,
    pub(super) coverage: f64,
}

#[derive(Debug, Clone)]
pub(super) struct MobileElementEvidence {
    pub(super) genome_name: String,
    pub(super) source: &'static str,
    pub(super) kind: &'static str,
    pub(super) label: String,
    pub(super) distance_bp: u64,
    pub(super) identity: Option<f64>,
    pub(super) coverage: Option<f64>,
}

#[derive(Debug, Clone)]
pub(super) struct MechanismCall {
    pub(super) mechanism: &'static str,
    pub(super) detail: String,
    pub(super) support: &'static str,
    pub(super) wording: &'static str,
    pub(super) mobile: Option<MobileElementEvidence>,
}

pub(super) struct EvolutionaryEventInputs<'a> {
    pub(super) pa: &'a crate::parser::PresenceAbsence,
    pub(super) matrix_rescue_map: &'a HashMap<(ClusterIdx, GenomeIdx), &'a RescueEvent>,
    pub(super) rescues: &'a [RescueEvent],
    pub(super) synteny_rejected: &'a [RescueEvent],
    pub(super) paralog_rejected: &'a [RescueEvent],
    pub(super) reference_cluster_events: &'a [ReferenceClusterEvent],
    pub(super) context_summaries: &'a [ClusterContextSummary],
    pub(super) gff_indices: &'a HashMap<GenomeIdx, GffIndex>,
    pub(super) mobile_element_hits: &'a [MobileElementHit],
    pub(super) mobile_element_window: u64,
    pub(super) phylo_tree: Option<&'a PhyloTree>,
    pub(super) reference_genome: Option<&'a str>,
    pub(super) nsc_result: &'a Option<nscluster::NsClusterResult>,
    pub(super) output_dir: &'a Path,
}

pub(super) fn nsc_presence_updates(
    pa: &crate::parser::PresenceAbsence,
    nsc_result: &Option<nscluster::NsClusterResult>,
    rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
) -> ahash::HashSet<(usize, usize)> {
    nsc_result
        .as_ref()
        .map(|r| {
            let (_, updates) = nscluster::apply_ns_merges(
                &r.merge_map,
                |ci, gi| pa.is_present(ci, gi) || rescue_map.contains_key(&(ci, gi)),
                pa.n_genomes(),
            );
            updates
                .into_iter()
                .map(|(ci, gi)| (ci.val(), gi.val()))
                .collect()
        })
        .unwrap_or_default()
}

pub(super) fn corrected_cluster_present(
    pa: &crate::parser::PresenceAbsence,
    ci: ClusterIdx,
    gi: GenomeIdx,
    matrix_rescue_map: &HashMap<(ClusterIdx, GenomeIdx), &RescueEvent>,
    nsc_presence: &ahash::HashSet<(usize, usize)>,
) -> bool {
    pa.is_present(ci, gi)
        || matrix_rescue_map.contains_key(&(ci, gi))
        || nsc_presence.contains(&(ci.val(), gi.val()))
}

pub(super) fn mobile_kind_from_labels(labels: &[&str]) -> Option<(&'static str, String)> {
    let joined = labels.join("|");
    let normalized = joined
        .to_ascii_lowercase()
        .replace("%20", " ")
        .replace('_', " ");
    let kind = if normalized.contains("is6110") {
        "IS6110"
    } else if normalized.contains("is1081") {
        "IS1081"
    } else if normalized.contains("insertion sequence")
        || normalized.contains("transposase")
        || normalized.contains("is element")
    {
        "insertion_sequence"
    } else if normalized.contains("integrase") || normalized.contains("phage") {
        "phage_or_integrase"
    } else if normalized.contains("recombinase") {
        "recombinase"
    } else {
        return None;
    };
    Some((kind, joined))
}

pub(super) struct MobileElementSearchInput<'a> {
    pub(super) path: &'a Path,
    pub(super) genome_files: &'a HashMap<String, PathBuf>,
    pub(super) pa: &'a crate::parser::PresenceAbsence,
    pub(super) rescues: &'a [RescueEvent],
    pub(super) synteny_rejected: &'a [RescueEvent],
    pub(super) paralog_rejected: &'a [RescueEvent],
    pub(super) threads: usize,
    pub(super) min_identity: f64,
    pub(super) min_coverage: f64,
    pub(super) preset: &'a str,
}

pub(super) fn search_mobile_element_reference(
    input: MobileElementSearchInput<'_>,
) -> Result<Vec<MobileElementHit>> {
    let mobile_ref = parse_pan_reference(input.path, true)
        .with_context(|| format!("Cannot parse mobile-element FASTA {}", input.path.display()))?;
    if mobile_ref.is_empty() {
        return Ok(Vec::new());
    }

    let mut mobile_names: Vec<String> = mobile_ref.keys().cloned().collect();
    mobile_names.sort();
    let mobile_indices: Vec<ClusterIdx> = (0..mobile_names.len()).map(ClusterIdx).collect();

    let mut genomes: ahash::HashSet<GenomeIdx> = ahash::HashSet::default();
    for rescue in input
        .rescues
        .iter()
        .chain(input.synteny_rejected.iter())
        .chain(input.paralog_rejected.iter())
    {
        genomes.insert(rescue.genome_idx);
    }
    if genomes.is_empty() {
        return Ok(Vec::new());
    }
    let mut searchable_genomes: Vec<(GenomeIdx, Vec<ClusterIdx>)> = genomes
        .into_iter()
        .map(|gi| (gi, mobile_indices.clone()))
        .collect();
    searchable_genomes.sort_by_key(|(gi, _)| gi.val());

    info!(
        "Running mobile-element sequence audit: {} queries × {} genomes",
        mobile_names.len(),
        searchable_genomes.len()
    );
    let hits = batch_search(
        &mobile_ref,
        input.genome_files,
        &searchable_genomes,
        input.pa.genome_names(),
        &mobile_names,
        input.threads,
        input.min_identity,
        input.preset,
        true,
        5000,
    )?;

    let mut mobile_hits = Vec::new();
    for ((gi, ci), cluster_hits) in hits {
        for hit in cluster_hits {
            if hit.query_coverage() < input.min_coverage {
                continue;
            }
            let identity = hit.identity();
            let coverage = hit.query_coverage();
            mobile_hits.push(MobileElementHit {
                genome_idx: gi,
                query_name: mobile_names[ci.val()].clone(),
                contig: hit.target,
                start: hit.target_start,
                end: hit.target_end,
                strand: hit.strand,
                identity,
                coverage,
            });
        }
    }
    mobile_hits.sort_by(|a, b| {
        a.genome_idx
            .cmp(&b.genome_idx)
            .then(a.contig.cmp(&b.contig))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
            .then(a.query_name.cmp(&b.query_name))
    });
    info!(
        "Mobile-element sequence audit retained {} hits",
        mobile_hits.len()
    );
    Ok(mobile_hits)
}

pub(super) fn interval_distance(start_a: u64, end_a: u64, start_b: u64, end_b: u64) -> u64 {
    if end_a < start_b {
        start_b - end_a
    } else {
        start_a.saturating_sub(end_b)
    }
}

pub(super) fn nearest_mobile_element_for_rescue(
    rescue: &RescueEvent,
    pa: &crate::parser::PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    mobile_hits: &[MobileElementHit],
    mobile_window: u64,
) -> Option<MobileElementEvidence> {
    let mut best: Option<MobileElementEvidence> = None;

    for hit in mobile_hits {
        if hit.genome_idx != rescue.genome_idx || hit.contig != rescue.contig {
            continue;
        }
        let distance = interval_distance(rescue.start, rescue.end, hit.start, hit.end);
        if distance > mobile_window {
            continue;
        }
        let candidate = MobileElementEvidence {
            genome_name: pa.genome_name(rescue.genome_idx).to_string(),
            source: "sequence",
            kind: "mobile_element_sequence",
            label: format!(
                "{}:{}-{}({};id={:.1};cov={:.1})",
                hit.query_name,
                hit.start + 1,
                hit.end,
                hit.strand,
                hit.identity,
                hit.coverage
            ),
            distance_bp: distance,
            identity: Some(hit.identity),
            coverage: Some(hit.coverage),
        };
        if best
            .as_ref()
            .is_none_or(|current| mobile_evidence_is_better(&candidate, current))
        {
            best = Some(candidate);
        }
    }

    let Some(gff) = gff_indices.get(&rescue.genome_idx) else {
        return best;
    };
    for gene in &gff.genes {
        if gene.contig != rescue.contig {
            continue;
        }
        let distance = interval_distance(rescue.start, rescue.end, gene.start, gene.end);
        if distance > mobile_window {
            continue;
        }
        let labels = gff.labels_for_gene(gene);
        let Some((kind, label)) = mobile_kind_from_labels(&labels) else {
            continue;
        };
        let candidate = MobileElementEvidence {
            genome_name: pa.genome_name(rescue.genome_idx).to_string(),
            source: "annotation",
            kind,
            label,
            distance_bp: distance,
            identity: None,
            coverage: None,
        };
        if best
            .as_ref()
            .is_none_or(|current| mobile_evidence_is_better(&candidate, current))
        {
            best = Some(candidate);
        }
    }
    best
}

pub(super) fn nearest_mobile_element_for_cluster(
    rescues: &[&RescueEvent],
    pa: &crate::parser::PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    mobile_hits: &[MobileElementHit],
    mobile_window: u64,
) -> Option<MobileElementEvidence> {
    rescues
        .iter()
        .filter_map(|rescue| {
            nearest_mobile_element_for_rescue(rescue, pa, gff_indices, mobile_hits, mobile_window)
        })
        .min_by(|a, b| {
            a.distance_bp
                .cmp(&b.distance_bp)
                .then(mobile_source_priority(a.source).cmp(&mobile_source_priority(b.source)))
                .then(a.genome_name.cmp(&b.genome_name))
                .then(a.label.cmp(&b.label))
        })
}

pub(super) fn mobile_source_priority(source: &str) -> u8 {
    if source == "sequence" {
        0
    } else {
        1
    }
}

pub(super) fn mobile_evidence_is_better(
    candidate: &MobileElementEvidence,
    current: &MobileElementEvidence,
) -> bool {
    candidate
        .distance_bp
        .cmp(&current.distance_bp)
        .then(mobile_source_priority(candidate.source).cmp(&mobile_source_priority(current.source)))
        .then(candidate.label.cmp(&current.label))
        .is_lt()
}

pub(super) fn nsc_role_for_cluster(
    ci: ClusterIdx,
    nsc_result: &Option<nscluster::NsClusterResult>,
) -> String {
    let Some(nsc) = nsc_result else {
        return String::new();
    };
    if nsc.merge_map.contains_key(&ci) {
        return "nsc_member_merged_into_representative".to_string();
    }
    for cluster in &nsc.merged_clusters {
        if cluster.members.contains(&ci) {
            return "nsc_representative".to_string();
        }
    }
    String::new()
}

pub(super) fn rescue_type_counts(rescues: &[&RescueEvent]) -> String {
    if rescues.is_empty() {
        return String::new();
    }
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for rescue in rescues {
        *counts.entry(rescue.rescue_type.to_string()).or_default() += 1;
    }
    counts
        .into_iter()
        .map(|(kind, count)| format!("{}:{}", kind, count))
        .collect::<Vec<_>>()
        .join(";")
}

pub(super) fn orf_signal_counts(rescues: &[&RescueEvent]) -> String {
    if rescues.is_empty() {
        return String::new();
    }
    let mut counts: BTreeMap<&'static str, usize> = BTreeMap::new();
    for rescue in rescues {
        *counts.entry(rescue_orf_signal(rescue)).or_default() += 1;
    }
    counts
        .into_iter()
        .map(|(signal, count)| format!("{}:{}", signal, count))
        .collect::<Vec<_>>()
        .join(";")
}

pub(super) fn synteny_mean_score(rescues: &[&RescueEvent]) -> String {
    let scores: Vec<f64> = rescues
        .iter()
        .filter_map(|rescue| rescue.synteny_evidence.score)
        .collect();
    if scores.is_empty() {
        String::new()
    } else {
        format!("{:.3}", scores.iter().sum::<f64>() / scores.len() as f64)
    }
}

pub(super) fn synteny_pass_fraction(rescues: &[&RescueEvent]) -> String {
    let assessed = rescues
        .iter()
        .filter(|rescue| rescue.synteny_ok.is_some())
        .count();
    if assessed == 0 {
        return String::new();
    }
    let passed = rescues
        .iter()
        .filter(|rescue| rescue.synteny_ok == Some(true))
        .count();
    format!("{:.3}", passed as f64 / assessed as f64)
}

pub(super) fn paralog_locus_risk(
    context: Option<&ClusterContextSummary>,
    accepted: &[&RescueEvent],
    synteny_rejected: &[&RescueEvent],
    paralog_rejected: &[&RescueEvent],
) -> &'static str {
    let has_paralog = accepted
        .iter()
        .any(|rescue| matches!(rescue.rescue_type, RescueType::Paralog { .. }))
        || !paralog_rejected.is_empty();
    if has_paralog || context.is_some_and(|c| c.interpretation == "possible_multilocus_or_paralog")
    {
        "high"
    } else if !synteny_rejected.is_empty()
        || context.is_some_and(|c| c.interpretation == "mostly_stable_review_minor_contexts")
    {
        "moderate"
    } else if context.is_some_and(|c| {
        matches!(
            c.interpretation,
            "stable_single_locus" | "stable_single_locus_minor_noise"
        )
    }) {
        "low"
    } else {
        "unknown"
    }
}

pub(super) fn phylogenetic_reliability(phylo: &PhyloPatternSummary) -> &'static str {
    if !phylo.assessed {
        return "not_assessed";
    }
    if matches!(
        phylo.pattern,
        "singleton_or_noise" | "homoplastic_gain_loss" | "dispersed_gain_loss"
    ) {
        return "low";
    }
    if phylo.parsimony_changes >= 6 {
        return "low";
    }
    if matches!(
        phylo.pattern,
        "lineage_specific_gain"
            | "clade_associated_presence"
            | "ancestral_retained"
            | "ancestral_loss"
    ) && phylo.parsimony_changes <= 2
    {
        return "high";
    }
    "moderate"
}

pub(super) struct MechanismInput<'a> {
    pub(super) pa: &'a crate::parser::PresenceAbsence,
    pub(super) phylo: &'a PhyloPatternSummary,
    pub(super) reference_state: &'a str,
    pub(super) reference_event: Option<&'a ReferenceClusterEvent>,
    pub(super) accepted: &'a [&'a RescueEvent],
    pub(super) synteny_rejected: &'a [&'a RescueEvent],
    pub(super) paralog_rejected: &'a [&'a RescueEvent],
    pub(super) context: Option<&'a ClusterContextSummary>,
    pub(super) mobile: Option<MobileElementEvidence>,
    pub(super) nsc_role: &'a str,
}

pub(super) fn classify_mechanism(input: MechanismInput<'_>) -> MechanismCall {
    let MechanismInput {
        pa,
        phylo,
        reference_state,
        reference_event,
        accepted,
        synteny_rejected,
        paralog_rejected,
        context,
        mobile,
        nsc_role,
    } = input;
    let has_pseudogene = accepted
        .iter()
        .any(|r| matches!(r.rescue_type, RescueType::Pseudogene { .. }));
    let has_fragment = accepted.iter().any(|r| {
        matches!(
            r.rescue_type,
            RescueType::Fragmented { .. } | RescueType::AlignmentSplit { .. }
        )
    });
    let has_split_or_boundary = accepted.iter().any(|r| {
        matches!(
            r.rescue_type,
            RescueType::Split { .. } | RescueType::StartVariant { .. }
        )
    });
    let has_paralog = accepted
        .iter()
        .any(|r| matches!(r.rescue_type, RescueType::Paralog { .. }))
        || !paralog_rejected.is_empty();
    let unstable_context =
        context.is_some_and(|c| c.interpretation == "possible_multilocus_or_paralog");
    let mobile_is = mobile.as_ref().is_some_and(|m| {
        (m.source == "sequence" || matches!(m.kind, "IS6110" | "IS1081" | "insertion_sequence"))
            && m.distance_bp <= 250
    });

    if let Some(event) = reference_event {
        if (has_paralog || unstable_context)
            && matches!(
                event.event_class,
                "partial_reference_homology" | "same_gene_full_length"
            )
        {
            return MechanismCall {
                mechanism: "duplication_or_paralog_expansion",
                detail: format!(
                    "reference_event={};target_reference_cluster={};paralog_or_multilocus_context",
                    event.event_class,
                    pa.cluster_name(event.target_cluster_idx)
                ),
                support: "low",
                wording: "Reference-like sequence is present in a paralogous or unstable context; interpret as copy-number/paralog or cluster-definition complexity before claiming gain/loss.",
                mobile,
            };
        }
        let mechanism = match event.event_class {
            "short_in_frame_allele" | "long_in_frame_allele" | "same_gene_full_length" => {
                "reference_gene_allelic_derivative"
            }
            "frameshift_or_truncated_reference_homology"
            | "frameshift_or_extension_reference_homology" => {
                "reference_gene_disruption_or_truncation"
            }
            _ => "partial_reference_homology",
        };
        return MechanismCall {
            mechanism,
            detail: format!(
                "{};target_reference_cluster={};target_reference_gene={}",
                event.event_class,
                pa.cluster_name(event.target_cluster_idx),
                event.target_gene_id
            ),
            support: if event.identity >= 95.0 { "high" } else { "moderate" },
            wording: "Keep as a separate matrix cluster, but describe as an event-derived reference-gene variant rather than a clean novel gene.",
            mobile,
        };
    }

    if mobile_is && (has_pseudogene || has_fragment) {
        return MechanismCall {
            mechanism: "IS_interruption_supported",
            detail: "fragmented_or_disrupted_rescue_with_nearby_insertion_sequence".to_string(),
            support: "high",
            wording: "Interpret as an IS-supported gene interruption/disruption candidate; keep matrix action governed by the orthology evidence.",
            mobile,
        };
    }

    if has_pseudogene {
        return MechanismCall {
            mechanism: "pseudogenization",
            detail: "pseudogene_rescue_evidence".to_string(),
            support: "moderate",
            wording: "Describe as disrupted homologous sequence, not as intact gene presence.",
            mobile,
        };
    }

    if has_fragment {
        return MechanismCall {
            mechanism: "gene_fragmentation_or_structural_break",
            detail: "fragmented_or_alignment_split_rescue_evidence".to_string(),
            support: "moderate",
            wording: "Describe as partial/fragmented homology and review the local structure before making a gene gain/loss claim.",
            mobile,
        };
    }

    if has_paralog || unstable_context {
        return MechanismCall {
            mechanism: "duplication_or_paralog_expansion",
            detail: "paralog_evidence_or_multilocus_context".to_string(),
            support: if has_paralog { "moderate" } else { "low" },
            wording: "Treat as a copy-number/paralog or multi-locus signal; avoid interpreting as a single clean gain/loss without locus review.",
            mobile,
        };
    }

    if !synteny_rejected.is_empty() {
        return MechanismCall {
            mechanism: "translocation_or_context_shift_candidate",
            detail: "sequence_evidence_rejected_by_synteny".to_string(),
            support: "low",
            wording: "Sequence similarity exists in an incompatible context; review as possible translocation, paralog, or annotation artefact.",
            mobile,
        };
    }

    if has_split_or_boundary || !nsc_role.is_empty() {
        return MechanismCall {
            mechanism: "annotation_or_clustering_artifact",
            detail: if nsc_role.is_empty() {
                "split_or_start_boundary_evidence".to_string()
            } else {
                nsc_role.to_string()
            },
            support: "moderate",
            wording: "Treat as a cluster split, boundary shift, or redundant annotation signal rather than an independent evolutionary gene event.",
            mobile,
        };
    }

    if reference_state == "absent"
        && matches!(
            phylo.pattern,
            "lineage_specific_gain"
                | "clade_associated_presence"
                | "single_gain_with_incomplete_clade_support"
        )
    {
        return MechanismCall {
            mechanism: "clean_insertion_or_acquisition_candidate",
            detail: "reference_absent_clade_associated_presence".to_string(),
            support: if phylo.carrier_monophyly { "moderate" } else { "low" },
            wording: "Candidate accessory acquisition/insertion relative to the audited reference; validate flanking sequence before claiming a novel locus.",
            mobile,
        };
    }

    if reference_state == "absent" && phylo.root_state == "present" {
        return MechanismCall {
            mechanism: "reference_or_lineage_loss_candidate",
            detail: "inferred_ancestral_presence_with_reference_absence".to_string(),
            support: "moderate",
            wording: "Candidate loss in the reference/background lineage rather than a de novo gain in all other genomes.",
            mobile,
        };
    }

    if matches!(phylo.pattern, "ancestral_loss" | "recurrent_loss") {
        return MechanismCall {
            mechanism: "deletion_loss_candidate",
            detail: phylo.pattern.to_string(),
            support: "moderate",
            wording: "Presence is inferred ancestral with one or more losses; inspect flanks for deletion scars or mobile-element mediation.",
            mobile,
        };
    }

    if phylo.pattern == "homoplastic_gain_loss" {
        return MechanismCall {
            mechanism: "unstable_or_homoplastic_accessory_signal",
            detail: "many_parsimony_changes".to_string(),
            support: "review",
            wording: "Highly recurrent presence/absence can reflect real mobility, recombination, or cluster/annotation instability; review before highlighting.",
            mobile,
        };
    }

    if phylo.pattern == "ancestral_retained" {
        return MechanismCall {
            mechanism: "ancestral_retained_no_disruption",
            detail: "present_across_matched_tree".to_string(),
            support: "high",
            wording: "Consistent with conserved ancestral presence in the analysed tree.",
            mobile,
        };
    }

    MechanismCall {
        mechanism: "unresolved_accessory_variation",
        detail: phylo.pattern.to_string(),
        support: "review",
        wording: "Evolutionary process is not resolved from matrix, synteny, reference-event, and optional tree evidence alone.",
        mobile,
    }
}

pub(super) fn write_evolutionary_events(inputs: EvolutionaryEventInputs<'_>) -> Result<()> {
    let nsc_presence = nsc_presence_updates(inputs.pa, inputs.nsc_result, inputs.matrix_rescue_map);
    let reference_gi = inputs
        .reference_genome
        .and_then(|name| resolve_reference_genome(inputs.pa, name));

    let mut accepted_by_cluster: Vec<Vec<&RescueEvent>> = vec![Vec::new(); inputs.pa.n_clusters()];
    for rescue in inputs.rescues {
        accepted_by_cluster[rescue.cluster_idx.val()].push(rescue);
    }
    let mut synteny_by_cluster: Vec<Vec<&RescueEvent>> = vec![Vec::new(); inputs.pa.n_clusters()];
    for rescue in inputs.synteny_rejected {
        synteny_by_cluster[rescue.cluster_idx.val()].push(rescue);
    }
    let mut paralog_by_cluster: Vec<Vec<&RescueEvent>> = vec![Vec::new(); inputs.pa.n_clusters()];
    for rescue in inputs.paralog_rejected {
        paralog_by_cluster[rescue.cluster_idx.val()].push(rescue);
    }
    let mut reference_events_by_cluster: Vec<Option<&ReferenceClusterEvent>> =
        vec![None; inputs.pa.n_clusters()];
    for event in inputs.reference_cluster_events {
        reference_events_by_cluster[event.query_cluster_idx.val()] = Some(event);
    }
    let mut contexts_by_cluster: Vec<Option<&ClusterContextSummary>> =
        vec![None; inputs.pa.n_clusters()];
    for context in inputs.context_summaries {
        contexts_by_cluster[context.cluster_idx.val()] = Some(context);
    }

    let path = inputs.output_dir.join("evolutionary_events.tsv");
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "gene_cluster\tn_present_corrected\tprevalence_corrected\ttree_matched_tips\treference_state\tphylogenetic_pattern\tphylogenetic_reliability\tancestral_state\tphylogenetic_changes\tphylogenetic_gains\tphylogenetic_losses\tcarrier_monophyly\tmax_carrier_clade_fraction\tmolecular_mechanism\tmechanism_detail\tsupport_level\tis_interruption_supported\tmobile_element_source\tmobile_element_kind\tmobile_element_label\tmobile_element_genome\tmobile_element_distance_bp\tmobile_element_identity\tmobile_element_coverage\treference_event_class\ttarget_reference_cluster\ttarget_reference_gene\tcontext_interpretation\tparalog_locus_risk\tsynteny_mean_score\tsynteny_pass_fraction\torf_integrity_summary\trescue_evidence_counts\trecommended_wording"
    )?;

    for ci_val in 0..inputs.pa.n_clusters() {
        let ci = ClusterIdx(ci_val);
        let present_by_genome: Vec<bool> = (0..inputs.pa.n_genomes())
            .map(|gi_val| {
                corrected_cluster_present(
                    inputs.pa,
                    ci,
                    GenomeIdx(gi_val),
                    inputs.matrix_rescue_map,
                    &nsc_presence,
                )
            })
            .collect();
        let n_present = present_by_genome.iter().filter(|present| **present).count();
        let prevalence = if inputs.pa.n_genomes() > 0 {
            n_present as f64 / inputs.pa.n_genomes() as f64
        } else {
            0.0
        };
        let reference_state = reference_gi
            .map(|gi| {
                if corrected_cluster_present(
                    inputs.pa,
                    ci,
                    gi,
                    inputs.matrix_rescue_map,
                    &nsc_presence,
                ) {
                    "present"
                } else {
                    "absent"
                }
            })
            .unwrap_or("not_assessed");
        let phylo = summarize_phylogenetic_pattern(inputs.phylo_tree, &present_by_genome);
        let accepted = &accepted_by_cluster[ci_val];
        let synteny_rejected = &synteny_by_cluster[ci_val];
        let paralog_rejected = &paralog_by_cluster[ci_val];
        let all_rescue_evidence: Vec<&RescueEvent> = accepted
            .iter()
            .copied()
            .chain(synteny_rejected.iter().copied())
            .chain(paralog_rejected.iter().copied())
            .collect();
        let mobile = nearest_mobile_element_for_cluster(
            &all_rescue_evidence,
            inputs.pa,
            inputs.gff_indices,
            inputs.mobile_element_hits,
            inputs.mobile_element_window,
        );
        let nsc_role = nsc_role_for_cluster(ci, inputs.nsc_result);
        let mechanism = classify_mechanism(MechanismInput {
            pa: inputs.pa,
            phylo: &phylo,
            reference_state,
            reference_event: reference_events_by_cluster[ci_val],
            accepted,
            synteny_rejected,
            paralog_rejected,
            context: contexts_by_cluster[ci_val],
            mobile,
            nsc_role: &nsc_role,
        });
        let is_interruption = mechanism.mechanism == "IS_interruption_supported";
        let (
            mobile_source,
            mobile_kind,
            mobile_label,
            mobile_genome,
            mobile_distance,
            mobile_identity,
            mobile_coverage,
        ) = if let Some(mobile) = mechanism.mobile.as_ref() {
            (
                mobile.source.to_string(),
                mobile.kind.to_string(),
                mobile.label.clone(),
                mobile.genome_name.clone(),
                mobile.distance_bp.to_string(),
                mobile
                    .identity
                    .map(|value| format!("{:.1}", value))
                    .unwrap_or_default(),
                mobile
                    .coverage
                    .map(|value| format!("{:.1}", value))
                    .unwrap_or_default(),
            )
        } else {
            (
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
            )
        };
        let reference_event = reference_events_by_cluster[ci_val];
        let context_interpretation = contexts_by_cluster[ci_val]
            .map(|context| context.interpretation)
            .unwrap_or("");
        let paralog_risk = paralog_locus_risk(
            contexts_by_cluster[ci_val],
            accepted,
            synteny_rejected,
            paralog_rejected,
        );
        write_tsv_fields(
            &mut w,
            &[
                inputs.pa.cluster_name(ci).to_string(),
                n_present.to_string(),
                format!("{:.4}", prevalence),
                if phylo.assessed {
                    phylo.matched_tips.to_string()
                } else {
                    String::new()
                },
                reference_state.to_string(),
                phylo.pattern.to_string(),
                phylogenetic_reliability(&phylo).to_string(),
                phylo.root_state.to_string(),
                if phylo.assessed {
                    phylo.parsimony_changes.to_string()
                } else {
                    String::new()
                },
                if phylo.assessed {
                    phylo.gains.to_string()
                } else {
                    String::new()
                },
                if phylo.assessed {
                    phylo.losses.to_string()
                } else {
                    String::new()
                },
                if phylo.assessed {
                    phylo.carrier_monophyly.to_string()
                } else {
                    String::new()
                },
                if phylo.assessed {
                    format!("{:.3}", phylo.max_carrier_clade_fraction)
                } else {
                    String::new()
                },
                mechanism.mechanism.to_string(),
                mechanism.detail,
                mechanism.support.to_string(),
                is_interruption.to_string(),
                mobile_source,
                mobile_kind,
                mobile_label,
                mobile_genome,
                mobile_distance,
                mobile_identity,
                mobile_coverage,
                reference_event
                    .map(|event| event.event_class.to_string())
                    .unwrap_or_default(),
                reference_event
                    .map(|event| inputs.pa.cluster_name(event.target_cluster_idx).to_string())
                    .unwrap_or_default(),
                reference_event
                    .map(|event| event.target_gene_id.clone())
                    .unwrap_or_default(),
                context_interpretation.to_string(),
                paralog_risk.to_string(),
                synteny_mean_score(&all_rescue_evidence),
                synteny_pass_fraction(&all_rescue_evidence),
                orf_signal_counts(&all_rescue_evidence),
                rescue_type_counts(accepted),
                mechanism.wording.to_string(),
            ],
        )?;
    }
    w.flush()?;
    info!("Evolutionary event audit written to {}", path.display());
    Ok(())
}

#[derive(Debug, Clone, Default)]
pub(super) struct ContextAccumulator {
    pub(super) count: usize,
    pub(super) examples: Vec<String>,
}

#[derive(Debug, Clone)]
pub(super) struct ClusterContextSummary {
    pub(super) cluster_idx: ClusterIdx,
    pub(super) n_present: usize,
    pub(super) n_mapped: usize,
    pub(super) n_unresolved: usize,
    pub(super) n_unique_contexts: usize,
    pub(super) dominant_count: usize,
    pub(super) dominant_fraction: f64,
    pub(super) interpretation: &'static str,
}

pub(super) fn context_signature(
    gff: &GffIndex,
    genes: &[&GffGene],
    gene_to_cluster: &HashMap<String, ClusterIdx>,
    pa: &crate::parser::PresenceAbsence,
) -> String {
    let names: Vec<&str> = genes
        .iter()
        .filter_map(|gene| gff.cluster_for_gene(gene, gene_to_cluster))
        .map(|ci| pa.cluster_name(ci))
        .collect();
    if names.is_empty() {
        "NA".to_string()
    } else {
        names.join("|")
    }
}

pub(super) fn context_interpretation(
    n_mapped: usize,
    n_unique_contexts: usize,
    dominant_fraction: f64,
) -> &'static str {
    if n_mapped < 3 {
        "insufficient_context"
    } else if n_unique_contexts <= 1 {
        "stable_single_locus"
    } else if dominant_fraction >= 0.95 {
        "stable_single_locus_minor_noise"
    } else if dominant_fraction >= 0.80 {
        "mostly_stable_review_minor_contexts"
    } else {
        "possible_multilocus_or_paralog"
    }
}

pub(super) fn write_cluster_locus_contexts(
    pa: &crate::parser::PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    gene_to_cluster_maps: &HashMap<GenomeIdx, HashMap<String, ClusterIdx>>,
    output_dir: &std::path::Path,
) -> Result<(usize, Vec<ClusterContextSummary>)> {
    let mut context_rows: BTreeMap<(usize, String, String), ContextAccumulator> = BTreeMap::new();
    let mut summaries = Vec::new();

    for ci_val in 0..pa.n_clusters() {
        let ci = ClusterIdx(ci_val);
        let mut n_mapped = 0usize;
        let mut n_unresolved = 0usize;
        let mut context_counts: HashMap<(String, String), usize> = HashMap::new();

        let mut genomes: Vec<GenomeIdx> = pa.present_genomes(ci).iter().copied().collect();
        genomes.sort_by_key(|gi| gi.val());

        for gi in genomes {
            let Some(gff) = gff_indices.get(&gi) else {
                continue;
            };
            let Some(g2c) = gene_to_cluster_maps.get(&gi) else {
                continue;
            };
            let cell = pa.cell(ci, gi);
            let Some(gene) = gff_gene_for_cell(gff, cell) else {
                n_unresolved += 1;
                continue;
            };

            let (upstream, downstream) = gff.context_window(
                &gene.contig,
                gene.start,
                gene.end,
                Some(gene.gene_id.as_str()),
                2,
            );
            let upstream_sig = context_signature(gff, &upstream, g2c, pa);
            let downstream_sig = context_signature(gff, &downstream, g2c, pa);
            let key = (upstream_sig.clone(), downstream_sig.clone());
            *context_counts.entry(key).or_insert(0) += 1;
            n_mapped += 1;

            let row_key = (ci.val(), upstream_sig, downstream_sig);
            let row = context_rows.entry(row_key).or_default();
            row.count += 1;
            if row.examples.len() < 5 {
                row.examples
                    .push(format!("{}:{}", pa.genome_name(gi), gene.gene_id));
            }
        }

        let n_unique_contexts = context_counts.len();
        let dominant_count = context_counts.values().copied().max().unwrap_or(0);
        let dominant_fraction = if n_mapped > 0 {
            dominant_count as f64 / n_mapped as f64
        } else {
            0.0
        };
        let interpretation = context_interpretation(n_mapped, n_unique_contexts, dominant_fraction);
        summaries.push(ClusterContextSummary {
            cluster_idx: ci,
            n_present: pa.present_count(ci),
            n_mapped,
            n_unresolved,
            n_unique_contexts,
            dominant_count,
            dominant_fraction,
            interpretation,
        });
    }

    let context_path = output_dir.join("cluster_locus_context.tsv");
    let context_file = fs::File::create(&context_path)?;
    let mut context_w = BufWriter::new(context_file);
    writeln!(
        context_w,
        "gene_cluster\tcontext_rank\tn_genomes\tfraction_of_mapped_contexts\tupstream_context\tdownstream_context\texample_genomes"
    )?;
    let mut written_context_rows = 0usize;
    let mut rows_by_cluster: BTreeMap<usize, Vec<(String, String, ContextAccumulator)>> =
        BTreeMap::new();
    for ((ci, upstream, downstream), acc) in context_rows {
        rows_by_cluster
            .entry(ci)
            .or_default()
            .push((upstream, downstream, acc));
    }
    for (ci_val, mut rows) in rows_by_cluster {
        rows.sort_by(|a, b| {
            b.2.count
                .cmp(&a.2.count)
                .then(a.0.cmp(&b.0))
                .then(a.1.cmp(&b.1))
        });
        let mapped = summaries[ci_val].n_mapped.max(1);
        for (rank, (upstream, downstream, acc)) in rows.into_iter().enumerate() {
            write_tsv_fields(
                &mut context_w,
                &[
                    pa.cluster_name(ClusterIdx(ci_val)).to_string(),
                    (rank + 1).to_string(),
                    acc.count.to_string(),
                    format!("{:.3}", acc.count as f64 / mapped as f64),
                    upstream,
                    downstream,
                    acc.examples.join(";"),
                ],
            )?;
            written_context_rows += 1;
        }
    }
    context_w.flush()?;

    let summary_path = output_dir.join("cluster_locus_summary.tsv");
    let summary_file = fs::File::create(&summary_path)?;
    let mut summary_w = BufWriter::new(summary_file);
    writeln!(
        summary_w,
        "gene_cluster\tn_present\tn_context_mapped\tn_unresolved_gene_ids\tn_unique_contexts\tdominant_context_count\tdominant_context_fraction\tinterpretation"
    )?;
    for summary in &summaries {
        write_tsv_fields(
            &mut summary_w,
            &[
                pa.cluster_name(summary.cluster_idx).to_string(),
                summary.n_present.to_string(),
                summary.n_mapped.to_string(),
                summary.n_unresolved.to_string(),
                summary.n_unique_contexts.to_string(),
                summary.dominant_count.to_string(),
                format!("{:.3}", summary.dominant_fraction),
                summary.interpretation.to_string(),
            ],
        )?;
    }
    summary_w.flush()?;

    Ok((written_context_rows, summaries))
}
