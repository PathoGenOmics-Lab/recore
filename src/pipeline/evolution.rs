use super::reference_audit::{display_reference_gene_id, ReferenceClusterEvent};
use super::*;
use crate::relationships::ClusterComponentInfo;

#[derive(Debug, Clone)]
pub(super) struct PhyloNode {
    pub(super) label: Option<String>,
    pub(super) genome_idx: Option<GenomeIdx>,
    pub(super) branch_length: Option<f64>,
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
    pub(super) root_state_ambiguous: bool,
    pub(super) carrier_monophyly: bool,
    pub(super) max_carrier_clade_fraction: f64,
    pub(super) carrier_mrca_tips: usize,
    pub(super) carrier_mrca_purity: f64,
    pub(super) dollo_loss_events_in_mrca: usize,
    pub(super) pattern: &'static str,
}

#[derive(Debug, Clone, Default)]
pub(super) struct PenalizedPhyloSummary {
    pub(super) model: &'static str,
    pub(super) dollo_loss_count: usize,
    pub(super) dollo_violation_count: usize,
    pub(super) dollo_origin_tips: usize,
    pub(super) dollo_origin_purity: f64,
    pub(super) score: f64,
}

#[derive(Debug, Clone, Default)]
pub(super) struct FitchNodeSummary {
    pub(super) state_set: u8,
    pub(super) matched_tips: usize,
    pub(super) present_tips: usize,
    pub(super) changes: usize,
    pub(super) children: Vec<FitchNodeSummary>,
    /// Fewest changes in this subtree with the node absent ([0]) or present
    /// ([1]).
    pub(super) cost: [usize; 2],
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
            let branch_length = self.parse_branch_length()?;
            Ok(PhyloNode {
                label,
                genome_idx: None,
                branch_length,
                children,
            })
        } else {
            let label = self.parse_optional_label();
            if label.as_deref().unwrap_or("").is_empty() {
                anyhow::bail!("Expected Newick leaf label at byte {}", self.pos);
            }
            let branch_length = self.parse_branch_length()?;
            Ok(PhyloNode {
                label,
                genome_idx: None,
                branch_length,
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

    fn parse_branch_length(&mut self) -> Result<Option<f64>> {
        self.skip_ignored();
        if self.peek() != Some(':') {
            return Ok(None);
        }
        self.pos += 1;
        let start = self.pos;
        while let Some(ch) = self.peek() {
            if matches!(ch, ',' | ')' | ';' | '[') || ch.is_whitespace() {
                break;
            }
            self.pos += ch.len_utf8();
        }
        let raw = self.input[start..self.pos].trim();
        if raw.is_empty() {
            anyhow::bail!("Expected Newick branch length at byte {}", start);
        }
        let branch_length = raw
            .parse::<f64>()
            .with_context(|| format!("Invalid Newick branch length '{raw}'"))?;
        if !branch_length.is_finite() {
            anyhow::bail!("Invalid non-finite Newick branch length '{raw}'");
        }
        Ok(Some(branch_length))
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

/// Match tree tips to matrix genomes, one tip per genome.
///
/// A second tip that resolves to a genome another tip already took, by the
/// same name or through a shared AS identifier, used to be matched as well:
/// the genome was counted twice, the genome the tip really stood for was
/// hidden, and the coverage check, which compares counts, saw nothing wrong.
/// An exact name beats an AS-identifier match, and otherwise the first tip in
/// the tree keeps the genome; the others are left unmatched, by name, in a
/// warning.
pub(super) fn resolve_phylo_tips(
    node: &mut PhyloNode,
    pa: &crate::parser::PresenceAbsence,
) -> (usize, usize) {
    let mut candidates: Vec<(Option<(GenomeIdx, bool)>, String)> = Vec::new();
    collect_tip_candidates(node, pa, &mut candidates);

    let mut winner: HashMap<GenomeIdx, usize> = HashMap::new();
    for (i, (candidate, _)) in candidates.iter().enumerate() {
        let Some((gi, exact)) = candidate else {
            continue;
        };
        match winner.get(gi) {
            None => {
                winner.insert(*gi, i);
            }
            Some(&j) => {
                let held_exact = candidates[j].0.is_some_and(|(_, e)| e);
                if *exact && !held_exact {
                    winner.insert(*gi, i);
                }
            }
        }
    }

    let mut duplicates = Vec::new();
    let assigned: Vec<Option<GenomeIdx>> = candidates
        .iter()
        .enumerate()
        .map(|(i, (candidate, label))| {
            let (gi, _) = (*candidate)?;
            if winner[&gi] == i {
                Some(gi)
            } else {
                duplicates.push(format!("{} (as {})", label, pa.genome_name(gi)));
                None
            }
        })
        .collect();
    if !duplicates.is_empty() {
        warn!(
            "{} tree tips resolve to a genome another tip already matched, and are left unmatched: {}",
            duplicates.len(),
            duplicates.join(", ")
        );
    }

    let matched = assigned.iter().filter(|gi| gi.is_some()).count();
    let unmatched = assigned.len() - matched;
    let mut assigned = assigned.into_iter();
    assign_tip_genomes(node, &mut assigned);
    (matched, unmatched)
}

fn collect_tip_candidates(
    node: &PhyloNode,
    pa: &crate::parser::PresenceAbsence,
    out: &mut Vec<(Option<(GenomeIdx, bool)>, String)>,
) {
    if node.children.is_empty() {
        let label = node.label.clone().unwrap_or_default();
        let candidate = node
            .label
            .as_deref()
            .and_then(|label| resolve_tree_tip_with_kind(label, pa));
        out.push((candidate, label));
        return;
    }
    for child in &node.children {
        collect_tip_candidates(child, pa, out);
    }
}

fn assign_tip_genomes(
    node: &mut PhyloNode,
    assigned: &mut impl Iterator<Item = Option<GenomeIdx>>,
) {
    if node.children.is_empty() {
        node.genome_idx = assigned.next().flatten();
        return;
    }
    for child in &mut node.children {
        assign_tip_genomes(child, assigned);
    }
}

/// The genome a tip label stands for, and whether it matched by exact name
/// (true) or only through an embedded AS identifier (false).
fn resolve_tree_tip_with_kind(
    label: &str,
    pa: &crate::parser::PresenceAbsence,
) -> Option<(GenomeIdx, bool)> {
    resolve_tree_tip(label, pa).map(|gi| (gi, pa.genome_name(gi) == label))
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

/// Parsimony over presence and absence, as Sankoff costs per state.
///
/// It used to be Fitch's rule applied to children one pair at a time, which is
/// only right for two children: at a node with three or more, the result
/// depended on the order the siblings were written in the Newick file, and the
/// counts could contradict themselves (two changes reported as two gains and a
/// loss). Per-state costs summed over all children give the same answer for
/// any order and any number of children.
pub(super) fn fitch_postorder(node: &PhyloNode, present_by_genome: &[bool]) -> FitchNodeSummary {
    const ABSENT: u8 = 0b01;
    const PRESENT: u8 = 0b10;
    const UNREACHABLE: usize = usize::MAX / 4;

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
            cost: if present {
                [UNREACHABLE, 0]
            } else {
                [0, UNREACHABLE]
            },
        };
    }

    let mut children = Vec::new();
    let mut matched_tips = 0usize;
    let mut present_tips = 0usize;
    let mut cost = [0usize; 2];
    for child in &node.children {
        let summary = fitch_postorder(child, present_by_genome);
        if summary.matched_tips > 0 {
            matched_tips += summary.matched_tips;
            present_tips += summary.present_tips;
            for (state, total) in cost.iter_mut().enumerate() {
                *total += summary.cost[state].min(summary.cost[1 - state] + 1);
            }
        }
        children.push(summary);
    }
    if matched_tips == 0 {
        return FitchNodeSummary {
            children,
            ..FitchNodeSummary::default()
        };
    }

    let best = cost[0].min(cost[1]);
    let state_set =
        (if cost[0] == best { ABSENT } else { 0 }) | (if cost[1] == best { PRESENT } else { 0 });
    FitchNodeSummary {
        state_set,
        matched_tips,
        present_tips,
        changes: best,
        children,
        cost,
    }
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
    // Keep the parent's state whenever that costs no more than switching.
    let state = if let Some(parent) = parent_state {
        let (keep, switch) = (usize::from(parent), usize::from(!parent));
        if summary.cost[keep] <= summary.cost[switch] + 1 {
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

pub(super) fn dollo_loss_events_inside(summary: &FitchNodeSummary) -> usize {
    if summary.matched_tips == 0 || summary.present_tips == summary.matched_tips {
        return 0;
    }
    if summary.present_tips == 0 {
        return 1;
    }
    summary
        .children
        .iter()
        .filter(|child| child.matched_tips > 0)
        .map(|child| {
            if child.present_tips == 0 {
                1
            } else {
                dollo_loss_events_inside(child)
            }
        })
        .sum()
}

pub(super) fn carrier_mrca_stats(
    summary: &FitchNodeSummary,
    total_present: usize,
) -> (usize, f64, usize) {
    if total_present == 0 || summary.matched_tips == 0 {
        return (0, 0.0, 0);
    }

    fn visit(summary: &FitchNodeSummary, total_present: usize, best: &mut Option<(usize, usize)>) {
        if summary.matched_tips > 0 && summary.present_tips == total_present {
            let losses = dollo_loss_events_inside(summary);
            match best {
                Some((best_tips, best_losses))
                    if summary.matched_tips > *best_tips
                        || (summary.matched_tips == *best_tips && losses >= *best_losses) => {}
                _ => *best = Some((summary.matched_tips, losses)),
            }
        }
        for child in &summary.children {
            visit(child, total_present, best);
        }
    }

    let mut best = None;
    visit(summary, total_present, &mut best);
    let (mrca_tips, losses) =
        best.unwrap_or((summary.matched_tips, dollo_loss_events_inside(summary)));
    let purity = total_present as f64 / mrca_tips.max(1) as f64;
    (mrca_tips, purity, losses)
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
    if summary.root_state_ambiguous {
        if summary.carrier_mrca_purity >= 0.80 && summary.dollo_loss_events_in_mrca <= 3 {
            return "clade_associated_presence_with_losses";
        }
        return "ambiguous_gain_loss";
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
    let root_state_ambiguous = fitch.state_set == 0b11;
    let root_present = assign_fitch_states(&fitch, None, &mut gains, &mut losses).unwrap_or(false);
    let (carrier_monophyly, max_carrier_clade_fraction) =
        carrier_clade_stats(&fitch, fitch.present_tips);
    let (carrier_mrca_tips, carrier_mrca_purity, dollo_loss_events_in_mrca) =
        carrier_mrca_stats(&fitch, fitch.present_tips);
    let mut summary = PhyloPatternSummary {
        assessed: true,
        matched_tips: fitch.matched_tips,
        present_tips: fitch.present_tips,
        absent_tips: fitch.matched_tips.saturating_sub(fitch.present_tips),
        parsimony_changes: fitch.changes,
        gains,
        losses,
        root_state: if root_state_ambiguous {
            "ambiguous"
        } else if root_present {
            "present"
        } else {
            "absent"
        },
        root_state_ambiguous,
        carrier_monophyly,
        max_carrier_clade_fraction,
        carrier_mrca_tips,
        carrier_mrca_purity,
        dollo_loss_events_in_mrca,
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
    pub(super) component_info: &'a [Option<ClusterComponentInfo>],
    pub(super) coverage: &'a SearchCoverage,
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
    pub(super) search_cache: Option<&'a Path>,
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
        input.search_cache,
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
    for gene in gff.genes() {
        if gene.contig != rescue.contig {
            continue;
        }
        let distance = interval_distance(rescue.start, rescue.end, gene.start, gene.end);
        if distance > mobile_window {
            continue;
        }
        // Worked out from the gene's names and descriptions when the GFF was
        // read; see GffIndex::mobile_annotation.
        let Some((kind, label)) = gff.mobile_annotation(&gene) else {
            continue;
        };
        let label = label.to_string();
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
    // No carrier on the tree is no reading at all, not a moderate one: it used
    // to score above a cluster with no tree.
    if !phylo.assessed || phylo.pattern == "absent_from_tree" {
        return "not_assessed";
    }
    if matches!(
        phylo.pattern,
        "singleton_or_noise"
            | "homoplastic_gain_loss"
            | "dispersed_gain_loss"
            | "ambiguous_gain_loss"
    ) {
        return "low";
    }
    if phylo.parsimony_changes >= 6 {
        return "low";
    }
    if phylo.root_state_ambiguous && phylo.carrier_mrca_purity < 0.80 {
        return "low";
    }
    if matches!(
        phylo.pattern,
        "lineage_specific_gain"
            | "clade_associated_presence"
            | "clade_associated_presence_with_losses"
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
                | "clade_associated_presence_with_losses"
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

/// How much of a cluster's absence the search actually tested.
#[derive(Debug, Clone, Copy, Default)]
pub(super) struct AbsenceCoverage {
    /// Genomes where the cluster is absent from the corrected matrix.
    pub(super) absent: usize,
    /// Of those, the ones whose gap was never searched.
    pub(super) not_searched: usize,
    /// Why, when any were not searched: the cluster's own reason when it was
    /// left out of the search, else the genomes without a FASTA.
    pub(super) reason: Option<NotSearched>,
    /// The reference genome lacks the cluster and its gap was never searched.
    pub(super) reference_untested: bool,
}

/// Mechanisms whose call rests on the cluster being absent somewhere: a gain
/// needs the gene missing elsewhere, a loss needs it missing where it was lost,
/// and a homoplastic or unresolved pattern is made of absences too.
const ABSENCE_BASED_MECHANISMS: [&str; 5] = [
    "clean_insertion_or_acquisition_candidate",
    "reference_or_lineage_loss_candidate",
    "deletion_loss_candidate",
    "unstable_or_homoplastic_accessory_signal",
    "unresolved_accessory_variation",
];

fn downgraded_support(support: &'static str) -> &'static str {
    match support {
        "high" => "moderate",
        "moderate" => "low",
        other => other,
    }
}

/// Qualify a mechanism call that rests on absences Recore never searched.
///
/// Gaps outside the search keep the input matrix's call, so a gain inferred
/// from a cluster missing everywhere else, or a loss inferred where it is
/// missing, is only as good as that call when the absences behind it were not
/// tested. Every such call is tagged in `mechanism_detail` and explained in the
/// returned sentence for `recommended_wording`; when none of its absences were
/// searched (for a reference loss, when the reference genome's was not), its
/// support also drops one level, which lowers its score and priority with it.
/// Under the default `--min-prevalence` this is every gain call on a cluster in
/// fewer than one genome in ten, and it used to read as an audited result.
pub(super) fn qualify_unsearched_absences(
    mechanism: &mut MechanismCall,
    coverage: &AbsenceCoverage,
) -> Option<String> {
    if !ABSENCE_BASED_MECHANISMS.contains(&mechanism.mechanism) {
        return None;
    }
    let reason = coverage.reason.filter(|_| coverage.not_searched > 0)?;
    let one = coverage.not_searched == 1;
    let why = match reason {
        NotSearched::LowPrevalence => "the cluster is below --min-prevalence",
        NotSearched::NoClusterSequence => "the cluster has no sequence to search with",
        NotSearched::NoFasta if one => "that genome has no FASTA",
        NotSearched::NoFasta => "those genomes have no FASTA",
    };
    let reference_loss = mechanism.mechanism == "reference_or_lineage_loss_candidate";
    let none_searched = coverage.not_searched == coverage.absent;
    if none_searched || (reference_loss && coverage.reference_untested) {
        mechanism.support = downgraded_support(mechanism.support);
        mechanism.detail.push_str(";absences_not_searched");
        return Some(if !none_searched {
            // Only reachable without a cluster-level reason: the reference
            // genome is one of the genomes without a FASTA.
            "The reference genome has no FASTA, so its absence was never searched and is the input matrix's call, not Recore's.".to_string()
        } else if coverage.absent == 1 {
            format!("The only absence behind this call was never searched ({why}); it is the input matrix's call, not Recore's.")
        } else {
            format!(
                "None of the {} absences behind this call were searched ({why}); they are the input matrix's calls, not Recore's.",
                coverage.absent
            )
        });
    }
    mechanism.detail.push_str(&format!(
        ";absences_not_searched={}/{}",
        coverage.not_searched, coverage.absent
    ));
    Some(if one {
        format!(
            "1 of the {} absences was never searched ({why}); that one is the input matrix's call.",
            coverage.absent
        )
    } else {
        format!(
            "{} of the {} absences were never searched ({why}); those are the input matrix's calls.",
            coverage.not_searched, coverage.absent
        )
    })
}

fn phylogenetic_blink_flag(phylo: &PhyloPatternSummary) -> &'static str {
    if !phylo.assessed || phylo.pattern == "absent_from_tree" {
        "not_assessed"
    } else if matches!(
        phylo.pattern,
        "homoplastic_gain_loss"
            | "dispersed_gain_loss"
            | "singleton_or_noise"
            | "ambiguous_gain_loss"
    ) || (phylo.parsimony_changes >= 6 && phylo.max_carrier_clade_fraction < 0.75)
        || (phylo.carrier_mrca_purity < 0.35 && phylo.parsimony_changes >= 4)
    {
        "high_blink"
    } else if phylo.parsimony_changes >= 4 {
        "moderate_blink"
    } else {
        "low_blink"
    }
}

pub(super) fn penalized_phylo_summary(phylo: &PhyloPatternSummary) -> PenalizedPhyloSummary {
    if !phylo.assessed {
        return PenalizedPhyloSummary {
            model: "not_assessed",
            score: 0.5,
            ..PenalizedPhyloSummary::default()
        };
    }
    if phylo.present_tips == 0 {
        return PenalizedPhyloSummary {
            model: "absent_from_tree",
            score: 0.5,
            ..PenalizedPhyloSummary::default()
        };
    }
    if phylo.absent_tips == 0 {
        return PenalizedPhyloSummary {
            model: "ancestral_retained",
            dollo_origin_tips: phylo.matched_tips,
            dollo_origin_purity: 1.0,
            score: 1.0,
            ..PenalizedPhyloSummary::default()
        };
    }

    let recurrent_gain_violation = phylo.gains.saturating_sub(1);
    let ambiguous_span_violation =
        usize::from(phylo.root_state_ambiguous && phylo.carrier_mrca_purity < 0.80);
    let dollo_violation_count = recurrent_gain_violation.max(ambiguous_span_violation);
    let dollo_loss_count = phylo.dollo_loss_events_in_mrca;
    let blink_penalty = match phylogenetic_blink_flag(phylo) {
        "high_blink" => 0.25,
        "moderate_blink" => 0.10,
        _ => 0.0,
    };
    let span_penalty = (1.0 - phylo.carrier_mrca_purity).clamp(0.0, 1.0) * 0.35;
    let score = (1.0
        - dollo_violation_count as f64 * 0.14
        - dollo_loss_count as f64 * 0.05
        - span_penalty
        - blink_penalty)
        .clamp(0.0, 1.0);
    let model = if phylo.carrier_mrca_purity >= 0.95 && dollo_loss_count == 0 {
        "single_gain_clade"
    } else if phylo.carrier_mrca_purity >= 0.80 && dollo_loss_count <= 3 {
        "single_gain_with_losses"
    } else if phylo.carrier_mrca_purity >= 0.50 && dollo_loss_count <= 6 {
        "broad_origin_many_losses"
    } else if dollo_violation_count <= 2 && score >= 0.45 {
        "weak_dollo_fit"
    } else {
        "poor_single_gain_fit"
    };

    PenalizedPhyloSummary {
        model,
        dollo_loss_count,
        dollo_violation_count,
        dollo_origin_tips: phylo.carrier_mrca_tips,
        dollo_origin_purity: phylo.carrier_mrca_purity,
        score,
    }
}

fn evolutionary_signal_score(
    phylo: &PhyloPatternSummary,
    penalized_phylo: &PenalizedPhyloSummary,
    mechanism: &MechanismCall,
    context: Option<&ClusterContextSummary>,
    component: Option<&ClusterComponentInfo>,
    rescue_evidence: &[&RescueEvent],
) -> f64 {
    let mut score: f64 = 0.35;
    score += match phylogenetic_reliability(phylo) {
        "high" => 0.25,
        "moderate" => 0.10,
        "low" => -0.20,
        _ => 0.0,
    };
    score += match mechanism.support {
        "high" => 0.25,
        "moderate" => 0.12,
        "low" => -0.10,
        "review" => -0.05,
        _ => 0.0,
    };
    score += match context.map(|c| c.interpretation) {
        Some("stable_single_locus") => 0.12,
        Some("stable_single_locus_minor_noise") => 0.08,
        Some("mostly_stable_review_minor_contexts") => 0.03,
        Some("possible_multilocus_or_paralog") => -0.22,
        _ => 0.0,
    };
    if rescue_evidence.iter().any(|r| r.synteny_ok == Some(true)) {
        score += 0.08;
    }
    if rescue_evidence.iter().any(|r| r.synteny_ok == Some(false)) {
        score -= 0.18;
    }
    if let Some(component) = component {
        if component.recommendation.contains("paralog") {
            score -= 0.22;
        } else if component.recommendation.contains("collapse_aliases") {
            score -= 0.08;
        }
    }
    match phylogenetic_blink_flag(phylo) {
        "high_blink" => score -= 0.25,
        "moderate_blink" => score -= 0.10,
        _ => {}
    }
    score += (penalized_phylo.score - 0.5) * 0.20;
    score.clamp(0.0, 1.0)
}

fn highlight_priority(
    score: f64,
    phylo: &PhyloPatternSummary,
    component: Option<&ClusterComponentInfo>,
    mechanism: &MechanismCall,
) -> &'static str {
    if component.is_some_and(|c| c.recommendation.contains("paralog")) {
        return "review_paralog_component";
    }
    if phylogenetic_blink_flag(phylo) == "high_blink" {
        return "filter_or_review_blinking";
    }
    if mechanism.mechanism == "annotation_or_clustering_artifact" {
        return "deprioritize_cluster_artifact";
    }
    if score >= 0.75 {
        "highlight"
    } else if score >= 0.55 {
        "candidate"
    } else {
        "review"
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
        "gene_cluster\tn_present_corrected\tprevalence_corrected\tabsent_not_searched\tnot_searched_reason\ttree_matched_tips\treference_state\tphylogenetic_pattern\tphylogenetic_reliability\tphylogenetic_blink_flag\tpenalized_phylo_model\tpenalized_phylo_score\tdollo_loss_count\tdollo_violation_count\tdollo_origin_tips\tdollo_origin_purity\tancestral_state\troot_state_confidence\tphylogenetic_changes\tphylogenetic_gains\tphylogenetic_losses\tcarrier_monophyly\tmax_carrier_clade_fraction\tcarrier_mrca_tips\tcarrier_mrca_purity\tmolecular_mechanism\tmechanism_detail\tsupport_level\tevolutionary_signal_score\thighlight_priority\tcluster_component_id\tcluster_component_size\tcluster_component_role\tcluster_component_relations\tcluster_component_recommendation\tis_interruption_supported\tmobile_element_source\tmobile_element_kind\tmobile_element_label\tmobile_element_genome\tmobile_element_distance_bp\tmobile_element_identity\tmobile_element_coverage\treference_event_class\ttarget_reference_cluster\ttarget_reference_gene\tcontext_interpretation\tparalog_locus_risk\tsynteny_mean_score\tsynteny_pass_fraction\torf_integrity_summary\trescue_evidence_counts\trecommended_wording"
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
        let penalized_phylo = penalized_phylo_summary(&phylo);
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
        // Which of the cluster's absences the search tested; the cluster's own
        // reason wins, as it does in the states file.
        let mut absence = AbsenceCoverage {
            reason: inputs.coverage.cluster_reason(ci),
            ..AbsenceCoverage::default()
        };
        for (gi_val, &present) in present_by_genome.iter().enumerate() {
            if present {
                continue;
            }
            absence.absent += 1;
            let gi = GenomeIdx(gi_val);
            if let Some(reason) = inputs.coverage.reason(ci, gi) {
                absence.not_searched += 1;
                absence.reason.get_or_insert(reason);
                absence.reference_untested |= Some(gi) == reference_gi;
            }
        }
        if absence.not_searched == 0 {
            absence.reason = None;
        }
        let mut mechanism = classify_mechanism(MechanismInput {
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
        // Before the score, so a call that rests on untested absences is scored
        // and prioritised as the weaker call it is.
        let not_searched_note = qualify_unsearched_absences(&mut mechanism, &absence);
        let component = inputs
            .component_info
            .get(ci_val)
            .and_then(|info| info.as_ref());
        let signal_score = evolutionary_signal_score(
            &phylo,
            &penalized_phylo,
            &mechanism,
            contexts_by_cluster[ci_val],
            component,
            &all_rescue_evidence,
        );
        let priority = highlight_priority(signal_score, &phylo, component, &mechanism);
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
                absence.not_searched.to_string(),
                absence
                    .reason
                    .map(NotSearched::tag)
                    .unwrap_or("")
                    .to_string(),
                if phylo.assessed {
                    phylo.matched_tips.to_string()
                } else {
                    String::new()
                },
                reference_state.to_string(),
                phylo.pattern.to_string(),
                phylogenetic_reliability(&phylo).to_string(),
                phylogenetic_blink_flag(&phylo).to_string(),
                penalized_phylo.model.to_string(),
                if phylo.assessed {
                    format!("{:.3}", penalized_phylo.score)
                } else {
                    String::new()
                },
                if phylo.assessed {
                    penalized_phylo.dollo_loss_count.to_string()
                } else {
                    String::new()
                },
                if phylo.assessed {
                    penalized_phylo.dollo_violation_count.to_string()
                } else {
                    String::new()
                },
                if phylo.assessed {
                    penalized_phylo.dollo_origin_tips.to_string()
                } else {
                    String::new()
                },
                if phylo.assessed {
                    format!("{:.3}", penalized_phylo.dollo_origin_purity)
                } else {
                    String::new()
                },
                phylo.root_state.to_string(),
                if phylo.assessed {
                    if phylo.root_state_ambiguous {
                        "ambiguous".to_string()
                    } else {
                        "unambiguous".to_string()
                    }
                } else {
                    String::new()
                },
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
                if phylo.assessed {
                    phylo.carrier_mrca_tips.to_string()
                } else {
                    String::new()
                },
                if phylo.assessed {
                    format!("{:.3}", phylo.carrier_mrca_purity)
                } else {
                    String::new()
                },
                mechanism.mechanism.to_string(),
                mechanism.detail,
                mechanism.support.to_string(),
                format!("{:.3}", signal_score),
                priority.to_string(),
                component
                    .map(|info| info.component_id.clone())
                    .unwrap_or_default(),
                component
                    .map(|info| info.component_size.to_string())
                    .unwrap_or_default(),
                component.map(|info| info.role.clone()).unwrap_or_default(),
                component
                    .map(|info| info.relation_types.clone())
                    .unwrap_or_default(),
                component
                    .map(|info| info.recommendation.clone())
                    .unwrap_or_default(),
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
                match not_searched_note {
                    Some(note) => format!("{} {}", mechanism.wording, note),
                    None => mechanism.wording.to_string(),
                },
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

#[derive(Debug, Clone, Default)]
struct SyntenyBlockAccumulator {
    count: usize,
    examples: Vec<String>,
    carrier_genomes: Vec<String>,
    upstream_signature: String,
    downstream_signature: String,
    representative: Option<SyntenyBlockRepresentative>,
}

#[derive(Debug, Clone)]
struct SyntenyGeneRecord {
    cluster: String,
    gene_id: String,
    contig: String,
    start: u64,
    end: u64,
    strand: char,
    length: u64,
}

#[derive(Debug, Clone)]
struct SyntenyBlockRepresentative {
    genome: String,
    gene_id: String,
    contig: String,
    start: u64,
    end: u64,
    strand: char,
    length: u64,
    genes: Vec<SyntenyGeneRecord>,
}

#[derive(Debug, Clone)]
struct ReferenceBoundaryContext {
    genes: Vec<SyntenyGeneRecord>,
    cluster_ids: Vec<ClusterIdx>,
}

#[derive(Debug, Clone, Default)]
struct ReferenceBoundaryIndex {
    contexts: Vec<ReferenceBoundaryContext>,
    by_cluster: HashMap<ClusterIdx, Vec<usize>>,
    by_focal: HashMap<ClusterIdx, Vec<SyntenyGeneRecord>>,
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
    genes: &[GffGene],
    gene_to_cluster: &GeneClusters,
    pa: &crate::parser::PresenceAbsence,
) -> String {
    let names: Vec<&str> = genes
        .iter()
        .filter_map(|gene| gene_to_cluster.cluster_for_gene(gene))
        .map(|ci| pa.cluster_name(ci))
        .collect();
    if names.is_empty() {
        "NA".to_string()
    } else {
        names.join("|")
    }
}

/// Whether two (upstream, downstream) contexts name a neighbouring cluster in
/// common. A context with no named neighbour carries no evidence either way.
pub(super) fn shares_neighbour(a: &(String, String), b: &(String, String)) -> bool {
    let names = |c: &(String, String)| -> Vec<String> {
        [&c.0, &c.1]
            .iter()
            .flat_map(|side| side.split('|'))
            .filter(|n| !n.is_empty() && *n != "NA")
            .map(str::to_string)
            .collect()
    };
    let (na, nb) = (names(a), names(b));
    na.is_empty() || nb.is_empty() || na.iter().any(|n| nb.contains(n))
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

fn synteny_gene_record(
    gff: &GffIndex,
    gene: &GffGene,
    gene_to_cluster: &GeneClusters,
    pa: &crate::parser::PresenceAbsence,
) -> SyntenyGeneRecord {
    let cluster = gene_to_cluster
        .cluster_for_gene(gene)
        .map(|ci| pa.cluster_name(ci).to_string())
        .unwrap_or_else(|| "NA".to_string());
    let display_gene_id = gff
        .aliases_for_gene(gene)
        .into_iter()
        .map(display_reference_gene_id)
        .find(|id| id.starts_with("Rv"))
        .unwrap_or_else(|| display_reference_gene_id(gene.gene_id));
    SyntenyGeneRecord {
        cluster,
        gene_id: display_gene_id,
        contig: gene.contig.to_string(),
        start: gene.start + 1,
        end: gene.end,
        strand: gene.strand,
        length: gene.end.saturating_sub(gene.start),
    }
}

fn synteny_gap_record(contig: &str, start_0based: u64, end_0based: u64) -> SyntenyGeneRecord {
    let low = start_0based.min(end_0based);
    let high = start_0based.max(end_0based);
    let start = low.saturating_add(1);
    let (end, length) = if low == high {
        (start, 0)
    } else {
        (high, high.saturating_sub(low))
    };
    SyntenyGeneRecord {
        cluster: "__GAP__".to_string(),
        gene_id: "H37Rv_context_slot".to_string(),
        contig: contig.to_string(),
        start,
        end,
        strand: '.',
        length,
    }
}

fn ordered_synteny_gene_records(
    gff: &GffIndex,
    upstream: &[GffGene],
    center: Option<&GffGene>,
    downstream: &[GffGene],
    gene_to_cluster: &GeneClusters,
    pa: &crate::parser::PresenceAbsence,
) -> Vec<SyntenyGeneRecord> {
    let mut records =
        Vec::with_capacity(upstream.len() + downstream.len() + usize::from(center.is_some()));
    for gene in upstream.iter().rev() {
        records.push(synteny_gene_record(gff, gene, gene_to_cluster, pa));
    }
    if let Some(gene) = center {
        records.push(synteny_gene_record(gff, gene, gene_to_cluster, pa));
    }
    for gene in downstream {
        records.push(synteny_gene_record(gff, gene, gene_to_cluster, pa));
    }
    records
}

fn encode_synteny_gene_records(records: &[SyntenyGeneRecord]) -> String {
    records
        .iter()
        .map(|record| {
            format!(
                "{},{},{},{},{},{},{}",
                record.cluster,
                record.gene_id,
                record.strand,
                record.contig,
                record.start,
                record.end,
                record.length
            )
        })
        .collect::<Vec<_>>()
        .join(";")
}

fn synteny_gene_identity(record: &SyntenyGeneRecord) -> &str {
    let gene_id = record.gene_id.trim();
    if !gene_id.is_empty() && gene_id != "H37Rv_context_slot" {
        gene_id
    } else {
        record.cluster.trim()
    }
}

fn synteny_gene_identity_signature(records: &[SyntenyGeneRecord]) -> String {
    let names: Vec<&str> = records
        .iter()
        .filter_map(|record| match synteny_gene_identity(record) {
            "" | "NA" => None,
            identity => Some(identity),
        })
        .collect();
    if names.is_empty() {
        "NA".to_string()
    } else {
        names.join("|")
    }
}

fn build_reference_boundary_index(
    gff: &GffIndex,
    gene_to_cluster: &GeneClusters,
    pa: &crate::parser::PresenceAbsence,
    n_neighbors: usize,
) -> ReferenceBoundaryIndex {
    let mut index = ReferenceBoundaryIndex::default();
    let all: Vec<GffGene> = gff.genes().collect();
    if all.is_empty() {
        return index;
    }

    let mut range_start = 0usize;
    while range_start < all.len() {
        let contig = all[range_start].contig.to_string();
        let mut range_end = range_start + 1;
        while range_end < all.len() && all[range_end].contig == contig {
            range_end += 1;
        }
        let slice = &all[range_start..range_end];
        for (gene_idx, gene) in slice.iter().enumerate() {
            let Some(focal_ci) = gene_to_cluster.cluster_for_gene(gene) else {
                continue;
            };
            index.by_focal.entry(focal_ci).or_insert_with(|| {
                let upstream: Vec<GffGene> = slice[..gene_idx]
                    .iter()
                    .rev()
                    .take(n_neighbors)
                    .copied()
                    .collect();
                let downstream: Vec<GffGene> = slice[gene_idx + 1..]
                    .iter()
                    .take(n_neighbors)
                    .copied()
                    .collect();
                ordered_synteny_gene_records(
                    gff,
                    &upstream,
                    Some(gene),
                    &downstream,
                    gene_to_cluster,
                    pa,
                )
            });
        }
        for boundary in 0..=slice.len() {
            let upstream: Vec<GffGene> = slice[..boundary]
                .iter()
                .rev()
                .take(n_neighbors)
                .copied()
                .collect();
            let downstream: Vec<GffGene> = slice[boundary..]
                .iter()
                .take(n_neighbors)
                .copied()
                .collect();
            let gap_start = if boundary > 0 {
                slice[boundary - 1].end
            } else {
                slice.first().map(|gene| gene.start).unwrap_or(0)
            };
            let gap_end = if boundary < slice.len() {
                slice[boundary].start
            } else {
                slice.last().map(|gene| gene.end).unwrap_or(gap_start)
            };
            let mut genes = ordered_synteny_gene_records(
                gff,
                &upstream,
                None,
                &downstream,
                gene_to_cluster,
                pa,
            );
            let gap_idx = upstream.len();
            genes.insert(gap_idx, synteny_gap_record(&contig, gap_start, gap_end));
            let mut cluster_ids: Vec<ClusterIdx> = upstream
                .iter()
                .chain(downstream.iter())
                .filter_map(|gene| gene_to_cluster.cluster_for_gene(gene))
                .collect();
            cluster_ids.sort_unstable();
            cluster_ids.dedup();
            let context_idx = index.contexts.len();
            for ci in &cluster_ids {
                index.by_cluster.entry(*ci).or_default().push(context_idx);
            }
            index
                .contexts
                .push(ReferenceBoundaryContext { genes, cluster_ids });
        }
        range_start = range_end;
    }
    index
}

fn parse_context_cluster_indices(
    upstream: &str,
    downstream: &str,
    name_to_idx: &HashMap<String, ClusterIdx>,
) -> Vec<ClusterIdx> {
    let mut out: Vec<ClusterIdx> = upstream
        .split('|')
        .chain(downstream.split('|'))
        .map(str::trim)
        .filter(|value| !value.is_empty() && *value != "NA")
        .filter_map(|name| name_to_idx.get(name).copied())
        .collect();
    out.sort_unstable();
    out.dedup();
    out
}

fn best_reference_context(
    index: &ReferenceBoundaryIndex,
    focal_cluster: ClusterIdx,
    target_clusters: &[ClusterIdx],
) -> Vec<SyntenyGeneRecord> {
    if let Some(centered_context) = index.by_focal.get(&focal_cluster) {
        return centered_context.clone();
    }
    if target_clusters.is_empty() || index.contexts.is_empty() {
        return Vec::new();
    }
    let mut candidate_counts: HashMap<usize, usize> = HashMap::default();
    for ci in target_clusters {
        if let Some(boundaries) = index.by_cluster.get(ci) {
            for boundary_idx in boundaries {
                *candidate_counts.entry(*boundary_idx).or_default() += 1;
            }
        }
    }
    let Some((best_idx, _)) = candidate_counts.into_iter().max_by(|a, b| {
        a.1.cmp(&b.1).then_with(|| {
            index.contexts[a.0]
                .cluster_ids
                .len()
                .cmp(&index.contexts[b.0].cluster_ids.len())
                .reverse()
        })
    }) else {
        return Vec::new();
    };
    index.contexts[best_idx].genes.clone()
}

/// Runs of consecutive clusters to build in parallel and write in order. A run
/// is cut by the genomes its clusters are present in, not by a count of
/// clusters, so what one holds does not grow with the number of genomes: a
/// cluster for every thread at least, then clusters until the run's
/// genome-cluster pairs reach the budget.
fn cluster_chunks(
    pa: &crate::parser::PresenceAbsence,
    threads: usize,
) -> Vec<std::ops::Range<usize>> {
    const PAIRS_PER_CHUNK: usize = 32 * 1024;
    let mut chunks = Vec::new();
    let mut chunk_start = 0usize;
    while chunk_start < pa.n_clusters() {
        let mut chunk_end = chunk_start;
        let mut pairs = 0usize;
        while chunk_end < pa.n_clusters()
            && (chunk_end - chunk_start < threads || pairs < PAIRS_PER_CHUNK)
        {
            pairs += pa.present_count(ClusterIdx(chunk_end));
            chunk_end += 1;
        }
        chunks.push(chunk_start..chunk_end);
        chunk_start = chunk_end;
    }
    chunks
}

/// One cluster's locus contexts: its summary, and its context rows ranked by
/// carriers.
fn cluster_locus_contexts(
    ci: ClusterIdx,
    pa: &crate::parser::PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    gene_to_cluster_maps: &HashMap<GenomeIdx, GeneClusters>,
) -> (
    ClusterContextSummary,
    Vec<(String, String, ContextAccumulator)>,
) {
    let mut context_rows: BTreeMap<(String, String), ContextAccumulator> = BTreeMap::new();
    let mut n_mapped = 0usize;
    let mut n_unresolved = 0usize;

    let mut genomes = pa.present_genomes(ci).to_vec();
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

        let (upstream, downstream) =
            gff.context_window(gene.contig, gene.start, gene.end, Some(gene.gene_id), 2);
        let upstream_sig = context_signature(&upstream, g2c, pa);
        let downstream_sig = context_signature(&downstream, g2c, pa);
        n_mapped += 1;

        let row = context_rows
            .entry((upstream_sig, downstream_sig))
            .or_default();
        row.count += 1;
        if row.examples.len() < 5 {
            row.examples
                .push(format!("{}:{}", pa.genome_name(gi), gene.gene_id));
        }
    }

    let n_unique_contexts = context_rows.len();
    let dominant_count = context_rows
        .values()
        .map(|row| row.count)
        .max()
        .unwrap_or(0);
    let dominant_fraction = if n_mapped > 0 {
        dominant_count as f64 / n_mapped as f64
    } else {
        0.0
    };
    // Judged by neighbourhood, not by exact signature. A context that shares
    // even one neighbouring cluster with the dominant one is the same locus
    // with something gained, lost or unannotated next to it: a pseudogene,
    // a missing CDS or a true loss beside a single-copy gene changes its
    // signature without making it a second copy. Only a context with no
    // neighbour in common is evidence of another locus. Signature
    // differences alone used to label single-copy genes as paralogs.
    let dominant = context_rows
        .iter()
        .max_by(|a, b| a.1.count.cmp(&b.1.count).then_with(|| b.0.cmp(a.0)))
        .map(|(key, _)| key.clone());
    let compatible: usize = context_rows
        .iter()
        .filter(|(key, _)| dominant.as_ref().is_none_or(|d| shares_neighbour(d, key)))
        .map(|(_, row)| row.count)
        .sum();
    let compatible_fraction = if n_mapped > 0 {
        compatible as f64 / n_mapped as f64
    } else {
        0.0
    };
    let interpretation = context_interpretation(n_mapped, n_unique_contexts, compatible_fraction);
    let summary = ClusterContextSummary {
        cluster_idx: ci,
        n_present: pa.present_count(ci),
        n_mapped,
        n_unresolved,
        n_unique_contexts,
        dominant_count,
        dominant_fraction,
        interpretation,
    };

    let mut rows: Vec<(String, String, ContextAccumulator)> = context_rows
        .into_iter()
        .map(|((upstream, downstream), acc)| (upstream, downstream, acc))
        .collect();
    rows.sort_by(|a, b| {
        b.2.count
            .cmp(&a.2.count)
            .then(a.0.cmp(&b.0))
            .then(a.1.cmp(&b.1))
    });
    (summary, rows)
}

/// Each cluster's neighbourhoods, in `cluster_locus_context.tsv`, and a summary
/// per cluster, in `cluster_locus_summary.tsv`.
///
/// One cluster at a time in one thread, this was a sixth of a run on 400
/// genomes and grew with every genome. Clusters do not depend on each other,
/// so they are built in parallel a run at a time and written in cluster order,
/// the same rows in the same order.
pub(super) fn write_cluster_locus_contexts(
    pa: &crate::parser::PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    gene_to_cluster_maps: &HashMap<GenomeIdx, GeneClusters>,
    threads: usize,
    output_dir: &std::path::Path,
) -> Result<(usize, Vec<ClusterContextSummary>)> {
    let threads = threads.max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("Failed to build the locus context thread pool")?;

    let context_path = output_dir.join("cluster_locus_context.tsv");
    let context_file = fs::File::create(&context_path)?;
    let mut context_w = BufWriter::new(context_file);
    writeln!(
        context_w,
        "gene_cluster\tcontext_rank\tn_genomes\tfraction_of_mapped_contexts\tupstream_context\tdownstream_context\texample_genomes"
    )?;
    let mut summaries = Vec::with_capacity(pa.n_clusters());
    let mut written_context_rows = 0usize;
    for chunk in cluster_chunks(pa, threads) {
        let built: Vec<_> = pool.install(|| {
            use rayon::prelude::*;
            chunk
                .into_par_iter()
                .map(|ci_val| {
                    cluster_locus_contexts(
                        ClusterIdx(ci_val),
                        pa,
                        gff_indices,
                        gene_to_cluster_maps,
                    )
                })
                .collect()
        });
        for (summary, rows) in built {
            let mapped = summary.n_mapped.max(1);
            for (rank, (upstream, downstream, acc)) in rows.into_iter().enumerate() {
                write_tsv_fields(
                    &mut context_w,
                    &[
                        pa.cluster_name(summary.cluster_idx).to_string(),
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
            summaries.push(summary);
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
                // 0 of 0 is no reading, not a fraction of 0: blank, so a filter
                // for multi-locus clusters (fraction < 0.5) does not pick it up.
                if summary.n_mapped == 0 {
                    String::new()
                } else {
                    format!("{:.3}", summary.dominant_fraction)
                },
                summary.interpretation.to_string(),
            ],
        )?;
    }
    summary_w.flush()?;

    Ok((written_context_rows, summaries))
}

/// Unique flanking-gene blocks per cluster, in `cluster_synteny_blocks.tsv`.
///
/// Every cluster's blocks used to be gathered before the first row was written,
/// in one thread: on 200 genomes of 4,000 genes that held 1.2 GB and took half
/// the run. Clusters do not depend on each other, so they are now built in
/// parallel a chunk at a time and written in cluster order, the same rows in
/// the same order, holding one chunk's rows at a time.
///
/// A block is a genome's neighbourhood, so where annotations share no gene IDs
/// every genome is a block of its own: chunks of 256 clusters held 200 MB of
/// rows on 400 genomes, more with every genome added. Chunks are now cut by
/// genome-cluster pairs; see cluster_chunks.
pub(super) fn write_cluster_synteny_blocks(
    pa: &crate::parser::PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    gene_to_cluster_maps: &HashMap<GenomeIdx, GeneClusters>,
    reference_gi: Option<GenomeIdx>,
    threads: usize,
    output_dir: &std::path::Path,
) -> Result<usize> {
    let reference_boundary_index = reference_gi
        .and_then(|gi| gff_indices.get(&gi).zip(gene_to_cluster_maps.get(&gi)))
        .map(|(gff, g2c)| build_reference_boundary_index(gff, g2c, pa, 5))
        .unwrap_or_default();
    let cluster_index = pa.cluster_index();

    let path = output_dir.join("cluster_synteny_blocks.tsv");
    let file = fs::File::create(&path)?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "gene_cluster\tblock_rank\tn_genomes\tfraction_of_present\tupstream_block_5\tdownstream_block_5\tblock_interpretation\texample_genomes\tcarrier_genomes\trepresentative_genome\trepresentative_gene\trepresentative_contig\trepresentative_start\trepresentative_end\trepresentative_strand\trepresentative_length_bp\trepresentative_block_genes\treference_block_genes"
    )?;

    let threads = threads.max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("Failed to build the synteny block thread pool")?;
    let mut written = 0usize;
    for chunk in cluster_chunks(pa, threads) {
        let rows: Vec<Result<(Vec<u8>, usize)>> = pool.install(|| {
            use rayon::prelude::*;
            chunk
                .into_par_iter()
                .map(|ci_val| {
                    cluster_synteny_block_rows(
                        ClusterIdx(ci_val),
                        pa,
                        gff_indices,
                        gene_to_cluster_maps,
                        &reference_boundary_index,
                        &cluster_index,
                    )
                })
                .collect()
        });
        for cluster_rows in rows {
            let (lines, n_rows) = cluster_rows?;
            w.write_all(&lines)?;
            written += n_rows;
        }
    }
    w.flush()?;
    Ok(written)
}

/// One cluster's synteny block rows, ranked by carriers, as TSV lines, and how
/// many there are.
fn cluster_synteny_block_rows(
    ci: ClusterIdx,
    pa: &crate::parser::PresenceAbsence,
    gff_indices: &HashMap<GenomeIdx, GffIndex>,
    gene_to_cluster_maps: &HashMap<GenomeIdx, GeneClusters>,
    reference_boundary_index: &ReferenceBoundaryIndex,
    cluster_index: &HashMap<String, ClusterIdx>,
) -> Result<(Vec<u8>, usize)> {
    let mut blocks: BTreeMap<String, SyntenyBlockAccumulator> = BTreeMap::new();
    let mut genomes = pa.present_genomes(ci).to_vec();
    genomes.sort_by_key(|gi| gi.val());
    for gi in genomes {
        let Some(gff) = gff_indices.get(&gi) else {
            continue;
        };
        let Some(g2c) = gene_to_cluster_maps.get(&gi) else {
            continue;
        };
        let Some(gene) = gff_gene_for_cell(gff, pa.cell(ci, gi)) else {
            continue;
        };

        let (upstream, downstream) =
            gff.context_window(gene.contig, gene.start, gene.end, Some(gene.gene_id), 5);
        let upstream_sig = context_signature(&upstream, g2c, pa);
        let downstream_sig = context_signature(&downstream, g2c, pa);
        let representative_genes =
            ordered_synteny_gene_records(gff, &upstream, Some(&gene), &downstream, g2c, pa);
        let block_signature = synteny_gene_identity_signature(&representative_genes);
        let row = blocks.entry(block_signature).or_default();
        row.count += 1;
        if row.upstream_signature.is_empty() {
            row.upstream_signature = upstream_sig;
            row.downstream_signature = downstream_sig;
        }
        row.carrier_genomes.push(pa.genome_name(gi).to_string());
        if row.examples.len() < 5 {
            row.examples
                .push(format!("{}:{}", pa.genome_name(gi), gene.gene_id));
        }
        if row.representative.is_none() {
            row.representative = Some(SyntenyBlockRepresentative {
                genome: pa.genome_name(gi).to_string(),
                gene_id: gene.gene_id.to_string(),
                contig: gene.contig.to_string(),
                start: gene.start + 1,
                end: gene.end,
                strand: gene.strand,
                length: gene.end.saturating_sub(gene.start),
                genes: representative_genes,
            });
        }
    }

    let mut rows: Vec<(String, SyntenyBlockAccumulator)> = blocks.into_iter().collect();
    rows.sort_by(|a, b| b.1.count.cmp(&a.1.count).then(a.0.cmp(&b.0)));
    let present = pa.present_count(ci).max(1);
    let n_blocks = rows.len();
    // Each row is formatted and written as it is made.
    let row_fields = rows
        .into_iter()
        .enumerate()
        .map(|(rank, (_signature, acc))| {
            let fraction = acc.count as f64 / present as f64;
            let interpretation = if n_blocks <= 1 {
                "single_synteny_block"
            } else if rank == 0 && fraction >= 0.80 {
                "dominant_block_minor_alternatives"
            } else {
                "alternative_synteny_block"
            };
            let representative = acc.representative.as_ref();
            let target_clusters = parse_context_cluster_indices(
                &acc.upstream_signature,
                &acc.downstream_signature,
                cluster_index,
            );
            let reference_genes =
                best_reference_context(reference_boundary_index, ci, &target_clusters);
            vec![
                pa.cluster_name(ci).to_string(),
                (rank + 1).to_string(),
                acc.count.to_string(),
                format!("{:.3}", fraction),
                acc.upstream_signature.clone(),
                acc.downstream_signature.clone(),
                interpretation.to_string(),
                acc.examples.join(";"),
                acc.carrier_genomes.join(";"),
                representative
                    .map(|rep| rep.genome.clone())
                    .unwrap_or_default(),
                representative
                    .map(|rep| rep.gene_id.clone())
                    .unwrap_or_default(),
                representative
                    .map(|rep| rep.contig.clone())
                    .unwrap_or_default(),
                representative
                    .map(|rep| rep.start.to_string())
                    .unwrap_or_default(),
                representative
                    .map(|rep| rep.end.to_string())
                    .unwrap_or_default(),
                representative
                    .map(|rep| rep.strand.to_string())
                    .unwrap_or_default(),
                representative
                    .map(|rep| rep.length.to_string())
                    .unwrap_or_default(),
                representative
                    .map(|rep| encode_synteny_gene_records(&rep.genes))
                    .unwrap_or_default(),
                encode_synteny_gene_records(&reference_genes),
            ]
        });
    let mut lines = Vec::new();
    for fields in row_fields {
        write_tsv_fields(&mut lines, &fields)?;
    }
    Ok((lines, n_blocks))
}
