//! Scale harness for the HTML report. Not a correctness test: it synthesises a
//! realistic run of N genomes x M clusters in memory, the way the rich-sample
//! test builds its inputs, and writes the report, so the browser side can be
//! measured at the sizes the report promises to handle.
//!
//!   RECORE_SCALE_GENOMES=10000 RECORE_SCALE_CLUSTERS=30000 \
//!   RECORE_SCALE_OUT=/abs/report_10000x30000.html \
//!   cargo test --release --lib scale_report -- --ignored --nocapture
//!
//! Optional knobs:
//!   RECORE_SCALE_SEED      PRNG seed (default 42)
//!   RECORE_SCALE_SYNTENY   "some" (default: a quarter of clusters, 1-4 blocks),
//!                          "full" (every cluster with carriers, the way the real
//!                          writer emits one row per distinct context), or "none"
//!   RECORE_SCALE_TREE      "mtbc" (default) or "caterpillar" (worst-case depth)
//!   RECORE_SCALE_KEEP_WORK keep the side TSVs next to the HTML
//!
//! It prints one `SCALE {...}` JSON line with phase timings, peak RSS before and
//! after the report call, and the composition of the matrix it generated.

use super::*;
use crate::parser::PresenceAbsence;
use std::time::Instant;

struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        // splitmix64
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn below(&mut self, n: usize) -> usize {
        if n == 0 {
            0
        } else {
            (self.next_u64() % n as u64) as usize
        }
    }
    fn range(&mut self, lo: usize, hi: usize) -> usize {
        lo + self.below(hi.saturating_sub(lo) + 1)
    }
    fn chance(&mut self, p: f64) -> bool {
        self.f64() < p
    }
    fn exp(&mut self, mean: f64) -> f64 {
        -mean * (1.0 - self.f64()).ln()
    }
    fn shuffle<T>(&mut self, items: &mut [T]) {
        for i in (1..items.len()).rev() {
            let j = self.below(i + 1);
            items.swap(i, j);
        }
    }
}

/// Resident set size of this process right now, in bytes. The crate forbids
/// unsafe code, so this asks `ps` rather than calling getrusage.
fn current_rss_bytes() -> u64 {
    let pid = std::process::id().to_string();
    std::process::Command::new("ps")
        .args(["-o", "rss=", "-p", &pid])
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .and_then(|s| s.trim().parse::<u64>().ok())
        .map(|kib| kib * 1024)
        .unwrap_or(0)
}

/// Polls RSS every 50 ms and keeps the maximum since the last `take`.
struct RssSampler {
    max: std::sync::Arc<std::sync::atomic::AtomicU64>,
    stop: std::sync::Arc<std::sync::atomic::AtomicBool>,
}

impl RssSampler {
    fn start() -> Self {
        use std::sync::atomic::Ordering;
        let max = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
        let stop = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
        let (m, s) = (max.clone(), stop.clone());
        std::thread::spawn(move || {
            while !s.load(Ordering::Relaxed) {
                m.fetch_max(current_rss_bytes(), Ordering::Relaxed);
                std::thread::sleep(std::time::Duration::from_millis(50));
            }
        });
        RssSampler { max, stop }
    }
    /// Maximum seen since the last call, then restart from the current value.
    fn take(&self) -> u64 {
        use std::sync::atomic::Ordering;
        let now = current_rss_bytes();
        self.max.fetch_max(now, Ordering::Relaxed);
        self.max.swap(now, Ordering::Relaxed)
    }
}

impl Drop for RssSampler {
    fn drop(&mut self) {
        self.stop.store(true, std::sync::atomic::Ordering::Relaxed);
    }
}

/// Global MTBC collections are dominated by L4 and L2; every lineage gets at
/// least one genome so the lineage colours and the L8 rooting are exercised.
const LINEAGES: &[(&str, f64)] = &[
    ("L4", 0.44),
    ("L2", 0.22),
    ("L1", 0.09),
    ("L3", 0.08),
    ("L5", 0.03),
    ("L6", 0.03),
    ("L7", 0.012),
    ("L8", 0.004),
    ("L9", 0.004),
    ("L10", 0.002),
    ("A1", 0.02),
    ("A2", 0.015),
    ("A3", 0.01),
    ("A4", 0.033),
];

const ANNOTATIONS: &[&str] = &[
    "hypothetical protein",
    "conserved hypothetical protein",
    "PE family protein",
    "PPE family protein PPE38",
    "PE-PGRS family protein PE_PGRS33",
    "IS6110 transposase",
    "IS1081 transposase",
    "ABC transporter ATP-binding protein",
    "conserved membrane protein",
    "toxin VapC",
    "antitoxin VapB",
    "lipoprotein LpqB",
    "ESX-1 secretion-associated protein EspA",
    "transcriptional regulator, TetR family",
    "two-component sensor histidine kinase",
    "polyketide synthase Pks12",
    "mycocerosic acid synthase",
    "MCE-family protein Mce1A",
    "cytochrome P450 Cyp125",
    "probable oxidoreductase",
];

#[derive(Clone, Copy, PartialEq)]
enum Category {
    Core,
    Deletion,
    Shell,
    Cloud,
    NscMember,
}

struct ClusterDraft {
    name: String,
    annotation: String,
    category: Category,
    present: Vec<u32>,
    rescued: Vec<u32>,
    pseudogene: Vec<u32>,
    /// Absent cells inside the carrying clade that nothing explains, the pool an
    /// NSC member is drawn from.
    unexplained: Vec<u32>,
}

fn leaf(name: &str, gi: usize, rng: &mut Rng) -> PhyloNode {
    PhyloNode {
        label: Some(name.to_string()),
        genome_idx: Some(GenomeIdx(gi)),
        branch_length: Some(0.00002 + rng.exp(0.00015)),
        children: Vec::new(),
    }
}

fn clade(children: Vec<PhyloNode>, length: f64) -> PhyloNode {
    PhyloNode {
        label: None,
        genome_idx: None,
        branch_length: Some(length),
        children,
    }
}

/// Random-split binary tree over the given tips (a random-BST shape, depth
/// roughly 4 ln n, about 40 at 10,000 tips).
fn random_subtree(tips: &[usize], genomes: &[String], rng: &mut Rng) -> PhyloNode {
    if tips.len() == 1 {
        return leaf(&genomes[tips[0]], tips[0], rng);
    }
    let split = 1 + rng.below(tips.len() - 1);
    let left = random_subtree(&tips[..split], genomes, rng);
    let right = random_subtree(&tips[split..], genomes, rng);
    clade(vec![left, right], 0.00003 + rng.exp(0.0002))
}

fn caterpillar(tips: &[usize], genomes: &[String], rng: &mut Rng) -> PhyloNode {
    let mut node = leaf(&genomes[tips[tips.len() - 1]], tips[tips.len() - 1], rng);
    for &gi in tips[..tips.len() - 1].iter().rev() {
        node = clade(vec![leaf(&genomes[gi], gi, rng), node], 0.00005);
    }
    node
}

/// Tip ranges of every internal node, in the DFS leaf order also returned.
fn collect_clades(node: &PhyloNode, order: &mut Vec<u32>, clades: &mut Vec<(u32, u32)>) {
    // Iterative post-order so a caterpillar does not overflow the stack.
    let mut stack: Vec<(&PhyloNode, usize, u32)> = vec![(node, 0, 0)];
    while let Some((n, child, start)) = stack.pop() {
        if n.children.is_empty() {
            if let Some(gi) = n.genome_idx {
                order.push(gi.val() as u32);
            }
            continue;
        }
        let start = if child == 0 {
            order.len() as u32
        } else {
            start
        };
        if child < n.children.len() {
            stack.push((n, child + 1, start));
            stack.push((&n.children[child], 0, 0));
        } else {
            clades.push((start, order.len() as u32 - start));
        }
    }
}

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

#[test]
#[ignore = "scale harness: run explicitly with RECORE_SCALE_* set"]
fn scale_report() {
    let n_genomes = env_usize("RECORE_SCALE_GENOMES", 100).max(20);
    let n_clusters = env_usize("RECORE_SCALE_CLUSTERS", 500).max(20);
    let seed = env_usize("RECORE_SCALE_SEED", 42) as u64;
    let synteny_mode = std::env::var("RECORE_SCALE_SYNTENY").unwrap_or_else(|_| "some".into());
    let tree_mode = std::env::var("RECORE_SCALE_TREE").unwrap_or_else(|_| "mtbc".into());
    let out_html = std::path::PathBuf::from(
        std::env::var("RECORE_SCALE_OUT")
            .unwrap_or_else(|_| format!("report_{n_genomes}x{n_clusters}.html")),
    );
    let work = out_html.with_extension("work");
    std::fs::create_dir_all(&work).unwrap();
    let mut rng = Rng(seed);
    let sampler = RssSampler::start();
    let t_start = Instant::now();

    // ---- Genomes, named the way MTBC collections are ----
    let weight_sum: f64 = LINEAGES.iter().map(|(_, w)| w).sum();
    let mut counts: Vec<usize> = LINEAGES
        .iter()
        .map(|(_, w)| ((w / weight_sum) * n_genomes as f64).round().max(1.0) as usize)
        .collect();
    let assigned: usize = counts.iter().sum();
    if assigned > n_genomes {
        counts[0] -= assigned - n_genomes;
    } else {
        counts[0] += n_genomes - assigned;
    }
    let mut accession: Vec<usize> = (1..=n_genomes).collect();
    rng.shuffle(&mut accession);
    let mut named: Vec<(String, usize)> = Vec::with_capacity(n_genomes);
    let mut k = 0;
    for (li, (lineage, _)) in LINEAGES.iter().enumerate() {
        for _ in 0..counts[li] {
            named.push((format!("{lineage}_AS{:05}", accession[k]), li));
            k += 1;
        }
    }
    // Matrix columns in name order, the way an input directory lists them.
    named.sort();
    let genomes: Vec<String> = named.iter().map(|(n, _)| n.clone()).collect();
    let mut by_lineage: Vec<Vec<usize>> = vec![Vec::new(); LINEAGES.len()];
    for (gi, (_, li)) in named.iter().enumerate() {
        by_lineage[*li].push(gi);
    }

    // ---- Tree: MTBC backbone, random clades inside each lineage ----
    let tree_root = if tree_mode == "caterpillar" {
        let mut tips: Vec<usize> = (0..n_genomes).collect();
        rng.shuffle(&mut tips);
        caterpillar(&tips, &genomes, &mut rng)
    } else {
        let mut sub: Vec<Option<PhyloNode>> = Vec::new();
        for tips in by_lineage.iter_mut() {
            rng.shuffle(tips);
            let mut node = random_subtree(tips, &genomes, &mut rng);
            node.branch_length = Some(0.001 + rng.exp(0.0015));
            sub.push(Some(node));
        }
        let mut take = |lineage: &str| -> PhyloNode {
            let li = LINEAGES.iter().position(|(l, _)| *l == lineage).unwrap();
            sub[li].take().unwrap()
        };
        let animal = clade(
            vec![
                clade(vec![take("A1"), take("A2")], 0.0008),
                clade(vec![take("A3"), take("A4")], 0.0009),
            ],
            0.0011,
        );
        let modern = clade(
            vec![
                take("L7"),
                clade(
                    vec![take("L4"), clade(vec![take("L2"), take("L3")], 0.0006)],
                    0.0007,
                ),
            ],
            0.0009,
        );
        let africanum = clade(
            vec![
                take("L5"),
                clade(
                    vec![
                        take("L6"),
                        clade(
                            vec![take("L9"), clade(vec![take("L10"), animal], 0.0005)],
                            0.0006,
                        ),
                    ],
                    0.0007,
                ),
            ],
            0.0010,
        );
        PhyloNode {
            label: None,
            genome_idx: None,
            branch_length: None,
            children: vec![
                take("L8"),
                clade(
                    vec![take("L1"), clade(vec![modern, africanum], 0.0012)],
                    0.0015,
                ),
            ],
        }
    };
    let mut leaf_order = Vec::with_capacity(n_genomes);
    let mut clades = Vec::new();
    collect_clades(&tree_root, &mut leaf_order, &mut clades);
    clades.sort_by_key(|&(_, size)| size);
    let tree = PhyloTree { root: tree_root };
    let pick_clade = |rng: &mut Rng, lo: usize, hi: usize| -> (u32, u32) {
        let from = clades.partition_point(|&(_, s)| (s as usize) < lo);
        let to = clades.partition_point(|&(_, s)| (s as usize) <= hi);
        if from < to {
            clades[from + rng.below(to - from)]
        } else if from < clades.len() {
            clades[from]
        } else {
            *clades.last().unwrap()
        }
    };

    // ---- Clusters: a U-shaped pangenome ----
    let mut drafts: Vec<ClusterDraft> = Vec::with_capacity(n_clusters);
    let mut rv = 1usize;
    let mut group = 1usize;
    let mut nsc_pairs: Vec<(usize, usize)> = Vec::new();
    let all: Vec<u32> = (0..n_genomes as u32).collect();
    while drafts.len() < n_clusters {
        let roll = rng.f64();
        let category = if roll < 0.12 {
            Category::Core
        } else if roll < 0.20 {
            Category::Deletion
        } else if roll < 0.35 {
            Category::Shell
        } else {
            Category::Cloud
        };
        let name = if matches!(category, Category::Core | Category::Deletion) && rv < 3990 {
            let suffix = if rng.chance(0.3) { "c" } else { "" };
            let n = format!("Rv{rv:04}{suffix}");
            rv += 1;
            n
        } else {
            let n = format!("group_{group}");
            group += 1;
            n
        };
        let annotation = ANNOTATIONS[rng.below(ANNOTATIONS.len())].to_string();
        let mut draft = ClusterDraft {
            name,
            annotation,
            category,
            present: Vec::new(),
            rescued: Vec::new(),
            pseudogene: Vec::new(),
            unexplained: Vec::new(),
        };
        // Which genomes carry the locus at all, before annotation loses some.
        let mut carriers: Vec<u32> = match category {
            Category::Core => all.clone(),
            Category::Deletion => {
                let (start, size) = pick_clade(&mut rng, 2, (n_genomes / 5).max(2));
                let mut lost = vec![false; n_genomes];
                for &g in &leaf_order[start as usize..(start + size) as usize] {
                    lost[g as usize] = true;
                }
                all.iter().copied().filter(|&g| !lost[g as usize]).collect()
            }
            Category::Shell => {
                let lo = (n_genomes / 20).max(2);
                let hi = (n_genomes * 3 / 5).max(lo);
                let (start, size) = pick_clade(&mut rng, lo, hi);
                leaf_order[start as usize..(start + size) as usize].to_vec()
            }
            _ => {
                if rng.chance(0.35) {
                    // Homoplastic: scattered across the tree.
                    let k = rng.range(1, 5);
                    (0..k).map(|_| rng.below(n_genomes) as u32).collect()
                } else {
                    let (start, size) = pick_clade(&mut rng, 1, 50.min(n_genomes));
                    leaf_order[start as usize..(start + size) as usize].to_vec()
                }
            }
        };
        carriers.sort_unstable();
        carriers.dedup();
        // False absences: annotation or clustering lost the gene in some carriers.
        let loss = match category {
            Category::Core | Category::Deletion => 0.004,
            Category::Shell => 0.03,
            _ => 0.02,
        };
        for g in carriers {
            if rng.chance(loss) {
                let r = rng.f64();
                if r < 0.6 {
                    draft.rescued.push(g);
                } else if r < 0.7 {
                    draft.pseudogene.push(g);
                } else {
                    draft.unexplained.push(g);
                }
            } else {
                draft.present.push(g);
            }
        }
        if draft.present.is_empty() {
            continue;
        }
        let idx = drafts.len();
        let wants_nsc = category == Category::Shell
            && !draft.unexplained.is_empty()
            && rng.chance(0.15)
            && drafts.len() + 1 < n_clusters;
        let member_present: Vec<u32> = if wants_nsc {
            draft
                .unexplained
                .iter()
                .copied()
                .filter(|_| rng.chance(0.7))
                .collect()
        } else {
            Vec::new()
        };
        drafts.push(draft);
        if !member_present.is_empty() {
            let name = format!("group_{group}");
            group += 1;
            drafts.push(ClusterDraft {
                name,
                annotation: drafts[idx].annotation.clone(),
                category: Category::NscMember,
                present: member_present,
                rescued: Vec::new(),
                pseudogene: Vec::new(),
                unexplained: Vec::new(),
            });
            nsc_pairs.push((drafts.len() - 1, idx));
        }
    }
    drafts.truncate(n_clusters);
    nsc_pairs.retain(|&(m, r)| m < drafts.len() && r < drafts.len());

    // ---- The matrix, as the parser would hold it ----
    let short_ids: Vec<String> = genomes
        .iter()
        .map(|g| g.rsplit('_').next().unwrap_or(g).to_string())
        .collect();
    let mut cell_values: Vec<Vec<String>> = Vec::with_capacity(drafts.len());
    let mut present_cells = 0usize;
    for (ci, draft) in drafts.iter().enumerate() {
        let mut row = vec![String::new(); n_genomes];
        for &g in &draft.present {
            // Same allocation class as a real locus tag (16 bytes).
            row[g as usize] = format!("{}_{:05}", short_ids[g as usize], ci % 100_000);
            present_cells += 1;
        }
        cell_values.push(row);
    }
    let cluster_names: Vec<String> = drafts.iter().map(|d| d.name.clone()).collect();
    let annotations: Vec<String> = drafts.iter().map(|d| d.annotation.clone()).collect();
    let pa = PresenceAbsence::new_for_test_with_annotations(
        cluster_names.clone(),
        annotations,
        genomes.clone(),
        cell_values,
    );

    // ---- Rescue events: matrix-correcting types go into the map; pseudogenes
    // stay out of it, exactly as the balanced policy keeps them out ----
    let mut events: Vec<RescueEvent> = Vec::new();
    let mut by_type: BTreeMap<&'static str, usize> = BTreeMap::new();
    let mk = |ci: usize, g: u32, rescue_type: RescueType, rng: &mut Rng| -> RescueEvent {
        let start = rng.below(4_400_000) as u64;
        RescueEvent {
            cluster_idx: ClusterIdx(ci),
            genome_idx: GenomeIdx(g as usize),
            rescue_type,
            identity: 95.0 + rng.f64() * 5.0,
            coverage: 85.0 + rng.f64() * 15.0,
            contig: format!("{}_contig{}", short_ids[g as usize], 1 + rng.below(120)),
            start,
            end: start + 300 + rng.below(2400) as u64,
            strand: if rng.chance(0.5) { '+' } else { '-' },
            annotation_match: None,
            confidence: 0.7 + rng.f64() * 0.3,
            synteny_ok: Some(true),
            synteny_evidence: SyntenyEvidence {
                score: Some(0.8 + rng.f64() * 0.2),
                support: 10,
                reference_contexts: 12,
                n_neighbors: 5,
                matched_upstream: 5,
                matched_downstream: 4,
            },
            fragment_evidence: crate::rescue::FragmentEvidence::default(),
            query_start: None,
        }
    };
    for (ci, draft) in drafts.iter().enumerate() {
        for &g in &draft.rescued {
            let r = rng.f64();
            let (label, rescue_type) = if r < 0.35 {
                ("dropped", RescueType::Dropped)
            } else if r < 0.60 {
                ("missing_cds", RescueType::MissingCds)
            } else if r < 0.70 {
                (
                    "split",
                    RescueType::Split {
                        other_cluster: format!("group_{}", rng.range(1, group)),
                    },
                )
            } else if r < 0.85 {
                (
                    "fragmented",
                    RescueType::Fragmented {
                        n_fragments: 2,
                        total_coverage: 88.0,
                    },
                )
            } else if r < 0.90 {
                (
                    "alignment_split",
                    RescueType::AlignmentSplit {
                        n_fragments: 2,
                        total_coverage: 91.0,
                    },
                )
            } else {
                (
                    "start_variant",
                    RescueType::StartVariant { length_diff_bp: 12 },
                )
            };
            *by_type.entry(label).or_default() += 1;
            events.push(mk(ci, g, rescue_type, &mut rng));
        }
    }
    let matrix_events = events.len();
    for (ci, draft) in drafts.iter().enumerate() {
        for &g in &draft.pseudogene {
            *by_type.entry("pseudogene").or_default() += 1;
            events.push(mk(
                ci,
                g,
                RescueType::Pseudogene {
                    reason: PseudogeneReason::Frameshift,
                },
                &mut rng,
            ));
        }
    }
    let rescue_map: HashMap<(ClusterIdx, GenomeIdx), &RescueEvent> = events[..matrix_events]
        .iter()
        .map(|r| ((r.cluster_idx, r.genome_idx), r))
        .collect();

    let nsc = if nsc_pairs.is_empty() {
        None
    } else {
        let mut merge_map = HashMap::new();
        let mut merged_clusters = Vec::new();
        for &(member, rep) in &nsc_pairs {
            merge_map.insert(ClusterIdx(member), ClusterIdx(rep));
            merged_clusters.push(crate::nscluster::NsCluster {
                members: vec![ClusterIdx(rep), ClusterIdx(member)],
                representative: drafts[rep].name.clone(),
            });
        }
        Some(crate::nscluster::NsClusterResult {
            total_merged_genes: merged_clusters.len(),
            n_clusters: merged_clusters.len(),
            merged_clusters,
            merge_map,
        })
    };
    let nsc_cells = nsc_presence_updates(&pa, &nsc, &rescue_map).len();
    let t_inputs = t_start.elapsed().as_secs_f64();

    // ---- Side files the report reads ----
    let t_side = Instant::now();
    let reference_name = by_lineage[0].first().map(|&gi| genomes[gi].clone());
    let write =
        |name: &str| std::io::BufWriter::new(std::fs::File::create(work.join(name)).unwrap());

    {
        let mut w = write("genome_qc.tsv");
        writeln!(w, "genome\tqc_status").unwrap();
        for g in &genomes {
            let status = match rng.below(40) {
                0 => "fail",
                1..=3 => "review",
                _ => "pass",
            };
            writeln!(w, "{g}\t{status}").unwrap();
        }
    }
    {
        let mut w = write("orthology_evidence.tsv");
        writeln!(w, "matrix_action\torthology_confidence_tier").unwrap();
        for (i, _) in events.iter().enumerate() {
            if i < matrix_events {
                let tier = if rng.chance(0.8) {
                    "high_confidence"
                } else {
                    "moderate"
                };
                writeln!(w, "correct_presence\t{tier}").unwrap();
            } else {
                writeln!(w, "audit_only_quality\treview").unwrap();
            }
        }
        for _ in 0..(events.len() / 20) {
            writeln!(w, "rejected_synteny\tlow").unwrap();
        }
    }
    {
        let mut w = write("evolutionary_events.tsv");
        writeln!(w, "gene_cluster\thighlight_priority\tevolutionary_signal_score\tmolecular_mechanism\tsupport_level\tphylogenetic_pattern\trecommended_wording\tmechanism_detail\treference_event_class\ttarget_reference_cluster\ttarget_reference_gene").unwrap();
        for draft in &drafts {
            if draft.category == Category::NscMember {
                continue;
            }
            let options: &[(&str, &str, &str)] = match draft.category {
                Category::Core if !draft.pseudogene.is_empty() => &[
                    ("pseudogenization", "loss", "review"),
                    (
                        "gene_fragmentation_or_structural_break",
                        "fragmented",
                        "candidate",
                    ),
                ],
                Category::Core => &[(
                    "ancestral_retained_no_disruption",
                    "ancestral_retained",
                    "filter",
                )],
                Category::Deletion => &[
                    (
                        "deletion_loss_candidate",
                        "ancestral_with_loss",
                        "candidate",
                    ),
                    (
                        "reference_or_lineage_loss_candidate",
                        "ancestral_with_loss",
                        "highlight",
                    ),
                ],
                Category::Shell => &[
                    (
                        "clean_insertion_or_acquisition_candidate",
                        "clade_gain",
                        "highlight",
                    ),
                    ("IS_interruption_supported", "homoplastic", "candidate"),
                    (
                        "duplication_or_paralog_expansion",
                        "paralog_complex",
                        "review",
                    ),
                    (
                        "translocation_or_context_shift_candidate",
                        "mixed_context",
                        "review",
                    ),
                ],
                _ => &[
                    ("annotation_or_clustering_artifact", "singleton", "filter"),
                    (
                        "unstable_or_homoplastic_accessory_signal",
                        "blinking",
                        "review",
                    ),
                    ("partial_reference_homology", "noise", "filter"),
                    (
                        "clean_insertion_or_acquisition_candidate",
                        "clade_gain",
                        "candidate",
                    ),
                ],
            };
            let (mechanism, phylo, priority) = options[rng.below(options.len())];
            let support = ["high", "moderate", "low"][rng.below(3)];
            let (event_class, target_cluster, target_gene) =
                if mechanism == "partial_reference_homology" {
                    (
                        "partial_reference_homology",
                        format!("Rv{:04}", rng.range(1, 3900)),
                        format!("Rv{:04}", rng.range(1, 3900)),
                    )
                } else {
                    ("", String::new(), String::new())
                };
            writeln!(
                w,
                "{}\t{}\t{:.3}\t{}\t{}\t{}\tReport as {} after checking the flanks.\trescue_evidence:{}\t{}\t{}\t{}",
                draft.name,
                priority,
                rng.f64(),
                mechanism,
                support,
                phylo,
                mechanism.replace('_', " "),
                mechanism,
                event_class,
                target_cluster,
                target_gene
            )
            .unwrap();
        }
    }
    if let Some(reference) = &reference_name {
        let mut w = write("reference_genome_audit.tsv");
        writeln!(w, "gene_cluster\treference_genome\treference_absence_class\tfinal_publication_status\treference_status\tmatrix_action\tevidence_state\tevidence_strength\tpublication_confidence\trecommended_wording").unwrap();
        for draft in &drafts {
            writeln!(w, "{}\t{reference}\treference_like_or_present\treference_matrix_present\treference_matrix\talready_present\tpresent_in_matrix\thigh\tpublication_grade\tConsistent with corrected matrix call.", draft.name).unwrap();
        }
    }
    let mut synteny_rows = 0usize;
    let mut synteny_bytes = 0u64;
    if synteny_mode != "none" {
        let mut w = write("cluster_synteny_blocks.tsv");
        writeln!(w, "gene_cluster\tblock_rank\tn_genomes\tfraction_of_present\tupstream_block_5\tdownstream_block_5\tblock_interpretation\texample_genomes\tcarrier_genomes\trepresentative_genome\trepresentative_gene\trepresentative_contig\trepresentative_start\trepresentative_end\trepresentative_strand\trepresentative_length_bp\trepresentative_block_genes\treference_block_genes").unwrap();
        let n = drafts.len();
        for (ci, draft) in drafts.iter().enumerate() {
            if draft.category == Category::NscMember {
                continue;
            }
            let carriers = draft.present.len();
            let n_blocks = if synteny_mode == "full" {
                // A real run writes one row per distinct context; with draft
                // assemblies a core gene near contig ends has dozens.
                if carriers >= 50 {
                    12
                } else if carriers >= 5 {
                    rng.range(2, 4)
                } else {
                    1
                }
            } else {
                if !rng.chance(0.25) {
                    continue;
                }
                let mut b = 1;
                while b < 4 && rng.chance(0.5) {
                    b += 1;
                }
                b
            }
            .min(carriers);
            // Partition the carriers: a dominant block and a geometric tail.
            let mut shuffled = draft.present.clone();
            rng.shuffle(&mut shuffled);
            let mut sizes = vec![1usize; n_blocks];
            let mut rest = carriers - n_blocks;
            for size in sizes.iter_mut() {
                let take = if rest == 0 {
                    0
                } else {
                    ((rest as f64) * 0.8).ceil() as usize
                };
                *size += take.min(rest);
                rest -= take.min(rest);
            }
            sizes[0] += rest;
            let neighbour = |offset: i64| -> String {
                let idx = (ci as i64 + offset).rem_euclid(n as i64) as usize;
                cluster_names[idx].clone()
            };
            let mut cursor = 0;
            for (rank, size) in sizes.iter().enumerate() {
                let members = &shuffled[cursor..cursor + size];
                cursor += size;
                let rep = members[0] as usize;
                let rep_contig = format!("{}_contig{}", short_ids[rep], 1 + rng.below(120));
                let start = 1000 + rng.below(4_400_000) as u64;
                let len = 300 + rng.below(2400) as u64;
                let mut genes = Vec::with_capacity(11);
                let mut pos = start.saturating_sub(5 * 1100);
                for offset in -5i64..=5 {
                    let (cluster, gene) = if offset == 0 {
                        (draft.name.clone(), draft.name.clone())
                    } else if rank > 0 && (offset + rank as i64) % 4 == 0 {
                        // Alternative contexts swap or lose a neighbour.
                        (
                            "NA".to_string(),
                            format!("{}_{:05}", short_ids[rep], rng.below(99_999)),
                        )
                    } else {
                        let c = neighbour(offset);
                        (c.clone(), c)
                    };
                    let glen = if offset == 0 {
                        len
                    } else {
                        300 + rng.below(1800) as u64
                    };
                    genes.push(format!(
                        "{cluster},{gene},+,{rep_contig},{pos},{},{glen}",
                        pos + glen
                    ));
                    pos += glen + 120;
                }
                let reference_genes: Vec<String> = if reference_name.is_some() {
                    let mut p = 1_500_000u64 + (ci as u64 % 3000) * 1000;
                    (-5i64..=5)
                        .map(|offset| {
                            let c = if offset == 0 {
                                draft.name.clone()
                            } else {
                                neighbour(offset)
                            };
                            let s = format!("{c},{c},+,AS00001_contig1,{p},{},900", p + 900);
                            p += 1020;
                            s
                        })
                        .collect()
                } else {
                    Vec::new()
                };
                let up: Vec<String> = (1..=5).map(|o| neighbour(-o)).collect();
                let down: Vec<String> = (1..=5).map(neighbour).collect();
                let interpretation = if n_blocks == 1 {
                    "single_synteny_block"
                } else if rank == 0 {
                    "dominant_block_minor_alternatives"
                } else {
                    "alternative_synteny_block"
                };
                let examples: Vec<String> = members
                    .iter()
                    .take(5)
                    .map(|&g| {
                        format!(
                            "{}:{}_{:05}",
                            genomes[g as usize],
                            short_ids[g as usize],
                            ci % 100_000
                        )
                    })
                    .collect();
                let carrier_list: Vec<&str> = members
                    .iter()
                    .map(|&g| genomes[g as usize].as_str())
                    .collect();
                let rep_gene = format!("{}_{:05}", short_ids[rep], ci % 100_000);
                let line = format!(
                    "{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t+\t{}\t{}\t{}\n",
                    draft.name,
                    rank + 1,
                    size,
                    *size as f64 / carriers as f64,
                    up.join("|"),
                    down.join("|"),
                    interpretation,
                    examples.join(";"),
                    carrier_list.join(";"),
                    genomes[rep],
                    rep_gene,
                    rep_contig,
                    start,
                    start + len,
                    len,
                    genes.join(";"),
                    reference_genes.join(";"),
                );
                synteny_bytes += line.len() as u64;
                w.write_all(line.as_bytes()).unwrap();
                synteny_rows += 1;
            }
        }
    }
    let t_side_files = t_side.elapsed().as_secs_f64();

    // ---- The report ----
    let tree_path = work.join("tree.nwk");
    std::fs::write(&tree_path, b"(placeholder);\n").unwrap();
    let args = Args {
        presence_absence: "pa.csv".into(),
        pan_reference: "pan_ref.fa".into(),
        genomes: "genomes".into(),
        annotations: None,
        output: work.clone(),
        min_identity: 90.0,
        min_coverage: 80.0,
        threads: 4,
        preset: "asm5".into(),
        min_prevalence: 0.1,
        multi_allele: false,
        synteny: synteny_mode != "none",
        write_corrected_gff: false,
        write_html_report: true,
        matrix_policy: MatrixPolicy::Balanced,
        reference_genome: reference_name.clone(),
        tree: Some(tree_path),
        mobile_element_reference: None,
        mobile_element_min_identity: 85.0,
        mobile_element_min_coverage: 50.0,
        mobile_element_window: 1000,
        nscluster: nsc.is_some(),
        nsc_kmer_size: 21,
        nsc_threshold: 0.8,
        nsc_engine: "minhash".into(),
        verbose: false,
        nsc_min_seq_id: 0.9,
        nsc_coverage: 0.8,
    };
    let count = |label: &str| by_type.get(label).copied().unwrap_or(0);
    let gaps = n_genomes * drafts.len() - present_cells;
    let stats = SummaryStats {
        core_before: pa.core_size(),
        total_gaps: gaps,
        filtered_gaps: gaps / 3,
        skipped_by_prevalence: gaps - gaps / 3,
        searchable_gaps: gaps / 3,
        missing_fasta_gaps: 0,
        skipped_no_reference: 0,
        remaining_gaps: gaps / 3 - rescue_map.len().min(gaps / 3),
        dropped: count("dropped"),
        missing_cds: count("missing_cds"),
        split: count("split"),
        fragmented: count("fragmented"),
        alignment_split: count("alignment_split"),
        start_variant: count("start_variant"),
        paralog: 0,
        pseudogene: count("pseudogene"),
        synteny_rejected: events.len() / 20,
        paralog_rejected: 0,
        matrix_policy: MatrixPolicy::Balanced,
        matrix_corrected: rescue_map.len(),
        nsc_matrix_corrected: nsc_cells,
        synteny_enabled: synteny_mode != "none",
        nsc_clusters: nsc_pairs.len(),
        nsc_merged_genes: nsc_pairs.len(),
    };
    // Evidence rows the way the pipeline builds them, so the per-cell
    // evidence panel and the kept-out cells are exercised at this size too.
    let no_annotations: HashMap<GenomeIdx, HashMap<String, ClusterIdx>> = HashMap::new();
    let prevalences: Vec<f64> = (0..pa.n_clusters())
        .map(|ci| pa.prevalence(ClusterIdx(ci)))
        .collect();
    let seq_lengths = vec![1000usize; pa.n_clusters()];
    let evidence = build_orthology_evidence(OrthologyEvidenceInputs {
        rescues: &events,
        synteny_rejected: &[],
        paralog_rejected: &[],
        pa: &pa,
        gene_to_cluster_maps: &no_annotations,
        cluster_prevalences: &prevalences,
        cluster_seq_lengths: &seq_lengths,
        matrix_policy: MatrixPolicy::Balanced,
    });
    // Clusters under the prevalence floor go unsearched, as in a real run, so
    // the not-searched cells are drawn at this size too. The synthetic rescues
    // ignore the floor, so a cluster with any stays searched.
    let rescued_clusters: ahash::HashSet<usize> =
        events.iter().map(|e| e.cluster_idx.val()).collect();
    let coverage = SearchCoverage {
        clusters: prevalences
            .iter()
            .enumerate()
            .map(|(ci, &p)| {
                (p < args.min_prevalence && !rescued_clusters.contains(&ci))
                    .then_some(NotSearched::LowPrevalence)
            })
            .collect(),
        genomes_without_fasta: Vec::new(),
    };
    let rss_before = sampler.take();
    let t_report = Instant::now();
    write_html_report(
        &args,
        &pa,
        &stats,
        &rescue_map,
        &nsc,
        Some(&tree),
        &evidence,
        &coverage,
        &work,
    )
    .unwrap();
    let t_report = t_report.elapsed().as_secs_f64();
    let rss_after = sampler.take();

    let produced = work.join("recore_report.html");
    let html_bytes = std::fs::metadata(&produced).unwrap().len();
    std::fs::rename(&produced, &out_html).unwrap();
    if std::env::var_os("RECORE_SCALE_KEEP_WORK").is_none() {
        std::fs::remove_dir_all(&work).ok();
    }
    let core = drafts
        .iter()
        .filter(|d| d.category == Category::Core)
        .count();
    let deletion = drafts
        .iter()
        .filter(|d| d.category == Category::Deletion)
        .count();
    let shell = drafts
        .iter()
        .filter(|d| d.category == Category::Shell)
        .count();
    let cloud = drafts
        .iter()
        .filter(|d| d.category == Category::Cloud)
        .count();
    println!(
        "SCALE {{\"genomes\":{n_genomes},\"clusters\":{},\"rendered_clusters\":{},\"tree\":\"{tree_mode}\",\"synteny\":\"{synteny_mode}\",\"present_cells\":{present_cells},\"rescued_cells\":{},\"pseudogene_events\":{},\"nsc_cells\":{nsc_cells},\"nsc_pairs\":{},\"core\":{core},\"deletion\":{deletion},\"shell\":{shell},\"cloud\":{cloud},\"synteny_rows\":{synteny_rows},\"synteny_tsv_bytes\":{synteny_bytes},\"inputs_s\":{t_inputs:.2},\"side_files_s\":{t_side_files:.2},\"report_s\":{t_report:.2},\"peak_rss_inputs\":{rss_before},\"peak_rss_during_report\":{rss_after},\"html_bytes\":{html_bytes},\"out\":\"{}\"}}",
        drafts.len(),
        drafts.len() - nsc_pairs.len(),
        rescue_map.len(),
        count("pseudogene"),
        nsc_pairs.len(),
        out_html.display()
    );
}
