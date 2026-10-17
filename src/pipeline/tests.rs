use super::*;
use crate::parser::PresenceAbsence;

fn default_args_for_validation() -> Args {
    Args {
        presence_absence: "pa.csv".into(),
        pan_reference: "pan_ref.fa".into(),
        genomes: "genomes".into(),
        annotations: None,
        output: "recore_output".into(),
        min_identity: 90.0,
        min_coverage: 80.0,
        threads: 4,
        preset: "asm5".into(),
        min_prevalence: 0.1,
        multi_allele: false,
        synteny: false,
        write_corrected_gff: false,
        write_html_report: false,
        matrix_policy: MatrixPolicy::Balanced,
        reference_genome: None,
        tree: None,
        mobile_element_reference: None,
        mobile_element_min_identity: 85.0,
        mobile_element_min_coverage: 50.0,
        mobile_element_window: 1000,
        nscluster: false,
        nsc_kmer_size: 21,
        nsc_threshold: 0.8,
        nsc_engine: "minhash".into(),
        verbose: false,
        nsc_min_seq_id: 0.9,
        nsc_coverage: 0.8,
    }
}

fn rescue_with_type(rescue_type: RescueType, synteny_ok: Option<bool>) -> RescueEvent {
    RescueEvent {
        cluster_idx: ClusterIdx(0),
        genome_idx: GenomeIdx(0),
        rescue_type,
        identity: 90.0,
        coverage: 90.0,
        contig: "chr".into(),
        start: 0,
        end: 100,
        strand: '+',
        annotation_match: None,
        confidence: 0.8,
        synteny_ok,
        synteny_evidence: SyntenyEvidence::default(),
        fragment_evidence: crate::rescue::FragmentEvidence::default(),
    }
}

#[test]
fn invalid_args_are_rejected_before_running_pipeline() {
    let mut args = default_args_for_validation();
    args.threads = 0;
    assert!(validate_args(&args).is_err());

    let mut args = default_args_for_validation();
    args.min_identity = 101.0;
    assert!(validate_args(&args).is_err());

    let mut args = default_args_for_validation();
    args.min_prevalence = -0.1;
    assert!(validate_args(&args).is_err());

    let mut args = default_args_for_validation();
    args.nsc_kmer_size = 0;
    assert!(validate_args(&args).is_err());

    let mut args = default_args_for_validation();
    args.nsc_engine = "minhsah".into();
    assert!(validate_args(&args).is_err());
}

/// Real input paths in a temporary directory, for the path checks.
fn validation_fixture() -> (tempfile::TempDir, Args) {
    let dir = tempfile::tempdir().unwrap();
    std::fs::write(dir.path().join("pa.csv"), "Gene\n").unwrap();
    std::fs::write(dir.path().join("pan_ref.fa"), ">g\nACGT\n").unwrap();
    std::fs::create_dir(dir.path().join("genomes")).unwrap();
    let mut args = default_args_for_validation();
    args.presence_absence = dir.path().join("pa.csv");
    args.pan_reference = dir.path().join("pan_ref.fa");
    args.genomes = dir.path().join("genomes");
    args.output = dir.path().join("out");
    (dir, args)
}

#[test]
fn default_args_pass_validation() {
    let (_dir, args) = validation_fixture();
    validate_args(&args).unwrap();
    assert!(
        args.output.is_dir(),
        "the output directory is created and probed before the search"
    );
}

/// Every one of these used to be discovered late or not at all: a misspelled
/// --annotations ran to exit 0 without annotation logic, --synteny without -a
/// was silently ignored, and a bad --mobile-element-reference failed after the
/// whole search, having already deleted the previous run's outputs.
#[test]
fn input_paths_are_checked_before_any_work() {
    let message = |args: &Args| format!("{:#}", validate_args(args).unwrap_err());

    let (dir, mut args) = validation_fixture();
    args.annotations = Some(dir.path().join("annotationz"));
    assert!(message(&args).contains("-a/--annotations"));

    let (_dir, mut args) = validation_fixture();
    args.synteny = true;
    assert!(message(&args).contains("--synteny needs --annotations"));

    let (dir, mut args) = validation_fixture();
    args.mobile_element_reference = Some(dir.path().join("missing_IS.fa"));
    assert!(message(&args).contains("--mobile-element-reference"));

    let (dir, mut args) = validation_fixture();
    args.tree = Some(dir.path().join("missing.nwk"));
    assert!(message(&args).contains("--tree"));

    let (_dir, mut args) = validation_fixture();
    args.genomes = args.presence_absence.clone();
    assert!(message(&args).contains("is not a directory"));

    let (dir, mut args) = validation_fixture();
    std::fs::write(dir.path().join("afile"), "x").unwrap();
    args.output = dir.path().join("afile");
    assert!(message(&args).contains("exists and is not a directory"));
}

#[test]
fn reference_genome_must_match_exactly_one_column() {
    let pa = PresenceAbsence::new_for_test(
        vec!["c".into()],
        vec![
            "L4_AS00401".into(),
            "L2_AS00401".into(),
            "L1_AS00101".into(),
        ],
        vec![vec!["x".into(), "".into(), "".into()]],
    );
    check_reference_genome(&pa, "L1_AS00101").unwrap();
    // A unique AS identifier still resolves, with a warning.
    check_reference_genome(&pa, "AS00101").unwrap();

    let shared = format!("{:#}", check_reference_genome(&pa, "AS00401").unwrap_err());
    assert!(shared.contains("several"), "{shared}");
    assert!(
        shared.contains("L4_AS00401") && shared.contains("L2_AS00401"),
        "{shared}"
    );

    let unknown = format!("{:#}", check_reference_genome(&pa, "H37Rv").unwrap_err());
    assert!(unknown.contains("matches no column"), "{unknown}");
}

#[test]
fn rescue_orf_signal_flags_disrupted_and_frame_review_spans() {
    let mut rescue = rescue_with_type(RescueType::Dropped, Some(true));
    rescue.start = 0;
    rescue.end = 99;
    assert_eq!(rescue_orf_signal(&rescue), "in_frame_span");

    rescue.end = 100;
    assert_eq!(rescue_orf_signal(&rescue), "out_of_frame_span_review");

    rescue.rescue_type = RescueType::Pseudogene {
        reason: PseudogeneReason::Frameshift,
    };
    assert_eq!(rescue_orf_signal(&rescue), "disrupted_orf");
}

#[test]
fn corrected_csv_recomputes_multicopy_metadata_counts() {
    let mut input = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        input,
        "Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,Avg sequences per isolate,g1,g2,g3"
    )
    .unwrap();
    writeln!(input, "clA,clA,prod,2,3,1.5,a1;a2,b1,").unwrap();
    let pa = parse_presence_absence(input.path()).unwrap();
    let rescue = RescueEvent {
        cluster_idx: ClusterIdx(0),
        genome_idx: GenomeIdx(2),
        contig: "ctg".into(),
        start: 9,
        end: 99,
        ..rescue_with_type(RescueType::Dropped, Some(true))
    };
    let mut rescue_map = HashMap::new();
    rescue_map.insert((rescue.cluster_idx, rescue.genome_idx), &rescue);
    let out = tempfile::tempdir().unwrap();

    write_corrected_csv(&pa, &rescue_map, out.path(), &None).unwrap();

    let mut reader =
        csv::Reader::from_path(out.path().join("gene_presence_absence_recore.csv")).unwrap();
    let row = reader.records().next().unwrap().unwrap();
    let headers = reader.headers().unwrap().clone();
    let idx = |name: &str| headers.iter().position(|h| h == name).unwrap();
    assert_eq!(&row[idx("No. isolates")], "3");
    assert_eq!(&row[idx("No. sequences")], "4");
    assert_eq!(&row[idx("Avg sequences per isolate")], "1.333");
}

#[test]
fn html_report_summarizes_run_outputs() {
    let out = tempfile::tempdir().unwrap();
    std::fs::write(
        out.path().join("orthology_evidence.tsv"),
        b"matrix_action\torthology_confidence_tier\ncorrect_presence\thigh_confidence\naudit_only_quality\treview\n",
    )
    .unwrap();
    std::fs::write(
        out.path().join("evolutionary_events.tsv"),
        b"gene_cluster\thighlight_priority\tevolutionary_signal_score\tmolecular_mechanism\tsupport_level\tphylogenetic_pattern\trecommended_wording\ngroup_1\thighlight\t0.910\tclade_associated_gain\thigh\tclade_gain\tDescribe as clade-associated.\n",
    )
    .unwrap();
    std::fs::write(
        out.path().join("genome_qc.tsv"),
        b"genome\tqc_status\ng1\tpass\ng2\treview\n",
    )
    .unwrap();
    let pa = PresenceAbsence::new_for_test(
        vec!["group_1".into()],
        vec!["g1".into(), "g2".into()],
        vec![vec!["gene1".into(), "".into()]],
    );
    let mut args = default_args_for_validation();
    args.output = out.path().to_path_buf();
    args.write_html_report = true;
    args.nscluster = true;
    let stats = SummaryStats {
        core_before: 0,
        total_gaps: 1,
        filtered_gaps: 1,
        skipped_by_prevalence: 0,
        searchable_gaps: 1,
        missing_fasta_gaps: 0,
        dropped: 1,
        missing_cds: 0,
        split: 0,
        fragmented: 0,
        alignment_split: 0,
        start_variant: 0,
        paralog: 0,
        pseudogene: 0,
        synteny_rejected: 0,
        paralog_rejected: 0,
        matrix_policy: MatrixPolicy::Balanced,
        matrix_corrected: 0,
        nsc_matrix_corrected: 0,
        synteny_enabled: false,
        nsc_clusters: 0,
        nsc_merged_genes: 0,
        skipped_no_reference: 0,
        remaining_gaps: 0,
    };

    write_html_report(
        &args,
        &pa,
        &stats,
        &HashMap::new(),
        &None,
        None,
        &[],
        out.path(),
    )
    .unwrap();

    let text = std::fs::read_to_string(out.path().join("recore_report.html")).unwrap();
    assert!(text.contains("<title>Recore run report</title>"));
    assert!(text.contains("group_1"));
    assert!(text.contains("Orthology actions"));
    assert!(text.contains("gene_presence_absence_recore.csv"));
    // Without a tree the workspace still opens on a star tree in matrix order,
    // and says so, instead of leaving a page with one sentence on it.
    assert!(text.contains("No --tree was given"));
    assert!(text.contains("id='workspaceNotices'"));
    for control in [
        "id='clusterOmnibox'",
        "id='prevCluster'",
        "id='nextCluster'",
        "id='toggleSynteny'",
        "id='exportView'",
        "id='palette'",
    ] {
        assert!(text.contains(control), "tree-less report lacks {control}");
    }
    assert!(text.contains("window.recoreExplorerFailed"));
    // What does work without a tree stays.
    assert!(text.contains("id='themeToggle'"));
    assert!(text.contains("id='openDetails'"));
    assert_no_remote_assets(&text);
}

/// The report is opened from result directories on offline machines, so it must
/// not pull stylesheets, fonts, scripts, or images over the network. Navigational
/// links to the project are fine: they only resolve if the reader clicks them.
fn assert_no_remote_assets(html: &str) {
    for tag in ["<link", "<script", "<img", "<iframe", "<use"] {
        let mut rest = html;
        while let Some(at) = rest.find(tag) {
            rest = &rest[at + tag.len()..];
            let end = rest.find('>').unwrap_or(rest.len());
            let element = &rest[..end];
            assert!(
                !element.contains("http://") && !element.contains("https://"),
                "report loads a remote resource: {tag}{element}>"
            );
        }
    }
    assert!(
        !html.contains("@import"),
        "report pulls a stylesheet with @import"
    );
    assert!(
        !html.contains("url(http"),
        "report references a remote url() in CSS"
    );
}

#[test]
fn html_report_embeds_small_phylogeny_matrix_explorer() {
    let out = tempfile::tempdir().unwrap();
    std::fs::write(
        out.path().join("orthology_evidence.tsv"),
        b"matrix_action\torthology_confidence_tier\ncorrect_presence\thigh_confidence\n",
    )
    .unwrap();
    std::fs::write(
        out.path().join("evolutionary_events.tsv"),
        b"gene_cluster\thighlight_priority\tevolutionary_signal_score\tmolecular_mechanism\tsupport_level\tphylogenetic_pattern\trecommended_wording\ngroup_1\thighlight\t0.910\tclade_associated_gain\thigh\tclade_gain\tDescribe as clade-associated.\n",
    )
    .unwrap();
    std::fs::write(
        out.path().join("genome_qc.tsv"),
        b"genome\tqc_status\ng1\tpass\n",
    )
    .unwrap();
    let pa = PresenceAbsence::new_for_test(
        vec!["group_1".into(), "group_2".into()],
        vec!["g1".into(), "g2".into()],
        vec![
            vec!["gene1".into(), "".into()],
            vec!["".into(), "gene2".into()],
        ],
    );
    let mut args = default_args_for_validation();
    args.output = out.path().to_path_buf();
    args.write_html_report = true;
    args.tree = Some(out.path().join("tree.nwk"));
    let stats = SummaryStats {
        core_before: 0,
        total_gaps: 2,
        filtered_gaps: 2,
        skipped_by_prevalence: 0,
        searchable_gaps: 2,
        missing_fasta_gaps: 0,
        dropped: 0,
        missing_cds: 0,
        split: 0,
        fragmented: 0,
        alignment_split: 0,
        start_variant: 0,
        paralog: 0,
        pseudogene: 0,
        synteny_rejected: 0,
        paralog_rejected: 0,
        matrix_policy: MatrixPolicy::Balanced,
        matrix_corrected: 0,
        nsc_matrix_corrected: 0,
        synteny_enabled: false,
        nsc_clusters: 0,
        nsc_merged_genes: 0,
        skipped_no_reference: 0,
        remaining_gaps: 0,
    };
    let tree = PhyloTree {
        root: PhyloNode {
            label: None,
            genome_idx: None,
            branch_length: None,
            children: vec![
                PhyloNode {
                    label: Some("g1".into()),
                    genome_idx: Some(GenomeIdx(0)),
                    branch_length: None,
                    children: Vec::new(),
                },
                PhyloNode {
                    label: Some("g2".into()),
                    genome_idx: Some(GenomeIdx(1)),
                    branch_length: None,
                    children: Vec::new(),
                },
            ],
        },
    };

    write_html_report(
        &args,
        &pa,
        &stats,
        &HashMap::new(),
        &None,
        Some(&tree),
        &[],
        out.path(),
    )
    .unwrap();

    let text = std::fs::read_to_string(out.path().join("recore_report.html")).unwrap();
    assert!(text.contains("Corrected matrix"));
    assert!(text.contains("window.recoreExplorerData"));
    assert!(text.contains("clusterDetails"));
    assert!(text.contains("columnHeaderCanvas"));
    assert!(text.contains("syntenyPanel"));
    assert!(text.contains("syntenySvg"));
    assert!(text.contains("toggleSynteny"));
    assert!(text.contains("Show synteny"));
    assert!(text.contains("class='workspace synteny-hidden'"));
    assert!(text.contains("applySyntenyVisibility"));
    assert!(text.contains("applyLayout"));
    assert!(text.contains("setupDragHandle"));
    assert!(text.contains("treeHandle"));
    assert!(text.contains("syntenyHandle"));
    assert!(text.contains("syntenyZoom"));
    assert!(text.contains("let syntenyHidden = true"));
    assert!(text.contains("selectedSyntenyReference"));
    assert!(text.contains("id='syntenySvg'"));
    assert!(text.contains("treeEdgeColor"));
    assert!(text.contains("subtreeLineage"));
    assert!(text.contains("\"synteny\""));
    assert!(text.contains("selectedSyntenyBlock"));
    assert!(text.contains("function tidyTree"));
    assert!(text.contains("genomeBits"));
    assert!(text.contains("syntenyCarrierIndices"));
    assert!(text.contains("referenceGenes"));
    assert!(text.contains("representativeLength"));
    assert!(text.contains("H37Rv context"));
    assert!(text.contains("renderSyntenyPanel"));
    assert!(text.contains("cellSize"));
    assert!(text.contains("rowHeight"));
    assert!(text.contains("id='cellSize' type='range' min='0.25'"));
    assert!(text.contains("id='rowHeight' type='range' min='1'"));
    assert!(text.contains("step='0.25'"));
    assert!(text.contains("const headerH = 24"));
    assert!(text.contains("tipColorMode"));
    assert!(text.contains("matrixColorMode"));
    assert!(text.contains("processLegend"));
    assert!(text.contains("lineageLegend"));

    assert!(text.contains("data-mechanism"));
    assert!(text.contains("data-clear-mechanism"));
    assert!(text.contains("mechanismFilters = new Set"));
    assert!(text.contains("mechanismFilters.has"));
    assert!(text.contains("Click mechanism cards to add/remove"));
    assert!(text.contains("referenceCell"));
    assert!(text.contains("targetReferenceGene"));
    assert!(text.contains("Reference (audit)"));
    assert!(text.contains("value='lineage'"));
    assert!(text.contains("#C94434"));
    assert!(text.contains("Colour tips by"));
    assert!(text.contains("Colour cells by"));
    assert!(text.contains("clade_associated_gain"));
    assert!(text.contains("tipLabels"));
    assert!(text.contains("visibleColumnRange"));
    assert!(text.contains("visibleRowRange"));
    assert!(text.contains("matrixSpacer"));
    assert!(text.contains("phylogramMode"));
    assert!(text.contains("const tipX = 14 + usable"));
    // Both the branch layout and the scale bar measure against this.
    assert!(text.contains("function treeUsableWidth"));
    assert!(!text.contains("clusterHeaderLabel"));
    assert!(text.contains("referenceOnly"));
    assert!(text.contains("explorerFullscreen"));
    // The chrome the redesign introduced: one find field over both kinds of
    // thing, per-pane menus, a genome card and a scale for the branch lengths.
    assert!(text.contains("id='clusterOmnibox'"));
    assert!(text.contains("id='paletteList'"));
    assert!(text.contains("id='shortcutSheet'"));
    assert!(text.contains("id='treeMenu'"));
    assert!(text.contains("id='matrixMenu'"));
    assert!(text.contains("id='genomeDetails'"));
    assert!(text.contains("id='treeScaleBar'"));
    assert!(text.contains("id='overviewCanvas'"));
    assert!(text.contains("substitutions/site"));
    assert!(text.contains("Alt-click any tip"));
    assert!(text.contains("\"referenceGenome\":null"));
    // Cells are the compact per-cluster encoding, not one character per genome.
    // Two genomes fall inside a single byte, so both clusters use bitset mode:
    // group_1 is present in genome 0 and group_2 in genome 1.
    assert!(text.contains("\"cellEncoding\":\"recore.cells.v2\""));
    assert!(text.contains("\"cells\":[\"AgEAAAAAAAA=\",\"AgIAAAAAAAA=\"]"));
    assert!(text.contains("\"genomeStats\":[[1,0,0,0,1],[1,0,0,0,1]]"));
    assert!(text.contains("matrixCanvas"));
    // Canvases must be allocated at the device ratio, and wide columns must be named.
    assert!(text.contains("devicePixelRatio"));
    assert!(text.contains("function canvasScale"));
    assert!(text.contains("function columnLabelBand"));
    assert!(text.contains("fitLabel"));
}

/// One row of the sample synteny table:
/// (cluster, representative genome, contig, start, length, carrier genomes,
/// carrier fraction, interpretation, variant index).
type SampleSyntenyBlock = (
    &'static str,
    &'static str,
    &'static str,
    u64,
    u64,
    u64,
    f64,
    &'static str,
    usize,
);

#[test]
fn rich_sample_report_renders_full_explorer() {
    // Builds a realistic MTBC-style toy dataset (16 genomes spanning 12 lineages,
    // 12 clusters with varied mechanisms, synteny blocks, NSC merge, rescue
    // events, reference audit) and renders a full HTML report. When the
    // RECORE_SAMPLE_DUMP env var is set the rendered report is copied there
    // so it can be opened in a browser for visual review.
    let out = tempfile::tempdir().unwrap();

    let genomes: Vec<String> = vec![
        "L4_AS00146",
        "L4_AS00150",
        "L4_AS00200",
        "L2_AS00021",
        "L2_AS00060",
        "L3_AS00094",
        "L1_AS00033",
        "L1_AS00100",
        "L5_AS00161",
        "L6_AS00352",
        "L7_AS00044",
        "L8_AS00022",
        "A1_AS00349",
        "A2_AS00357",
        "A4_AS00245",
        "A4_AS00248",
    ]
    .into_iter()
    .map(String::from)
    .collect();
    let n_g = genomes.len();
    let gi = |name: &str| -> usize { genomes.iter().position(|g| g == name).unwrap() };

    let cluster_specs: Vec<(&str, &str, Vec<&str>)> = vec![
        (
            "group_001",
            "dnaA chromosomal replication initiator",
            vec![
                "L4_AS00146",
                "L4_AS00150",
                "L4_AS00200",
                "L2_AS00021",
                "L2_AS00060",
                "L3_AS00094",
                "L1_AS00033",
                "L1_AS00100",
                "L5_AS00161",
                "L6_AS00352",
                "L7_AS00044",
                "L8_AS00022",
                "A1_AS00349",
                "A2_AS00357",
                "A4_AS00245",
                "A4_AS00248",
            ],
        ),
        (
            "group_002",
            "PE_PGRS33 lineage-restricted surface protein",
            vec![
                "L4_AS00146",
                "L4_AS00150",
                "L4_AS00200",
                "L2_AS00021",
                "L2_AS00060",
            ],
        ),
        (
            "group_003",
            "plcA phospholipase C (RD5 locus)",
            vec![
                "L4_AS00146",
                "L4_AS00150",
                "L4_AS00200",
                "L2_AS00021",
                "L2_AS00060",
                "L3_AS00094",
                "L1_AS00033",
                "L1_AS00100",
                "L5_AS00161",
                "L6_AS00352",
                "L7_AS00044",
                "L8_AS00022",
                "A1_AS00349",
                "A2_AS00357",
            ],
        ),
        (
            "group_004",
            "Rv2820c hypothetical (IS6110 hotspot)",
            vec![
                "L4_AS00146",
                "L4_AS00150",
                "L4_AS00200",
                "L3_AS00094",
                "L1_AS00033",
                "L1_AS00100",
                "L5_AS00161",
                "L6_AS00352",
                "L7_AS00044",
                "L8_AS00022",
                "A1_AS00349",
                "A2_AS00357",
                "A4_AS00245",
                "A4_AS00248",
            ],
        ),
        (
            "group_005",
            "pks15-1 polyketide synthase (degraded)",
            vec![
                "L4_AS00146",
                "L4_AS00150",
                "L4_AS00200",
                "L2_AS00021",
                "L2_AS00060",
                "L3_AS00094",
                "L1_AS00033",
            ],
        ),
        (
            "group_006",
            "RD1 / esxA virulence (lost in BCG-like ref)",
            vec![
                "L2_AS00021",
                "L2_AS00060",
                "L3_AS00094",
                "L1_AS00033",
                "L1_AS00100",
                "L5_AS00161",
                "L6_AS00352",
                "L7_AS00044",
                "L8_AS00022",
                "A1_AS00349",
                "A2_AS00357",
                "A4_AS00245",
                "A4_AS00248",
            ],
        ),
        (
            "group_007",
            "rpsL ribosomal S12 (duplicated)",
            vec![
                "L4_AS00146",
                "L4_AS00150",
                "L4_AS00200",
                "L2_AS00021",
                "L2_AS00060",
                "L3_AS00094",
                "L1_AS00033",
                "L1_AS00100",
                "L5_AS00161",
                "L6_AS00352",
                "L7_AS00044",
                "L8_AS00022",
                "A1_AS00349",
                "A2_AS00357",
                "A4_AS00245",
                "A4_AS00248",
            ],
        ),
        (
            "group_007b",
            "rpsL paralog (NSC-merged)",
            vec!["L4_AS00146", "L4_AS00150", "L1_AS00033", "A4_AS00245"],
        ),
        (
            "group_008",
            "mmpL5 efflux pump (split/fragmented)",
            vec![
                "L4_AS00146",
                "L4_AS00150",
                "L4_AS00200",
                "L2_AS00021",
                "L2_AS00060",
                "L3_AS00094",
                "L1_AS00033",
                "L1_AS00100",
                "L8_AS00022",
            ],
        ),
        (
            "group_009",
            "TB7.7 antigen partial copy",
            vec!["A1_AS00349", "A2_AS00357", "A4_AS00245", "A4_AS00248"],
        ),
        (
            "group_010",
            "Rv0095c translocated paralog",
            vec!["L4_AS00146", "L2_AS00021", "L3_AS00094", "L1_AS00033"],
        ),
        (
            "group_011",
            "katG start-variant annotation",
            vec![
                "L4_AS00146",
                "L4_AS00150",
                "L4_AS00200",
                "L2_AS00021",
                "L2_AS00060",
                "L3_AS00094",
                "L1_AS00033",
                "L1_AS00100",
                "L5_AS00161",
                "L6_AS00352",
                "L7_AS00044",
                "L8_AS00022",
                "A1_AS00349",
                "A2_AS00357",
                "A4_AS00245",
                "A4_AS00248",
            ],
        ),
        (
            "group_012",
            "Rv1986 patchy lineage-specific accessory",
            vec![
                "L4_AS00146",
                "L2_AS00060",
                "L3_AS00094",
                "L7_AS00044",
                "A4_AS00248",
            ],
        ),
    ];

    let cluster_names: Vec<String> = cluster_specs
        .iter()
        .map(|(n, _, _)| (*n).to_string())
        .collect();
    let mut matrix: Vec<Vec<String>> = vec![vec![String::new(); n_g]; cluster_specs.len()];
    for (ci, (_name, _ann, present_names)) in cluster_specs.iter().enumerate() {
        for genome_name in present_names {
            let g = gi(genome_name);
            matrix[ci][g] = format!("{}_{}", cluster_names[ci], genomes[g]);
        }
    }

    let pa = PresenceAbsence::new_for_test(cluster_names.clone(), genomes.clone(), matrix);

    // Rescue events for a few cells that would otherwise be absent.
    let mk_rescue = |cluster: usize,
                     genome: usize,
                     rtype: RescueType,
                     contig: &str,
                     start: u64,
                     end: u64|
     -> RescueEvent {
        RescueEvent {
            cluster_idx: ClusterIdx(cluster),
            genome_idx: GenomeIdx(genome),
            rescue_type: rtype,
            identity: 96.5,
            coverage: 92.0,
            contig: contig.into(),
            start,
            end,
            strand: '+',
            annotation_match: Some("RecoreRescue".into()),
            confidence: 0.88,
            synteny_ok: Some(true),
            synteny_evidence: SyntenyEvidence {
                score: Some(0.91),
                support: 12,
                reference_contexts: 14,
                n_neighbors: 5,
                matched_upstream: 5,
                matched_downstream: 5,
            },
            fragment_evidence: crate::rescue::FragmentEvidence::default(),
        }
    };

    #[rustfmt::skip]
    let rescues = [
        mk_rescue(3, gi("L2_AS00021"), RescueType::Dropped, "AS00021", 819_540, 820_660),
        mk_rescue(3, gi("L2_AS00060"), RescueType::Dropped, "AS00060", 821_010, 822_130),
        mk_rescue(4, gi("L1_AS00033"), RescueType::MissingCds, "AS00033", 3_120_440, 3_121_310),
        mk_rescue(5, gi("L4_AS00200"), RescueType::Pseudogene { reason: PseudogeneReason::PartialDeletion { coverage: 0.71 } }, "AS00200", 778_900, 779_510),
        mk_rescue(7, gi("L1_AS00100"), RescueType::Split { other_cluster: "group_007b".into() }, "AS00100", 4_120_010, 4_121_220),
        mk_rescue(8, gi("L5_AS00161"), RescueType::Fragmented { n_fragments: 2, total_coverage: 84.5 }, "AS00161", 2_544_910, 2_546_330),
    ];
    let rescue_map: HashMap<(ClusterIdx, GenomeIdx), &RescueEvent> = rescues
        .iter()
        .map(|r| ((r.cluster_idx, r.genome_idx), r))
        .collect();

    // NSC: merge group_007b (idx 7) into group_007 (idx 6) as paralog redundancy.
    let nsc = crate::nscluster::NsClusterResult {
        merged_clusters: vec![crate::nscluster::NsCluster {
            members: vec![ClusterIdx(6), ClusterIdx(7)],
            representative: "group_007".into(),
        }],
        total_merged_genes: 2,
        n_clusters: 1,
        merge_map: {
            let mut m = HashMap::new();
            m.insert(ClusterIdx(7), ClusterIdx(6));
            m
        },
    };

    // MTBC-shaped phylogeny: animal clade + L8/L7/L5/L6 + L1 + L3 + L2 + L4.
    let leaf = |name: &str| PhyloNode {
        label: Some(name.into()),
        genome_idx: Some(GenomeIdx(gi(name))),
        branch_length: Some(0.002),
        children: Vec::new(),
    };
    let clade = |bl: f64, children: Vec<PhyloNode>| PhyloNode {
        label: None,
        genome_idx: None,
        branch_length: Some(bl),
        children,
    };
    let l4 = clade(
        0.001,
        vec![leaf("L4_AS00146"), leaf("L4_AS00150"), leaf("L4_AS00200")],
    );
    let l2 = clade(0.001, vec![leaf("L2_AS00021"), leaf("L2_AS00060")]);
    let l3 = leaf("L3_AS00094");
    let l1 = clade(0.0012, vec![leaf("L1_AS00033"), leaf("L1_AS00100")]);
    let modern = clade(
        0.002,
        vec![clade(0.001, vec![l4, l2]), clade(0.0015, vec![l3, l1])],
    );
    let africanum = clade(
        0.002,
        vec![
            leaf("L5_AS00161"),
            leaf("L6_AS00352"),
            leaf("L7_AS00044"),
            leaf("L8_AS00022"),
        ],
    );
    let animals = clade(
        0.0035,
        vec![
            leaf("A2_AS00357"),
            leaf("A4_AS00245"),
            leaf("A4_AS00248"),
            leaf("A1_AS00349"),
        ],
    );
    let ancient = clade(0.003, vec![africanum, animals]);
    let tree = PhyloTree {
        root: PhyloNode {
            label: None,
            genome_idx: None,
            branch_length: None,
            children: vec![modern, ancient],
        },
    };

    // ── TSV side-files ──
    // orthology evidence: spread across actions and tiers.
    std::fs::write(
        out.path().join("orthology_evidence.tsv"),
        "matrix_action\torthology_confidence_tier\n\
correct_presence\thigh_confidence\n\
correct_presence\thigh_confidence\n\
correct_presence\tmoderate\n\
audit_only_quality\treview\n\
audit_only_quality\treview\n\
audit_only_quality\treview\n\
rejected_paralog\tlow\n\
rejected_synteny\tlow\n\
correct_presence\thigh_confidence\n\
correct_presence\thigh_confidence\n",
    )
    .unwrap();

    // genome QC: mostly pass with a few reviews/fails.
    let mut qc = String::from("genome\tqc_status\n");
    for (i, name) in genomes.iter().enumerate() {
        let status = match i % 7 {
            6 => "fail",
            5 | 3 => "review",
            _ => "pass",
        };
        qc.push_str(&format!("{}\t{}\n", name, status));
    }
    std::fs::write(out.path().join("genome_qc.tsv"), qc).unwrap();

    // evolutionary events: one row per cluster with realistic mechanism distribution.
    let mech_rows: Vec<(&str, &str, f64, &str, &str, &str, &str)> = vec![
        ("group_001", "highlight", 0.98, "ancestral_retained_no_disruption", "high",  "ancestral_retained",      "Conserved core gene retained across the analysed tree."),
        ("group_002", "highlight", 0.91, "clean_insertion_or_acquisition_candidate", "high", "clade_gain",       "Reference-absent locus concentrated in a modern clade; validate flanks before claiming novelty."),
        ("group_003", "candidate", 0.85, "deletion_loss_candidate",          "moderate","ancestral_with_loss",   "Lost in the A4/M. bovis clade; consistent with RD5-style deletion."),
        ("group_004", "candidate", 0.72, "IS_interruption_supported",         "moderate","homoplastic",          "Disrupted near IS6110 in two L2 isolates; treat as interrupted, not absent."),
        ("group_005", "review",    0.55, "pseudogenization",                  "low",     "loss",                 "Frame-disrupting evidence; report as pseudogene rather than presence."),
        ("group_006", "highlight", 0.93, "reference_or_lineage_loss_candidate","high",   "ancestral_with_loss",  "Background suggests presence; reference appears to have lost the locus."),
        ("group_007", "review",    0.40, "duplication_or_paralog_expansion",  "low",     "paralog_complex",      "Paralog/duplication signal; group_007b merged via NSC."),
        ("group_008", "candidate", 0.66, "gene_fragmentation_or_structural_break","moderate","fragmented",       "Split/fragmented alignment; inspect breakpoints before treating as present."),
        ("group_009", "filter",    0.22, "partial_reference_homology",        "low",     "noise",                "Short partial match to a reference locus; not a novel gene."),
        ("group_010", "review",    0.48, "translocation_or_context_shift_candidate","low","mixed_context",       "Sequence in incompatible context; possible translocation or paralog."),
        ("group_011", "filter",    0.18, "annotation_or_clustering_artifact", "low",     "singleton",            "Likely annotation/clustering boundary; not a biological event."),
        ("group_012", "review",    0.51, "unstable_or_homoplastic_accessory_signal","low","blinking",            "Patchy presence across distant lineages; deprioritise for headlines."),
    ];
    let mut evo = String::from(
        "gene_cluster\thighlight_priority\tevolutionary_signal_score\tmolecular_mechanism\tsupport_level\tphylogenetic_pattern\trecommended_wording\tmechanism_detail\n",
    );
    for (cluster, priority, score, mechanism, support, phylo, wording) in &mech_rows {
        evo.push_str(&format!(
            "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\trescue_evidence:{}\n",
            cluster, priority, score, mechanism, support, phylo, wording, mechanism,
        ));
    }
    std::fs::write(out.path().join("evolutionary_events.tsv"), evo).unwrap();

    // reference audit: each cluster vs L4_AS00146 (H37Rv stand-in).
    let mut ref_audit = String::from(
        "gene_cluster\treference_genome\treference_absence_class\tfinal_publication_status\treference_status\tmatrix_action\tevidence_state\tevidence_strength\tpublication_confidence\trecommended_wording\trescue_decision\trescue_type\trejection_reason\tidentity\tcoverage\tcontig\tstart\tend\tstrand\tsynteny\tannotation_match\tannotation_cluster\tfragment_note\n",
    );
    let ref_specs: Vec<(&str, &str, &str, &str, &str)> = vec![
        (
            "group_001",
            "reference_like_or_present",
            "reference_matrix_present",
            "already_present",
            "STRICT_PRESENT",
        ),
        (
            "group_002",
            "reference_absent_partial",
            "reference_absent_clean",
            "absent_in_reference",
            "ABSENT",
        ),
        (
            "group_003",
            "reference_like_or_present",
            "reference_matrix_present",
            "already_present",
            "STRICT_PRESENT",
        ),
        (
            "group_004",
            "reference_like_or_present",
            "reference_matrix_present",
            "already_present",
            "STRICT_PRESENT",
        ),
        (
            "group_005",
            "reference_like_or_present",
            "reference_matrix_present",
            "already_present",
            "STRICT_PRESENT",
        ),
        (
            "group_006",
            "reference_absent_clean",
            "reference_absent_in_lineage",
            "absent_in_reference",
            "ABSENT",
        ),
        (
            "group_007",
            "reference_like_or_present",
            "reference_matrix_present",
            "already_present",
            "STRICT_PRESENT",
        ),
        (
            "group_008",
            "reference_like_or_present",
            "reference_matrix_present",
            "already_present",
            "STRICT_PRESENT",
        ),
        (
            "group_009",
            "reference_like_partial",
            "reference_like_ambiguous",
            "audit_only_quality",
            "PARTIAL_HOMOLOGY",
        ),
        (
            "group_010",
            "reference_like_or_present",
            "reference_matrix_present",
            "already_present",
            "STRICT_PRESENT",
        ),
        (
            "group_011",
            "reference_like_or_present",
            "reference_matrix_present",
            "already_present",
            "STRICT_PRESENT",
        ),
        (
            "group_012",
            "reference_absent_partial",
            "reference_absent_clean",
            "absent_in_reference",
            "ABSENT",
        ),
    ];
    for (cluster, class, status, action, state) in &ref_specs {
        ref_audit.push_str(&format!(
            "{}\tL4_AS00146\t{}\t{}\treference_matrix\t{}\t{}\thigh\tpublication_grade\tConsistent with corrected matrix call.\tmatrix_present\t\t\t\t\t\t\t\t\t\t\t\t\n",
            cluster, class, status, action, state,
        ));
    }
    std::fs::write(out.path().join("reference_genome_audit.tsv"), ref_audit).unwrap();

    // reference cluster events: one cluster identified as event-derived ref allele.
    std::fs::write(
        out.path().join("reference_cluster_events.tsv"),
        "gene_cluster\ttarget_reference_cluster\ttarget_reference_gene\treference_event_class\tidentity\tquery_coverage\treference_coverage\tquery_length\treference_length\tevent_size_nt\tevent_size_aa\tmatches\tmismatches\tgaps_in_query\tgaps_in_reference\tshared_kmers\tmatrix_recommendation\tinterpretation\n\
group_009\tgroup_001\tRv0440\tpartial_reference_homology\t91.2\t98.5\t42.6\t237\t561\t-324\t-108\t218\t19\t324\t0\t1\tkeep_separate_cluster_annotate_event\tShort partial reference-gene homology; keep separate but flag as event-derived allele.\n",
    )
    .unwrap();

    // cluster_synteny_blocks: one block per cluster (or two for group_003 to show variation).
    // The schema columns:
    let synteny_header = "gene_cluster\tblock_rank\tn_genomes\tfraction_of_present\tupstream_block_5\tdownstream_block_5\tblock_interpretation\texample_genomes\tcarrier_genomes\trepresentative_genome\trepresentative_gene\trepresentative_contig\trepresentative_start\trepresentative_end\trepresentative_strand\trepresentative_length_bp\trepresentative_block_genes\treference_block_genes\n";
    let mut syn = String::from(synteny_header);
    // helper to build a representative gene block string of 5 upstream + target + 5 downstream
    let block_genes =
        |target: &str, target_genome: &str, start: u64, len: u64, variant: usize| -> String {
            // Each neighbour cluster carries its own Rv number so block signatures
            // differ between variants (signatures are computed from gene_id, so the
            // dedup would collapse rows otherwise). Variants share several neighbours
            // for inter-row connection lines, but swap a few at the edges/middle to
            // show divergence.
            let upstream: Vec<(&str, u32)> = match variant {
                0 => vec![
                    ("group_101", 1200),
                    ("group_102", 1201),
                    ("group_103", 1202),
                    ("group_104", 1203),
                    ("group_105", 1204),
                ],
                1 => vec![
                    ("group_201", 2200),
                    ("group_102", 1201),
                    ("group_103", 1202),
                    ("group_104", 1203),
                    ("group_105", 1204),
                ],
                2 => vec![
                    ("group_101", 1200),
                    ("group_102", 1201),
                    ("group_211", 2110),
                    ("group_104", 1203),
                    ("group_105", 1204),
                ],
                3 => vec![
                    ("group_201", 2200),
                    ("group_212", 2120),
                    ("group_103", 1202),
                    ("group_104", 1203),
                    ("group_105", 1204),
                ],
                _ => vec![
                    ("group_101", 1200),
                    ("group_102", 1201),
                    ("group_103", 1202),
                    ("group_104", 1203),
                    ("group_105", 1204),
                ],
            };
            let downstream: Vec<(&str, u32)> = match variant {
                0 => vec![
                    ("group_106", 1305),
                    ("group_107", 1306),
                    ("group_108", 1307),
                    ("group_109", 1308),
                    ("group_110", 1309),
                ],
                1 => vec![
                    ("group_106", 1305),
                    ("group_107", 1306),
                    ("group_108", 1307),
                    ("group_109", 1308),
                    ("group_220", 2200),
                ],
                2 => vec![
                    ("group_106", 1305),
                    ("group_107", 1306),
                    ("group_108", 1307),
                    ("group_220", 2200),
                    ("group_110", 1309),
                ],
                3 => vec![
                    ("group_106", 1305),
                    ("group_221", 2210),
                    ("group_222", 2220),
                    ("group_109", 1308),
                    ("group_110", 1309),
                ],
                _ => vec![
                    ("group_106", 1305),
                    ("group_107", 1306),
                    ("group_108", 1307),
                    ("group_109", 1308),
                    ("group_110", 1309),
                ],
            };
            let mut parts: Vec<String> = Vec::new();
            let pos = start.saturating_sub(6 * 900);
            for (i, (name, rv)) in upstream.iter().enumerate() {
                let g_start = pos + (i as u64) * 900;
                let g_end = g_start + 720;
                parts.push(format!(
                    "{},Rv{:04},+,{},{},{},{}",
                    name, rv, target_genome, g_start, g_end, 720
                ));
            }
            parts.push(format!(
                "{},Rv{:04},+,{},{},{},{}",
                target,
                1300,
                target_genome,
                start,
                start + len,
                len
            ));
            for (i, (name, rv)) in downstream.iter().enumerate() {
                let g_start = start + len + 200 + (i as u64) * 900;
                let g_end = g_start + 720;
                parts.push(format!(
                    "{},Rv{:04},+,{},{},{},{}",
                    name, rv, target_genome, g_start, g_end, 720
                ));
            }
            parts.join(";")
        };
    let reference_block_genes = |target: &str| -> String {
        let upstream_names = [
            "group_101",
            "group_102",
            "group_103",
            "group_104",
            "group_105",
        ];
        let downstream_names = [
            "group_106",
            "group_107",
            "group_108",
            "group_109",
            "group_110",
        ];
        let mut parts: Vec<String> = Vec::new();
        let pos = 1_500_000u64;
        for (i, name) in upstream_names.iter().enumerate() {
            let g_start = pos + (i as u64) * 900;
            let g_end = g_start + 720;
            parts.push(format!(
                "{},Rv{:04},+,AS00146,{},{},{}",
                name,
                1200 + i,
                g_start,
                g_end,
                720
            ));
        }
        parts.push(format!(
            "{},Rv{:04},+,AS00146,{},{},{}",
            target, 1300, 1_504_500, 1_505_700, 1200
        ));
        for (i, name) in downstream_names.iter().enumerate() {
            let g_start = 1_506_000 + (i as u64) * 900;
            let g_end = g_start + 720;
            parts.push(format!(
                "{},Rv{:04},+,AS00146,{},{},{}",
                name,
                1310 + i,
                g_start,
                g_end,
                720
            ));
        }
        parts.join(";")
    };
    #[rustfmt::skip]
    // Extended block table: several clusters get 3-4 alternative synteny blocks
    // so the connection lines between rows are clearly visible and the panel is
    // worth interacting with.
    let syn_blocks: Vec<SampleSyntenyBlock> = vec![
        ("group_001", "L4_AS00146", "AS00146", 1_504_500, 1_500, 16, 1.000, "dominant_block",                          0),
        ("group_002", "L4_AS00146", "AS00146",   801_200,   980,  3, 0.600, "dominant_block_minor_alternatives",       0),
        ("group_002", "L2_AS00021", "AS00021",   812_440,   980,  2, 0.400, "alternative_synteny_block",               1),
        ("group_003", "L4_AS00146", "AS00146", 1_734_734, 1_200,  8, 0.571, "dominant_block_minor_alternatives",       0),
        ("group_003", "A2_AS00357", "AS00357", 1_690_330, 1_200,  3, 0.214, "alternative_synteny_block",               1),
        ("group_003", "L1_AS00033", "AS00033", 1_812_990, 1_180,  2, 0.143, "alternative_synteny_block",               2),
        ("group_003", "A1_AS00349", "AS00349", 1_940_120, 1_180,  1, 0.071, "alternative_synteny_block",               3),
        ("group_004", "L4_AS00146", "AS00146", 2_440_010,   870, 11, 0.786, "dominant_block_minor_alternatives",       0),
        ("group_004", "L2_AS00021", "AS00021", 2_502_440,   870,  2, 0.143, "alternative_synteny_block",               1),
        ("group_004", "A4_AS00245", "AS00245", 2_690_320,   870,  1, 0.071, "alternative_synteny_block",               2),
        ("group_005", "L4_AS00146", "AS00146", 3_120_400, 1_410,  5, 0.714, "dominant_block_minor_alternatives",       0),
        ("group_005", "L1_AS00033", "AS00033", 3_240_120, 1_410,  2, 0.286, "alternative_synteny_block",               2),
        ("group_006", "L2_AS00021", "AS00021", 4_200_990,   980,  9, 0.692, "dominant_block_reference_absent",         0),
        ("group_006", "A4_AS00245", "AS00245", 4_310_140,   980,  3, 0.231, "alternative_synteny_block",               1),
        ("group_006", "L1_AS00033", "AS00033", 4_420_080,   980,  1, 0.077, "alternative_synteny_block",               3),
        ("group_007", "L4_AS00146", "AS00146",   980_100,   375, 10, 0.625, "dominant_block_minor_alternatives",       0),
        ("group_007", "L1_AS00033", "AS00033", 1_120_330,   375,  4, 0.250, "alternative_synteny_block_paralog_locus", 2),
        ("group_007", "A4_AS00245", "AS00245", 1_240_500,   375,  2, 0.125, "alternative_synteny_block",               1),
        ("group_008", "L4_AS00146", "AS00146", 1_120_000, 2_410,  6, 0.667, "dominant_block_minor_alternatives",       0),
        ("group_008", "L1_AS00033", "AS00033", 1_180_440, 2_410,  3, 0.333, "alternative_synteny_block",               1),
        ("group_009", "A4_AS00245", "AS00245",   544_900,   237,  3, 0.750, "dominant_block_partial_only",             0),
        ("group_009", "A1_AS00349", "AS00349",   612_080,   237,  1, 0.250, "alternative_synteny_block",               1),
        ("group_010", "L4_AS00146", "AS00146",   220_010,   612,  2, 0.500, "dominant_block_minor_alternatives",       0),
        ("group_010", "L3_AS00094", "AS00094",   320_140,   612,  1, 0.250, "alternative_synteny_block",               2),
        ("group_010", "L1_AS00033", "AS00033",   430_220,   612,  1, 0.250, "alternative_synteny_block",               3),
        ("group_011", "L4_AS00146", "AS00146",   745_500, 2_220, 16, 1.000, "dominant_block",                          0),
        ("group_012", "L7_AS00044", "AS00044", 1_990_100,   615,  3, 0.600, "dominant_block_minor_alternatives",       0),
        ("group_012", "L4_AS00146", "AS00146", 2_120_220,   615,  1, 0.200, "alternative_synteny_block",               1),
        ("group_012", "A4_AS00248", "AS00248", 2_240_440,   615,  1, 0.200, "alternative_synteny_block",               3),
    ];
    let mut rank_counter: HashMap<&str, usize> = HashMap::new();
    for (cluster, rep_genome, rep_contig, start, len, n, frac, interp, variant) in &syn_blocks {
        let rank_entry = rank_counter.entry(cluster).or_insert(0);
        *rank_entry += 1;
        let rank = *rank_entry;
        let carriers_list = cluster_specs
            .iter()
            .find(|(c, _, _)| c == cluster)
            .map(|(_, _, names)| names.to_vec())
            .unwrap_or_default()
            .join(";");
        let example = format!("{}:{}_{}", rep_genome, cluster, rep_genome);
        syn.push_str(&format!(
            "{cluster}\t{rank}\t{n}\t{frac:.3}\tgroup_105|group_104|group_103|group_102|group_101\tgroup_106|group_107|group_108|group_109|group_110\t{interp}\t{example}\t{carriers_list}\t{rep_genome}\t{cluster}_{rep_genome}\t{rep_contig}\t{start}\t{end}\t+\t{len}\t{rep_block}\t{ref_block}\n",
            end = start + len,
            rep_block = block_genes(cluster, rep_contig, *start, *len, *variant),
            ref_block = reference_block_genes(cluster),
        ));
    }
    std::fs::write(out.path().join("cluster_synteny_blocks.tsv"), syn).unwrap();

    // Phylogeny file path is also embedded in args (the file doesn't need to exist
    // for the report renderer because we pass the tree object directly, but the
    // banner reports the tree as "interactive" when args.tree is Some).
    let tree_path = out.path().join("tree.nwk");
    std::fs::write(&tree_path, b"(L4_AS00146,(L4_AS00150,L4_AS00200));\n").unwrap();

    let mut args = default_args_for_validation();
    args.output = out.path().to_path_buf();
    args.write_html_report = true;
    args.tree = Some(tree_path);
    args.reference_genome = Some("L4_AS00146".into());
    args.synteny = true;
    args.nscluster = true;
    args.matrix_policy = MatrixPolicy::Balanced;

    let stats = SummaryStats {
        // Must come from the matrix itself. A hardcoded value made the sample
        // report claim the core shrank after correction, which cannot happen:
        // Recore only ever adds presence.
        core_before: pa.core_size(),
        total_gaps: 64,
        filtered_gaps: 58,
        skipped_by_prevalence: 6,
        searchable_gaps: 58,
        missing_fasta_gaps: 0,
        dropped: 2,
        missing_cds: 1,
        split: 1,
        fragmented: 1,
        alignment_split: 0,
        start_variant: 0,
        paralog: 0,
        pseudogene: 1,
        synteny_rejected: 3,
        paralog_rejected: 1,
        matrix_policy: MatrixPolicy::Balanced,
        matrix_corrected: rescues.len(),
        nsc_matrix_corrected: 4,
        synteny_enabled: true,
        nsc_clusters: 1,
        nsc_merged_genes: 1,
        skipped_no_reference: 0,
        remaining_gaps: 0,
    };

    write_html_report(
        &args,
        &pa,
        &stats,
        &rescue_map,
        &Some(nsc),
        Some(&tree),
        &[],
        out.path(),
    )
    .unwrap();

    if let Some(dump) = std::env::var_os("RECORE_SAMPLE_DUMP") {
        let dest = std::path::PathBuf::from(dump);
        if let Some(parent) = dest.parent() {
            std::fs::create_dir_all(parent).ok();
        }
        std::fs::copy(out.path().join("recore_report.html"), &dest).ok();
    }

    let text = std::fs::read_to_string(out.path().join("recore_report.html")).unwrap();
    // Spot-check key sections.
    assert!(text.contains("Corrected matrix"));
    assert!(text.contains("L4_AS00146"));
    assert!(text.contains("group_001"));
    assert!(text.contains("ancestral_retained_no_disruption"));
    assert!(text.contains("clean_insertion_or_acquisition_candidate"));
    assert!(text.contains("<span class='mech-badge'"));
    assert!(text.contains("id='candidatesTable'"));
    assert!(text.contains("id='syntenyPanel'"));
    // Sortable headers must stay reachable and announced for keyboard/screen readers.
    assert!(text.contains("aria-sort='none'"));
    assert!(text.contains("scope='col'"));
    assert_no_remote_assets(&text);
}

#[test]
fn clear_optional_outputs_removes_stale_optional_files() {
    let out = tempfile::tempdir().unwrap();
    let stale_file = out.path().join("synteny_rejected_report.tsv");
    std::fs::write(&stale_file, b"old\n").unwrap();
    let stale_html = out.path().join("recore_report.html");
    std::fs::write(&stale_html, b"old\n").unwrap();
    let corrected_dir = out.path().join("corrected_gff");
    std::fs::create_dir_all(&corrected_dir).unwrap();
    std::fs::write(corrected_dir.join("old.recore.gff3"), b"old\n").unwrap();

    clear_optional_outputs(out.path()).unwrap();

    assert!(!stale_file.exists());
    assert!(!stale_html.exists());
    assert!(!corrected_dir.exists());
}

#[test]
fn corrected_gff_appends_features_before_fasta_section() {
    let out = tempfile::tempdir().unwrap();
    let source = out.path().join("g1.gff3");
    std::fs::write(
        &source,
        b"##gff-version 3\nchr1\tAMAP\tCDS\t1\t90\t.\t+\t0\tID=gene1\n##FASTA\n>chr1\nATGCATGC\n",
    )
    .unwrap();
    let pa =
        PresenceAbsence::new_for_test(vec!["clA".into()], vec!["g1".into()], vec![vec!["".into()]]);
    let rescue = RescueEvent {
        cluster_idx: ClusterIdx(0),
        genome_idx: GenomeIdx(0),
        contig: "chr1".into(),
        start: 99,
        end: 180,
        ..rescue_with_type(RescueType::Dropped, Some(true))
    };
    let mut gff_paths = HashMap::new();
    gff_paths.insert(GenomeIdx(0), source);

    write_corrected_gffs(&[rescue], &pa, &gff_paths, out.path()).unwrap();

    let text = std::fs::read_to_string(out.path().join("corrected_gff/g1.recore.gff3")).unwrap();
    let rescue_pos = text.find("\tRecore\tCDS\t").unwrap();
    let fasta_pos = text.find("##FASTA").unwrap();
    assert!(rescue_pos < fasta_pos);
    assert!(text.contains(">chr1\nATGCATGC"));
}

#[test]
fn gff_copy_helper_splits_fasta_tail() {
    let mut input = tempfile::NamedTempFile::new().unwrap();
    writeln!(input, "##gff-version 3").unwrap();
    writeln!(
        input,
        "chr1\tAMAP\tCDS\t1\t90\t.\t+\t0\tID=gene1;Name=gene1"
    )
    .unwrap();
    writeln!(input, "##FASTA").unwrap();
    writeln!(input, ">chr1").unwrap();
    writeln!(input, "ATGC").unwrap();

    let mut output = Vec::new();
    let tail = copy_gff_records_before_fasta_with_header(input.path(), &mut output)
        .unwrap()
        .unwrap();
    let text = String::from_utf8(output).unwrap();
    let tail_text = String::from_utf8(tail).unwrap();

    assert!(text.contains("ID=gene1"));
    assert!(!text.contains("##FASTA"));
    assert!(tail_text.starts_with("##FASTA"));
}

/// The corrected GFF3 writer reads the user's source annotation itself, so it
/// truncated a bgzipped one exactly the way the parser did, and separately.
/// Every other test of this helper hands it a `NamedTempFile`, whose generated
/// name has no extension at all, so the compressed branch was never once
/// executed. Without this the production half of the fix can be reverted and
/// the suite stays green, which is the blindness that let the defect live.
#[test]
fn gff_copy_helper_reads_every_member_of_a_multi_member_gz() {
    use flate2::write::GzEncoder;
    use flate2::Compression;

    let dir = tempfile::tempdir().unwrap();
    let gz_path = dir.path().join("annotation.gff3.gz");
    let members = [
        "##gff-version 3\nchr1\tAMAP\tCDS\t1\t90\t.\t+\t0\tID=gene1;Name=gene1\n",
        "chr1\tAMAP\tCDS\t200\t290\t.\t+\t0\tID=gene2;Name=gene2\n",
    ];
    {
        let mut file = std::fs::File::create(&gz_path).unwrap();
        for member in members {
            let mut gz = GzEncoder::new(Vec::new(), Compression::fast());
            gz.write_all(member.as_bytes()).unwrap();
            file.write_all(&gz.finish().unwrap()).unwrap();
        }
    }

    let mut output = Vec::new();
    copy_gff_records_before_fasta_with_header(&gz_path, &mut output).unwrap();
    let text = String::from_utf8(output).unwrap();

    assert!(text.contains("ID=gene1"), "first member missing: {text}");
    assert!(
        text.contains("ID=gene2"),
        "second member dropped, the file was truncated: {text}"
    );
}

#[test]
fn reference_audit_quality_guard_keeps_weak_sequence_only_ambiguous() {
    let mut rescue = rescue_with_type(RescueType::Dropped, None);
    rescue.identity = 94.0;
    rescue.coverage = 86.0;
    rescue.end = 120;
    let state = "SEQUENCE_ONLY";
    let quality = orthology_quality(&rescue, state, true, 0.02, 120);

    assert_eq!(
        matrix_action_with_quality(MatrixPolicy::Balanced, true, state, &quality),
        "audit_only_quality"
    );
    assert_eq!(
        quality_aware_reference_publication_status(state, true, &quality),
        "flag_reference_like_ambiguous"
    );
    assert_eq!(
        quality_aware_reference_publication_confidence(
            state,
            reference_status_for(state, true),
            true,
            &quality,
        ),
        "unresolved"
    );
    assert!(quality_aware_reference_wording(state, true, &quality)
        .contains("below adaptive quality thresholds"));
}

#[test]
fn paralog_requires_explicit_synteny_support() {
    let no_synteny = rescue_with_type(
        RescueType::Paralog {
            other_cluster: "other".into(),
        },
        None,
    );
    assert!(is_paralog_without_synteny_support(&no_synteny));

    let failed_synteny = rescue_with_type(
        RescueType::Paralog {
            other_cluster: "other".into(),
        },
        Some(false),
    );
    assert!(is_paralog_without_synteny_support(&failed_synteny));

    let passed_synteny = rescue_with_type(
        RescueType::Paralog {
            other_cluster: "other".into(),
        },
        Some(true),
    );
    assert!(!is_paralog_without_synteny_support(&passed_synteny));

    let dropped = rescue_with_type(RescueType::Dropped, None);
    assert!(!is_paralog_without_synteny_support(&dropped));
}

#[test]
fn best_rescue_map_keeps_strongest_duplicate_event() {
    let weak = rescue_with_type(RescueType::Dropped, None);
    let mut strong = rescue_with_type(RescueType::MissingCds, Some(true));
    strong.identity = 99.0;
    strong.coverage = 100.0;
    strong.confidence = 0.99;
    strong.start = 200;
    strong.end = 500;

    let rescues = [strong.clone(), weak.clone()];
    let selected = build_best_rescue_map(rescues.iter());
    let rescue = selected.get(&(ClusterIdx(0), GenomeIdx(0))).unwrap();

    assert!(matches!(rescue.rescue_type, RescueType::MissingCds));
    assert_eq!(rescue.start, 200);
}

#[test]
fn best_rescue_map_uses_coordinate_tie_breaker() {
    let mut later = rescue_with_type(RescueType::Dropped, None);
    later.contig = "chr2".into();
    later.start = 200;
    later.end = 300;
    let mut earlier = rescue_with_type(RescueType::Dropped, None);
    earlier.contig = "chr1".into();
    earlier.start = 10;
    earlier.end = 110;

    let rescues = [later, earlier];
    let selected = build_best_rescue_map(rescues.iter());
    let rescue = selected.get(&(ClusterIdx(0), GenomeIdx(0))).unwrap();

    assert_eq!(rescue.contig, "chr1");
    assert_eq!(rescue.start, 10);
}

#[test]
fn tsv_fields_are_single_line_safe() {
    assert_eq!(clean_tsv_field("a\tb\nc\rd\0e"), "a b c d e");
}

#[test]
fn sensitive_policy_does_not_correct_paralog_state() {
    assert!(!matrix_policy_corrects(
        MatrixPolicy::Sensitive,
        "HOMOLOG_PARALOG"
    ));
}

#[test]
fn reference_genome_resolution_accepts_embedded_as_id() {
    let pa = crate::parser::PresenceAbsence::new_for_test(
        vec!["group_1".into()],
        vec!["L4_AS00146".into(), "AS00999".into()],
        vec![vec!["gene1".into(), "".into()]],
    );

    assert_eq!(resolve_reference_genome(&pa, "AS00146"), Some(GenomeIdx(0)));
    assert_eq!(
        resolve_reference_genome(&pa, "L4_AS00146"),
        Some(GenomeIdx(0))
    );
    assert_eq!(resolve_reference_genome(&pa, "AS00000"), None);
}

#[test]
fn evidence_states_separate_strict_alias_and_paralog() {
    let strict = rescue_with_type(RescueType::Dropped, Some(true));
    assert_eq!(
        evidence_state(&strict, true, "group_1", "", ""),
        "STRICT_PRESENT"
    );

    let alias = rescue_with_type(
        RescueType::Split {
            other_cluster: "group_2".into(),
        },
        Some(true),
    );
    assert_eq!(
        evidence_state(&alias, true, "group_1", "group_2", ""),
        "ALIAS_CLUSTER"
    );

    let paralog = rescue_with_type(
        RescueType::Paralog {
            other_cluster: "group_2".into(),
        },
        None,
    );
    assert_eq!(
        evidence_state(
            &paralog,
            false,
            "group_1",
            "group_2",
            "paralog_without_synteny"
        ),
        "HOMOLOG_PARALOG"
    );
}

#[test]
fn reference_shadow_candidate_flags_partial_paralog_evidence() {
    let partial_hit = RescueEvent {
        coverage: 91.0,
        identity: 100.0,
        annotation_match: Some("rv_duf222".into()),
        ..rescue_with_type(
            RescueType::Split {
                other_cluster: "group_537".into(),
            },
            None,
        )
    };

    let state = reference_shadow_state(&partial_hit, "group_1044", "group_537");
    assert_eq!(state, "PARTIAL_PARALOG_EVIDENCE");
    assert_eq!(
        reference_publication_status(state, false),
        "flag_reference_partial_paralog_evidence"
    );
    assert_eq!(
        reference_publication_confidence(state, reference_status_for(state, false)),
        "complex_paralog"
    );
}

#[test]
fn reference_hit_shadow_candidates_capture_partial_annotated_paralog() {
    let pa = crate::parser::PresenceAbsence::new_for_test(
        vec!["group_1044".into(), "group_537".into()],
        vec!["L4_AS00146".into()],
        vec![vec!["".into()], vec!["rv_duf222".into()]],
    );
    let gff = GffIndex::build(vec![GffGene {
        contig: "chr".into(),
        start: 0,
        end: 954,
        strand: '+',
        gene_id: "rv_duf222".into(),
    }]);
    let mut gff_indices = HashMap::new();
    gff_indices.insert(GenomeIdx(0), gff);
    let mut gene_map = HashMap::new();
    gene_map.insert("rv_duf222".to_string(), ClusterIdx(1));
    let mut gene_to_cluster_maps = HashMap::new();
    gene_to_cluster_maps.insert(GenomeIdx(0), gene_map);
    let mut hits = HashMap::new();
    hits.insert(
        (GenomeIdx(0), ClusterIdx(0)),
        vec![Hit {
            query: "group_1044".into(),
            cluster: "group_1044".into(),
            query_length: 972,
            query_start: 0,
            query_end: 885,
            strand: '+',
            target: "chr".into(),
            target_start: 0,
            target_end: 885,
            target_length: 0,
            matches: 885,
            alignment_length: 885,
            mapq: 0,
        }],
    );

    let candidates =
        reference_hit_shadow_candidates(&hits, &pa, &gff_indices, &gene_to_cluster_maps);

    assert_eq!(candidates.len(), 1);
    assert!(matches!(
        candidates[0].rescue_type,
        RescueType::Paralog { .. }
    ));
    assert_eq!(
        reference_shadow_state(&candidates[0], "group_1044", "group_537"),
        "PARTIAL_PARALOG_EVIDENCE"
    );
}

#[test]
fn reference_cluster_events_detect_short_in_frame_allele() {
    let left = b"ATGAAACCCGGGTTTAAACCCGGGTTTAAACCC";
    let middle = b"GCTGCTGCTGCTGCTGCTGCTGCTGCTGCT";
    let right = b"GGGAAATTTCCCGGGAAATTTCCCGGGAAATAA";
    let mut reference_seq = Vec::new();
    reference_seq.extend_from_slice(left);
    reference_seq.extend_from_slice(middle);
    reference_seq.extend_from_slice(right);
    let mut short_seq = Vec::new();
    short_seq.extend_from_slice(left);
    short_seq.extend_from_slice(right);

    assert_eq!((reference_seq.len() - short_seq.len()) % 3, 0);

    let pa = crate::parser::PresenceAbsence::new_for_test(
        vec!["group_short".into(), "group_ref".into()],
        vec!["L4_AS00146".into()],
        vec![vec!["".into()], vec!["Rv0601c".into()]],
    );
    let mut pan_ref = HashMap::new();
    pan_ref.insert("group_short".to_string(), vec![short_seq]);
    pan_ref.insert("group_ref".to_string(), vec![reference_seq]);

    let events = build_reference_cluster_events(&pa, &pan_ref, GenomeIdx(0));

    assert_eq!(events.len(), 1);
    let event = &events[0];
    assert_eq!(event.query_cluster_idx, ClusterIdx(0));
    assert_eq!(event.target_cluster_idx, ClusterIdx(1));
    assert_eq!(event.target_gene_id, "Rv0601c");
    assert_eq!(event.event_class, "short_in_frame_allele");
    assert_eq!(event.event_size_nt, -30);
    assert_eq!(event.event_size_aa, Some(-10));
    assert!(event.identity >= 99.0);
}

#[test]
fn reference_gene_display_prefers_h37rv_locus() {
    assert_eq!(display_reference_gene_id("L4_AS00146_Rv0601c"), "Rv0601c");
    assert_eq!(
        display_reference_gene_id("L4_AS00146_Rv0601c_pseudo"),
        "Rv0601c"
    );
    assert_eq!(
        display_reference_gene_id("domain_dup_family"),
        "domain_dup_family"
    );
    assert_eq!(display_reference_gene_id("custom_gene"), "custom_gene");
}

#[test]
fn newick_fitch_detects_clade_gain() {
    let pa = crate::parser::PresenceAbsence::new_for_test(
        vec!["g".into()],
        vec!["s1".into(), "s2".into(), "s3".into(), "s4".into()],
        vec![vec!["".into(), "g_s2".into(), "g_s3".into(), "".into()]],
    );
    let mut root = NewickParser::new("(s1,(s2:0.1,s3:0.1),s4);")
        .parse()
        .unwrap();
    let (matched, unmatched) = resolve_phylo_tips(&mut root, &pa);
    assert_eq!(matched, 4);
    assert_eq!(unmatched, 0);
    let tree = PhyloTree { root };
    let present = vec![false, true, true, false];

    let summary = summarize_phylogenetic_pattern(Some(&tree), &present);

    assert_eq!(summary.present_tips, 2);
    assert_eq!(summary.gains, 1);
    assert_eq!(summary.losses, 0);
    assert!(summary.carrier_monophyly);
    assert_eq!(summary.pattern, "lineage_specific_gain");
}

#[test]
fn newick_parser_preserves_branch_lengths() {
    let root = NewickParser::new("(s1:0.25,(s2:1e-3,s3:0.4):0.5);")
        .parse()
        .unwrap();

    assert_eq!(root.branch_length, None);
    assert!((root.children[0].branch_length.unwrap() - 0.25).abs() < 1e-12);
    assert!((root.children[1].branch_length.unwrap() - 0.5).abs() < 1e-12);
    assert!((root.children[1].children[0].branch_length.unwrap() - 0.001).abs() < 1e-12);
}

#[test]
fn ambiguous_root_is_not_called_ancestral_loss() {
    let pa = crate::parser::PresenceAbsence::new_for_test(
        vec!["g".into()],
        vec!["s1".into(), "s2".into(), "s3".into(), "s4".into()],
        vec![vec!["g_s1".into(), "".into(), "g_s3".into(), "".into()]],
    );
    let mut root = NewickParser::new("((s1:0.1,s2:0.1),(s3:0.1,s4:0.1));")
        .parse()
        .unwrap();
    let (matched, unmatched) = resolve_phylo_tips(&mut root, &pa);
    assert_eq!((matched, unmatched), (4, 0));
    let tree = PhyloTree { root };
    let present = vec![true, false, true, false];

    let summary = summarize_phylogenetic_pattern(Some(&tree), &present);
    let penalized = penalized_phylo_summary(&summary);

    assert_eq!(summary.root_state, "ambiguous");
    assert!(summary.root_state_ambiguous);
    assert_ne!(summary.pattern, "ancestral_loss");
    assert_eq!(summary.pattern, "ambiguous_gain_loss");
    assert_eq!(summary.carrier_mrca_tips, 4);
    assert_eq!(summary.dollo_loss_events_in_mrca, 2);
    assert!(penalized.dollo_violation_count >= 1);
}

#[test]
fn carrier_mrca_purity_tracks_single_gain_with_internal_loss() {
    let pa = crate::parser::PresenceAbsence::new_for_test(
        vec!["g".into()],
        vec!["s1".into(), "s2".into(), "s3".into(), "s4".into()],
        vec![vec!["".into(), "g_s2".into(), "".into(), "g_s4".into()]],
    );
    let mut root = NewickParser::new("(s1:0.1,(s2:0.1,(s3:0.1,s4:0.1):0.1):0.1);")
        .parse()
        .unwrap();
    let (matched, unmatched) = resolve_phylo_tips(&mut root, &pa);
    assert_eq!((matched, unmatched), (4, 0));
    let tree = PhyloTree { root };
    let present = vec![false, true, false, true];

    let summary = summarize_phylogenetic_pattern(Some(&tree), &present);
    let penalized = penalized_phylo_summary(&summary);

    assert_eq!(summary.carrier_mrca_tips, 3);
    assert!((summary.carrier_mrca_purity - 2.0 / 3.0).abs() < 1e-9);
    assert_eq!(summary.dollo_loss_events_in_mrca, 1);
    assert_eq!(penalized.dollo_loss_count, 1);
    assert_eq!(penalized.dollo_origin_tips, 3);
}

#[test]
fn mobile_element_detection_uses_gff_aliases() {
    let pa = crate::parser::PresenceAbsence::new_for_test(
        vec!["g".into()],
        vec!["s1".into()],
        vec![vec!["".into()]],
    );
    let gff = GffIndex::build_with_aliases(vec![(
        GffGene {
            contig: "chr".into(),
            start: 150,
            end: 300,
            strand: '+',
            gene_id: "mobile1".into(),
        },
        vec!["IS6110 transposase".into()],
    )]);
    let mut gff_indices = HashMap::new();
    gff_indices.insert(GenomeIdx(0), gff);
    let rescue = RescueEvent {
        cluster_idx: ClusterIdx(0),
        genome_idx: GenomeIdx(0),
        rescue_type: RescueType::Fragmented {
            n_fragments: 2,
            total_coverage: 91.0,
        },
        identity: 99.0,
        coverage: 91.0,
        contig: "chr".into(),
        start: 100,
        end: 220,
        strand: '+',
        annotation_match: None,
        confidence: 0.9,
        synteny_ok: None,
        synteny_evidence: SyntenyEvidence::default(),
        fragment_evidence: crate::rescue::FragmentEvidence::default(),
    };

    let mobile = nearest_mobile_element_for_rescue(&rescue, &pa, &gff_indices, &[], 1_000).unwrap();

    assert_eq!(mobile.source, "annotation");
    assert_eq!(mobile.kind, "IS6110");
    assert_eq!(mobile.distance_bp, 0);
}

#[test]
fn mobile_element_sequence_hits_are_name_agnostic_and_prioritized() {
    let pa = crate::parser::PresenceAbsence::new_for_test(
        vec!["g".into()],
        vec!["s1".into()],
        vec![vec!["".into()]],
    );
    let mut gff_indices = HashMap::new();
    gff_indices.insert(
        GenomeIdx(0),
        GffIndex::build_with_aliases(vec![(
            GffGene {
                contig: "chr".into(),
                start: 90,
                end: 210,
                strand: '+',
                gene_id: "anonymous_cds".into(),
            },
            vec!["hypothetical protein".into()],
        )]),
    );
    let rescue = RescueEvent {
        cluster_idx: ClusterIdx(0),
        genome_idx: GenomeIdx(0),
        rescue_type: RescueType::Fragmented {
            n_fragments: 2,
            total_coverage: 91.0,
        },
        identity: 99.0,
        coverage: 91.0,
        contig: "chr".into(),
        start: 100,
        end: 220,
        strand: '+',
        annotation_match: None,
        confidence: 0.9,
        synteny_ok: None,
        synteny_evidence: SyntenyEvidence::default(),
        fragment_evidence: crate::rescue::FragmentEvidence::default(),
    };
    let mobile_hits = vec![MobileElementHit {
        genome_idx: GenomeIdx(0),
        query_name: "ISfinder_anon_001".into(),
        contig: "chr".into(),
        start: 120,
        end: 260,
        strand: '+',
        identity: 97.5,
        coverage: 82.0,
    }];

    let mobile =
        nearest_mobile_element_for_rescue(&rescue, &pa, &gff_indices, &mobile_hits, 1_000).unwrap();

    assert_eq!(mobile.source, "sequence");
    assert_eq!(mobile.kind, "mobile_element_sequence");
    assert_eq!(mobile.distance_bp, 0);
    assert!(mobile.label.contains("ISfinder_anon_001"));
}

#[test]
fn nsc_role_labels_member_and_representative() {
    let mut merge_map = HashMap::new();
    merge_map.insert(ClusterIdx(2), ClusterIdx(1));
    let nsc = Some(nscluster::NsClusterResult {
        merged_clusters: vec![nscluster::NsCluster {
            members: vec![ClusterIdx(1), ClusterIdx(2)],
            representative: "group_1".into(),
        }],
        total_merged_genes: 2,
        n_clusters: 1,
        merge_map,
    });

    assert_eq!(
        nsc_role_for_cluster(ClusterIdx(2), &nsc),
        "nsc_member_merged_into_representative"
    );
    assert_eq!(
        nsc_role_for_cluster(ClusterIdx(1), &nsc),
        "nsc_representative"
    );
    assert_eq!(nsc_role_for_cluster(ClusterIdx(3), &nsc), "");
}

#[test]
fn nsc_presence_union_includes_matrix_rescued_members() {
    let pa = crate::parser::PresenceAbsence::new_for_test(
        vec!["rep".into(), "member".into()],
        vec!["g1".into(), "g2".into()],
        vec![vec!["rep_g1".into(), "".into()], vec!["".into(), "".into()]],
    );
    let mut merge_map = HashMap::new();
    merge_map.insert(ClusterIdx(1), ClusterIdx(0));
    let nsc = Some(nscluster::NsClusterResult {
        merged_clusters: vec![nscluster::NsCluster {
            members: vec![ClusterIdx(0), ClusterIdx(1)],
            representative: "rep".into(),
        }],
        total_merged_genes: 2,
        n_clusters: 1,
        merge_map,
    });
    let member_rescue = RescueEvent {
        cluster_idx: ClusterIdx(1),
        genome_idx: GenomeIdx(1),
        ..rescue_with_type(RescueType::Dropped, Some(true))
    };
    let mut rescue_map = HashMap::new();
    rescue_map.insert(
        (member_rescue.cluster_idx, member_rescue.genome_idx),
        &member_rescue,
    );

    let updates = nsc_presence_updates(&pa, &nsc, &rescue_map);

    assert!(updates.contains(&(0, 1)));
}

/// The summary's matrix-corrected count, its remaining gaps and the report's
/// share come from what write_corrected_csv writes. A rescue in a row that NSC
/// merges away is not written, and used to be counted as a correction anyway;
/// "Remaining gaps" used to be searched gaps minus corrections, which is not
/// what is left in the file.
#[test]
fn corrected_matrix_counts_only_what_is_written() {
    let dir = tempfile::tempdir().unwrap();
    let csv_path = dir.path().join("pa.csv");
    std::fs::write(
        &csv_path,
        "Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,\
         Avg sequences per isolate,Genome Fragment,Order within Fragment,\
         Accessory Order,QC,Min group size nuc,Max group size nuc,Avg group size nuc,g1,g2,g3\n\
         rep,,p,1,1,1,0,0,0,,1,1,1,rep_g1,,\n\
         member,,p,0,0,0,0,0,0,,1,1,1,,,\n\
         other,,p,1,1,1,0,0,0,,1,1,1,other_g1,,\n",
    )
    .unwrap();
    let pa = crate::parser::parse_presence_absence(&csv_path).unwrap();

    let mut merge_map = HashMap::new();
    merge_map.insert(ClusterIdx(1), ClusterIdx(0));
    let nsc = Some(nscluster::NsClusterResult {
        merged_clusters: vec![nscluster::NsCluster {
            members: vec![ClusterIdx(0), ClusterIdx(1)],
            representative: "rep".into(),
        }],
        total_merged_genes: 2,
        n_clusters: 1,
        merge_map,
    });
    let member_rescue = RescueEvent {
        cluster_idx: ClusterIdx(1),
        genome_idx: GenomeIdx(1),
        ..rescue_with_type(RescueType::Dropped, Some(true))
    };
    let other_rescue = RescueEvent {
        cluster_idx: ClusterIdx(2),
        genome_idx: GenomeIdx(1),
        ..rescue_with_type(RescueType::Dropped, Some(true))
    };
    let mut rescue_map = HashMap::new();
    rescue_map.insert((ClusterIdx(1), GenomeIdx(1)), &member_rescue);
    rescue_map.insert((ClusterIdx(2), GenomeIdx(1)), &other_rescue);

    let counts = write_corrected_csv(&pa, &rescue_map, dir.path(), &nsc).unwrap();
    assert_eq!(
        counts.rescued_cells, 1,
        "the member's rescue sits in a row that is not written"
    );
    assert_eq!(
        counts.nsc_cells, 1,
        "it reaches rep/g2 as NSC_MERGED instead"
    );

    // rep: rep_g1, NSC_MERGED, "" and other: other_g1, RESCUED, "".
    let written =
        std::fs::read_to_string(dir.path().join("gene_presence_absence_recore.csv")).unwrap();
    let empty: usize = written
        .lines()
        .skip(1)
        .map(|l| l.split(',').skip(13).filter(|c| c.is_empty()).count())
        .sum();
    assert_eq!(
        written.lines().count(),
        3,
        "header plus two rows: member is merged away"
    );
    assert_eq!(counts.absent_cells, empty);
    assert_eq!(counts.absent_cells, 2);
}

#[test]
fn reference_shadow_candidate_prefers_synteny_supported_evidence() {
    let weak = RescueEvent {
        cluster_idx: ClusterIdx(1),
        genome_idx: GenomeIdx(2),
        confidence: 1.5,
        synteny_ok: Some(false),
        ..rescue_with_type(RescueType::Dropped, Some(false))
    };
    let supported = RescueEvent {
        cluster_idx: ClusterIdx(1),
        genome_idx: GenomeIdx(2),
        confidence: 0.9,
        synteny_ok: Some(true),
        ..rescue_with_type(RescueType::Dropped, Some(true))
    };
    let candidates = vec![weak, supported];

    let best = best_reference_candidate(&candidates, ClusterIdx(1), GenomeIdx(2)).unwrap();
    assert_eq!(best.synteny_ok, Some(true));
}

#[test]
fn matrix_policy_controls_binary_correction() {
    assert!(matrix_policy_corrects(
        MatrixPolicy::Strict,
        "STRICT_PRESENT"
    ));
    assert!(matrix_policy_corrects(
        MatrixPolicy::Strict,
        "ALIAS_CLUSTER"
    ));
    assert!(!matrix_policy_corrects(
        MatrixPolicy::Strict,
        "SEQUENCE_ONLY"
    ));

    assert!(matrix_policy_corrects(
        MatrixPolicy::Balanced,
        "SEQUENCE_ONLY"
    ));
    assert!(!matrix_policy_corrects(
        MatrixPolicy::Balanced,
        "PSEUDOGENE"
    ));
    assert!(!matrix_policy_corrects(
        MatrixPolicy::Balanced,
        "PARTIAL_HOMOLOGY"
    ));
    assert!(!matrix_policy_corrects(
        MatrixPolicy::Balanced,
        "MISSING_CDS_UNVALIDATED"
    ));

    assert!(matrix_policy_corrects(
        MatrixPolicy::Sensitive,
        "PSEUDOGENE"
    ));
    assert!(matrix_policy_corrects(
        MatrixPolicy::Sensitive,
        "PARTIAL_HOMOLOGY"
    ));
    assert!(!matrix_policy_corrects(
        MatrixPolicy::Sensitive,
        "MISSING_CDS_UNVALIDATED"
    ));
}

#[test]
fn gff_stream_copy_adds_missing_header() {
    let mut input = tempfile::NamedTempFile::new().unwrap();
    writeln!(
        input,
        "chr1\tAMAP\tCDS\t1\t90\t.\t+\t0\tID=gene1;Name=gene1"
    )
    .unwrap();

    let mut output = Vec::new();
    copy_text_maybe_gzip_with_header(input.path(), &mut output).unwrap();
    let text = String::from_utf8(output).unwrap();

    assert!(text.starts_with("##gff-version 3\n"));
    assert!(text.contains("ID=gene1"));
}

#[test]
fn gff_stream_copy_preserves_existing_header_after_blank_prefix() {
    let mut input = tempfile::NamedTempFile::new().unwrap();
    writeln!(input).unwrap();
    writeln!(input, "##gff-version 3").unwrap();
    writeln!(
        input,
        "chr1\tAMAP\tCDS\t1\t90\t.\t+\t0\tID=gene1;Name=gene1"
    )
    .unwrap();

    let mut output = Vec::new();
    copy_text_maybe_gzip_with_header(input.path(), &mut output).unwrap();
    let text = String::from_utf8(output).unwrap();

    assert_eq!(text.matches("##gff-version 3").count(), 1);
    assert!(text.starts_with('\n'));
    assert!(text.contains("ID=gene1"));
}

/// A single-copy gene beside a pseudogene, a missing CDS or a true loss gets a
/// different context signature in those genomes, but the same neighbourhood.
/// Signature differences alone used to label it a duplication or paralog
/// expansion; only a neighbourhood with nothing in common is another locus.
#[test]
fn locus_context_is_judged_by_shared_neighbours() {
    let ctx = |up: &str, down: &str| (up.to_string(), down.to_string());
    let usual = ctx("Rv0006|Rv0005", "Rv0008c|Rv0009");
    assert!(evolution::shares_neighbour(
        &usual,
        &ctx("Rv0005|Rv0004", "Rv0008c")
    ));
    assert!(evolution::shares_neighbour(&usual, &ctx("NA", "NA")));
    assert!(!evolution::shares_neighbour(
        &usual,
        &ctx("Rv3001|Rv3002", "Rv3004|Rv3005")
    ));

    // Interpretation thresholds now read the compatible fraction.
    assert_eq!(
        evolution::context_interpretation(12, 4, 1.0),
        "stable_single_locus_minor_noise"
    );
    assert_eq!(
        evolution::context_interpretation(12, 2, 0.5),
        "possible_multilocus_or_paralog"
    );
}

/// Parsimony at a node with more than two children must not depend on the
/// order the children are written in. Pairwise Fitch did: one ordering gave
/// dispersed gain and loss, the other an ancestral loss.
#[test]
fn parsimony_is_independent_of_sibling_order() {
    let genomes: Vec<String> = ["a", "b", "c", "d", "e"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let pa = PresenceAbsence::new_for_test(
        vec!["gene".into()],
        genomes,
        vec![vec![
            "".into(),
            "".into(),
            "x".into(),
            "x".into(),
            "x".into(),
        ]],
    );
    let present: Vec<bool> = (0..5).map(|i| i >= 2).collect();
    let mut outcomes = Vec::new();
    for newick in ["(e,(a,b,c,d));", "(e,(c,a,d,b));", "(e,(d,c,b,a));"] {
        let mut root = evolution::NewickParser::new(newick).parse().unwrap();
        evolution::resolve_phylo_tips(&mut root, &pa);
        let summary = evolution::fitch_postorder(&root, &present);
        let (mut gains, mut losses) = (0, 0);
        evolution::assign_fitch_states(&summary, None, &mut gains, &mut losses);
        assert_eq!(summary.changes, gains + losses, "{newick}");
        outcomes.push((summary.changes, gains, losses));
    }
    assert!(outcomes.windows(2).all(|w| w[0] == w[1]), "{outcomes:?}");
    assert_eq!(outcomes[0].0, 2, "two losses, a and b");
}

/// "<!--" followed by "<script" in any name or annotation switched the HTML
/// parser into the script's escaped state and left the report blank. In JSON
/// these characters only occur inside strings, so escaping them keeps every
/// value as it was.
#[test]
fn report_payload_is_safe_inside_a_script_element() {
    let json = "{\"a\":\"<!--<script>\",\"b\":\"x & y > z\",\"c\":\"line\u{2028}sep\"}";
    let safe = output::script_safe_json(json);
    assert!(!safe.contains('<') && !safe.contains('>') && !safe.contains('&'));
    assert!(!safe.contains('\u{2028}'));
    assert!(safe.contains("\\u003c!--\\u003cscript\\u003e"));
}

fn evidence_row(cluster: usize, genome: usize, values: &[(&str, &str)]) -> OrthologyEvidenceRow {
    let fields = ORTHOLOGY_EVIDENCE_COLUMNS
        .iter()
        .map(|column| {
            values
                .iter()
                .find(|(c, _)| c == column)
                .map(|(_, v)| v.to_string())
                .unwrap_or_default()
        })
        .collect();
    OrthologyEvidenceRow {
        cluster_idx: ClusterIdx(cluster),
        genome_idx: GenomeIdx(genome),
        fields,
    }
}

fn test_base64_decode(text: &str) -> Vec<u8> {
    let alphabet = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut bits = 0u32;
    let mut n_bits = 0;
    let mut out = Vec::new();
    for c in text.bytes().filter(|&c| c != b'=') {
        bits = (bits << 6) | alphabet.iter().position(|&a| a == c).unwrap() as u32;
        n_bits += 6;
        if n_bits >= 8 {
            n_bits -= 8;
            out.push((bits >> n_bits) as u8);
            bits &= (1 << n_bits) - 1;
        }
    }
    out
}

/// The JSON string array that follows `"key":` in `json` (no escapes inside).
fn test_json_strings(json: &str, key: &str) -> Vec<String> {
    let start = json.find(&format!("\"{key}\":[")).unwrap() + key.len() + 4;
    let end = start + json[start..].find(']').unwrap();
    json[start..end]
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.trim_matches('"').to_string())
        .collect()
}

/// The explorer reads each cell's evidence from one base64 string per cluster.
/// Decoding it must give back the rows of orthology_evidence.tsv, field for
/// field, in genome order, with empty numbers as "not given".
#[test]
fn explorer_evidence_round_trips_per_cluster() {
    let rows = vec![
        evidence_row(
            0,
            2,
            &[
                ("decision", "accepted"),
                ("matrix_action", "audit_only_policy"),
                ("rescue_type", "PSEUDOGENE:1stops"),
                ("contig", "g3_ctg"),
                ("start", "101"),
                ("end", "400"),
                ("identity", "99.1"),
                ("coverage", "100.0"),
                ("strand", "-"),
            ],
        ),
        evidence_row(
            0,
            1,
            &[
                ("decision", "accepted"),
                ("matrix_action", "correct_presence"),
                ("rescue_type", "DROPPED"),
                ("contig", "g2_ctg"),
                ("start", "7"),
                ("end", "1206"),
                ("identity", "93.6"),
                ("coverage", ""),
                ("strand", "+"),
            ],
        ),
        // A cluster the explorer does not draw (merged away by NSC) is left out.
        evidence_row(5, 0, &[("decision", "accepted"), ("contig", "hidden")]),
    ];
    let rendered: HashMap<usize, usize> = [(0usize, 0usize), (1, 1)].into_iter().collect();
    let json = output::explorer_evidence_json(&rows, &rendered, 2);
    let strings = test_json_strings(&json, "strings");
    let by_cluster = test_json_strings(&json, "byCluster");
    assert_eq!(by_cluster.len(), 2);
    assert_eq!(
        by_cluster[1], "",
        "a cluster with no evidence costs nothing"
    );
    assert!(!strings.iter().any(|s| s == "hidden"));

    let bytes = test_base64_decode(&by_cluster[0]);
    let mut pos = 0;
    let mut next = || {
        let (mut value, mut shift) = (0u64, 0);
        loop {
            let byte = bytes[pos];
            pos += 1;
            value |= u64::from(byte & 0x7f) << shift;
            if byte & 0x80 == 0 {
                return value;
            }
            shift += 7;
        }
    };
    assert_eq!(next(), 2);
    let mut genome = 0;
    let mut decoded = Vec::new();
    for _ in 0..2 {
        genome += next();
        let _profile = next();
        let text: Vec<String> = (0..output::EXPLORER_EVIDENCE_TEXT.len())
            .map(|_| strings[next() as usize].clone())
            .collect();
        let numbers: Vec<u64> = (0..4).map(|_| next()).collect();
        decoded.push((genome, text, numbers));
    }
    // Genome order, not input order.
    assert_eq!(decoded[0].0, 1);
    assert_eq!(decoded[0].1[0], "DROPPED");
    assert_eq!(decoded[0].1[1], "g2_ctg");
    assert_eq!(decoded[0].2, vec![8, 1207, 937, 0]);
    assert_eq!(decoded[1].0, 2);
    assert_eq!(decoded[1].1[0], "PSEUDOGENE:1stops");
    assert_eq!(decoded[1].2, vec![102, 401, 992, 1001]);
    assert!(json.contains("\"audit_only_policy\"") && json.contains("\"correct_presence\""));
}

/// The cell tooltip was unstyled for a release: it sat static and transparent
/// under the scroller. Keep the rule, and keep legend keys as swatch+label pairs.
#[test]
fn report_styles_the_cell_tooltip_and_pairs_legend_keys() {
    let css = include_str!("report_assets/report.css");
    let rule = css
        .split(".matrix-tooltip {")
        .nth(1)
        .and_then(|rest| rest.split('}').next())
        .expect("a .matrix-tooltip rule");
    assert!(rule.contains("position: absolute"));
    assert!(rule.contains("pointer-events: none"));
    assert!(rule.contains("background: var(--ink)"));
    let writer = include_str!("output.rs");
    assert!(writer.contains("<span class='key'><i class='sw p'></i>input</span>"));
}

/// A row the matrix policy or a quality floor kept out was labelled
/// corrected_but_review_for_publication. Only a corrected cell may carry a
/// status that says corrected, and a span that is not a whole number of codons
/// is never publication-ready.
#[test]
fn publication_status_follows_the_matrix_action() {
    let states = [
        "STRICT_PRESENT",
        "ALIAS_CLUSTER",
        "SEQUENCE_ONLY",
        "PARTIAL_HOMOLOGY",
        "PSEUDOGENE",
        "HOMOLOG_PARALOG",
        "MISSING_CDS_UNVALIDATED",
    ];
    let actions = [
        "correct_presence",
        "audit_only_policy",
        "audit_only_quality",
        "audit_only",
    ];
    for state in states {
        for action in actions {
            for orf in ["in_frame_span", "out_of_frame_span_review", "disrupted_orf"] {
                for accepted in [true, false] {
                    let status = publication_status(accepted, state, action, orf);
                    assert!(PUBLICATION_STATUSES.contains(&status), "{status}");
                    if status.starts_with("corrected") || status.starts_with("publication") {
                        assert!(
                            accepted && action == "correct_presence",
                            "{state} {action} {status}"
                        );
                    }
                }
            }
        }
    }
    assert_eq!(
        publication_status(true, "PSEUDOGENE", "audit_only_policy", "disrupted_orf"),
        "held_out_policy"
    );
    assert_eq!(
        publication_status(
            true,
            "STRICT_PRESENT",
            "correct_presence",
            "out_of_frame_span_review"
        ),
        "corrected_review"
    );
    assert_eq!(
        publication_status(true, "STRICT_PRESENT", "correct_presence", "in_frame_span"),
        "publication_ready_correction"
    );
}

/// outputs.md lists the publication_status values; it listed four the code
/// never wrote. The documented set and the emitted set must be the same.
#[test]
fn documented_publication_statuses_match_the_code() {
    let docs = include_str!("../../docs/guide/outputs.md");
    let section = docs
        .split("`publication_status` takes one of these values:")
        .nth(1)
        .expect("the publication_status value table");
    let documented: Vec<&str> = section
        .lines()
        .skip_while(|line| !line.starts_with("| `"))
        .take_while(|line| line.starts_with("| `"))
        .map(|line| line.split('`').nth(1).unwrap())
        .collect();
    assert_eq!(documented, PUBLICATION_STATUSES.to_vec());
}

/// Gaps searched with nothing found are what true absences look like. A genome
/// with no correction must pass, and only a burden far above the cohort's asks
/// for review.
#[test]
fn genome_qc_does_not_review_genomes_without_corrections() {
    let out = tempfile::tempdir().unwrap();
    let n = 8;
    let genomes: Vec<String> = (0..n).map(|i| format!("g{i}")).collect();
    let clusters: Vec<String> = (0..40).map(|i| format!("c{i}")).collect();
    let cells = vec![vec![String::new(); n]; clusters.len()];
    let pa = PresenceAbsence::new_for_test(clusters, genomes.clone(), cells);
    let genome_files: HashMap<String, PathBuf> = genomes
        .iter()
        .map(|g| (g.clone(), PathBuf::from(format!("{g}.fa"))))
        .collect();
    // g0 has one correction, g7 has thirty, everyone else none.
    let mut rescues = Vec::new();
    for (genome, count) in [(0usize, 1usize), (7, 30)] {
        for ci in 0..count {
            let mut r = rescue_with_type(RescueType::Dropped, Some(true));
            r.cluster_idx = ClusterIdx(ci);
            r.genome_idx = GenomeIdx(genome);
            rescues.push(r);
        }
    }
    let matrix: HashMap<(ClusterIdx, GenomeIdx), &RescueEvent> = rescues
        .iter()
        .map(|r| ((r.cluster_idx, r.genome_idx), r))
        .collect();
    write_genome_qc(GenomeQcInputs {
        pa: &pa,
        genome_files: &genome_files,
        gff_indices: &HashMap::new(),
        gff_contig_remaps: &[],
        gaps_passing_filter_by_genome: &[40; 8],
        rescues: &rescues,
        synteny_rejected: &[],
        paralog_rejected: &[],
        matrix_rescue_map: &matrix,
        nsc_presence: &ahash::HashSet::default(),
        output_dir: out.path(),
    })
    .unwrap();
    let text = std::fs::read_to_string(out.path().join("genome_qc.tsv")).unwrap();
    let status: HashMap<&str, (&str, &str)> = text
        .lines()
        .skip(1)
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            (f[0], (f[10], f[11]))
        })
        .collect();
    for g in ["g1", "g2", "g3", "g4", "g5", "g6"] {
        assert_eq!(status[g].0, "pass", "{g}: {:?}", status[g]);
        assert!(status[g]
            .1
            .contains("no_rescue_evidence_for_searchable_gaps"));
    }
    assert_eq!(status["g0"].0, "pass");
    assert_eq!(status["g7"].0, "review");
    assert!(status["g7"].1.contains("correction_burden_above_cohort"));
}

fn tip(name: &str, genome: Option<usize>, length: f64) -> PhyloNode {
    PhyloNode {
        label: Some(name.into()),
        genome_idx: genome.map(GenomeIdx),
        branch_length: Some(length),
        children: Vec::new(),
    }
}

fn inner(children: Vec<PhyloNode>, length: Option<f64>) -> PhyloNode {
    PhyloNode {
        label: None,
        genome_idx: None,
        branch_length: length,
        children,
    }
}

/// A clade with no matrix genome is dropped, and a node left with one child
/// collapses into it with the lengths summed: a unary node became a nameless
/// tip, painted present in every core column, after a reroot in the browser.
#[test]
fn tree_json_drops_unmatched_tips_and_collapses_unary_nodes() {
    let tree = inner(
        vec![
            inner(vec![tip("a", Some(0), 1.0)], Some(2.0)),
            inner(vec![tip("unmatched", None, 1.0)], Some(3.0)),
            inner(
                vec![tip("b", Some(1), 0.5), tip("c", Some(2), 0.25)],
                Some(4.0),
            ),
        ],
        None,
    );
    let json = output::pruned_phylo_json(&tree).unwrap();
    assert!(!json.contains("unmatched"));
    assert!(json.contains("{\"label\":\"a\",\"genome\":0,\"length\":3,\"children\":[]}"));
    assert_eq!(json.matches("\"genome\":null").count(), 2, "{json}");
    // No node of the output has exactly one child.
    assert!(!json
        .contains("\"children\":[{\"label\":\"a\",\"genome\":0,\"length\":3,\"children\":[]}]}"));
    assert!(output::pruned_phylo_json(&tip("x", None, 1.0)).is_none());
}

/// A ladder tree is as deep as it has tips. The writer must not recurse: run it
/// on a 50,000-level ladder in a thread with a 2 MB stack.
#[test]
fn tree_json_handles_a_fifty_thousand_level_ladder() {
    let handle = std::thread::Builder::new()
        .stack_size(2 << 20)
        .spawn(|| {
            let depth = 50_000;
            let mut node = tip("t0", Some(0), 1.0);
            for i in 1..depth {
                node = inner(vec![node, tip(&format!("t{i}"), Some(i), 1.0)], Some(0.5));
            }
            let json = output::pruned_phylo_json(&node).unwrap();
            let tips = json.matches("\"children\":[]").count();
            // Take the tree apart by hand: dropping it would recurse too.
            let mut stack = vec![node];
            while let Some(mut n) = stack.pop() {
                stack.extend(std::mem::take(&mut n.children));
            }
            tips
        })
        .unwrap();
    assert_eq!(handle.join().unwrap(), 50_000);
}

/// A name holding "<!--" and then "<script" put the HTML parser into the
/// script's escaped state, so the viewer and the drawer ended up inside the JSON
/// and the page was blank with no error. The data block must hold neither, and
/// the boot and viewer scripts must still follow it.
#[test]
fn hostile_names_cannot_break_out_of_the_data_block() {
    let out = tempfile::tempdir().unwrap();
    let pa = PresenceAbsence::new_for_test(
        vec!["<!--<script>".into(), "a</script><b>&".into()],
        vec!["g1".into(), "<!--g2".into()],
        vec![
            vec!["gene1".into(), "".into()],
            vec!["".into(), "gene2".into()],
        ],
    );
    let mut args = default_args_for_validation();
    args.output = out.path().to_path_buf();
    args.write_html_report = true;
    let stats = SummaryStats {
        core_before: 0,
        total_gaps: 2,
        filtered_gaps: 2,
        skipped_by_prevalence: 0,
        searchable_gaps: 2,
        missing_fasta_gaps: 0,
        dropped: 0,
        missing_cds: 0,
        split: 0,
        fragmented: 0,
        alignment_split: 0,
        start_variant: 0,
        paralog: 0,
        pseudogene: 0,
        synteny_rejected: 0,
        paralog_rejected: 0,
        matrix_policy: MatrixPolicy::Balanced,
        matrix_corrected: 0,
        nsc_matrix_corrected: 0,
        synteny_enabled: false,
        nsc_clusters: 0,
        nsc_merged_genes: 0,
        skipped_no_reference: 0,
        remaining_gaps: 0,
    };
    write_html_report(
        &args,
        &pa,
        &stats,
        &HashMap::new(),
        &None,
        None,
        &[],
        out.path(),
    )
    .unwrap();
    let text = std::fs::read_to_string(out.path().join("recore_report.html")).unwrap();
    let start = text
        .find("<script id='explorerData' type='application/json'>")
        .unwrap();
    let body_start = start + "<script id='explorerData' type='application/json'>".len();
    let end = body_start + text[body_start..].find("</script>").unwrap();
    let block = &text[body_start..end];
    assert!(
        !block.contains('<') && !block.contains('>'),
        "raw angle bracket in the data block"
    );
    let after = &text[end..];
    assert!(after.contains("window.recoreExplorerFailed = function"));
    assert!(after.contains("const data = window.recoreExplorerData;"));
}
