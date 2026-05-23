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

#[test]
fn default_args_pass_validation() {
    let args = default_args_for_validation();
    validate_args(&args).unwrap();
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
    };

    write_html_report(&args, &pa, &stats, &HashMap::new(), &None, None, out.path()).unwrap();

    let text = std::fs::read_to_string(out.path().join("recore_report.html")).unwrap();
    assert!(text.contains("<title>Recore run report</title>"));
    assert!(text.contains("group_1"));
    assert!(text.contains("Orthology actions"));
    assert!(text.contains("Add --tree"));
    assert!(text.contains("gene_presence_absence_recore.csv"));
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
        out.path(),
    )
    .unwrap();

    let text = std::fs::read_to_string(out.path().join("recore_report.html")).unwrap();
    assert!(text.contains("Phylogeny and corrected matrix explorer"));
    assert!(text.contains("window.recoreExplorerData"));
    assert!(text.contains("clusterDetails"));
    assert!(text.contains("columnHeaderCanvas"));
    assert!(text.contains("cellSize"));
    assert!(text.contains("rowHeight"));
    assert!(text.contains("min='0.25'"));
    assert!(text.contains("step='0.25'"));
    assert!(text.contains("const headerH = 24"));
    assert!(text.contains("tipColorMode"));
    assert!(text.contains("matrixColorMode"));
    assert!(text.contains("processLegend"));
    assert!(text.contains("lineageLegend"));
    assert!(text.contains("viewer-controls"));
    assert!(text.contains("data-mechanism"));
    assert!(text.contains("data-clear-mechanism"));
    assert!(text.contains("referenceCell"));
    assert!(text.contains("targetReferenceGene"));
    assert!(text.contains("Reference equivalent"));
    assert!(text.contains("Color tips by lineage"));
    assert!(text.contains("#C94434"));
    assert!(text.contains("Color tips by selected process"));
    assert!(text.contains("Matrix by evolutionary process"));
    assert!(text.contains("clade_associated_gain"));
    assert!(text.contains("tipLabels"));
    assert!(text.contains("rebuildBaseMatrix"));
    assert!(text.contains("phylogramMode"));
    assert!(text.contains("const tipX = 14 + usable"));
    assert!(!text.contains("clusterHeaderLabel"));
    assert!(text.contains("referenceOnly"));
    assert!(text.contains("explorerFullscreen"));
    assert!(text.contains("\"referenceGenome\":null"));
    assert!(text.contains("\"states\":[\"PA\",\"AP\"]"));
    assert!(text.contains("matrixCanvas"));
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
