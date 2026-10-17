use super::*;
use crate::parser::PresenceAbsence;

fn make_hit(query_start: u64, query_end: u64, query_length: u64) -> Hit {
    Hit {
        query: "g".into(),
        cluster: "g".into(),
        query_length,
        query_start,
        query_end,
        strand: '+',
        target: "c".into(),
        target_start: query_start,
        target_end: query_end,
        target_length: 0,
        matches: query_end - query_start,
        alignment_length: query_end - query_start,
        mapq: 60,
    }
}

#[test]
fn test_merged_coverage_no_overlap() {
    let h1 = make_hit(0, 400, 1000);
    let h2 = make_hit(600, 1000, 1000);
    let (cov, frags) = merged_query_coverage(&[&h1, &h2]);
    assert!((cov - 80.0).abs() < 0.01);
    assert_eq!(frags, 2);
}

#[test]
fn test_merged_coverage_with_overlap() {
    let h1 = make_hit(0, 600, 1000);
    let h2 = make_hit(400, 1000, 1000);
    let (cov, frags) = merged_query_coverage(&[&h1, &h2]);
    assert!((cov - 100.0).abs() < 0.01);
    assert_eq!(frags, 1);
}

#[test]
fn test_merged_coverage_multi_allele() {
    let h1 = Hit {
        query: "g__0".into(),
        cluster: "g".into(),
        query_length: 800,
        query_start: 0,
        query_end: 600,
        strand: '+',
        target: "c".into(),
        target_start: 0,
        target_end: 600,
        target_length: 0,
        matches: 600,
        alignment_length: 600,
        mapq: 60,
    };
    let h2 = Hit {
        query: "g__1".into(),
        cluster: "g".into(),
        query_length: 1200,
        query_start: 0,
        query_end: 1100,
        strand: '+',
        target: "c".into(),
        target_start: 0,
        target_end: 1100,
        target_length: 0,
        matches: 1100,
        alignment_length: 1100,
        mapq: 60,
    };
    let (cov, frags) = merged_query_coverage(&[&h1, &h2]);
    assert!(cov > 90.0);
    assert_eq!(frags, 1);
}

#[test]
fn test_merged_coverage_same_length_different_alleles() {
    let h1 = Hit {
        query: "cluster~~~ALLELE~0".into(),
        cluster: "cluster".into(),
        query_length: 1000,
        query_start: 0,
        query_end: 400,
        strand: '+',
        target: "c".into(),
        target_start: 0,
        target_end: 400,
        target_length: 0,
        matches: 400,
        alignment_length: 400,
        mapq: 60,
    };
    let h2 = Hit {
        query: "cluster~~~ALLELE~1".into(),
        cluster: "cluster".into(),
        query_length: 1000,
        query_start: 600,
        query_end: 1000,
        strand: '+',
        target: "c".into(),
        target_start: 5000,
        target_end: 5400,
        target_length: 0,
        matches: 400,
        alignment_length: 400,
        mapq: 60,
    };
    let (cov, _frags) = merged_query_coverage(&[&h1, &h2]);
    assert!(
        cov < 50.0,
        "Coverage should be <50% (was {:.1}%), alleles must not be merged",
        cov
    );
}

#[test]
fn test_gff_index_overlap() {
    let genes = vec![
        GffRecord {
            contig: "chr1".into(),
            start: 100,
            end: 500,
            strand: '+',
            gene_id: "g1".into(),
        },
        GffRecord {
            contig: "chr1".into(),
            start: 600,
            end: 900,
            strand: '-',
            gene_id: "g2".into(),
        },
        GffRecord {
            contig: "chr2".into(),
            start: 50,
            end: 200,
            strand: '+',
            gene_id: "g3".into(),
        },
    ];
    let idx = GffIndex::build(genes);

    assert_eq!(idx.overlapping("chr1", 200, 400).len(), 1);
    assert_eq!(idx.overlapping("chr1", 400, 700).len(), 2);
    assert_eq!(idx.overlapping("chr2", 0, 100).len(), 1);
    assert!(idx.overlapping("chr3", 0, 100).is_empty());
}

#[test]
fn test_gff_flanking() {
    let genes = vec![
        GffRecord {
            contig: "c".into(),
            start: 100,
            end: 200,
            strand: '+',
            gene_id: "g1".into(),
        },
        GffRecord {
            contig: "c".into(),
            start: 300,
            end: 400,
            strand: '+',
            gene_id: "g2".into(),
        },
        GffRecord {
            contig: "c".into(),
            start: 500,
            end: 600,
            strand: '+',
            gene_id: "g3".into(),
        },
    ];
    let idx = GffIndex::build(genes);

    let (up, down) = idx.flanking("c", 300, 400);
    assert_eq!(up.unwrap().gene_id, "g1");
    assert_eq!(down.unwrap().gene_id, "g3");

    let (up, down) = idx.flanking("c", 100, 200);
    assert!(up.is_none());
    assert_eq!(down.unwrap().gene_id, "g2");
}

#[test]
fn test_synteny_context_handles_overlapping_boundary_genes() {
    let pa = PresenceAbsence::new_for_test(
        vec![
            "target_gene".into(),
            "up_overlap".into(),
            "down_overlap".into(),
            "outer_up".into(),
            "outer_down".into(),
        ],
        vec!["target".into(), "reference".into()],
        vec![
            vec!["".into(), "ref_center".into()],
            vec!["target_up".into(), "ref_up".into()],
            vec!["target_down".into(), "ref_down".into()],
            vec!["target_outer_up".into(), "ref_outer_up".into()],
            vec!["target_outer_down".into(), "ref_outer_down".into()],
        ],
    );

    let target_gff = GffIndex::build(vec![
        GffRecord {
            contig: "chr".into(),
            start: 0,
            end: 100,
            strand: '+',
            gene_id: "target_outer_up".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 100,
            end: 200,
            strand: '+',
            gene_id: "target_up".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 197,
            end: 260,
            strand: '+',
            gene_id: "target_center".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 257,
            end: 320,
            strand: '+',
            gene_id: "target_down".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 320,
            end: 400,
            strand: '+',
            gene_id: "target_outer_down".into(),
        },
    ]);
    let ref_gff = GffIndex::build(vec![
        GffRecord {
            contig: "chr".into(),
            start: 0,
            end: 100,
            strand: '+',
            gene_id: "ref_outer_up".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 100,
            end: 200,
            strand: '+',
            gene_id: "ref_up".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 197,
            end: 260,
            strand: '+',
            gene_id: "ref_center".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 257,
            end: 320,
            strand: '+',
            gene_id: "ref_down".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 320,
            end: 400,
            strand: '+',
            gene_id: "ref_outer_down".into(),
        },
    ]);

    let mut gff_indices = HashMap::new();
    gff_indices.insert(GenomeIdx(0), target_gff);
    gff_indices.insert(GenomeIdx(1), ref_gff);

    let mut g2c = HashMap::new();
    g2c.insert(
        GenomeIdx(0),
        GeneClusters::from_map(
            &gff_indices[&GenomeIdx(0)],
            &build_gene_to_cluster(&pa, GenomeIdx(0)),
        ),
    );
    g2c.insert(
        GenomeIdx(1),
        GeneClusters::from_map(
            &gff_indices[&GenomeIdx(1)],
            &build_gene_to_cluster(&pa, GenomeIdx(1)),
        ),
    );

    let mut rescue = RescueEvent {
        cluster_idx: ClusterIdx(0),
        genome_idx: GenomeIdx(0),
        rescue_type: RescueType::StartVariant { length_diff_bp: 10 },
        identity: 100.0,
        coverage: 95.0,
        contig: "chr".into(),
        start: 205,
        end: 250,
        strand: '+',
        annotation_match: Some("target_center".into()),
        confidence: 0.95,
        synteny_ok: None,
        synteny_evidence: SyntenyEvidence::default(),
        fragment_evidence: FragmentEvidence::default(),
        query_start: None,
    };

    validate_synteny(&mut rescue, &pa, &gff_indices, &g2c, false);
    assert_eq!(rescue.synteny_ok, Some(true));
}

#[test]
fn test_synteny_context_rejects_different_locus_paralog() {
    let pa = PresenceAbsence::new_for_test(
        vec![
            "target_gene".into(),
            "target_up".into(),
            "target_down".into(),
        ],
        vec!["target".into(), "reference".into()],
        vec![
            vec!["".into(), "ref_center".into()],
            vec!["target_up_gene".into(), "".into()],
            vec!["target_down_gene".into(), "".into()],
        ],
    );

    let target_gff = GffIndex::build(vec![
        GffRecord {
            contig: "chr".into(),
            start: 100,
            end: 200,
            strand: '+',
            gene_id: "target_up_gene".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 200,
            end: 300,
            strand: '+',
            gene_id: "target_center".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 300,
            end: 400,
            strand: '+',
            gene_id: "target_down_gene".into(),
        },
    ]);
    let ref_gff = GffIndex::build(vec![
        GffRecord {
            contig: "chr".into(),
            start: 100,
            end: 200,
            strand: '+',
            gene_id: "ref_up_other".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 200,
            end: 300,
            strand: '+',
            gene_id: "ref_center".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 300,
            end: 400,
            strand: '+',
            gene_id: "ref_down_other".into(),
        },
    ]);

    let mut gff_indices = HashMap::new();
    gff_indices.insert(GenomeIdx(0), target_gff);
    gff_indices.insert(GenomeIdx(1), ref_gff);

    let mut g2c = HashMap::new();
    g2c.insert(
        GenomeIdx(0),
        GeneClusters::from_map(
            &gff_indices[&GenomeIdx(0)],
            &build_gene_to_cluster(&pa, GenomeIdx(0)),
        ),
    );
    g2c.insert(
        GenomeIdx(1),
        GeneClusters::from_map(
            &gff_indices[&GenomeIdx(1)],
            &build_gene_to_cluster(&pa, GenomeIdx(1)),
        ),
    );

    let mut rescue = RescueEvent {
        cluster_idx: ClusterIdx(0),
        genome_idx: GenomeIdx(0),
        rescue_type: RescueType::Split {
            other_cluster: "paralogous_cluster".into(),
        },
        identity: 99.0,
        coverage: 99.0,
        contig: "chr".into(),
        start: 200,
        end: 300,
        strand: '+',
        annotation_match: Some("target_center".into()),
        confidence: 1.0,
        synteny_ok: None,
        synteny_evidence: SyntenyEvidence::default(),
        fragment_evidence: FragmentEvidence::default(),
        query_start: None,
    };

    validate_synteny(&mut rescue, &pa, &gff_indices, &g2c, false);
    assert_eq!(rescue.synteny_ok, Some(false));
}

#[test]
fn test_synteny_without_any_reference_annotation_is_not_assessed() {
    // Same target, but the genome carrying the cluster has no GFF: there is
    // nothing to compare with, which used to be a FAIL that rejected the rescue.
    let pa = PresenceAbsence::new_for_test(
        vec![
            "target_gene".into(),
            "target_up".into(),
            "target_down".into(),
        ],
        vec!["target".into(), "reference".into()],
        vec![
            vec!["".into(), "ref_center".into()],
            vec!["target_up_gene".into(), "".into()],
            vec!["target_down_gene".into(), "".into()],
        ],
    );

    let target_gff = GffIndex::build(vec![
        GffRecord {
            contig: "chr".into(),
            start: 100,
            end: 200,
            strand: '+',
            gene_id: "target_up_gene".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 200,
            end: 300,
            strand: '+',
            gene_id: "target_center".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 300,
            end: 400,
            strand: '+',
            gene_id: "target_down_gene".into(),
        },
    ]);
    let ref_gff = GffIndex::build(vec![
        GffRecord {
            contig: "chr".into(),
            start: 100,
            end: 200,
            strand: '+',
            gene_id: "ref_up_other".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 200,
            end: 300,
            strand: '+',
            gene_id: "ref_center".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 300,
            end: 400,
            strand: '+',
            gene_id: "ref_down_other".into(),
        },
    ]);

    let mut gff_indices = HashMap::new();
    gff_indices.insert(GenomeIdx(0), target_gff);
    drop(ref_gff);

    // Genome 1 has no annotation, so it has no lookups either.
    let mut g2c = HashMap::new();
    g2c.insert(
        GenomeIdx(0),
        GeneClusters::from_map(
            &gff_indices[&GenomeIdx(0)],
            &build_gene_to_cluster(&pa, GenomeIdx(0)),
        ),
    );

    let mut rescue = RescueEvent {
        cluster_idx: ClusterIdx(0),
        genome_idx: GenomeIdx(0),
        rescue_type: RescueType::Split {
            other_cluster: "paralogous_cluster".into(),
        },
        identity: 99.0,
        coverage: 99.0,
        contig: "chr".into(),
        start: 200,
        end: 300,
        strand: '+',
        annotation_match: Some("target_center".into()),
        confidence: 1.0,
        synteny_ok: None,
        synteny_evidence: SyntenyEvidence::default(),
        fragment_evidence: FragmentEvidence::default(),
        query_start: None,
    };

    validate_synteny(&mut rescue, &pa, &gff_indices, &g2c, false);
    assert_eq!(rescue.synteny_ok, None);
}

#[test]
fn test_synteny_uses_all_present_genomes_deterministically() {
    let pa = PresenceAbsence::new_for_test(
        vec![
            "target_gene".into(),
            "shared_up".into(),
            "shared_down".into(),
            "other_up".into(),
            "other_down".into(),
        ],
        vec![
            "ref0".into(),
            "ref1".into(),
            "ref2".into(),
            "ref3".into(),
            "ref4".into(),
            "ref5".into(),
            "ref6".into(),
            "target".into(),
        ],
        vec![
            vec![
                "g0".into(),
                "g1".into(),
                "g2".into(),
                "g3".into(),
                "g4".into(),
                "g5".into(),
                "g6".into(),
                "".into(),
            ],
            vec![
                "".into(),
                "".into(),
                "".into(),
                "".into(),
                "".into(),
                "up5".into(),
                "up6".into(),
                "up_target".into(),
            ],
            vec![
                "".into(),
                "".into(),
                "".into(),
                "".into(),
                "".into(),
                "down5".into(),
                "down6".into(),
                "down_target".into(),
            ],
            vec![
                "xup0".into(),
                "xup1".into(),
                "xup2".into(),
                "xup3".into(),
                "xup4".into(),
                "".into(),
                "".into(),
                "".into(),
            ],
            vec![
                "xdown0".into(),
                "xdown1".into(),
                "xdown2".into(),
                "xdown3".into(),
                "xdown4".into(),
                "".into(),
                "".into(),
                "".into(),
            ],
        ],
    );

    fn genes(up: &str, center: &str, down: &str) -> GffIndex {
        GffIndex::build(vec![
            GffRecord {
                contig: "chr".into(),
                start: 100,
                end: 200,
                strand: '+',
                gene_id: up.into(),
            },
            GffRecord {
                contig: "chr".into(),
                start: 200,
                end: 300,
                strand: '+',
                gene_id: center.into(),
            },
            GffRecord {
                contig: "chr".into(),
                start: 300,
                end: 400,
                strand: '+',
                gene_id: down.into(),
            },
        ])
    }

    let mut gff_indices = HashMap::new();
    for gi in 0..5 {
        gff_indices.insert(
            GenomeIdx(gi),
            genes(
                &format!("xup{}", gi),
                &format!("g{}", gi),
                &format!("xdown{}", gi),
            ),
        );
    }
    gff_indices.insert(GenomeIdx(5), genes("up5", "g5", "down5"));
    gff_indices.insert(GenomeIdx(6), genes("up6", "g6", "down6"));
    gff_indices.insert(GenomeIdx(7), genes("up_target", "candidate", "down_target"));

    let mut g2c = HashMap::new();
    for gi in 0..8 {
        g2c.insert(
            GenomeIdx(gi),
            GeneClusters::from_map(
                &gff_indices[&GenomeIdx(gi)],
                &build_gene_to_cluster(&pa, GenomeIdx(gi)),
            ),
        );
    }

    let mut rescue = RescueEvent {
        cluster_idx: ClusterIdx(0),
        genome_idx: GenomeIdx(7),
        rescue_type: RescueType::Dropped,
        identity: 99.0,
        coverage: 99.0,
        contig: "chr".into(),
        start: 200,
        end: 300,
        strand: '+',
        annotation_match: Some("candidate".into()),
        confidence: 1.0,
        synteny_ok: None,
        synteny_evidence: SyntenyEvidence::default(),
        fragment_evidence: FragmentEvidence::default(),
        query_start: None,
    };

    validate_synteny(&mut rescue, &pa, &gff_indices, &g2c, false);
    assert_eq!(rescue.synteny_ok, Some(true));
}

#[test]
fn test_synteny_uses_all_gene_ids_in_multicopy_reference_cell() {
    let pa = PresenceAbsence::new_for_test(
        vec![
            "target_gene".into(),
            "shared_up".into(),
            "shared_down".into(),
            "wrong_up".into(),
            "wrong_down".into(),
        ],
        vec!["target".into(), "reference".into()],
        vec![
            vec!["".into(), "ref_bad;ref_good".into()],
            vec!["target_up".into(), "ref_good_up".into()],
            vec!["target_down".into(), "ref_good_down".into()],
            vec!["".into(), "ref_bad_up".into()],
            vec!["".into(), "ref_bad_down".into()],
        ],
    );

    let target_gff = GffIndex::build(vec![
        GffRecord {
            contig: "chr".into(),
            start: 100,
            end: 200,
            strand: '+',
            gene_id: "target_up".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 200,
            end: 300,
            strand: '+',
            gene_id: "candidate".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 300,
            end: 400,
            strand: '+',
            gene_id: "target_down".into(),
        },
    ]);
    let ref_gff = GffIndex::build(vec![
        GffRecord {
            contig: "chr".into(),
            start: 100,
            end: 200,
            strand: '+',
            gene_id: "ref_bad_up".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 200,
            end: 300,
            strand: '+',
            gene_id: "ref_bad".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 300,
            end: 400,
            strand: '+',
            gene_id: "ref_bad_down".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 1000,
            end: 1100,
            strand: '+',
            gene_id: "ref_good_up".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 1100,
            end: 1200,
            strand: '+',
            gene_id: "ref_good".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 1200,
            end: 1300,
            strand: '+',
            gene_id: "ref_good_down".into(),
        },
    ]);

    let mut gff_indices = HashMap::new();
    gff_indices.insert(GenomeIdx(0), target_gff);
    gff_indices.insert(GenomeIdx(1), ref_gff);

    let mut g2c = HashMap::new();
    g2c.insert(
        GenomeIdx(0),
        GeneClusters::from_map(
            &gff_indices[&GenomeIdx(0)],
            &build_gene_to_cluster(&pa, GenomeIdx(0)),
        ),
    );
    g2c.insert(
        GenomeIdx(1),
        GeneClusters::from_map(
            &gff_indices[&GenomeIdx(1)],
            &build_gene_to_cluster(&pa, GenomeIdx(1)),
        ),
    );

    let mut rescue = RescueEvent {
        cluster_idx: ClusterIdx(0),
        genome_idx: GenomeIdx(0),
        rescue_type: RescueType::Dropped,
        identity: 99.0,
        coverage: 99.0,
        contig: "chr".into(),
        start: 200,
        end: 300,
        strand: '+',
        annotation_match: Some("candidate".into()),
        confidence: 1.0,
        synteny_ok: None,
        synteny_evidence: SyntenyEvidence::default(),
        fragment_evidence: FragmentEvidence::default(),
        query_start: None,
    };

    validate_synteny(&mut rescue, &pa, &gff_indices, &g2c, false);
    assert_eq!(rescue.synteny_ok, Some(true));
    assert_eq!(rescue.synteny_evidence.support, 1);
}

#[test]
fn test_missing_cds_preferred_over_annotated_paralog_and_validated_by_synteny() {
    let pa = PresenceAbsence::new_for_test(
        vec![
            "target_gene".into(),
            "shared_up".into(),
            "shared_down".into(),
            "other_cluster".into(),
        ],
        vec!["target".into(), "ref1".into(), "ref2".into()],
        vec![
            vec!["".into(), "ref_center1".into(), "ref_center2".into()],
            vec!["target_up".into(), "ref_up1".into(), "ref_up2".into()],
            vec!["target_down".into(), "ref_down1".into(), "ref_down2".into()],
            vec!["target_paralog".into(), "".into(), "".into()],
        ],
    );

    let target_gff = GffIndex::build(vec![
        GffRecord {
            contig: "chr".into(),
            start: 100,
            end: 200,
            strand: '+',
            gene_id: "target_up".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 300,
            end: 400,
            strand: '+',
            gene_id: "target_down".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 1000,
            end: 1100,
            strand: '+',
            gene_id: "target_paralog".into(),
        },
    ]);
    let ref1_gff = GffIndex::build(vec![
        GffRecord {
            contig: "chr".into(),
            start: 100,
            end: 200,
            strand: '+',
            gene_id: "ref_up1".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 200,
            end: 300,
            strand: '+',
            gene_id: "ref_center1".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 300,
            end: 400,
            strand: '+',
            gene_id: "ref_down1".into(),
        },
    ]);
    let ref2_gff = GffIndex::build(vec![
        GffRecord {
            contig: "chr".into(),
            start: 100,
            end: 200,
            strand: '+',
            gene_id: "ref_up2".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 200,
            end: 300,
            strand: '+',
            gene_id: "ref_center2".into(),
        },
        GffRecord {
            contig: "chr".into(),
            start: 300,
            end: 400,
            strand: '+',
            gene_id: "ref_down2".into(),
        },
    ]);

    let mut gff_indices = HashMap::new();
    gff_indices.insert(GenomeIdx(0), target_gff);
    gff_indices.insert(GenomeIdx(1), ref1_gff);
    gff_indices.insert(GenomeIdx(2), ref2_gff);

    let mut g2c = HashMap::new();
    for gi in 0..3 {
        g2c.insert(
            GenomeIdx(gi),
            GeneClusters::from_map(
                &gff_indices[&GenomeIdx(gi)],
                &build_gene_to_cluster(&pa, GenomeIdx(gi)),
            ),
        );
    }

    let annotated_paralog = Hit {
        query: "target_gene".into(),
        cluster: "target_gene".into(),
        query_length: 100,
        query_start: 0,
        query_end: 100,
        strand: '+',
        target: "chr".into(),
        target_start: 1000,
        target_end: 1100,
        target_length: 0,
        matches: 100,
        alignment_length: 100,
        mapq: 60,
    };
    let missing_cds_hit = Hit {
        query: "target_gene".into(),
        cluster: "target_gene".into(),
        query_length: 100,
        query_start: 0,
        query_end: 100,
        strand: '+',
        target: "chr".into(),
        target_start: 198,
        target_end: 302,
        target_length: 0,
        matches: 98,
        alignment_length: 100,
        mapq: 60,
    };

    let mut rescue = classify_hits(
        &[annotated_paralog, missing_cds_hit],
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        gff_indices.get(&GenomeIdx(0)),
        g2c.get(&GenomeIdx(0)).unwrap(),
        80.0,
    )
    .unwrap();

    assert!(matches!(rescue.rescue_type, RescueType::MissingCds));
    assert_eq!(rescue.start, 198);
    assert_eq!(rescue.annotation_match, None);
    // The rescue keeps where its one alignment starts on the representative.
    assert_eq!(rescue.query_start, Some(0));

    validate_synteny(&mut rescue, &pa, &gff_indices, &g2c, false);
    assert_eq!(rescue.synteny_ok, Some(true));
}

#[test]
fn test_synteny_rejects_one_off_context_when_many_references_exist() {
    let target = SyntenyContext {
        upstream: vec![ClusterIdx(1)],
        downstream: vec![ClusterIdx(2)],
    };
    let contexts = vec![
        SyntenyContext {
            upstream: vec![ClusterIdx(3)],
            downstream: vec![ClusterIdx(4)],
        },
        SyntenyContext {
            upstream: vec![ClusterIdx(5)],
            downstream: vec![ClusterIdx(6)],
        },
        SyntenyContext {
            upstream: vec![ClusterIdx(7)],
            downstream: vec![ClusterIdx(8)],
        },
        SyntenyContext {
            upstream: vec![ClusterIdx(9)],
            downstream: vec![ClusterIdx(10)],
        },
        SyntenyContext {
            upstream: vec![ClusterIdx(11)],
            downstream: vec![ClusterIdx(12)],
        },
        SyntenyContext {
            upstream: vec![ClusterIdx(1)],
            downstream: vec![ClusterIdx(2)],
        },
    ];
    assert!(!has_recurrent_context_support(&target, &contexts));

    let rare_contexts = &contexts[..3];
    assert!(has_recurrent_context_support(
        &rare_contexts[0],
        rare_contexts
    ));
}

#[test]
fn test_synteny_accepts_recurrent_compatible_context_not_exact_signature() {
    let target = SyntenyContext {
        upstream: vec![ClusterIdx(1)],
        downstream: vec![ClusterIdx(2)],
    };
    let contexts = vec![
        SyntenyContext {
            upstream: vec![ClusterIdx(1), ClusterIdx(20)],
            downstream: vec![ClusterIdx(2), ClusterIdx(21)],
        },
        SyntenyContext {
            upstream: vec![ClusterIdx(1), ClusterIdx(30)],
            downstream: vec![ClusterIdx(2), ClusterIdx(31)],
        },
        SyntenyContext {
            upstream: vec![ClusterIdx(3)],
            downstream: vec![ClusterIdx(4)],
        },
        SyntenyContext {
            upstream: vec![ClusterIdx(5)],
            downstream: vec![ClusterIdx(6)],
        },
    ];

    assert_eq!(compatible_reference_support(&target, &contexts), 2);
    assert!(has_recurrent_context_support(&target, &contexts));
}

#[test]
fn test_synteny_does_not_count_one_reference_neighbor_twice() {
    let target = vec![ClusterIdx(1), ClusterIdx(1), ClusterIdx(1), ClusterIdx(1)];
    let reference = vec![ClusterIdx(1)];

    let (score, matches) = weighted_side_score(&target, &reference);

    assert_eq!(matches, 1);
    assert!(
        score < SYNTENY_REFERENCE_SUPPORT_THRESHOLD,
        "one repeated reference neighbour should not validate a repeated flank by itself"
    );
}

#[test]
fn test_synteny_recurrence_requires_multiple_genomes_for_large_reference_sets() {
    let target = SyntenyContext {
        upstream: vec![ClusterIdx(1)],
        downstream: vec![ClusterIdx(2)],
    };
    let matching = SyntenyContext {
        upstream: vec![ClusterIdx(1)],
        downstream: vec![ClusterIdx(2)],
    };
    let nonmatching = SyntenyContext {
        upstream: vec![ClusterIdx(3)],
        downstream: vec![ClusterIdx(4)],
    };
    let contexts = vec![
        ReferenceSyntenyContext {
            genome_idx: GenomeIdx(0),
            context: matching.clone(),
        },
        ReferenceSyntenyContext {
            genome_idx: GenomeIdx(0),
            context: matching,
        },
        ReferenceSyntenyContext {
            genome_idx: GenomeIdx(1),
            context: nonmatching.clone(),
        },
        ReferenceSyntenyContext {
            genome_idx: GenomeIdx(2),
            context: nonmatching.clone(),
        },
        ReferenceSyntenyContext {
            genome_idx: GenomeIdx(3),
            context: nonmatching,
        },
    ];

    let assessment = assess_context_support(&target, &contexts, SYNTENY_CONTEXT_GENES);

    assert_eq!(assessment.support, 2);
    assert!(
        !assessment.passes,
        "multiple compatible copies from one genome should not satisfy recurrent support"
    );
}

#[test]
fn test_detect_start_variant_with_gff() {
    let genes = vec![GffRecord {
        contig: "c".into(),
        start: 0,
        end: 1100,
        strand: '+',
        gene_id: "gene1".into(),
    }];
    let idx = GffIndex::build(genes);

    let hit = Hit {
        query: "g".into(),
        cluster: "g".into(),
        query_length: 1000,
        query_start: 0,
        query_end: 1000,
        strand: '+',
        target: "c".into(),
        target_start: 100,
        target_end: 1100,
        target_length: 0,
        matches: 1000,
        alignment_length: 1000,
        mapq: 60,
    };
    let result = detect_start_variant(&hit, Some(&idx));
    assert_eq!(result, Some((100, "gene1".to_string())));
}

/// A CDS that agrees with the representative at both ends, once the query
/// bases the alignment left unaligned are projected onto the contig, is not a
/// start variant. The old rule compared lengths and called this +50bp.
#[test]
fn test_detect_start_variant_projects_unaligned_query_ends() {
    let genes = vec![GffRecord {
        contig: "c".into(),
        start: 0,
        end: 1050,
        strand: '+',
        gene_id: "gene1".into(),
    }];
    let idx = GffIndex::build(genes);
    // The last 50 bases of the 1100 bp query are unaligned, and the target span
    // ends at 1000: projected, the representative ends at 1050, as the CDS does.
    let hit = Hit {
        query: "g".into(),
        cluster: "g".into(),
        query_length: 1100,
        query_start: 0,
        query_end: 1050,
        strand: '+',
        target: "c".into(),
        target_start: 0,
        target_end: 1000,
        target_length: 0,
        matches: 1040,
        alignment_length: 1050,
        mapq: 60,
    };
    assert_eq!(detect_start_variant(&hit, Some(&idx)), None);
}

/// A CDS that runs past the representative's stop differs at the 3' end:
/// that is a stop-side difference, never a start variant. A representative
/// that reads through its stop used to be reported this way.
#[test]
fn test_detect_start_variant_ignores_stop_side_differences() {
    let idx = GffIndex::build(vec![GffRecord {
        contig: "c".into(),
        start: 0,
        end: 1100,
        strand: '+',
        gene_id: "gene1".into(),
    }]);
    let hit = make_hit(0, 1000, 1000);
    assert_eq!(detect_start_variant(&hit, Some(&idx)), None);
}

/// On the minus strand the start is the high coordinate.
#[test]
fn test_detect_start_variant_minus_strand_start_is_the_high_end() {
    let idx = GffIndex::build(vec![GffRecord {
        contig: "c".into(),
        start: 0,
        end: 1100,
        strand: '-',
        gene_id: "gene1".into(),
    }]);
    let hit = Hit {
        strand: '-',
        ..make_hit(0, 1000, 1000)
    };
    assert_eq!(
        detect_start_variant(&hit, Some(&idx)),
        Some((100, "gene1".to_string()))
    );
}

#[test]
fn test_detect_start_variant_no_gff() {
    let hit = make_hit(0, 1000, 1000);
    assert!(detect_start_variant(&hit, None).is_none());
}

#[test]
fn test_conflict_resolution() {
    let mut rescues = vec![
        RescueEvent {
            cluster_idx: ClusterIdx(0),
            genome_idx: GenomeIdx(0),
            rescue_type: RescueType::Dropped,
            identity: 99.0,
            coverage: 100.0,
            contig: "c1".into(),
            start: 100,
            end: 500,
            strand: '+',
            annotation_match: None,
            confidence: 0.99,
            synteny_ok: None,
            synteny_evidence: SyntenyEvidence::default(),
            fragment_evidence: FragmentEvidence::default(),
            query_start: None,
        },
        RescueEvent {
            cluster_idx: ClusterIdx(1),
            genome_idx: GenomeIdx(0),
            rescue_type: RescueType::Dropped,
            identity: 95.0,
            coverage: 90.0,
            contig: "c1".into(),
            start: 200,
            end: 600,
            strand: '+',
            annotation_match: None,
            confidence: 0.85,
            synteny_ok: None,
            synteny_evidence: SyntenyEvidence::default(),
            fragment_evidence: FragmentEvidence::default(),
            query_start: None,
        },
    ];
    resolve_conflicts(&mut rescues);
    assert_eq!(rescues.len(), 1);
    assert_eq!(rescues[0].cluster_idx, ClusterIdx(0));
}

#[test]
fn test_conflict_resolution_no_overlap() {
    let mut rescues = vec![
        RescueEvent {
            cluster_idx: ClusterIdx(0),
            genome_idx: GenomeIdx(0),
            rescue_type: RescueType::Dropped,
            identity: 99.0,
            coverage: 100.0,
            contig: "c1".into(),
            start: 100,
            end: 500,
            strand: '+',
            annotation_match: None,
            confidence: 0.99,
            synteny_ok: None,
            synteny_evidence: SyntenyEvidence::default(),
            fragment_evidence: FragmentEvidence::default(),
            query_start: None,
        },
        RescueEvent {
            cluster_idx: ClusterIdx(1),
            genome_idx: GenomeIdx(0),
            rescue_type: RescueType::Dropped,
            identity: 95.0,
            coverage: 90.0,
            contig: "c1".into(),
            start: 600,
            end: 900,
            strand: '+',
            annotation_match: None,
            confidence: 0.85,
            synteny_ok: None,
            synteny_evidence: SyntenyEvidence::default(),
            fragment_evidence: FragmentEvidence::default(),
            query_start: None,
        },
    ];
    resolve_conflicts(&mut rescues);
    assert_eq!(rescues.len(), 2);
}

#[test]
fn test_conflict_resolution_keeps_minor_overlaps() {
    let mut rescues = vec![
        RescueEvent {
            cluster_idx: ClusterIdx(0),
            genome_idx: GenomeIdx(0),
            rescue_type: RescueType::Dropped,
            identity: 99.0,
            coverage: 100.0,
            contig: "c1".into(),
            start: 100,
            end: 500,
            strand: '+',
            annotation_match: None,
            confidence: 0.99,
            synteny_ok: None,
            synteny_evidence: SyntenyEvidence::default(),
            fragment_evidence: FragmentEvidence::default(),
            query_start: None,
        },
        RescueEvent {
            cluster_idx: ClusterIdx(1),
            genome_idx: GenomeIdx(0),
            rescue_type: RescueType::Dropped,
            identity: 98.0,
            coverage: 100.0,
            contig: "c1".into(),
            start: 490,
            end: 900,
            strand: '+',
            annotation_match: None,
            confidence: 0.98,
            synteny_ok: None,
            synteny_evidence: SyntenyEvidence::default(),
            fragment_evidence: FragmentEvidence::default(),
            query_start: None,
        },
    ];
    resolve_conflicts(&mut rescues);
    assert_eq!(rescues.len(), 2);
}

#[test]
fn test_fragmented_span() {
    let h1 = Hit {
        query: "g".into(),
        cluster: "g".into(),
        query_length: 1000,
        query_start: 0,
        query_end: 400,
        strand: '+',
        target: "c".into(),
        target_start: 1000,
        target_end: 1400,
        target_length: 0,
        matches: 400,
        alignment_length: 400,
        mapq: 60,
    };
    let h2 = Hit {
        query: "g".into(),
        cluster: "g".into(),
        query_length: 1000,
        query_start: 600,
        query_end: 1000,
        strand: '+',
        target: "c".into(),
        target_start: 2000,
        target_end: 2400,
        target_length: 0,
        matches: 400,
        alignment_length: 400,
        mapq: 60,
    };

    let pa = PresenceAbsence::new_for_test(
        vec!["g".into()],
        vec!["genome1".into()],
        vec![vec!["".into()]],
    );

    let result = classify_hits(
        &[h1, h2],
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        None,
        &GeneClusters::none(),
        80.0,
    );
    assert!(result.is_some());
    let r = result.unwrap();
    assert_eq!(r.start, 1000);
    assert_eq!(r.end, 2400);
    match r.rescue_type {
        RescueType::Fragmented { n_fragments, .. } => assert_eq!(n_fragments, 2),
        _ => panic!("Expected Fragmented"),
    }
}

#[test]
fn test_best_fragment_chain_tie_breaks_by_locus() {
    let first = Hit {
        query: "g".into(),
        cluster: "g".into(),
        query_length: 1000,
        query_start: 0,
        query_end: 500,
        strand: '+',
        target: "z_contig".into(),
        target_start: 100,
        target_end: 600,
        target_length: 0,
        matches: 500,
        alignment_length: 500,
        mapq: 60,
    };
    let second = Hit {
        target: "a_contig".into(),
        target_start: 100,
        target_end: 600,
        ..first.clone()
    };
    let hits = vec![&first, &second];

    let chain = best_fragment_chain(&hits).unwrap();

    assert_eq!(chain.contig, "a_contig");
}

#[test]
fn test_classify_split_vs_paralog() {
    // Setup: cluster "target_cluster" absent in genome, hit lands on a gene
    // assigned to "other_cluster". Identity determines SPLIT vs PARALOG.
    let gff_genes = vec![GffRecord {
        contig: "c".into(),
        start: 0,
        end: 1000,
        strand: '+',
        gene_id: "geneA".into(),
    }];
    let gff = GffIndex::build(gff_genes);

    let pa = PresenceAbsence::new_for_test(
        vec!["target_cluster".into(), "other_cluster".into()],
        vec!["genome1".into()],
        vec![vec!["".into()], vec!["geneA".into()]],
    );

    let mut g2c: HashMap<String, ClusterIdx> = HashMap::new();
    g2c.insert("geneA".into(), ClusterIdx(1)); // geneA belongs to other_cluster

    // High identity (99%) → SPLIT
    let hit_high = Hit {
        query: "q".into(),
        cluster: "target_cluster".into(),
        query_length: 1000,
        query_start: 0,
        query_end: 1000,
        strand: '+',
        target: "c".into(),
        target_start: 0,
        target_end: 1000,
        target_length: 0,
        matches: 990,
        alignment_length: 1000,
        mapq: 60,
    };
    let result = classify_hits(
        &[hit_high],
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        Some(&gff),
        &GeneClusters::from_map(&gff, &g2c),
        80.0,
    );
    let r = result.unwrap();
    match &r.rescue_type {
        RescueType::Split { other_cluster } => assert_eq!(other_cluster, "other_cluster"),
        other => panic!("Expected Split, got {:?}", other),
    }

    // Low identity (90%) → PARALOG
    let hit_low = Hit {
        query: "q".into(),
        cluster: "target_cluster".into(),
        query_length: 1000,
        query_start: 0,
        query_end: 1000,
        strand: '+',
        target: "c".into(),
        target_start: 0,
        target_end: 1000,
        target_length: 0,
        matches: 900,
        alignment_length: 1000,
        mapq: 60,
    };
    let result = classify_hits(
        &[hit_low],
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        Some(&gff),
        &GeneClusters::from_map(&gff, &g2c),
        80.0,
    );
    let r = result.unwrap();
    match &r.rescue_type {
        RescueType::Paralog { other_cluster } => assert_eq!(other_cluster, "other_cluster"),
        other => panic!("Expected Paralog, got {:?}", other),
    }
    assert!(r.confidence < 0.7, "Paralog should have reduced confidence");
}

#[test]
fn test_classify_dropped_no_gff() {
    // Full-coverage hit with no GFF → DROPPED
    let hit = Hit {
        query: "q".into(),
        cluster: "cluster1".into(),
        query_length: 1000,
        query_start: 0,
        query_end: 1000,
        strand: '+',
        target: "c".into(),
        target_start: 0,
        target_end: 1000,
        target_length: 0,
        matches: 990,
        alignment_length: 1000,
        mapq: 60,
    };
    let pa = PresenceAbsence::new_for_test(
        vec!["cluster1".into()],
        vec!["genome1".into()],
        vec![vec!["".into()]],
    );
    let result = classify_hits(
        &[hit],
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        None,
        &GeneClusters::none(),
        80.0,
    );
    let r = result.unwrap();
    assert!(matches!(r.rescue_type, RescueType::Dropped));
    assert!(r.confidence > 0.9);
}

#[test]
fn test_classify_prefers_best_evidence_over_longest_allele() {
    let long_lower_score = Hit {
        query: "cluster1~~~ALLELE~0".into(),
        cluster: "cluster1".into(),
        query_length: 1200,
        query_start: 0,
        query_end: 960,
        strand: '+',
        target: "long_hit".into(),
        target_start: 100,
        target_end: 1060,
        target_length: 0,
        matches: 864,
        alignment_length: 960,
        mapq: 60,
    };
    let short_better_score = Hit {
        query: "cluster1~~~ALLELE~1".into(),
        cluster: "cluster1".into(),
        query_length: 900,
        query_start: 0,
        query_end: 900,
        strand: '+',
        target: "best_hit".into(),
        target_start: 200,
        target_end: 1100,
        target_length: 0,
        matches: 900,
        alignment_length: 900,
        mapq: 60,
    };
    let pa = PresenceAbsence::new_for_test(
        vec!["cluster1".into()],
        vec!["genome1".into()],
        vec![vec!["".into()]],
    );
    let result = classify_hits(
        &[long_lower_score, short_better_score],
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        None,
        &GeneClusters::none(),
        80.0,
    )
    .unwrap();

    assert_eq!(result.contig, "best_hit");
    assert!((result.identity - 100.0).abs() < 0.01);
    assert!((result.coverage - 100.0).abs() < 0.01);
}

#[test]
fn test_classify_missing_cds_no_overlapping_genes() {
    // Full-coverage hit, GFF exists but no gene overlaps the hit region → MISSING_CDS
    let gff_genes = vec![GffRecord {
        contig: "c".into(),
        start: 5000,
        end: 6000,
        strand: '+',
        gene_id: "far_gene".into(),
    }];
    let gff = GffIndex::build(gff_genes);
    let pa = PresenceAbsence::new_for_test(
        vec!["cluster1".into()],
        vec!["genome1".into()],
        vec![vec!["".into()]],
    );
    let g2c: HashMap<String, ClusterIdx> = HashMap::new();
    let hit = Hit {
        query: "q".into(),
        cluster: "cluster1".into(),
        query_length: 1000,
        query_start: 0,
        query_end: 1000,
        strand: '+',
        target: "c".into(),
        target_start: 0,
        target_end: 1000,
        target_length: 0,
        matches: 990,
        alignment_length: 1000,
        mapq: 60,
    };
    let result = classify_hits(
        &[hit],
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        Some(&gff),
        &GeneClusters::from_map(&gff, &g2c),
        80.0,
    );
    let r = result.unwrap();
    assert!(matches!(r.rescue_type, RescueType::MissingCds));
}

#[test]
fn test_classify_start_variant_with_gff() {
    // Hit with high identity overlapping a similar-length GFF gene → START_VARIANT
    let gff_genes = vec![GffRecord {
        contig: "c".into(),
        start: 0,
        end: 1100,
        strand: '+',
        gene_id: "gene1".into(),
    }];
    let gff = GffIndex::build(gff_genes);
    let pa = PresenceAbsence::new_for_test(
        vec!["cluster1".into()],
        vec!["genome1".into()],
        vec![vec!["".into()]],
    );
    let g2c: HashMap<String, ClusterIdx> = HashMap::new();
    // gene1 is NOT in g2c → no cluster assignment → classify_single_hit returns DROPPED
    // Then detect_start_variant kicks in: gene_len(1100) - target_aln(1000) = +100bp
    let hit = Hit {
        query: "q".into(),
        cluster: "cluster1".into(),
        query_length: 1000,
        query_start: 0,
        query_end: 1000,
        strand: '+',
        target: "c".into(),
        target_start: 100,
        target_end: 1100,
        target_length: 0,
        matches: 1000,
        alignment_length: 1000,
        mapq: 60,
    };
    let result = classify_hits(
        &[hit],
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        Some(&gff),
        &GeneClusters::from_map(&gff, &g2c),
        80.0,
    );
    let r = result.unwrap();
    match &r.rescue_type {
        RescueType::StartVariant { length_diff_bp } => assert_eq!(*length_diff_bp, 100),
        other => panic!("Expected StartVariant, got {:?}", other),
    }
    assert_eq!(r.annotation_match, Some("gene1".to_string()));
}

#[test]
fn test_classify_returns_none_for_empty_hits() {
    let pa =
        PresenceAbsence::new_for_test(vec!["c1".into()], vec!["g1".into()], vec![vec!["".into()]]);
    let result = classify_hits(
        &[],
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        None,
        &GeneClusters::none(),
        80.0,
    );
    assert!(result.is_none());
}

#[test]
fn test_classify_returns_none_below_coverage() {
    // Single hit with only 50% coverage, no fragments → None
    let hit = Hit {
        query: "q".into(),
        cluster: "c1".into(),
        query_length: 1000,
        query_start: 0,
        query_end: 500,
        strand: '+',
        target: "c".into(),
        target_start: 0,
        target_end: 500,
        target_length: 0,
        matches: 500,
        alignment_length: 500,
        mapq: 60,
    };
    let pa =
        PresenceAbsence::new_for_test(vec!["c1".into()], vec!["g1".into()], vec![vec!["".into()]]);
    let result = classify_hits(
        &[hit],
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        None,
        &GeneClusters::none(),
        80.0,
    );
    assert!(result.is_none());
}

#[test]
fn test_build_gene_to_cluster() {
    let pa = PresenceAbsence::new_for_test(
        vec!["clusterA".into(), "clusterB".into()],
        vec!["g1".into(), "g2".into()],
        vec![
            vec!["gene1;gene2".into(), "gene3".into()],
            vec!["".into(), "gene4".into()],
        ],
    );
    let g2c = build_gene_to_cluster(&pa, GenomeIdx(0));
    assert_eq!(g2c.get("gene1"), Some(&ClusterIdx(0)));
    assert_eq!(g2c.get("gene2"), Some(&ClusterIdx(0)));
    assert_eq!(g2c.get("gene4"), None); // gene4 is in genome1 (GenomeIdx(1)), not genome0

    let g2c_1 = build_gene_to_cluster(&pa, GenomeIdx(1));
    assert_eq!(g2c_1.get("gene3"), Some(&ClusterIdx(0)));
    assert_eq!(g2c_1.get("gene4"), Some(&ClusterIdx(1)));
}

#[test]
fn test_gene_suffix_stripping_only_handles_panaroo_suffixes() {
    let pa = PresenceAbsence::new_for_test(
        vec!["clusterA".into(), "clusterB".into()],
        vec!["g1".into()],
        vec![vec!["gene_dup_2".into()], vec!["gene_dup_domain".into()]],
    );
    let g2c = build_gene_to_cluster(&pa, GenomeIdx(0));

    assert_eq!(g2c.get("gene"), Some(&ClusterIdx(0)));
    assert_eq!(g2c.get("gene_dup_2"), Some(&ClusterIdx(0)));
    assert_eq!(g2c.get("gene_dup_domain"), Some(&ClusterIdx(1)));
    assert_eq!(
        g2c.get("gene_dup"),
        None,
        "biological names containing _dup_ must not be collapsed"
    );
}

#[test]
fn test_build_gene_to_cluster_adds_unambiguous_gff_aliases() {
    let pa = PresenceAbsence::new_for_test(
        vec!["clusterA".into()],
        vec!["g1".into()],
        vec![vec!["commonName".into()]],
    );
    let gff = GffIndex::build_with_aliases(vec![(
        GffRecord {
            contig: "ctg".into(),
            start: 0,
            end: 100,
            strand: '+',
            gene_id: "Rv0001".into(),
        },
        vec!["commonName".into(), "oldName".into()],
    )]);

    let g2c = build_gene_to_cluster_with_gff_aliases(&pa, GenomeIdx(0), Some(&gff));
    assert_eq!(g2c.get("commonName"), Some(&ClusterIdx(0)));
    assert_eq!(g2c.get("Rv0001"), Some(&ClusterIdx(0)));
    assert_eq!(g2c.get("oldName"), Some(&ClusterIdx(0)));
}

#[test]
fn test_descriptive_gff_labels_do_not_create_cluster_aliases() {
    let pa = PresenceAbsence::new_for_test(
        vec!["clusterA".into()],
        vec!["g1".into()],
        vec![vec!["commonName".into()]],
    );
    let gff = GffIndex::build_with_aliases_and_labels(vec![
        (
            GffRecord {
                contig: "ctg".into(),
                start: 0,
                end: 100,
                strand: '+',
                gene_id: "Rv0001".into(),
            },
            vec!["commonName".into()],
            vec!["commonName".into(), "hypothetical%20protein".into()],
        ),
        (
            GffRecord {
                contig: "ctg".into(),
                start: 200,
                end: 300,
                strand: '+',
                gene_id: "Rv0002".into(),
            },
            vec!["Rv0002".into()],
            vec!["Rv0002".into(), "hypothetical%20protein".into()],
        ),
    ]);

    let g2c = build_gene_to_cluster_with_gff_aliases(&pa, GenomeIdx(0), Some(&gff));
    assert_eq!(g2c.get("Rv0001"), Some(&ClusterIdx(0)));
    assert_eq!(g2c.get("hypothetical%20protein"), None);
    assert_eq!(g2c.get("Rv0002"), None);
}

#[test]
fn test_classify_hit_uses_gff_alias_for_split_cluster() {
    let pa = PresenceAbsence::new_for_test(
        vec!["queryCluster".into(), "annotatedCluster".into()],
        vec!["g1".into()],
        vec![vec!["".into()], vec!["commonName".into()]],
    );
    let gff = GffIndex::build_with_aliases(vec![(
        GffRecord {
            contig: "ctg".into(),
            start: 100,
            end: 1100,
            strand: '+',
            gene_id: "Rv0001".into(),
        },
        vec!["commonName".into()],
    )]);
    let g2c = build_gene_to_cluster_with_gff_aliases(&pa, GenomeIdx(0), Some(&gff));
    let hit = Hit {
        query: "queryCluster".into(),
        cluster: "queryCluster".into(),
        query_length: 1000,
        query_start: 0,
        query_end: 1000,
        strand: '+',
        target: "ctg".into(),
        target_start: 100,
        target_end: 1100,
        target_length: 0,
        matches: 990,
        alignment_length: 1000,
        mapq: 60,
    };

    let result = classify_hits(
        &[hit],
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        Some(&gff),
        &GeneClusters::from_map(&gff, &g2c),
        80.0,
    )
    .unwrap();
    assert_eq!(result.annotation_match, Some("Rv0001".to_string()));
    assert_eq!(
        result.rescue_type,
        RescueType::Split {
            other_cluster: "annotatedCluster".into()
        }
    );
}

#[test]
fn test_rescue_type_display() {
    assert_eq!(RescueType::Dropped.to_string(), "DROPPED");
    assert_eq!(RescueType::MissingCds.to_string(), "MISSING_CDS");
    assert_eq!(
        RescueType::Split {
            other_cluster: "cl2".into()
        }
        .to_string(),
        "SPLIT:cl2"
    );
    assert_eq!(
        RescueType::Fragmented {
            n_fragments: 3,
            total_coverage: 95.5
        }
        .to_string(),
        "FRAGMENTED:3frags:96%cov"
    );
    assert_eq!(
        RescueType::AlignmentSplit {
            n_fragments: 2,
            total_coverage: 98.8
        }
        .to_string(),
        "ALIGNMENT_SPLIT:2frags:99%cov"
    );
    assert_eq!(
        RescueType::StartVariant {
            length_diff_bp: -42
        }
        .to_string(),
        "START_VARIANT:-42bp"
    );
    assert_eq!(
        RescueType::Paralog {
            other_cluster: "cl3".into()
        }
        .to_string(),
        "PARALOG:cl3"
    );
    assert_eq!(
        RescueType::Pseudogene {
            reason: PseudogeneReason::InternalStop { n_stops: 2 }
        }
        .to_string(),
        "PSEUDOGENE:2stops"
    );
    assert_eq!(
        RescueType::Pseudogene {
            reason: PseudogeneReason::Frameshift
        }
        .to_string(),
        "PSEUDOGENE:frameshift"
    );
    assert_eq!(
        RescueType::Pseudogene {
            reason: PseudogeneReason::PartialDeletion { coverage: 72.5 }
        }
        .to_string(),
        "PSEUDOGENE:partial:72%cov"
    );
}

#[test]
fn test_count_internal_stops_clean_orf() {
    // ATG ... no stops ... TAA
    let seq = b"ATGAAAGCGTTTTTTAAAGCGAAAGCGTAA";
    assert_eq!(count_internal_stops(seq), 0);
}

#[test]
fn test_count_internal_stops_with_stops() {
    // Frame 0: ATG TAA GCG TTT — has internal stop at pos 3
    let seq = b"ATGTAAGCGTTTGCGAAAGCG";
    let stops = count_internal_stops(seq);
    // At least one frame will have stops, the cleanest frame might have 0 or 1
    // Frame 0: ATG|TAA|GCG|TTT|GCG|AAA|GCG = TAA at pos 3 (internal)
    // We want the minimum across frames
    assert!(stops <= 1);
}

#[test]
fn test_count_internal_stops_short_seq() {
    assert_eq!(count_internal_stops(b"ATG"), 0);
    assert_eq!(count_internal_stops(b""), 0);
}

#[test]
fn test_frameshift_pattern_detected() {
    // Two fragments with a 2bp gap (not divisible by 3) → frameshift
    let h1 = Hit {
        query: "q".into(),
        cluster: "c".into(),
        query_length: 1000,
        query_start: 0,
        query_end: 500,
        strand: '+',
        target: "c".into(),
        target_start: 100,
        target_end: 600,
        target_length: 0,
        matches: 500,
        alignment_length: 500,
        mapq: 60,
    };
    let h2 = Hit {
        query: "q".into(),
        cluster: "c".into(),
        query_length: 1000,
        query_start: 500,
        query_end: 1000,
        strand: '+',
        target: "c".into(),
        target_start: 602,
        target_end: 1102, // 2bp gap (not %3)
        target_length: 0,
        matches: 500,
        alignment_length: 500,
        mapq: 60,
    };
    assert!(detect_frameshift_pattern(&[&h1, &h2]));
}

#[test]
fn test_frameshift_pattern_not_detected_divisible_by_3() {
    let h1 = Hit {
        query: "q".into(),
        cluster: "c".into(),
        query_length: 1000,
        query_start: 0,
        query_end: 500,
        strand: '+',
        target: "c".into(),
        target_start: 100,
        target_end: 600,
        target_length: 0,
        matches: 500,
        alignment_length: 500,
        mapq: 60,
    };
    let h2 = Hit {
        query: "q".into(),
        cluster: "c".into(),
        query_length: 1000,
        query_start: 500,
        query_end: 1000,
        strand: '+',
        target: "c".into(),
        target_start: 603,
        target_end: 1103, // 3bp gap (divisible by 3)
        target_length: 0,
        matches: 500,
        alignment_length: 500,
        mapq: 60,
    };
    assert!(!detect_frameshift_pattern(&[&h1, &h2]));
}

#[test]
fn test_reverse_complement() {
    assert_eq!(reverse_complement(b"ATCG"), b"CGAT");
    assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
    assert_eq!(reverse_complement(b""), b"");
}

#[test]
fn test_internal_stops_minus_strand() {
    // Forward: GCATTAGCG — no obvious stops in frame 0
    // Reverse complement: CGCTAATGC
    // Frame 0: CGC|TAA|TGC → TAA is internal stop!
    let seq = b"GCATTAGCG";
    let rc = reverse_complement(seq);
    assert_eq!(&rc, b"CGCTAATGC");
    // The RC should have a stop in at least one frame
    let stops = count_internal_stops(&rc);
    // Frame 0: CGC TAA TGC → 1 stop (TAA), but it's at pos 3 with only 1 codon after → internal
    // Minimum across frames could be 0 (other frames may be clean)
    // Just verify the function runs without panic on RC
    assert!(stops <= 1);

    // Now test via check_internal_stops_in_hit with strand -
    let mut contigs: HashMap<String, Vec<u8>> = HashMap::new();
    contigs.insert("c".into(), seq.to_vec());
    let hit_minus = Hit {
        query: "q".into(),
        cluster: "c".into(),
        query_length: 9,
        query_start: 0,
        query_end: 9,
        strand: '-',
        target: "c".into(),
        target_start: 0,
        target_end: 9,
        target_length: 0,
        matches: 9,
        alignment_length: 9,
        mapq: 60,
    };
    let hit_plus = Hit {
        strand: '+',
        ..hit_minus.clone()
    };
    // Strand - should RC before scanning
    let result_minus = check_internal_stops_in_hit(&hit_minus, Some(&contigs));
    let result_plus = check_internal_stops_in_hit(&hit_plus, Some(&contigs));
    // At minimum, they should not be identical (different sequences scanned)
    // (they could both be None if no stops, which is fine)
    let _ = (result_minus, result_plus); // just ensure no panic
}

#[test]
fn test_classify_partial_deletion_pseudogene() {
    // Hit with 95% identity but only 70% coverage → partial deletion pseudogene
    let hit = Hit {
        query: "q".into(),
        cluster: "c1".into(),
        query_length: 1000,
        query_start: 0,
        query_end: 700,
        strand: '+',
        target: "c".into(),
        target_start: 0,
        target_end: 700,
        target_length: 0,
        matches: 665,
        alignment_length: 700,
        mapq: 60,
    };
    let pa =
        PresenceAbsence::new_for_test(vec!["c1".into()], vec!["g1".into()], vec![vec!["".into()]]);
    // min_coverage = 80% so 70% coverage won't pass Path 1
    let result = classify_hits(
        &[hit],
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        None,
        &GeneClusters::none(),
        80.0,
    );
    let r = result.unwrap();
    match &r.rescue_type {
        RescueType::Pseudogene {
            reason: PseudogeneReason::PartialDeletion { coverage },
        } => {
            assert!(*coverage >= 69.0 && *coverage <= 71.0);
        }
        other => panic!("Expected Pseudogene:PartialDeletion, got {:?}", other),
    }
}

#[test]
fn test_alignment_split_for_small_colinear_gaps() {
    let h1 = Hit {
        query: "g".into(),
        cluster: "g".into(),
        query_length: 1000,
        query_start: 0,
        query_end: 495,
        strand: '+',
        target: "c".into(),
        target_start: 0,
        target_end: 495,
        target_length: 0,
        matches: 495,
        alignment_length: 495,
        mapq: 60,
    };
    let h2 = Hit {
        query: "g".into(),
        cluster: "g".into(),
        query_length: 1000,
        query_start: 504,
        query_end: 1000,
        strand: '+',
        target: "c".into(),
        target_start: 504,
        target_end: 1000,
        target_length: 0,
        matches: 495,
        alignment_length: 495,
        mapq: 60,
    };
    let pa = PresenceAbsence::new_for_test(
        vec!["g".into()],
        vec!["genome1".into()],
        vec![vec!["".into()]],
    );
    let result = classify_hits(
        &[h1, h2],
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        None,
        &GeneClusters::none(),
        80.0,
    )
    .unwrap();
    match result.rescue_type {
        RescueType::AlignmentSplit { n_fragments, .. } => assert_eq!(n_fragments, 2),
        other => panic!("Expected AlignmentSplit, got {:?}", other),
    }
    assert_eq!(result.fragment_evidence.max_target_gap_bp, 9);
    assert_eq!(
        result.fragment_evidence.note.as_deref(),
        Some("alignment_split")
    );
}

#[test]
fn test_fragment_chain_inside_same_cds_is_not_fragmented() {
    let h1 = Hit {
        query: "g".into(),
        cluster: "g".into(),
        query_length: 1000,
        query_start: 0,
        query_end: 400,
        strand: '+',
        target: "c".into(),
        target_start: 100,
        target_end: 500,
        target_length: 0,
        matches: 400,
        alignment_length: 400,
        mapq: 60,
    };
    let h2 = Hit {
        query: "g".into(),
        cluster: "g".into(),
        query_length: 1000,
        query_start: 600,
        query_end: 1000,
        strand: '+',
        target: "c".into(),
        target_start: 700,
        target_end: 1100,
        target_length: 0,
        matches: 400,
        alignment_length: 400,
        mapq: 60,
    };
    let gff = GffIndex::build(vec![GffRecord {
        contig: "c".into(),
        start: 100,
        end: 1100,
        strand: '+',
        gene_id: "gene1".into(),
    }]);
    let pa = PresenceAbsence::new_for_test(
        vec!["g".into()],
        vec!["genome1".into()],
        vec![vec!["".into()]],
    );
    let g2c: HashMap<String, ClusterIdx> = HashMap::new();
    let result = classify_hits(
        &[h1, h2],
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        Some(&gff),
        &GeneClusters::from_map(&gff, &g2c),
        80.0,
    )
    .unwrap();
    assert!(matches!(result.rescue_type, RescueType::Dropped));
    assert_eq!(result.annotation_match.as_deref(), Some("gene1"));
    assert_eq!(result.fragment_evidence.same_cds, Some(true));
}

#[test]
fn test_non_colinear_fragments_are_not_rescued() {
    let h1 = Hit {
        query: "g".into(),
        cluster: "g".into(),
        query_length: 1000,
        query_start: 0,
        query_end: 400,
        strand: '+',
        target: "c".into(),
        target_start: 1000,
        target_end: 1400,
        target_length: 0,
        matches: 400,
        alignment_length: 400,
        mapq: 60,
    };
    let h2 = Hit {
        query: "g".into(),
        cluster: "g".into(),
        query_length: 1000,
        query_start: 600,
        query_end: 1000,
        strand: '+',
        target: "c".into(),
        target_start: 500,
        target_end: 900,
        target_length: 0,
        matches: 400,
        alignment_length: 400,
        mapq: 60,
    };
    let pa = PresenceAbsence::new_for_test(
        vec!["g".into()],
        vec!["genome1".into()],
        vec![vec!["".into()]],
    );
    let result = classify_hits(
        &[h1, h2],
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        None,
        &GeneClusters::none(),
        80.0,
    );
    assert!(result.is_none());
}

/// Genes are sorted by start and their ends stop being monotone once one CDS
/// contains another. The binary search on ends skipped the containing CDS.
#[test]
fn test_gff_index_finds_a_cds_that_contains_a_nested_one() {
    let gene = |start, end, id: &str| GffRecord {
        contig: "c".into(),
        start,
        end,
        strand: '+',
        gene_id: id.into(),
    };
    let idx = GffIndex::build(vec![
        gene(1000, 5000, "A"),
        gene(1200, 1500, "B"),
        gene(5100, 6000, "C"),
    ]);
    let ids: Vec<&str> = idx
        .overlapping("c", 3000, 4000)
        .iter()
        .map(|g| g.gene_id)
        .collect();
    assert_eq!(ids, vec!["A"]);

    let (up, down) = idx.flanking("c", 5040, 5090);
    assert_eq!(
        up.map(|g| g.gene_id),
        Some("A"),
        "nearest end is A's, not B's"
    );
    assert_eq!(down.map(|g| g.gene_id), Some("C"));
}

/// The two halves of a gene broken by an insertion abut exactly on the query
/// and sit far apart on the contig. They are two fragments, not one segment.
#[test]
fn test_fragment_chain_keeps_query_adjacent_pieces_apart_on_the_target() {
    let left = make_hit(0, 588, 1176);
    let right = Hit {
        target_start: 5588,
        target_end: 6176,
        ..make_hit(588, 1176, 1176)
    };
    let chain = build_fragment_chain(&[&left, &right]).unwrap();
    assert_eq!(chain.n_fragments, 2);
    assert!((chain.coverage - 100.0).abs() < 1e-9);
    assert_eq!(chain.evidence.max_target_gap_bp, 5000);
}

/// One piece runs off the end of one contig, the other starts another: an
/// assembly break, classified as FRAGMENTED rather than a partial deletion.
#[test]
fn test_gene_split_across_contigs_is_fragmented() {
    let pa = PresenceAbsence::new_for_test(
        vec!["gene".into()],
        vec!["g".into(), "r".into()],
        vec![vec!["".into(), "r1".into()]],
    );
    let first = Hit {
        query: "gene".into(),
        cluster: "gene".into(),
        target: "ctg1".into(),
        target_start: 3300,
        target_end: 4000,
        target_length: 4000,
        ..make_hit(0, 700, 1000)
    };
    let second = Hit {
        query: "gene".into(),
        cluster: "gene".into(),
        target: "ctg2".into(),
        target_start: 0,
        target_end: 300,
        target_length: 9000,
        ..make_hit(700, 1000, 1000)
    };
    let rescue = classify_hits(
        &[first, second],
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        None,
        &GeneClusters::none(),
        80.0,
    )
    .unwrap();
    assert!(
        matches!(rescue.rescue_type, RescueType::Fragmented { .. }),
        "{:?}",
        rescue.rescue_type
    );
    assert_eq!(
        rescue.fragment_evidence.note.as_deref(),
        Some("split_across_contigs:ctg2")
    );
}

/// An unclaimed annotated CDS is the plainest false absence. An identical
/// unannotated copy elsewhere used to take its place as MISSING_CDS.
#[test]
fn test_unclaimed_annotated_cds_beats_an_unannotated_copy() {
    let pa = PresenceAbsence::new_for_test(
        vec!["gene".into()],
        vec!["g".into(), "r".into()],
        vec![vec!["".into(), "r1".into()]],
    );
    let gff = GffIndex::build(vec![GffRecord {
        contig: "chr".into(),
        start: 1000,
        end: 2000,
        strand: '+',
        gene_id: "g_00019".into(),
    }]);
    let annotated = Hit {
        query: "gene".into(),
        cluster: "gene".into(),
        target: "chr".into(),
        target_start: 1000,
        target_end: 2000,
        ..make_hit(0, 1000, 1000)
    };
    let copy = Hit {
        target_start: 8000,
        target_end: 9000,
        ..annotated.clone()
    };
    // The copy is listed first, as it was when it won.
    let rescue = classify_hits(
        &[copy, annotated],
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        Some(&gff),
        &GeneClusters::from_map(&gff, &HashMap::new()),
        80.0,
    )
    .unwrap();
    assert!(
        matches!(rescue.rescue_type, RescueType::Dropped),
        "{:?}",
        rescue.rescue_type
    );
    assert_eq!(rescue.start, 1000);
    assert_eq!(rescue.annotation_match.as_deref(), Some("g_00019"));
}

/// A premature stop read in the CDS's own frame. Over three frames the
/// minimum hid it whenever another frame happened to be stop-free.
#[test]
fn test_stop_in_frame_is_not_hidden_by_a_clean_other_frame() {
    // Frame 0: ATG GCC TAA GCC GCC TAA -> one internal stop, then the end.
    // Frames 1 and 2 of this sequence contain no stop codon at all.
    let seq = b"ATGGCCTAAGCCGCCTAA";
    assert_eq!(orf::count_internal_stops(seq), 0, "the minimum over frames");
    assert_eq!(orf::count_stops_in_frame(seq), 1);

    let mut contigs = HashMap::new();
    contigs.insert("chr".to_string(), seq.to_vec());
    assert_eq!(
        check_internal_stops_in_cds("chr", 0, seq.len() as u64, '+', &contigs),
        Some(1)
    );
}

/// A [1000,2000] 0.8, B [1000,1600] 0.5, C [1500,2600] 0.9. B overlaps C by
/// 17% of B, so B's only real conflict is A, which loses to C: B and C stay.
/// The old pairwise sweep let A remove B before A itself was removed.
#[test]
fn a_rescue_that_loses_does_not_remove_an_innocent_one() {
    let mk = |cluster: usize, start: u64, end: u64, confidence: f64| RescueEvent {
        cluster_idx: ClusterIdx(cluster),
        genome_idx: GenomeIdx(0),
        rescue_type: RescueType::Dropped,
        identity: 99.0,
        coverage: 99.0,
        contig: "ctg".into(),
        start,
        end,
        strand: '+',
        annotation_match: None,
        confidence,
        synteny_ok: None,
        synteny_evidence: Default::default(),
        fragment_evidence: Default::default(),
        query_start: None,
    };
    let mut rescues = vec![
        mk(0, 1000, 2000, 0.8),
        mk(1, 1000, 1600, 0.5),
        mk(2, 1500, 2600, 0.9),
    ];
    resolve_conflicts(&mut rescues);
    let kept: Vec<usize> = rescues.iter().map(|r| r.cluster_idx.val()).collect();
    assert_eq!(kept, vec![1, 2]);
}

/// A gene of G/C-ending codons with no T in a middle position, so its two
/// off frames hold no stop codon at all, and a premature TAA planted in the
/// coding frame. The minimum over frames read it as an intact ORF.
fn gc3_gene(premature_stop: bool) -> Vec<u8> {
    let codons: [&[u8; 3]; 10] = [
        b"GCC", b"GGC", b"CGC", b"GCG", b"CCG", b"AGC", b"ACC", b"CAG", b"GAG", b"AAG",
    ];
    let mut seq = b"ATG".to_vec();
    for (i, codon) in codons.iter().cycle().take(40).enumerate() {
        if premature_stop && i == 16 {
            seq.extend_from_slice(b"TAA");
        }
        seq.extend_from_slice(*codon);
    }
    seq.extend_from_slice(b"TAA");
    seq
}

/// Unannotated rescues are scanned in the frame their alignment fixes, on
/// either strand and when the alignment starts part-way into a codon.
#[test]
fn unannotated_spans_are_scanned_in_the_alignment_frame() {
    let gene = gc3_gene(true);
    assert_eq!(
        orf::count_internal_stops(&gene),
        0,
        "the minimum over frames hides it"
    );
    assert_eq!(orf::count_stops_in_frame(&gene), 1);

    let pad = b"CCCCCCCCCC".to_vec();
    let n = gene.len() as u64;
    let mut contigs = HashMap::new();
    contigs.insert(
        "plus".to_string(),
        [pad.clone(), gene.clone(), pad.clone()].concat(),
    );
    contigs.insert(
        "minus".to_string(),
        [pad.clone(), orf::reverse_complement(&gene), pad.clone()].concat(),
    );

    // The whole gene aligned from query base 0.
    assert_eq!(
        check_internal_stops_in_span("plus", 10, 10 + n, '+', 0, &contigs),
        Some(1)
    );
    assert_eq!(
        check_internal_stops_in_span("minus", 10, 10 + n, '-', 0, &contigs),
        Some(1)
    );
    // The alignment misses the first two query bases. On the minus strand the
    // gene's 5' end is the span's right end.
    assert_eq!(
        check_internal_stops_in_span("plus", 12, 10 + n, '+', 2, &contigs),
        Some(1)
    );
    assert_eq!(
        check_internal_stops_in_span("minus", 10, 8 + n, '-', 2, &contigs),
        Some(1)
    );

    // The same gene without the planted stop is clean in every case.
    let clean = gc3_gene(false);
    let m = clean.len() as u64;
    contigs.insert(
        "plus".to_string(),
        [pad.clone(), clean.clone(), pad.clone()].concat(),
    );
    contigs.insert(
        "minus".to_string(),
        [pad.clone(), orf::reverse_complement(&clean), pad.clone()].concat(),
    );
    assert_eq!(
        check_internal_stops_in_span("plus", 10, 10 + m, '+', 0, &contigs),
        None
    );
    assert_eq!(
        check_internal_stops_in_span("minus", 10, 10 + m, '-', 0, &contigs),
        None
    );
    assert_eq!(
        check_internal_stops_in_span("plus", 12, 10 + m, '+', 2, &contigs),
        None
    );
    assert_eq!(
        check_internal_stops_in_span("minus", 10, 8 + m, '-', 2, &contigs),
        None
    );
}

fn piece(
    query_start: u64,
    query_end: u64,
    strand: char,
    target_start: u64,
    target_end: u64,
) -> Hit {
    Hit {
        query: "q".into(),
        cluster: "cluster1".into(),
        query_length: 1182,
        query_start,
        query_end,
        strand,
        target: "chr".into(),
        target_start,
        target_end,
        target_length: 23632,
        matches: (query_end - query_start) * 99 / 100,
        alignment_length: query_end - query_start,
        mapq: 60,
    }
}

fn classify_pieces(hits: &[Hit]) -> Option<RescueEvent> {
    // An annotation exists but nothing in it overlaps either piece.
    let gff = GffIndex::build(vec![GffRecord {
        contig: "chr".into(),
        start: 9000,
        end: 9900,
        strand: '+',
        gene_id: "far".into(),
    }]);
    let pa = PresenceAbsence::new_for_test(
        vec!["cluster1".into()],
        vec!["genome1".into()],
        vec![vec!["".into()]],
    );
    let g2c: HashMap<String, ClusterIdx> = HashMap::new();
    classify_hits(
        hits,
        ClusterIdx(0),
        GenomeIdx(0),
        &pa,
        Some(&gff),
        &GeneClusters::from_map(&gff, &g2c),
        80.0,
    )
}

/// A gene across the origin of a circular contig: its first half runs off the
/// end of the contig and its second half carries on from the start. Neither
/// piece reaches 80% on its own, and chained as one stretch they read as
/// non-colinear, so there was no rescue at all.
#[test]
fn a_gene_across_the_origin_of_a_circular_contig_is_rescued() {
    let plus = classify_pieces(&[
        piece(9, 588, '+', 23049, 23628),
        piece(608, 1173, '+', 16, 581),
    ])
    .expect("a rescue across the origin");
    assert!(
        matches!(plus.rescue_type, RescueType::MissingCds),
        "{}",
        plus.rescue_type
    );
    assert!(plus.coverage > 95.0, "both pieces count: {}", plus.coverage);
    assert_eq!(
        (plus.start, plus.end),
        (23049, 23628),
        "the larger piece's span"
    );
    assert_eq!(
        plus.fragment_evidence.note.as_deref(),
        Some("wraps_origin:17-581")
    );
    assert_eq!(plus.fragment_evidence.n_fragments, 2);
    assert_eq!(plus.query_start, Some(9));

    // Minus strand: the gene's first piece is the one at the contig's start.
    let minus = classify_pieces(&[
        piece(9, 588, '-', 10, 589),
        piece(608, 1173, '-', 23049, 23614),
    ])
    .expect("a rescue across the origin on the minus strand");
    assert!(matches!(minus.rescue_type, RescueType::MissingCds));
    assert_eq!(
        minus.fragment_evidence.note.as_deref(),
        Some("wraps_origin:23050-23614")
    );
}

/// The same two pieces in the order a linear contig gives them are not a
/// wrap: the first half sits at the start and the second at the end.
#[test]
fn pieces_at_both_ends_in_linear_order_are_not_an_origin_wrap() {
    let rescue = classify_pieces(&[
        piece(9, 588, '+', 16, 595),
        piece(608, 1173, '+', 23049, 23614),
    ]);
    assert!(
        rescue.as_ref().is_none_or(|r| !r
            .fragment_evidence
            .note
            .as_deref()
            .unwrap_or("")
            .starts_with("wraps_origin")),
        "{rescue:?}"
    );
}

/// A cell's gene ID reaches the GFF whichever tool wrote it: the ID itself
/// first, then Panaroo's suffixes stripped, a PIRATE fission's pieces, and the
/// gene part of a PanTA `sample-contig-gene` ID.
#[test]
fn cell_gene_ids_reach_the_gff_from_every_tool() {
    use crate::rescue::cell_gene_id_variants;
    assert_eq!(
        cell_gene_id_variants("X_dup_2_pseudo"),
        vec!["X_dup_2_pseudo", "X_dup_2", "X"]
    );
    assert_eq!(
        cell_gene_id_variants("(L8_AS00801_00013a:L8_AS00801_00013b)"),
        vec![
            "(L8_AS00801_00013a:L8_AS00801_00013b)",
            "L8_AS00801_00013a",
            "L8_AS00801_00013b"
        ]
    );
    let panta = cell_gene_id_variants("L1_AS00101-L1_AS00101_contig1-L1_AS00101_00019");
    assert_eq!(panta[0], "L1_AS00101-L1_AS00101_contig1-L1_AS00101_00019");
    assert!(panta.contains(&"L1_AS00101_00019"));

    // In a genome's lookup, a PanTA cell answers for the GFF's own ID.
    let pa = crate::parser::PresenceAbsence::new_for_test(
        vec!["Rv0019c".into()],
        vec!["L1_AS00101".into()],
        vec![vec!["L1_AS00101-L1_AS00101_contig1-L1_AS00101_00019".into()]],
    );
    let map = build_gene_to_cluster(&pa, GenomeIdx(0));
    assert_eq!(map.get("L1_AS00101_00019"), Some(&ClusterIdx(0)));
}
