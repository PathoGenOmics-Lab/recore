"""Tests for Panaroo output parsing."""

import pandas as pd
import pytest
from recore.parser import find_gaps, parse_pan_reference, parse_presence_absence


@pytest.fixture
def sample_pa_csv(tmp_path):
    """Create a minimal Panaroo gene_presence_absence.csv."""
    csv_path = tmp_path / "gene_presence_absence.csv"
    csv_path.write_text(
        "Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,"
        "Avg sequences per isolate,Genome Fragment,Order within Fragment,"
        "Accessory Order,QC,Min group size nuc,Max group size nuc,"
        "Avg group size nuc,genome1,genome2,genome3\n"
        "geneA,geneA,product A,3,3,1.0,0,0,0,,100,100,100,g1_001,g2_001,g3_001\n"
        "geneB,geneB,product B,2,2,1.0,0,0,0,,200,200,200,g1_002,,g3_002\n"
        "geneC,geneC,product C,1,1,1.0,0,0,0,,150,150,150,g1_003,,\n"
    )
    return csv_path


@pytest.fixture
def sample_pan_ref(tmp_path):
    """Create a minimal pan_genome_reference.fa."""
    fasta_path = tmp_path / "pan_genome_reference.fa"
    fasta_path.write_text(
        ">geneA~~~g1_001 product A\nATGCATGCATGC\n"
        ">geneB~~~g1_002 product B\nATGCATGCATGCATGCATGC\n"
        ">geneC~~~g1_003 product C\nATGCATGCATGCATGC\n"
    )
    return fasta_path


def test_parse_presence_absence(sample_pa_csv):
    pa = parse_presence_absence(sample_pa_csv)
    assert pa.shape == (3, 3)
    assert list(pa.columns) == ["genome1", "genome2", "genome3"]
    assert list(pa.index) == ["geneA", "geneB", "geneC"]
    assert pa.at["geneA", "genome1"] == "g1_001"
    assert pd.isna(pa.at["geneB", "genome2"])


def test_find_gaps(sample_pa_csv):
    pa = parse_presence_absence(sample_pa_csv)
    gaps = find_gaps(pa)
    # geneB missing in genome2, geneC missing in genome2 and genome3
    assert len(gaps) == 3
    gene_genome_pairs = {(g, gn) for g, gn in gaps}
    assert ("geneB", "genome2") in gene_genome_pairs
    assert ("geneC", "genome2") in gene_genome_pairs
    assert ("geneC", "genome3") in gene_genome_pairs


def test_parse_pan_reference(sample_pan_ref):
    refs = parse_pan_reference(sample_pan_ref)
    assert len(refs) == 3
    assert "geneA" in refs
    assert refs["geneA"] == "ATGCATGCATGC"
    assert refs["geneB"] == "ATGCATGCATGCATGCATGC"
