"""Sequence search via minimap2."""

import logging
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

log = logging.getLogger("recore")


@dataclass
class SearchHit:
    """A single minimap2 hit (PAF format)."""

    query: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target: str  # genome contig
    target_length: int
    target_start: int
    target_end: int
    matches: int
    alignment_length: int
    mapq: int

    @property
    def identity(self) -> float:
        """Alignment identity as percentage."""
        if self.alignment_length == 0:
            return 0.0
        return (self.matches / self.alignment_length) * 100.0

    @property
    def query_coverage(self) -> float:
        """Query coverage as percentage."""
        if self.query_length == 0:
            return 0.0
        return ((self.query_end - self.query_start) / self.query_length) * 100.0


def run_minimap2(
    query_fasta: Path,
    target_fasta: Path,
    threads: int = 4,
    preset: str = "asm20",
) -> list[SearchHit]:
    """Run minimap2 and parse PAF output.

    preset 'asm20' is for cross-species with ≤20% divergence (fine for MTBC).
    For intra-species (>95% identity), 'asm5' would also work.

    Returns list of SearchHit objects.
    """
    cmd = [
        "minimap2",
        "-x", preset,
        "-t", str(threads),
        "--secondary=yes",  # report secondary alignments for split detection
        "-N", "5",          # max secondary alignments per query
        str(target_fasta),
        str(query_fasta),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"minimap2 failed: {result.stderr}")

    hits = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        hits.append(SearchHit(
            query=fields[0],
            query_length=int(fields[1]),
            query_start=int(fields[2]),
            query_end=int(fields[3]),
            strand=fields[4],
            target=fields[5],
            target_length=int(fields[6]),
            target_start=int(fields[7]),
            target_end=int(fields[8]),
            matches=int(fields[9]),
            alignment_length=int(fields[10]),
            mapq=int(fields[11]),
        ))

    return hits


def search_genome(
    queries: dict[str, str],
    genome_fasta: Path,
    threads: int = 4,
    min_identity: float = 90.0,
    min_coverage: float = 80.0,
) -> list[SearchHit]:
    """Search query sequences against a genome with minimap2.

    Args:
        queries: {gene_cluster: nucleotide_sequence}
        genome_fasta: path to genome FASTA
        threads: minimap2 threads
        min_identity: minimum % identity to keep a hit
        min_coverage: minimum % query coverage to keep a hit

    Returns:
        List of SearchHit objects passing thresholds (best hit per query).
    """
    if not queries:
        return []

    with tempfile.TemporaryDirectory(prefix="recore_") as tmpdir:
        query_fasta = Path(tmpdir) / "queries.fasta"
        records = [
            SeqRecord(Seq(seq), id=name, description="")
            for name, seq in queries.items()
        ]
        SeqIO.write(records, str(query_fasta), "fasta")

        hits = run_minimap2(query_fasta, genome_fasta, threads=threads)

    # Filter by thresholds and keep best hit per query (highest identity × coverage)
    best_per_query: dict[str, SearchHit] = {}
    for hit in hits:
        if hit.identity < min_identity:
            continue
        if hit.query_coverage < min_coverage:
            continue
        score = hit.identity * hit.query_coverage
        if hit.query not in best_per_query:
            best_per_query[hit.query] = hit
        else:
            prev_score = best_per_query[hit.query].identity * best_per_query[hit.query].query_coverage
            if score > prev_score:
                best_per_query[hit.query] = hit

    return list(best_per_query.values())
