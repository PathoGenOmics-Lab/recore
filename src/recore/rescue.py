"""Core rescue logic: classify and rescue falsely absent genes."""

import logging
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from recore.search import SearchHit, search_genome

log = logging.getLogger("recore")


@dataclass
class RescueEvent:
    """A single gene rescue."""

    gene_cluster: str
    genome: str
    rescue_type: str  # SPLIT, DROPPED, OVERLAP
    identity: float
    coverage: float
    contig: str
    start: int
    end: int
    strand: str
    original_cluster: str | None = None  # if SPLIT: which cluster had this gene
    annotation_match: str | None = None  # if --annotations: overlapping CDS


def classify_rescue(
    hit: SearchHit,
    genome: str,
    gene_cluster: str,
    pa: pd.DataFrame,
    gene_coords: dict[str, list[tuple[int, int, str]]] | None = None,
) -> RescueEvent:
    """Classify a rescue event.

    Args:
        hit: minimap2 hit for the absent gene
        genome: genome name
        gene_cluster: cluster where the gene was absent
        pa: presence/absence matrix
        gene_coords: optional {genome: [(start, end, gene_id)]} from GFF
    """
    start = hit.target_start
    end = hit.target_end

    rescue_type = "DROPPED"
    original_cluster = None
    annotation_match = None

    # Check annotation overlap if available
    if gene_coords and genome in gene_coords:
        for gstart, gend, gene_id in gene_coords[genome]:
            if start < gend and end > gstart:
                overlap = min(end, gend) - max(start, gstart)
                hit_len = end - start
                if hit_len > 0 and overlap / hit_len > 0.5:
                    annotation_match = gene_id
                    break

    return RescueEvent(
        gene_cluster=gene_cluster,
        genome=genome,
        rescue_type=rescue_type,
        identity=hit.identity,
        coverage=hit.query_coverage,
        contig=hit.target,
        start=start,
        end=end,
        strand=hit.strand,
        original_cluster=original_cluster,
        annotation_match=annotation_match,
    )


def rescue_one_genome(
    genome: str,
    genome_fasta: Path,
    gaps: list[str],
    pan_ref: dict[str, str],
    pa: pd.DataFrame,
    gene_coords: dict[str, list[tuple[int, int, str]]] | None,
    threads: int,
    min_identity: float,
    min_coverage: float,
) -> list[RescueEvent]:
    """Run rescue for all gaps in a single genome.

    Args:
        genome: genome name
        genome_fasta: path to genome FASTA
        gaps: list of gene cluster names absent in this genome
        pan_ref: {cluster: representative_sequence}
        pa: full presence/absence matrix
        gene_coords: annotation coords if available
        threads: minimap2 threads
        min_identity: minimum % identity
        min_coverage: minimum % query coverage

    Returns:
        List of RescueEvent for genes successfully rescued.
    """
    # Build query set
    queries = {}
    skipped = 0
    for gene in gaps:
        if gene in pan_ref:
            queries[gene] = pan_ref[gene]
        else:
            skipped += 1

    if skipped > 0:
        log.debug(f"{genome}: skipped {skipped} clusters with no representative sequence")

    if not queries:
        return []

    log.debug(f"{genome}: searching {len(queries)} absent genes")

    # minimap2 search
    hits = search_genome(
        queries=queries,
        genome_fasta=genome_fasta,
        threads=threads,
        min_identity=min_identity,
        min_coverage=min_coverage,
    )

    # Classify
    rescues = []
    for hit in hits:
        rescue = classify_rescue(
            hit=hit,
            genome=genome,
            gene_cluster=hit.query,
            pa=pa,
            gene_coords=gene_coords,
        )
        rescues.append(rescue)

    return rescues
