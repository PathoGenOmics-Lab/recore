"""Main pipeline: orchestrate parsing, BLAST, rescue, and output."""

import logging
from collections import defaultdict
from pathlib import Path

import pandas as pd

from recore.parser import (
    find_gaps,
    find_genome_files,
    parse_gff_genes,
    parse_pan_reference,
    parse_presence_absence,
)
from recore.rescue import rescue_one_genome

log = logging.getLogger("recore")


def run_recore(args) -> dict:
    """Run the full Recore pipeline.

    Returns dict with summary statistics.
    """
    # ── Parse inputs ──
    pa = parse_presence_absence(args.presence_absence)
    pan_ref = parse_pan_reference(args.pan_reference)
    genome_files = find_genome_files(args.genomes)
    gaps = find_gaps(pa)

    # Validate: genome columns in PA must match available FASTA files
    pa_genomes = set(pa.columns)
    file_genomes = set(genome_files.keys())
    missing_files = pa_genomes - file_genomes
    if missing_files:
        log.warning(
            f"{len(missing_files)} genomes in PA matrix have no FASTA file: "
            f"{', '.join(sorted(missing_files)[:5])}..."
        )

    # Group gaps by genome
    gaps_by_genome: dict[str, list[str]] = defaultdict(list)
    for gene, genome in gaps:
        if genome in genome_files:
            gaps_by_genome[genome].append(gene)

    log.info(
        f"Will search {sum(len(g) for g in gaps_by_genome.values())} gaps "
        f"across {len(gaps_by_genome)} genomes"
    )

    # ── Parse annotations (optional) ──
    genome_gene_coords = None
    if args.annotations:
        genome_gene_coords = {}
        for genome_name in gaps_by_genome:
            gff_path = args.annotations / f"{genome_name}.gff"
            if not gff_path.exists():
                gff_path = args.annotations / f"{genome_name}.gff3"
            if gff_path.exists():
                genes = parse_gff_genes(gff_path)
                genome_gene_coords[genome_name] = [
                    (g["start"], g["end"], g["gene_id"]) for g in genes
                ]
        log.info(f"Loaded annotations for {len(genome_gene_coords)} genomes")

    # ── Run rescue per genome ──
    all_rescues = []
    for i, (genome, gene_gaps) in enumerate(sorted(gaps_by_genome.items()), 1):
        if i % 50 == 0 or i == len(gaps_by_genome):
            log.info(f"Processing genome {i}/{len(gaps_by_genome)}: {genome}")

        rescues = rescue_one_genome(
            genome=genome,
            genome_fasta=genome_files[genome],
            gaps=gene_gaps,
            pan_ref=pan_ref,
            pa=pa,
            gene_coords=genome_gene_coords,
            threads=args.threads,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
        )
        all_rescues.extend(rescues)

    log.info(f"Total rescues: {len(all_rescues)}")

    # ── Classify rescue types ──
    type_counts = defaultdict(int)
    for r in all_rescues:
        type_counts[r.rescue_type] += 1
    for rtype, count in sorted(type_counts.items()):
        log.info(f"  {rtype}: {count}")

    # ── Generate outputs ──
    args.output.mkdir(parents=True, exist_ok=True)

    # 1. Corrected presence/absence matrix
    pa_corrected = pa.copy()
    for rescue in all_rescues:
        pa_corrected.at[rescue.gene_cluster, rescue.genome] = (
            f"RESCUED:{rescue.contig}:{rescue.start}-{rescue.end}"
        )
    pa_corrected.to_csv(args.output / "gene_presence_absence_recore.csv")

    # 2. Rescue report
    rescue_rows = []
    for r in all_rescues:
        rescue_rows.append({
            "gene_cluster": r.gene_cluster,
            "genome": r.genome,
            "rescue_type": r.rescue_type,
            "identity": f"{r.identity:.1f}",
            "coverage": f"{r.coverage:.1f}",
            "contig": r.contig,
            "start": r.start,
            "end": r.end,
            "strand": r.strand,
            "original_cluster": r.original_cluster or "",
            "annotation_match": r.annotation_match or "",
        })
    rescue_df = pd.DataFrame(rescue_rows)
    rescue_df.to_csv(args.output / "rescue_report.tsv", sep="\t", index=False)

    # 3. Summary statistics
    n_clusters = len(pa)
    n_genomes = len(pa.columns)
    n_gaps_total = len(gaps)
    n_rescued = len(all_rescues)
    n_remaining = n_gaps_total - n_rescued

    # Core genome sizes: before and after
    def core_size(matrix):
        return int((matrix.notna() & (matrix != "")).all(axis=1).sum())

    core_before = core_size(pa)
    core_after = core_size(pa_corrected)

    summary_lines = [
        "Recore Summary",
        "=" * 40,
        f"Clusters:           {n_clusters}",
        f"Genomes:            {n_genomes}",
        f"Total gaps:         {n_gaps_total}",
        f"Rescued:            {n_rescued} ({n_rescued/n_gaps_total*100:.1f}%)" if n_gaps_total > 0 else f"Rescued:            {n_rescued}",
        f"Remaining gaps:     {n_remaining}",
        f"Core size (before): {core_before}",
        f"Core size (after):  {core_after} (+{core_after - core_before})",
        "",
        "Rescue types:",
    ]
    for rtype, count in sorted(type_counts.items()):
        summary_lines.append(f"  {rtype}: {count}")

    summary_text = "\n".join(summary_lines) + "\n"
    (args.output / "summary.txt").write_text(summary_text)
    log.info("\n" + summary_text)

    return {
        "rescued": n_rescued,
        "genomes": len(gaps_by_genome),
        "core_before": core_before,
        "core_after": core_after,
        "type_counts": dict(type_counts),
    }
