"""Command-line interface for Recore."""

import argparse
import logging
import sys
from pathlib import Path

from recore import __version__
from recore.pipeline import run_recore


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="recore",
        description="Post-hoc rescue of falsely absent genes in pangenome matrices.",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )

    # Required inputs
    required = parser.add_argument_group("required inputs")
    required.add_argument(
        "--presence-absence",
        required=True,
        type=Path,
        help="Panaroo gene_presence_absence.csv (or Roary-format CSV)",
    )
    required.add_argument(
        "--pan-reference",
        required=True,
        type=Path,
        help="Panaroo pan_genome_reference.fa (representative sequences per cluster)",
    )
    required.add_argument(
        "--genomes",
        required=True,
        type=Path,
        help="Directory of genome FASTA files (one per genome)",
    )

    # Optional inputs
    optional = parser.add_argument_group("optional inputs")
    optional.add_argument(
        "--annotations",
        type=Path,
        default=None,
        help="Directory of GFF files for validation (one per genome, matching names)",
    )

    # BLAST thresholds
    blast = parser.add_argument_group("BLAST thresholds")
    blast.add_argument(
        "--min-identity",
        type=float,
        default=90.0,
        help="Minimum %% identity for a BLAST hit (default: 90.0)",
    )
    blast.add_argument(
        "--min-coverage",
        type=float,
        default=80.0,
        help="Minimum %% query coverage for a BLAST hit (default: 80.0)",
    )

    # Performance
    perf = parser.add_argument_group("performance")
    perf.add_argument(
        "--threads",
        type=int,
        default=4,
        help="BLAST threads (default: 4)",
    )

    # Output
    output = parser.add_argument_group("output")
    output.add_argument(
        "--output",
        type=Path,
        default=Path("recore_output"),
        help="Output directory (default: recore_output)",
    )
    output.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging",
    )

    return parser.parse_args(argv)


def validate_args(args):
    """Validate inputs exist and are accessible."""
    errors = []

    if not args.presence_absence.exists():
        errors.append(f"Presence/absence file not found: {args.presence_absence}")
    if not args.pan_reference.exists():
        errors.append(f"Pan-genome reference not found: {args.pan_reference}")
    if not args.genomes.is_dir():
        errors.append(f"Genomes directory not found: {args.genomes}")
    if args.annotations and not args.annotations.is_dir():
        errors.append(f"Annotations directory not found: {args.annotations}")

    if args.min_identity < 0 or args.min_identity > 100:
        errors.append(f"--min-identity must be 0-100, got {args.min_identity}")
    if args.min_coverage < 0 or args.min_coverage > 100:
        errors.append(f"--min-coverage must be 0-100, got {args.min_coverage}")

    return errors


def main(argv=None):
    args = parse_args(argv)

    # Setup logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="[%(asctime)s] %(levelname)s — %(message)s",
        datefmt="%H:%M:%S",
    )
    log = logging.getLogger("recore")

    # Validate
    errors = validate_args(args)
    if errors:
        for e in errors:
            log.error(e)
        sys.exit(1)

    # Run
    log.info(f"Recore v{__version__}")
    log.info(f"Presence/absence: {args.presence_absence}")
    log.info(f"Pan-genome reference: {args.pan_reference}")
    log.info(f"Genomes: {args.genomes}")
    log.info(f"Thresholds: identity≥{args.min_identity}%, coverage≥{args.min_coverage}%")

    try:
        stats = run_recore(args)
        log.info(f"Done. {stats['rescued']} genes rescued across {stats['genomes']} genomes.")
        log.info(f"Output: {args.output}")
    except Exception:
        log.exception("Recore failed")
        sys.exit(1)
