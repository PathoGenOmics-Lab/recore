#!/usr/bin/env python3
"""Compare two Recore output directories.

The script compares matrix/state calls cell-by-cell and evolutionary mechanism
calls cluster-by-cluster. It is designed for old-vs-new regression review.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path


METADATA_COLS = {
    "Gene",
    "Non-unique Gene name",
    "Annotation",
    "No. isolates",
    "No. sequences",
    "Avg sequences per isolate",
    "Genome Fragment",
    "Order within Fragment",
    "Accessory Fragment",
    "Accessory Order with Fragment",
    "Accessory Order",
    "QC",
    "Min group size nuc",
    "Max group size nuc",
    "Avg group size nuc",
}


def matrix_path(run_dir: Path) -> Path:
    states = run_dir / "gene_presence_absence_recore_states.csv"
    if states.exists():
        return states
    matrix = run_dir / "gene_presence_absence_recore.csv"
    if matrix.exists():
        return matrix
    raise SystemExit(f"ERROR: no Recore matrix found in {run_dir}")


def read_states(run_dir: Path) -> tuple[list[str], dict[tuple[str, str], str]]:
    path = matrix_path(run_dir)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        headers = reader.fieldnames or []
        genomes = [h for h in headers if h not in METADATA_COLS]
        states = {}
        is_state_matrix = path.name.endswith("_states.csv")
        for row in reader:
            cluster = row["Gene"]
            for genome in genomes:
                raw = row.get(genome, "").strip()
                if is_state_matrix:
                    state = raw or "ABSENT"
                else:
                    state = "PRESENT" if raw else "ABSENT"
                states[(cluster, genome)] = state
    return genomes, states


def read_evolution(run_dir: Path) -> dict[str, dict[str, str]]:
    path = run_dir / "evolutionary_events.tsv"
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["gene_cluster"]: row for row in reader}


def is_absent(state: str) -> bool:
    """No presence and no evidence: a gap searched and found empty, or one that
    was never searched. A NOT_SEARCHED_* state is an absence, and counting it as
    a call made every unsearched cell a "new correction" against older runs,
    which wrote ABSENT for both."""
    return state in {"", "ABSENT"} or state.startswith("NOT_SEARCHED")


def classify_change(old_state: str, new_state: str) -> str:
    if old_state == "MISSING_CELL":
        return "new_only_cell"
    if new_state == "MISSING_CELL":
        return "old_only_cell"
    if old_state == new_state:
        return "same"
    old_present = not is_absent(old_state)
    new_present = not is_absent(new_state)
    if not old_present and new_present:
        return "new_correction"
    if old_present and not new_present:
        return "lost_correction"
    if not old_present and not new_present:
        # Absent in both runs, but searched in one and not the other, or left
        # out of the search for a different reason: the search changed, the
        # matrix did not.
        return "search_coverage_changed"
    return "state_changed"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old", required=True, type=Path)
    parser.add_argument("--new", required=True, type=Path)
    parser.add_argument("-o", "--outdir", required=True, type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    old_genomes, old_states = read_states(args.old)
    new_genomes, new_states = read_states(args.new)
    common = set(old_states) & set(new_states)
    all_cells = sorted(set(old_states) | set(new_states))

    cell_counts = Counter()
    cluster_counts: dict[str, Counter] = defaultdict(Counter)
    cell_path = args.outdir / "run_comparison_cells.tsv"
    with cell_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_cluster", "genome", "old_state", "new_state", "change_class"])
        for cluster, genome in all_cells:
            old_state = old_states.get((cluster, genome), "MISSING_CELL")
            new_state = new_states.get((cluster, genome), "MISSING_CELL")
            change = classify_change(old_state, new_state)
            if change != "same":
                writer.writerow([cluster, genome, old_state, new_state, change])
            cell_counts[change] += 1
            cluster_counts[cluster][change] += 1

    cluster_path = args.outdir / "run_comparison_clusters.tsv"
    with cluster_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "gene_cluster",
            "same",
            "new_correction",
            "lost_correction",
            "state_changed",
            "old_only_cell",
            "new_only_cell",
            "changed_cells",
            "search_coverage_changed",
        ])
        for cluster in sorted(cluster_counts):
            counts = cluster_counts[cluster]
            changed = (
                counts["new_correction"]
                + counts["lost_correction"]
                + counts["state_changed"]
                + counts["old_only_cell"]
                + counts["new_only_cell"]
            )
            if changed or counts["search_coverage_changed"]:
                writer.writerow([
                    cluster,
                    counts["same"],
                    counts["new_correction"],
                    counts["lost_correction"],
                    counts["state_changed"],
                    counts["old_only_cell"],
                    counts["new_only_cell"],
                    changed,
                    counts["search_coverage_changed"],
                ])

    old_evo = read_evolution(args.old)
    new_evo = read_evolution(args.new)
    evo_path = args.outdir / "run_comparison_evolution.tsv"
    with evo_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "gene_cluster",
            "old_mechanism",
            "new_mechanism",
            "old_priority",
            "new_priority",
            "old_signal_score",
            "new_signal_score",
        ])
        for cluster in sorted(set(old_evo) & set(new_evo)):
            old = old_evo[cluster]
            new = new_evo[cluster]
            fields = (
                old.get("molecular_mechanism", ""),
                new.get("molecular_mechanism", ""),
                old.get("highlight_priority", ""),
                new.get("highlight_priority", ""),
                old.get("evolutionary_signal_score", ""),
                new.get("evolutionary_signal_score", ""),
            )
            if fields[0] != fields[1] or fields[2] != fields[3] or fields[4] != fields[5]:
                writer.writerow([cluster, *fields])

    summary_path = args.outdir / "run_comparison_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["old_dir", str(args.old)])
        writer.writerow(["new_dir", str(args.new)])
        writer.writerow(["old_genomes", len(old_genomes)])
        writer.writerow(["new_genomes", len(new_genomes)])
        writer.writerow(["old_only_genomes", ";".join(sorted(set(old_genomes) - set(new_genomes)))])
        writer.writerow(["new_only_genomes", ";".join(sorted(set(new_genomes) - set(old_genomes)))])
        writer.writerow(["old_cells", len(old_states)])
        writer.writerow(["new_cells", len(new_states)])
        writer.writerow(["common_cells", len(common)])
        writer.writerow(["union_cells", len(all_cells)])
        for key in [
            "same",
            "new_correction",
            "lost_correction",
            "state_changed",
            "old_only_cell",
            "new_only_cell",
            "search_coverage_changed",
        ]:
            writer.writerow([key, cell_counts[key]])

    print(f"Wrote {cell_path}")
    print(f"Wrote {cluster_path}")
    print(f"Wrote {evo_path}")
    print(f"Wrote {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
