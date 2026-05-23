#!/usr/bin/env python3
"""Dataset-specific pseudo-absence benchmark for Recore.

This benchmark masks real present cells from a user-provided pangenome matrix,
runs Recore on the masked matrix, and measures how often the masked positives
are recovered. It also samples unchanged absent cells as a negative-control
proxy; those cells are useful for review, but they are not a biological gold
standard because some original absences may be genuine false absences.
"""

from __future__ import annotations

import argparse
import csv
import os
import random
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORKDIR = Path("/private/tmp/recore_self_benchmark")
DEFAULT_METRICS = ROOT / "benchmarks" / "self_benchmark_metrics.tsv"

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


def read_matrix(path: Path) -> tuple[list[str], list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        headers = reader.fieldnames or []
        genomes = [h for h in headers if h not in METADATA_COLS]
        rows = list(reader)
    return headers, genomes, rows


def update_metadata(row: dict[str, str], genomes: list[str]) -> None:
    present = sum(1 for genome in genomes if row.get(genome, "").strip())
    if "No. isolates" in row:
        row["No. isolates"] = str(present)
    if "No. sequences" in row:
        row["No. sequences"] = str(present)
    if "Avg sequences per isolate" in row:
        row["Avg sequences per isolate"] = f"{present:.3f}" if present else "0"


def write_matrix(path: Path, headers: list[str], genomes: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            update_metadata(row, genomes)
            writer.writerow(row)


def ensure_binary(skip_build: bool) -> Path:
    binary = ROOT / "target" / "release" / "recore"
    if not skip_build:
        subprocess.run(["cargo", "build", "--release", "--locked"], cwd=ROOT, check=True)
    if not binary.exists():
        raise SystemExit(f"ERROR: Recore binary not found: {binary}")
    return binary


def select_truth(
    rows: list[dict[str, str]],
    genomes: list[str],
    n_positive: int,
    n_negative: int,
    min_prevalence: float,
    seed: int,
) -> tuple[list[tuple[str, str, str]], list[tuple[str, str, str]]]:
    rng = random.Random(seed)
    positives: list[tuple[str, str, str]] = []
    negatives: list[tuple[str, str, str]] = []
    n_genomes = max(len(genomes), 1)
    for row in rows:
        cluster = row["Gene"]
        present = [g for g in genomes if row.get(g, "").strip()]
        prevalence = len(present) / n_genomes
        if prevalence < min_prevalence:
            continue
        for genome in present:
            positives.append((cluster, genome, row[genome]))
        for genome in genomes:
            if not row.get(genome, "").strip():
                negatives.append((cluster, genome, ""))
    rng.shuffle(positives)
    rng.shuffle(negatives)
    return positives[:n_positive], negatives[:n_negative]


def run_recore(args: argparse.Namespace, masked_matrix: Path, output_dir: Path) -> None:
    binary = Path(args.recore) if args.recore else ensure_binary(args.skip_build)
    cmd = [
        str(binary),
        "-p",
        str(masked_matrix),
        "-r",
        str(args.pan_reference),
        "-g",
        str(args.genomes),
        "--min-prevalence",
        str(args.min_prevalence),
        "--matrix-policy",
        args.matrix_policy,
        "--preset",
        args.preset,
        "-t",
        str(args.threads),
        "-o",
        str(output_dir),
    ]
    if args.annotations:
        cmd.extend(["-a", str(args.annotations), "--synteny"])
    if args.multi_allele:
        cmd.append("--multi-allele")
    if args.nscluster:
        cmd.extend(["--nscluster", "--nsc-engine", args.nsc_engine])
    if args.reference_genome:
        cmd.extend(["--reference-genome", args.reference_genome])
    if args.tree:
        cmd.extend(["--tree", str(args.tree)])
    if args.mobile_element_reference:
        cmd.extend(["--mobile-element-reference", str(args.mobile_element_reference)])

    env = os.environ.copy()
    if args.minimap2:
        env["PATH"] = f"{Path(args.minimap2).parent}{os.pathsep}{env.get('PATH', '')}"
    print("Command:", " ".join(cmd))
    subprocess.run(cmd, cwd=ROOT, env=env, check=True)


def read_corrected_presence(path: Path) -> dict[tuple[str, str], bool]:
    _headers, genomes, rows = read_matrix(path)
    presence: dict[tuple[str, str], bool] = {}
    for row in rows:
        cluster = row["Gene"]
        for genome in genomes:
            presence[(cluster, genome)] = bool(row.get(genome, "").strip())
    return presence


def read_orthology_scores(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {(r["gene_cluster"], r["genome"]): r for r in reader}


def score_bin(score_text: str) -> str:
    try:
        score = float(score_text)
    except (TypeError, ValueError):
        return "missing"
    low = int(score * 10) / 10
    high = min(low + 0.1, 1.0)
    return f"{low:.1f}-{high:.1f}"


def write_truth(path: Path, positives: list[tuple[str, str, str]], negatives: list[tuple[str, str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["truth_class", "gene_cluster", "genome", "original_cell"])
        for cluster, genome, value in positives:
            writer.writerow(["masked_present_positive", cluster, genome, value])
        for cluster, genome, value in negatives:
            writer.writerow(["unchanged_absent_control", cluster, genome, value])


def write_calibration(
    path: Path,
    positives: list[tuple[str, str, str]],
    negatives: list[tuple[str, str, str]],
    corrected: dict[tuple[str, str], bool],
    scores: dict[tuple[str, str], dict[str, str]],
) -> None:
    rows = []
    for label, records in [
        ("masked_present_positive", positives),
        ("unchanged_absent_control", negatives),
    ]:
        for cluster, genome, _ in records:
            key = (cluster, genome)
            score_row = scores.get(key, {})
            recovered = corrected.get(key, False)
            rows.append(
                {
                    "truth_class": label,
                    "gene_cluster": cluster,
                    "genome": genome,
                    "recovered_or_corrected": str(recovered),
                    "orthology_confidence_score": score_row.get("orthology_confidence_score", ""),
                    "orthology_confidence_tier": score_row.get("orthology_confidence_tier", ""),
                    "matrix_action": score_row.get("matrix_action", ""),
                    "risk_flags": score_row.get("risk_flags", ""),
                }
            )

    with path.open("w", newline="") as handle:
        fieldnames = list(rows[0].keys()) if rows else [
            "truth_class",
            "gene_cluster",
            "genome",
            "recovered_or_corrected",
            "orthology_confidence_score",
            "orthology_confidence_tier",
            "matrix_action",
            "risk_flags",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_score_bins(path: Path, calibration_path: Path) -> None:
    bins: dict[tuple[str, str], Counter] = defaultdict(Counter)
    with calibration_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["truth_class"], score_bin(row["orthology_confidence_score"]))
            bins[key]["n"] += 1
            if row["recovered_or_corrected"] == "True":
                bins[key]["recovered"] += 1
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["truth_class", "score_bin", "n", "recovered_or_corrected", "rate"])
        for (truth_class, bin_name), counts in sorted(bins.items()):
            n = counts["n"]
            rec = counts["recovered"]
            writer.writerow([truth_class, bin_name, n, rec, f"{rec / n:.4f}" if n else ""])


def write_metrics(path: Path, row: dict[str, str]) -> None:
    exists = path.exists()
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def prepare_workdir(path: Path, keep_workdir: bool) -> None:
    path = path.resolve()
    protected = {ROOT.resolve(), ROOT.parent.resolve(), Path("/").resolve(), Path.home().resolve()}
    if path.exists() and not keep_workdir:
        if path in protected:
            raise SystemExit(f"ERROR: refusing to clean protected workdir: {path}")
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-p", "--presence-absence", required=True, type=Path)
    parser.add_argument("-r", "--pan-reference", required=True, type=Path)
    parser.add_argument("-g", "--genomes", required=True, type=Path)
    parser.add_argument("-a", "--annotations", type=Path)
    parser.add_argument("-o", "--workdir", default=DEFAULT_WORKDIR, type=Path)
    parser.add_argument("--metrics", default=DEFAULT_METRICS, type=Path)
    parser.add_argument("--recore", type=Path)
    parser.add_argument("--minimap2", type=Path)
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--n-positive", type=int, default=200)
    parser.add_argument("--n-negative", type=int, default=200)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--min-prevalence", type=float, default=0.05)
    parser.add_argument("--matrix-policy", default="balanced", choices=["strict", "balanced", "sensitive"])
    parser.add_argument("--preset", default="asm5")
    parser.add_argument("--multi-allele", action="store_true")
    parser.add_argument("--nscluster", action="store_true")
    parser.add_argument("--nsc-engine", default="minhash")
    parser.add_argument("--reference-genome")
    parser.add_argument("--tree", type=Path)
    parser.add_argument("--mobile-element-reference", type=Path)
    parser.add_argument("--keep-workdir", action="store_true")
    args = parser.parse_args()

    prepare_workdir(args.workdir, args.keep_workdir)
    input_dir = args.workdir / "input"
    output_dir = args.workdir / "output"
    input_dir.mkdir(exist_ok=True)
    output_dir.mkdir(exist_ok=True)

    headers, genome_cols, rows = read_matrix(args.presence_absence)
    positives, negatives = select_truth(
        rows,
        genome_cols,
        args.n_positive,
        args.n_negative,
        args.min_prevalence,
        args.seed,
    )
    if not positives:
        raise SystemExit("ERROR: no positive cells selected; lower --min-prevalence or --n-positive.")
    masked = {(cluster, genome) for cluster, genome, _ in positives}
    for row in rows:
        cluster = row["Gene"]
        for genome in genome_cols:
            if (cluster, genome) in masked:
                row[genome] = ""

    masked_matrix = input_dir / "gene_presence_absence_masked.csv"
    truth_path = args.workdir / "self_benchmark_truth.tsv"
    write_matrix(masked_matrix, headers, genome_cols, rows)
    write_truth(truth_path, positives, negatives)

    started = time.monotonic()
    run_recore(args, masked_matrix, output_dir)
    runtime = time.monotonic() - started

    corrected = read_corrected_presence(output_dir / "gene_presence_absence_recore.csv")
    scores = read_orthology_scores(output_dir / "orthology_evidence.tsv")
    tp = sum(1 for cluster, genome, _ in positives if corrected.get((cluster, genome), False))
    fn = len(positives) - tp
    negative_corrected = sum(1 for cluster, genome, _ in negatives if corrected.get((cluster, genome), False))
    recall = tp / len(positives) if positives else 0.0
    precision_proxy = (
        f"{tp / (tp + negative_corrected):.4f}"
        if negatives and (tp + negative_corrected)
        else "NA"
    )

    calibration_path = args.workdir / "self_benchmark_calibration.tsv"
    bins_path = args.workdir / "self_benchmark_score_bins.tsv"
    write_calibration(calibration_path, positives, negatives, corrected, scores)
    write_score_bins(bins_path, calibration_path)

    row = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "matrix_policy": args.matrix_policy,
        "threads": str(args.threads),
        "masked_positives": str(len(positives)),
        "negative_controls": str(len(negatives)),
        "true_positive_recovered": str(tp),
        "false_negative_not_recovered": str(fn),
        "negative_controls_corrected_proxy": str(negative_corrected),
        "masked_positive_recall": f"{recall:.4f}",
        "precision_proxy": precision_proxy,
        "runtime_seconds": f"{runtime:.2f}",
        "workdir": str(args.workdir),
        "output_dir": str(output_dir),
    }
    write_metrics(args.metrics, row)
    for key, value in row.items():
        print(f"{key}: {value}")
    print(f"Wrote truth: {truth_path}")
    print(f"Wrote calibration: {calibration_path}")
    print(f"Wrote score bins: {bins_path}")
    print(f"Wrote metrics: {args.metrics}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
