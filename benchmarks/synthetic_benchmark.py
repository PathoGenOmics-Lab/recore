#!/usr/bin/env python3
"""Synthetic precision/recall benchmark for Recore matrix policies.

The benchmark creates a small pangenome matrix with known false absences:

- intact hidden genes that should be corrected in the binary matrix
- genuinely absent genes that should remain absent
- partial pseudogenes that should be audited but not treated as clean binary
  presence under the balanced policy

It then runs Recore with strict, balanced, and sensitive matrix policies,
computes precision/recall/F1 from the corrected matrix, writes a metrics TSV,
and renders a dependency-free SVG plot for the README.
"""

from __future__ import annotations

import argparse
import csv
import html
import os
import random
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORKDIR = Path(tempfile.gettempdir()) / "recore_synthetic_benchmark"
DEFAULT_METRICS = ROOT / "benchmarks" / "synthetic_benchmark_metrics.tsv"
DEFAULT_PLOT = ROOT / "assets" / "synthetic_benchmark_precision_recall.svg"
POLICIES = ("strict", "balanced", "sensitive")
GENOMES = ("s1", "s2", "s3", "s4")

METADATA_COLS = [
    "Gene",
    "Non-unique Gene name",
    "Annotation",
    "No. isolates",
    "No. sequences",
    "Avg sequences per isolate",
    "Genome Fragment",
    "Order within Fragment",
    "Accessory Order",
    "QC",
    "Min group size nuc",
    "Max group size nuc",
    "Avg group size nuc",
]


PRESENT_IN_MATRIX = {
    "gCore": set(GENOMES),
    "gClean1": {"s1"},
    "gClean2": {"s1"},
    "gClean3": {"s1"},
    "gPseudo": {"s1"},
    "gAbsent": {"s1"},
}

TRUE_BINARY_POSITIVES = {
    ("gClean1", "s2"),
    ("gClean1", "s3"),
    ("gClean2", "s2"),
    ("gClean2", "s4"),
    ("gClean3", "s3"),
}

PARTIAL_PSEUDOGENES = {
    ("gPseudo", "s2"),
    ("gPseudo", "s3"),
}


CODONS = [
    "GCT",
    "GCC",
    "GCA",
    "GCG",
    "CGT",
    "CGC",
    "CGA",
    "CGG",
    "AAT",
    "AAC",
    "GAT",
    "GAC",
    "TGT",
    "TGC",
    "CAA",
    "CAG",
    "GAA",
    "GAG",
    "GGT",
    "GGC",
    "GGA",
    "GGG",
    "CAT",
    "CAC",
    "ATT",
    "ATC",
    "ATA",
    "CTT",
    "CTC",
    "CTA",
    "CTG",
    "AAA",
    "AAG",
    "TTT",
    "TTC",
    "CCT",
    "CCC",
    "CCA",
    "CCG",
    "TCT",
    "TCC",
    "TCA",
    "TCG",
    "ACT",
    "ACC",
    "ACA",
    "ACG",
    "TAT",
    "TAC",
    "GTT",
    "GTC",
    "GTA",
    "GTG",
]


def coding_sequence(seed: str, codons: int = 78) -> str:
    rng = random.Random(seed)
    body = "".join(rng.choice(CODONS) for _ in range(codons))
    return "ATG" + body + "TAA"


def wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def spacer(seed: str, length: int = 80) -> str:
    rng = random.Random(f"spacer-{seed}")
    return "".join(rng.choice("ACGT") for _ in range(length))


def metric_fmt(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{value:.3f}"


def require_tool(path_or_name: str, label: str) -> str:
    path = Path(path_or_name)
    if path.exists():
        return str(path)
    found = shutil.which(path_or_name)
    if found:
        return found
    raise SystemExit(f"ERROR: {label} not found: {path_or_name}")


def find_minimap2(explicit: str | None) -> str:
    if explicit:
        return require_tool(explicit, "minimap2")
    found = shutil.which("minimap2")
    if found:
        return found
    candidates = [
        Path.home() / "miniconda3" / "envs" / "bioconda" / "bin" / "minimap2",
        Path.home() / "miniconda3" / "bin" / "minimap2",
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    raise SystemExit(
        "ERROR: minimap2 not found in PATH. Pass --minimap2 /path/to/minimap2."
    )


def ensure_recore_binary(recore: str, skip_build: bool) -> str:
    default_binary = ROOT / "target" / "debug" / "recore"
    path = Path(recore)
    if path.exists() and (skip_build or path != default_binary):
        return str(path)
    if path != default_binary and not path.exists():
        raise SystemExit(f"ERROR: Recore binary not found: {path}")
    if not skip_build:
        subprocess.run(["cargo", "build"], cwd=ROOT, check=True)
    if not default_binary.exists():
        raise SystemExit(f"ERROR: Recore binary not found: {default_binary}")
    return str(default_binary)


def write_inputs(input_dir: Path) -> tuple[set[tuple[str, str]], set[tuple[str, str]]]:
    input_dir.mkdir(parents=True, exist_ok=True)
    genomes_dir = input_dir / "genomes"
    genomes_dir.mkdir(parents=True, exist_ok=True)

    sequences = {cluster: coding_sequence(cluster) for cluster in PRESENT_IN_MATRIX}

    with (input_dir / "pan_ref.fa").open("w") as handle:
        for cluster, seq in sequences.items():
            handle.write(f">{cluster}\n{wrap_fasta(seq)}\n")

    header = METADATA_COLS + list(GENOMES)
    with (input_dir / "pa.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for cluster, present in PRESENT_IN_MATRIX.items():
            length = len(sequences[cluster])
            meta = [
                cluster,
                cluster,
                f"synthetic {cluster}",
                str(len(present)),
                str(len(present)),
                "1.0",
                "0",
                "0",
                "0",
                "",
                str(length),
                str(length),
                str(length),
            ]
            cells = [
                f"{cluster}_{genome}" if genome in present else "" for genome in GENOMES
            ]
            writer.writerow(meta + cells)

    for genome in GENOMES:
        chunks = []
        for cluster, seq in sequences.items():
            if genome in PRESENT_IN_MATRIX[cluster]:
                chunks.extend([spacer(f"{genome}-{cluster}-left"), seq])
            elif (cluster, genome) in TRUE_BINARY_POSITIVES:
                chunks.extend([spacer(f"{genome}-{cluster}-hidden"), seq])
            elif (cluster, genome) in PARTIAL_PSEUDOGENES:
                partial_len = int(len(seq) * 0.70)
                partial_len -= partial_len % 3
                chunks.extend([spacer(f"{genome}-{cluster}-pseudo"), seq[:partial_len]])
        genome_seq = "".join(chunks) + spacer(f"{genome}-tail", 120)
        with (genomes_dir / f"{genome}.fasta").open("w") as handle:
            handle.write(f">{genome}_contig\n{wrap_fasta(genome_seq)}\n")

    evaluated_negatives = set()
    for cluster, present in PRESENT_IN_MATRIX.items():
        if cluster == "gCore":
            continue
        for genome in GENOMES:
            if genome in present:
                continue
            cell = (cluster, genome)
            if cell not in TRUE_BINARY_POSITIVES:
                evaluated_negatives.add(cell)

    return set(TRUE_BINARY_POSITIVES), evaluated_negatives


def run_recore(
    recore: str,
    input_dir: Path,
    output_dir: Path,
    policy: str,
    env: dict[str, str],
) -> None:
    cmd = [
        recore,
        "-p",
        str(input_dir / "pa.csv"),
        "-r",
        str(input_dir / "pan_ref.fa"),
        "-g",
        str(input_dir / "genomes"),
        "--matrix-policy",
        policy,
        "--min-prevalence",
        "0.0",
        "-t",
        "1",
        "-o",
        str(output_dir),
    ]
    subprocess.run(cmd, cwd=ROOT, check=True, env=env)


def read_matrix(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        return {
            row["Gene"]: {genome: row.get(genome, "") for genome in GENOMES}
            for row in reader
        }


def compute_metrics(
    output_dir: Path,
    positives: set[tuple[str, str]],
    negatives: set[tuple[str, str]],
    policy: str,
) -> dict[str, object]:
    corrected = read_matrix(output_dir / "gene_presence_absence_recore.csv")
    states = read_matrix(output_dir / "gene_presence_absence_recore_states.csv")

    evaluated = positives | negatives
    predicted = {
        cell
        for cell in evaluated
        if corrected[cell[0]][cell[1]].startswith("RESCUED:")
    }

    tp = len(predicted & positives)
    fp = len(predicted & negatives)
    fn = len(positives - predicted)
    tn = len(negatives - predicted)

    precision = tp / (tp + fp) if tp + fp else None
    recall = tp / (tp + fn) if tp + fn else None
    if precision is None:
        f1 = 0.0 if recall == 0 else None
    elif precision + (recall or 0.0) == 0:
        f1 = 0.0
    else:
        f1 = 2 * precision * (recall or 0.0) / (precision + (recall or 0.0))

    state_counts = Counter(states[cluster][genome] or "ABSENT" for cluster, genome in evaluated)

    return {
        "policy": policy,
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "evaluated_positives": len(positives),
        "evaluated_negatives": len(negatives),
        "sequence_only_states": state_counts.get("SEQUENCE_ONLY", 0),
        "pseudogene_states": state_counts.get("PSEUDOGENE", 0),
        "absent_states": state_counts.get("ABSENT", 0),
    }


def write_metrics(metrics: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "policy",
        "tp",
        "fp",
        "fn",
        "tn",
        "precision",
        "recall",
        "f1",
        "evaluated_positives",
        "evaluated_negatives",
        "sequence_only_states",
        "pseudogene_states",
        "absent_states",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in metrics:
            out = dict(row)
            for key in ("precision", "recall", "f1"):
                out[key] = metric_fmt(out[key])  # type: ignore[arg-type]
            writer.writerow(out)


def plot_value(value: object) -> float:
    if value is None:
        return 0.0
    return float(value)


def render_svg(metrics: list[dict[str, object]], path: Path) -> None:
    width = 760
    height = 430
    left = 82
    right = 28
    top = 72
    plot_height = 235
    bottom = top + plot_height
    group_width = 190
    bar_width = 28
    bar_gap = 9
    colors = {
        "precision": "#276EF1",
        "recall": "#12B886",
        "f1": "#7C3AED",
    }
    labels = {
        "precision": "Precision",
        "recall": "Recall",
        "f1": "F1",
    }

    def y(value: float) -> float:
        return bottom - value * plot_height

    lines = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        "<style>text{font-family:DejaVu Sans,Arial,sans-serif;fill:#1f2933}.axis{stroke:#344054;stroke-width:1}.grid{stroke:#d8dee9;stroke-width:1}.label{font-size:12px}.small{font-size:11px}.title{font-size:17px;font-weight:700}</style>",
        f"<text x='{left}' y='28' class='title'>Synthetic benchmark: binary matrix correction</text>",
        f"<text x='{left}' y='49' class='small'>Intact hidden genes are positives; partial pseudogenes are negatives for binary presence.</text>",
    ]

    for tick in (0.0, 0.25, 0.5, 0.75, 1.0):
        yy = y(tick)
        lines.append(f"<line x1='{left}' x2='{width - right}' y1='{yy:.1f}' y2='{yy:.1f}' class='grid'/>")
        lines.append(
            f"<text x='{left - 12}' y='{yy + 4:.1f}' text-anchor='end' class='small'>{tick:.2f}</text>"
        )

    lines.append(f"<line x1='{left}' x2='{left}' y1='{top}' y2='{bottom}' class='axis'/>")
    lines.append(
        f"<line x1='{left}' x2='{width - right}' y1='{bottom}' y2='{bottom}' class='axis'/>"
    )

    for gi, row in enumerate(metrics):
        group_x = left + 55 + gi * group_width
        for mi, key in enumerate(("precision", "recall", "f1")):
            value = plot_value(row[key])
            x = group_x + mi * (bar_width + bar_gap)
            yy = y(value)
            height_bar = bottom - yy
            lines.append(
                f"<rect x='{x}' y='{yy:.1f}' width='{bar_width}' height='{height_bar:.1f}' rx='3' fill='{colors[key]}'/>"
            )
            label = "n/a" if row[key] is None else f"{float(row[key]):.2f}"
            lines.append(
                f"<text x='{x + bar_width / 2:.1f}' y='{yy - 7:.1f}' text-anchor='middle' class='small'>{label}</text>"
            )
        policy = html.escape(str(row["policy"]))
        center = group_x + (3 * bar_width + 2 * bar_gap) / 2
        lines.append(
            f"<text x='{center:.1f}' y='{bottom + 25}' text-anchor='middle' class='label'>{policy}</text>"
        )
        confusion = f"TP {row['tp']}  FP {row['fp']}  FN {row['fn']}"
        lines.append(
            f"<text x='{center:.1f}' y='{bottom + 44}' text-anchor='middle' class='small'>{confusion}</text>"
        )

    legend_x = width - 290
    for i, key in enumerate(("precision", "recall", "f1")):
        x = legend_x + i * 92
        lines.append(f"<rect x='{x}' y='24' width='13' height='13' rx='2' fill='{colors[key]}'/>")
        lines.append(f"<text x='{x + 19}' y='35' class='small'>{labels[key]}</text>")

    footer = (
        "Generated by benchmarks/synthetic_benchmark.py; balanced keeps "
        "pseudogenes in the multistate audit instead of the binary matrix."
    )
    lines.append(
        f"<text x='{left}' y='{height - 24}' class='small'>{html.escape(footer)}</text>"
    )
    lines.append("</svg>")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--recore",
        default=str(ROOT / "target" / "debug" / "recore"),
        help="Path to the Recore binary. Built with cargo build if missing.",
    )
    parser.add_argument(
        "--minimap2",
        default=None,
        help="Path to minimap2. If omitted, PATH and common conda paths are checked.",
    )
    parser.add_argument(
        "--workdir",
        default=str(DEFAULT_WORKDIR),
        help="Temporary working directory for generated inputs and Recore outputs.",
    )
    parser.add_argument(
        "--metrics",
        default=str(DEFAULT_METRICS),
        help="Output TSV with precision/recall/F1 metrics.",
    )
    parser.add_argument(
        "--plot",
        default=str(DEFAULT_PLOT),
        help="Output SVG plot for README embedding.",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Use the existing target/debug/recore binary without running cargo build.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    recore = ensure_recore_binary(args.recore, args.skip_build)
    minimap2 = find_minimap2(args.minimap2)

    workdir = Path(args.workdir)
    if workdir.exists():
        shutil.rmtree(workdir)
    input_dir = workdir / "input"
    output_root = workdir / "outputs"

    positives, negatives = write_inputs(input_dir)
    env = os.environ.copy()
    env["PATH"] = f"{Path(minimap2).parent}{os.pathsep}{env.get('PATH', '')}"

    metrics = []
    for policy in POLICIES:
        output_dir = output_root / policy
        run_recore(recore, input_dir, output_dir, policy, env)
        metrics.append(compute_metrics(output_dir, positives, negatives, policy))

    write_metrics(metrics, Path(args.metrics))
    render_svg(metrics, Path(args.plot))

    print(f"Wrote metrics: {args.metrics}")
    print(f"Wrote plot: {args.plot}")
    print(f"Workdir: {workdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
