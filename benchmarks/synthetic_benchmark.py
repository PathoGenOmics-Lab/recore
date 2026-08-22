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
DEFAULT_PLOT_DARK = ROOT / "assets" / "synthetic_benchmark_precision_recall_dark.svg"

# The plot is drawn once on paper colours and then recoloured for dark pages.
# GitHub picks between the two with <picture media="(prefers-color-scheme: dark)">,
# so the figure stops being a white slab on a dark README. Neutrals invert;
# the semantic accents keep their hue and are lifted just enough to clear the
# 4.5:1 they need against #0d1117.
DARK_PALETTE = {
    "white": "#0d1117",   # page
    "#f8fafc": "#151b23",  # panel
    "#ecfdf3": "#132d21",  # highlighted row
    "#eaecf0": "#2a313c",  # hairlines and bar tracks
    "#101828": "#e6edf3",  # body text
    "#475467": "#a8b0ba",  # subtitle
    "#667085": "#9aa4b0",  # muted labels
    "#d0d5dd": "#4d5866",  # "not applicable" bar
    "#98a2b3": "#8b95a3",  # true-absence swatch
    "#079455": "#3fb950",  # good
    "#7f56d9": "#a879f5",  # partial
    "#dc6803": "#e8912a",  # weak
    "#d92d20": "#f85149",  # bad
}


def recolour_for_dark(svg: str) -> str:
    """Swap the paper palette for the dark one, longest keys first so that a
    colour which is a prefix of another cannot be half-replaced."""
    for old, new in sorted(DARK_PALETTE.items(), key=lambda kv: -len(kv[0])):
        svg = svg.replace(f"'{old}'", f"'{new}'").replace(f":{old}", f":{new}")
    return svg
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
    width = 1120
    height = 520
    panel_x = 42
    panel_w = width - 84

    def score_label(value: object) -> str:
        if value is None:
            return "NA"
        return f"{float(value):.2f}"

    def score_fill(value: object) -> str:
        if value is None:
            return "#d0d5dd"
        value_f = float(value)
        if value_f >= 0.95:
            return "#079455"
        if value_f >= 0.8:
            return "#7f56d9"
        if value_f >= 0.5:
            return "#dc6803"
        return "#d92d20"

    def metric_bar(x: int, y: int, value: object) -> list[str]:
        track_w = 126
        bar_h = 14
        fill_w = track_w * plot_value(value)
        value_text = score_label(value)
        fill = score_fill(value)
        return [
            f"<rect x='{x}' y='{y}' width='{track_w}' height='{bar_h}' rx='6' fill='#eaecf0'/>",
            f"<rect x='{x}' y='{y}' width='{fill_w:.1f}' height='{bar_h}' rx='6' fill='{fill}'/>",
            f"<text x='{x + track_w + 10}' y='{y + 11}' class='tiny strong'>{value_text}</text>",
        ]

    lines = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        (
            "<style>"
            "text{font-family:DejaVu Sans,Arial,sans-serif;fill:#101828}"
            ".title{font-size:22px;font-weight:700}"
            ".subtitle{font-size:13px;fill:#475467}"
            ".section{font-size:12px;font-weight:700;letter-spacing:.08em;fill:#667085}"
            ".label{font-size:13px;font-weight:700}"
            ".small{font-size:12px}"
            ".tiny{font-size:11px}"
            ".muted{fill:#667085}"
            ".white{fill:white}"
            ".strong{font-weight:700}"
            ".row{fill:#f8fafc;stroke:#eaecf0;stroke-width:1}"
            ".balanced{fill:#ecfdf3;stroke:#079455;stroke-width:1.5}"
            "</style>"
        ),
        "<text x='42' y='38' class='title'>Synthetic benchmark: matrix policy trade-off</text>",
        (
            "<text x='42' y='61' class='subtitle'>Full-length false absences should be rescued; "
            "partial pseudogenes should remain audit-only in the binary matrix.</text>"
        ),
    ]

    lines.extend(
        [
            "<text x='42' y='101' class='section'>GROUND TRUTH</text>",
            f"<rect x='{panel_x}' y='115' width='{panel_w}' height='82' rx='12' fill='#f8fafc' stroke='#eaecf0'/>",
            "<rect x='70' y='139' width='16' height='16' rx='4' fill='#079455'/>",
            "<text x='96' y='148' class='label'>5 intact false absences</text>",
            "<text x='96' y='170' class='small muted'>binary positives: should become present</text>",
            "<rect x='420' y='139' width='16' height='16' rx='4' fill='#dc6803'/>",
            "<text x='446' y='148' class='label'>2 partial pseudogene traps</text>",
            "<text x='446' y='170' class='small muted'>biological evidence, not clean binary presence</text>",
            "<rect x='802' y='139' width='16' height='16' rx='4' fill='#98a2b3'/>",
            "<text x='828' y='148' class='label'>8 true absences</text>",
            "<text x='828' y='170' class='small muted'>binary negatives: should stay absent</text>",
            "<text x='42' y='235' class='section'>POLICY PERFORMANCE</text>",
            "<text x='70' y='262' class='tiny muted'>Policy</text>",
            "<text x='350' y='262' class='tiny muted'>Precision</text>",
            "<text x='545' y='262' class='tiny muted'>Recall</text>",
            "<text x='740' y='262' class='tiny muted'>F1</text>",
            "<text x='940' y='262' class='tiny muted'>Confusion</text>",
        ]
    )

    policy_notes = {
        "strict": "annotation/synteny-grade only",
        "balanced": "paper-ready binary matrix",
        "sensitive": "exploratory, accepts partial evidence",
    }
    row_y = {"strict": 278, "balanced": 346, "sensitive": 414}
    for row in metrics:
        policy = str(row["policy"])
        y = row_y[policy]
        klass = "balanced" if policy == "balanced" else "row"
        lines.append(
            f"<rect x='{panel_x}' y='{y}' width='{panel_w}' height='54' rx='12' class='{klass}'/>"
        )
        lines.append(
            f"<text x='70' y='{y + 23}' class='label'>{html.escape(policy)}</text>"
        )
        lines.append(
            f"<text x='70' y='{y + 42}' class='tiny muted'>{html.escape(policy_notes[policy])}</text>"
        )
        if policy == "balanced":
            lines.append(
                f"<rect x='190' y='{y + 12}' width='104' height='22' rx='11' fill='#079455'/>"
            )
            lines.append(
                f"<text x='242' y='{y + 27}' text-anchor='middle' class='tiny white'>recommended</text>"
            )
        lines.extend(metric_bar(350, y + 21, row["precision"]))
        lines.extend(metric_bar(545, y + 21, row["recall"]))
        lines.extend(metric_bar(740, y + 21, row["f1"]))
        confusion = f"TP {row['tp']}  FP {row['fp']}  FN {row['fn']}"
        lines.append(
            f"<text x='940' y='{y + 26}' class='small strong'>{html.escape(confusion)}</text>"
        )
        lines.append(
            f"<text x='940' y='{y + 44}' class='tiny muted'>TN {row['tn']}</text>"
        )

    lines.extend(
        [
            "<line x1='42' x2='1078' y1='494' y2='494' stroke='#eaecf0'/>",
            (
                "<text x='42' y='512' class='tiny muted'>Generated by benchmarks/synthetic_benchmark.py. "
                "The balanced policy preserves disrupted homologues in the multistate audit instead of the binary matrix.</text>"
            ),
        ]
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
        "--plot-dark",
        default=str(DEFAULT_PLOT_DARK),
        help="Output SVG plot recoloured for dark backgrounds.",
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
    Path(args.plot_dark).write_text(recolour_for_dark(Path(args.plot).read_text()))

    print(f"Wrote metrics: {args.metrics}")
    print(f"Wrote plot: {args.plot}")
    print(f"Wrote dark plot: {args.plot_dark}")
    print(f"Workdir: {workdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
