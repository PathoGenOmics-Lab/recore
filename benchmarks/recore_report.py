#!/usr/bin/env python3
"""Generate a compact Markdown/HTML report from a Recore output directory."""

from __future__ import annotations

import argparse
import csv
import html
import json
import sys
from collections import Counter
from pathlib import Path


def read_summary(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    out = {}
    for line in path.read_text().splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            out[key.strip().lower().replace(" ", "_")] = value.strip()
    return out


def count_column(path: Path, column: str) -> Counter:
    counts = Counter()
    if not path.exists():
        return counts
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            counts[row.get(column, "") or "NA"] += 1
    return counts


def top_rows(path: Path, sort_column: str, n: int = 12) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    def score(row: dict[str, str]) -> float:
        try:
            return float(row.get(sort_column, "0") or 0)
        except ValueError:
            return 0.0
    return sorted(rows, key=score, reverse=True)[:n]


def md_table(rows: list[list[str]]) -> str:
    if not rows:
        return ""
    header = rows[0]
    body = rows[1:]
    out = ["| " + " | ".join(header) + " |"]
    out.append("| " + " | ".join(["---"] * len(header)) + " |")
    for row in body:
        out.append("| " + " | ".join(str(x) for x in row) + " |")
    return "\n".join(out)


def write_html(markdown: str, path: Path) -> None:
    escaped = html.escape(markdown)
    body = escaped.replace("\n", "<br>\n")
    path.write_text(
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<title>Recore report</title>"
        "<style>body{font-family:Arial,sans-serif;max-width:980px;margin:40px auto;"
        "line-height:1.45}code{background:#f3f3f3;padding:2px 4px}</style>"
        "</head><body><pre style='white-space:pre-wrap'>"
        f"{body}</pre></body></html>\n"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument("--html", type=Path)
    parser.add_argument("--top", type=int, default=12)
    args = parser.parse_args()

    run_dir = args.run_dir
    output = args.output or run_dir / "recore_report.md"
    html_output = args.html or run_dir / "recore_report.html"

    summary = read_summary(run_dir / "summary.txt")
    manifest = {}
    manifest_path = run_dir / "run_manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())

    ortho_actions = count_column(run_dir / "orthology_evidence.tsv", "matrix_action")
    ortho_tiers = count_column(run_dir / "orthology_evidence.tsv", "orthology_confidence_tier")
    priorities = count_column(run_dir / "evolutionary_events.tsv", "highlight_priority")
    mechanisms = count_column(run_dir / "evolutionary_events.tsv", "molecular_mechanism")
    qc = count_column(run_dir / "genome_qc.tsv", "qc_status")

    lines = [
        "# Recore Report",
        "",
        f"Run directory: `{run_dir}`",
        f"Tool version: `{manifest.get('tool_version', 'unknown')}`",
        f"Schema: `{manifest.get('schema_version', 'unknown')}`",
        "",
        "## Headline Counts",
        "",
        md_table(
            [
                ["Metric", "Value"],
                ["Clusters", summary.get("clusters", manifest.get("counts", {}).get("clusters", ""))],
                ["Genomes", summary.get("genomes", manifest.get("counts", {}).get("genomes", ""))],
                ["Rescued", summary.get("rescued", manifest.get("counts", {}).get("rescued", ""))],
                ["Matrix corrected", summary.get("matrix-corrected", manifest.get("counts", {}).get("matrix_corrected_total", ""))],
                ["Core genome", summary.get("core_genome_(strict)", "")],
            ]
        ),
        "",
        "## Orthology Evidence",
        "",
        md_table([["Matrix action", "Rows"], *[[k, v] for k, v in sorted(ortho_actions.items())]]),
        "",
        md_table([["Confidence tier", "Rows"], *[[k, v] for k, v in sorted(ortho_tiers.items())]]),
        "",
        "## Evolutionary Prioritization",
        "",
        md_table([["Highlight priority", "Clusters"], *[[k, v] for k, v in sorted(priorities.items())]]),
        "",
        md_table([["Mechanism", "Clusters"], *[[k, v] for k, v in mechanisms.most_common(15)]]),
        "",
        "## Genome QC",
        "",
        md_table([["QC status", "Genomes"], *[[k, v] for k, v in sorted(qc.items())]]),
        "",
        "## Top Candidate Clusters",
        "",
    ]

    top = top_rows(run_dir / "evolutionary_events.tsv", "evolutionary_signal_score", args.top)
    lines.append(
        md_table(
            [["Gene cluster", "Priority", "Score", "Mechanism", "Support", "Recommendation"]]
            + [
                [
                    row.get("gene_cluster", ""),
                    row.get("highlight_priority", ""),
                    row.get("evolutionary_signal_score", ""),
                    row.get("molecular_mechanism", ""),
                    row.get("support_level", ""),
                    row.get("recommended_wording", ""),
                ]
                for row in top
            ]
        )
    )
    lines.extend(
        [
            "",
            "## Files To Inspect",
            "",
            "- `orthology_evidence.tsv`: per-rescue decision, adaptive score, and risk flags.",
            "- `evolutionary_events.tsv`: cluster-level mechanism, priority, and phylogenetic flags.",
            "- `genome_qc.tsv`: genome-level assembly/annotation/rescue burden flags.",
            "- `cluster_components.tsv`: related cluster components.",
        ]
    )

    markdown = "\n".join(lines) + "\n"
    output.write_text(markdown)
    write_html(markdown, html_output)
    print(f"Wrote {output}")
    print(f"Wrote {html_output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
