#!/usr/bin/env python3
"""Biological event benchmark for Recore synteny and paralog safeguards."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORKDIR = Path(tempfile.gettempdir()) / "recore_biological_event_benchmark"
DEFAULT_METRICS = ROOT / "benchmarks" / "biological_event_benchmark_metrics.tsv"
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


def coding(seed: str, codons: int = 80) -> str:
    alphabet = ["GCT", "GAA", "TTC", "CGT", "AAC", "GGT", "CTG", "ACC"]
    body = "".join(alphabet[(i + len(seed)) % len(alphabet)] for i in range(codons))
    return "ATG" + body + "TAA"


SEQS = {
    "target": coding("target"),
    "up": coding("up", 40),
    "down": coding("down", 40),
    "other_up": coding("other_up", 40),
    "other_down": coding("other_down", 40),
    "paralog_cluster": coding("target"),
}

PRESENT = {
    "target": {"s1", "s2"},
    "up": {"s1", "s2", "s3", "s4"},
    "down": {"s1", "s2", "s3", "s4"},
    "other_up": {"s4"},
    "other_down": {"s4"},
    "paralog_cluster": {"s4"},
}


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def spacer(length: int = 60) -> str:
    return "ACGT" * (length // 4)


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
    for candidate in [
        Path.home() / "miniconda3" / "envs" / "bioconda" / "bin" / "minimap2",
        Path.home() / "miniconda3" / "bin" / "minimap2",
    ]:
        if candidate.exists():
            return str(candidate)
    raise SystemExit("ERROR: minimap2 not found. Pass --minimap2 /path/to/minimap2.")


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


def write_pa(path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(METADATA_COLS + list(GENOMES))
        for cluster, genomes in PRESENT.items():
            length = len(SEQS[cluster])
            meta = [
                cluster,
                cluster,
                f"benchmark {cluster}",
                str(len(genomes)),
                str(len(genomes)),
                "1.0",
                "0",
                "0",
                "0",
                "",
                str(length),
                str(length),
                str(length),
            ]
            cells = [f"{genome}_{cluster}" if genome in genomes else "" for genome in GENOMES]
            writer.writerow(meta + cells)


def add_feature(rows: list[str], contig: str, start: int, end: int, gene_id: str) -> None:
    rows.append(
        f"{contig}\tbenchmark\tCDS\t{start}\t{end}\t.\t+\t0\tID={gene_id};locus_tag={gene_id}"
    )


def build_genome(
    genome: str,
    include_target: bool,
    annotate_target: bool,
    paralog_locus: bool,
) -> tuple[str, list[str]]:
    contig = f"{genome}_contig"
    seq = spacer()
    gff_rows = ["##gff-version 3"]

    def append_gene(cluster: str, gene_id: str, annotate: bool = True) -> None:
        nonlocal seq
        start = len(seq) + 1
        seq += SEQS[cluster]
        end = len(seq)
        if annotate:
            add_feature(gff_rows, contig, start, end, gene_id)
        seq += spacer()

    append_gene("up", f"{genome}_up")
    if include_target:
        append_gene("target", f"{genome}_target", annotate_target)
    append_gene("down", f"{genome}_down")
    if paralog_locus:
        append_gene("other_up", f"{genome}_other_up")
        append_gene("paralog_cluster", f"{genome}_paralog_cluster")
        append_gene("other_down", f"{genome}_other_down")

    return seq, gff_rows


def write_inputs(input_dir: Path) -> None:
    genomes_dir = input_dir / "genomes"
    annotations_dir = input_dir / "annotations"
    genomes_dir.mkdir(parents=True, exist_ok=True)
    annotations_dir.mkdir(parents=True, exist_ok=True)

    with (input_dir / "pan_ref.fa").open("w") as handle:
        for cluster, seq in SEQS.items():
            handle.write(f">{cluster}\n{wrap(seq)}\n")

    write_pa(input_dir / "pa.csv")

    config = {
        "s1": (True, True, False),
        "s2": (True, True, False),
        "s3": (True, False, False),
        "s4": (False, False, True),
    }
    for genome, (include_target, annotate_target, paralog_locus) in config.items():
        seq, gff_rows = build_genome(genome, include_target, annotate_target, paralog_locus)
        (genomes_dir / f"{genome}.fasta").write_text(f">{genome}_contig\n{wrap(seq)}\n")
        (annotations_dir / f"{genome}.gff3").write_text("\n".join(gff_rows) + "\n")

    (input_dir / "tree.nwk").write_text("((s1:0.1,s2:0.1):0.1,s3:0.2,s4:0.3);\n")


def run_recore(recore: str, input_dir: Path, output_dir: Path, env: dict[str, str]) -> None:
    cmd = [
        recore,
        "-p",
        str(input_dir / "pa.csv"),
        "-r",
        str(input_dir / "pan_ref.fa"),
        "-g",
        str(input_dir / "genomes"),
        "-a",
        str(input_dir / "annotations"),
        "--synteny",
        "--tree",
        str(input_dir / "tree.nwk"),
        "--matrix-policy",
        "balanced",
        "--min-prevalence",
        "0.0",
        "-t",
        "1",
        "-o",
        str(output_dir),
    ]
    subprocess.run(cmd, cwd=ROOT, check=True, env=env)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def add_check(rows: list[dict[str, str]], name: str, ok: bool, detail: str) -> None:
    rows.append({"check": name, "status": "PASS" if ok else "FAIL", "detail": detail})


def evaluate(output_dir: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    rescue_rows = read_tsv(output_dir / "rescue_report.tsv")
    rejected_path = output_dir / "synteny_rejected_report.tsv"
    rejected_rows = read_tsv(rejected_path) if rejected_path.exists() else []
    event_rows = read_tsv(output_dir / "evolutionary_events.tsv")

    same_locus = [
        row for row in rescue_rows if row["gene_cluster"] == "target" and row["genome"] == "s3"
    ]
    add_check(
        rows,
        "same_locus_missing_cds_rescued",
        bool(same_locus)
        and same_locus[0]["rescue_type"] == "MISSING_CDS"
        and same_locus[0]["synteny"] == "PASS"
        and float(same_locus[0]["synteny_score"]) >= 0.5,
        str(same_locus[0] if same_locus else "missing"),
    )

    wrong_locus = [
        row for row in rejected_rows if row["gene_cluster"] == "target" and row["genome"] == "s4"
    ]
    add_check(
        rows,
        "wrong_context_paralog_rejected",
        bool(wrong_locus)
        and wrong_locus[0]["synteny"] == "FAIL"
        and wrong_locus[0]["rescue_type"].startswith(("SPLIT", "PARALOG"))
        and wrong_locus[0]["annotation_match"] == "s4_paralog_cluster",
        str(wrong_locus[0] if wrong_locus else "missing"),
    )

    target_event = next(row for row in event_rows if row["gene_cluster"] == "target")
    required_columns = {
        "phylogenetic_reliability",
        "paralog_locus_risk",
        "synteny_mean_score",
        "synteny_pass_fraction",
        "orf_integrity_summary",
    }
    add_check(
        rows,
        "evolutionary_audit_columns_present",
        required_columns.issubset(target_event),
        ",".join(sorted(required_columns)),
    )
    add_check(
        rows,
        "target_cluster_has_quantitative_synteny",
        target_event["synteny_mean_score"] != ""
        and target_event["synteny_pass_fraction"] != "",
        str(target_event),
    )
    return rows


def write_metrics(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["check", "status", "detail"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recore", default=str(ROOT / "target" / "debug" / "recore"))
    parser.add_argument("--minimap2", default=None)
    parser.add_argument("--workdir", default=str(DEFAULT_WORKDIR))
    parser.add_argument("--metrics", default=str(DEFAULT_METRICS))
    parser.add_argument("--skip-build", action="store_true")
    return parser.parse_args()



WORKDIR_MARKER = ".recore-benchmark-workdir"


def clean_workdir(workdir: Path) -> None:
    """Delete a previous benchmark workdir, and only that.

    The path comes from --workdir, and an unguarded rmtree would take the
    repository, a home directory or any folder passed by mistake with it. Only a
    directory this script created (it carries a marker file), an empty one, or
    an older run's (it holds an input/ folder) is removed."""
    workdir = workdir.resolve()
    if not workdir.exists():
        return
    root = Path(__file__).resolve().parents[1]
    protected = {root, root.parent, Path("/").resolve(), Path.home().resolve()}
    ours = (workdir / WORKDIR_MARKER).exists() or (workdir / "input").is_dir()
    empty = not any(workdir.iterdir())
    if workdir in protected or not (ours or empty):
        raise SystemExit(
            f"ERROR: refusing to delete {workdir}: it is not a workdir this benchmark "
            f"created (no {WORKDIR_MARKER}). Remove it yourself or pick another --workdir."
        )
    shutil.rmtree(workdir)


def mark_workdir(workdir: Path) -> None:
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / WORKDIR_MARKER).touch()

def main() -> int:
    args = parse_args()
    recore = ensure_recore_binary(args.recore, args.skip_build)
    minimap2 = find_minimap2(args.minimap2)
    workdir = Path(args.workdir)
    clean_workdir(workdir)
    mark_workdir(workdir)
    input_dir = workdir / "input"
    output_dir = workdir / "output"
    input_dir.mkdir(parents=True)
    write_inputs(input_dir)
    env = dict(os.environ)
    env["PATH"] = str(Path(minimap2).parent) + ":" + env.get("PATH", "")
    run_recore(recore, input_dir, output_dir, env)
    rows = evaluate(output_dir)
    write_metrics(rows, Path(args.metrics))
    for row in rows:
        print(f"{row['status']}\t{row['check']}")
    return 0 if all(row["status"] == "PASS" for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
