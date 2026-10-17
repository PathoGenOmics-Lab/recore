#!/usr/bin/env python3
"""Run a reproducible Recore performance benchmark on the local AMAP dataset.

The benchmark uses the PANTERA 459-genome working set by filtering the AMAP
presence/absence matrix with `pantera/data/phylogeny/pantera_459.treefile`.
It creates only symlinks to FASTA/GFF inputs, runs the release Recore binary,
and samples the process tree to estimate wall time, peak RSS, and temporary
disk use.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = REPO_ROOT.parent
DEFAULT_WORKDIR = Path("/private/tmp/recore_amap_performance")
DEFAULT_METRICS = REPO_ROOT / "benchmarks" / "amap_performance_metrics.tsv"
AS_RE = re.compile(r"AS\d{5}")

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


def as_id(value: str) -> str | None:
    match = AS_RE.search(value)
    return match.group(0) if match else None


def read_tree_ids(tree_path: Path) -> tuple[set[str], set[str]]:
    text = tree_path.read_text()
    tips = set(re.findall(r"([A-Za-z0-9]+_AS\d{5}|AS\d{5})", text))
    asids = {m.group(0) for m in AS_RE.finditer(text)}
    if not tips or not asids:
        raise SystemExit(f"ERROR: no AS IDs found in {tree_path}")
    return tips, asids


def ensure_recore_binary(skip_build: bool) -> Path:
    binary = REPO_ROOT / "target" / "release" / "recore"
    if not skip_build:
        subprocess.run(["cargo", "build", "--release"], cwd=REPO_ROOT, check=True)
    if not binary.exists():
        raise SystemExit(f"ERROR: release binary not found: {binary}")
    return binary


def find_minimap2(explicit: str | None) -> Path:
    candidates = []
    if explicit:
        candidates.append(Path(explicit))
    found = shutil.which("minimap2")
    if found:
        candidates.append(Path(found))
    candidates.extend(
        [
            Path.home() / "miniconda3" / "envs" / "bioconda" / "bin" / "minimap2",
            Path.home() / "miniconda3" / "bin" / "minimap2",
        ]
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise SystemExit("ERROR: minimap2 not found. Pass --minimap2 /path/to/minimap2.")


def find_existing(paths: list[Path], label: str) -> Path:
    for path in paths:
        if path.exists():
            return path
    raise SystemExit(f"ERROR: cannot find {label}: " + ", ".join(map(str, paths)))


def source_paths(project_root: Path) -> dict[str, Path]:
    return {
        "pa": project_root / "03pangenome-cds" / "AMAP-gene_presence_absence.csv",
        "pan_ref": project_root / "03pangenome-cds" / "AMAP-pan_genome_reference.fa",
        "tree": project_root / "pantera" / "data" / "phylogeny" / "pantera_459.treefile",
        "genomes": project_root / "01data" / "genomes" / "all_genomes",
        "annotations": project_root / "02annot" / "AMAP",
    }


def filter_pa_matrix(
    pa_path: Path,
    tree_path: Path,
    output_path: Path,
    max_clusters: int | None,
) -> tuple[list[str], int]:
    tips, asids = read_tree_ids(tree_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with pa_path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        keep_indices = []
        kept_genomes = []
        for i, col in enumerate(header):
            if col in METADATA_COLS:
                keep_indices.append(i)
                continue
            col_as = as_id(col)
            if col in tips or (col_as and col_as in asids):
                keep_indices.append(i)
                kept_genomes.append(col)

        if len(kept_genomes) != len(asids):
            raise SystemExit(
                f"ERROR: kept {len(kept_genomes)} genomes but tree has {len(asids)} AS IDs"
            )

        isolates_col = header.index("No. isolates") if "No. isolates" in header else None
        sequences_col = header.index("No. sequences") if "No. sequences" in header else None
        avg_col = (
            header.index("Avg sequences per isolate")
            if "Avg sequences per isolate" in header
            else None
        )
        kept_genome_indices = [i for i in keep_indices if header[i] not in METADATA_COLS]

        n_rows = 0
        with output_path.open("w", newline="") as out_handle:
            writer = csv.writer(out_handle, lineterminator="\n")
            writer.writerow([header[i] for i in keep_indices])
            for row in reader:
                if max_clusters is not None and n_rows >= max_clusters:
                    break
                kept_row = [row[i] if i < len(row) else "" for i in keep_indices]
                present = sum(
                    1
                    for i in kept_genome_indices
                    if i < len(row) and row[i].strip()
                )
                for source_col in (isolates_col, sequences_col):
                    if source_col is not None and source_col in keep_indices:
                        kept_row[keep_indices.index(source_col)] = str(present)
                if avg_col is not None and avg_col in keep_indices:
                    kept_row[keep_indices.index(avg_col)] = "1.00" if present else "0.00"
                writer.writerow(kept_row)
                n_rows += 1

    return kept_genomes, n_rows


def symlink_inputs(
    kept_genomes: list[str],
    source_genomes: Path,
    source_annotations: Path,
    input_root: Path,
) -> tuple[int, int]:
    genome_dir = input_root / "genomes"
    ann_dir = input_root / "annotations"
    genome_dir.mkdir(parents=True, exist_ok=True)
    ann_dir.mkdir(parents=True, exist_ok=True)
    linked_genomes = 0
    linked_annotations = 0

    for genome in kept_genomes:
        genome_as = as_id(genome)
        if not genome_as:
            continue
        genome_src = find_existing(
            [
                source_genomes / f"{genome_as}.fna",
                source_genomes / f"{genome_as}.fa",
                source_genomes / f"{genome_as}.fasta",
                source_genomes / f"{genome_as}.fna.gz",
            ],
            f"FASTA for {genome}",
        )
        genome_dst = genome_dir / genome_src.name
        if not genome_dst.exists():
            genome_dst.symlink_to(genome_src)
        linked_genomes += 1

        ann_src = find_existing(
            [
                source_annotations / f"{genome}.gff",
                source_annotations / f"{genome}.gff3",
                source_annotations / f"{genome_as}.gff",
                source_annotations / f"{genome_as}.gff3",
                source_annotations / f"{genome}.gff.gz",
                source_annotations / f"{genome}.gff3.gz",
            ],
            f"GFF for {genome}",
        )
        ann_dst = ann_dir / ann_src.name
        if not ann_dst.exists():
            ann_dst.symlink_to(ann_src)
        linked_annotations += 1

    return linked_genomes, linked_annotations


def dir_size(path: Path) -> int:
    if not path.exists():
        return 0
    total = 0
    for root, _dirs, files in os.walk(path):
        for file_name in files:
            file_path = Path(root) / file_name
            try:
                if not file_path.is_symlink():
                    total += file_path.stat().st_size
            except FileNotFoundError:
                pass
    return total


def process_tree_rss_kb(root_pid: int) -> int:
    try:
        output = subprocess.run(
            ["ps", "-axo", "pid=,ppid=,rss="],
            check=True,
            text=True,
            capture_output=True,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return 0
    children: dict[int, list[int]] = {}
    rss: dict[int, int] = {}
    for line in output.splitlines():
        parts = line.split()
        if len(parts) != 3:
            continue
        pid, ppid, rss_kb = map(int, parts)
        children.setdefault(ppid, []).append(pid)
        rss[pid] = rss_kb

    stack = [root_pid]
    total = 0
    seen = set()
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        total += rss.get(pid, 0)
        stack.extend(children.get(pid, []))
    return total


def parse_summary(path: Path) -> dict[str, str]:
    summary = {}
    if not path.exists():
        return summary
    for line in path.read_text(errors="replace").splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        key = key.strip().lower().replace(" ", "_")
        value = value.strip()
        if key:
            summary[key] = value
    return summary


def run_monitored(
    cmd: list[str],
    env: dict[str, str],
    cwd: Path,
    log_dir: Path,
    output_dir: Path,
) -> dict[str, float | int]:
    log_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = log_dir / "stdout.log"
    stderr_path = log_dir / "stderr.log"
    start = time.monotonic()
    peak_rss_kb = 0
    peak_tmp_bytes = 0
    peak_output_bytes = 0
    tmpdir = Path(env["TMPDIR"])

    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=stdout, stderr=stderr)
        last_disk_sample = 0.0
        while proc.poll() is None:
            try:
                peak_rss_kb = max(peak_rss_kb, process_tree_rss_kb(proc.pid))
            except subprocess.SubprocessError:
                pass
            now = time.monotonic()
            if now - last_disk_sample >= 2.0:
                peak_tmp_bytes = max(peak_tmp_bytes, dir_size(tmpdir))
                peak_output_bytes = max(peak_output_bytes, dir_size(output_dir))
                last_disk_sample = now
            time.sleep(0.5)
        return_code = proc.wait()

    peak_tmp_bytes = max(peak_tmp_bytes, dir_size(tmpdir))
    peak_output_bytes = max(peak_output_bytes, dir_size(output_dir))
    runtime = time.monotonic() - start
    if return_code != 0:
        raise SystemExit(
            f"ERROR: Recore failed with exit code {return_code}. See {stderr_path}"
        )
    return {
        "runtime_seconds": runtime,
        "peak_rss_mb": peak_rss_kb / 1024.0,
        "peak_tmp_mb": peak_tmp_bytes / (1024.0 * 1024.0),
        "peak_output_mb": peak_output_bytes / (1024.0 * 1024.0),
    }


def write_metrics(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "timestamp",
        "profile",
        "threads",
        "clusters",
        "genomes",
        "runtime_seconds",
        "peak_rss_mb",
        "peak_tmp_mb",
        "peak_output_mb",
        "rescued",
        "matrix_corrected",
        "core_genome",
        "output_dir",
    ]
    exists = path.exists()
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        if not exists:
            writer.writeheader()
        writer.writerow({field: row.get(field, "") for field in fields})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=str(PROJECT_ROOT))
    parser.add_argument("--workdir", default=str(DEFAULT_WORKDIR))
    parser.add_argument("--metrics", default=str(DEFAULT_METRICS))
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--minimap2")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--max-clusters", type=int)
    parser.add_argument(
        "--profile",
        choices=("publication", "rescue_only"),
        default="publication",
    )
    parser.add_argument("--keep-workdir", action="store_true")
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
    project_root = Path(args.project_root).resolve()
    paths = source_paths(project_root)
    for label, path in paths.items():
        if not path.exists():
            raise SystemExit(f"ERROR: missing {label}: {path}")

    recore = ensure_recore_binary(args.skip_build)
    minimap2 = find_minimap2(args.minimap2)
    workdir = Path(args.workdir).resolve()
    if not args.keep_workdir:
        clean_workdir(workdir)
    mark_workdir(workdir)
    input_dir = workdir / "input"
    output_dir = workdir / f"output_{args.profile}"
    tmpdir = workdir / "tmp"
    log_dir = workdir / "logs"
    tmpdir.mkdir(parents=True, exist_ok=True)

    filtered_pa = input_dir / "AMAP-gene_presence_absence_459.csv"
    kept_genomes, n_clusters = filter_pa_matrix(
        paths["pa"],
        paths["tree"],
        filtered_pa,
        args.max_clusters,
    )
    linked_genomes, linked_annotations = symlink_inputs(
        kept_genomes,
        paths["genomes"],
        paths["annotations"],
        input_dir,
    )

    cmd = [
        str(recore),
        "-p",
        str(filtered_pa),
        "-r",
        str(paths["pan_ref"]),
        "-g",
        str(input_dir / "genomes"),
        "--preset",
        "asm5",
        "--min-prevalence",
        "0.05",
        "--matrix-policy",
        "balanced",
        "--reference-genome",
        "AS00146",
        "-t",
        str(args.threads),
        "-o",
        str(output_dir),
    ]
    if args.profile == "publication":
        cmd.extend(
            [
                "-a",
                str(input_dir / "annotations"),
                "--synteny",
                "--multi-allele",
                "--nscluster",
                "--nsc-engine",
                "minhash",
            ]
        )

    env = os.environ.copy()
    env["PATH"] = f"{minimap2.parent}{os.pathsep}{env.get('PATH', '')}"
    env["TMPDIR"] = str(tmpdir)

    print(f"Running AMAP benchmark: profile={args.profile}, clusters={n_clusters}, genomes={len(kept_genomes)}")
    print(f"Linked inputs: genomes={linked_genomes}, annotations={linked_annotations}")
    print("Command:", " ".join(cmd))

    monitored = run_monitored(cmd, env, REPO_ROOT, log_dir, output_dir)
    summary = parse_summary(output_dir / "summary.txt")
    row = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "profile": args.profile,
        "threads": args.threads,
        "clusters": n_clusters,
        "genomes": len(kept_genomes),
        "runtime_seconds": f"{monitored['runtime_seconds']:.2f}",
        "peak_rss_mb": f"{monitored['peak_rss_mb']:.1f}",
        "peak_tmp_mb": f"{monitored['peak_tmp_mb']:.1f}",
        "peak_output_mb": f"{monitored['peak_output_mb']:.1f}",
        "rescued": summary.get("rescued", ""),
        "matrix_corrected": summary.get("matrix-corrected", ""),
        "core_genome": summary.get("core_genome_(strict)", ""),
        "output_dir": str(output_dir),
    }
    write_metrics(Path(args.metrics), row)

    print("Metrics:")
    for key, value in row.items():
        print(f"  {key}: {value}")
    print(f"Wrote metrics: {args.metrics}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
