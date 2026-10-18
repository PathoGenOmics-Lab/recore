#!/usr/bin/env python3
"""Measure how long Recore takes, and how much memory and disk it uses, on a
pangenome you have.

Give it what your pangenome tool wrote and the assemblies, the way you would
give them to Recore. A Newick tree keeps only the genomes that are its tips,
and --max-clusters keeps the first N clusters, to time a subset; both need a
Panaroo or Roary gene_presence_absence.csv. The inputs are linked, not copied,
and the run's wall time, peak memory of the whole process tree, temporary disk
and output size are appended to a TSV.

It is a resource benchmark, not a truth set: synthetic_benchmark.py and
self_benchmark.py measure whether the answers are right.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WORKDIR = Path(tempfile.gettempdir()) / "recore_performance_benchmark"
DEFAULT_METRICS = REPO_ROOT / "benchmarks" / "performance_metrics.tsv"
FASTA_EXTENSIONS = (".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz")
GFF_EXTENSIONS = (".gff", ".gff3", ".gff.gz", ".gff3.gz")

# The columns of a Panaroo or Roary gene_presence_absence.csv that are not
# genomes.
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


def read_tips(tree_path: Path) -> set[str]:
    """The tip labels of a Newick tree: the names right after '(' or ','."""
    text = tree_path.read_text()
    tips = {label.strip("'\"") for label in re.findall(r"[(,]\s*([^():,;\s]+)", text)}
    if not tips:
        raise SystemExit(f"ERROR: no tip labels found in {tree_path}")
    return tips


def matrix_genomes(matrix: Path) -> tuple[list[str], list[str]]:
    with matrix.open(newline="") as handle:
        header = next(csv.reader(handle))
    return header, [col for col in header if col not in METADATA_COLS]


def filter_matrix(
    matrix: Path,
    keep: set[str] | None,
    output: Path,
    max_clusters: int | None,
) -> tuple[list[str], int]:
    """Write the matrix with only the genomes in `keep` (all of them when None)
    and at most `max_clusters` rows, recounting the isolate columns."""
    header, genomes = matrix_genomes(matrix)
    if header[:1] != ["Gene"]:
        raise SystemExit(
            f"ERROR: {matrix} is not a Panaroo or Roary gene_presence_absence.csv; "
            "--tree and --max-clusters only work on that format"
        )
    kept = [g for g in genomes if keep is None or g in keep]
    if not kept:
        raise SystemExit("ERROR: no matrix column is a tip of the tree")
    if keep is not None and len(kept) < len(keep):
        missing = sorted(keep - set(kept))
        print(f"Note: {len(missing)} tree tips are not matrix columns, e.g. {', '.join(missing[:3])}")
    keep_indices = [i for i, col in enumerate(header) if col in METADATA_COLS or col in kept]
    genome_indices = [i for i in keep_indices if header[i] not in METADATA_COLS]
    recount = [header.index(c) for c in ("No. isolates", "No. sequences") if c in header]
    average = header.index("Avg sequences per isolate") if "Avg sequences per isolate" in header else None

    output.parent.mkdir(parents=True, exist_ok=True)
    n_rows = 0
    with matrix.open(newline="") as handle, output.open("w", newline="") as out:
        reader = csv.reader(handle)
        next(reader)
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow([header[i] for i in keep_indices])
        for row in reader:
            if max_clusters is not None and n_rows >= max_clusters:
                break
            kept_row = [row[i] if i < len(row) else "" for i in keep_indices]
            present = sum(1 for i in genome_indices if i < len(row) and row[i].strip())
            for column in recount:
                kept_row[keep_indices.index(column)] = str(present)
            if average is not None:
                kept_row[keep_indices.index(average)] = "1.00" if present else "0.00"
            writer.writerow(kept_row)
            n_rows += 1
    return kept, n_rows


def find_input(directory: Path, genome: str, extensions: tuple[str, ...], label: str) -> Path:
    for extension in extensions:
        path = directory / f"{genome}{extension}"
        if path.exists():
            return path
    raise SystemExit(f"ERROR: no {label} for {genome} in {directory} (tried {', '.join(extensions)})")


def link_inputs(genomes: list[str], fasta_dir: Path, gff_dir: Path | None, input_root: Path) -> None:
    """Link each genome's FASTA, and its GFF when there is one, named after its
    matrix column the way Recore expects."""
    for name, source, extensions in (("genomes", fasta_dir, FASTA_EXTENSIONS), ("annotations", gff_dir, GFF_EXTENSIONS)):
        if source is None:
            continue
        target = input_root / name
        target.mkdir(parents=True, exist_ok=True)
        for genome in genomes:
            found = find_input(source, genome, extensions, "FASTA" if name == "genomes" else "GFF")
            link = target / found.name
            if not link.exists():
                link.symlink_to(found.resolve())


def recore_binary(explicit: str | None, skip_build: bool) -> Path:
    if explicit:
        binary = Path(explicit).resolve()
    else:
        binary = REPO_ROOT / "target" / "release" / "recore"
        if not skip_build:
            subprocess.run(["cargo", "build", "--release"], cwd=REPO_ROOT, check=True)
    if not binary.exists():
        raise SystemExit(f"ERROR: Recore binary not found: {binary}")
    return binary


def find_minimap2(explicit: str | None) -> Path:
    found = explicit or shutil.which("minimap2")
    if not found or not Path(found).exists():
        raise SystemExit("ERROR: minimap2 not found on PATH. Pass --minimap2 /path/to/minimap2.")
    return Path(found).resolve()


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-p", "--matrix", required=True, help="the presence/absence matrix your pangenome tool wrote")
    parser.add_argument("-r", "--reference", required=True, help="its representative sequences, pan_genome_reference.fa")
    parser.add_argument("-g", "--genomes", required=True, help="one assembly per genome, named like the matrix columns")
    parser.add_argument("-a", "--annotations", help="one GFF per genome; the publication profile needs it")
    parser.add_argument("--tree", help="a Newick tree: keep only the genomes that are its tips")
    parser.add_argument("--max-clusters", type=int, help="keep only the first N clusters of the matrix")
    parser.add_argument("--profile", choices=("publication", "rescue_only"), default="publication",
                        help="publication adds --synteny, --multi-allele and --nscluster (default)")
    parser.add_argument("--preset", help="passed to Recore; its default when omitted")
    parser.add_argument("--min-prevalence", help="passed to Recore; its default when omitted")
    parser.add_argument("--reference-genome", help="passed to Recore, to audit a reference column")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--recore", help="the Recore binary to time (default: build target/release/recore)")
    parser.add_argument("--skip-build", action="store_true", help="use target/release/recore as it is")
    parser.add_argument("--minimap2", help="the minimap2 to put first on PATH (default: the one on PATH)")
    parser.add_argument("--workdir", default=str(DEFAULT_WORKDIR))
    parser.add_argument("--metrics", default=str(DEFAULT_METRICS))
    parser.add_argument("--keep-workdir", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    matrix = Path(args.matrix).resolve()
    reference = Path(args.reference).resolve()
    fasta_dir = Path(args.genomes).resolve()
    gff_dir = Path(args.annotations).resolve() if args.annotations else None
    for label, path in (("matrix", matrix), ("reference", reference), ("genomes", fasta_dir), ("annotations", gff_dir)):
        if path is not None and not path.exists():
            raise SystemExit(f"ERROR: missing {label}: {path}")
    if args.profile == "publication" and gff_dir is None:
        raise SystemExit("ERROR: the publication profile checks neighbours and needs -a/--annotations")

    recore = recore_binary(args.recore, args.skip_build)
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

    if args.tree or args.max_clusters is not None:
        keep = read_tips(Path(args.tree)) if args.tree else None
        run_matrix = input_dir / matrix.name
        genomes, n_clusters = filter_matrix(matrix, keep, run_matrix, args.max_clusters)
    else:
        run_matrix = matrix
        _header, genomes = matrix_genomes(matrix) if matrix.suffix == ".csv" else ([], [])
        with matrix.open() as handle:
            n_clusters = max(0, sum(1 for _ in handle) - 1)
    if genomes:
        link_inputs(genomes, fasta_dir, gff_dir, input_dir)
        genomes_arg, annotations_arg = input_dir / "genomes", input_dir / "annotations"
    else:
        # A format whose columns this script does not parse: hand Recore the
        # folders as they are.
        genomes_arg, annotations_arg = fasta_dir, gff_dir

    cmd = [str(recore), "-p", str(run_matrix), "-r", str(reference), "-g", str(genomes_arg)]
    for flag, value in (("--preset", args.preset), ("--min-prevalence", args.min_prevalence),
                        ("--reference-genome", args.reference_genome)):
        if value:
            cmd += [flag, value]
    if args.profile == "publication":
        cmd += ["-a", str(annotations_arg), "--synteny", "--multi-allele", "--nscluster", "--nsc-engine", "minhash"]
    cmd += ["--matrix-policy", "balanced", "-t", str(args.threads), "-o", str(output_dir)]

    env = os.environ.copy()
    env["PATH"] = f"{minimap2.parent}{os.pathsep}{env.get('PATH', '')}"
    env["TMPDIR"] = str(tmpdir)

    print(f"Timing Recore: profile={args.profile}, clusters={n_clusters}, genomes={len(genomes) or 'as given'}")
    print("Command:", " ".join(cmd))

    monitored = run_monitored(cmd, env, workdir, log_dir, output_dir)
    summary = parse_summary(output_dir / "summary.txt")
    row = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "profile": args.profile,
        "threads": args.threads,
        "clusters": n_clusters,
        "genomes": len(genomes) or "",
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
