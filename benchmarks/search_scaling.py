#!/usr/bin/env python3
"""How long the search takes as the number of genomes grows.

Generates MTBC-sized synthetic datasets (one 4.4 Mb contig and about 4,000
genes per genome by default), with a Panaroo-style matrix in which every genome
has some false absences (the gene is in the assembly but not in the matrix and
not annotated) and some true ones (the gene is gone from the assembly). Each
size is run three times with the release binary: without a search cache, with
a cold cache and with the warm one, and wall time, peak memory and the rescue
count are written to a TSV.

The numbers are for the machine they ran on; the shape (time per genome, and
what the cache saves) is what carries over.

    python3 benchmarks/search_scaling.py --genomes 100 200 400 --threads 8
"""

from __future__ import annotations

import argparse
import csv
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METRICS = REPO_ROOT / "benchmarks" / "search_scaling_metrics.tsv"
DEFAULT_RICH_METRICS = REPO_ROOT / "benchmarks" / "search_scaling_rich_metrics.tsv"
BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
STOPS = {b"TAA", b"TAG", b"TGA"}


def random_orf(rng: np.random.Generator, codons: int) -> bytes:
    """ATG, `codons` sense codons, TAA."""
    body = BASES[rng.integers(0, 4, size=(codons + 16) * 3)].tobytes()
    sense = [body[i:i + 3] for i in range(0, len(body), 3) if body[i:i + 3] not in STOPS]
    return b"ATG" + b"".join(sense[:codons]) + b"TAA"


PRODUCTS = [
    "DNA-directed RNA polymerase subunit beta",
    "ATP-dependent Clp protease proteolytic subunit",
    "PE family protein",
    "hypothetical protein",
    "ABC transporter ATP-binding protein",
    "transcriptional regulator, TetR family",
]


def gff_attributes(genome: str, i: int, rich: bool) -> str:
    """The CDS attributes: an ID and a locus tag, or a PGAP-style record with a
    name, a product, a long note, cross-references and a protein ID, which is
    what most real annotations carry and what an annotation index holds."""
    tag = f"{genome}_{i:05d}"
    if not rich:
        return f"ID={tag};locus_tag={tag}"
    return (
        f"ID=cds-{tag};Parent=gene-{tag};Name=WP_{i:09d}.1;gene=gen{i:05d};locus_tag={tag};"
        f"old_locus_tag=Rv{i:04d};protein_id=WP_{i:09d}.1;Dbxref=GenBank:WP_{i:09d}.1,GeneID:{1000000 + i};"
        f"product={PRODUCTS[i % len(PRODUCTS)]};"
        "Note=Derived by automated computational analysis using gene prediction method%3A Protein Homology. "
        "GO_function%3A GO%3A0003677 - DNA binding [Evidence IEA]"
    )


def build_dataset(out: Path, n_genomes: int, n_genes: int, gaps: int, seed: int, rich_gff: bool = False) -> dict:
    """Write pa.csv, pan_ref.fa, genomes/ and annotations/ under `out`."""
    rng = np.random.default_rng(seed)
    genes = [random_orf(rng, int(rng.integers(250, 500))) for _ in range(n_genes)]
    names = [f"cl{i:05d}" for i in range(n_genes)]
    spacer_len = 120
    (out / "genomes").mkdir(parents=True, exist_ok=True)
    (out / "annotations").mkdir(parents=True, exist_ok=True)
    with open(out / "pan_ref.fa", "w") as fa:
        for name, seq in zip(names, genes):
            fa.write(f">{name}\n{seq.decode()}\n")

    genomes = [f"G{g:05d}" for g in range(n_genomes)]
    # absent[g] = genes missing from the matrix in genome g; deleted[g] = the
    # subset also missing from the assembly. Half of each genome's gaps are
    # false absences, the other half true ones.
    absent: list[set[int]] = []
    deleted: list[set[int]] = []
    for _ in genomes:
        chosen = rng.choice(n_genes, size=gaps, replace=False)
        absent.append(set(int(i) for i in chosen))
        deleted.append(set(int(i) for i in chosen[: gaps // 2]))

    genome_bp = 0
    for g, genome in enumerate(genomes):
        seq_parts = []
        gff_lines = ["##gff-version 3"]
        pos = 0
        for i, gene in enumerate(genes):
            spacer = BASES[rng.integers(0, 4, size=spacer_len)].tobytes()
            seq_parts.append(spacer)
            pos += spacer_len
            if i in deleted[g]:
                continue
            seq_parts.append(gene)
            if i not in absent[g]:
                gff_lines.append(
                    f"{genome}_1\tsynthetic\tCDS\t{pos + 1}\t{pos + len(gene)}\t.\t+\t0\t"
                    + gff_attributes(genome, i, rich_gff)
                )
            pos += len(gene)
        contig = bytearray(b"".join(seq_parts))
        # About one difference in a thousand, so alignments are not trivial.
        sites = rng.integers(0, len(contig), size=len(contig) // 1000)
        for s in sites:
            contig[s] = BASES[(BASES.tolist().index(contig[s]) + 1) % 4]
        genome_bp = len(contig)
        with open(out / "genomes" / f"{genome}.fasta", "w") as fh:
            fh.write(f">{genome}_1\n")
            text = contig.decode()
            for k in range(0, len(text), 80):
                fh.write(text[k:k + 80] + "\n")
        (out / "annotations" / f"{genome}.gff").write_text("\n".join(gff_lines) + "\n")

    header = (
        "Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,"
        "Avg sequences per isolate,Genome Fragment,Order within Fragment,"
        "Accessory Order,QC,Min group size nuc,Max group size nuc,Avg group size nuc,"
        + ",".join(genomes)
    )
    with open(out / "pa.csv", "w") as fh:
        fh.write(header + "\n")
        for i, name in enumerate(names):
            cells = ["" if i in absent[g] else f"{genomes[g]}_{i:05d}" for g in range(n_genomes)]
            n = sum(1 for c in cells if c)
            fh.write(f"{name},{name},synthetic,{n},{n},1,0,0,0,,0,0,0," + ",".join(cells) + "\n")
    return {
        "genome_mb": round(genome_bp / 1e6, 2),
        "false_absences": sum(len(a) - len(d) for a, d in zip(absent, deleted)),
        "true_absences": sum(len(d) for d in deleted),
    }


def find_minimap2(explicit: str | None) -> Path:
    for candidate in [explicit, shutil.which("minimap2")]:
        if candidate and Path(candidate).exists():
            return Path(candidate)
    raise SystemExit("ERROR: minimap2 not found; pass --minimap2")


def run_recore(binary: Path, data: Path, out: Path, threads: int, cache: Path | None, env: dict) -> dict:
    """One run under /usr/bin/time; wall seconds, peak RSS and rescues."""
    cmd = [
        str(binary), "-p", str(data / "pa.csv"), "-r", str(data / "pan_ref.fa"),
        "-g", str(data / "genomes"), "-a", str(data / "annotations"), "--synteny",
        "-t", str(threads), "-o", str(out),
    ]
    if cache is not None:
        cmd += ["--search-cache", str(cache)]
    timer = ["/usr/bin/time", "-l"] if platform.system() == "Darwin" else ["/usr/bin/time", "-v"]
    start = time.monotonic()
    proc = subprocess.run(timer + cmd, env=env, capture_output=True, text=True)
    seconds = time.monotonic() - start
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr[-3000:])
        raise SystemExit(f"ERROR: recore failed on {data}")
    rss = 0
    m = re.search(r"(\d+)\s+maximum resident set size", proc.stderr)
    if m:
        rss = int(m.group(1)) // (1024 * 1024)
    else:
        m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", proc.stderr)
        if m:
            rss = int(m.group(1)) // 1024
    summary = (out / "summary.txt").read_text()
    rescued = int(re.search(r"Rescued:\s+(\d+)", summary).group(1))
    searched = int(re.search(r"Gaps searched:\s+(\d+)", summary).group(1))
    return {"seconds": round(seconds, 1), "peak_rss_mb": rss, "rescued": rescued, "searched": searched}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genomes", type=int, nargs="+", default=[50, 100, 200])
    parser.add_argument("--genes", type=int, default=4000, help="genes per genome")
    parser.add_argument("--gaps", type=int, default=80, help="matrix gaps per genome, half of them false")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument(
        "--rich-gff",
        action="store_true",
        help="PGAP-style annotations (names, products, notes, cross-references) instead of an ID and a locus tag",
    )
    parser.add_argument("--minimap2")
    parser.add_argument("--recore", help="recore binary (default: target/release/recore, built if missing)")
    parser.add_argument("--workdir", help="where datasets and outputs go (default: a temporary directory)")
    parser.add_argument("--keep-workdir", action="store_true")
    parser.add_argument(
        "--metrics",
        help="output TSV (default: search_scaling_metrics.tsv, or search_scaling_rich_metrics.tsv with --rich-gff)",
    )
    args = parser.parse_args()

    minimap2 = find_minimap2(args.minimap2)
    env = os.environ.copy()
    env["PATH"] = f"{minimap2.parent}{os.pathsep}{env.get('PATH', '')}"
    binary = Path(args.recore) if args.recore else REPO_ROOT / "target" / "release" / "recore"
    if not binary.exists():
        subprocess.run(["cargo", "build", "--release"], cwd=REPO_ROOT, check=True)
    version = subprocess.run([str(binary), "--version"], capture_output=True, text=True).stdout.strip()
    mm2_version = subprocess.run([str(minimap2), "--version"], capture_output=True, text=True).stdout.strip()

    workdir = Path(args.workdir) if args.workdir else Path(tempfile.mkdtemp(prefix="recore_scaling_"))
    workdir.mkdir(parents=True, exist_ok=True)
    machine = f"{platform.machine()} {platform.system()} {os.cpu_count()} cpus"
    rows = []
    try:
        for n in args.genomes:
            data = workdir / f"data_{n}"
            if data.exists():
                shutil.rmtree(data)
            t0 = time.monotonic()
            info = build_dataset(data, n, args.genes, args.gaps, args.seed, args.rich_gff)
            print(f"{n} genomes: dataset in {time.monotonic() - t0:.0f}s ({info['genome_mb']} Mb each)", flush=True)
            cache = workdir / f"cache_{n}"
            if cache.exists():
                shutil.rmtree(cache)
            plain = run_recore(binary, data, workdir / f"out_{n}_plain", args.threads, None, env)
            cold = run_recore(binary, data, workdir / f"out_{n}_cold", args.threads, cache, env)
            warm = run_recore(binary, data, workdir / f"out_{n}_warm", args.threads, cache, env)
            row = {
                "genomes": n,
                "genome_mb": info["genome_mb"],
                "genes_per_genome": args.genes,
                "annotation": "rich" if args.rich_gff else "minimal",
                "gaps_searched": plain["searched"],
                "false_absences": info["false_absences"],
                "rescued": plain["rescued"],
                "threads": args.threads,
                "seconds_no_cache": plain["seconds"],
                "seconds_cold_cache": cold["seconds"],
                "seconds_warm_cache": warm["seconds"],
                "seconds_per_genome": round(plain["seconds"] / n, 3),
                "peak_rss_mb": max(plain["peak_rss_mb"], cold["peak_rss_mb"], warm["peak_rss_mb"]),
                "recore": version,
                "minimap2": mm2_version,
                "machine": machine,
            }
            rows.append(row)
            print("  " + ", ".join(f"{k}={v}" for k, v in row.items() if k not in {"recore", "minimap2", "machine"}), flush=True)
            if not args.keep_workdir:
                shutil.rmtree(data, ignore_errors=True)
                for suffix in ("plain", "cold", "warm"):
                    shutil.rmtree(workdir / f"out_{n}_{suffix}", ignore_errors=True)
                shutil.rmtree(cache, ignore_errors=True)
    finally:
        if not args.keep_workdir and not args.workdir:
            shutil.rmtree(workdir, ignore_errors=True)

    metrics = Path(args.metrics) if args.metrics else DEFAULT_RICH_METRICS if args.rich_gff else DEFAULT_METRICS
    metrics.parent.mkdir(parents=True, exist_ok=True)
    with open(metrics, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {metrics}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
