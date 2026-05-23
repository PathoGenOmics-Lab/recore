# Benchmarks

## Synthetic benchmark

This benchmark is a reproducible sanity check for Recore's binary matrix
policies. It generates a small pangenome with known ground truth:

- intact genes hidden in assemblies but marked absent in the matrix
- truly absent genes
- partial pseudogenes that should be audited but not counted as clean binary
  presence under the balanced policy

Run it from the repository root:

```bash
python3 benchmarks/synthetic_benchmark.py
```

If `minimap2` is installed outside `PATH`, pass it explicitly:

```bash
python3 benchmarks/synthetic_benchmark.py --minimap2 /path/to/minimap2
```

Outputs:

- `benchmarks/synthetic_benchmark_metrics.tsv`
- `assets/synthetic_benchmark_precision_recall.svg`

Interpretation:

- `strict` requires annotation/synteny-grade evidence, so this annotation-free
  benchmark has no binary corrections.
- `balanced` rescues full-length sequence evidence while keeping partial
  pseudogenes audit-only.
- `sensitive` allows partial/disrupted evidence at the policy layer, but the
  adaptive quality guard still blocks the weak pseudogene-like controls in the
  current benchmark.

This is a regression benchmark for Recore behavior, not a species-specific
biological truth set.

## Biological event benchmark

`biological_event_benchmark.py` is a small GFF/synteny truth set for the
guardrails that matter most in real pangenomes:

- a same-locus sequence with no CDS annotation should be rescued as
  `MISSING_CDS` when synteny supports it
- an identical sequence at a wrong/paralogous context should be rejected to
  `synteny_rejected_report.tsv`
- `evolutionary_events.tsv` should expose quantitative synteny, ORF integrity,
  phylogenetic reliability, and paralog-risk columns

Run it from the repository root:

```bash
python3 benchmarks/biological_event_benchmark.py
```

If `minimap2` is installed outside `PATH`, pass it explicitly:

```bash
python3 benchmarks/biological_event_benchmark.py --minimap2 /path/to/minimap2
```

Output:

- `benchmarks/biological_event_benchmark_metrics.tsv`

This is a regression benchmark for Recore's locus-context and audit logic, not
a species-specific biological truth set.

## Dataset self benchmark

`self_benchmark.py` masks real present cells from a user-provided matrix,
runs Recore on the masked matrix, and estimates how often those known-present
cells are recovered. It also writes score-calibration bins from
`orthology_confidence_score`.

```bash
python3 benchmarks/self_benchmark.py \
  -p /path/to/gene_presence_absence.csv \
  -r /path/to/pan_genome_reference.fa \
  -g /path/to/genomes \
  -a /path/to/annotations \
  --threads 8 \
  --multi-allele \
  --nscluster \
  --minimap2 /path/to/minimap2
```

Outputs in the workdir:

- `self_benchmark_truth.tsv`
- `self_benchmark_calibration.tsv`
- `self_benchmark_score_bins.tsv`
- Recore output directory for the masked run

The negative-control cells are a proxy only: original matrix absences can still
include true false absences. `precision_proxy` is reported as `NA` when no
negative-control cells are available. Existing workdirs are cleaned by default
to avoid stale files; add `--keep-workdir` for manual inspection/reuse.

## Run comparison and reporting

Compare two Recore output directories:

```bash
python3 benchmarks/compare_runs.py --old old_recore_output --new new_recore_output -o compare_out
```

The comparison uses the union of matrix cells, so clusters or genomes present
in only one run are reported as `old_only_cell` / `new_only_cell` instead of
being silently ignored.

Generate a compact Markdown/HTML report:

```bash
python3 benchmarks/recore_report.py recore_output
```

## AMAP performance benchmark

`amap_performance_benchmark.py` profiles Recore on the full local AMAP
pangenome. It is a resource benchmark, not a precision/recall truth set.

The script:

- filters `AMAP-gene_presence_absence.csv` with the canonical 459-tip tree
- links FASTA and GFF inputs instead of copying them
- runs the release binary with the publication settings
- records wall time, peak RSS, temporary disk use, output size, rescues, and
  corrected matrix entries

Run from the repository root:

```bash
python3 benchmarks/amap_performance_benchmark.py \
  --threads 8 \
  --profile publication \
  --minimap2 /path/to/minimap2
```

For a quick smoke test:

```bash
python3 benchmarks/amap_performance_benchmark.py \
  --max-clusters 50 \
  --threads 4 \
  --profile rescue_only \
  --minimap2 /path/to/minimap2 \
  --metrics /tmp/recore_amap_smoke_metrics.tsv
```

The default AMAP output metrics file is
`benchmarks/amap_performance_metrics.tsv`. Because it records local hardware
and temporary paths, treat it as a local run log rather than a portable
benchmark table.
