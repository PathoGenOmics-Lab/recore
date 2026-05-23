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
- `sensitive` maximizes recall but also turns partial pseudogene evidence into
  binary presence, which creates false positives in this benchmark.

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
