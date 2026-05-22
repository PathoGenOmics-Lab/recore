# Synthetic benchmark

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
