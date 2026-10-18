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
being silently ignored. A cell absent in both runs but searched in only one of
them (an `ABSENT` against a `NOT_SEARCHED_*`, as older runs wrote `ABSENT` for
both) is `search_coverage_changed`, not a correction, and is counted apart from
`changed_cells`.

Generate a compact Markdown/HTML report:

```bash
python3 benchmarks/recore_report.py recore_output
```

For new runs, prefer the integrated report flag:

```bash
recore ... --write-html-report -o recore_output
```

## Search scaling benchmark

`search_scaling.py` measures how a run grows with the number of genomes, on
synthetic MTBC-sized datasets: one 4.9 Mb contig and 4,000 genes per genome,
and per genome 40 false absences (gene in the assembly, not in the matrix or
the GFF) and 40 true ones (gene gone). Each size runs three times with the
release binary and `--synteny`: without a search cache, with a cold one and
with the warm one.

```bash
python3 benchmarks/search_scaling.py --genomes 100 200 400 800 --threads 8 \
  --minimap2 /path/to/minimap2
```

It writes `benchmarks/search_scaling_metrics.tsv`: wall seconds for each of the
three runs, seconds per genome, peak memory and the rescue count (which should
equal the false absences planted). The datasets are generated afresh and
removed afterwards unless `--keep-workdir` is given. The numbers belong to the
machine they ran on, which the table records; the shape is what carries over.

The GFFs carry an ID and a locus tag. `--rich-gff` writes PGAP-style records
instead (a name, a gene symbol, an old locus tag, a protein ID,
cross-references, a product and a long note), six names per CDS, which is what
an annotation index has to hold for most real assemblies; that run goes to
`benchmarks/search_scaling_rich_metrics.tsv`.

## Performance on your own pangenome

`performance_benchmark.py` times Recore on a pangenome you have: wall time,
peak memory of the whole process tree, temporary disk, output size, and the
rescue counts, appended to `benchmarks/performance_metrics.tsv`. It is a
resource benchmark, not a truth set.

Give it the same inputs as Recore. A Newick tree keeps only the genomes that are
its tips, and `--max-clusters` keeps the first clusters, to time a subset; both
need a Panaroo or Roary `gene_presence_absence.csv`. FASTA and GFF files are
linked, not copied. The `publication` profile adds `--synteny`,
`--multi-allele` and `--nscluster`; `rescue_only` searches without them.

```bash
python3 benchmarks/performance_benchmark.py \
  -p gene_presence_absence.csv -r pan_genome_reference.fa \
  -g genomes/ -a annotations/ --tree species.nwk --threads 8
```

For a quick smoke test on the tutorial dataset (below):

```bash
python3 benchmarks/performance_benchmark.py \
  -p tutorial_data/gene_presence_absence.csv \
  -r tutorial_data/pan_genome_reference.fa \
  -g tutorial_data/genomes --max-clusters 20 --profile rescue_only \
  --threads 4 --metrics recore_smoke_metrics.tsv
```

It builds `target/release/recore` first unless you pass `--recore` with a
binary or `--skip-build`. The metrics file records the machine and temporary
paths, so treat it as a local run log rather than a portable benchmark table.

## Tutorial dataset

`make_tutorial_dataset.py` writes the small input set the documentation
[tutorial](../docs/getting-started/tutorial.md)
walks through: twelve MTBC-style genomes, twenty-four clusters, GFF3
annotations, a Newick tree, and eleven false absences planted on purpose in four
flavours (dropped gene, missing CDS, pseudogene, and a diverged copy that must
not be rescued).

```bash
python3 benchmarks/make_tutorial_dataset.py --out tutorial_data
```

It is deterministic: the same `--seed` gives the same dataset, so the numbers
quoted in the tutorial stay true.
