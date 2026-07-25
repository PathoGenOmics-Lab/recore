---
hide:
  - navigation
---

<div class="hero" markdown>

![Recore](assets/recore_logo.svg){ .hero-logo }

# Recore

Post-hoc rescue, audit, and quality control of pangenome presence/absence
matrices.

</div>

Pangenome pipelines can mark a gene as absent even when the sequence is still
present in the assembly. These false absences inflate the accessory genome,
shrink the core genome, and create lineage-specific signals that are artefacts
of graph cleaning rather than biology.

Recore searches every absent cell against the genome assemblies with minimap2,
classifies the evidence, validates orthology against annotation and synteny when
they are available, and writes both:

- a **Panaroo-compatible corrected binary matrix** for downstream analyses;
- a **multistate evidence matrix** for audit, plotting, and manuscript review.

The central design choice is conservative. Recore separates *sequence evidence*
from *clean binary presence*: full-length evidence can be corrected in the binary
matrix, while partial homology, pseudogenes, paralog traps, and
incompatible-context hits stay visible in audit outputs instead of silently
becoming presence calls.

<div class="grid cards" markdown>

-   :material-download: **[Install](getting-started/installation.md)**

    Build from source or take a tagged binary. Needs Rust and minimap2.

-   :material-rocket-launch: **[Quick start](getting-started/quickstart.md)**

    One command for a publication-grade correction pass, and what to read first.

-   :material-tune: **[Matrix policies](guide/matrix-policies.md)**

    Which evidence reaches the binary matrix, and how to word it in a paper.

-   :material-chart-box: **[Run report](guide/report.md)**

    The self-contained HTML viewer: phylogeny, matrix and synteny in one page.

-   :material-table: **[Output files](guide/outputs.md)**

    Every table Recore writes, column by column, with coordinate conventions.

-   :material-dna: **[How it works](how-it-works/algorithm.md)**

    Each pipeline phase, the data structures behind it, and the complexity.

</div>

## What a run gives you

```bash
recore \
  -p gene_presence_absence.csv \
  -r pan_genome_reference.fa \
  -g genomes/ \
  -a annotations/ \
  --synteny --multi-allele --nscluster \
  --matrix-policy balanced \
  --reference-genome AS00146 \
  --write-html-report \
  -t 8 -o recore_output/
```

| Read this first | For |
|-----------------|-----|
| `gene_presence_absence_recore.csv` | the corrected binary matrix, for downstream tools |
| `gene_presence_absence_recore_states.csv` | audit and figures, one biological state per cell |
| `orthology_evidence.tsv` | publication wording and manual review |
| `recore_report.html` | a self-contained viewer for the whole run |
| `evolutionary_events.tsv` | acquisition, loss, IS interruption, pseudogene or paralog per cluster |

The full list, with every column, is in [Output files](guide/outputs.md).

## What Recore does not claim

Recore is not a pangenome caller and it does not replace one. It cannot invent
evidence that is absent from the assembly sequence, it does not decide that a
partial homologue is a functional gene, and it does not force an ancestral story
onto an ambiguous phylogenetic pattern. Everything it declines to correct is
still written down. See [Limits and comparisons](how-it-works/limits.md).

## Citation

> Ruiz-Rodriguez P, et al. (2026). Recore: post-hoc rescue and quality control
> of pangenome presence/absence matrices. *In preparation.*
