# Quick start

If you would rather see a full run on data with known answers before pointing
Recore at your own, start with the [tutorial](tutorial.md).

## The run most people want

A publication-grade correction pass, with FASTA assemblies and GFF annotations:

```bash
recore \
  -p gene_presence_absence.csv \
  -r pan_genome_reference.fa \
  -g genomes/ \
  -a annotations/ \
  --synteny \
  --multi-allele \
  --nscluster \
  --nsc-engine mmseqs2 \
  --nsc-min-seq-id 0.9 \
  --nsc-coverage 0.8 \
  --preset asm5 \
  --min-prevalence 0.05 \
  --matrix-policy balanced \
  --reference-genome AS00146 \
  --write-html-report \
  -t 8 \
  -o recore_output/
```

Every flag is described in [Parameters](../guide/parameters.md). The four that
change what a run means, rather than how fast it is:

`--matrix-policy`

:   Which accepted evidence reaches the binary matrix. `balanced` is the default
    recommendation. See [Matrix policies](../guide/matrix-policies.md).

`--synteny`

:   Validates candidate loci against their flanking genes. Turn it on whenever
    paralogy or multi-locus clusters are plausible, which is most of the time.

`--reference-genome`

:   Asks whether each cluster is genuinely absent from a chosen reference, or
    whether it looks like a reference-derived allele or fragment.

`--write-html-report`

:   Writes `recore_report.html`, a self-contained viewer for the whole run. See
    [Run report](../guide/report.md).

## What to read first

| Read this | For |
|-----------|-----|
| `gene_presence_absence_recore.csv` | the corrected binary matrix, for downstream tools |
| `gene_presence_absence_recore.Rtab` | the same matrix as 1/0, for pyseer and Rtab readers |
| `gene_presence_absence_recore_states.csv` | audit and figures, one biological state per cell |
| `orthology_evidence.tsv` | publication wording and manual review |
| `recore_report.html` | the whole run in one page, offline |
| `reference_genome_audit.tsv` | whether clusters are truly absent from a reference |
| `evolutionary_events.tsv` | what kind of event best explains each cluster |

## Other run profiles

=== "Evolutionary and mobile-element audit"

    When the question is not only "is the gene present?" but "what kind of event
    explains the signal?".

    ```bash
    recore \
      -p gene_presence_absence.csv \
      -r pan_genome_reference.fa \
      -g genomes/ \
      -a annotations/ \
      --synteny --multi-allele --nscluster \
      --matrix-policy balanced \
      --reference-genome AS00146 \
      --tree species.tree \
      --mobile-element-reference mobile_elements.fa \
      --mobile-element-min-identity 90 \
      --mobile-element-min-coverage 70 \
      --mobile-element-window 500 \
      -t 8 \
      -o recore_evolutionary_audit/
    ```

    `evolutionary_events.tsv` then keeps matrix correction, phylogenetic pattern,
    reference similarity, local context, NSC redundancy, and nearby
    mobile-element support in one row per cluster.

    `--tree` also unlocks the phylogeny pane of the HTML report.

=== "Without GFF annotations"

    ```bash
    recore \
      -p gene_presence_absence.csv \
      -r pan_genome_reference.fa \
      -g genomes/ \
      --matrix-policy balanced \
      -t 8 \
      -o recore_output/
    ```

    Full-length sequence evidence is still rescued, but without annotations
    Recore cannot reliably tell `SPLIT`, `START_VARIANT` and `PARALOG` apart by
    annotation overlap, and synteny validation is unavailable.

=== "Exhaustive reference audit"

    ```bash
    recore \
      -p gene_presence_absence.csv \
      -r pan_genome_reference.fa \
      -g genomes/ \
      -a annotations/ \
      --synteny --multi-allele \
      --reference-genome AS00146 \
      --min-prevalence 0.0 \
      --matrix-policy balanced \
      -t 8 \
      -o recore_reference_audit/
    ```

    `--min-prevalence 0.0` stops rare clusters being filtered out, which is what
    you want when the point is to classify every cluster against a reference
    background.

=== "Built-in NSC engine"

    Without MMseqs2 installed, use the built-in MinHash engine:

    ```bash
    recore \
      -p gene_presence_absence.csv \
      -r pan_genome_reference.fa \
      -g genomes/ \
      --nscluster \
      --nsc-engine minhash \
      -o recore_output/
    ```

## Reporting the run

The wording that matches what `balanced` actually does:

> We used Recore with the balanced matrix policy, which rescues full-length
> false absences while retaining partial or disrupted homologues as multistate
> audit-only evidence.

`orthology_evidence.tsv` carries per-cluster wording guidance for the cases that
need a sentence of their own.
