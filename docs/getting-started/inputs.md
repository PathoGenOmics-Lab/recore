# Inputs

| Input | Required | Description |
|-------|----------|-------------|
| `gene_presence_absence.csv` | yes | Panaroo-style matrix: metadata columns followed by one column per genome. Roary, PPanGGOLiN and PanTA write this layout, and PIRATE's own table is read too (see below) |
| `pan_genome_reference.fa` | yes | Cluster representative sequences; FASTA IDs must match matrix cluster names |
| genome FASTAs | yes | One assembly per matrix genome column |
| GFF annotations | optional, recommended | Enables annotation-aware `SPLIT`, `START_VARIANT`, `PARALOG`, pseudogene and synteny classification |
| reference genome name | optional | Used by `--reference-genome` to audit whether clusters are absent from a chosen reference |
| phylogenetic tree | optional | Newick, for `--tree`: parsimony on the corrected state, and the phylogeny pane of the HTML report |
| mobile-element FASTA | optional | For `--mobile-element-reference`: sequence-based IS support instead of product names |

## Which files each pangenome tool gives you

Pass `-p` the matrix and `-r` the representative sequences, and `-a` the same
GFF files the tool was run on, so the gene IDs in the matrix cells are the ones
in the annotations. Each pair was checked on real output of the version shown.

| Tool | `-p` | `-r` |
|------|------|------|
| Panaroo | `gene_presence_absence.csv` | `pan_genome_reference.fa` |
| Roary | `gene_presence_absence.csv` | `pan_genome_reference.fa` |
| PPanGGOLiN (2.3.2) | `matrix.csv`, from `ppanggolin all` or `ppanggolin write_pangenome --csv` | `all_nucleotide_families.fasta`, from `ppanggolin fasta --gene_families all` |
| PanTA (1.0.1) | `gene_presence_absence.csv` | `representative_clusters_nucl.fasta` |
| PIRATE (1.0.5) | `PIRATE.gene_families.tsv` | `representative_sequences.ffn` |

For PIRATE, give recore its own table, not the output of `PIRATE_to_roary.pl`:
that script names every row `allele__family`, which matches no sequence in
PIRATE's representative FASTA. The table is read into Panaroo's layout, so the
corrected files are the same as for a Panaroo run, and a fission written
`(locus1:locus2)` counts its pieces as one gene. PanTA writes its gene IDs as
`sample-contig-gene`; recore matches them to the GFF by the gene part.

## Name matching

This is where most failed runs come from.

- Matrix genome columns should match FASTA and GFF file stems, or contain the
  same embedded `ASxxxxx` identifier.
- Pan-reference FASTA cluster names must match the `Gene` column of the matrix.
- GFF3 attributes are percent-decoded before matching, so `ID`, `Name`, `gene`,
  `locus_tag`, aliases, products and notes all reach the matcher with their real
  values. An encoded locus tag no longer silently misses its matrix entry.
- Tree tips are matched to matrix genomes by name. Genomes with no matching tip
  are ignored, and Recore reports that explicitly rather than dropping them
  quietly.

## What is rejected before the pipeline starts

Recore validates the inputs before any expensive search, so a malformed run
fails in seconds rather than after an hour of minimap2:

- empty presence/absence matrices;
- empty pan-reference FASTA files;
- duplicate genome names;
- duplicate cluster names.

## Choosing the optional inputs

> **Use GFFs whenever you have them.** Without annotations Recore can still
> detect full-length sequence evidence and some partial or deletion evidence,
> but it cannot reliably distinguish `SPLIT`, `START_VARIANT` and `PARALOG` by
> annotation overlap, and synteny validation is not available at all.

- `--synteny`: Use it with GFFs whenever paralogy or multi-locus clusters are plausible. It
  is the check that keeps a paralogous hit from becoming a presence call.

- `--multi-allele`: Use it when clusters have divergent alleles, so a diverged true allele is not
  mistaken for a different gene.

- `--min-prevalence`: `0.0` for exhaustive audits where every cluster matters. A higher value such
  as `0.05` cuts rare-cluster noise in large publication runs.

- `--preset`: `asm5` for assemblies within about 5% divergence, which is faster. `asm20`
  for up to 20%, which is safer across a broader species range.
