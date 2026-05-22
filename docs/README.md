# Recore Documentation

## Contents

| Document | Description |
|----------|-------------|
| [algorithm.md](algorithm.md) | Detailed description of every pipeline phase, data structures, and complexity |
| [rescue_types.md](rescue_types.md) | All 8 rescue types with examples and biological interpretation |
| [output_files.md](output_files.md) | Format of every output file, coordinate system, log messages |
| [parameters.md](parameters.md) | All CLI parameters with defaults, effects, and organism-specific guides |

## Quick reference

```
recore -p gene_presence_absence.csv \
       -r pan_genome_reference.fa \
       -g genomes/ \
       -a annotations/ \
       --synteny --multi-allele --nscluster \
       --matrix-policy balanced \
       --reference-genome AS00146 \
       --preset asm5 --min-prevalence 0.0 \
       -t 8 -o recore_output/
```

Outputs: `rescue_report.tsv`, `orthology_evidence.tsv`, `reference_genome_audit.tsv` when `--reference-genome` is used, `gene_presence_absence_recore.csv`, `gene_presence_absence_recore_states.csv`, `summary.txt`, `cluster_relationships.tsv`, `cluster_locus_context.tsv` / `cluster_locus_summary.tsv` when GFFs are loaded, `nscluster_report.tsv`, `synteny_rejected_report.tsv` when explicit synteny failures are found, and `paralog_rejected_report.tsv` when PARALOG candidates lack synteny support.
