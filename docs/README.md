# Recore Documentation

## Contents

| Document | Description |
|----------|-------------|
| [algorithm.md](algorithm.md) | Detailed description of every pipeline phase, data structures, and complexity |
| [rescue_types.md](rescue_types.md) | All 7 rescue types with examples and biological interpretation |
| [output_files.md](output_files.md) | Format of every output file, coordinate system, log messages |
| [parameters.md](parameters.md) | All CLI parameters with defaults, effects, and organism-specific guides |

## Quick reference

```
recore -p gene_presence_absence.csv \
       -r pan_genome_reference.fa \
       -g genomes/ \
       -a annotations/ \
       --synteny --multi-allele --nscluster \
       --preset asm5 --min-prevalence 0.0 \
       -t 8 -o recore_output/
```

Outputs: `rescue_report.tsv`, `gene_presence_absence_recore.csv`, `summary.txt`, `nscluster_report.tsv`
