# Recore

**Post-hoc rescue of falsely absent genes in pangenome matrices.**

Pangenome tools (Panaroo, Roary, PPanGGOLiN) sometimes report genes as absent in a genome when they are actually present — inflating the accessory genome and shrinking the core. Recore detects and corrects these false negatives by searching each "absent" gene against the genome assembly with minimap2.

Written in Rust for speed and low memory. Parallelised across genomes with rayon.

## How it works

```
Panaroo output                    Genome assemblies
(gene_presence_absence.csv)       (FASTA files)
         │                                │
         └──────────┬─────────────────────┘
                    ▼
           For each gap (gene × genome):
           minimap2 representative seq → genome
                    │
            ┌───────┼───────┐
            ▼       ▼       ▼
         No hit   Hit in    Hit not in
         (true    another   any cluster
         absence) cluster   (dropped)
                  (SPLIT)   (RESCUE)
                    │
                    ▼
           Corrected presence/absence matrix
           + rescue report + statistics
```

## Installation

### From source

```bash
git clone https://github.com/PathoGenOmics-Lab/recore.git
cd recore
cargo build --release
# Binary: target/release/recore
```

**Requirements:** Rust ≥1.70, [minimap2](https://github.com/lh3/minimap2) in PATH

### Pre-built binaries

See [Releases](https://github.com/PathoGenOmics-Lab/recore/releases).

## Usage

```bash
recore \
  -p gene_presence_absence.csv \
  -r pan_genome_reference.fa \
  -g genomes/ \
  -o recore_output/
```

### Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min-identity` | 90.0 | Minimum % identity for a hit |
| `--min-coverage` | 80.0 | Minimum % query coverage |
| `-t, --threads` | 4 | Parallelism (one minimap2 per genome) |

### With curated annotations (optional validation)

```bash
recore \
  -p gene_presence_absence.csv \
  -r pan_genome_reference.fa \
  -g genomes/ \
  -a annotations/ \
  -o recore_output/
```

When `-a` is provided (directory of GFF files), Recore validates each rescue against curated annotations: does the rescued region overlap a known CDS?

## Output

```
recore_output/
├── gene_presence_absence_recore.csv   # Corrected matrix
├── rescue_report.tsv                  # Per-rescue details
└── summary.txt                        # Statistics
```

### rescue_report.tsv

| Column | Description |
|--------|-------------|
| `gene_cluster` | Panaroo cluster name |
| `genome` | Genome where the gene was rescued |
| `rescue_type` | DROPPED or SPLIT:cluster_name |
| `identity` | % identity of the hit |
| `coverage` | % query coverage |
| `contig` | Target contig name |
| `start` | Start position on contig |
| `end` | End position on contig |
| `strand` | Strand (+/-) |
| `annotation_match` | Overlapping CDS from GFF (if -a provided) |

### summary.txt

```
Recore Summary
========================================
Clusters:           4200
Genomes:            480
Total gaps:         12450
Rescued:            312
Remaining gaps:     12138
Core size (before): 3421
Core size (after):  3587 (+166)

Rescue types:
  DROPPED: 245
  SPLIT:   67
```

## How does this differ from Panaroo's `--refind_prop_match`?

Panaroo has a built-in refind step, but:
- It uses a fixed proportion threshold, not coverage+identity
- It does not classify *why* a gene was missed (SPLIT vs DROPPED)
- It runs during the graph cleaning phase, where aggressive settings may re-remove refound genes
- It uses BLAST internally — slow for large datasets
- Recore runs **post-hoc**: it takes the final output and corrects it
- Recore uses minimap2 — orders of magnitude faster than BLAST
- Recore parallelises across genomes with rayon

## Citation

If you use Recore, please cite:

> Ruiz-Rodriguez P, et al. (2026). Recore: post-hoc rescue of falsely absent genes in pangenome matrices. *In preparation.*

## License

GNU General Public License v3.0
