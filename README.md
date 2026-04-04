# Recore

**Post-hoc rescue of falsely absent genes in pangenome matrices.**

Pangenome tools (Panaroo, Roary, PPanGGOLiN) sometimes report genes as absent when they are actually present — inflating the accessory genome and shrinking the core. Recore detects and corrects these false negatives by searching each "absent" gene against the genome assembly with minimap2.

Written in Rust for speed and low memory. Single-invocation batch minimap2 with pre-indexing and streamed PAF parsing.

## How it works

```
Panaroo output                    Genome assemblies       GFF annotations
(gene_presence_absence.csv        (FASTA directory)       (optional)
 + pan_genome_reference.fa)
         │                                │                      │
         └────────────────┬───────────────┘                      │
                          ▼                                      │
              Prevalence filter (skip truly                      │
              accessory genes, default ≥10%)                     │
                          │                                      │
                          ▼                                      │
              Batch minimap2 (all queries ×         ◄────────────┘
              all genomes in ONE invocation)         GFF-based classification
                          │
           ┌──────────┬───┼──────────┬──────────┐
           ▼          ▼   ▼          ▼          ▼
        DROPPED    SPLIT  FRAGMENTED START_     No hit
        (gene      (gene  (IS/frame- VARIANT   (true
        removed    in a   shift,     (alt.     absence)
        by graph   diff.  multiple   start
        cleaning)  clust) fragments) codon)
           │          │       │          │
           └──────────┴───────┴──────────┘
                          │
                    Conflict resolution
                   (overlapping rescues
                    → highest confidence wins)
                          │
                   Synteny validation
                   (optional, --synteny)
                          │
                          ▼
              Corrected presence/absence matrix
              + rescue report + summary statistics
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

### Basic

```bash
recore \
  -p gene_presence_absence.csv \
  -r pan_genome_reference.fa \
  -g genomes/ \
  -o recore_output/
```

### With annotations and synteny validation

```bash
recore \
  -p gene_presence_absence.csv \
  -r pan_genome_reference.fa \
  -g genomes/ \
  -a annotations/ \
  --synteny \
  --multi-allele \
  -t 8 \
  -o recore_output/
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-p, --presence-absence` | — | Panaroo `gene_presence_absence.csv` |
| `-r, --pan-reference` | — | Panaroo `pan_genome_reference.fa` |
| `-g, --genomes` | — | Directory of genome FASTA files |
| `-a, --annotations` | — | Directory of GFF files (enables rescue classification) |
| `-o, --output` | `recore_output/` | Output directory |
| `--min-identity` | 90.0 | Minimum % identity for a hit |
| `--min-coverage` | 80.0 | Minimum % query coverage (single hit or merged fragments) |
| `-t, --threads` | 4 | Number of threads |
| `--preset` | `asm5` | minimap2 preset (`asm5` for ≤5% divergence, `asm20` for ≤20%) |
| `--min-prevalence` | 0.1 | Skip clusters present in fewer genomes (truly accessory) |
| `--multi-allele` | off | Query with ALL alleles per cluster (slower, more sensitive) |
| `--synteny` | off | Validate rescues by flanking gene conservation (requires `-a`) |

## Rescue types

| Type | Description | Confidence |
|------|-------------|------------|
| **DROPPED** | Gene is present but was removed by Panaroo graph cleaning | Base |
| **SPLIT** | Gene locus is assigned to a different Panaroo cluster | Base × 1.1 |
| **FRAGMENTED** | Gene broken into ≥2 fragments (IS insertion, frameshift, assembly break) | Base × 0.8 |
| **START_VARIANT** | CDS present but with different start codon (length differs <300bp) | Identity × 0.95 |

## Output

```
recore_output/
├── gene_presence_absence_recore.csv   # Corrected matrix (Panaroo-compatible)
├── rescue_report.tsv                  # Per-rescue details
└── summary.txt                        # Statistics
```

### rescue_report.tsv

| Column | Description |
|--------|-------------|
| `gene_cluster` | Panaroo cluster name |
| `genome` | Genome where the gene was rescued |
| `rescue_type` | DROPPED, SPLIT:cluster, FRAGMENTED:Nfrags:cov%, START_VARIANT:±Nbp |
| `identity` | % identity of the minimap2 hit |
| `coverage` | % query coverage (merged for fragments) |
| `contig` | Target contig |
| `start` / `end` | Coordinates on contig |
| `strand` | +/- |
| `confidence` | Score for conflict resolution (0–1+) |
| `synteny` | PASS/FAIL (if `--synteny` enabled) |
| `annotation_match` | Overlapping CDS gene ID from GFF (if `-a` provided) |

### gene_presence_absence_recore.csv

Drop-in replacement for Panaroo's output. Rescued cells contain `RESCUED:contig:start-end`. Metadata columns (`No. isolates`, `No. sequences`, `Avg sequences per isolate`) are recalculated.

## Performance

Recore uses a single-invocation batch search strategy:

1. All genome FASTAs are concatenated with prefixed contig names
2. A minimap2 index (`.mmi`) is built once
3. All queries are searched in one minimap2 run with streamed PAF output
4. Hit classification is parallelised with rayon

This avoids spawning N separate minimap2 processes (one per genome), which would dominate runtime with process overhead and repeated indexing.

## How does this differ from Panaroo's `--refind_prop_match`?

| Feature | Panaroo refind | Recore |
|---------|---------------|--------|
| Timing | During graph cleaning | Post-hoc (on final output) |
| Aligner | BLAST | minimap2 (batch) |
| Classification | Binary (found/not found) | 4 types (DROPPED, SPLIT, FRAGMENTED, START_VARIANT) |
| Thresholds | Fixed proportion | Identity + coverage |
| Validation | None | Synteny + conflict resolution |
| Multi-allele | No | Yes (`--multi-allele`) |
| Speed | Slow (BLAST per genome) | Fast (single minimap2, streamed, parallel) |

## Citation

If you use Recore, please cite:

> Ruiz-Rodriguez P, et al. (2026). Recore: post-hoc rescue of falsely absent genes in pangenome matrices. *In preparation.*

## Contact

Paula Ruiz-Rodriguez — paula.ruiz.rodriguez@csic.es

## License

GNU General Public License v3.0
