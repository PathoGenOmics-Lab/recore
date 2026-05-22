# Parameters Reference

## Required

| Flag | Description |
|------|-------------|
| `-p, --presence-absence` | Panaroo `gene_presence_absence.csv` |
| `-r, --pan-reference` | Panaroo `pan_genome_reference.fa` |
| `-g, --genomes` | Directory of genome FASTA files (`.fasta`, `.fa`, `.fna`, also `.gz`) |

## Optional

| Flag | Default | Description |
|------|---------|-------------|
| `-a, --annotations` | — | Directory of GFF3 annotation files. Required for SPLIT/START_VARIANT classification and synteny validation. |
| `-o, --output` | `recore_output/` | Output directory |
| `-t, --threads` | `4` | Number of parallel threads (for rayon classification and minimap2) |

Recore matches FASTA and GFF files by exact matrix column name first, then by an
embedded `ASxxxxx` identifier. For single-contig genomes, GFF seqids are
normalized internally to the FASTA contig name so `--synteny` can work when a
GFF uses `contig_1` or `L4_AS00146` but the FASTA contig is `AS00146`.

## Alignment

| Flag | Default | Description |
|------|---------|-------------|
| `--preset` | `asm5` | minimap2 preset. `asm5` for sequences ≤5% divergence (MTBC, closely related strains). `asm20` for up to 20% divergence (more divergent species). |
| `--min-identity` | `90.0` | Minimum % nucleotide identity to consider a hit. For MTBC, the default (90%) is appropriate — real genes have >99%. Lowering this increases PARALOG noise. |
| `--min-coverage` | `80.0` | Minimum % of reference sequence that must be covered by a hit (single hit or merged fragments). |
| `--multi-allele` | `false` | Use all alleles per cluster as queries. More sensitive for clusters with high intra-cluster diversity, but slower. Panaroo produces multiple alleles in `pan_genome_reference.fa`. |

## Filtering

| Flag | Default | Description |
|------|---------|-------------|
| `--min-prevalence` | `0.1` | Minimum fraction of genomes where a cluster must be present to attempt rescue. Clusters present in fewer genomes are likely truly accessory and skipped. Set to `0.0` to attempt rescue of all gaps. |

## Rescue options

| Flag | Default | Description |
|------|---------|-------------|
| `--synteny` | `false` | Validate each rescue by comparing flanking gene context to reference genomes. Requires `--annotations`. Boosts synteny-consistent rescues (×1.2) and rejects explicit synteny failures from the corrected matrix while writing them to `synteny_rejected_report.tsv`. Also enables the exact missing-CDS pass for unannotated loci. Independently, PARALOG candidates without explicit synteny support are rejected to `paralog_rejected_report.tsv`. |
| `--write-corrected-gff` | `false` | Write per-genome `corrected_gff/*.recore.gff3` files with rescued CDS features appended. Original GFF records are preserved when `--annotations` is provided. |
| `--matrix-policy` | `balanced` | Controls which accepted rescues enter the binary corrected matrix. `strict` corrects only `STRICT_PRESENT` and `ALIAS_CLUSTER`; `balanced` also corrects full-length `SEQUENCE_ONLY`; `sensitive` also corrects partial/disrupted accepted evidence for exploratory analyses. |
| `--reference-genome` | — | Audit one genome column/reference background by exact matrix name or embedded `ASxxxxx` identifier. Runs a reference-only shadow search and writes `reference_genome_audit.tsv` with `new_to_reference`, `reference_like_or_present`, and `reference_like_ambiguous` classifications, plus `publication_confidence` labels for clean, complex/paralogous, fragmented, pseudogene, and unresolved reference evidence. |

## Nucleotide Similarity Clustering

| Flag | Default | Description |
|------|---------|-------------|
| `--nscluster` | `false` | Enable NSC. Detects and merges Panaroo clusters that are nearly identical at the DNA level but split at the protein level. |
| `--nsc-engine` | `minhash` | Clustering engine: `minhash` (built-in, approximate, no dependencies) or `mmseqs2` (exact, requires MMseqs2 in PATH). |
| `--nsc-kmer-size` | `21` | k-mer size for MinHash sketching. Larger values are more specific (fewer false positives). |
| `--nsc-threshold` | `0.8` | Maximum containment threshold for MinHash. Pairs above this are merged. |
| `--nsc-min-seq-id` | `0.9` | Minimum sequence identity for MMseqs2 linclust (0.0–1.0). |
| `--nsc-coverage` | `0.8` | Minimum bidirectional coverage for MMseqs2 linclust (0.0–1.0). |

## Logging

| Flag | Default | Description |
|------|---------|-------------|
| `-v, --verbose` | `false` | Enable DEBUG-level log output. Shows GFF loading details, synteny decisions (including why individual rescues pass/fail), per-genome missing FASTA warnings, and the top 10 rescues. |

---

## Parameter Guide by Organism

### *M. tuberculosis* complex (MTBC)

MTBC has >99.95% average nucleotide identity across strains. This has two consequences:

1. Real genes have very high identity — PARALOGs (90–94%) are nearly always noise
2. PE/PPE genes (~170 families with high mutual similarity) are a major source of cross-matches

**Recommended command:**

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
  --min-prevalence 0.0 \
  -t 8
```

**Why `--min-prevalence 0.0`:** In MTBC, even rare accessory clusters can have real false absences due to annotation differences. `0.0` ensures nothing is skipped.

**Why `asm5`:** Strains are <5% divergent. `asm5` is faster and more accurate at this range.

**Why `mmseqs2` for NSC:** More precise than MinHash for detecting the redundant annotation clusters common in MTBC pangenome analyses.

---

### More divergent bacteria / pan-genome with 95% ANI

```bash
recore \
  -p gene_presence_absence.csv \
  -r pan_genome_reference.fa \
  -g genomes/ \
  -a annotations/ \
  --synteny \
  --multi-allele \
  --nscluster \
  --preset asm20 \
  --min-prevalence 0.05 \
  -t 8
```

**Why `asm20`:** Allows up to 20% divergence — needed when comparing strains from the same species but different lineages.

**Why `--min-prevalence 0.05`:** Rare accessory clusters (present in <5% of genomes) generate excessive cross-match noise across the matrix.

---

## Threshold Effects

### `--min-identity` (alignment filter)

Affects PARALOG counts most strongly. If you see an unexpectedly large number of PARALOGs:

- Check if `--min-identity` was set below 90% by mistake
- For MTBC: real genes should be >99% identical; PARALOGs at 90–94% are noise
- Increasing to 95% eliminates most PE/PPE cross-matches but may miss a few rare variants

### `--min-coverage` (coverage filter)

Controls the boundary between PSEUDOGENE:PartialDeletion (60–80%) and no rescue (<60%):

- Increasing to 90%: only rescues nearly complete genes; misses short pseudogenes
- Decreasing to 60%: rescues more partial genes but increases false positives

### `--min-prevalence` (gap filter)

Affects total rescues and PARALOG rate:

- `0.0`: search all gaps (maximum sensitivity, higher PARALOG noise)
- `0.05`: skip clusters present in <5% of genomes
- `0.1` (default): skip clusters present in <10% of genomes

For core/soft-core analysis, `--min-prevalence 0.0` is recommended. For reducing runtime on large datasets, `0.05` is a reasonable compromise.

### `--nsc-threshold` (NSC similarity)

Controls how aggressively redundant clusters are merged:

- `0.8` (default): merge if 80% of one cluster's k-mers are contained in the other
- `0.9`: stricter, fewer merges
- `0.7`: more aggressive, may merge genuinely distinct paralogs

For MTBC, prefer MMseqs2 with `--nsc-min-seq-id 0.9` over MinHash, as it uses exact alignment identity.
