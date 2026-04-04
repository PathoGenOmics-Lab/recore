# Output Files

All output files are written to the directory specified by `-o` (default: `recore_output/`).

---

## rescue_report.tsv

Primary output. One row per rescued gene. Tab-separated.

### Columns

| Column | Type | Description |
|--------|------|-------------|
| `gene_cluster` | string | Panaroo cluster name (from CSV `Gene` column) |
| `genome` | string | Genome name (matches CSV column header and FASTA filename stem) |
| `rescue_type` | string | See [rescue types](rescue_types.md) |
| `identity` | float | % nucleotide identity of the best hit (0–100) |
| `coverage` | float | % of reference sequence covered by hit(s) (0–100) |
| `contig` | string | Contig name in the genome FASTA |
| `start` | int | Start coordinate, **1-based** |
| `end` | int | End coordinate, 1-based inclusive |
| `strand` | char | `+` or `-` |
| `confidence` | float | Combined quality score (0–1.2+); see [algorithm](algorithm.md#confidence-score) |
| `synteny` | string | `PASS`, `FAIL`, or empty (not assessed) |
| `annotation_match` | string | GFF gene ID that overlapped the rescue region, if any |

### Example

```tsv
gene_cluster	genome	rescue_type	identity	coverage	contig	start	end	strand	confidence	synteny	annotation_match
group_2044	L1_AS00349	DROPPED	100.0	100.0	L1_AS00349	245301	246100	+	1.198	PASS	Rv0210c
group_0187	A4_AS00248	SPLIT:group_0188	99.2	99.8	A4_AS00248	1203401	1204300	+	1.087	PASS	Rv1099c
group_3940	L1_AS00366	PSEUDOGENE:2stops	99.0	100.0	L1_AS00366	1783748	1784254	-	0.495	PASS	Rv2651c
```

### Filtering suggestions

```bash
# High-confidence rescues only
awk '$10 > 0.8' rescue_report.tsv

# Synteny-validated only
awk '$11 == "PASS"' rescue_report.tsv

# Non-pseudogene rescues
grep -v "PSEUDOGENE" rescue_report.tsv

# DROPPED and SPLIT only (most confident true rescues)
grep -E "DROPPED|SPLIT" rescue_report.tsv | awk '$10 > 0.9'
```

---

## gene_presence_absence_recore.csv

Drop-in replacement for Panaroo's `gene_presence_absence.csv`. Identical format and column order.

### Changes from the original

1. **Rescued cells**: Empty cells (gene absent) become `RESCUED:contig:start-end` where coordinates are 1-based.
   ```
   Before: ""
   After:  "RESCUED:L1_AS00349:245301-246100"
   ```

2. **NSC-merged cells** (only if `--nscluster`): Cells where a merged cluster would have added presence become `NSC_MERGED`.

3. **Metadata columns updated**: `No. isolates`, `No. sequences`, and `Avg sequences per isolate` are updated to reflect the new presence count.

4. **NSC clusters removed**: Panaroo clusters that were merged by NSC into a representative cluster are dropped entirely (their rows are absent from the output).

### Compatibility

This file can be used directly as input to downstream tools that accept Panaroo's format:
- `distree` for phylogenetic distance analysis
- `fstic` for allele frequency-based FST calculations
- Custom R/Python scripts expecting Panaroo output

---

## summary.txt

Human-readable run summary. Example:

```
Recore Summary
==================================================
Clusters:              4103
Genomes:               479
Total gaps:            91509
Gaps searched:         91509
Skipped (prevalence):  0
Rescued:               947
Remaining gaps:        90562

Core genome (strict):  1539 → 1546 (+7)
Soft-core (≥95%):      3781 → 3776 (-5)

Rescue types:
  DROPPED:       378
  SPLIT:         206
  FRAGMENTED:    0
  START_VARIANT: 60
  PARALOG:       113
  PSEUDOGENE:    190

Synteny validation:
  PASS:  747
  FAIL:  195

Nucleotide similarity clustering:
  NS clusters:      19
  Genes merged:     39
  Clusters removed: 20
```

### Interpreting core genome changes

- **Core +7**: 7 clusters that were missing from at least one genome are now found in all — true core genes that Panaroo had incorrectly excluded.
- **Soft-core -5**: 5 clusters that appeared to be soft-core were actually redundant with other clusters (merged by NSC). The true soft-core is slightly smaller.

---

## nscluster_report.tsv

Written only if `--nscluster` is enabled and at least one merge was found.

| Column | Description |
|--------|-------------|
| `ns_cluster` | NS cluster ID (NS_0, NS_1, ...) |
| `representative` | Name of the representative Panaroo cluster (highest prevalence) |
| `n_members` | Total number of clusters in the group (including representative) |
| `members` | Semicolon-separated list of all member cluster names |

### Example

```tsv
ns_cluster	representative	n_members	members
NS_0	group_2891	3	group_2891;group_2892;group_2893
NS_1	group_0441	2	group_0441;group_0442
```

These clusters had nearly identical nucleotide sequences but were split by Panaroo due to amino acid-level differences (typically from upstream annotation artifacts).

---

## Coordinate System

All coordinates in Recore output are **1-based** (standard for bioinformatics formats like GFF3, VCF, GenBank):

- `rescue_report.tsv`: `start` column is 1-based
- `gene_presence_absence_recore.csv`: `RESCUED:contig:start-end` coordinates are 1-based

minimap2 internally uses 0-based half-open coordinates (PAF format). Recore converts these to 1-based for all output.

---

## Log Output

Recore writes structured log output to stderr. Use `-v` / `--verbose` for debug-level detail.

### Standard output (INFO level)

```
[INFO recore::pipeline] Parsed presence/absence: 4103 clusters × 479 genomes (core: 1539, soft-core 95%: 3781)
[INFO recore::pipeline] Gaps: 91509 total, 91509 after prevalence filter (≥0%), 0 skipped
[INFO recore::pipeline] Core genome (before): 1539 / 4103 clusters
[INFO recore::pipeline] Loaded GFF annotations for 479 genomes
[INFO recore::pipeline] Prepared batch search: 4103 queries × 479 genomes
[INFO recore::pipeline] Built minimap2 index
[INFO recore::pipeline] Validating synteny...
[INFO recore::pipeline] Synteny: 747 validated, 195 failed, 5 skipped (no flanking genes found)
[INFO recore::pipeline] Final rescues: 947 (DROPPED: 378, SPLIT: 206, ...)
[INFO recore::pipeline] Core genome: 1539 → 1546 (+7)
```

### Verbose output (DEBUG level, with -v)

```
[DEBUG recore::rescue] Synteny skip: no flanking genes mapped to clusters for rescue cluster=3940 genome=75
[DEBUG recore::rescue] Upstream gene 'L1_AS00366_Rv2650c' not found in any Panaroo cluster
[DEBUG recore::pipeline] Top rescues (first 10):
[DEBUG recore::pipeline]   group_2044 in A1_AS00349 → Dropped (id=100.0%, cov=100.0%, conf=1.20, synteny=Some(true))
```

### Common debug messages and their meaning

| Message | Meaning |
|---------|---------|
| `ref gene 'xxx_pseudo' not found in GFF of genome idx=N` | Panaroo refound gene ID has no counterpart in GFF (expected for `refound_*` genes) |
| `Upstream gene 'Rv2650c' not found in any Panaroo cluster` | Flanking gene was cleaned from the pangenome; synteny cannot be assessed |
| `Synteny skip: no flanking genes mapped` | Neither upstream nor downstream gene has a cluster assignment; rescue is kept but synteny_ok = None |
| `X genomes with gaps have no FASTA file (skipped)` | Genome names in CSV don't match FASTA filenames; check naming |
