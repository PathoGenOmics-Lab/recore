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
| `synteny` | string | `PASS` or empty (not assessed). `FAIL` rows are written to `synteny_rejected_report.tsv` when `--synteny` is enabled |
| `annotation_match` | string | GFF gene ID that overlapped the rescue region, if any |
| `fragment_n` | int | Number of fragments in the selected multi-hit chain, if applicable |
| `fragment_query_gap_bp` | int | Total query-space gap across the selected fragment chain |
| `fragment_target_gap_bp` | int | Total target-space gap across the selected fragment chain |
| `fragment_colinear` | bool | Whether fragments are co-linear on the target |
| `fragment_same_cds` | bool | Whether the chain is contained in one compatible CDS |
| `fragment_note` | string | Fragment interpretation, such as `alignment_split`, `biological_gap`, `frameshift_gap`, or `same_cds_alignment_split` |

### Example

```tsv
gene_cluster	genome	rescue_type	identity	coverage	contig	start	end	strand	confidence	synteny	annotation_match	fragment_n	fragment_query_gap_bp	fragment_target_gap_bp	fragment_colinear	fragment_same_cds	fragment_note
group_2044	L1_AS00349	DROPPED	100.0	100.0	L1_AS00349	245301	246100	+	1.198	PASS	Rv0210c	NA	NA	NA	NA	NA	NA
group_0187	A4_AS00248	SPLIT:group_0188	99.2	99.8	A4_AS00248	1203401	1204300	+	1.087	PASS	Rv1099c	NA	NA	NA	NA	NA	NA
group_3940	L1_AS00366	PSEUDOGENE:2stops	99.0	100.0	L1_AS00366	1783748	1784254	-	0.495	PASS	Rv2651c	NA	NA	NA	NA	NA	NA
group_2101	L2_AS00401	ALIGNMENT_SPLIT:2frags:98%cov	99.6	98.4	L2_AS00401	552100	553020	+	0.882	PASS		2	9	9	true	false	alignment_split
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

## synteny_rejected_report.tsv

Written only when `--synteny` is enabled and one or more candidate rescues have
explicitly incompatible flanking context. The schema is identical to
`rescue_report.tsv`, but these rows are audit records: they are not added to
`gene_presence_absence_recore.csv` and are not written to corrected GFF output.

This file is useful for reviewing high-identity paralog/cross-family hits that
look like sequence rescues but fail the local genomic-context check.

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

## cluster_relationships.tsv

Written for every run. The file may contain only the header when no inter-cluster relationships are detected.

This table summarizes reproducible evidence that one Panaroo cluster is related
to another. It is intended as an audit layer: Recore reports relationships and
recommendations, but only nucleotide similarity clustering (`--nscluster`) changes
the output matrix automatically.

### Columns

| Column | Description |
|--------|-------------|
| `query_cluster` | Cluster being rescued or reviewed |
| `related_cluster` | Cluster linked by SPLIT/PARALOG or NSC evidence |
| `relationship_type` | `split`, `paralog`, `annotation_match_cluster`, or `nucleotide_similarity_cluster` |
| `evidence_type` | Specific rescue type or `nscluster` |
| `n_genomes_supporting` | Number of genomes supporting the relationship |
| `synteny_support` | `pass`, `fail`, `mixed`, or `not_tested` |
| `identity_*` | Minimum, median and maximum nucleotide identity across supporting rescues |
| `coverage_*` | Minimum, median and maximum query coverage across supporting rescues |
| `recommendation` | Conservative interpretation: `collapse_alias`, `review_alias_without_full_synteny`, `review_split`, `review_annotation_match`, or `review_paralog_keep_distinct_by_default` |
| `example_regions` | Up to five genome/region examples |

### Interpretation

- `collapse_alias`: near-complete, high-identity, synteny-supported evidence that the query is not an independent cluster.
- `review_alias_without_full_synteny`: high sequence support, but synteny is absent, mixed, or not tested.
- `review_annotation_match`: the rescued locus overlaps a feature already assigned to another cluster; inspect as a possible naming or boundary issue.
- `review_paralog_keep_distinct_by_default`: related sequence family member; keep distinct unless additional evidence supports collapse.

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
Rejected (synteny):    4
Remaining gaps:        90562

Core genome (strict):  1539 → 1546 (+7)
Soft-core (≥95%):      3781 → 3776 (-5)

Rescue types:
  DROPPED:       378
  SPLIT:         206
  FRAGMENTED:    0
  ALIGNMENT_SPLIT: 12
  START_VARIANT: 60
  PARALOG:       113
  PSEUDOGENE:    190

Synteny validation:
  PASS:  747
  FAIL:  0
  REJECTED_FAIL: 4

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

## corrected_gff/*.recore.gff3

Written only when `--write-corrected-gff` is enabled.

For each genome with at least one rescue, Recore writes a GFF3 file under
`corrected_gff/`. If an original GFF was supplied with `--annotations`, its
records are preserved and Recore appends rescued CDS features. If no original
GFF is available for a rescued genome, the file contains only the appended
rescue features.

Appended rescue features use source `Recore`, type `CDS`, and include:

| Attribute | Description |
|-----------|-------------|
| `gene_cluster` | Panaroo/Recore cluster name |
| `rescue_type` | Rescue classification |
| `identity` | Nucleotide identity (%) |
| `coverage` | Query coverage (%) |
| `confidence` | Recore confidence score |
| `synteny` | `PASS`, `FAIL`, or `NA` |
| `annotation_match` | Original GFF feature ID overlapped by the rescue, if any |
| `fragment_n` | Number of fragments in the selected multi-hit chain, if applicable |
| `fragment_query_gap_bp` | Total query-space gap across the fragment chain |
| `fragment_target_gap_bp` | Total target-space gap across the fragment chain |
| `fragment_note` | Fragment interpretation, if applicable |

These GFFs are designed as auditable corrected annotations. If a rescue overlaps
an existing feature, the original feature is preserved and the appended Recore
feature records the rescue evidence explicitly.

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
[INFO recore::pipeline] Synteny rejection: 947 → 752 rescues (195 explicit FAIL candidates written to synteny_rejected_report.tsv)
[INFO recore::pipeline] Final rescues: 752 (DROPPED: 378, SPLIT: 206, ...)
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
