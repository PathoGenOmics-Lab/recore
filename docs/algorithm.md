# Recore Algorithm

## Overview

Recore rescues falsely absent genes from pangenome matrices produced by Panaroo (or Roary). It takes a gene presence/absence CSV, the pan-genome reference FASTA, individual genome FASTAs, and (optionally) GFF annotations, then identifies genes that are physically present in genomes but missing from the matrix.

The pipeline runs in 12 sequential phases. Each phase is deterministic given the same inputs and thread count.

---

## Phase 1: Input Parsing

### Presence/absence matrix

Parses Panaroo's `gene_presence_absence.csv`. The CSV has 13–15 metadata columns (Gene, Annotation, No. isolates, etc.) followed by one column per genome. Each cell is either empty (gene absent) or contains semicolon-separated gene IDs.

The matrix is stored as a `PresenceAbsence` struct with encapsulated access via newtype indices (`ClusterIdx`, `GenomeIdx`) to prevent index confusion.

### Pan-genome reference

Parses Panaroo's `pan_genome_reference.fa`. Each entry has format `>clusterName~~~genomeName description`. The `~~~` separator is stripped to extract the cluster name.

**Single-allele mode** (default): keeps only the longest allele per cluster.  
**Multi-allele mode** (`--multi-allele`): keeps all alleles, increasing sensitivity at the cost of more minimap2 queries.

### Genome FASTAs

Scans the genomes directory for files with extensions `.fasta`, `.fa`, `.fna` (also `.gz` variants). File stems become genome names matched against CSV column headers.

### GFF annotations (optional)

Parses GFF3 files matching genome names. Extracts CDS features only: contig, start, end, strand, gene ID. Coordinates are converted from GFF3's 1-based closed to 0-based half-open (matching minimap2 PAF output). Supports `.gff`, `.gff3`, and gzipped variants.

### Cross-validation

Before proceeding, validates:
- Genome and cluster names are safe for path construction (no `../`, `/`, `\0`)
- Warns about clusters in CSV but not in reference (and vice versa)
- Warns about genomes without matching FASTA files

---

## Phase 2: Gap Identification

Iterates over every (cluster, genome) pair. A **gap** is any pair where the gene is absent.

Gaps are filtered by **prevalence**: if a cluster is present in fewer than `--min-prevalence` fraction of genomes (default 10%), its gaps are skipped. The rationale is that ultra-rare clusters are likely truly accessory, and searching for them in hundreds of genomes generates cross-match noise.

The remaining gaps are grouped by genome for batch processing.

---

## Phase 3: GFF Index Construction

For each genome with a GFF file, builds a `GffIndex`:

1. Genes sorted by (contig, start position)
2. Per-contig range index for O(log n) overlap queries via binary search (`partition_point`)
3. Gene ID index for O(1) lookup by name

The index supports:
- `overlapping(contig, start, end)` → all CDS features overlapping a query region
- `flanking(contig, start, end)` → nearest upstream and downstream genes
- `get_by_id(gene_id)` → direct lookup

---

## Phase 4: Gene-to-Cluster Maps

For each genome involved in the search, builds a HashMap from gene ID → ClusterIdx. This maps every gene ID found in the CSV cells to its Panaroo cluster assignment.

Panaroo appends suffixes to gene IDs that don't exist in the original GFF:
- `_pseudo` for pseudogene annotations
- `_dup_N` for duplicate gene copies

The map indexes both the original ID and stripped versions to maximize GFF matching success.

---

## Phase 5: Batch Search (minimap2)

This is the computational core. Instead of running minimap2 once per genome (which would be O(n) process spawns), Recore uses batched minimap2 searches over concatenated targets.

### 5a. Query preparation

All unique cluster sequences are written to one temporary FASTA. In multi-allele mode, alleles are named `clusterName~~~ALLELE~0`, `clusterName~~~ALLELE~1`, etc.

### 5b. Target concatenation

All genome FASTAs are concatenated into one target file. Contig names are prefixed with `genomeName|||originalContig` to prevent collisions between genomes. The `|||` separator is chosen because it never appears in real contig names.

### 5c. Indexing

minimap2 builds an index (`.mmi`) of the concatenated target. This amortizes the O(n) indexing cost across all queries.

### 5d. Alignment

The primary minimap2 pass runs with:
- `--secondary=yes -N 50 -p 0.5`: report secondary alignments (needed for split/paralog detection)
- Preset `asm5` (default) or `asm20`: controls alignment sensitivity
- Streamed PAF output (stdout piped, no intermediate file)
- Stderr captured in a background thread to prevent deadlock

### 5e. Hit filtering

PAF lines are parsed in streaming fashion (no full output in RAM). Each hit is:
1. Demuxed back to its genome by parsing the `|||` prefix
2. Checked against the gap set — hits for non-gap (cluster, genome) pairs are discarded
3. Filtered by minimum identity (`--min-identity`, default 90%)

Result: a HashMap from `(genome_idx, cluster_name)` → `Vec<Hit>`.

### 5f. Exact missing-CDS pass

When `--synteny` and `--annotations` are enabled, Recore runs an additional exact-base alignment pass with minimap2 `-c` in small genome chunks. This pass addresses a specific failure mode in repetitive gene families: a true unannotated CDS-locus hit can be hidden when annotated paralogs consume minimap2's secondary-hit budget in a large concatenated target.

Only candidates classified as `MISSING_CDS` are imported from this second pass. All other extra hits are ignored, and a `MISSING_CDS` candidate is corrected in the binary matrix only if later synteny validation confirms the locus.

---

## Phase 6: Hit Classification

Each (cluster, genome) gap with hits is classified into one of 9 rescue types. Classification follows three paths:

### Path 1: Single full-coverage hit (≥ min_coverage)

The best hit (by query length, then score) covers ≥80% of the reference sequence (default `--min-coverage`). Classification depends on GFF annotation overlap:

| Condition | Type | Confidence |
|-----------|------|------------|
| No GFF provided | DROPPED | identity × coverage |
| GFF loaded, no substantial CDS overlaps hit, identity ≥95% | MISSING_CDS | identity × coverage |
| No GFF gene overlaps hit | DROPPED | identity × coverage |
| Overlapping gene belongs to a different cluster, identity ≥ 95% | SPLIT | base × 1.1 |
| Overlapping gene belongs to a different cluster, identity < 95% | PARALOG | base × 0.7 |
| Overlapping gene not in any cluster | DROPPED | identity × coverage |

`MISSING_CDS` is used only when the alignment is full-length and no CDS overlaps more than half of both the hit span and the annotated CDS. Short boundary overlaps with adjacent genes are ignored so a real missing annotation is not swallowed by a neighboring feature.

After initial classification, DROPPED hits undergo **start variant detection**: if the hit has ≥98% identity and a GFF gene within 300bp length difference overlaps ≥80% of the hit, it's reclassified as START_VARIANT.

### Path 2: Conservative fragment-chain assembly (same allele, contig, and strand)

If no single hit has sufficient coverage, hits are first grouped by query allele, query length, target contig, and strand. Recore builds the best co-linear chain within each group and uses only the best chain. This prevents independent alleles or paralogous hits from being merged into a false fragmented gene.

If the best chain reaches `--min-coverage` with ≥2 fragments:
- Non-co-linear chains are not rescued as fragments.
- If the whole chain lies inside one compatible CDS, Recore classifies it through the single-CDS annotation logic (DROPPED, START_VARIANT, SPLIT, or PARALOG) instead of calling FRAGMENTED.
- If average identity is ≥95% and a small gap is not divisible by 3, the event is classified as PSEUDOGENE:Frameshift.
- If the chain has a large query or target gap (≥50 bp), the event is classified as FRAGMENTED.
- Otherwise, the event is classified as ALIGNMENT_SPLIT: the gene is reconstructed by multiple alignments, but the evidence does not support a biological break.

Fragment-chain evidence is written to `rescue_report.tsv` for auditability.

### Path 3: Partial deletion pseudogene

If the best hit has ≥90% identity but only 60–80% coverage (below min_coverage), the gene is classified as PSEUDOGENE:PartialDeletion. Multi-fragment coverage is also considered.

### If none of the paths match → no rescue (gene is truly absent)

---

## Phase 7: Internal Stop Codon Scanning

DROPPED and MISSING_CDS rescues are post-processed for pseudogene detection. For each genome with at least one sequence-only rescue:

1. Load genome FASTA contigs on-demand (not preloaded for all genomes)
2. Extract the hit region from the contig
3. If the hit is on the minus strand, reverse-complement the sequence
4. Scan all 3 reading frames for internal stop codons (TAA, TAG, TGA)
5. Report the frame with the fewest stops (the "real" reading frame)
6. If stops found → reclassify as PSEUDOGENE:InternalStop, confidence × 0.5

---

## Phase 8: Synteny Validation

When `--synteny` is enabled and GFF files are loaded, each rescue is validated against the genomic context of reference genomes where the gene IS present.

### Algorithm

For the rescued gene at position (contig, start, end) in the target genome:
1. Build a small side-specific context window around the rescued locus: up to two upstream and two downstream CDSs mapped to Panaroo clusters. If the rescue substantially overlaps an annotated CDS, that CDS is used as the centre so short start-site differences and overlapping CDS boundaries do not shift the context. Otherwise the locus midpoint is used, which is important for missing-CDS calls with only tiny overlaps into adjacent genes.
2. For every reference genome where this cluster IS present, in deterministic genome order:
   a. Find the reference gene in the GFF (try exact ID, then strip `_pseudo`/`_dup_N` suffixes)
   b. Build the same side-specific context window around the reference gene
   c. Count recurrent reference contexts
3. If three or fewer reference contexts are available, any compatible reference context can validate the rescue. If more than three are available, a compatible context must be observed at least twice. This prevents a single noisy/paralogous context from validating a rescue in a large dataset.
4. If both sides are comparable, require at least one upstream and one downstream match on the same side. If only one side is comparable, require a same-side match on that side.
5. Otherwise → FAIL

Explicit FAIL candidates are written to `synteny_rejected_report.tsv` and are not used to correct the output matrix. Rescues with no mappable flanking context remain in `rescue_report.tsv` with `synteny` empty.

PARALOG candidates are more conservative: if they do not have explicit `synteny=PASS`, they are written to `paralog_rejected_report.tsv` and are not inserted into the corrected matrix.

MISSING_CDS candidates follow the same conservative matrix rule: `synteny=PASS` upgrades them to `STRICT_PRESENT`; without synteny support they are kept as `MISSING_CDS_UNVALIDATED` audit evidence and are not corrected by any matrix policy.

### Confidence adjustment

- Synteny PASS: confidence × 1.2 (boost)
- Synteny FAIL: confidence × 0.6, then reject from the corrected matrix
- No flanking genes found: no adjustment (synteny_ok = None)

---

## Phase 9: Conflict Resolution

Multiple rescues can overlap on the genome (e.g., a DROPPED and a PARALOG rescue at the same locus). Conflicts are resolved by:

1. Sort rescues by (genome, contig, start position)
2. For substantial overlaps on the same contig (≥50% of the shorter interval): keep the one with higher confidence
3. Single-pass O(n) sweep with a removal bitmap

---

## Phase 10: Nucleotide Similarity Clustering (NSC)

Optional (`--nscluster`). Detects Panaroo clusters that are distinct at the protein level but nearly identical at the nucleotide level (annotation-level redundancy).

### MinHash engine (default, built-in)

1. **Sketch**: Bottom-k MinHash with 128 hashes per cluster. Uses canonical k-mers (min of forward and reverse-complement hash) with FNV-1a hashing. Multi-allele clusters are concatenated with N-spacers to prevent cross-allele k-mers.

2. **LSH banding**: 16 bands of 8 rows each. Two clusters are candidates if they share any identical band hash. This reduces O(n²) comparisons to O(n × candidates).

3. **Containment**: For candidate pairs, compute max Jaccard containment (max of containment(A,B) and containment(B,A)). This is more appropriate than Jaccard similarity when genes differ in length.

4. **Clustering**: Union-Find with path compression and union-by-rank. Pairs above threshold (default 0.8) are merged.

5. **Representative selection**: Within each merged group, the cluster with highest prevalence becomes the representative. Others are removed from the output, and their presences are unioned into the representative.

### MMseqs2 engine (optional, exact)

Uses MMseqs2 `linclust` for linear-time clustering with alignment-based identity:
- `createdb` → `linclust` → `createtsv`
- Connected-component clustering mode
- Bidirectional coverage requirement
- Automatic temp directory cleanup

---

## Phase 11: Output Generation

### rescue_report.tsv

One row per rescue event. Columns: gene_cluster, genome, rescue_type, identity, coverage, contig, start (1-based), end, strand, confidence, synteny, annotation_match, plus fragment-chain evidence columns for multi-hit rescues.

### synteny_rejected_report.tsv

Written only when `--synteny` rejects at least one candidate. It uses the same schema as `rescue_report.tsv`, but these rows are audit records only and are not inserted into `gene_presence_absence_recore.csv`.

### paralog_rejected_report.tsv

Written only when PARALOG candidates lack explicit synteny support. It uses the same schema as `rescue_report.tsv`, but these rows are audit records only and are not inserted into `gene_presence_absence_recore.csv`.

### orthology_evidence.tsv

One row per accepted rescue plus rejected synteny or paralog candidate. This table makes the decision layer explicit: accepted/rejected status, orthology status, evidence state (`STRICT_PRESENT`, `ALIAS_CLUSTER`, `SEQUENCE_ONLY`, `PARTIAL_HOMOLOGY`, `PSEUDOGENE`, `HOMOLOG_PARALOG`, `REJECTED_CONTEXT`), matrix action, publication status, annotation-matched cluster, synteny state, confidence, coordinates, and rejection reason.

### reference_genome_audit.tsv

Written when `--reference-genome` is provided. The reference can be an exact matrix column or an embedded `ASxxxxx` identifier. For every cluster, Recore reports whether the reference matrix already contains it, whether Recore rescued reference evidence, whether only rejected/paralog/partial evidence exists, or whether no evidence was found. The `reference_absence_class` column is deliberately simple:
- `new_to_reference`: no reference evidence under the current thresholds.
- `reference_like_or_present`: already present or rescued in the audited reference.
- `reference_like_ambiguous`: partial, pseudogene, paralog, or rejected-context evidence; review before calling novel to the reference.

The reference audit also runs an exact reference-only shadow search and keeps a
shadow copy of candidates before conflict resolution. This prevents
partial/paralogous evidence from disappearing simply because it was correctly
excluded from the binary matrix. Such rows are reported as
`flag_reference_partial_paralog_evidence` with
`publication_confidence=complex_paralog`: the exact cluster can remain absent,
but the absence should not be described as a clean absence of the broader
sequence family.

The reference audit also compares reference-absent cluster representatives
against reference-present cluster representatives. Recore seeds candidate pairs
with canonical 15-mers, globally aligns the best candidates, and classifies
high-confidence matches as full-length same-gene sequence, short/long in-frame
alleles, frameshift-like length changes, or partial reference homology. These
events are written to `reference_cluster_events.tsv` and mirrored in
`reference_genome_audit.tsv` as `EVENT_DERIVED_REFERENCE_CLUSTER` with
`publication_confidence=event_derived`. They are deliberately audit-only:
Recore keeps the cluster as a separate matrix column but prevents downstream
tables from calling it a clean novel gene relative to the reference.

### cluster_locus_context.tsv / cluster_locus_summary.tsv

Written when GFF annotations are loaded. Recore builds a two-gene upstream and two-gene downstream context window for each present gene ID that can be resolved in the GFF. The context table lists every observed cluster-side context; the summary table flags clusters as stable single-locus, mostly stable with minor context noise, possible multi-locus/paralog, or insufficient context.

### gene_presence_absence_recore.csv

Drop-in replacement for Panaroo's CSV. Changes:
- Rescued genes: cell value `RESCUED:contig:start-end` (1-based coordinates)
- Binary correction is controlled by `--matrix-policy` (`strict`, `balanced`, `sensitive`)
- NSC-merged genes: cell value `NSC_MERGED`
- Clusters absorbed by NSC are removed (rows dropped)
- Metadata columns (No. isolates, No. sequences) updated to reflect rescues

### gene_presence_absence_recore_states.csv

Same shape as the corrected matrix, but genome cells are states rather than
binary calls: `PRESENT_INPUT`, `ABSENT`, `STRICT_PRESENT`, `ALIAS_CLUSTER`,
`SEQUENCE_ONLY`, `PARTIAL_HOMOLOGY`, `PSEUDOGENE`, `HOMOLOG_PARALOG`, or
`REJECTED_CONTEXT`. This is the preferred output for plotting or reviewing
ambiguous biology without forcing it into a clean 0/1 matrix.

### summary.txt

Human-readable summary with gap counts, accepted rescues, matrix-corrected
rescues under the selected policy, rescue type breakdown, core/soft-core genome
changes, synteny stats, and NSC stats.

### nscluster_report.tsv (if NSC enabled)

One row per NS cluster: cluster ID, representative name, member count, semicolon-separated member names.

---

## Phase 12: Core Genome Recalculation

After all rescues and NSC merges, the core and soft-core genome sizes are recalculated:

- **Core (strict)**: clusters present in ALL genomes (including rescues and NSC unions)
- **Soft-core (≥95%)**: clusters present in ≥95% of genomes

NSC-merged clusters are excluded from the count (their presences are absorbed into representatives).

---

## Rescue Type Reference

| Type | When | Confidence modifier |
|------|------|-------------------|
| DROPPED | Gene present in genome but cleaned by Panaroo graph cleaning | 1.0× (base) |
| MISSING_CDS | Gene sequence present at a locus with no substantial overlapping CDS in the GFF | 1.0× (base) |
| SPLIT | Gene assigned to wrong Panaroo cluster (≥95% identity to another cluster's gene) | 1.1× |
| FRAGMENTED | Gene broken into ≥2 fragments with combined coverage ≥80% and a large biological gap | 0.8× |
| ALIGNMENT_SPLIT | Gene reconstructed by ≥2 small co-linear alignments without biological break evidence | 0.9× |
| START_VARIANT | Same gene with different start codon (≤300bp length difference) | 0.95× |
| PARALOG | Cross-match to a related gene family (<95% identity) | 0.7× |
| PSEUDOGENE:InternalStop | Internal stop codon(s) in the best reading frame | 0.5× |
| PSEUDOGENE:Frameshift | Fragments with gaps not divisible by 3 | 0.5× |
| PSEUDOGENE:PartialDeletion | High identity (≥90%) but only 60–80% coverage | 0.4× |

---

## Computational Complexity

| Phase | Time | Space |
|-------|------|-------|
| Parsing | O(clusters × genomes) | O(clusters × genomes) |
| Gap filtering | O(clusters × genomes) | O(gaps) |
| GFF indexing | O(genes × log(genes)) per genome | O(genes) per genome |
| minimap2 | O(queries × target_size) | O(streaming) — PAF not buffered |
| Classification | O(hits) with rayon parallelism | O(rescues) |
| Internal stops | O(DROPPED_rescues × hit_length) | O(genome) per genome |
| Synteny | O(rescues × 5_references) | O(1) per rescue |
| Conflict resolution | O(rescues × log(rescues)) | O(rescues) |
| NSC MinHash | O(clusters × seq_length) + O(candidates) | O(clusters × sketch_size) |
| NSC MMseqs2 | O(clusters) linear time | O(mmseqs2 internal) |
| Output | O(clusters × genomes) | O(streaming) |
