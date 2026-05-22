# Rescue Types

Recore classifies each rescued gene into one of 8 types based on the alignment evidence and GFF annotation context.

Examples below show the core `rescue_report.tsv` columns; see [output_files.md](output_files.md) for the full schema, including fragment evidence columns.

---

## DROPPED

**What it means:** The gene is physically present in the genome but was removed by Panaroo during graph cleaning. Panaroo's graph-based cleaning can over-aggressively trim genes that are divergent or connected to fragmented paths in the de Bruijn graph.

**Evidence:** Full-coverage hit (≥80%) with no overlapping GFF gene assigned to a different cluster.

**Confidence:** `identity × coverage` (base confidence, unmodified)

**Example in rescue_report.tsv:**
```
geneX   genome1   DROPPED   99.8   100.0   contig_1   230401   231234   +   0.998   PASS   Rv0210c
```

---

## SPLIT

**What it means:** The gene locus is present and annotated, but Panaroo assigned it to a different cluster. This typically happens when a gene has subtle sequence variation that caused it to cluster separately from the canonical reference.

**Evidence:** Full-coverage hit (≥80%) with a GFF gene at that locus that belongs to a different Panaroo cluster, **with ≥95% nucleotide identity**.

**Confidence:** `identity × coverage × 1.1` (slight boost — evidence is strong)

**Example:**
```
geneY   genome2   SPLIT:group_1847   99.1   99.5   contig_1   845201   846100   +   1.08   PASS   Rv0756c
```

The value after `SPLIT:` is the cluster that actually contains the annotation.

---

## FRAGMENTED

**What it means:** The gene is broken into multiple biological pieces in this genome. Common causes are assembly breaks (contig boundaries), IS element insertions, or large indels that split the reading frame.

**Evidence:** A co-linear chain of ≥2 hits from the same query allele, target contig, and strand whose merged coverage reaches ≥80% of the reference, with a large query or target gap (≥50 bp). Non-co-linear chains are not rescued as FRAGMENTED, and chains contained inside one compatible CDS are classified through the single-CDS annotation logic instead.

**Confidence:** `avg_identity × total_coverage × 0.8`

**Example:**
```
geneZ   genome3   FRAGMENTED:2frags:87%cov   97.2   87.4   contig_3   112000   113400   +   0.672   FAIL   
```

The `start` and `end` coordinates span from the first to the last fragment.

---

## ALIGNMENT_SPLIT

**What it means:** The gene is reconstructed by multiple minimap2 alignments, but the evidence does not support a biological gene break. This is usually an alignment or representation split rather than a functional disruption.

**Evidence:** A co-linear multi-hit chain reaches the coverage threshold, but the query and target gaps are small (<50 bp), no frameshift-like gap is detected, and the chain is not better explained as a single-CDS annotation match.

**Confidence:** `avg_identity × total_coverage × 0.9`

**Example:**
```
geneZ   genome3   ALIGNMENT_SPLIT:2frags:98%cov   99.6   98.4   contig_3   112000   113400   +   0.882   PASS
```

These events are reported separately so downstream analyses can avoid treating technical split alignments as biological fragmentation.

---

## START_VARIANT

**What it means:** The gene is present with the correct sequence but starts at a different position than the reference annotation. Alternative start codons (GTG, TTG) and upstream ATGs can lead to start site disagreements between annotators.

**Evidence:** Full-coverage hit (≥80%) with ≥98% identity, overlapping a GFF gene of similar length (within ±300 bp), where the length difference is non-zero.

**Confidence:** `identity × coverage × 0.95`

**Example:**
```
geneA   genome4   START_VARIANT:+57bp   99.6   100.0   contig_1   34100   34956   +   0.941   PASS   Rv0033
```

`+57bp` means the reference annotation extends 57bp upstream of where this genome's annotator placed the start.

---

## PARALOG

**What it means:** The hit is a related gene from a different family (a paralog), not the missing gene itself. Common in high-GC organisms with large multi-copy gene families (PE/PPE in *M. tuberculosis*, IS elements).

**Evidence:** Full-coverage hit with a GFF gene at that locus belonging to a different cluster, **with <95% identity**.

**Confidence:** `identity × coverage × 0.7` (reduced — this is a cross-match, not a true rescue)

**Example:**
```
geneB   genome5   PARALOG:PE_PGRS_group   91.3   98.2   contig_2   1920400   1921800   +   0.501   FAIL   PE_PGRS47
```

> **Note for MTBC analyses:** A large fraction of PARALOGs in tuberculosis data come from the PE/PPE superfamily, which comprises ~170 genes with high nucleotide similarity. These are not rescued genes — they are correctly absent.

---

## PSEUDOGENE:InternalStop

**What it means:** The gene sequence is present with good coverage but contains one or more premature stop codons (TAA, TAG, TGA) in the reading frame. This is a functional pseudogene.

**Evidence:** Full-coverage hit reclassified after scanning the genome FASTA sequence. The region is extracted, reverse-complemented if on the minus strand, and scanned in all 3 reading frames. The frame with the fewest stops is reported.

**Confidence:** `identity × coverage × 0.5` (halved)

**Example:**
```
geneC   genome6   PSEUDOGENE:2stops   99.1   100.0   contig_1   2103400   2104200   -   0.496   PASS   Rv1884c
```

---

## PSEUDOGENE:Frameshift

**What it means:** The gene is fragmented in a pattern consistent with a frameshift mutation — fragments have high identity but are separated by small gaps whose size is not divisible by 3.

**Evidence:** ≥2 fragments with average identity ≥95%, gaps between fragments ≤50bp, and at least one gap whose size mod 3 ≠ 0.

**Confidence:** `avg_identity × total_coverage × 0.5`

**Example:**
```
geneD   genome7   PSEUDOGENE:frameshift   97.8   91.0   contig_2   654000   655100   +   0.446   PASS   
```

---

## PSEUDOGENE:PartialDeletion

**What it means:** The gene is truncated — only 60–80% of the reference sequence is present with high identity. This is consistent with a partial deletion of one end of the gene.

**Evidence:** Hit with ≥90% identity but coverage 60–80% (below `--min-coverage`, above 60%).

**Confidence:** `identity × coverage × 0.4` (lowest — truncated gene, evidence is weakest)

**Example:**
```
geneE   genome8   PSEUDOGENE:partial:72%cov   94.8   72.1   contig_1   890000   890900   +   0.273   FAIL   
```

---

## Confidence Score

The confidence score combines alignment quality and rescue type into a single value:

```
confidence = identity_fraction × coverage_fraction × type_modifier × [synteny_modifier]
```

Where:
- `identity_fraction` = identity / 100
- `coverage_fraction` = coverage / 100
- `type_modifier` = see table above (1.0 for DROPPED, down to 0.4 for PartialDeletion)
- `synteny_modifier` = 1.2 if PASS, 0.6 if FAIL, 1.0 if not assessed

A value close to 1.0 indicates a high-confidence rescue. Values below ~0.5 should be interpreted with caution.

---

## Rescue Type Priority in Conflict Resolution

When two rescues overlap on the genome, the one with higher confidence wins. Because of the type modifiers, the natural priority order is:

```
SPLIT > DROPPED ≈ START_VARIANT > ALIGNMENT_SPLIT > FRAGMENTED > PARALOG > PSEUDOGENE:InternalStop ≈ PSEUDOGENE:Frameshift > PSEUDOGENE:PartialDeletion
```

This prevents a PARALOG cross-match from overwriting a genuine DROPPED rescue at the same locus.
