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
| `synteny_score` | float | Best adaptive context score across 1-to-5-neighbour windows; blank when not assessed |
| `synteny_support_n` | int | Number of compatible reference contexts supporting the best synteny assessment |
| `synteny_reference_contexts` | int | Number of reference locus contexts compared for the selected neighbour window |
| `synteny_neighbors` | int | Selected upstream/downstream neighbour window size |
| `synteny_upstream_matches`, `synteny_downstream_matches` | int | Best upstream/downstream cluster matches observed in the selected assessment |
| `orf_integrity_signal` | string | `in_frame_span`, `out_of_frame_span_review`, `fragmented_orf`, `disrupted_orf`, or `unknown` |
| `annotation_match` | string | GFF gene ID that overlapped the rescue region, if any |
| `fragment_n` | int | Number of fragments in the selected multi-hit chain, if applicable |
| `fragment_query_gap_bp` | int | Total query-space gap across the selected fragment chain |
| `fragment_target_gap_bp` | int | Total target-space gap across the selected fragment chain |
| `fragment_colinear` | bool | Whether fragments are co-linear on the target |
| `fragment_same_cds` | bool | Whether the chain is contained in one compatible CDS |
| `fragment_note` | string | Fragment interpretation, such as `alignment_split`, `biological_gap`, `frameshift_gap`, or `same_cds_alignment_split` |

### Example

```tsv
gene_cluster	genome	rescue_type	identity	coverage	contig	start	end	strand	confidence	synteny	synteny_score	synteny_support_n	synteny_reference_contexts	synteny_neighbors	synteny_upstream_matches	synteny_downstream_matches	orf_integrity_signal	annotation_match	fragment_n	fragment_query_gap_bp	fragment_target_gap_bp	fragment_colinear	fragment_same_cds	fragment_note
group_2044	L1_AS00349	DROPPED	100.0	100.0	L1_AS00349	245301	246100	+	1.198	PASS	1.000	12	18	1	1	1	in_frame_span	Rv0210c	NA	NA	NA	NA	NA	NA
group_681	L4_AS00146	MISSING_CDS	99.2	89.5	L4_AS00146	2195924	2197220	+	1.066	PASS	0.833	4	5	2	1	1	in_frame_span		NA	NA	NA	NA	NA	missing_cds_no_substantial_gff_overlap
group_0187	A4_AS00248	SPLIT:group_0188	99.2	99.8	A4_AS00248	1203401	1204300	+	1.087	PASS	1.000	8	9	1	1	1	in_frame_span	Rv1099c	NA	NA	NA	NA	NA	NA
group_3940	L1_AS00366	PSEUDOGENE:2stops	99.0	100.0	L1_AS00366	1783748	1784254	-	0.495	PASS	0.750	3	4	2	1	1	disrupted_orf	Rv2651c	NA	NA	NA	NA	NA	NA
group_2101	L2_AS00401	ALIGNMENT_SPLIT:2frags:98%cov	99.6	98.4	L2_AS00401	552100	553020	+	0.882	PASS	1.000	6	6	1	1	1	fragmented_orf		2	9	9	true	false	alignment_split
```

### Filtering suggestions

```bash
# High-confidence rescues only
awk '$10 > 0.8' rescue_report.tsv

# Synteny-validated only
awk '$11 == "PASS"' rescue_report.tsv

# Non-pseudogene rescues
grep -v "PSEUDOGENE" rescue_report.tsv

# Same-locus full-length rescues with strong support
grep -E "DROPPED|MISSING_CDS|SPLIT" rescue_report.tsv | awk '$10 > 0.9 && $11 == "PASS"'
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

## paralog_rejected_report.tsv

Written only when one or more PARALOG candidates lack explicit synteny support.
The schema is identical to `rescue_report.tsv`, but these rows are audit
records: they are not added to `gene_presence_absence_recore.csv` and are not
written to corrected GFF output.

This is the conservative guardrail for multi-copy families: sequence similarity
alone is not allowed to rescue a missing gene when the annotation says the hit
belongs to a different cluster.

---

## orthology_evidence.tsv

Written for every run. One row per accepted rescue plus every rejected synteny
or paralog candidate.

| Column | Description |
|--------|-------------|
| `decision` | `accepted` or `rejected` |
| `orthology_status` | Conservative interpretation, e.g. `accepted_synteny`, `accepted_boundary_variant`, `rejected_incompatible_context`, or `rejected_paralog_without_synteny` |
| `evidence_state` | Explicit evidence class: `STRICT_PRESENT`, `ALIAS_CLUSTER`, `SEQUENCE_ONLY`, `MISSING_CDS_UNVALIDATED`, `PARTIAL_HOMOLOGY`, `PSEUDOGENE`, `HOMOLOG_PARALOG`, or `REJECTED_CONTEXT` |
| `matrix_action` | `correct_presence` for accepted rescues, `audit_only_quality` for accepted evidence that fails the adaptive guard, `audit_only_policy` for accepted evidence excluded by policy, or `audit_only` for rejected candidates |
| `publication_status` | Conservative publication-level interpretation of the evidence |
| `evidence_strength` | `high`, `moderate`, `review`, or `low` |
| `recommended_wording` | Short text for how to describe the event in downstream tables/manuscripts |
| `gene_cluster`, `genome`, `rescue_type`, `identity`, `coverage`, `contig`, `start`, `end`, `strand`, `confidence`, `synteny`, `annotation_match` | Same meaning as `rescue_report.tsv` |
| `adaptive_identity_floor`, `adaptive_coverage_floor`, `adaptive_sequence_pass` | Per-cluster sequence guard used before a rescue can enter the binary matrix; rare or short clusters require stronger evidence unless same-locus context supports the rescue |
| `orthology_confidence_score`, `orthology_confidence_tier`, `risk_flags` | Integrated confidence and review flags combining sequence, synteny, ORF, state, prevalence, and cluster length |
| `synteny_score`, `synteny_support_n`, `synteny_reference_contexts`, `synteny_neighbors`, `synteny_upstream_matches`, `synteny_downstream_matches` | Quantitative adaptive synteny evidence, same meaning as `rescue_report.tsv` |
| `orf_integrity_signal` | ORF-level signal, same meaning as `rescue_report.tsv` |
| `annotation_cluster` | Cluster assigned to the overlapping annotation, when known |
| `rejection_reason` | Reason rejected, e.g. `synteny_fail` or `paralog_without_synteny` |
| `fragment_note` | Fragment interpretation, if applicable |

Use this table as the primary audit trail for why a candidate entered, or did
not enter, the corrected matrix.

---

## reference_genome_audit.tsv

Written only when `--reference-genome` is provided. The value can be an exact
matrix column name or an embedded `ASxxxxx` identifier.

| Column | Description |
|--------|-------------|
| `gene_cluster` | Panaroo cluster name |
| `reference_genome` | Resolved matrix column used as the audited reference |
| `reference_absence_class` | `new_to_reference`, `reference_like_or_present`, or `reference_like_ambiguous` |
| `final_publication_status` | Whether to keep in a strict reference-absent set, remove, or review |
| `reference_status` | Detailed reason: matrix-present, clean absent, strict rescue, alias, partial-like, pseudogene-like, paralog-like, or incompatible-context evidence |
| `matrix_action`, `evidence_state`, `evidence_strength`, `publication_confidence`, `recommended_wording` | Same interpretation layer used in `orthology_evidence.tsv`, plus a publication-facing confidence label |
| `rescue_decision`, `rescue_type`, `rejection_reason`, `identity`, `coverage`, `contig`, `start`, `end`, `strand`, `synteny`, `annotation_match`, `annotation_cluster`, `fragment_note` | Evidence fields when reference-like sequence evidence exists |

This file is meant for reference-absence claims: keep `new_to_reference` for
strictly absent clusters, remove `reference_like_or_present`, and manually
review `reference_like_ambiguous`. The `publication_confidence` field uses
`clean`, `complex_paralog`, `event_derived`, `fragmented`, `pseudogene`,
`unresolved`, or `not_reference_absent`. `complex_paralog` is used when the
exact cluster remains absent from the reference but discarded/subthreshold
sequence evidence points to a partial or paralogous reference-like copy.
`event_derived` is used when the exact cluster remains separate in the matrix
but its representative sequence is best explained as a reference-gene allele or
fragment. Accepted sequence-only evidence that fails the adaptive
orthology-quality guard is kept as `reference_like_ambiguous` with
`publication_confidence=unresolved`.

---

## reference_cluster_events.tsv

Written only when `--reference-genome` is provided and Recore finds
reference-absent clusters with strong sequence evidence to a reference-present
cluster. This file is audit-only: it never merges clusters and never changes
the binary corrected matrix.

| Column | Description |
|--------|-------------|
| `gene_cluster` | Cluster absent from the audited reference |
| `target_reference_cluster` | Best reference-present cluster matched by k-mer seed plus global alignment |
| `target_reference_gene` | Gene ID carried by the audited reference for the target cluster |
| `reference_event_class` | `same_gene_full_length`, `short_in_frame_allele`, `long_in_frame_allele`, `frameshift_or_truncated_reference_homology`, `frameshift_or_extension_reference_homology`, or `partial_reference_homology` |
| `identity`, `query_coverage`, `reference_coverage` | Global-alignment summary between the cluster representative and the matched reference-cluster representative |
| `query_length`, `reference_length`, `event_size_nt`, `event_size_aa` | Length comparison. Amino-acid delta is filled only when the nucleotide length difference is divisible by three |
| `matches`, `mismatches`, `gaps_in_query`, `gaps_in_reference`, `shared_kmers` | Alignment and seed evidence used to choose the best event |
| `matrix_recommendation` | Currently `keep_separate_cluster_annotate_event` |
| `interpretation` | Publication-facing wording for the event |

Rows in this file also appear in `reference_genome_audit.tsv` as
`EVENT_DERIVED_REFERENCE_CLUSTER` with
`final_publication_status=flag_reference_event_derived_cluster` and
`publication_confidence=event_derived`.

---

## evolutionary_events.tsv

Written for every run. This file is audit-only: it does not merge clusters and
does not change `gene_presence_absence_recore.csv`. When `--tree` is provided,
Recore parses the rooted Newick tree, matches tips to matrix columns by exact
name or embedded `ASxxxxx` identifier, and applies Fitch parsimony to the
corrected binary state for each cluster.

| Column | Description |
|--------|-------------|
| `gene_cluster` | Panaroo cluster name |
| `n_present_corrected`, `prevalence_corrected` | Presence after accepted matrix corrections and NSC presence updates |
| `tree_matched_tips` | Number of tree tips matched to matrix genomes; blank when no tree is provided |
| `reference_state` | `present`, `absent`, or `not_assessed` for the requested `--reference-genome` |
| `phylogenetic_pattern` | `ancestral_retained`, `ancestral_loss`, `lineage_specific_gain`, `clade_associated_presence`, `clade_associated_presence_with_losses`, `recurrent_gain`, `recurrent_loss`, `homoplastic_gain_loss`, `ambiguous_gain_loss`, `singleton_or_noise`, `dispersed_gain_loss`, or `not_assessed_no_tree` |
| `phylogenetic_reliability` | `high`, `moderate`, `low`, or `not_assessed`, based on tree pattern, root-state ambiguity, MRCA purity, and parsimony changes |
| `phylogenetic_blink_flag` | `low_blink`, `moderate_blink`, `high_blink`, or `not_assessed`; high blink marks clusters with singleton/dispersed/homoplastic behaviour that should be filtered or reviewed before highlighting |
| `penalized_phylo_model`, `penalized_phylo_score`, `dollo_loss_count`, `dollo_violation_count`, `dollo_origin_tips`, `dollo_origin_purity` | Dollo-like penalized model using the smallest carrier MRCA as the single-gain origin; broad origins, many losses, root ambiguity, and repeated independent gains reduce the score |
| `ancestral_state`, `root_state_confidence` | Parsimony root state: `present`, `absent`, `ambiguous`, or `unknown`, plus whether the Fitch root assignment is unambiguous |
| `phylogenetic_changes`, `phylogenetic_gains`, `phylogenetic_losses` | Parsimony counts on the matched tree |
| `carrier_monophyly`, `max_carrier_clade_fraction` | Whether all carrier tips form one pure carrier-only clade, and the largest carrier-only clade fraction |
| `carrier_mrca_tips`, `carrier_mrca_purity` | Size and carrier purity of the smallest tree clade containing all carrier tips; low purity means a single-gain story would require many internal losses |
| `molecular_mechanism` | Best conservative mechanism call: reference-derived allele, IS interruption, pseudogenization, fragmentation, duplication/paralog expansion, context shift, acquisition/insertion candidate, deletion/loss candidate, homoplastic/unstable signal, ancestral retention, or unresolved |
| `mechanism_detail`, `support_level`, `evolutionary_signal_score`, `highlight_priority`, `recommended_wording` | Evidence summary, integrated prioritization, and publication-facing interpretation |
| `cluster_component_id`, `cluster_component_size`, `cluster_component_role`, `cluster_component_relations`, `cluster_component_recommendation` | Connected-component context for split, annotation-matched, NSC, and paralog-linked clusters |
| `is_interruption_supported`, `mobile_element_*` | Nearby mobile element signal. `mobile_element_source=sequence` comes from `--mobile-element-reference` alignments and is name-agnostic; `mobile_element_source=annotation` is the fallback from GFF labels/products/notes. IS calls are conservative and require sequence support or an insertion-sequence-like label within 250 bp of a fragmented/disrupted rescue |
| `reference_event_class`, `target_reference_cluster`, `target_reference_gene` | Reference-derived event evidence when available |
| `context_interpretation`, `paralog_locus_risk` | Locus-context stability and paralog/multilocus risk summary |
| `synteny_mean_score`, `synteny_pass_fraction` | Cluster-level summary of quantitative synteny support across assessed rescue evidence |
| `orf_integrity_summary`, `rescue_evidence_counts` | ORF signals and per-cluster accepted rescue-type summary |

Use this file as the top-level evolutionary interpretation table. The most
important guardrail is that `phylogenetic_pattern` and `molecular_mechanism`
are deliberately separate: a clade-associated presence can be an acquisition,
a reference loss, a duplicated/paralogous copy, or an event-derived allele
depending on the sequence, synteny, GFF, and reference evidence.

---

## cluster_locus_context.tsv

Written when GFF annotations are loaded. One row per observed locus context for
each cluster.

| Column | Description |
|--------|-------------|
| `gene_cluster` | Panaroo cluster name |
| `context_rank` | Rank within the cluster by number of genomes |
| `n_genomes` | Number of genomes supporting this context |
| `fraction_of_mapped_contexts` | Fraction among genomes where the cluster gene ID was resolved in the GFF |
| `upstream_context` | Nearest upstream cluster names, joined by `|`; `NA` if no mapped context |
| `downstream_context` | Nearest downstream cluster names, joined by `|`; `NA` if no mapped context |
| `example_genomes` | Up to five genome:gene examples |

---

## cluster_locus_summary.tsv

Written when GFF annotations are loaded. One row per cluster.

| Column | Description |
|--------|-------------|
| `gene_cluster` | Panaroo cluster name |
| `n_present` | Number of genomes where the cluster is present in the input matrix |
| `n_context_mapped` | Present genomes where the cell gene ID was found in the GFF |
| `n_unresolved_gene_ids` | Present genomes where the cell gene ID could not be resolved in the GFF |
| `n_unique_contexts` | Number of distinct upstream/downstream context signatures |
| `dominant_context_count` | Number of genomes in the largest context |
| `dominant_context_fraction` | Dominant context count divided by mapped contexts |
| `interpretation` | `stable_single_locus`, `stable_single_locus_minor_noise`, `mostly_stable_review_minor_contexts`, `possible_multilocus_or_paralog`, or `insufficient_context` |

Clusters flagged as `possible_multilocus_or_paralog` are not automatically
changed, but they should be reviewed before treating accessory-gene calls as
biological.

---

## cluster_synteny_blocks.tsv

Written when GFF annotations are loaded. This is an audit-only synteny-block
summary using up to five upstream and five downstream cluster neighbours.

| Column | Description |
|--------|-------------|
| `gene_cluster` | Panaroo cluster name |
| `block_rank` | Rank by number of supporting genomes |
| `n_genomes`, `fraction_of_present` | Support for this block among genomes where the cluster is present |
| `upstream_block_5`, `downstream_block_5` | Cluster-name signatures for the wider flanking block |
| `block_interpretation` | `single_synteny_block`, `dominant_block_minor_alternatives`, or `alternative_synteny_block` |
| `example_genomes` | Up to five genome:gene examples |

---

## gene_presence_absence_recore.csv

Drop-in replacement for Panaroo's `gene_presence_absence.csv`. Identical format and column order.

### Changes from the original

1. **Rescued cells**: Empty cells (gene absent) become `RESCUED:contig:start-end` where coordinates are 1-based.
   ```
   Before: ""
   After:  "RESCUED:L1_AS00349:245301-246100"
   ```

2. **Matrix-policy filter**: `--matrix-policy` controls which accepted rescues become binary presence. `strict` corrects only `STRICT_PRESENT`/`ALIAS_CLUSTER`; `balanced` also corrects full-length `SEQUENCE_ONLY`; `sensitive` also permits partial/disrupted accepted evidence. In all modes, the adaptive orthology-quality guard can keep weak accepted evidence out of the binary matrix.

3. **NSC-merged cells** (only if `--nscluster`): Cells where a merged cluster would have added presence become `NSC_MERGED`.

4. **Metadata columns updated**: `No. isolates`, `No. sequences`, and `Avg sequences per isolate` are recomputed from the final cells. Multicopy cells separated by `;` count as multiple sequences but one isolate.

5. **NSC clusters removed**: Panaroo clusters that were merged by NSC into a representative cluster are dropped entirely (their rows are absent from the output).

### Compatibility

This file can be used directly as input to downstream tools that accept Panaroo's format:
- `distree` for phylogenetic distance analysis
- `fstic` for allele frequency-based FST calculations
- Custom R/Python scripts expecting Panaroo output.

---

## gene_presence_absence_recore_states.csv

Written for every run. Same row/column layout as the corrected matrix, but each
genome cell stores the biological interpretation rather than a binary presence
call.

Possible states:
- `PRESENT_INPUT`: present in the original matrix.
- `ABSENT`: no evidence under current thresholds.
- `STRICT_PRESENT`: accepted same-locus rescue.
- `ALIAS_CLUSTER`: accepted split/alias or NSC-compatible cluster evidence.
- `SEQUENCE_ONLY`: full-length sequence evidence without context/annotation support.
- `PARTIAL_HOMOLOGY`: accepted fragmented or alignment-split evidence.
- `PSEUDOGENE`: accepted disrupted sequence evidence.
- `HOMOLOG_PARALOG`: rejected or ambiguous paralog-like evidence.
- `REJECTED_CONTEXT`: rejected incompatible-context evidence.

Use this file when a binary matrix would overstate biology, especially for
pseudogenes, partial genes, PE/PPE-like families, or publication figures that
need to show uncertainty explicitly.

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

## cluster_components.tsv / cluster_component_members.tsv

Written for every run. These tables convert pairwise relationships into
connected components, so split clusters, annotation-matched clusters, NSC
members, and paralog-linked clusters can be reviewed as one unit.

`cluster_components.tsv` contains one row per non-singleton component:

| Column | Description |
|--------|-------------|
| `component_id` | Stable component label for the run |
| `n_clusters` | Number of clusters in the component |
| `relation_types` | Relationship classes connecting the component |
| `recommendation` | Component-level action, such as `review_or_collapse_aliases` or `review_component_keep_paralogs_distinct` |
| `member_clusters` | Semicolon-separated cluster names |

`cluster_component_members.tsv` contains one row per member cluster:

| Column | Description |
|--------|-------------|
| `gene_cluster` | Panaroo cluster name |
| `component_id`, `component_size` | Component membership |
| `component_role` | `nsc_representative`, `nsc_member`, `related_cluster_member`, or `paralog_component_member` |
| `relation_types`, `recommendation` | Same component evidence and review action |

These components are also reported in `evolutionary_events.tsv` so manuscript
highlighting can deprioritize redundant aliases or paralog components.

---

## genome_qc.tsv

Written for every run. One row per genome with rescue and annotation burden
metrics.

| Column | Description |
|--------|-------------|
| `genome` | Matrix genome column |
| `fasta_found`, `gff_found` | Whether matching FASTA/GFF inputs were available |
| `gff_contig_remapped_records` | Number of GFF records remapped to FASTA contig names |
| `gaps_passing_filter` | Missing cells searched or eligible for search after prevalence filtering |
| `accepted_rescues`, `matrix_corrected_rescues`, `nsc_added_presence` | Evidence and matrix-correction burden by genome |
| `synteny_rejected`, `paralog_rejected` | Rejected evidence by genome |
| `qc_status`, `qc_flags` | `pass`, `review`, or `fail` plus semicolon-separated flags |

`gff_contig_seqid_remapped` is informational when remapping succeeds; it is
kept in `qc_flags` but does not by itself trigger `review`.

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
  MISSING_CDS:   0
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

## run_manifest.json

Machine-readable run manifest written for every run. It records:

- `tool` and `tool_version`;
- `schema_version` (`recore.outputs.v6`);
- input paths used for the run;
- key parameters, including matrix policy, prevalence threshold, synteny, NSC,
  reference genome, tree, and mobile-element settings;
- headline counts such as clusters, genomes, searched gaps, direct rescue
  corrections, NSC-merged corrections, rejected synteny/paralog candidates, and
  NSC groups.

Use this file when archiving outputs for a manuscript or comparing runs across
parameter settings.

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
records are preserved and Recore appends rescued CDS features. If the original
GFF contains a `##FASTA` section, rescue features are inserted before `##FASTA`
and the FASTA tail is preserved after them. If no original GFF is available for
a rescued genome, the file contains only the appended rescue features.

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

## recore_report.html

Written only when `--write-html-report` is enabled.

This is a self-contained browser report for run review and manuscript-facing
QC. It summarizes:

- input matrix size, searchable gaps, accepted evidence, corrected matrix cells,
  rejected synteny/paralog evidence, and core/soft-core changes;
- orthology actions and confidence tiers from `orthology_evidence.tsv`;
- evolutionary priorities and dominant mechanisms from `evolutionary_events.tsv`;
- per-genome QC status from `genome_qc.tsv`;
- the top candidate clusters ranked by `evolutionary_signal_score`;
- direct links to the exact output files used to generate the summary.

When `--tree` is supplied and the matched tree/matrix is small enough for a
responsive browser view, the report also embeds a phylogeny plus corrected
matrix explorer. The explorer prunes the tree to matched matrix tips, keeps
matrix rows in tree order, lets the user search/select clusters, recolors tree
tips dynamically by the selected cluster state or evolutionary mechanism, and
lets users click matrix columns to change the selected cluster. The matrix can
also be recolored by evolutionary mechanism, with a legend explaining each
observed process from `evolutionary_events.tsv`. MTBC-style genome names can be
colored by lineage with the mycolorsTB Nature palette. Clicking a mechanism
card filters the matrix to clusters with that mechanism, and selected-cluster
details show the reference cell/gene or the reference-derived event source when
available. Dense views can be collapsed or zoomed with the Y-axis slider, column
zoom slider, mouse wheel over the tree or column header, and fullscreen mode.

The embedded explorer is intentionally size-limited. Larger runs keep the
static audit summary and report why the explorer was omitted, rather than
writing a sluggish HTML file.

The HTML report is an audit dashboard, not a replacement for the TSV outputs.
Use the TSV files for exact filtering, statistics, and downstream figures.

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
