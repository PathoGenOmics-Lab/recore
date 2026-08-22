# Evolutionary events

`evolutionary_events.tsv` is an interpretation layer, not a new correction
policy. Its job is to answer a different question from the corrected matrix:

- The corrected matrix asks: **is there enough evidence to mark this cluster as
  present in this genome?**
- The evolutionary audit asks: **what biological or technical process best
  explains the cluster-level pattern after correction?**

The two answers are deliberately kept separate. A cluster can remain a real
matrix cluster while the recommended wording says "reference-gene derivative",
"paralog/duplication signal", "fragmented homology", or "review before claiming
novel acquisition".

## Evidence used

Each `molecular_mechanism` call combines:

1. corrected cluster state after accepted rescues and optional NSC support;
2. reference-genome audit and reference-derived cluster events;
3. local locus-context stability and synteny/paralog rejection evidence;
4. optional tree parsimony, including explicit ambiguous-root reporting;
5. optional sequence-based mobile-element support;
6. connected components of related clusters, so paralogous, redundant, split,
   or annotation-linked clusters are interpreted together.

Tree tips are matched to matrix columns by exact name or embedded `ASxxxxx`
identifier. If only a subset of tips matches, evolutionary counts and patterns
are calculated on the matched subset and the run summary/log warns about the
incomplete match.

The tree is used for topology-based interpretation. Recore first scores the
corrected presence/absence pattern on the matched tree using binary parsimony
and a Dollo-like single-origin model. This produces fields such as
`phylogenetic_pattern`, `phylogenetic_blink_flag`,
`carrier_mrca_purity`, `dollo_origin_purity`, and
`penalized_phylo_score`. Branch lengths are preserved for the HTML viewer, but
the mechanism call itself is conservative and is driven primarily by evidence
type, locus context, and topology.

## Decision order

The mechanism call is a deterministic evidence ladder. Direct molecular/locus
evidence has priority over a visually attractive tree pattern, because a clean
clade can still be caused by a reference allele, paralog collapse, gene
fragmentation, or annotation boundaries.

1. **Reference-derived event check**. If `reference_cluster_events.tsv` links
   the cluster to a gene in the audited reference, Recore first decides whether
   the cluster is a full-length/same-gene allele, a short/long in-frame allele,
   a frameshift/truncation/extension, or partial reference homology. If the same
   signal is also paralogous or multi-locus, it is treated as
   `duplication_or_paralog_expansion`.
2. **Mobile-element disruption check**. If a disrupted or fragmented rescue has
   nearby mobile-element evidence within the conservative window, Recore calls
   `IS_interruption_supported`. With `--mobile-element-reference`, this is
   sequence-based and name-agnostic; without it, Recore falls back to annotation
   labels such as insertion sequence, transposase, integrase, phage, or
   recombinase terms.
3. **ORF integrity check**. Pseudogene evidence is called
   `pseudogenization`; fragmented or split homology is called
   `gene_fragmentation_or_structural_break`.
4. **Paralog and locus-context check**. Explicit paralog evidence, unstable
   locus context, or multi-locus behaviour becomes
   `duplication_or_paralog_expansion`.
5. **Synteny rejection check**. Strong sequence similarity in an incompatible
   local context becomes `translocation_or_context_shift_candidate`, because it
   may represent a real sequence at the wrong locus, a translocation, or a
   paralog/annotation problem.
6. **Annotation/redundancy check**. Split genes, start-boundary shifts, or NSC
   representative/member relationships become
   `annotation_or_clustering_artifact`.
7. **Reference-absent clade check**. Only after the previous traps are ruled
   out, a cluster absent from the reference and concentrated in a lineage/clade
   can become `clean_insertion_or_acquisition_candidate`.
8. **Reference loss check**. If the reference is absent but the inferred root
   state is present, Recore prefers `reference_or_lineage_loss_candidate` over
   "gain everywhere else".
9. **Tree loss/homoplasy/retention check**. Ancestral or recurrent losses
   become `deletion_loss_candidate`; highly recurrent gain/loss patterns become
   `unstable_or_homoplastic_accessory_signal`; conserved ancestral presence
   becomes `ancestral_retained_no_disruption`.
10. **Fallback**. If the available matrix, reference, synteny, mobile-element,
    and tree evidence do not resolve the process, Recore reports
    `unresolved_accessory_variation`.

## Mechanism meanings

| `molecular_mechanism` | When Recore assigns it | How to interpret it |
|-----------------------|------------------------|---------------------|
| `ancestral_retained_no_disruption` | The corrected cluster is broadly present and the tree pattern supports ancestral retention. | Conserved ancestral presence in the analysed set; not a novel acquisition. |
| `clean_insertion_or_acquisition_candidate` | The cluster is absent from the audited reference, present in a coherent lineage/clade, and lacks stronger reference/paralog/fragmentation/artifact explanations. | Best candidate for reference-absent accessory acquisition or insertion; validate flanks before claiming a novel locus. |
| `reference_or_lineage_loss_candidate` | The reference is absent, but the inferred root/background state is present. | More consistent with loss in the reference or reference lineage than independent gain in all other genomes. |
| `deletion_loss_candidate` | Tree pattern suggests ancestral presence with one or more losses. | Candidate deletion/loss event; inspect flanking genes, breakpoints, and mobile elements. |
| `IS_interruption_supported` | Fragmented/disrupted evidence is close to mobile-element evidence. | Candidate IS/transposase/mobile-element-mediated interruption; not automatically intact gene presence. |
| `pseudogenization` | Rescue evidence indicates frameshift, internal stop, or pseudogene-like disruption. | Homologous sequence exists, but should be described as disrupted unless additional evidence proves function. |
| `gene_fragmentation_or_structural_break` | Evidence is split, fragmented, or alignment-split. | Partial homology or local structural break; review before using as gain/loss evidence. |
| `duplication_or_paralog_expansion` | Paralog evidence, multi-locus context, or copy-number-like behaviour is detected. | Copy-number/paralog/cluster-definition complexity; avoid treating as one clean locus without manual locus review. |
| `reference_gene_allelic_derivative` | The cluster is best explained as same-gene, short in-frame, or long in-frame reference-derived sequence. | Event-derived allele/version of a reference gene; not a completely new gene. |
| `reference_gene_disruption_or_truncation` | Reference-like sequence has frameshift, truncation, or extension evidence. | Reference-gene-derived disrupted/truncated form; interpret with ORF integrity. |
| `partial_reference_homology` | The representative sequence partially matches a reference gene. | Some homology to reference exists; do not call "absent from reference" as if it were wholly novel. |
| `translocation_or_context_shift_candidate` | Sequence evidence exists but synteny/local context is incompatible. | Real sequence may be in a different context, translocated, paralogous, or annotation-driven; review locus manually. |
| `annotation_or_clustering_artifact` | Split/start-boundary evidence, NSC redundancy, or cluster-component evidence explains the signal. | Likely cluster definition, boundary, or redundancy issue rather than a biological gene gain/loss. |
| `unstable_or_homoplastic_accessory_signal` | The tree shows repeated appearances/disappearances or high blinking. | Could be real mobility/recombination, but also annotation/assembly/clustering noise; deprioritize for headline claims. |
| `unresolved_accessory_variation` | Evidence is insufficient or conflicting after all checks. | Keep visible for audit; do not use as a strong mechanistic claim without external validation. |

## How to use the calls

Use `phylogenetic_pattern` to describe **where** the signal sits on the tree.
Use `molecular_mechanism` to describe **what process is the safest
interpretation**. Use `recommended_wording` for manuscript text, because it is
written to avoid overclaiming.

For manuscript highlighting, the safest candidates are usually clusters with:

- `clean_insertion_or_acquisition_candidate`,
  `reference_or_lineage_loss_candidate`, or `deletion_loss_candidate`;
- coherent `phylogenetic_pattern` with low/moderate blink flag;
- high `carrier_mrca_purity` or `dollo_origin_purity`;
- no paralog/component warning;
- compatible synteny and clear ORF integrity.

Clusters with `duplication_or_paralog_expansion`,
`annotation_or_clustering_artifact`,
`translocation_or_context_shift_candidate`, or
`unstable_or_homoplastic_accessory_signal` can still be biologically
interesting, but they should be presented as review/complex-locus signals until
the locus is inspected manually. Clusters with `reference_gene_*` or
`partial_reference_homology` should not be described as completely novel
reference-absent genes, even if the original pangenome matrix marked them as
absent from the reference.

## Where the fields live

Every field named on this page is a column of `evolutionary_events.tsv`; its
schema is in [Output files](../guide/outputs.md#evolutionary_eventstsv). The
report colours the matrix by `molecular_mechanism` when you ask it to, and its
mechanism cards double as column filters: see [Run report](../guide/report.md).

