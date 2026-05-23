# Changelog

All notable changes to Recore are documented in this file.

The format follows the spirit of [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and Recore uses semantic versioning.

## [Unreleased]

No changes yet.

## [0.2.0] - 2026-05-23

### Added

- Publication-grade matrix policies: `strict`, `balanced`, and `sensitive`.
- Multistate presence/absence output in `gene_presence_absence_recore_states.csv`.
- `orthology_evidence.tsv` with explicit accepted/rejected evidence, matrix action,
  publication status, evidence strength, recommended wording, coordinates, synteny,
  and annotation-cluster interpretation.
- `reference_genome_audit.tsv` for auditing whether clusters are truly absent from
  a chosen reference genome.
- Reference-only shadow search for `--reference-genome`, so partial or paralogous
  reference-like evidence is preserved even when it is excluded from the binary
  matrix.
- `publication_confidence` labels in the reference audit, including `clean`,
  `complex_paralog`, `fragmented`, `pseudogene`, `unresolved`, and
  `not_reference_absent`.
- Exact missing-CDS rescue pass in GFF/synteny mode to recover unannotated CDS
  loci that can be hidden by nearby annotated paralogs.
- Optional corrected GFF output in `corrected_gff/*.recore.gff3`.
- Cluster relationship auditing for split, alias-like, annotation-matched,
  nucleotide-similar, and paralog-like cluster relationships.
- Nucleotide similarity clustering (`--nscluster`) with built-in MinHash and
  optional mmseqs2 engines.
- Synthetic precision/recall benchmark and AMAP-scale performance benchmark.
- README logo, workflow diagram, benchmark plot, and expanded documentation for
  algorithms, parameters, rescue types, and output files.
- Event-aware reference cluster audit for `--reference-genome`.
- `reference_cluster_events.tsv`, an audit-only table that links
  reference-absent clusters to the best reference-present cluster when the
  sequence supports a same-gene, in-frame allele, frameshift-like, or partial
  reference-homology event.
- `EVENT_DERIVED_REFERENCE_CLUSTER` and
  `publication_confidence=event_derived` in `reference_genome_audit.tsv`, so
  event-derived reference-gene variants are not reported as clean novel genes.
- `--tree` and `evolutionary_events.tsv` for audit-only evolutionary
  interpretation of each cluster, including ancestral state, parsimony
  gain/loss counts, carrier monophyly, reference-derived alleles, IS-supported
  interruptions, pseudogenization, fragmentation, duplication/paralogy,
  context-shift candidates, acquisition/insertion candidates, and deletion/loss
  candidates.
- `--mobile-element-reference` plus identity/coverage/window controls for
  name-agnostic mobile-element support. Recore can now support IS/transposase
  interruption calls by sequence proximity instead of relying only on GFF
  product names.
- `run_manifest.json` with tool version, output schema version, input paths,
  key parameters, and headline run counts for reproducible archiving.
- Quantitative adaptive synteny evidence in `rescue_report.tsv` and
  `orthology_evidence.tsv`, including score, support count, reference-context
  count, selected neighbour window, and upstream/downstream matches.
- ORF integrity signals in rescue and orthology reports.
- Cluster-level `phylogenetic_reliability`, `paralog_locus_risk`,
  `synteny_mean_score`, `synteny_pass_fraction`, and `orf_integrity_summary`
  fields in `evolutionary_events.tsv`.
- Biological event benchmark for same-locus missing-CDS rescue, wrong-context
  paralog rejection, and evolutionary audit columns.
- Output manifest schema version `recore.outputs.v2` for the expanded TSV
  report schema.

### Changed

- Batched minimap2 search now uses streamed PAF parsing and indexed concatenated
  targets to reduce IO overhead.
- Rescue decisions are more conservative for fragments, pseudogenes, paralogs,
  and incompatible synteny contexts.
- Synteny validation now uses adaptive 1-to-5-neighbour quantitative context
  scoring rather than a single exact-context match.
- GFF seqids can be normalized against FASTA contig names to handle common
  annotation/assembly naming mismatches.
- Gene-to-cluster matching now uses unambiguous GFF aliases, so locus tags and
  gene names from different annotation layers are treated consistently.
- Multi-copy reference cells now contribute all resolved gene IDs to synteny
  validation instead of only the first semicolon-separated ID.
- Output reports are sorted deterministically for reproducible downstream review.
- TSV report fields are normalized to single-line values so annotation-derived
  tabs or line breaks cannot corrupt tabular outputs.

### Fixed

- Explicit synteny failures are no longer inserted into the corrected binary
  matrix; they are retained as audit-only evidence.
- Paralog-only sequence hits require explicit synteny support before being
  treated as orthology evidence.
- Missing CDS rescues can now be detected when an annotated paralog would
  otherwise absorb the strongest minimap2 hit.
- Reference-absence claims no longer classify partial/paralogous reference-like
  evidence as `NO_EVIDENCE`; these cases are flagged as
  `PARTIAL_PARALOG_EVIDENCE` with `publication_confidence=complex_paralog`.
- NSC representative clusters now inherit matrix-corrected rescue evidence
  from redundant member clusters, not only original matrix presence.
- `summary.txt` now counts NSC-merged matrix cells separately from direct
  rescue cells and reports the combined corrected-cell total.
- `sensitive` matrix policy no longer treats homologous/paralogous evidence
  states as binary-correctable unless they have already passed same-locus
  support checks.
- Tree-based evolutionary audits now warn when only a subset of matrix genomes
  is matched by the supplied Newick tree.
- Invalid numeric options and mistyped NSC engines are rejected before the
  pipeline starts, avoiding silent fallback to unintended behaviour.

## [0.1.0] - 2026-04-04

### Added

- Initial public release of Recore.
- Panaroo/Roary-style `gene_presence_absence.csv` parsing.
- Pan-genome reference FASTA parsing.
- Genome FASTA discovery and minimap2-based search for matrix gaps.
- Basic rescue classification and `rescue_report.tsv`.
- Panaroo-compatible corrected binary matrix output.
- Summary report with total gaps, rescued gaps, remaining gaps, and core-genome
  change.
- Command-line interface with configurable identity, coverage, prevalence,
  preset, thread count, and output directory.
- CI and release workflow scaffolding.

[Unreleased]: https://github.com/PathoGenOmics-Lab/recore/compare/0.2.0...HEAD
[0.2.0]: https://github.com/PathoGenOmics-Lab/recore/compare/v0.1.0...0.2.0
[0.1.0]: https://github.com/PathoGenOmics-Lab/recore/releases/tag/v0.1.0
