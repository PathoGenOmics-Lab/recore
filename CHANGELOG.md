# Changelog

All notable changes to Recore are documented in this file.

The format follows the spirit of [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and Recore uses semantic versioning.

## [Unreleased]

### Changed

- Presence/absence storage now uses compact indexed vectors instead of a
  per-cluster hash set, reducing memory overhead and making matrix lookups
  cache-friendly.
- Gap filtering and PAF hit filtering now use dense indexed lookup tables
  instead of hash maps/sets for genome/cluster pairs.
- Release builds use a single codegen unit for stronger whole-binary
  optimization.
- Matrix correction now applies an adaptive orthology-quality guard on top of
  the selected matrix policy. Rare and short clusters require stronger
  sequence evidence unless same-locus/synteny evidence supports the rescue.
- `orthology_evidence.tsv` now reports adaptive identity/coverage floors,
  whether the sequence evidence passes them, an integrated orthology confidence
  score, confidence tier, and risk flags.
- Related clusters are now grouped into explicit components in
  `cluster_components.tsv` and `cluster_component_members.tsv`, so split,
  annotation-matched, NSC, and paralog-linked clusters can be reviewed as one
  biological/technical unit.
- `evolutionary_events.tsv` now includes phylogenetic blink flags,
  evolutionary signal scores, highlight priorities, and cluster-component
  context for each cluster.
- `evolutionary_events.tsv` now includes Dollo-like penalized phylogenetic
  scores, violation counts, root-state confidence, carrier-MRCA size/purity,
  and Dollo origin size/purity.
- `cluster_synteny_blocks.tsv` reports wider five-neighbour synteny blocks for
  audit-only structural/context review.
- `genome_qc.tsv` reports per-genome FASTA/GFF availability, remapping,
  rescue burden, matrix-correction burden, and rejection burden.
- Added reproducible utilities: `benchmarks/self_benchmark.py`,
  `benchmarks/compare_runs.py`, and `benchmarks/recore_report.py`.
- Added `--write-html-report` to write a self-contained HTML run report with
  headline QC, correction counts, mechanisms, top candidate clusters, and
  links to the underlying audit tables.
- The HTML report now embeds an interactive phylogeny plus matrix explorer
  for large `--tree` runs, with searchable clusters, clickable matrix columns,
  and dynamic tip coloring by corrected gene state.
- The HTML report layout was redesigned as a more mature audit document, with
  a sober header, executive summary, clearer review metrics, and richer
  selected-cluster details in the explorer.
- The embedded tree/matrix viewer now uses an aligned vertical viewport,
  synchronized matrix column header, adjustable cell size, row/column hover
  highlighting, and per-cell tooltips.
- The embedded tree/matrix viewer now caches the matrix canvas, auto-compacts
  dense trees, supports mouse inspection on the tree and column header, and can
  recolor tips by selected cluster state, rescue burden, absence burden, or any
  correction burden.
- The embedded viewer now preserves Newick branch lengths, roots the displayed
  tree on the L8 outgroup when present, switches between phylogram and
  cladogram modes, supports fullscreen exploration, and can restrict columns to
  clusters absent from the selected reference genome.
- The embedded viewer can now zoom/collapse matrix rows and columns with
  sliders or mouse wheel modifiers, and can color both tree tips and matrix
  cells by evolutionary mechanism with an in-report explanation for each
  observed process.
- Output manifest schema version is now `recore.outputs.v6`.
- README now documents the current v6 audit guarantees, input guardrails,
  multicopy-corrected matrix metadata, and evolutionary interpretation fields.

### Fixed

- Dataset self-benchmark reports `precision_proxy=NA` when no negative controls
  are sampled instead of implying perfect precision.
- Run comparison now includes cells present in only one run, so cluster/genome
  set changes are not silently ignored.
- Successful GFF contig remapping is treated as informational genome QC rather
  than forcing every remapped genome into `review`.
- Dollo-like phylogenetic scoring no longer rewards clusters absent from all
  matched tree tips.
- Evolutionary event calls now preserve ambiguous Fitch root states instead of
  forcing tie-breaks into overconfident ancestral-loss calls.
- Corrected matrices now recompute `No. isolates`, `No. sequences`, and
  `Avg sequences per isolate` from the final multicopy-aware cell values.
- Optional output files and the `corrected_gff/` directory are cleared before
  writing a new run, preventing stale reports when an output directory is
  reused.
- Corrected GFF3 export now inserts Recore rescue features before any original
  `##FASTA` section, preserving valid GFF3 ordering.
- `reference_genome_audit.tsv` now applies the adaptive orthology-quality guard
  before classifying accepted weak sequence-only evidence as reference-like or
  present.
- Synthetic and biological benchmark metrics/plots were regenerated against
  the current adaptive quality guard and evolutionary-event schema.
- Input validation now rejects empty matrices, empty reference FASTA files, and
  duplicate genome or cluster names before running expensive searches.
- GFF3 percent-encoded attributes are decoded before gene/alias/product matching.
- The integration test now checks the current output schema version.

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
- Synteny validation now caches reference locus contexts per cluster/window
  during a run, avoiding repeated GFF context reconstruction for every rescue.
- Internal stop-codon scanning now groups rescues by genome, avoiding repeated
  full rescue-list scans after each FASTA load.
- Output reports are sorted deterministically for reproducible downstream review.
- TSV report fields are normalized to single-line values so annotation-derived
  tabs or line breaks cannot corrupt tabular outputs.
- Large pipeline and rescue implementation files were split into smaller
  section files, separating tests, reference audits, evolutionary summaries,
  and output writers into real Rust modules without changing runtime
  behaviour.

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
- Adaptive synteny scoring now matches each reference neighbour at most once,
  avoiding inflated scores in tandem/repeated flanks.
- Recurrent synteny support in large reference sets now requires compatible
  contexts in at least two genomes, so multiple copies from one genome cannot
  validate a rescue by themselves.
- NSC representative selection is now deterministic on prevalence ties, and
  NSC component reporting is stable for both MinHash and mmseqs2 engines.
- Duplicate rescue events for the same cluster/genome are now collapsed by an
  explicit evidence ranking instead of relying on hash-map overwrite order.
- Gene ID suffix normalization now strips only Panaroo-style `_dup_<number>`
  and `_pseudo` suffixes, avoiding false alias collapse for biological names
  that contain `_dup_`.
- Equal-scoring fragment-chain candidates now use deterministic locus
  tie-breakers, avoiding arbitrary fragment evidence in repeated regions.

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
