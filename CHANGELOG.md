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
- The embedded viewer can now color tree tips by MTBC lineage using the
  mycolorsTB Nature palette, and places the interactive viewer controls next to
  the tree/matrix viewport.
- Evolutionary mechanism cards in the embedded viewer now support multi-select
  column filtering, and selected-cluster details report the reference gene/cell
  or reference-derived event source when available.
- The embedded viewer now uses `cluster_synteny_blocks.tsv` to add a split
  synteny panel for the selected gene cluster, drawing ranked unique
  five-neighbour blocks with the focal cluster highlighted and shared neighbours
  connected between adjacent rows.
- Synteny blocks now include carrier genomes, representative gene coordinates,
  strand, length, gene-level block records, and an inferred reference/H37Rv
  context from shared neighbours; selecting a synteny block marks its carrier
  tips on the embedded tree.
- The synteny panel can now be hidden from the embedded viewer controls so the
  tree plus matrix expand to the full available width during dense inspection.
- The embedded synteny panel is hidden by default, includes matrix/synteny
  split and synteny x-axis zoom sliders, draws proportional gene arrows, lets
  the H37Rv/reference context row mark the reference tip, and lets gene-arrow
  clicks jump back to the corresponding matrix column.
- `cluster_synteny_blocks.tsv` now collapses synteny blocks by ordered displayed
  gene/locus identity, ignoring internal cluster gaps and representative
  length/coordinate differences, while keeping one representative locus and all
  carrier genomes per unique architecture.
- The embedded viewer now colors monophyletic tree branches by lineage when
  lineage coloring is selected, and places the matrix/synteny split slider as
  a vertical control between the two panes.
- Output manifest schema version is now `recore.outputs.v6`.
- README now documents the current v6 audit guarantees, input guardrails,
  multicopy-corrected matrix metadata, and evolutionary interpretation fields.
- The HTML report is now a full-viewport workspace built around the phylogeny,
  the corrected matrix and synteny, with the run summary, breakdowns, candidate
  table and file list moved into a run-details drawer. Synteny is a full-width
  band under the matrix, and the zoom controls sit on the axis they drive.
- The displayed tree can be rerooted: on the L8 outgroup, as supplied, at the
  midpoint, or on any genome by clicking its tip.
- Chrome controls are icons with immediate tooltips rather than a row of text
  buttons.
- Tree tips are aligned row for row with the matrix, and the top bar carries the
  version with links to the documentation and the repository.
- The workspace can order matrix columns by prevalence, mechanism, corrected
  cells or name instead of only the input matrix order, and can find a genome by
  name, ringing the matching tips and scrolling to the first hit.
- The report can export the current view as a PNG, and the synteny panel as SVG.
- The embedded matrix ships as compact per-cluster state sets rather than one
  character per cell, and each pane renders only its visible window, so the
  viewer handles runs of 10,000 genomes by 30,000 clusters. The tip limit rises
  from 5,000 to 20,000.
- The report chrome is one row. Everything that shapes a single pane moved into
  that pane's own menu, which leaves identity, finding, the run readout and the
  global actions on top and gives the data back the second row.
- Finding a cluster is a palette on `Cmd-K` or `/` rather than a `<select>`, and
  it answers for genomes as well, ranking name prefixes over name hits over
  annotation hits. At 30,000 clusters the old control rebuilt 30,000 option
  elements on every keystroke.
- The run figures are a readout rather than a row of pills, and they give way in
  a deliberate order as the window narrows instead of widening the page.
- Clicking a tip opens a genome card in the inspector instead of rerooting the
  tree; rerooting moved to Alt-click and to a button on the card. The inspector
  is tabbed, so a genome question no longer evicts the selected cluster.
- The phylogeny has a scale bar, and the matrix has a whole-matrix prevalence
  strip for runs above 40 columns.
- Matrix column names are drawn across the column once it is wide enough to hold
  them, and the opening zoom is fitted to the pane that exists rather than to a
  fixed target height.
- The synteny band opens at the height of its own content until the reader drags
  the seam.
- Keyboard shortcuts throughout, with a reference sheet on `?`.
- The report logo is cropped to the wordmark and stripped of the opaque plate the
  export wrapped it in, so it fills the box it is given rather than two thirds
  of it.

### Added

- The documentation site takes the design the other PathoGenOmics tools use: a
  neutral masthead with a hairline instead of a saturated band, one accent colour
  that only ever means "this is a link or a control", hairlines instead of boxes
  on code, tables and callouts, headings that carry the structure, and the four
  nav groups always open as small-caps landmarks. IBM Plex, breadcrumbs, and edit
  and view actions on every page.
- Dark-background variants of the wordmark and the benchmark figure. Three of
  the mark's six inks sit between 1.8:1 and 2.3:1 against GitHub's dark canvas,
  under the 3:1 WCAG asks of a non-text graphic, and the figure was a white slab.
  The README picks between them with `<picture>` and `prefers-color-scheme`; the
  documentation keys the swap on its own theme toggle instead, because that is
  what the reader actually set; and the report, which inlines the mark, repaints
  its inks in CSS rather than carrying a second copy of a 30 KB asset.
- `tools/make_logo_variants.py` derives all of that from one source file: it
  crops the export to the artwork, strips the editor metadata and the fixed size
  and `preserveAspectRatio="none"` that let any CSS box stretch the mark, lifts
  each ink to a contrast floor with its hue fixed, and writes the report's
  override rules between markers so they cannot drift from the dark file.

- Documentation site built with MkDocs Material, published from `docs/` to
  GitHub Pages. The four reference documents are now the deep half of a site
  with a landing page, installation, quick start, inputs, matrix policies, the
  run report, evidence states, limits, citation and contributing.
- `docs/getting-started/tutorial.md`, a full run on a dataset with known answers.
  `benchmarks/make_tutorial_dataset.py` generates that dataset deterministically:
  twelve MTBC-style genomes, twenty-four clusters, GFF3 annotations, a Newick
  tree, and eleven false absences planted on purpose.

### Fixed

- The HTML report no longer requests a remote webfont stylesheet, so it renders
  identically on offline and air-gapped machines and makes no network request
  when opened.
- The report logo keeps its aspect ratio instead of being stretched into a
  square, and the theme is resolved before the first paint so dark-mode readers
  no longer get a white flash.
- The embedded viewer draws its matrix, column header, and lineage strip at the
  device pixel ratio, so cells and selection outlines are sharp on HiDPI
  screens.
- The embedded tree caption moved to the pane header; it used to be drawn over
  the first tip labels.
- The embedded viewer now opens at a row/column zoom fitted to the tip and
  cluster counts, sizes its panes to the data instead of always reserving the
  maximum height, and labels matrix columns with the cluster name when the
  columns are wide enough.
- Showing the synteny panel no longer stretches the whole report to the width of
  the synteny drawing and gives the page a horizontal scrollbar.
- Candidate table headers are keyboard-operable and report `aria-sort`, section
  anchors clear the sticky navigation bar, and transitions honour
  `prefers-reduced-motion`.
- The report prints correctly: the light palette is forced regardless of the
  selected theme, bar and matrix colours survive printing, and the desktop
  layout is kept instead of collapsing to the narrow-screen layout on paper.
- Bar breakdowns report each row's share of the panel total, so a full-width bar
  is no longer ambiguous between "most clusters" and "one of twelve".
- Mouse-wheel column zoom in the embedded viewer no longer snaps the zoom back
  down when it is above the old ceiling.
- Leaving fullscreen wrote text into the fullscreen button, which replaced its
  icon with a blank box.
- The chrome could grow wider than the window, which took the workspace with it:
  the shell's single implicit grid column sized itself to the widest child.
- The workspace slid into the wrong grid row whenever the lineage legend was
  absent, which is whenever tips are not coloured by lineage, and the free
  vertical space went to the status bar.
- The reference audit block rendered as a bare heading with an indented body,
  because its `dt` and `dd` sat outside the definition list.
- The run-details drawer and the shortcut sheet can be closed with a button, not
  only with `Esc`.
- The tree scale bar chose its step from a quarter of the tree depth, which in a
  narrow pane rounded to a bar just under the display threshold, so no bar
  appeared at all.
- The wordmark's counters are real holes now rather than white paint, so the
  mark sits on any background without a plate behind it.
- The lineage strip beside the tree was blanked on every frame: `render()` called
  its draw function a second time with no arguments, after `renderTree()` had
  already drawn it with the correct window, which sized its canvas to 0 by 0.
- Without `--tree` the topbar drew the find field, cluster stepping, synteny,
  export, fullscreen and the shortcut sheet, none of which have handlers in a
  report that does not embed the viewer. They are no longer written, and a test
  asserts it.
- Clicking a candidate row in the run-details drawer did nothing: it looked for
  the `<select>` the redesign removed. The drawer now asks the viewer by cluster
  name and closes only if the viewer answers.
- The narrow-screen breakpoint targeted the pane container by its pre-rewrite
  class name, so below 820px it set the tree pane to full width inside a still
  horizontal row and pushed the matrix off screen.
- Removed 300 lines of stylesheet whose selectors no longer match anything the
  report emits, four rule blocks that were declared twice, a `:fullscreen` rule
  aimed at a class the fullscreen element does not carry, and the in-page
  navigation observer left over from the document layout.
- Six report assertions were matched by strings inside the inlined stylesheet
  and script rather than by the markup they claimed to test, so they could not
  fail. They are anchored on emitted markup now.
- `icon()` leaked a `String` per call to satisfy a `&'static str` signature.
- The matrix and its column header each carried two scroll listeners: an early
  unguarded pair that bounced `scrollLeft` between the two elements, and the
  guarded pair that actually keeps them in step. The unguarded pair is gone.
- `treeDirty` was assigned in fourteen places and read in none. It was the
  cache-invalidation flag of the renderer that preceded the virtual one.
- Documentation: the output-files page still described the pre-rewrite synteny
  layout and a `summary.txt` without four of the lines the writer always emits;
  the contributing page listed a benchmark command CI does not run; the
  START_VARIANT confidence formula included a coverage factor the code does not
  apply.
- `cargo fmt --check` and `cargo clippy -D warnings` pass again.
- CI installs minimap2 and fails when it is missing, so the end-to-end test
  actually runs instead of skipping silently. Its assertions were updated to the
  current adaptive-guard behaviour: only evidence accepted into the binary matrix
  reaches a corrected GFF, and audit-only evidence stays out of it.
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
