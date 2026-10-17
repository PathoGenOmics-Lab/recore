# Run report

`--write-html-report` writes `recore_report.html` into the output directory. It
is a single file with no external requests: it opens from a results directory on
a cluster node with no network, and it keeps working when you email it to a
co-author.

```bash
recore ... --tree species.tree --write-html-report -o recore_output/
```

[Open an example report](../assets/example-report.html){ .md-button .md-button--primary target=_blank }

That one is the output of the [tutorial](../getting-started/tutorial.md) run:
twelve genomes, twenty-four clusters, eleven planted false absences. It is the
whole file, opened straight from this page, which is the point of it being
self-contained.

Without `--tree` the workspace still opens: the rows follow the matrix order, the
tree pane only lists the genomes, and a notice above the workspace says so.
Ancestral readings, tree colouring by subtree and rerooting need a real tree.

If the workspace cannot start (a browser that cannot hold the embedded data, for
instance), the report says why in place of the workspace, and ++d++ still opens
the run details.

## The workspace

The viewer *is* the report. Four regions share the window:

| Region | Shows |
|--------|-------|
| Phylogeny | the matched tips, in the order the matrix rows use |
| Corrected matrix | one column per cluster, one row per genome, coloured by evidence state or by mechanism |
| Inspector | the selected cluster, or the selected genome |
| Unique syntenies | the flanking-gene context for the selected cluster, full width along the bottom |

Every seam between them is draggable, and the sizes are remembered in the
browser for the next report you open.

## Finding things

Press ++cmd+k++ (or ++ctrl+k++, or ++slash++) to open the find field. It searches
clusters **and** genomes at once, ranking prefix matches on the name above name
matches above annotation matches. Enter opens the result: a cluster becomes the
selected column, a genome scrolls the panes to its tip and opens its card.

The arrow buttons beside the field step through the visible columns one at a
time, and ++left++ / ++right++ do the same from the keyboard. ++shift++ with them
jumps ten.

## The tree

The sliders icon in the phylogeny header opens everything that shapes the tree.

**Shape.** Phylogram draws branch lengths to scale; cladogram ignores them and
shows topology only. A scale bar under the pane gives the phylogram its units.

**Root.** Every ancestral reading in the report depends on where the tree is
rooted, so it is an explicit choice: the L8 outgroup when one is present, the
tree as supplied, or the midpoint of the longest tip-to-tip path. ++alt++ and a
click on any tip roots on that genome instead.

**Tip colour.** By the selected cluster's state, by MTBC lineage, by mechanism,
by rescue burden, by absence burden, or by any correction. For MTBC-style names
such as `L4_AS00146`, lineage colouring uses the mycolorsTB palette and colours
monophyletic branches as well as tips.

A plain click on a tip opens that genome in the inspector: its lineage, its
state in the selected cluster, its counts across all clusters, and its distance
from the root.

The burden modes draw a ramp from the line colour to the state colour, with a
hollow ring for zero, and a key under the top bar gives the ramp's maximum.

## The matrix

Each cell is one cluster in one genome, in one of five states:

| State | Colour | Meaning |
|-------|--------|---------|
| input | grey | present in the input matrix |
| corrected | blue | absent in the input; Recore found the gene and corrected the cell |
| NSC | orange | absent in the input; near-identical clusters were merged and the merged cluster is present (only when NSC ran) |
| kept out | crimson | Recore found sequence evidence, but the matrix policy, a quality floor or a synteny or paralog guard kept it out of the corrected matrix |
| absent | pale | absent, with no evidence |

Input presence is the quiet background on purpose, so what Recore found is what
draws the eye at any zoom. The five colours stay distinguishable under deutan and
protan colour vision. The legend in the matrix header only lists the states the
run has.

**Cell evidence.** Hover a cell for its genome, cluster and state; for a
corrected or kept-out cell the tooltip adds the rescue type, identity and
coverage. Click a cell for the full row from `orthology_evidence.tsv` at the top
of the inspector: the locus, the synteny and ORF checks, whether it went into
the matrix and why not, and its publication status. The cluster card lists every
evidence row of the cluster, corrections first; clicking one opens that cell and
brings its tip into view.

**Mechanism colour.** In mechanism mode the input cells take the colour of their
cluster's mechanism family (retained, gain, loss, disruption, paralog or context
change, artifact or unresolved), while corrected and kept-out cells keep their
own colours. The legend key opens the full list, with counts and filters.

The sliders icon in the matrix header holds cell colour, column order, the
reference-absent filter and the mechanism filter.

**Column order.** Input matrix order, prevalence, mechanism, number of corrected
cells, or cluster name.

**Mechanism filter.** The mechanism cards double as filters: clicking one adds or
removes that mechanism from the visible columns. The header button stays lit
while any filter is on, so a filter left running two minutes ago does not
silently explain a matrix that looks wrong.

**Zoom.** The slider in each pane header drives its own axis. ++shift++ and the
wheel zooms columns; ++cmd++ and the wheel zooms rows. Press ++z++ to fit the
data back to the window.

Above 40 columns a prevalence strip appears under the matrix: one bar per pixel
showing how prevalent the clusters at that position are, with the current
viewport marked. At 30,000 columns the scrollbar tells you where you are and
nothing about what is there; the strip tells you what is there, and clicking it
jumps.

## Synteny

Press ++s++ or use the layers button to show the synteny band. It draws the top
unique five-neighbour blocks for the selected cluster, at most
10 per cluster, with gene arrows keeping strand and proportional length.

Uniqueness is defined by the ordered displayed gene architecture, ignoring
internal cluster gaps, coordinate differences and length differences. Each row
shows one representative locus.

- Selecting a row marks its carrier tips on the tree.
- Row 2 is the inferred reference context from shared neighbours, and it is
  clickable.
- Clicking any gene arrow jumps to that cluster in the matrix.
- Shared neighbouring clusters are joined only between adjacent rows.

This needs `cluster_synteny_blocks.tsv`, which means GFFs were loaded.

## Getting figures out

The download button writes `recore_matrix.png` (the tree, lineage strip, column
header and matrix composited at 2x) and `recore_synteny.svg` (vector, so the gene
arrows stay crisp at print resolution).

++cmd+p++ prints the whole thing: the workspace as page one, then the run details
as a document. The canvases are redrawn for the print layout rather than being
clipped to whatever the screen happened to be.

## Keyboard shortcuts

Press ++question++ for this list in the report itself.

| Keys | Action |
|------|--------|
| ++cmd+k++ / ++ctrl+k++ / ++slash++ | find a cluster or genome |
| ++left++ ++right++ | previous / next cluster |
| ++shift++ + ++left++ ++right++ | jump ten clusters |
| ++s++ | show or hide the synteny panel |
| ++d++ | open the run details |
| ++e++ | export the view |
| ++f++ | fullscreen |
| ++t++ | switch theme |
| ++c++ | phylogram or cladogram |
| ++r++ | reroot at the midpoint |
| ++z++ | fit the matrix to the window |
| ++shift++ + wheel | zoom columns |
| ++cmd++ + wheel | zoom rows |
| ++alt++ + click a tip | root the tree on that genome |
| ++esc++ | close whatever is open |

## Run details

The info button opens the drawer: the run summary, the QC breakdowns, the top
candidate clusters, and links to the exact TSV files behind every number.
Clicking a candidate takes you to it in the matrix.

## Scale

The matrix ships as compact per-cluster state sets rather than one character per
cell, and each pane draws only its visible window onto a viewport-sized canvas.
A 10,000 genome by 30,000 cluster matrix, 300 million cells, is a 17 MB file that
reaches interactive in about 130 ms and re-renders in single-digit milliseconds
at every zoom level.

The workspace is embedded up to **20,000 matched tips** and **1.2 billion cells**.
Above either limit the report says so and points at
`gene_presence_absence_recore_states.csv` for the full matrix, rather than
shipping something a browser cannot open.

## If you did not use `--write-html-report`

`benchmarks/recore_report.py` builds a compact Markdown or HTML report from an
existing output directory.
