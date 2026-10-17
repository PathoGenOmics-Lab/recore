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

Each cell is one cluster in one genome, in one of six states:

| State | Colour | Meaning |
|-------|--------|---------|
| input | grey | present in the input matrix |
| corrected | blue | absent in the input; Recore found the gene and corrected the cell |
| NSC | orange | absent in the input; near-identical clusters were merged and the merged cluster is present (only when NSC ran) |
| kept out | crimson | Recore found sequence evidence, but the matrix policy, a quality floor or a synteny or paralog guard kept it out of the corrected matrix |
| absent | pale | absent; Recore searched the gap and found no evidence |
| not searched | pale, hatched | absent in the input and never searched: the cluster is below `--min-prevalence` or has no sequence to search with, or the genome has no FASTA |

Input presence is the quiet background on purpose, so what Recore found is what
draws the eye at any zoom. The five colours stay distinguishable under deutan and
protan colour vision, and not-searched cells add a hatch rather than a sixth
colour. The legend in the matrix header only lists the states the run has, and
drops to swatches when the pane is too narrow for the labels; each swatch still
names its state on hover.

**Absent or not searched.** An absent cell is a gap Recore searched and found
nothing in. A hatched cell is a gap it never searched, so its absence is the
input's call and not a finding. Under the default `--min-prevalence` of 0.1 every
cluster in fewer than one genome in ten is hatched. Its tooltip and cell card say
why it was not searched, the cluster card says so for a whole cluster, and in the
tree a tip with a dashed ring is a genome the selected cluster was not searched in.

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
reference-absent filter, the searched-clusters filter and the mechanism filter.
**Searched clusters only** hides the clusters whose gaps were never searched and
refits the columns, leaving the part of the matrix Recore actually audited.

**Column order.** Input matrix order, prevalence, mechanism, number of corrected
cells, or cluster name.

**Mechanism filter.** The matrix menu lists the mechanisms in view as a
checklist, with each one's description in its tip. Ticking one shows only its
clusters and refits the columns to the pane; active filters appear as chips
beside the matrix count, and a click on a chip removes it. **Download TSV**
saves the clusters in view with their counts of input, corrected, NSC,
kept-out, absent and not-searched cells.

**Corrections by lineage.** For MTBC-style genome names the cluster card shows,
per lineage, how many genomes were corrected or had evidence kept out, out of
the lineage's genomes.

**Zoom.** The slider in each pane header drives its own axis, in equal ratios
from end to end. A plain wheel always scrolls. ++cmd++ (or ++ctrl++) and the
wheel zoom columns, which is also what a trackpad pinch does; ++alt++ and the
wheel zoom rows. Each notch is a factor of 1.25 and the cell under the pointer
stays under the pointer. Press ++z++ to fit the data back to the window.

Press ++z++ to see the whole matrix at once. When that makes cells smaller than
a pixel, the matrix is drawn as density: each pixel shows the share of present
cells it covers, between the absent and present colours, and is tinted blue
where any cell in it was corrected (crimson where evidence was kept out), more
strongly the larger that share. The hatch carries on in density, as strong as
the share of not-searched cells in the pixel. Zoom in on a region to get back to
single cells.

Above 40 columns a strip under the matrix shows the mean prevalence of the
clusters at each position, with a blue mark where any of them has a changed
cell, and the current viewport outlined. Clicking it jumps.

The open genome's row stays outlined in the tree and the matrix after the
pointer moves on, so a genome found from the find field is still easy to spot.

## Comparing two genomes

Two genomes far apart in the tree are far apart in the matrix. **Pin to
compare** in the genome card, or ++p++, keeps a genome's row in a band under the
column header, aligned with every column whatever the scroll; its name sits
level with it in the tree pane's header, marked **A** or **B**, and its tip gets
a blue ring. Pinning a third genome unpins the oldest. Hover a cell of the band
for its state; click it to open that cell.

With two genomes pinned, the genome card compares them over every cluster of the
corrected matrix: how many are in both, only in A, only in B and in neither,
and the clusters that differ, each one a click away. When a difference rests on
a gap that was never searched, the card says so, because then the input matrix,
not Recore, says the gene is missing. **Show only the clusters that differ**
filters the matrix to them and refits the columns; a chip beside the matrix
count takes the filter off again.

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

## Run details

++d++ or the info button opens the run details. They start with what you need
to write the run up:

- **Using these results**: a methods paragraph built from the run's own
  settings and counts, with a button that copies it; the caveats the evidence
  calls for (cells kept out of the matrix, corrections whose span is not whole
  codons, rejected candidates, synteny or the tree being off); and the
  parameters, the command line and the time the report was written.
- **Changed cells**: every evidence row, corrections first, with a button that
  copies them as tab-separated text and one that downloads
  `recore_changed_cells.tsv`. Clicking a row opens that cell in the matrix.

Below come the run statistics, the dominant mechanisms, the top candidate
clusters (the first ten, then **Show all**) and every file in the output
directory with its row count. The links open the files beside the report, so
they work while the report sits in its output directory.

The genome card shows each genome's QC status and flags from `genome_qc.tsv`,
and tips can be coloured by QC status.

## Getting figures out

The download button writes `recore_matrix.png` and, when the synteny band is
open, `recore_synteny.svg`.

The PNG is drawn for the page, not copied from the screen: the tree, the lineage
strip, the column names and the matrix share one geometry, so every tip sits on
its row, at three times the on-screen size, always in the light palette and
without the hover crosshair or the selection marks. Below the matrix come a key
for the cell states, the tree's scale bar and a caption with what is shown, the
root, the matrix policy and the Recore version. When the whole matrix fits in
4,000 by 4,000 screen pixels it is drawn in full; otherwise the window on screen
is drawn and the caption names its rows and columns.

The SVG carries every colour and font as a plain attribute, with no stylesheet
or CSS variables, so illustration software and PDF tools draw it as the page
does.

++cmd+p++ prints the figure described above for the whole matrix (in the light
palette, whatever the screen is in), then the run details as a document with
every candidate. The live panes are swapped for the figure only while printing.

## Keyboard shortcuts

Press ++question++ for this list in the report itself. A switch at its top turns
the one-letter shortcuts off, for example when they fire by accident with speech
input; ++cmd+k++ and the arrows keep working, and the choice is remembered in
that browser.

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
| ++p++ | pin or unpin the open genome, to compare two |
| ++cmd++ / ++ctrl++ + wheel, or pinch | zoom columns under the pointer |
| ++alt++ + wheel | zoom rows under the pointer |
| ++shift++ + wheel | scroll sideways |
| ++alt++ + click a tip | root the tree on that genome |
| ++esc++ | close whatever is open |

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

`benchmarks/recore_report.py` builds a compact Markdown and HTML report
(`recore_report.md` and `recore_compact_report.html`) from an existing output
directory. It refuses to replace an interactive `recore_report.html` unless
given `--force`.
