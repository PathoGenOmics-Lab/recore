# Contributing

## Development checks

Exactly what CI runs, in the order it runs, is three commands:

```bash
cargo fmt --all -- --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --all-targets
```

The end-to-end test needs minimap2 on `PATH`. It skips silently without one,
which is not what you want in CI, so CI exports the guard that turns a skip into
a failure, and you should too:

```bash
RECORE_REQUIRE_MINIMAP2=1 cargo test --locked --all-targets
```

The benchmarks are not part of CI. Run them by hand when a change could move the
numbers:

```bash
python3 benchmarks/synthetic_benchmark.py
```

## Working on the HTML report

The report's assets are inlined into the binary at compile time from
`src/pipeline/report_assets/`:

| File | Holds |
|------|-------|
| `report.css` | the whole stylesheet |
| `explorer.js` | the phylogeny, matrix and synteny viewer |
| `report-ui.js` | theme, drawer, sortable tables, tooltips |

Two rules the report has to keep:

1. **No remote assets.** It has to open from a results directory with no
   network. The integration test walks every `link`, `script`, `img`, `iframe`
   and `use` element and fails if any of them loads from a URL. Links to the
   project are fine; loading a resource is not.
2. **Nothing may assume the whole matrix fits.** Every pane draws only its
   visible window onto a viewport-sized canvas. A change that iterates all cells
   per frame will be fine on the test fixture and unusable at 30,000 clusters.

To look at a change, generate a report and open it over HTTP. Opening the file
directly with `file://` works, but some browsers block the inline script and you
get an empty workspace with no error.

```bash
recore ... --write-html-report -o out/
python3 -m http.server 8000 --directory out/
```

## Documentation

The website is `docs/index.html`, `docs/how-it-works.html` and
`docs/reference.html`, styled by `docs/assets/site.css` and driven by
`docs/assets/site.js`: plain files, so opening `docs/index.html` in a browser
shows a change. Keep the front page short. It is there to get someone running in
two minutes; the details belong in the documentation.

The drawings of How it works are `docs/assets/how.js`. They decide with a copy
of `evidence_state` and `matrix_policy_corrects` from `src/pipeline.rs`, so a
change to those rules is a change there too; `tools/check_site.py` fails if the
policy table, a state, a rescue type, a matrix action or a column no longer
matches the code.

The Docs tab is these markdown pages. `tools/build_docs.py` turns each one into
an `.html` page next to it, with the site's header, sidebar, contents and
search, and writes `docs/documentation.html` and the search index. The pages
are committed, so after editing a `.md` file, rebuild them:

```bash
python3 tools/build_docs.py
```

A new page is a line in `PAGES` at the top of that script.

Check it the way CI does:

```bash
cargo build
python3 tools/check_site.py --recore target/debug/recore --run-tutorial
```

That fails on a documentation page that no longer matches its markdown, a
dead link or anchor, an option the site mentions that
`recore --help` does not have, an option missing from the reference, a default
that differs from the binary's, an em-dash, and an example run whose counts
differ from the ones the page shows.

`docs/assets/example-report.html` is the report the website and the
[Run report](../guide/report.md) page link to. It is the tutorial run's own output, so regenerate it whenever the
viewer changes:

```bash
python3 benchmarks/make_tutorial_dataset.py --out tutorial_data
recore -p tutorial_data/gene_presence_absence.csv \
       -r tutorial_data/pan_genome_reference.fa \
       -g tutorial_data/genomes -a tutorial_data/annotations \
       --synteny --multi-allele --matrix-policy balanced \
       --reference-genome L4_AS00401 --tree tutorial_data/species.tree \
       --write-html-report -o tutorial_out/
cp tutorial_out/recore_report.html docs/assets/example-report.html
```

The two logo variants and the report's dark-mode colour overrides all come from
one source file, so regenerate those together rather than editing either by hand:

```bash
python3 tools/make_logo_variants.py
```

## Reporting a problem

Open an issue at
[github.com/PathoGenOmics-Lab/recore/issues](https://github.com/PathoGenOmics-Lab/recore/issues).
For a wrong call, the most useful thing to include is the row from
`orthology_evidence.tsv`: it carries the evidence state, the adaptive guard
floors, the confidence tier and the risk flags, which together explain why
Recore decided what it did.
