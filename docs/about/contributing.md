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

```bash
pip install -r requirements-docs.txt
mkdocs serve
```

CI builds with `mkdocs build --strict`, so a broken internal link or a page
missing from the nav fails the build rather than shipping a dead link.

## Reporting a problem

Open an issue at
[github.com/PathoGenOmics-Lab/recore/issues](https://github.com/PathoGenOmics-Lab/recore/issues).
For a wrong call, the most useful thing to include is the row from
`orthology_evidence.tsv`: it carries the evidence state, the adaptive guard
floors, the confidence tier and the risk flags, which together explain why
Recore decided what it did.
