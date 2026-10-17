# Contributing to Recore

Thank you for looking. Recore is a small tool with a narrow job, so most useful
contributions are small too: a parser that chokes on a real-world file, a wrong
call you can reproduce, a paragraph of documentation that lied to you.

## Where the real guide lives

There are two contributing documents in this repository and they are not
duplicates:

- **[`docs/about/contributing.md`](../docs/about/contributing.md)** is the
  development guide. It is part of the documentation site, it is the one that
  gets updated when the build changes, and it is the source of truth for how to
  run the checks, work on the HTML report and regenerate the example report.
- **This file** is the GitHub-facing front door. GitHub links it from the "New
  issue" and "New pull request" screens, which is the moment somebody actually
  reads it, so it covers the social side (what to open, what to expect) and
  repeats only the handful of commands you need before pressing the button.

If the two ever disagree, the docs page wins and this file is the bug.

## Table of contents

- [What Recore is, and what it is not](#what-recore-is-and-what-it-is-not)
- [Ways to contribute](#ways-to-contribute)
- [The four checks](#the-four-checks)
- [The two report invariants](#the-two-report-invariants)
- [Where things live](#where-things-live)
- [Rust version, and one trap worth naming](#rust-version-and-one-trap-worth-naming)
- [Python in this repository](#python-in-this-repository)
- [Documentation](#documentation)
- [Pull requests](#pull-requests)
- [Licence and code of conduct](#licence-and-code-of-conduct)

## What Recore is, and what it is not

Recore is a post-hoc rescue, audit and quality-control layer for pangenome
presence/absence matrices. It takes a matrix that some other pipeline already
produced (Panaroo-style), searches every absent cell against the genome
assemblies with minimap2, classifies the evidence it finds, validates orthology
against annotation and synteny when those are available, and writes a corrected
binary matrix next to an audit-grade multistate matrix and evidence tables.

Knowing the boundary saves everyone time on issues and pull requests:

- Recore is **not a pangenome builder**. It does not cluster genes from scratch,
  it does not rebuild a graph, and it will not fix a matrix whose clustering was
  wrong to begin with. It reannotates cells that are already in the matrix.
- Recore is **not an oracle**. Homologous sequence is not proof of a functional
  gene, which is exactly why the default recommendation keeps partial
  homologues, pseudogenes and paralog traps out of the clean binary matrix and
  visible in the audit outputs instead. A patch that makes the binary matrix
  more generous needs an argument about false positives, not only about recall.
- Recore **cannot recover what the assembly does not contain**. If the sequence
  is not there, no policy will find it.

[`docs/how-it-works/limits.md`](../docs/how-it-works/limits.md) is the current,
maintained statement of scope. Read it before proposing a feature that widens
the boundary; the answer may already be there.

## Ways to contribute

**Report a problem.** Open an issue at
[github.com/PathoGenOmics-Lab/recore/issues](https://github.com/PathoGenOmics-Lab/recore/issues).
For a wrong call, the single most useful attachment is the offending row from
`orthology_evidence.tsv`. That row carries the evidence state, the adaptive
identity and coverage floors, the confidence tier and the risk flags, which
together are Recore's own explanation of why it decided what it did. Without it
a maintainer is guessing; with it the diagnosis is usually one read.

Also include the Recore version (`recore --version`), the exact command line,
the minimap2 version, and the operating system. "It gave the wrong answer" with
no command line is not reproducible, and an issue nobody can reproduce ends up
closed rather than fixed.

**Ask a question or float an idea.** Discussions are enabled on this
repository. A question about interpreting an evidence state, or a "would it make
sense if Recore did X" that you are not sure is a feature request yet, belongs
there rather than in the issue tracker. Issues are for things with a definite
outcome: this is broken, or this should exist.

**Request a feature.** Open an issue and describe the analysis you are trying to
finish, not only the flag you want. Recore has a lot of surface already, and the
most common good outcome of a feature request is discovering that an existing
parameter combination does the job, which only happens if the goal is on the
table.

**Improve the documentation.** Documentation fixes are contributions and are
reviewed like any other change. If a page misled you, that is a bug report with
a known fix.

## The four checks

These are exactly what continuous integration runs, and running them locally
first is the difference between one review round and three. The first three are
the CI workflow, in the order it runs them; the fourth belongs to the
documentation workflow and only fires when the website or the code it describes
changed:

```bash
cargo fmt --all -- --check
cargo clippy --locked --all-targets -- -D warnings
RECORE_REQUIRE_MINIMAP2=1 cargo test --locked --all-targets
python3 tools/check_site.py --recore target/debug/recore --run-tutorial  # docs/ or the CLI
```

Three things about that list are worth explaining, because they are not the
defaults you would guess:

- **Clippy runs with `-D warnings`.** A warning fails the build. This is
  deliberate: a warning that nobody has to fix is a warning everybody learns to
  scroll past, and then the one that mattered scrolls past too.
- **`RECORE_REQUIRE_MINIMAP2=1` is not optional.** The end-to-end test needs
  minimap2 on `PATH`, and without it the test prints `SKIPPED` and passes. That
  is the right behaviour on a laptop that does not have minimap2 installed, and
  exactly the wrong behaviour in CI, where a green run would then mean "the
  pipeline was never executed". The environment variable turns the skip into a
  failure. CI sets it; set it locally too, otherwise you are testing less than
  you think you are.
- **`--locked` is on purpose.** It makes the build fail rather than quietly
  update `Cargo.lock`, so a dependency bump is always a visible, reviewable
  change instead of a side effect of somebody's test run.

The benchmarks under `benchmarks/` are **not** part of CI, because they need
data and minutes that a pull request check should not spend. Run them by hand
when a change could move the numbers, and put the before and after in the pull
request.

## The two report invariants

The HTML report is a self-contained page written by
`src/pipeline/output.rs` from assets that live in
`src/pipeline/report_assets/` (`report.css`, `explorer.js`, `report-ui.js`) and
are inlined into the binary at compile time with `include_str!`. Roughly 3300
lines of that are JavaScript: the phylogeny, matrix and synteny viewer. Two
rules hold that design together, and both are easy to break by accident.

**1. No remote assets.** The report has to open from a results directory on a
machine with no network: a cluster login node, an air-gapped analysis box, a
reviewer's laptop on a plane, a supplementary file attached to a manuscript.
Every font, style, script and image is embedded. The integration test walks the
resource-loading elements of a generated report and fails if any of them points
at `http://` or `https://`. A hyperlink to the project is fine; *loading* a
resource over the network is not. If you add a CDN link because it is faster to
develop that way, the test will catch it, which is the point.

**2. Nothing may assume the whole matrix fits.** Every pane draws only its
visible window onto a viewport-sized canvas. Code that iterates all cells once
per frame, or builds one DOM node per cell, will look perfectly fine on the tiny
test fixture and will be unusable at 30,000 clusters by several hundred genomes,
which is an ordinary real input. The same rule applies on the Rust side: the
presence/absence storage is deliberately compact and indexed rather than a
hash set per cluster. When you change the viewer, say in the pull request what
size input you reasoned about.

One related habit: every value that reaches the report from an input file goes
through `html_escape` or `json_escape` in `src/pipeline/output.rs`. Cluster
names, genome names and annotation text come from files Recore did not write,
and they end up inside a page that a human opens in a browser. New fields must
go through the same helpers.

To look at a change, generate a report and serve it over HTTP rather than
opening the file directly. Some browsers refuse to run the inline script from a
`file://` URL and you get an empty workspace with no error message, which costs
an hour if you do not know it:

```bash
recore ... --write-html-report -o out/
python3 -m http.server 8000 --directory out/
```

## Where things live

| Path | What is in it |
|------|---------------|
| `src/` | the whole tool: `cli.rs`, `parser.rs`, `search.rs`, `rescue/`, `nscluster.rs`, `relationships.rs`, `pipeline/` |
| `src/pipeline/report_assets/` | the report stylesheet and viewer JavaScript, inlined with `include_str!` |
| `src/rescue/tests.rs` | the largest unit suite: rescue classification |
| `src/pipeline/tests.rs` | pipeline, output and report unit tests |
| `tests/integration.rs` | the end-to-end run, the one that needs minimap2 |
| `benchmarks/` | Python benchmarks and the tutorial dataset generator |
| `tools/` | `make_logo_variants.py`, which generates the logo variants and the report's dark-mode colour overrides |
| `docs/` | the website (plain HTML, no build step) and the long-form markdown pages |
| `.github/workflows/` | CI, documentation build, release builds |

Recore is a Rust command-line tool and nothing else. There is no desktop
application, no web service and no npm build step. The JavaScript in this
repository is the report viewer, and it is edited as plain files with no
bundler; the Python is benchmarks and asset generation, run by hand.

## Rust version, and one trap worth naming

`Cargo.toml` declares `rust-version = "1.85"`, the oldest compiler the
dependency set builds on (clap 4.6 sets it). The `msrv` job in `ci.yml` builds
the crate with that toolchain and checks that `Cargo.toml` still names it, so
an accidental use of a newer standard-library API fails the pull request.

The job passes the version through the action's `toolchain:` input:

```yaml
- uses: dtolnay/rust-toolchain@master
  with:
    toolchain: "1.85"
```

Never encode it in the action's git ref as `dtolnay/rust-toolchain@1.85`.
Dependabot treats that ref as a version to keep current, bumps it like any other
action pin, and from then on the job that is supposed to prove Recore builds on
1.85 quietly proves it builds on something else. The failure mode is a green
check that means nothing, which is worse than no check at all. When the minimum
moves, change `rust-version`, the job's `toolchain:` and its check together.

## Python in this repository

The Python is tooling, not a shipped package, and it is not covered by CI:

- `benchmarks/synthetic_benchmark.py` builds known false absences, true
  absences and partial pseudogenes, then computes precision, recall and F1 from
  the corrected matrices that Recore actually wrote.
- `benchmarks/biological_event_benchmark.py` and `benchmarks/self_benchmark.py`
  cover the evolutionary-event labelling and internal consistency.
- `benchmarks/compare_runs.py` diffs two runs, which is the fastest way to see
  what a change did to real output.
- `benchmarks/make_tutorial_dataset.py` generates the dataset used by the
  tutorial and by the example report.
- `benchmarks/amap_performance_benchmark.py` profiles a run on a full local
  dataset; `benchmarks/README.md` documents its invocation and its output.
- `benchmarks/recore_report.py` and `tools/make_logo_variants.py` generate
  assets.

Keep them dependency-light and runnable with a plain `python3`. They are
diagnostics, and a diagnostic that needs its own environment before it can tell
you anything tends not to get run.

## Documentation

The website is plain HTML in `docs/` (`index.html`, `reference.html` and
`assets/`), with no build step: open `docs/index.html` to see a change. The
documentation workflow runs `tools/check_site.py`, which fails on a dead link or
anchor, an option `recore --help` does not have, an option missing from the
reference, a default that differs from the binary's, or an example run whose
counts differ from the page's, instead of a stale page somebody finds six
months later. The long-form pages are markdown under `docs/` and read on GitHub.

`docs/assets/example-report.html` is the tutorial run's own output, linked from
the [Run report](../docs/guide/report.md) page. If you change the viewer, that
file is stale until you regenerate it; the exact commands are in
[`docs/about/contributing.md`](../docs/about/contributing.md).

## Pull requests

1. Work on a branch, not on `main`, and open a pull request even for your own
   changes. CI runs on both pushes and pull requests, but the pull request is
   where the reasoning is recorded, and in six months the reasoning is the part
   you will want.
2. Keep one change per pull request. A fix bundled with a refactor is hard to
   review and impossible to revert cleanly.
3. Fill in the pull request template. The section that matters is **how it was
   verified**: not "it should work", but what you ran and what it printed. That
   is the one thing a reviewer cannot reconstruct from the diff.
4. Commit messages: a short subject line, then two or three sentences of
   explanation. Say why, not what; the diff already says what.
5. Add an entry under `## [Unreleased]` in `CHANGELOG.md` for anything a user
   would notice. Behaviour changes to the matrix policies, the evidence states
   or the report output always count.
6. New behaviour needs a test, or a sentence in the pull request saying why it
   cannot have one. "Tested by hand" is an acceptable answer for a viewer
   change; it is not one for a rescue-classification change.

## Licence and code of conduct

By contributing you agree that your contribution is licensed under the same
terms as the project, which are the ones in [LICENSE](../LICENSE) at the
repository root.

> [!WARNING]
> The repository currently disagrees with itself about which licence that is.
> `Cargo.toml`, the README and the documentation all say GNU General Public
> License v3.0 (`GPL-3.0-only`), while the committed `LICENSE` file contains the
> text of the GNU **Affero** General Public License v3. Those are different
> licences with different obligations for anyone running a modified copy over a
> network. This needs resolving by the maintainer, one way or the other, before
> the first outside contribution lands.

Participation is covered by the [Code of Conduct](CODE_OF_CONDUCT.md).
Security issues do not go in the issue tracker: see [SECURITY.md](SECURITY.md).

Maintainer: Paula Ruiz-Rodriguez, paula.ruiz.rodriguez@csic.es
