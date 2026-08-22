<!--
Fill in what applies and delete what does not. A short pull request needs a
short description; nobody is asking for an essay to fix a typo.
-->

## What this changes

<!-- One or two sentences. What is different after this is merged? -->

## Why

<!-- The problem it solves. Link an issue with "Closes #123" if there is one. -->

## How it was verified

<!--
The important part, and the one a reviewer cannot reconstruct from the diff.

Not "it should work", but what you ran and what it printed. For example: the
test you added and the fact that it fails without the fix, the matrix you ran
through it and the cells whose state changed, the `orthology_evidence.tsv` rows
before and after, or the benchmark numbers from
`python3 benchmarks/synthetic_benchmark.py` on both sides of the change.

If CI already covers it, say which check covers it.
-->

## Scale it was reasoned about

<!--
Only for changes to the pipeline or to the report viewer; delete otherwise.

The test fixture is tiny and everything looks fast on it. Say what size input
you had in mind: number of clusters, number of genomes, and whether anything
here iterates the whole matrix rather than the visible window.
-->

## Checklist

- [ ] `cargo fmt --all -- --check` and `cargo clippy --locked --all-targets -- -D warnings` are clean.
- [ ] `RECORE_REQUIRE_MINIMAP2=1 cargo test --locked --all-targets` passes. <!-- The variable turns the "minimap2 not in PATH" skip into a failure, so a green run means the end-to-end test actually ran. If you do not have minimap2 locally, say so and let CI be the check. -->
- [ ] New behaviour has a test, or I have said above why it does not.
- [ ] `CHANGELOG.md` has an entry under `## [Unreleased]`, if a user would notice this.
- [ ] Documentation under `docs/` is updated if the change is user-visible, and `mkdocs build --strict` passes.
- [ ] If I touched the report: it still loads nothing from the network, it still draws only the visible window, and I regenerated `docs/assets/example-report.html` if the viewer changed.
