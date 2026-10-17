<div align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/recore_logo_dark.svg">
    <img src="assets/recore_logo.svg" alt="Recore" width="450">
  </picture>

**Post-hoc rescue, audit, and quality control of pangenome presence/absence matrices.**
**Pure Rust · one external dependency · Linux and macOS**

[Documentation](#documentation) · [Install](#installation) · [Quick start](#quick-start) · [Tutorial](docs/getting-started/tutorial.md) · [Citation](#citation)

<a href="https://github.com/PathoGenOmics-Lab/recore/actions/workflows/ci.yml"><img src="https://github.com/PathoGenOmics-Lab/recore/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
<a href="LICENSE"><img src="https://img.shields.io/badge/license-GPLv3-blue.svg" alt="License: GPLv3"></a>
<img src="https://img.shields.io/badge/Rust-1.85%2B-orange.svg" alt="Rust 1.85+">

</div>

## What is Recore?

Pangenome pipelines can mark a gene as absent even when the sequence is still
present in the assembly. These false absences inflate the accessory genome,
shrink the core genome, and create lineage-specific signals that are artefacts
of graph cleaning rather than biology.

Recore searches every absent cell against the genome assemblies with minimap2,
classifies the evidence, validates orthology against annotation and synteny when
they are available, and writes both:

- a **Panaroo-compatible corrected binary matrix** for downstream analyses,
  as CSV and as the 1/0 `.Rtab` that pyseer and most R code read
- a **multistate evidence matrix** for audit, plotting, and manuscript review

It reads the output of Panaroo, Roary, PPanGGOLiN, PanTA and PIRATE, can write
the rescued genes as nucleotide and protein FASTA for functional annotation,
and with `--search-cache` reruns with another policy or threshold without
aligning again.

The central design choice is conservative. Recore separates *sequence evidence*
from *clean binary presence*: full-length evidence can be corrected in the binary
matrix, while partial homology, pseudogenes, paralog traps and
incompatible-context hits stay visible in audit outputs instead of silently
becoming presence calls.

## Installation

Needs [Rust 1.85+](https://www.rust-lang.org/tools/install) to build and
[minimap2](https://github.com/lh3/minimap2) on `PATH` to run.

```bash
conda install -c bioconda minimap2     # or: brew install minimap2 / sudo apt install minimap2
cargo install --git https://github.com/PathoGenOmics-Lab/recore --locked
recore --version
```

Or build a clone with `cargo build --release` (the binary lands at
`target/release/recore`). `--nsc-engine mmseqs2` also needs
[MMseqs2](https://github.com/soedinglab/MMseqs2) on `PATH`; the default engine
is built in.

## Quick start

```bash
recore \
  -p gene_presence_absence.csv \
  -r pan_genome_reference.fa \
  -g genomes/ \
  -a annotations/ \
  --synteny --multi-allele \
  --matrix-policy balanced \
  --reference-genome AS00146 \
  --tree species.tree \
  --write-html-report \
  -t 8 -o recore_output/
```

Four outputs answer most questions:

| File | For |
|---|---|
| `gene_presence_absence_recore.csv` | the corrected binary matrix, for downstream tools |
| `gene_presence_absence_recore.Rtab` | the same matrix as 1/0, for pyseer and Rtab readers |
| `gene_presence_absence_recore_states.csv` | audit and figures, one biological state per cell |
| `orthology_evidence.tsv` | publication wording and manual review |
| `recore_report.html` | the whole run in one self-contained page |

`--matrix-policy balanced` is the recommendation for anything that will be
published: it rescues full-length false absences while keeping partial or
disrupted homologues as audit-only evidence. The other two policies, and the
adaptive guard that sits on top of all three, are in
[Matrix policies](docs/guide/matrix-policies.md).

To see it work on data with known answers before pointing it at your own, the
[tutorial](docs/getting-started/tutorial.md) generates a twelve-genome dataset
with eleven false absences planted in it and walks through the run.

## Documentation

The quickest way in is the website in [`docs/`](docs/index.html): an install
line, a command builder for your pangenome tool, a one-minute example, the
example report to click through, a reference page with every option, output
file and cell state, and a Docs tab with every page below, searchable. It is plain HTML with no build step, so open
`docs/index.html` from a clone; once the repository is public, GitHub Pages can
serve it straight from `main` and `/docs`. The Docs pages are built from the
markdown by `tools/build_docs.py`, and `tools/check_site.py` keeps the site
honest in CI: the pages match their markdown, every option it shows exists, the
defaults match `recore --help`, and the example run gives the numbers the page
promises.

The same pages are markdown files in this repository and read here on GitHub.

| Start here | |
|---|---|
| [Tutorial](docs/getting-started/tutorial.md) | a full run on a generated dataset with known answers |
| [Inputs](docs/getting-started/inputs.md) | what each input is, and how the names have to match |

| Reference | |
|---|---|
| [Parameters](docs/guide/parameters.md) | every option, with its default and what it changes |
| [Matrix policies](docs/guide/matrix-policies.md) | which evidence reaches the binary matrix |
| [Output files](docs/guide/outputs.md) | every table, column by column |
| [Run report](docs/guide/report.md) | the self-contained HTML viewer |

| How it works | |
|---|---|
| [Algorithm](docs/how-it-works/algorithm.md) | each pipeline phase, its data structures and complexity |
| [Evidence states](docs/how-it-works/evidence-states.md) | the ten states of the multistate matrix |
| [Rescue types](docs/how-it-works/rescue-types.md) | rescue classes and biological interpretation |
| [Evolutionary events](docs/how-it-works/evolutionary-events.md) | how a mechanism is called, and how to word it in a paper |
| [Limits and comparisons](docs/how-it-works/limits.md) | what Recore does not claim, benchmark results, Panaroo refinding |

Version history is in [CHANGELOG.md](CHANGELOG.md). Benchmark design and
reproduction is in [benchmarks/README.md](benchmarks/README.md).

If the documentation does not answer it, ask in
[Discussions](https://github.com/PathoGenOmics-Lab/recore/discussions). Issues
are for something that crashed, refused to run, or is plainly wrong.

## For developers

```bash
cargo fmt --all -- --check
cargo clippy --locked --all-targets -- -D warnings
RECORE_REQUIRE_MINIMAP2=1 cargo test --locked --all-targets
```

The end-to-end test skips silently without minimap2 on `PATH`, which is why CI
sets that variable. [Contributing](docs/about/contributing.md) covers the report
assets, the two rules they have to keep, and how to regenerate the generated
files.

## Citation

If you use Recore, please cite:

> Ruiz-Rodriguez P, et al. (2026). Recore: post-hoc rescue and quality control
> of pangenome presence/absence matrices. *In preparation.*

## Contact

Paula Ruiz-Rodriguez - paula.ruiz.rodriguez@csic.es

## License

GNU General Public License v3.0. See [LICENSE](LICENSE).
