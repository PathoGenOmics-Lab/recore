<div align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/recore_logo_dark.svg">
    <img src="assets/recore_logo.svg" alt="Recore" width="450">
  </picture>

**Post-hoc rescue, audit, and quality control of pangenome presence/absence matrices.**
**Pure Rust · one external dependency · Linux and macOS**

[Documentation](docs/) · [Install](#installation) · [Quick start](#quick-start) · [Tutorial](docs/getting-started/tutorial.md) · [Citation](#citation)

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

- a **Panaroo-compatible corrected binary matrix** for downstream analyses
- a **multistate evidence matrix** for audit, plotting, and manuscript review

The central design choice is conservative. Recore separates *sequence evidence*
from *clean binary presence*: full-length evidence can be corrected in the binary
matrix, while partial homology, pseudogenes, paralog traps and
incompatible-context hits stay visible in audit outputs instead of silently
becoming presence calls.

## Installation

Needs [Rust 1.85+](https://www.rust-lang.org/tools/install) to build and
[minimap2](https://github.com/lh3/minimap2) on `PATH` to run.

```bash
git clone https://github.com/PathoGenOmics-Lab/recore.git
cd recore
cargo build --release        # binary at target/release/recore
```

Tagged releases also ship Linux and macOS binaries on the
[Releases](https://github.com/PathoGenOmics-Lab/recore/releases) page.
[Installation](docs/getting-started/installation.md) covers conda, MMseqs2 and
checking the install.

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

Every page is a markdown file in this repository and reads here on GitHub.

| Start here | |
|---|---|
| [Tutorial](docs/getting-started/tutorial.md) | a full run on a generated dataset with known answers |
| [Installation](docs/getting-started/installation.md) | requirements, from source, pre-built binaries |
| [Quick start](docs/getting-started/quickstart.md) | the run most people want, and what to read first |
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

> [!NOTE]
> The same pages build into a searchable site with an example report you can
> open (`pip install -r requirements-docs.txt && mkdocs serve`). It is deployed
> to `gh-pages` on every push to `main`, but GitHub Pages needs the repository
> to be public, or the organisation to be on a paid plan, so
> <https://pathogenomics-lab.github.io/recore/> is not served yet.

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
