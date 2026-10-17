# Installation

## Requirements

| Requirement | Version | Why |
|-------------|---------|-----|
| [Rust](https://www.rust-lang.org/tools/install) | 1.70 or newer | building from source |
| [minimap2](https://github.com/lh3/minimap2) | any recent release, on `PATH` | every sequence search |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | optional | alignment-based nucleotide similarity clustering with `--nsc-engine mmseqs2` |

minimap2 is not optional. Recore shells out to it for the search, and a
run fails immediately if it is not on `PATH`.

## From source

```bash
git clone https://github.com/PathoGenOmics-Lab/recore.git
cd recore
cargo build --release
```

The binary lands at `target/release/recore`. Put it somewhere on your `PATH`, or
call it by path:

```bash
target/release/recore --help
```

## Pre-built binaries

Tagged releases build Linux and macOS binaries through GitHub Actions. Download
them from [Releases](https://github.com/PathoGenOmics-Lab/recore/releases).

## With conda for the dependencies

Recore itself is a static binary, but conda is a convenient way to get minimap2
and MMseqs2 without touching the system:

```bash
conda create -n recore -c bioconda -c conda-forge minimap2 mmseqs2
conda activate recore
```

## Check the install

```bash
recore --version
minimap2 --version
```

If `--nsc-engine mmseqs2` is going to be used, check that too:

```bash
mmseqs version
```

Asking for `--nsc-engine mmseqs2` without MMseqs2 on `PATH` is an error, not a
silent downgrade. The default engine is `minhash`, which is built in and needs
nothing installed.

## Building the documentation

This site is MkDocs Material. To preview it while editing:

```bash
pip install -r requirements-docs.txt
mkdocs serve
```

`mkdocs build --strict` is what CI runs, and it turns a broken internal link into
a failed build.
