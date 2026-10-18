# Security Policy

## Reporting a vulnerability

Email **paula.ruiz.rodriguez@csic.es** with `recore security` in the subject
line. Please do not open an issue for a security problem: an issue is visible to
everybody who can read this repository and is indexed the moment it exists,
which hands the details to anyone watching before there is anything to update
to.

A useful report contains:

- what an attacker can make Recore do, stated as an outcome (write outside the
  output directory, execute a command, run a script in the reader's browser,
  exhaust memory on an input that looks ordinary);
- the Recore version (`recore --version`) or commit, the operating system, and
  the minimap2 version;
- the exact command line;
- the smallest input that triggers it. A three-record FASTA that reproduces the
  bug is worth more than a real dataset that also reproduces it, and it avoids
  sending data that may not be yours to send. Please do not attach patient or
  otherwise identifiable sequence data to a security report.

You should get a first reply within about a week. This is a small academic
project with no on-call rotation, so that is an honest expectation rather than a
service level; if a week passes with silence, send the email again, because the
most likely explanation is that it was buried, not ignored.

Please give the fix a reasonable window before publishing details. If you tell
us you intend to disclose on a given date, we will work to it and will credit
you in the release notes unless you would rather stay anonymous.

## Supported versions

| Version | Status |
|---------|--------|
| `1.0.x` | supported, fixes land here first |
| older than 1.0 | not supported: upgrade to the newest release |

There is deliberately no backport branch. A security fix goes onto `main` and
out in the next release, and the answer to "which version should I run" is the
newest one.

Nothing else is supported. In particular there is no Bioconda package, no
container image and no third-party build maintained by this project, so a binary
you did not build yourself or download from this repository's
[Releases](https://github.com/PathoGenOmics-Lab/recore/releases) page is outside
what we can speak for.

## What counts as a vulnerability here

Recore is a command-line tool that reads files, shells out to minimap2 (and
optionally MMseqs2), and writes a self-contained HTML report that a human then
opens in a browser. That shape decides what is interesting.

**In scope:**

- Anything reachable from a **malformed or hostile input file**. Recore parses
  Panaroo-style CSV, FASTA, GFF3, Newick and minimap2 PAF, and in practice those
  files arrive from collaborators, public archives and other pipelines. Treat
  every one of them as untrusted input. A panic on garbage is a bug and worth
  reporting; a panic that is really an out-of-bounds write, an unbounded
  allocation driven by a length field, or a path escaping the output directory
  is a vulnerability.
- **Injection into the HTML report.** Cluster names, genome names and
  annotation text come from files Recore did not write and end up inside a page
  someone opens in a browser. Values are escaped on the way in by `html_escape`
  and `json_escape` in `src/pipeline/output.rs`; if you find a field that gets
  past them and can run script, that is a real finding, because the report is
  routinely emailed around and attached to manuscripts.
- **Argument or path handling when invoking external tools.** Recore builds
  minimap2 and MMseqs2 command lines from user-supplied paths. A file name that
  turns into an extra argument, or into a command, is in scope.
- **Writing outside the output directory**, including through symlinks or
  crafted identifiers used to build output file names.
- **Temporary file handling**, for example a predictable path in a shared
  temporary directory that another user on the same machine can pre-create.

**Out of scope:**

- A wrong biological call. Recore rescuing a cell it should not have, or missing
  one it should have found, is a correctness bug: open a normal issue and attach
  the `orthology_evidence.tsv` row.
- Vulnerabilities in minimap2, MMseqs2 or the Rust crates themselves. Report
  those upstream. If a crate advisory affects Recore, an issue here is welcome
  so the dependency can be bumped.
- Denial of service that only comes from genuinely large legitimate input.
  Recore is expected to work hard on a real pangenome; slow is not a
  vulnerability, though a pathological blow-up on a small file is.
- Anything that requires an attacker who already runs code as your user.

## Running Recore safely

- Run it as an ordinary user. Recore never needs elevated privileges, and it
  writes only into the output directory you give it.
- Point `-o` at a fresh directory you own rather than a shared scratch path that
  other users can write to.
- Treat inputs from outside your group as untrusted, the same way you would
  treat any file you did not produce.
- Keep dependencies current. `cargo update` plus the automated dependency
  updates configured for this repository are how crate advisories reach a build.
- The HTML report loads nothing from the network by design and can be opened
  offline. If you serve it, serve it from a directory you control.
