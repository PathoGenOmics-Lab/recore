#!/usr/bin/env python3
"""Check the website in docs/ against the code it describes.

The site is plain HTML with no build step, so nothing renders it and nothing
fails when it drifts. This does, in four parts:

  links    every internal link and #anchor resolves, on the front page, the
           reference and every documentation page, and no page loads a
           stylesheet, script or image from another host (the site has to open
           from a clone, offline);
  docs     the documentation pages are what tools/build_docs.py makes of the
           markdown today, so an edited .md cannot ship with a stale page;
  flags    every option the pages or the command builder mention exists in
           `recore --help`, the reference lists every option there is, and
           its defaults are the ones the binary prints;
  dashes   no em-dashes;
  try      with --run-tutorial, the "Try it" commands are run as written and
           summary.txt must give the numbers the page promises.

Run from the repository root:

    python3 tools/check_site.py --recore target/debug/recore --run-tutorial

Exits 1 and lists every problem it found.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from html.parser import HTMLParser
from pathlib import Path

SITE = Path("docs")
PAGES = ["index.html", "reference.html"]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_docs  # noqa: E402  (tools/build_docs.py, next to this file)
SCRIPTS = [SITE / "assets" / "site.js", SITE / "assets" / "site.css"]
DATASET_SCRIPT = Path("benchmarks/make_tutorial_dataset.py")


class Page(HTMLParser):
    """The parts of a page the checks read."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.ids: set[str] = set()
        self.links: list[tuple[str, str]] = []  # (tag, url)
        self.remote_assets: list[str] = []
        self.code_text: list[str] = []
        self.code_elems: list[tuple[str, str]] = []  # (tag, text) of each <code> and <pre>
        self.blocks: dict[str, str] = {}
        self.option_rows: list[list[str]] = []
        self._code_depth = 0
        self._code_open: list[tuple[str, list[str]]] = []
        self._capture: list[tuple[str, int, list[str]]] = []
        self._depth = 0
        self._row: list[str] | None = None
        self._cell: list[str] | None = None
        self._in_options = False

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        self._depth += 1
        if "id" in a:
            self.ids.add(a["id"])
        for key in ("href", "src"):
            if key in a and a[key] is not None:
                self.links.append((tag, a[key]))
        remote = lambda url: url and re.match(r"^(https?:)?//", url)
        if tag == "link" and "stylesheet" in (a.get("rel") or "") and remote(a.get("href")):
            self.remote_assets.append(a["href"])
        if tag in ("script", "img", "iframe") and remote(a.get("src")):
            self.remote_assets.append(a["src"])
        if tag in ("code", "pre"):
            self._code_depth += 1
            self._code_open.append((tag, []))
        if tag == "pre" and a.get("data-copy-source"):
            self._capture.append((a["data-copy-source"], self._depth, []))
        if tag == "div" and "out" in (a.get("class") or "").split():
            self._capture.append(("out", self._depth, []))
        if tag == "table" and a.get("id") == "options-table":
            self._in_options = True
        if self._in_options and tag == "tr":
            self._row = []
        if self._in_options and tag in ("td", "th") and self._row is not None:
            self._cell = []
        if tag == "br" and self._capture:
            self._capture[-1][2].append("\n")

    def handle_endtag(self, tag):
        if tag in ("code", "pre"):
            self._code_depth = max(0, self._code_depth - 1)
            if self._code_open and self._code_open[-1][0] == tag:
                kind, parts = self._code_open.pop()
                self.code_elems.append((kind, "".join(parts)))
        if self._capture and self._capture[-1][1] == self._depth and tag in ("pre", "div"):
            name, _, parts = self._capture.pop()
            self.blocks[name] = "".join(parts)
        if self._in_options and tag in ("td", "th") and self._cell is not None and self._row is not None:
            self._row.append("".join(self._cell).strip())
            self._cell = None
        if self._in_options and tag == "tr" and self._row is not None:
            self.option_rows.append(self._row)
            self._row = None
        if tag == "table":
            self._in_options = False
        self._depth -= 1

    def handle_data(self, data):
        if self._code_depth:
            self.code_text.append(data)
            for _, parts in self._code_open:
                parts.append(data)
        for _, _, parts in self._capture:
            parts.append(data)
        if self._cell is not None:
            self._cell.append(data)


def parse(path: Path) -> Page:
    page = Page()
    page.feed(path.read_text(encoding="utf-8"))
    return page


def help_options(recore: str) -> tuple[set[str], set[str], dict[str, str]]:
    """Long options, short options and defaults that `recore --help` prints."""
    out = subprocess.run([recore, "--help"], check=True, capture_output=True, text=True).stdout
    longs = set(re.findall(r"(?<![\w-])(--[a-z][a-z0-9-]*)", out))
    shorts = set(re.findall(r"^\s+(-[a-zA-Z]),", out, flags=re.M))
    defaults: dict[str, str] = {}
    current = None
    for line in out.splitlines():
        m = re.match(r"^\s+(?:-[a-zA-Z], )?(--[a-z][a-z0-9-]*)", line)
        if m:
            current = m.group(1)
            continue
        d = re.search(r"\[default: ([^\]]+)\]", line)
        if d and current:
            defaults[current] = d.group(1)
    return longs, shorts, defaults


def all_pages() -> dict[str, Page]:
    """Every page of the site, by its path under docs/; the example report is
    the report's own output and is checked by the report's tests."""
    out = {}
    for path in sorted(SITE.rglob("*.html")):
        rel = path.relative_to(SITE).as_posix()
        if rel != "assets/example-report.html":
            out[rel] = parse(path)
    return out


def check_links(pages: dict[str, Page]) -> list[str]:
    problems = []
    for name, page in pages.items():
        for url in page.remote_assets:
            problems.append(f"{name}: loads {url} from another host; the site must work offline")
        for tag, url in page.links:
            if re.match(r"^(https?:|mailto:)", url):
                continue
            target, _, anchor = url.partition("#")
            if not target:
                if anchor and anchor not in page.ids:
                    problems.append(f"{name}: link to #{anchor}, which no element on the page has")
                continue
            resolved = os.path.normpath(os.path.join(os.path.dirname(name), target)).replace(os.sep, "/")
            if resolved.startswith("..") or not (SITE / resolved).exists():
                problems.append(f"{name}: <{tag}> points at {target}, which does not exist")
                continue
            if anchor and resolved in pages and anchor not in pages[resolved].ids:
                problems.append(f"{name}: link to {target}#{anchor}, which no element there has")
    return problems


def check_docs_built() -> list[str]:
    stale = [
        path.as_posix()
        for path, text in build_docs.build().items()
        if not path.exists() or path.read_text(encoding="utf-8") != text
    ]
    return [f"{p}: does not match its markdown; run python3 tools/build_docs.py" for p in stale]


def check_flags(pages: dict[str, Page], recore: str) -> list[str]:
    problems = []
    longs, shorts, defaults = help_options(recore)
    mentioned_long: dict[str, str] = {}
    mentioned_short: dict[str, str] = {}

    def note(text: str, where: str) -> None:
        for flag in re.findall(r"(?<![\w-])(--[a-z][a-z0-9-]*)", text):
            mentioned_long.setdefault(flag, where)
        for flag in re.findall(r"(?<![\w./-])(-[a-zA-Z])(?![\w-])", text):
            mentioned_short.setdefault(flag, where)

    # Only what is recore's: a command line that starts with recore, or a
    # snippet that is an option on its own. `cargo install --git`, `conda
    # install -c` and `make_tutorial_dataset.py --out` are other programs'.
    for name, page in pages.items():
        for tag, text in page.code_elems:
            if tag == "pre":
                for command in text.replace("\\\n", " ").splitlines():
                    command = command.strip()
                    if command.startswith("$ "):
                        command = command[2:]
                    if command.split()[:1] == ["recore"]:
                        note(command, name)
            elif text.startswith("recore ") or text.startswith("-"):
                note(text, name)
    # The command builder's lines are written as line(["-p", tool.p]).
    js = (SITE / "assets" / "site.js").read_text(encoding="utf-8")
    for args in re.findall(r"line\(\[([^\]]*)\]", js):
        for flag in re.findall(r"\"(-{1,2}[a-zA-Z][a-z0-9-]*)\"", args):
            note(flag, "site.js")
    for flag, where in sorted(mentioned_long.items()):
        if flag not in longs:
            problems.append(f"{where}: mentions {flag}, which recore --help does not have")
    for flag, where in sorted(mentioned_short.items()):
        if flag not in shorts:
            problems.append(f"{where}: mentions {flag}, which recore --help does not have")

    reference = pages["reference.html"]
    listed = set(re.findall(r"(--[a-z][a-z0-9-]*)", " ".join(reference.code_text)))
    for flag in sorted(longs - listed):
        problems.append(f"reference.html: does not list {flag}")
    for row in reference.option_rows:
        if len(row) != 3:
            continue
        option, default, _ = row
        flag = next(iter(re.findall(r"(--[a-z][a-z0-9-]*)", option)), None)
        if flag in defaults and default != defaults[flag]:
            problems.append(
                f"reference.html: {flag} default reads {default!r}, recore --help says {defaults[flag]!r}"
            )
    return problems


def check_dashes(pages: dict[str, Page]) -> list[str]:
    problems = []
    for path in [SITE / p for p in pages] + SCRIPTS:
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if "\u2014" in line:
                problems.append(f"{path}:{number}: em-dash")
    return problems


def check_tutorial(pages: dict[str, Page], recore: str) -> list[str]:
    """Run the Try-it commands as the page shows them and compare the counts."""
    index = pages["index.html"]
    script = index.blocks.get("try")
    promised = index.blocks.get("out")
    if not script or not promised:
        return ["index.html: no Try-it block or no expected output block to check"]
    want = {
        key: re.search(rf"^{key}:\s+(\d+)", promised, flags=re.M)
        for key in ("Rescued", "Matrix-corrected")
    }
    if not all(want.values()):
        return ["index.html: the expected output does not show Rescued and Matrix-corrected counts"]
    commands = []
    for line in script.replace("\\\n", " ").splitlines():
        line = line.strip()
        if line.startswith("$ "):
            line = line[2:]
        if not line:
            continue
        if line.startswith("curl ") and DATASET_SCRIPT.name in line:
            commands.append(["cp", str(DATASET_SCRIPT.resolve()), "."])
            continue
        words = line.split()
        if words[0] == "recore":
            words[0] = str(Path(recore).resolve())
        commands.append(words)
    with tempfile.TemporaryDirectory() as work:
        for words in commands:
            run = subprocess.run(words, cwd=work, capture_output=True, text=True)
            if run.returncode != 0:
                return [f"Try it: `{' '.join(words)}` failed:\n{run.stderr[-2000:]}"]
        summary = Path(work, "tutorial_out", "summary.txt")
        if not summary.exists():
            return ["Try it: the commands ran but wrote no tutorial_out/summary.txt"]
        text = summary.read_text()
    problems = []
    for key, match in want.items():
        got = re.search(rf"^{key}:\s+(\d+)", text, flags=re.M)
        if not got or got.group(1) != match.group(1):
            problems.append(
                f"Try it: the page promises {key}: {match.group(1)}, the run gave "
                f"{got.group(1) if got else 'nothing'}"
            )
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--recore", default="target/debug/recore", help="the recore binary to read --help from")
    parser.add_argument("--run-tutorial", action="store_true", help="also run the Try-it commands (needs minimap2)")
    args = parser.parse_args()

    pages = all_pages()
    front = {name: pages[name] for name in PAGES}
    problems = check_docs_built() + check_links(pages) + check_flags(front, args.recore) + check_dashes(pages)
    if args.run_tutorial:
        if not shutil.which("minimap2"):
            problems.append("Try it: minimap2 is not on PATH")
        else:
            problems += check_tutorial(front, args.recore)
    if problems:
        print("The website does not match the code:", file=sys.stderr)
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        return 1
    checked = "generated pages, links, options, defaults, dashes" + (", and the Try-it run" if args.run_tutorial else "")
    print(f"docs/: {checked} all agree with the code")
    return 0


if __name__ == "__main__":
    sys.exit(main())
