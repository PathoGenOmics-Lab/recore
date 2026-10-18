#!/usr/bin/env python3
"""Mark the example report as made from synthetic data.

docs/assets/example-report.html is the tutorial run's own output, and the
tutorial's genomes are made up. After copying a regenerated report into place,
run this once: it puts a "Synthetic example data" note next to the report's
title, where it shows in the website's frame, full screen and on a phone.
Running it again replaces the note instead of adding a second one.

    python3 tools/label_example_report.py docs/assets/example-report.html
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

MARKER = "data-example-note"
NOTE = (
    "<span class='synthetic-note' data-example-note "
    "title='Made by the tutorial run from twelve made-up MTBC-style genomes with eleven planted absences'>"
    "<span class='synthetic-note-long'>Synthetic example data</span><span class='synthetic-note-short'>Synthetic</span></span>"
)
# Below 900 px the report's top bar has no room to spare, so the note shortens
# and the search box may narrow to make up the rest.
STYLE = (
    "<style data-example-note>"
    ".synthetic-note{flex:0 0 auto;display:inline-flex;align-items:center;height:24px;padding:0 9px;"
    "border-radius:999px;border:1px solid color-mix(in srgb,#e0a030 60%,transparent);"
    "background:color-mix(in srgb,#e0a030 16%,transparent);color:var(--ink);"
    "font:600 11.5px/1 var(--font-sans);white-space:nowrap}"
    ".synthetic-note-short{display:none}"
    "@media (max-width:900px){.synthetic-note-long{display:none}.synthetic-note-short{display:inline}"
    ".topbar-find{min-width:100px}}"
    "</style>"
)
TITLE = re.compile(r"<span class='topbar-title'>[^<]*</span>")
OLD = re.compile(
    r"<style data-example-note>.*?</style>"
    r"|<span class='synthetic-note' data-example-note[^>]*>(?:<span class='(?:synthetic-note-long|synthetic-note-short|long|short)'>[^<]*</span>)+</span>"
    r"|<span class='synthetic-note' data-example-note[^>]*>[^<]*</span>",
    re.S,
)


def label(html: str) -> str:
    """The report with exactly one note and its style, whatever an older run of
    this script left in it."""
    html = OLD.sub("", html)
    match = TITLE.search(html)
    if not match or "</head>" not in html:
        raise SystemExit("ERROR: no topbar title to put the note next to; is this a Recore report?")
    html = html[: match.end()] + NOTE + html[match.end():]
    return html.replace("</head>", STYLE + "</head>", 1)


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    path = Path(sys.argv[1])
    html = path.read_text(encoding="utf-8")
    labelled = label(html)
    if labelled != html:
        path.write_text(labelled, encoding="utf-8")
        print(f"{path}: marked as synthetic example data")
    else:
        print(f"{path}: already marked")
    return 0


if __name__ == "__main__":
    sys.exit(main())
