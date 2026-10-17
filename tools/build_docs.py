#!/usr/bin/env python3
"""Build the website's documentation pages from the markdown in docs/.

The markdown stays the source: GitHub renders it, and this writes an .html page
next to each file with the site's header, theme, sidebar, page contents and
search, plus docs/documentation.html (the Docs tab) and the search index. The
pages are committed, so the site keeps opening from a clone with no server and
no build step; --check fails when they no longer match the markdown, which is
what CI runs.

    python3 tools/build_docs.py           # write the pages
    python3 tools/build_docs.py --check   # exit 1 if any page is stale

Only the markdown these pages use is supported: headings, paragraphs, lists
(nested, loose or tight, with code inside), tables, block quotes, fenced code
(a `title="..."` info string becomes a caption, `mermaid` flowcharts become a
list of steps), rules, inline code, emphasis, links and a few inline tags.
Anything else shows up as text, and the link check in tools/check_site.py
catches what that breaks.
"""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
from pathlib import Path

DOCS = Path("docs")
REPO = "https://github.com/PathoGenOmics-Lab/recore"

# The Docs tab, in reading order: (section, markdown path, title, one line).
PAGES = [
    ("Getting started", "getting-started/tutorial.md", "Tutorial", "The example run, number by number."),
    ("Getting started", "getting-started/inputs.md", "Inputs", "What each tool writes, and how names match."),
    ("Guide", "guide/parameters.md", "Parameters", "Every option, with what it changes."),
    ("Guide", "guide/matrix-policies.md", "Matrix policies", "The three policies and the quality guard."),
    ("Guide", "guide/outputs.md", "Output files", "Every table, column by column."),
    ("Guide", "guide/report.md", "The report", "Every pane and control of the viewer."),
    ("How it works", "how-it-works/algorithm.md", "Algorithm", "Each step of the pipeline."),
    ("How it works", "how-it-works/evidence-states.md", "Evidence states", "What each cell of the states matrix means."),
    ("How it works", "how-it-works/rescue-types.md", "Rescue types", "What Recore saw where a gene was found."),
    ("How it works", "how-it-works/evolutionary-events.md", "Evolutionary events", "How an event is called, and how to word it."),
    ("How it works", "how-it-works/limits.md", "Limits and speed", "What Recore does not claim; run time and memory."),
    ("About", "about/citation.md", "Citation", "How to cite Recore."),
    ("About", "about/contributing.md", "Contributing", "Working on the code, the report and these pages."),
]
INDEX = "documentation.html"
SEARCH_INDEX = DOCS / "assets" / "search-index.js"

GH_PATH = (
    "M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94"
    "-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07"
    "-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0"
    " 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25"
    ".54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"
)
GH_SVG = f'<svg viewBox="0 0 16 16" fill="currentColor" aria-hidden="true"><path d="{GH_PATH}"/></svg>'
SUN_SVG = (
    '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" aria-hidden="true">'
    '<circle cx="12" cy="12" r="4.5"/><path d="M12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1'
    'l1.4-1.4M17.7 6.3l1.4-1.4"/></svg>'
)
THEME_BOOT = (
    '<script>try { document.documentElement.dataset.theme = localStorage.getItem("recore.theme") === "light" ? '
    '"light" : "dark"; } catch (e) { document.documentElement.dataset.theme = "dark"; }</script>'
)
INLINE_TAGS = re.compile(r"</?(kbd|br|sup|sub|strong|em|b|i)\s*/?>", re.I)


# ── Markdown ──


def slugify(text: str, seen: dict[str, int]) -> str:
    """GitHub's anchor for a heading, so links written for GitHub land here too."""
    slug = re.sub(r"[^\w\- ]", "", text.strip().lower()).replace(" ", "-")
    if slug in seen:
        seen[slug] += 1
        slug = f"{slug}-{seen[slug]}"
    else:
        seen[slug] = 0
    return slug


def plain(markdown_inline: str) -> str:
    """The text of an inline span, without its markup."""
    text = re.sub(r"`([^`]*)`", r"\1", markdown_inline)
    text = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", text)
    text = INLINE_TAGS.sub("", text)
    # The emphasis inline() renders, and only that: an underscore inside a word
    # (synteny_rejected_report.tsv) is part of the word, and of its anchor.
    text = re.sub(r"\*\*(.+?)\*\*", r"\1", text)
    text = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"\1", text)
    return re.sub(r"(?<!\w)_(?!\s)(.+?)(?<!\s)_(?!\w)", r"\1", text)


class Renderer:
    def __init__(self, source: str):
        self.source = source  # markdown path relative to docs/, e.g. guide/outputs.md
        self.seen: dict[str, int] = {}
        self.headings: list[tuple[int, str, str]] = []  # (level, text, slug)
        self.sections: list[dict] = []  # for the search index

    # Links

    def link(self, target: str) -> str:
        if re.match(r"^(https?:|mailto:|#)", target):
            return target
        path, _, anchor = target.partition("#")
        resolved = (DOCS / Path(self.source).parent / path).resolve()
        try:
            inside = resolved.relative_to(DOCS.resolve())
        except ValueError:
            inside = None
        if inside is not None and inside.suffix == ".md" and str(inside) in {p[1] for p in PAGES}:
            new = path[: -len(".md")] + ".html"
        elif inside is not None and inside.suffix != ".md":
            new = path
        else:
            repo_path = resolved.relative_to(Path.cwd().resolve())
            new = f"{REPO}/blob/main/{repo_path.as_posix()}"
        return new + (f"#{anchor}" if anchor else "")

    # Inline

    def inline(self, text: str) -> str:
        keep: list[str] = []

        def stash(fragment: str) -> str:
            keep.append(fragment)
            return f"\x00{len(keep) - 1}\x00"

        text = re.sub(r"(`+)(.+?)\1", lambda m: stash(f"<code>{html.escape(m.group(2).strip())}</code>"), text)
        text = re.sub(r"\\([\\`*_{}\[\]()#+\-.!|<>])", lambda m: stash(html.escape(m.group(1))), text)
        text = re.sub(r"<(https?://[^>\s]+)>", lambda m: stash(f'<a href="{html.escape(m.group(1))}">{html.escape(m.group(1))}</a>'), text)
        text = INLINE_TAGS.sub(lambda m: stash(m.group(0)), text)

        def link(m: re.Match) -> str:
            label, target = m.group(1), m.group(2).strip().split(" ")[0]
            return stash(f'<a href="{html.escape(self.link(target))}">') + label + stash("</a>")

        text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", link, text)
        text = re.sub(r"&(?!#?\w+;)", "&amp;", text).replace("<", "&lt;").replace(">", "&gt;")
        text = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", text)
        text = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"<em>\1</em>", text)
        text = re.sub(r"(?<![\w])_(?!\s)(.+?)(?<!\s)_(?![\w])", r"<em>\1</em>", text)
        while "\x00" in text:
            text = re.sub(r"\x00(\d+)\x00", lambda m: keep[int(m.group(1))], text)
        return text

    # Blocks

    def blocks(self, lines: list[str]) -> str:
        out: list[str] = []
        i = 0
        while i < len(lines):
            line = lines[i]
            stripped = line.strip()
            if not stripped:
                i += 1
                continue
            fence = re.match(r"^(\s*)(```+|~~~+)\s*(.*)$", line)
            if fence:
                indent, marker, info = len(fence.group(1)), fence.group(2), fence.group(3).strip()
                body = []
                i += 1
                while i < len(lines) and not lines[i].strip().startswith(marker):
                    body.append(lines[i][indent:] if lines[i][:indent].strip() == "" else lines[i])
                    i += 1
                i += 1
                out.append(self.code(info, body))
                continue
            heading = re.match(r"^(#{1,6})\s+(.*?)\s*#*$", line)
            if heading:
                level, text = len(heading.group(1)), heading.group(2)
                slug = slugify(plain(text), self.seen)
                self.headings.append((level, plain(text), slug))
                self.sections.append({"heading": plain(text), "slug": slug, "level": level, "text": []})
                anchor = f'<a class="anchor" href="#{slug}" aria-label="Link to this section">#</a>' if level > 1 else ""
                out.append(f'<h{level} id="{slug}">{self.inline(text)}{anchor}</h{level}>')
                i += 1
                continue
            if re.match(r"^(-{3,}|\*{3,}|_{3,})$", stripped):
                out.append("<hr>")
                i += 1
                continue
            if stripped.startswith("|") and i + 1 < len(lines) and re.match(r"^\|?\s*:?-{2,}", lines[i + 1].strip()):
                rows = []
                while i < len(lines) and lines[i].strip().startswith("|"):
                    rows.append(lines[i].strip())
                    i += 1
                out.append(self.table(rows))
                continue
            if stripped.startswith(">"):
                quoted = []
                while i < len(lines) and lines[i].strip().startswith(">"):
                    quoted.append(re.sub(r"^\s*> ?", "", lines[i]))
                    i += 1
                out.append(f"<blockquote>{self.blocks(quoted)}</blockquote>")
                continue
            if re.match(r"^<(picture|div|details|table|p|figure)\b", stripped):
                chunk = []
                while i < len(lines) and lines[i].strip():
                    chunk.append(lines[i])
                    i += 1
                out.append(self.raw_html("\n".join(chunk)))
                continue
            item = re.match(r"^(\s*)([-*+]|\d+\.)\s+", line)
            if item:
                i, rendered = self.list_block(lines, i)
                out.append(rendered)
                continue
            para = []
            while i < len(lines) and lines[i].strip() and not self.starts_block(lines, i):
                para.append(lines[i].strip())
                i += 1
            text = " ".join(para)
            self.note(text)
            out.append(f"<p>{self.inline(text)}</p>")
        return "\n".join(out)

    def starts_block(self, lines: list[str], i: int) -> bool:
        line = lines[i]
        s = line.strip()
        return bool(
            re.match(r"^(\s*)(```|~~~)", line)
            or re.match(r"^#{1,6}\s", line)
            or re.match(r"^(\s*)([-*+]|\d+\.)\s+", line)
            or s.startswith(">")
            or (s.startswith("|") and i + 1 < len(lines) and re.match(r"^\|?\s*:?-{2,}", lines[i + 1].strip()))
        )

    def note(self, text: str) -> None:
        if self.sections:
            self.sections[-1]["text"].append(plain(text))

    def list_block(self, lines: list[str], i: int) -> tuple[int, str]:
        first = re.match(r"^(\s*)([-*+]|\d+\.)\s+", lines[i])
        base = len(first.group(1))
        ordered = first.group(2)[0].isdigit()
        start = int(first.group(2)[:-1]) if ordered else 1
        items: list[list[str]] = []
        loose = False
        while i < len(lines):
            line = lines[i]
            m = re.match(r"^(\s*)([-*+]|\d+\.)(\s+)", line)
            if m and len(m.group(1)) == base and (m.group(2)[0].isdigit()) == ordered:
                content_indent = len(m.group(1)) + len(m.group(2)) + len(m.group(3))
                items.append([line[content_indent:]])
                i += 1
                while i < len(lines):
                    nxt = lines[i]
                    if not nxt.strip():
                        # A blank line ends the list unless the next line carries on in it.
                        j = i + 1
                        while j < len(lines) and not lines[j].strip():
                            j += 1
                        if j < len(lines) and (len(lines[j]) - len(lines[j].lstrip())) > base:
                            items[-1].append("")
                            loose = loose or not re.match(r"^\s*([-*+]|\d+\.)\s", lines[j])
                            i += 1
                            continue
                        if j < len(lines) and re.match(rf"^\s{{{base}}}([-*+]|\d+\.)\s", lines[j]):
                            loose = True
                            i = j
                        break
                    indent = len(nxt) - len(nxt.lstrip())
                    if indent > base:
                        items[-1].append(nxt[min(indent, content_indent):])
                        i += 1
                        continue
                    if re.match(r"^(\s*)([-*+]|\d+\.)\s+", nxt):
                        break
                    items[-1].append(nxt.strip())  # lazy continuation
                    i += 1
                continue
            break
        tag = "ol" if ordered else "ul"
        attrs = f' start="{start}"' if ordered and start != 1 else ""
        rendered = []
        for content in items:
            body = self.blocks(content)
            if not loose and body.startswith("<p>") and body.count("<p>") == 1:
                body = body[3:].replace("</p>", "", 1)
            rendered.append(f"<li>{body}</li>")
        return i, f"<{tag}{attrs}>" + "".join(rendered) + f"</{tag}>"

    def table(self, rows: list[str]) -> str:
        def cells(row: str) -> list[str]:
            row = row.strip().strip("|")
            return [c.strip().replace("\\|", "|") for c in re.split(r"(?<!\\)\|", row)]

        head = cells(rows[0])
        body = [cells(r) for r in rows[2:]]
        for row in body:
            self.note(" ".join(row))
        out = ['<div class="table-wrap"><table><thead><tr>']
        out += [f"<th>{self.inline(c)}</th>" for c in head]
        out.append("</tr></thead><tbody>")
        for row in body:
            out.append("<tr>" + "".join(f"<td>{self.inline(c)}</td>" for c in row) + "</tr>")
        out.append("</tbody></table></div>")
        return "".join(out)

    def code(self, info: str, body: list[str]) -> str:
        if info == "mermaid":
            return self.flow(body)
        title = re.match(r'^title="([^"]*)"', info)
        caption = f'<div class="code-title">{html.escape(title.group(1))}</div>' if title else ""
        lang = "" if title or not info else f' data-lang="{html.escape(info.split()[0])}"'
        text = html.escape("\n".join(body))
        return f'<div class="code"{lang}>{caption}<pre>{text}</pre></div>'

    def flow(self, body: list[str]) -> str:
        """A mermaid flowchart as the list of its steps, each with where it leads."""
        labels: dict[str, str] = {}
        decisions: set[str] = set()
        order: list[str] = []
        edges: list[tuple[str, str, str]] = []
        node = r'(\w+)(?:\["([^"]*)"\]|\{"([^"]*)"\})?'
        for line in body[1:]:
            m = re.match(rf"^\s*{node}\s*-->(?:\|([^|]*)\|)?\s*{node}\s*$", line)
            if not m:
                continue
            a, a_box, a_dec, label, b, b_box, b_dec = m.groups()
            for name, box, dec in ((a, a_box, a_dec), (b, b_box, b_dec)):
                if name not in order:
                    order.append(name)
                if box or dec:
                    labels[name] = box or dec
                if dec:
                    decisions.add(name)
            edges.append((a, b, label or ""))
        incoming = {n: 0 for n in order}
        for _, b, _ in edges:
            incoming[b] += 1
        ready = [n for n in order if incoming[n] == 0]
        sequence = []
        while ready:
            n = ready.pop(0)
            sequence.append(n)
            for a, b, _ in edges:
                if a == n:
                    incoming[b] -= 1
                    if incoming[b] == 0:
                        ready.append(b)
        sequence += [n for n in order if n not in sequence]

        def parts(name: str) -> list[str]:
            return re.split(r"<br\s*/?>", labels.get(name, name))

        def label_html(name: str) -> str:
            head, *rest = parts(name)
            return f"<b>{html.escape(head)}</b>" + "".join(f"<span>{html.escape(p)}</span>" for p in rest)

        def next_html(target: str, label: str) -> str:
            note = f" <i>{html.escape(label.strip())}</i>" if label.strip() else ""
            return f"<li>&rarr; {html.escape(parts(target)[0])}{note}</li>"

        out = ['<ol class="flow">']
        for n in sequence:
            kind = " flow-decision" if n in decisions else ""
            nexts = "".join(next_html(b, lab) for a, b, lab in edges if a == n)
            out.append(
                f'<li class="flow-step{kind}"><div class="flow-box">{label_html(n)}</div>'
                + (f'<ul class="flow-next">{nexts}</ul>' if nexts else "")
                + "</li>"
            )
        out.append("</ol>")
        self.note(" ".join(plain(labels.get(n, n)).replace("<br/>", " ") for n in sequence))
        return "".join(out)

    def raw_html(self, chunk: str) -> str:
        """Raw HTML passes through. A <picture> that swaps on the system colour
        scheme follows the site's own theme instead."""
        picture = re.search(
            r'<source media="\(prefers-color-scheme: dark\)" srcset="([^"]+)">\s*<img src="([^"]+)" alt="([^"]*)">', chunk
        )
        if picture:
            dark, light, alt = picture.groups()
            return (
                f'<p class="figure"><img class="only-dark" src="{dark}" alt="{alt}">'
                f'<img class="only-light" src="{light}" alt="{alt}"></p>'
            )
        return chunk


# ── Pages ──


def header(root: str, current: str) -> str:
    def nav(href: str, label: str, key: str, optional: bool = True) -> str:
        cls = ' class="nav-optional"' if optional else ""
        here = ' aria-current="page"' if key == current else ""
        return f'      <a{cls} href="{href}"{here}>{label}</a>\n'

    return (
        '<header class="top">\n  <div class="wrap">\n'
        f'    <a class="brand" href="{root}index.html" aria-label="Recore home">\n'
        f'      <img class="logo-light" src="{root}assets/recore_logo.svg" alt="Recore" width="74" height="34">\n'
        f'      <img class="logo-dark" src="{root}assets/recore_logo_dark.svg" alt="Recore" width="74" height="34">\n'
        "    </a>\n"
        '    <nav class="nav" aria-label="Site">\n'
        + nav(f"{root}index.html#start", "Install", "install")
        + nav(f"{root}index.html#command", "Your command", "command")
        + nav(f"{root}index.html#try", "Try it", "try")
        + nav(f"{root}{INDEX}", "Docs", "docs", optional=False)
        + nav(f"{root}reference.html", "Reference", "reference")
        + f'      <a class="gh-link" href="{REPO}">{GH_SVG}<span>GitHub</span></a>\n'
        + f'      <button class="icon-btn" type="button" data-theme-toggle aria-label="Switch between light and dark">\n        {SUN_SVG}\n      </button>\n'
        "    </nav>\n  </div>\n</header>\n"
    )


def footer(root: str) -> str:
    return (
        '<footer class="foot">\n  <div class="wrap">\n    <div>\n      <h3>Cite</h3>\n'
        "      <blockquote>Ruiz-Rodriguez P, et al. (2026). Recore: post-hoc rescue and quality control of pangenome "
        "presence/absence matrices. <i>In preparation.</i></blockquote>\n    </div>\n"
        "    <div>\n      <h3>Learn more</h3>\n      <ul>\n"
        f'        <li><a href="{root}index.html">Get started</a></li>\n'
        f'        <li><a href="{root}{INDEX}">Documentation</a></li>\n'
        f'        <li><a href="{root}reference.html">Reference: options, files, states</a></li>\n'
        f'        <li><a href="{REPO}/blob/main/CHANGELOG.md">Changelog</a></li>\n'
        "      </ul>\n    </div>\n    <div>\n      <h3>Project</h3>\n      <ul>\n"
        f'        <li><a href="{REPO}">GitHub</a></li>\n'
        f'        <li><a href="{REPO}/discussions">Questions</a></li>\n'
        f'        <li><a href="{REPO}/issues">Report a bug</a></li>\n'
        "        <li>GPL-3.0 &middot; PathoGenOmics Lab</li>\n      </ul>\n    </div>\n  </div>\n</footer>\n"
    )


def sidebar(root: str, current: str | None) -> str:
    out = [
        '<aside class="docs-nav" aria-label="Documentation">\n',
        '  <div class="docs-search">\n'
        '    <input type="search" placeholder="Search the docs" aria-label="Search the documentation" data-docs-search>\n'
        '    <div class="docs-results" data-docs-results hidden></div>\n  </div>\n',
        '  <details class="docs-menu" open>\n    <summary>Pages</summary>\n',
    ]
    section = None
    for sec, path, title, _ in PAGES:
        if sec != section:
            if section is not None:
                out.append("    </ul>\n")
            out.append(f'    <p class="docs-section">{sec}</p>\n    <ul>\n')
            section = sec
        here = ' aria-current="page"' if path == current else ""
        out.append(f'      <li><a href="{root}{path[:-3]}.html"{here}>{title}</a></li>\n')
    out.append("    </ul>\n  </details>\n</aside>\n")
    return "".join(out)


def page(root: str, title: str, description: str, body: str) -> str:
    return (
        "<!doctype html>\n"
        '<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"<title>{html.escape(title)}</title>\n"
        f'<meta name="description" content="{html.escape(description)}">\n'
        f'<link rel="icon" href="{root}assets/recore_logo.svg" type="image/svg+xml">\n'
        f'<link rel="stylesheet" href="{root}assets/site.css">\n'
        f"{THEME_BOOT}\n</head>\n<body>\n"
        '<a class="skip" href="#content">Skip to the content</a>\n'
        "<!-- Generated by tools/build_docs.py from the markdown next to this file. Edit that, then run it. -->\n"
        + body
        + f'<script src="{root}assets/search-index.js"></script>\n'
        + f'<script src="{root}assets/site.js"></script>\n'
        + "</body>\n</html>\n"
    )


def doc_page(position: int) -> tuple[Path, str, list[dict]]:
    section, path, title, description = PAGES[position]
    renderer = Renderer(path)
    content = renderer.blocks((DOCS / path).read_text(encoding="utf-8").splitlines())
    root = "../" * path.count("/")
    toc = "".join(
        f'<li class="toc-{level}"><a href="#{slug}">{html.escape(text)}</a></li>'
        for level, text, slug in renderer.headings
        if level in (2, 3)
    )
    prev_page = PAGES[position - 1] if position > 0 else None
    next_page = PAGES[position + 1] if position + 1 < len(PAGES) else None
    pager = '<nav class="pager" aria-label="Previous and next page">'
    if prev_page:
        pager += f'<a class="prev" href="{root}{prev_page[1][:-3]}.html"><span>Previous</span>{prev_page[2]}</a>'
    if next_page:
        pager += f'<a class="next" href="{root}{next_page[1][:-3]}.html"><span>Next</span>{next_page[2]}</a>'
    pager += "</nav>"
    body = (
        header(root, "docs")
        + '<div class="wrap docs-layout">\n'
        + sidebar(root, path)
        + f'<main class="prose" id="content">\n<p class="crumb"><a href="{root}{INDEX}">Docs</a> / {section}</p>\n'
        + content
        + f'\n<p class="edit"><a href="{REPO}/blob/main/docs/{path}">{GH_SVG}Edit this page on GitHub</a></p>\n'
        + pager
        + "\n</main>\n"
        + (f'<aside class="docs-toc" aria-label="On this page"><p>On this page</p><ul>{toc}</ul></aside>\n' if toc else "")
        + "</div>\n"
        + footer(root)
    )
    entries = [
        {
            "page": title,
            "section": s["heading"] if s["level"] > 1 else "",
            "url": f"{path[:-3]}.html" + (f"#{s['slug']}" if s["level"] > 1 else ""),
            "text": " ".join(s["text"])[:600],
        }
        for s in renderer.sections
    ]
    return DOCS / f"{path[:-3]}.html", page(root, f"{title} | Recore docs", description, body), entries


def index_page() -> tuple[Path, str]:
    cards = []
    section = None
    for sec, path, title, line in PAGES:
        if sec != section:
            if section is not None:
                cards.append("</div>\n")
            cards.append(f'<h2 class="docs-group">{sec}</h2>\n<div class="deep">\n')
            section = sec
        cards.append(f'<a href="{path[:-3]}.html"><b>{title}</b><span>{line}</span></a>\n')
    cards.append("</div>\n")
    body = (
        header("", "docs")
        + '<div class="wrap docs-layout docs-home">\n'
        + sidebar("", None)
        + '<main class="prose" id="content">\n<h1>Documentation</h1>\n'
        + "<p>The long versions: every option, every file and column, and how each decision is made. "
        + 'New to Recore? <a href="index.html">Start on the front page</a>, it takes two minutes.</p>\n'
        + "".join(cards)
        + "</main>\n</div>\n"
        + footer("")
    )
    return DOCS / INDEX, page("", "Documentation | Recore", "Recore's documentation, page by page.", body)


def build() -> dict[Path, str]:
    files: dict[Path, str] = {}
    entries: list[dict] = []
    for position in range(len(PAGES)):
        path, text, page_entries = doc_page(position)
        files[path] = text
        entries += page_entries
    index_path, index_text = index_page()
    files[index_path] = index_text
    files[SEARCH_INDEX] = (
        "/* Generated by tools/build_docs.py: the documentation, section by section, for the search box. */\n"
        "window.RECORE_DOCS_INDEX = " + json.dumps(entries, ensure_ascii=False, indent=0) + ";\n"
    )
    return files


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true", help="exit 1 if a generated page differs from what is committed")
    args = parser.parse_args()
    files = build()
    stale = [p for p, text in files.items() if not p.exists() or p.read_text(encoding="utf-8") != text]
    if args.check:
        if stale:
            print("These pages no longer match their markdown; run python3 tools/build_docs.py:", file=sys.stderr)
            for p in stale:
                print(f"  {p}", file=sys.stderr)
            return 1
        print(f"docs/: {len(files)} generated files match their markdown")
        return 0
    for p in stale:
        p.write_text(files[p], encoding="utf-8")
    print(f"docs/: wrote {len(stale)} of {len(files)} generated files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
