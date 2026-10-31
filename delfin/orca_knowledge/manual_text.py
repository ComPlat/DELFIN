"""Convert one ORCA-manual HTML page (Sphinx/furo) to plain text.

The previous conversion of the manual to plain text lost table content
(rows truncated) and flattened code blocks.  This module re-converts each
page's main article only, so agents of this swarm read a faithful text:

* the page's ``<article role="main">`` only -- no navigation, no footer;
* headings keep their section number (``2.6. Self-Consistent-Field (SCF)``);
* code/pre blocks are fenced with `` ``` `` and keep their line breaks;
* tables are rendered row by row as ``cell | cell | cell`` with every
  cell's full text on one line (a cell spanning several columns is a
  single entry on its row);
* list items render as ``- item``, nested lists indented under the item;
* the first line is ``SOURCE: ORCA <version> manual, page <relative path>``.

CLI::

    python -m delfin.orca_knowledge.manual_text SRC_DIR DST_DIR --version 6.1

converts every ``.html`` page under ``SRC_DIR``, preserving the relative
directory layout, into ``.txt`` files under ``DST_DIR``.

Only the standard library and BeautifulSoup are used.  No manual text
itself lives in this repository: the manual HTML is read from an external
directory the operator points the CLI at.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterator

from bs4 import BeautifulSoup, Tag

__all__ = ["convert_html", "convert_page", "convert_tree", "main"]

_HEADINGS = ("h1", "h2", "h3", "h4", "h5", "h6")


def _clean(text: str) -> str:
    """Collapse whitespace runs to single spaces and trim the ends."""
    return " ".join(text.split())


def _inline(tag: Tag) -> str:
    """Flatten a tag's content to one line of text.

    Raw concatenation (no separator): the markup carries its own spaces,
    and a separator would inject spaces around punctuation -- get_text(' ')
    turns ``(!SloppySCF)`` into ``( !SloppySCF )``.
    """
    return _clean(tag.get_text(""))


def _is_heading(tag: Tag) -> bool:
    return tag.name in _HEADINGS


def _skip(tag: Tag) -> bool:
    """Content that never belongs in the text (permalinks, edit links)."""
    return "headerlink" in (tag.get("class") or [])


def _find_article(soup: BeautifulSoup) -> Tag:
    article = soup.find("article")
    if isinstance(article, Tag):
        return article
    raise ValueError(
        "no <article> element found -- is this a Sphinx/furo manual page?"
    )


def _table_lines(tag: Tag) -> Iterator[str]:
    """Render a table row by row as 'cell | cell | cell'.

    Every cell's full text is kept on one line -- this is the part the old
    conversion truncated.  A cell spanning several columns is a single
    entry on its row.
    """
    for tr in tag.find_all("tr"):
        cells = [_inline(td) for td in tr.find_all(("th", "td"), recursive=False)]
        yield " | ".join(cell for cell in cells if cell)


def _code_lines(tag: Tag) -> Iterator[str]:
    """Render a code/pre block as a fenced block, keeping line breaks.

    Sphinx's Pygments highlighting splits tokens into spans; the block's
    own text is the code, so ``get_text()`` on the ``<pre>`` is enough.
    """
    pre = tag if tag.name == "pre" else tag.find("pre")
    if pre is None:
        return
    text = pre.get_text().strip("\n")
    if _clean(text):
        yield "```"
        yield from text.split("\n")
        yield "```"


def _list_lines(tag: Tag, indent: str) -> Iterator[str]:
    """Render a ul/ol as '- item' lines, nested lists indented further."""
    for li in tag.find_all("li", recursive=False):
        # Nested lists live either directly in the li or inside one of its
        # direct block children (a <p> wrapper); lift them out so the item's
        # own text is not doubled and they indent under the item.
        nested: list[Tag] = []
        for holder in (li, *[c for c in li.children if isinstance(c, Tag)]):
            for sub in holder.find_all(("ul", "ol"), recursive=False):
                nested.append(sub)
                sub.extract()
        text = _inline(li)
        yield f"{indent}- {text}" if text else f"{indent}-"
        for sub in nested:
            yield from _list_lines(sub, indent + "  ")


def _walk(tag: Tag, out: list[str]) -> None:
    """Append the text lines of a block-level walk of ``tag`` to ``out``.

    Document order, one entry per block element; container elements
    (div, admonitions, figures ...) are recursed into so nothing inside
    them is lost, while leaves render themselves.
    """
    for child in tag.children:
        if not isinstance(child, Tag) or child.name in ("style", "script"):
            continue
        if _skip(child):
            continue
        # Permalink pilcrows on headings and table/figure captions: their
        # <a class="headerlink"> text is the pilcrow, not real content.
        for link in child.find_all("a", class_="headerlink"):
            link.extract()
        if _is_heading(child):
            out.append(_inline(child))
        elif child.name == "table":
            out.extend(line for line in _table_lines(child) if line)
        elif child.name == "pre" or (
            child.name == "div" and "highlight" in (child.get("class") or [])
        ):
            out.extend(_code_lines(child))
        elif child.name in ("ul", "ol"):
            out.extend(_list_lines(child, ""))
        elif child.name == "p":
            text = _inline(child)
            if text:
                out.append(text)
        else:
            # Containers and anything else: recurse if it holds blocks,
            # otherwise render it as one inline line.
            if child.find(
                ("table", "pre", "ul", "ol", *_HEADINGS), recursive=False
            ) or child.name in ("div", "blockquote", "section", "figure"):
                _walk(child, out)
            else:
                text = _inline(child)
                if text:
                    out.append(text)


def convert_html(html: str, *, version: str, relpath: str) -> str:
    """Convert one manual page's HTML to text. See module docstring."""
    soup = BeautifulSoup(html, "html.parser")
    article = _find_article(soup)
    body: list[str] = []
    _walk(article, body)
    lines = [f"SOURCE: ORCA {version} manual, page {relpath}", ""]
    lines.extend(body)
    text = "\n".join(lines)
    # collapse accidental blank runs, trim trailing whitespace
    collapsed: list[str] = []
    for line in text.split("\n"):
        if line == "" and (not collapsed or collapsed[-1] == ""):
            continue
        collapsed.append(line)
    return "\n".join(collapsed).rstrip() + "\n"


def convert_page(
    page: Path, dst_dir: Path, *, version: str, root: Path | None = None
) -> str:
    """Convert one HTML file, write ``<relpath>.txt`` under ``dst_dir``,
    return the relative page path used in the SOURCE line.  ``root`` is
    the manual root the page's path is relative to (default: its parent)."""
    root = page.parent.resolve() if root is None else root.resolve()
    rel = page.resolve().relative_to(root).as_posix()
    text = convert_html(
        page.read_text(encoding="utf-8", errors="replace"),
        version=version,
        relpath=rel,
    )
    out = (dst_dir / rel).with_suffix(".txt")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    return rel


def convert_tree(src_dir: Path, dst_dir: Path, *, version: str) -> int:
    """Convert every ``.html`` page under ``src_dir``. Returns the count."""
    src_dir = src_dir.resolve()
    count = 0
    for page in sorted(src_dir.rglob("*.html")):
        convert_page(page, dst_dir, version=version, root=src_dir)
        count += 1
    return count


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m delfin.orca_knowledge.manual_text",
        description="Convert ORCA manual HTML pages to plain text.",
    )
    parser.add_argument("src_dir", type=Path, help="manual HTML root")
    parser.add_argument("dst_dir", type=Path, help="output root for .txt files")
    parser.add_argument(
        "--version", required=True,
        help="ORCA version for the SOURCE line, e.g. 6.1",
    )
    args = parser.parse_args(argv)
    if not args.src_dir.is_dir():
        parser.error(f"not a directory: {args.src_dir}")
    count = convert_tree(args.src_dir, args.dst_dir, version=args.version)
    print(f"converted {count} pages -> {args.dst_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
