"""Tests for delfin.orca_knowledge.manual_text: HTML -> text conversion.

Fixtures are small standalone Sphinx/furo-style pages (no ORCA manual text in
the repository). They check the conversions the old plain-text conversion
lost: full table rows, fenced code blocks, nested lists, and headings that
keep their section number.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from delfin.orca_knowledge.manual_text import convert_html, convert_page

PAGE = """\
<!doctype html><html lang="en"><head><title>2.6. SCF - ORCA 6.1 Manual</title></head>
<body>
<nav><a href="index.html">Index</a><ul><li>nav item</li></ul></nav>
<article role="main">
<h1><span class="section-number">2.6. </span>Self-Consistent-Field (SCF)</h1>
<p>Intro paragraph about <code>!StrongSCF</code> convergence.</p>
<h2><span class="section-number">2.6.1. </span>Convergence Tolerances</h2>
<div class="table-wrapper docutils container">
<table>
<thead><tr><th><p>Variable</p></th><th><p>Value</p></th></tr></thead>
<tbody>
<tr><td colspan="2"><p><em>Sloppy</em> (<code>!SloppySCF</code>)</p></td></tr>
<tr><td><p><code>TolE</code></p></td>
    <td><p>3e-5, the energy change between two successive iterations
        is smaller than this threshold</p></td></tr>
</tbody>
</table>
</div>
<div class="highlight-orca notranslate"><div class="highlight"><pre><span></span>%scf
  maxiter 200
end
</pre></div></div>
<ul>
<li>first item
<ul>
<li>nested item A</li>
<li>nested item B</li>
</ul>
</li>
<li>second item</li>
</ul>
</article>
</body></html>
"""


def _text() -> str:
    return convert_html(PAGE, version="6.1", relpath="contents/essentialelements/scf.html")


def test_source_line_first():
    text = _text()
    assert text.splitlines()[0] == "SOURCE: ORCA 6.1 manual, page contents/essentialelements/scf.html"


def test_navigation_excluded():
    text = _text()
    assert "nav item" not in text
    assert text.count("Index") == 0


def test_headings_keep_section_number():
    text = _text()
    assert "2.6. Self-Consistent-Field (SCF)" in text
    assert "2.6.1. Convergence Tolerances" in text


def test_table_full_cell_text():
    text = _text()
    assert "Variable | Value" in text
    row = [ln for ln in text.splitlines() if "TolE" in ln and "|" in ln]
    assert len(row) == 1
    # every cell's full text, including a long wrapped value cell, on one line
    assert row[0] == "TolE | 3e-5, the energy change between two successive iterations is smaller than this threshold"
    # spanned rows keep their single (merged) cell
    assert "Sloppy (!SloppySCF)" in text


def test_code_block_fenced():
    text = _text()
    assert "```\n%scf\n  maxiter 200\nend\n```" in text


def test_nested_lists():
    text = _text()
    lines = text.splitlines()
    assert "- first item" in lines
    assert "  - nested item A" in lines
    assert "  - nested item B" in lines
    assert "- second item" in lines


def test_paragraph_inline_code_kept():
    assert "!StrongSCF" in _text()


def test_convert_page_reports_relative_path(tmp_path: Path):
    src = tmp_path / "page.html"
    src.write_text(PAGE, encoding="utf-8")
    dst = tmp_path / "out"
    rel = convert_page(src, dst, version="6.1")
    assert rel == "page.html"
    text = (dst / "page.txt").read_text(encoding="utf-8")
    assert text.splitlines()[0] == "SOURCE: ORCA 6.1 manual, page page.html"


def test_cli_converts_every_page(tmp_path: Path):
    from conftest import child_env
    src = tmp_path / "manual"
    (src / "contents" / "sub").mkdir(parents=True)
    (src / "contents" / "a.html").write_text(PAGE, encoding="utf-8")
    (src / "contents" / "sub" / "b.html").write_text(PAGE, encoding="utf-8")
    (src / "style.css").write_text("body {}", encoding="utf-8")
    dst = tmp_path / "text"
    result = subprocess.run(
        [sys.executable, "-m", "delfin.orca_knowledge.manual_text",
         str(src), str(dst), "--version", "6.1"],
        capture_output=True, text=True, env=child_env(tmp_path),
    )
    assert result.returncode == 0, result.stderr
    a = (dst / "contents" / "a.txt").read_text(encoding="utf-8")
    assert a.splitlines()[0] == "SOURCE: ORCA 6.1 manual, page contents/a.html"
    b = (dst / "contents" / "sub" / "b.txt").read_text(encoding="utf-8")
    assert b.splitlines()[0] == "SOURCE: ORCA 6.1 manual, page contents/sub/b.html"
    assert not (dst / "style.txt").exists()
