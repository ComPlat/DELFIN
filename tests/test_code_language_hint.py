"""The write path must flag non-English text inside freshly written code.

The shared work cycle requires English INSIDE code (comments, docstrings,
identifiers, log/error strings) while the agent replies to the user in the
user's own language. Two consecutive field runs produced modules with German
docstrings and comments although the prompt already carried the rule, so the
rule now lives in the write path itself: a successful write/edit result gains
a short advisory note. The write is never blocked and the file is never
touched by the check.

Covered here:
- the two verbatim field strings are flagged;
- correct English code, MIT/URL headers and Greek science symbols are not;
- prose files (.md/.rst/.txt/.json) are exempt and written byte-identical;
- edits only report text the agent just inserted, not legacy lines;
- notebook code cells are checked, markdown cells are not;
- caps, binary-ish and malformed content never raise.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent.api_client import (
    KitToolPermissions,
    _doc_executor,
    _lang_has_marker,
    _lang_scan_text,
    _language_hint_for_write,
)

NOTE_MARKER = "non-English text in code"

# The two strings that survived into generated code in the field.
FIELD_DOCSTRING = '"""Tetris-Spiel mit Voilà-Integration."""'
FIELD_COMMENT = "# Project root zur Path hinzufügen"


@pytest.fixture
def workspace(tmp_path) -> Path:
    ws = tmp_path / "ws"
    ws.mkdir()
    return ws


@pytest.fixture
def perms(workspace) -> KitToolPermissions:
    return KitToolPermissions(workspace=workspace, mode="bypassPermissions")


def _write(perms, name: str, content: str) -> str:
    return _doc_executor.execute(
        "write_file", {"path": name, "content": content}, perms)


def _read_baseline(perms, name: str) -> None:
    _doc_executor.execute("read_file", {"path": name}, perms)


def _note(out: str) -> str:
    """Just the advisory note, so assertions never match the echoed diff."""
    idx = out.find("\n\nNote: " + NOTE_MARKER)
    return out[idx:] if idx >= 0 else ""


# ---------------------------------------------------------------------------
# The field failures
# ---------------------------------------------------------------------------


def test_field_case_german_docstring_is_flagged(workspace, perms):
    out = _write(perms, "tetris.py", FIELD_DOCSTRING + "\n\nimport sys\n")
    assert NOTE_MARKER in out
    assert "Tetris-Spiel" in _note(out)
    assert "line 1" in _note(out)


def test_field_case_german_comment_is_flagged(workspace, perms):
    content = (
        '"""Entry point."""\n'
        "import sys\n"
        "\n"
        + FIELD_COMMENT + "\n"
        "sys.path.insert(0, '.')\n"
    )
    out = _write(perms, "run.py", content)
    assert NOTE_MARKER in out
    assert "line 4" in _note(out)
    assert "hinzuf" in _note(out)


def test_note_states_the_rule(workspace, perms):
    out = _write(perms, "g.py", FIELD_COMMENT + "\n")
    assert "Code stays English" in _note(out)
    assert "user's language" in _note(out)


# ---------------------------------------------------------------------------
# What must NOT be flagged
# ---------------------------------------------------------------------------


ENGLISH_MODULE = '''"""Utilities for the workflow runner.

Reads a config file, validates it and returns the parsed steps. The parser is
strict: an unknown key is an error, and a missing file raises immediately.
"""

import logging

logger = logging.getLogger(__name__)

# The default timeout is intentionally short so a stuck step is noticed fast.
DEFAULT_TIMEOUT = 30


def run(steps, timeout=DEFAULT_TIMEOUT):
    """Run every step in order and return the collected results."""
    results = []
    for step in steps:  # keep the original ordering
        logger.info("running step %s", step)
        if not step:
            raise ValueError("empty step in the pipeline")
        results.append(step)
    print("done, {} steps".format(len(results)))
    return results
'''


def test_correct_english_module_is_not_flagged(workspace, perms):
    out = _write(perms, "runner.py", ENGLISH_MODULE)
    assert NOTE_MARKER not in out
    assert "File created" in out


def test_license_headers_urls_and_science_symbols_are_not_flagged(
        workspace, perms):
    content = (
        "# MIT License - see LICENSE for the full text\n"
        "# Docs: https://example.org/guide#install\n"
        "# Author: Jürgen Weber\n"
        "# Total DAS capacity, ΔE in kcal/mol, α orbitals only.\n"
        "TAG = '#not-a-comment'\n"
    )
    out = _write(perms, "header.py", content)
    assert NOTE_MARKER not in out


def test_german_prose_in_markdown_is_not_flagged(workspace, perms):
    content = (
        "# Anleitung\n\n"
        "Dieses Dokument beschreibt, wie das Programm gestartet wird.\n"
    )
    out = _write(perms, "README.md", content)
    assert NOTE_MARKER not in out
    assert (workspace / "README.md").read_text(encoding="utf-8") == content


@pytest.mark.parametrize("name", ["notes.txt", "guide.rst", "data.json"])
def test_prose_and_data_files_are_exempt(workspace, perms, name):
    content = '{"titel": "Fehler beim Lesen der Datei"}\n'
    out = _write(perms, name, content)
    assert NOTE_MARKER not in out
    assert (workspace / name).read_text(encoding="utf-8") == content


def test_german_data_literal_is_not_flagged(workspace, perms):
    # Translation tables and template blobs are data the user asked for.
    content = (
        '"""Localised labels."""\n'
        "LABELS = {'de': 'Fehler beim Lesen der Datei'}\n"
        'TEMPLATE = """Sehr geehrte Damen und Herren,\n'
        'die Auswertung ist beigefuegt.\n'
        '"""\n'
    )
    out = _write(perms, "labels.py", content)
    assert NOTE_MARKER not in out


def test_german_in_a_print_call_is_flagged(workspace, perms):
    content = (
        '"""Entry point."""\n'
        'print("Berechnung wurde abgebrochen")\n'
    )
    out = _write(perms, "cli.py", content)
    assert NOTE_MARKER in out
    assert "line 2" in _note(out)


def test_german_in_a_raise_message_is_flagged(workspace, perms):
    content = (
        '"""Loader."""\n'
        'def load(p):\n'
        '    raise ValueError("Datei wurde nicht gefunden")\n'
    )
    out = _write(perms, "loader.py", content)
    assert NOTE_MARKER in out
    assert "line 3" in _note(out)


# ---------------------------------------------------------------------------
# The write itself stays untouched
# ---------------------------------------------------------------------------


def test_flagged_write_leaves_the_file_byte_identical(workspace, perms):
    content = FIELD_DOCSTRING + "\n" + FIELD_COMMENT + "\nx = 1\n"
    out = _write(perms, "app.py", content)
    assert NOTE_MARKER in out
    assert (workspace / "app.py").read_bytes() == content.encode("utf-8")


def test_flagged_write_still_reports_success(workspace, perms):
    out = _write(perms, "app2.py", FIELD_COMMENT + "\n")
    assert out.startswith("File created")
    assert "app2.py" in out


def test_note_is_short(workspace, perms):
    out = _write(perms, "app3.py", FIELD_COMMENT + "\n")
    assert 0 < len(_note(out)) < 700


# ---------------------------------------------------------------------------
# Edits: only newly inserted text is reported
# ---------------------------------------------------------------------------


def test_edit_file_flags_newly_inserted_german(workspace, perms):
    _write(perms, "e1.py", '"""Module."""\nx = 1\n')
    _read_baseline(perms, "e1.py")
    out = _doc_executor.execute("edit_file", {
        "path": "e1.py", "old_string": "x = 1",
        "new_string": "# Wert wird neu berechnet\nx = 2",
    }, perms)
    assert NOTE_MARKER in out
    assert "Wert wird" in _note(out)


def test_edit_file_stays_silent_for_english_insert_into_german_file(
        workspace, perms):
    _write(perms, "e2.py", '"""Module."""\n' + FIELD_COMMENT + "\nx = 1\n")
    _read_baseline(perms, "e2.py")
    out = _doc_executor.execute("edit_file", {
        "path": "e2.py", "old_string": "x = 1",
        "new_string": "# recompute the cached value\nx = 2",
    }, perms)
    assert NOTE_MARKER not in out


def test_multi_edit_flags_new_german_comment(workspace, perms):
    _write(perms, "m1.py", '"""Module."""\na = 1\nb = 2\n')
    _read_baseline(perms, "m1.py")
    out = _doc_executor.execute("multi_edit", {
        "path": "m1.py",
        "edits": [
            {"old_string": "a = 1", "new_string": "# keep the first value\na = 1"},
            {"old_string": "b = 2", "new_string": "# Zeile wird spaeter entfernt\nb = 2"},
        ],
    }, perms)
    assert NOTE_MARKER in out
    assert "Zeile wird" in _note(out)
    assert "keep the first value" not in _note(out)


def test_multi_edit_stays_silent_for_english_edits(workspace, perms):
    _write(perms, "m2.py", '"""Module."""\na = 1\nb = 2\n')
    _read_baseline(perms, "m2.py")
    out = _doc_executor.execute("multi_edit", {
        "path": "m2.py",
        "edits": [
            {"old_string": "a = 1", "new_string": "# first value\na = 1"},
            {"old_string": "b = 2", "new_string": "# second value\nb = 2"},
        ],
    }, perms)
    assert NOTE_MARKER not in out


def test_fuzzy_edit_path_reports_reindented_insert(workspace, perms):
    _write(perms, "f1.py", '"""Module."""\ndef go():\n    return 1\n')
    _read_baseline(perms, "f1.py")
    out = _doc_executor.execute("edit_file", {
        "path": "f1.py",
        # Indentation deliberately wrong -> fuzzy fallback.
        "old_string": "return 1",
        "new_string": "# Rueckgabe wird nicht mehr benoetigt\n    return 2",
    }, perms)
    assert "Edited" in out
    assert NOTE_MARKER in out


# ---------------------------------------------------------------------------
# Notebooks
# ---------------------------------------------------------------------------


def _nb(cells: list[dict]) -> str:
    return json.dumps(
        {"cells": cells, "metadata": {}, "nbformat": 4, "nbformat_minor": 5})


def _code_cell(src: str) -> dict:
    return {"cell_type": "code", "source": src, "outputs": [],
            "execution_count": None, "metadata": {}}


def test_write_notebook_flags_code_cell_only(workspace, perms):
    content = _nb([
        {"cell_type": "markdown", "source": "# Anleitung\nDas ist die Doku.\n",
         "metadata": {}},
        _code_cell(FIELD_COMMENT + "\nimport sys\n"),
    ])
    out = _write(perms, "nb.ipynb", content)
    assert NOTE_MARKER in out
    assert "cell 1" in _note(out)
    assert "Anleitung" not in _note(out)


def test_write_notebook_with_english_code_is_not_flagged(workspace, perms):
    content = _nb([
        {"cell_type": "markdown", "source": "# Anleitung\nDas ist die Doku.\n",
         "metadata": {}},
        _code_cell("# load the dataset\nimport sys\n"),
    ])
    out = _write(perms, "nb2.ipynb", content)
    assert NOTE_MARKER not in out


def test_notebook_edit_notes_german_code_cell(workspace, perms):
    (workspace / "e.ipynb").write_text(_nb([_code_cell("x = 1\n")]))
    _doc_executor.execute("notebook_read", {"path": "e.ipynb"}, perms)
    out = json.loads(_doc_executor.execute("notebook_edit", {
        "path": "e.ipynb", "cell_idx": 0, "mode": "replace",
        "cell_type": "code", "source": FIELD_COMMENT + "\nx = 2\n",
    }, perms))
    assert out["status"] == "ok"
    assert NOTE_MARKER in out["note"]


def test_notebook_edit_markdown_cell_is_not_flagged(workspace, perms):
    (workspace / "md.ipynb").write_text(_nb([_code_cell("x = 1\n")]))
    _doc_executor.execute("notebook_read", {"path": "md.ipynb"}, perms)
    out = json.loads(_doc_executor.execute("notebook_edit", {
        "path": "md.ipynb", "cell_idx": 0, "mode": "insert_after",
        "cell_type": "markdown",
        "source": "## Anleitung\nDieser Abschnitt erklaert die Bedienung.\n",
    }, perms))
    assert out["status"] == "ok"
    assert "note" not in out


def test_notebook_edit_english_code_cell_has_no_note(workspace, perms):
    (workspace / "en.ipynb").write_text(_nb([_code_cell("x = 1\n")]))
    _doc_executor.execute("notebook_read", {"path": "en.ipynb"}, perms)
    out = json.loads(_doc_executor.execute("notebook_edit", {
        "path": "en.ipynb", "cell_idx": 0, "mode": "replace",
        "cell_type": "code", "source": "# recompute the value\nx = 2\n",
    }, perms))
    assert "note" not in out


# ---------------------------------------------------------------------------
# Cost caps and robustness
# ---------------------------------------------------------------------------


def test_oversized_content_is_capped_without_raising(workspace, perms):
    lines = ["x{} = {}  # keep this value".format(i, i) for i in range(6000)]
    lines.append(FIELD_COMMENT)
    content = "\n".join(lines) + "\n"
    out = _write(perms, "big.py", content)
    assert "File created" in out
    # Beyond the scan cap -> silently not reported, never raised.
    assert NOTE_MARKER not in out
    assert (workspace / "big.py").read_text(encoding="utf-8") == content


def test_at_most_three_findings_are_listed(workspace, perms):
    content = "\n".join(
        "# Zeile {} wird nicht mehr gebraucht".format(i) for i in range(20))
    out = _write(perms, "many.py", content + "\n")
    assert 0 < _note(out).count("line ") <= 3


def test_binary_ish_content_does_not_raise(workspace, perms):
    content = "\x00\x01\x02 � '''\" # \x00\n" * 5
    out = _write(perms, "weird.py", content)
    assert "File created" in out


@pytest.mark.parametrize("text", [
    "",
    '"""',
    "'''unterminated docstring mit Umlaut ä",
    "/* unbalanced",
    "# ok\n\"\"\"\n" * 3,
    "﻿# leading bom\n",
    "x = '\\\\' # trailing\n",
])
def test_scanner_never_raises_on_malformed_fragments(text):
    for suffix in (".py", ".js", ".sh", ".sql", ".css"):
        assert isinstance(_lang_scan_text(text, suffix), list)


def test_malformed_notebook_json_does_not_raise(workspace, perms):
    out = _write(perms, "broken.ipynb", "{not json at all")
    assert "File created" in out
    assert NOTE_MARKER not in out


def test_helper_returns_empty_string_on_bad_input():
    assert _language_hint_for_write(Path("a.py"), None) == ""
    assert _language_hint_for_write(Path("a.py"), 12345) == ""
    assert _language_hint_for_write(Path("no_suffix"), FIELD_COMMENT) == ""


# ---------------------------------------------------------------------------
# Marker precision
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("segment", [
    " Tetris-Spiel mit Voilà-Integration.",
    " Project root zur Path hinzufügen",
    " Wert wird neu berechnet",
    " groesse der Datei pruefen",
])
def test_marker_hits_on_german(segment):
    assert _lang_has_marker(segment) is True


@pytest.mark.parametrize("segment", [
    " Add the project root to sys.path",
    " MIT License, see LICENSE",
    " DAS and NAS storage back ends",
    " ΔE in kcal/mol, α and β orbitals",
    " see https://example.org/über-uns for details",
    " Author: Jürgen Weber",
    " the die() helper aborts on a fatal error",
    " a fast path for the war-room dashboard",
    " handle the rare case where a tag is missing",
])
def test_no_marker_on_english(segment):
    assert _lang_has_marker(segment) is False


@pytest.mark.parametrize("segment", [
    # Measured against this repository's own English sources: these two
    # phrasings produced every false positive of the first draft.
    " ``None`` falls back to the module defaults below.",
    " Gate ON: the parser falls through to the legacy branch.",
    " Bondi van der Waals radii for clash detection.",
    " Relieves van-der-Waals clashes by rotating distal sub-trees.",
    " Covalent radii (Pyykkoe 2009, single-bond values).",
    " Covalent radii (Pyykkö 2009, single-bond values in Å).",
    " Hückel MO theory next to the Schrödinger equation.",
    " Löwdin population analysis, Möller-Plesset correction.",
])
def test_no_marker_on_measured_false_positive_classes(segment):
    assert _lang_has_marker(segment) is False


def test_shell_and_js_comments_are_covered(workspace, perms):
    sh = _write(perms, "s.sh", "#!/bin/sh\n# Datei wird geloescht\nrm -f x\n")
    js = _write(perms, "s.js", "// Wert wird gesetzt\nconst x = 1;\n")
    css = _write(perms, "s.css", "/* Farbe für den Hintergrund */\na{color:red}\n")
    assert NOTE_MARKER in sh
    assert NOTE_MARKER in js
    assert NOTE_MARKER in css
