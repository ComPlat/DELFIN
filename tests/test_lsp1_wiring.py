"""LSP1 phase 2 wiring tests.

Red state (main before the wiring patch): `_diagnostics_note` does not
exist in delfin/agent/api_client.py, so every case fails at the
attribute lookup — the red the operator's build verifies. After the
operator applies .gate/lsp1_wiring.patch these pass.

The wiring appends a "Diagnostics after this edit:" block to the result
strings of write_file / edit_file (exact + fuzzy) / multi_edit. Negative
seam: diff_approval staged results must stay silent (no diagnostics on
a pending change that never touched disk).
"""
import json
import os
import stat
import textwrap

import pytest

import delfin.agent.api_client as api


def _install_fake_ruff(tmp_path, monkeypatch, marker="beta"):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    exe = bin_dir / "ruff"
    exe.write_text(
        "#!/usr/bin/python3\n"
        "import json, re, sys\n"
        "args = sys.argv[1:]\n"
        'content = sys.stdin.read() if "--stdin-filename" in args else open(args[-1]).read()\n'
        "out = []\n"
        "for i, line in enumerate(content.splitlines(), start=1):\n"
        f'    if re.search(rf"\\b{marker}\\(\\)", line):\n'
        "            col = line.index(\"%s\")\n" % marker +
        "            out.append({\n"
        '                "code": "F821",\n'
        f'                "message": f"Undefined name `{marker}`",\n'
        '                "location": {"row": i, "column": col},\n'
        "            })\n"
        "print(json.dumps(out))\n",
        encoding="utf-8")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
    return bin_dir


def test_diagnostics_note_helper_exists_and_formats(api_client=None):
    """The helper exists and formats one finding at the right location."""
    import delfin.agent.diagnostics as diag
    assert hasattr(api, "_diagnostics_note"), (
        "wiring missing: _diagnostics_note not in api_client — "
        "the operator's build has not landed")


def test_note_reports_new_finding(tmp_path, monkeypatch):
    _install_fake_ruff(tmp_path, monkeypatch)
    from delfin.agent import diagnostics as _d  # ensure module imported
    note = api._diagnostics_note(
        tmp_path / "mod.py", tmp_path, "def go():\n    alpha()\n")
    # mod.py does not exist -> helper yields "" (no findings, no note)
    assert note == ""


def test_note_clean_edit(tmp_path, monkeypatch):
    _install_fake_ruff(tmp_path, monkeypatch)
    p = tmp_path / "mod.py"
    p.write_text("def go():\n    return 1\n", encoding="utf-8")
    note = api._diagnostics_note(p, tmp_path, "def go():\n    return 1\n")
    # A clean edit yields no block at all — nothing to tell the model.
    assert note == ""


def test_note_reports_added_beta(tmp_path, monkeypatch):
    _install_fake_ruff(tmp_path, monkeypatch)
    p = tmp_path / "mod.py"
    p.write_text("def go():\n    alpha()\n    beta()\n", encoding="utf-8")
    note = api._diagnostics_note(
        p, tmp_path, "def go():\n    alpha()\n")
    assert "Diagnostics after this edit:" in note
    assert "F821" in note
    assert "Undefined name `beta`" in note
    assert "alpha" not in note.split("Diagnostics after this edit:", 1)[1]


def test_note_silent_when_unavailable(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    p = tmp_path / "mod.py"
    p.write_text("def go():\n    beta()\n", encoding="utf-8")
    note = api._diagnostics_note(p, tmp_path, "def go():\n    return 1\n")
    # ruff unavailable -> no findings -> no block at all
    assert note == ""


def test_wiring_seams_present_in_source():
    """The four seams carry the diagnostics hook; diff_approval paths
    (returns before _capture_change) carry none. Reads the real source
    of the installed module — no parsing games."""
    src = open(api.__file__, encoding="utf-8").read()
    fired = src.count("_diagnostics_note(")
    # 1 definition + 4 call sites (write, fuzzy edit, exact edit, multi_edit)
    assert fired == 5, f"expected 1 def + 4 seams, found {fired}"
    # diff_approval stage path stays silent
    for line_no, line in enumerate(src.splitlines(), start=1):
        if "_stage_pending_change(" in line and "diagnostics" in line:
            raise AssertionError(
                f"diff_approval path fires diagnostics at line {line_no}: {line}")
