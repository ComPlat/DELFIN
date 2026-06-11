"""Regression: slash-command arguments must preserve case.

Jerome's bug (2026-06-11): ``/grant /pfs/.../Porpoise always`` failed
with "path does not exist: /pfs/.../porpoise" — the dashboard slash
dispatcher lowercased the ENTIRE command line (``cmd = text.lower()``)
and parsed the path argument out of that, so a case-sensitive Linux
path was destroyed. Path/command arguments must come from the
original-case text; only the command token is matched case-insensitively.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions


def test_add_extra_dir_is_case_sensitive(tmp_path):
    """Proves WHY case must be preserved: a lowercased path won't resolve."""
    real = tmp_path / "Porpoise"
    real.mkdir()
    perms = KitToolPermissions(workspace=tmp_path, mode="default")

    # Correct case resolves.
    resolved = perms.add_extra_dir(str(real))
    assert resolved == real.resolve()

    # Lowercased case does not exist -> ValueError (Jerome's symptom).
    with pytest.raises(ValueError, match="does not exist"):
        perms.add_extra_dir(str(tmp_path / "porpoise"))


def _dispatcher_source() -> str:
    src = (Path(__file__).resolve().parent.parent
           / "delfin" / "dashboard" / "tab_agent.py").read_text(encoding="utf-8")
    i = src.find("def _handle_slash_command")
    assert i != -1, "slash dispatcher not found"
    return src[i:i + 80_000]


def test_grant_parses_argument_from_original_text():
    """/grant must slice its path argument from ``text``, not ``cmd``."""
    src = _dispatcher_source()
    g = src.find('if cmd == "/grant"')
    assert g != -1
    block = src[g:g + 700]
    assert 'text.strip()[len("/grant")' in block, (
        "/grant must parse the path from original-case text"
    )
    assert 'arg = cmd[len("/grant")' not in block, (
        "/grant must NOT parse the path from the lowercased cmd"
    )


def test_watch_parses_folder_from_original_text():
    """/watch add <jobid> <folder> — folder is a case-sensitive path."""
    src = _dispatcher_source()
    w = src.find('arg = text.strip()[len("/watch")')
    assert w != -1, "/watch must parse its argument from original-case text"
    # subcommand still matched case-insensitively
    assert "parts[0].lower()" in src[w:w + 200]


# --- Error-detection self-correction (Jerome's "where is the error?" bug) ---

def test_missing_file_hint_lists_real_outputs(tmp_path):
    """A guessed name like orca.out must surface the REAL output files."""
    from delfin.dashboard.tab_agent import _calc_missing_file_hint
    calc = tmp_path
    folder = calc / "MarcDy"
    folder.mkdir()
    (folder / "MarcDy_opt.out").write_text("some output without a normal end")
    (folder / "CONTROL.txt").write_text("input_file=input.txt\n")

    hint = _calc_missing_file_hint("orca.out", folder, calc)
    assert "Not a file: orca.out" in hint
    assert "MarcDy_opt.out" in hint                 # real name surfaced
    assert "/analyze errors" in hint                # points at the right tool
    assert "calc/MarcDy" in hint


def test_missing_file_hint_points_into_subfolders(tmp_path):
    """At a parent level with only subfolders, steer cd into them."""
    from delfin.dashboard.tab_agent import _calc_missing_file_hint
    calc = tmp_path
    (calc / "JobA").mkdir()
    (calc / "JobB").mkdir()
    hint = _calc_missing_file_hint("orca.out", calc, calc)
    assert "only subfolders" in hint
    assert "JobA" in hint and "JobB" in hint
    assert "/calc cd" in hint


def test_missing_file_hint_marks_ok_vs_error(tmp_path):
    """.out files are tagged so the agent sees which one failed."""
    from delfin.dashboard.tab_agent import _calc_missing_file_hint
    calc = tmp_path
    folder = calc / "run"
    folder.mkdir()
    # has_ok_marker keys off ORCA's normal-termination line.
    (folder / "good.out").write_text("...\n****ORCA TERMINATED NORMALLY****\n")
    (folder / "bad.out").write_text("...\nABORTING THE RUN\n")
    hint = _calc_missing_file_hint("orca.out", folder, calc)
    assert "good.out [OK]" in hint
    assert "bad.out [INCOMPLETE/ERROR]" in hint
