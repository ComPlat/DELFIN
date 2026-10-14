"""The citation check judges what an answer ASSERTS exists, where it exists.

Field data (33 sessions, 2026-09-16..10-14): 47 "[Verify]" rounds on
citations, 26 of them "path does not exist", and about half of those were
false: a file the agent said it WOULD create ("schreibe ich die Tests,
lasse `pytest tests/test_delfin_doctor.py` laufen"), a file it said does
NOT exist ("`delfin/doctor.py` existiert weiterhin nirgends"), a module
named as an example ("z.B. `delfin/chemdarwin_space.py`"), and files that
live in a sibling session's worktree of the same repository
(`report_costs.py`, `.delfin/worktrees/delfin-wt-...`). Each false alarm
is a forced correction round and a "Korrektur der Zitate" in front of the
user. A real miss, `gpu-cpu-slurm.md` for the real
`infrastruktur-gpu-cpu-slurm.md`, must still be caught.
"""

from __future__ import annotations

import subprocess

import pytest

from delfin.agent.verify_guard import (
    scan_for_ungrounded_code_claims,
    scan_for_ungrounded_location_claims,
)
from tests.conftest import child_env


def _hard(text, root, **kw):
    return [f.path for f in scan_for_ungrounded_code_claims(
        text, repo_root=root, **kw) if f.kind == "nonexistent"]


@pytest.mark.parametrize("text", [
    "Sobald Session A antwortet, schreibe ich die Tests, lasse "
    "`pytest tests/test_delfin_doctor.py -q` grün laufen.",
    "Next I will create `delfin/agent/report_costs.py` with collect().",
    "Soll die Funktion als eigenes Modul (z.B. `delfin/chemdarwin_space.py`) "
    "angelegt werden?",
    "Noch nichts: `delfin/doctor.py` existiert weiterhin nirgends.",
    "Tasks 1-3 are marked completed but `delfin/fitness_oracle.py` does not "
    "exist in the workspace.",
    "The batch failed because `/tmp/covered_tests.txt` doesn't exist inside "
    "the job's sandbox.",
    "Sobald du grünes Licht gibst, legt Session B `delfin/fitness_oracle.py` "
    "an.",
    # Asked about, not asserted; offered, not claimed.
    "Gibt es `delfin/doctor.py` schon auf einem Branch, und was steht in "
    "der Schnittstelle?Session A hat noch nichts committet.",
    "Wenn du magst, fülle ich inzwischen das leere "
    "`02_research/frontier-methods.md` mit der Recherche.",
])
def test_a_file_named_as_future_or_absent_is_not_called_invented(
        text, tmp_path):
    (tmp_path / "delfin").mkdir()
    assert _hard(text, tmp_path) == []


@pytest.mark.parametrize("text, path", [
    ("The SLURM notes are in `04_infrastructure/gpu-cpu-slurm.md`.",
     "04_infrastructure/gpu-cpu-slurm.md"),
    ("The loop is at `delfin/agent/runner_loop.py:120`.",
     "delfin/agent/runner_loop.py"),
    ("I created `delfin/agent/report_costs.py` and it passes.",
     "delfin/agent/report_costs.py"),
])
def test_a_file_asserted_to_exist_that_does_not_is_still_caught(
        text, path, tmp_path):
    (tmp_path / "delfin" / "agent").mkdir(parents=True)
    (tmp_path / "04_infrastructure").mkdir()
    (tmp_path / "04_infrastructure" / "infrastruktur-gpu-cpu-slurm.md"
     ).write_text("x")
    assert _hard(text, tmp_path) == [path]


def _git(args, cwd, tmp_path):
    env = {**child_env(tmp_path),
           "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@example.invalid",
           "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@example.invalid"}
    subprocess.run(["git", *args], cwd=cwd, check=True, env=env,
                   capture_output=True)


def test_a_file_in_a_sibling_worktree_is_not_called_invented(tmp_path):
    main = tmp_path / "repo"
    main.mkdir()
    (main / "README.md").write_text("r\n")
    _git(["init", "-q"], main, tmp_path)
    _git(["add", "."], main, tmp_path)
    _git(["commit", "-qm", "base"], main, tmp_path)
    mine = tmp_path / "wt-mine"
    theirs = tmp_path / "wt-theirs"
    _git(["worktree", "add", "-q", str(mine)], main, tmp_path)
    _git(["worktree", "add", "-q", str(theirs)], main, tmp_path)
    (theirs / "delfin" / "agent").mkdir(parents=True)
    (theirs / "delfin" / "agent" / "report_costs.py").write_text("x = 1\n")
    text = ("Session B built `delfin/agent/report_costs.py` with the same "
            "collect() contract.")
    flags = scan_for_ungrounded_code_claims(text, repo_root=mine)
    assert [f.kind for f in flags] == ["unread"], flags
    # A file in no worktree at all is still invented.
    flags = scan_for_ungrounded_code_claims(
        "See `delfin/agent/nowhere.py` for it.", repo_root=mine)
    assert [f.kind for f in flags] == ["nonexistent"]


def test_a_read_from_an_earlier_turn_still_grounds_a_location(tmp_path):
    """The engine hands the scanners the session's ledger, not the turn's:
    pinned here at the scanner, which is where a per-turn ledger would
    reappear as the 'not backed by a read this turn' false alarm."""
    observed = {str(tmp_path / "delfin" / "dashboard" / "tab_agent.py")}
    assert scan_for_ungrounded_location_claims(
        "The save hook is at `tab_agent.py:7222`.",
        observed_files=observed) == []
