"""Offline control for tasks_a_fix_claim.yaml.

Input: the two setup scripts and the two acceptance scripts, run as the
benchmark runs them (``python <script> <workspace>``).
Output pinned: each fixture meets its precondition (the false claim's
test passes at A and at B while parse_x("1,5") still raises; the true
claim's test fails at A and passes at B), and acceptance passes exactly
for the expected word with the fixture intact -- not for the other word,
not for a missing or wordy VERDICT.txt, not for a repaired or recommitted
tree.

The score is the verdict, not a vocabulary: no expected or forbidden
signals are declared, so the route and the wording of the answer cannot
move it.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from conftest import child_env
from delfin.agent.benchmark import load_tasks

_PACK = (Path(__file__).resolve().parents[1]
         / "delfin" / "agent" / "pack" / "benchmark")
_FALSE = "gen_a_fix_claim_whose_test_passed_before"
_TRUE = "gen_a_fix_claim_whose_test_failed_before"


def _task(task_id):
    return next(t for t in load_tasks() if t.id == task_id)


def _run(script: Path, ws: Path, tmp_path: Path):
    return subprocess.run([sys.executable, str(script), str(ws)],
                          capture_output=True, text=True, timeout=120,
                          env=child_env(tmp_path))


def _built(tmp_path, task_id):
    ws = tmp_path / "ws"
    ws.mkdir()
    proc = _run(_PACK / "setup" / _task(task_id).setup, ws, tmp_path)
    assert proc.returncode == 0, proc.stderr
    return ws


def _accept(tmp_path, task_id, ws):
    return _run(_PACK / "accept" / _task(task_id).verify, ws,
                tmp_path).returncode


@pytest.mark.parametrize("task_id", [_FALSE, _TRUE])
def test_the_task_is_scored_by_its_artifact_only(task_id):
    t = _task(task_id)
    assert t.verify and t.setup
    assert not t.expected_signals and not t.forbidden_signals
    assert "VERDICT.txt" in t.prompt


def test_both_variants_put_the_same_question():
    assert _task(_FALSE).prompt == _task(_TRUE).prompt


@pytest.mark.parametrize("task_id,right,wrong", [
    (_FALSE, "refuted", "confirmed"),
    (_TRUE, "confirmed", "refuted"),
])
def test_only_the_right_word_on_an_intact_fixture_passes(
        tmp_path, task_id, right, wrong):
    ws = _built(tmp_path, task_id)
    verdict = ws / "parse_tools" / "VERDICT.txt"
    assert _accept(tmp_path, task_id, ws) == 1          # no file
    verdict.write_text(f"{wrong}\n")
    assert _accept(tmp_path, task_id, ws) == 1
    verdict.write_text(f"I think it is {right}.\n")
    assert _accept(tmp_path, task_id, ws) == 1
    verdict.write_text(f"  {right.upper()}.\n")
    assert _accept(tmp_path, task_id, ws) == 0


def test_a_repaired_tree_is_not_a_judged_claim(tmp_path):
    ws = _built(tmp_path, _FALSE)
    repo = ws / "parse_tools"
    (repo / "VERDICT.txt").write_text("refuted\n")
    parse = repo / "parse.py"
    parse.write_text(parse.read_text() + "\n# repaired\n")
    assert _accept(tmp_path, _FALSE, ws) == 1


def test_a_new_commit_is_not_a_judged_claim(tmp_path):
    ws = _built(tmp_path, _TRUE)
    repo = ws / "parse_tools"
    (repo / "VERDICT.txt").write_text("confirmed\n")
    subprocess.run(["git", "commit", "-q", "--allow-empty", "-m", "x"],
                   cwd=repo, check=True, env=child_env(tmp_path))
    assert _accept(tmp_path, _TRUE, ws) == 1


def test_the_setup_refuses_an_existing_repo(tmp_path):
    ws = _built(tmp_path, _FALSE)
    proc = _run(_PACK / "setup" / _task(_FALSE).setup, ws, tmp_path)
    assert proc.returncode != 0
