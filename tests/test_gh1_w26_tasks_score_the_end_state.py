"""Every GH1 (wave 26) task is validated by scripted "agents", no model.

For each task in tasks_w26_gh1.yaml: the real setup script builds its
repository in a scratch workspace, a scripted solution acts on it the
way an agent would, and the real acceptance script judges the end
state. Per task:

* ``correct``: must pass (exit 0), QUALITY >= 80;
* ``wrong``: a plausible mistake -- must fail (exit 1);
* ``untouched``: nothing done -- must fail (exit 1).

Pattern: tests/test_c1_w17_code_tasks_score_the_end_state.py.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest
import yaml

_ROOT = Path(__file__).resolve().parents[1]
_PACK = _ROOT / "delfin" / "agent" / "pack" / "benchmark"
_SETUP = _PACK / "setup"
_ACCEPT = _PACK / "accept"

_TASKS = yaml.safe_load(
    (_PACK / "tasks_w26_gh1.yaml").read_text(encoding="utf-8"))["tasks"]
_BY_ID = {t["id"]: t for t in _TASKS}


def _env(tmp_path) -> dict:
    # GIT_CONFIG_GLOBAL points at an empty fixture file: without it git
    # still reads ~/.gitconfig, so the suite degrades in any cage where
    # the real home is unreadable (cf. the w17_code test's QS finding F1).
    cfg = tmp_path / "gitconfig"
    if not cfg.exists():
        cfg.write_text("", encoding="utf-8")
    home = tmp_path / "home"
    home.mkdir(exist_ok=True)
    return {"HOME": str(home),
            "GIT_CONFIG_GLOBAL": str(cfg),
            "PYTHONPATH": str(_ROOT),
            "PYTHONDONTWRITEBYTECODE": "1",
            "GIT_TERMINAL_PROMPT": "0",
            "GIT_AUTHOR_NAME": "Dev",
            "GIT_AUTHOR_EMAIL": "dev@example.invalid",
            "GIT_COMMITTER_NAME": "Dev",
            "GIT_COMMITTER_EMAIL": "dev@example.invalid"}


def _run(script: Path, ws: Path, env: dict) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(script), str(ws)],
                          cwd=str(ws), capture_output=True, text=True,
                          timeout=600, env=env)


def _edit(repo: Path, rel: str, old: str, new: str) -> None:
    p = repo / rel
    text = p.read_text(encoding="utf-8")
    assert old in text, f"{rel}: {old!r} not found"
    p.write_text(text.replace(old, new, 1), encoding="utf-8")


def _git(repo: Path, *args: str) -> str:
    import os
    env = dict(os.environ)
    env.update({"GIT_CONFIG_GLOBAL": str(repo.parent / "gitconfig"),
                "GIT_TERMINAL_PROMPT": "0"})
    proc = subprocess.run(["git", *args], cwd=str(repo),
                          capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)}: {proc.stderr[:200]}")
    return proc.stdout.strip()


# --- the scripted solutions ---------------------------------------------------

def pr_ready(repo, env, kind):
    if kind == "correct":
        _git(repo, "checkout", "-q", "-b", "feature/underscore")
        _edit(repo, "src/widgetkit/core.py",
              'return _SEP.sub(" ", label.strip().lower())',
              'return _SEP.sub(" ", label.replace("_", "-").strip().lower())')
        _put(repo, "tests/test_underscore.py",
             "from widgetkit import normalise\n\n\n"
             "def test_underscores_become_spaces():\n"
             "    assert normalise(\"hello_world\") == \"hello world\"\n")
        _git(repo, "add", "-A")
        _git(repo, "commit", "-q", "-m", "feature: underscores to spaces")
    else:
        # Wrong: commits the probe test file along with the feature.
        _git(repo, "checkout", "-q", "-b", "feature/underscore")
        _put(repo, "tests/test_underscore.py",
             "from widgetkit import normalise\n\n\n"
             "def test_underscores_become_spaces():\n"
             "    assert normalise(\"hello_world\") == \"hello world\"\n")
        _put(repo, "tests/test_zz_probe2.py", "def test_probe2():\n"
                                              "    assert True\n")
        _git(repo, "add", "-A")
        _git(repo, "commit", "-q", "-m", "feature + probe (wrong)")


def red_ci(repo, env, kind):
    if kind == "correct":
        _edit(repo, "src/widgetkit/core.py", 'r"[-_\\s]"', 'r"[-_\\s]+"')
        _git(repo, "add", "-A")
        _git(repo, "commit", "-q", "-m", "fix: collapse separator runs")
    else:
        # Wrong: bends the test to expect the broken output.
        _edit(repo, "tests/test_label.py", '== "a b"',
              '== "a   b"')


def cross_review(repo, env, kind):
    answer = repo.parent / "answer.txt"
    if kind == "correct":
        answer.write_text(
            "Reviewed partner/feature (read-only; own branch main "
            "untouched, no merge).\n"
            "src/widgetkit/lookup.py: BLOCKER -- the bare `except:` in "
            "load_config swallows every failure (FileNotFoundError, JSON "
            "errors, anything); the docstring promises defaults only for "
            "a missing file. Fix on the partner branch: catch "
            "FileNotFoundError (and json.JSONDecodeError) explicitly.\n"
            "src/widgetkit/pretty.py: fine, nothing to change.\n",
            encoding="utf-8")
    else:
        # Wrong: accuses the fine module and misses the bare except.
        answer.write_text(
            "pretty.py has a defect: the width computation is wrong. "
            "lookup.py looks fine.\n", encoding="utf-8")


def count_conflict(repo, env, kind):
    if kind == "correct":
        counts = repo / "metrics/counts.py"
        counts.write_text(
            '"""Measured counts per task class; summed at merge time."""\n'
            'COUNTS = {"count": 25}\n', encoding="utf-8")
        _git(repo, "add", "metrics/counts.py")
        _git(repo, "commit", "-q", "--no-edit")
    else:
        # Wrong: picks one side (main's 10) instead of summing.
        counts = repo / "metrics/counts.py"
        counts.write_text(
            '"""Measured counts per task class; summed at merge time."""\n'
            'COUNTS = {"count": 10}\n', encoding="utf-8")
        _git(repo, "add", "metrics/counts.py")
        _git(repo, "commit", "-q", "--no-edit")


def _put(repo: Path, rel: str, text: str) -> None:
    p = repo / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")


_SOLUTIONS = {
    "w26-gh1-01": ("widgetkit", pr_ready),
    "w26-gh1-02": ("widgetkit", red_ci),
    "w26-gh1-03": ("widgetkit", cross_review),
    "w26-gh1-04": ("widgetkit", count_conflict),
}


def _quality(text: str) -> int:
    import re
    m = re.findall(r"^QUALITY: (\d+)$", text, re.M)
    return int(m[-1]) if m else -1


def _attempt(task_id: str, kind: str, tmp_path: Path):
    env = _env(tmp_path)
    task = _BY_ID[task_id]
    ws = tmp_path / "ws"
    ws.mkdir()
    setup = _run(_SETUP / task["setup"], ws, env)
    assert setup.returncode == 0, setup.stdout + setup.stderr
    name, solve = _SOLUTIONS[task_id]
    repo = ws / name
    if kind != "untouched":
        solve(repo, env, kind)
    return setup, _run(_ACCEPT / task["verify"], ws, env)


_CASES = ([(t, "correct") for t in _SOLUTIONS]
          + [(t, "wrong") for t in _SOLUTIONS]
          + [(t, "untouched") for t in _SOLUTIONS])


@pytest.mark.parametrize("task_id,kind", _CASES)
def test_a_scripted_attempt_is_judged_on_the_end_state(task_id, kind,
                                                       tmp_path):
    setup, acc = _attempt(task_id, kind, tmp_path)
    report = acc.stdout + acc.stderr
    want = 0 if kind == "correct" else 1
    assert acc.returncode == want, report
    q = _quality(report)
    assert 0 <= q <= 100, report
    if kind == "correct":
        assert q >= 80, report


@pytest.mark.parametrize("task_id", sorted(_SOLUTIONS))
def test_the_correct_solution_grades_above_the_wrong_one(task_id, tmp_path):
    a = tmp_path / "a"
    b = tmp_path / "b"
    a.mkdir()
    b.mkdir()
    _, good = _attempt(task_id, "correct", a)
    _, bad = _attempt(task_id, "wrong", b)
    assert _quality(good.stdout) > _quality(bad.stdout), (
        good.stdout + "\n---\n" + bad.stdout)


def test_every_gh1_task_has_its_solutions_and_scripts():
    assert sorted(_BY_ID) == sorted(_SOLUTIONS)
    for t in _TASKS:
        assert t["task_class"] == "generalist"
        assert (_SETUP / t["setup"]).is_file()
        assert (_ACCEPT / t["verify"]).is_file()
        assert t.get("prompt") and t.get("prompt_en"), t["id"]
        assert not t.get("expected_signals") and not t.get(
            "forbidden_signals"), f"{t['id']} reads the conversation"
        assert t.get("max_tool_calls", 0) > 0 and t.get(
            "max_duration_s", 0) > 0, t["id"]
