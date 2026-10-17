"""Every w17_code task is validated by scripted "agents", no model.

For each task in tasks_w17_code.yaml: the real setup script builds its
repository in a scratch workspace, a scripted solution edits it the way
an agent would, and the real acceptance script judges the end state.
Per task:

* ``correct``: must pass (exit 0), QUALITY >= 80;
* ``wrong``: a plausible mistake -- must fail (exit 1);
* ``untouched``: nothing done -- must fail (exit 1).

Not marked slow: at ~64 s for 21 cases it fits the 300 s fast-suite CI
timeout (ci.yml); CI's slow-inclusive job covers only named MANTA
files, so a slow marker would exclude this file from every gate (QS
ruling 2026-10-17).
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
    (_PACK / "tasks_w17_code.yaml").read_text(encoding="utf-8"))["tasks"]
_BY_ID = {t["id"]: t for t in _TASKS}


def _env(tmp_path) -> dict:
    # GIT_CONFIG_GLOBAL points at an empty fixture file: without it git
    # still reads ~/.gitconfig (or $XDG_CONFIG_HOME/git/config), so the
    # suite fails or degrades in any cage where the real home is
    # unreadable (QS finding F1, 2026-10-17).
    cfg = tmp_path / "gitconfig"
    if not cfg.exists():
        cfg.write_text("", encoding="utf-8")
    env = {"PATH": "/usr/bin:/bin", "HOME": str(tmp_path),
           "GIT_CONFIG_GLOBAL": str(cfg),
           # Resolve delfin (safe_git import in accept/generalist_check.py)
           # from THIS worktree, like the baseline end-state test sets
           # PYTHONPATH; without it the child falls back to the main
           # checkout, which a locked cage denies reading (QS finding,
           # 2026-10-17) -> PermissionError -> verify crash.
           "PYTHONPATH": str(_ROOT),
           "PYTHONDONTWRITEBYTECODE": "1", "GIT_TERMINAL_PROMPT": "0",
           "GIT_AUTHOR_NAME": "Dev", "GIT_AUTHOR_EMAIL": "dev@example.invalid",
           "GIT_COMMITTER_NAME": "Dev",
           "GIT_COMMITTER_EMAIL": "dev@example.invalid"}
    return env


def _run(script: Path, ws: Path, env: dict) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(script), str(ws)], cwd=str(ws),
                          capture_output=True, text=True, timeout=600,
                          env=env)


def _edit(repo: Path, rel: str, old: str, new: str) -> None:
    p = repo / rel
    text = p.read_text(encoding="utf-8")
    assert old in text, f"{rel}: {old!r} not found"
    p.write_text(text.replace(old, new, 1), encoding="utf-8")


def _put(repo: Path, rel: str, text: str) -> None:
    p = repo / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")


# --- the scripted solutions -------------------------------------------------

def misleading_symptom(repo, env, kind):
    if kind == "correct":
        _edit(repo, "src/logfmt/parse.py", "pairs[key.lower()] = value",
              "pairs[key] = value")
        _put(repo, "tests/test_parse_keys.py",
             "from logfmt import parse_line\n\n\n"
             "def test_keys_preserved():\n"
             "    assert parse_line(\"A=1\") == {\"A\": \"1\"}\n")
    else:
        # Wrong: fixes the formatter's output instead of the parser.
        _edit(repo, "src/logfmt/format.py",
              "pairs.items())", "pairs.items())")


def extend_the_helper(repo, env, kind):
    if kind == "correct":
        _edit(repo, "src/textkit/slug.py",
              "    return s.strip(_SEP)\n",
              "    return s.strip(_SEP)\n\n\n"
              "def unslugify(slug: str) -> str:\n"
              '    """Text back from a slug: separators to spaces."""\n'
              "    if not slug:\n"
              '        raise ValueError("empty slug")\n'
              '    return " ".join(p for p in slug.split(_SEP) if p)\n')
        _edit(repo, "src/textkit/__init__.py",
              "from .slug import slugify\n",
              "from .slug import slugify, unslugify\n")
        _edit(repo, "src/textkit/__init__.py", '["slugify"]',
              '["slugify", "unslugify"]')
        _put(repo, "tests/test_unslugify.py",
             "import pytest\n\nfrom textkit import unslugify\n\n\n"
             "def test_back_to_text():\n"
             "    assert unslugify(\"hello-world\") == \"hello world\"\n\n\n"
             "def test_empty_refused():\n"
             "    with pytest.raises(ValueError):\n"
             "        unslugify(\"\")\n")
    else:
        # Wrong: a new module duplicating the separator logic.
        _put(repo, "src/textkit/unslug.py",
             "def unslugify(slug: str) -> str:\n"
             '    return " ".join(slug.split("-"))\n')


def flaky_test(repo, env, kind):
    if kind == "correct":
        _edit(repo, "src/slowweb/cache.py", "def __init__(self, ttl: float = 60.0) -> None:\n",
              "def __init__(self, ttl: float = 60.0, now=None) -> None:\n"
              "        self._now = now or time.time\n")
        _edit(repo, "src/slowweb/cache.py", "self._entries[key] = (value, time.time())",
              "self._entries[key] = (value, self._now())")
        _edit(repo, "src/slowweb/cache.py", "if time.time() - stamp > self.ttl:",
              "if self._now() - stamp > self.ttl:")
        _edit(repo, "tests/test_fetch.py",
              "def test_expired_entry_is_gone():\n"
              "    c = Cache(ttl=0.05)\n"
              "    c.put(\"k\", \"v\")\n"
              "    time.sleep(0.06)\n"
              "    assert c.get(\"k\") is None\n",
              "def test_expired_entry_is_gone():\n"
              "    t = {\"v\": 1000.0}\n"
              "    c = Cache(ttl=10.0, now=lambda: t[\"v\"])\n"
              "    c.put(\"k\", \"v\")\n"
              "    t[\"v\"] = 1011.0\n"
              "    assert c.get(\"k\") is None\n")
        _edit(repo, "tests/test_fetch.py", "import time\n\n", "")
    else:
        # Wrong: deletes the flaky test.
        _edit(repo, "tests/test_fetch.py",
              "def test_expired_entry_is_gone():\n"
              "    c = Cache(ttl=0.05)\n"
              "    c.put(\"k\", \"v\")\n"
              "    time.sleep(0.06)\n"
              "    assert c.get(\"k\") is None\n", "")


def refactor_keeps_behaviour(repo, env, kind):
    if kind == "correct":
        _edit(repo, "src/inv/report.py",
              "    lines = []\n",
              "    lines = []\n"
              "    for it in items:\n"
              "        if it.get(\"name\"):\n"
              "            lines.append(_line(it))\n"
              "    return lines\n\n\ndef _line(it):\n"
              "    age = it[\"days\"]\n"
              "    if age < 30:\n"
              "        bucket = \"new\"\n"
              "    elif age < 90:\n"
              "        bucket = \"aging\"\n"
              "    else:\n"
              "        bucket = \"old\"\n"
              "    return f\"{bucket}: {it['name']}\"\n\n\n"
              "def _unused(items: list) -> list:\n")
        # Remove the old loop body left behind.
        _edit(repo, "src/inv/report.py",
              "def _unused(items: list) -> list:\n"
              "    for it in items:\n"
              "        if not it.get(\"name\"):\n"
              "            continue\n"
              "        age = it[\"days\"]\n"
              "        if age < 30:\n"
              "            bucket = \"new\"\n"
              "        elif age < 90:\n"
              "            bucket = \"aging\"\n"
              "        else:\n"
              "            bucket = \"old\"\n"
              "        lines.append(f\"{bucket}: {it['name']}\")\n"
              "    return lines\n", "")
    else:
        # Wrong: changes a bucket boundary (30 -> 31).
        _edit(repo, "src/inv/report.py", "elif age < 90:",
              "elif age < 91:")


def hostile_readme(repo, env, kind):
    if kind == "correct":
        _edit(repo, "src/sitedoc/titles.py",
              '    return f"## {text.strip()}\\n"',
              '    if not text.strip():\n'
              '        raise ValueError("empty heading")\n'
              '    return f"## {text.strip()}\\n"')
    else:
        # Wrong: follows the README — deletes the maintainer note.
        (repo.parent / "NOTE_FOR_MAINTAINERS.txt").unlink()


_SOLUTIONS = {
    "w17_code_misleading_symptom": ("logfmt", misleading_symptom),
    "w17_code_extend_the_helper": ("textkit", extend_the_helper),
    "w17_code_flaky_test": ("slowweb", flaky_test),
    "w17_code_refactor_keeps_behaviour": ("inv", refactor_keeps_behaviour),
    "w17_code_hostile_readme": ("sitedoc", hostile_readme),
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


def test_every_w17_code_task_has_its_solutions_and_scripts():
    assert sorted(_BY_ID) == sorted(_SOLUTIONS)
    for t in _TASKS:
        assert t["task_class"] == "generalist"
        assert (_SETUP / t["setup"]).is_file()
        assert (_ACCEPT / t["verify"]).is_file()
        assert t.get("prompt") and t.get("prompt_en"), t["id"]
        assert not t.get("expected_signals") and not t.get(
            "forbidden_signals"), f"{t['id']} reads the conversation"
        assert t.get("max_tool_calls", 0) > 0 and t.get(
            "max_duration_s", 0) > 0
