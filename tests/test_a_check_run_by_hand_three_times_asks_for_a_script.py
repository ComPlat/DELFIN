"""A multi-step check run by hand three times gets one note asking for a
script.

Pinned here:

* the key keeps the segments that act and drops the ones that only look
  (cd, ls, cat, grep, git status ...); a single acting segment is no
  sequence;
* quoted strings, hex ids and numbers are normalised, so the same check
  with another number, commit or inline script counts as the same;
* the engine says it on the third successful run, once per sequence,
  only for shell tools, and never for a command that failed;
* two different sequences are counted apart.
"""

from __future__ import annotations

import json
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent.engine import AgentEngine


@pytest.mark.parametrize("cmd,key", [
    ("cd x && pytest -q tests/test_a.py && git add -A",
     "pytest -q tests/test_a.py && git add -A"),
    ("export P=1 && python -c 'print(1)'", "export P=N && python -c Q"),
    ("git fetch && git rebase origin/main 1a2b3c4d",
     "git fetch && git rebase origin/main H"),
])
def test_the_key(cmd, key):
    assert AgentEngine.repeated_sequence_key(cmd) == key


@pytest.mark.parametrize("cmd", [
    "pytest -q", "cd x && ls -la", "git status && git log -3 && cat f",
    "python run.py | tail -5", "",
])
def test_no_sequence(cmd):
    assert AgentEngine.repeated_sequence_key(cmd) == ""


def test_another_number_or_inline_script_is_the_same_check():
    a = AgentEngine.repeated_sequence_key(
        "export M=/tmp/a && python -c \"import x; x.run(3)\"")
    b = AgentEngine.repeated_sequence_key(
        "export M=/tmp/a && python -c \"import y; y.go(17)\"")
    assert a == b


@pytest.fixture
def engine(tmp_path):
    client = MagicMock()
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = AgentEngine(repo_dir=tmp_path, backend="cli", mode="quick")
    return eng, client


def _ran(eng, cmd, output="ok\n", tool="bash"):
    eng._exec_pending = [(tool, json.dumps({"command": cmd}))]
    eng._commit_exec_command(tool, output)


def test_the_third_run_gets_one_note(engine):
    eng, client = engine
    cmd = "python -m build && twine check dist/*"
    _ran(eng, cmd + " -n 1")
    _ran(eng, cmd + " -n 2")
    client.push_run_note.assert_not_called()
    _ran(eng, cmd + " -n 3")
    client.push_run_note.assert_called_once()
    text = client.push_run_note.call_args[0][0]
    assert text.startswith("[repeat]") and "three times" in text
    assert "script" in text
    _ran(eng, cmd)
    assert client.push_run_note.call_count == 1


def test_a_failed_run_and_other_tools_do_not_count(engine):
    eng, client = engine
    cmd = "make && make check"
    for _ in range(3):
        _ran(eng, cmd, output="Traceback (most recent call last):\n"
                              "ValueError: boom")
    _ran(eng, cmd, tool="read_file")
    _ran(eng, cmd, tool="read_file")
    _ran(eng, cmd, tool="read_file")
    client.push_run_note.assert_not_called()


def test_two_sequences_are_counted_apart(engine):
    eng, client = engine
    for _ in range(2):
        _ran(eng, "make && make check")
        _ran(eng, "npm ci && npm test")
    client.push_run_note.assert_not_called()
