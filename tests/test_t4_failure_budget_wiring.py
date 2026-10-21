"""The third identical failing call is told to stop and ask.

The executor feeds every call's outcome into a process-wide, memory-bounded
FailureBudget. A call that fails the same way a third time in a row gets a
stop_and_ask note in its JSON error; a success in between resets the streak.
"""
from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as A


@pytest.fixture(autouse=True)
def _fresh_budget(monkeypatch):
    monkeypatch.setattr(A._DocToolExecutor, "_T4_FAIL_BUDGET", None, raising=False)


def _read(tmp_path, name):
    perms = A.KitToolPermissions(workspace=tmp_path, mode="acceptEdits")
    out = A._DocToolExecutor().execute("read_file", {"path": str(tmp_path / name)}, perms)
    try:
        return json.loads(out)
    except ValueError:
        return {"raw": out}


def test_the_third_identical_failure_says_stop_and_ask(tmp_path):
    first = _read(tmp_path, "missing.txt")
    second = _read(tmp_path, "missing.txt")
    third = _read(tmp_path, "missing.txt")
    assert "error" in first and "stop_and_ask" not in first
    assert "stop_and_ask" not in second
    assert "stop_and_ask" in third, third


def test_a_success_in_between_resets_the_streak(tmp_path):
    _read(tmp_path, "late.txt")
    _read(tmp_path, "late.txt")
    (tmp_path / "late.txt").write_text("now here", encoding="utf-8")
    ok = _read(tmp_path, "late.txt")
    assert "error" not in ok
    (tmp_path / "late.txt").unlink()
    again = _read(tmp_path, "late.txt")
    assert "stop_and_ask" not in again
