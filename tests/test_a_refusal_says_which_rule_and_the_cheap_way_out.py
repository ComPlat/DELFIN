"""A refusal that is right and unhelpful costs the same as a wrong one.

The write gate refused `> suite_full.log` for want of a prior read_file.
The way out was named — and for a file the command is about to replace
whole, "call read_file on it first" reads as "read this entire log", so
a session on 2026-09-18 hit it twice rather than once.

The rule is unchanged. The message now names the cheap baseline: for a
file being replaced, read_file with limit=1 establishes it.
"""

from __future__ import annotations

import json
import pathlib

import pytest

from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


@pytest.fixture()
def engine(tmp_path):
    perms = KitToolPermissions(workspace=tmp_path)
    perms.mode = "bypassPermissions"
    eng = _DocToolExecutor.__new__(_DocToolExecutor)
    eng._permissions = perms
    return eng, perms, tmp_path


def test_the_write_gate_names_the_cheap_baseline(engine):
    eng, perms, ws = engine
    log = ws / "suite_full.log"
    log.write_text("x" * 5000, encoding="utf-8")

    out = json.loads(eng._execute_bash(
        {"command": f"echo run > {log}", "description": "rerun"}, perms))
    err = str(out.get("error") or "")
    assert "without a prior read_file" in err
    assert "limit=1" in err, (
        "for a file about to be replaced whole, the message read as "
        "'read this entire log first': " + err)


def test_the_write_gate_still_refuses_an_unread_file(engine):
    eng, perms, ws = engine
    existing = ws / "someone_elses.txt"
    existing.write_text("theirs\n", encoding="utf-8")
    out = json.loads(eng._execute_bash(
        {"command": f"echo mine > {existing}", "description": "x"}, perms))
    assert "without a prior read_file" in str(out.get("error"))


def test_a_new_file_needs_no_baseline(engine):
    """Creating one has nothing to clobber — that half must not change."""
    eng, perms, ws = engine
    out = json.loads(eng._execute_bash(
        {"command": f"echo fresh > {ws / 'new.txt'}", "description": "x"},
        perms))
    assert out.get("exit_code") == 0
