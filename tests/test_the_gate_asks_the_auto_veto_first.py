"""The tool gate asks an armed auto veto before anything runs.

Input: ``KitToolPermissions.auto_veto`` -- a callable ``(name, args) ->
refusal text | None`` that auto mode sets for the length of an iteration
-- and a tool call through the single native entry
(``_DocToolExecutor.execute``) or the MCP entry (``_mcp_call_refusal``).
Output: the call's result. A refusal comes back as ``{"error": ...}`` and
the tool does not run; ``None`` (or no veto) leaves the gate exactly as it
was. Sub-agents get a copy of the permissions with the veto in it.
"""

from __future__ import annotations

import dataclasses
import json
import subprocess

import pytest

from delfin.agent import api_client as ac
from delfin.agent.api_client import KitToolPermissions


@pytest.fixture()
def perms(tmp_path):
    (tmp_path / "a.txt").write_text("hello\n", encoding="utf-8")
    return KitToolPermissions(workspace=str(tmp_path), mode="acceptEdits")


def test_no_veto_changes_nothing(perms):
    assert getattr(perms, "auto_veto", "missing") is None
    out = ac._doc_executor.execute("read_file", {"path": "a.txt"}, perms)
    assert "hello" in out


def test_a_refusal_stops_the_call_before_it_runs(perms, monkeypatch):
    calls = []
    perms.auto_veto = lambda name, args: (calls.append((name, args))
                                          or "blocked: queued for the user")

    def _no(*a, **k):
        raise AssertionError("a vetoed command must not start a process")

    monkeypatch.setattr(subprocess, "Popen", _no)
    monkeypatch.setattr(subprocess, "run", _no)
    out = ac._doc_executor.execute("bash", {"command": "git push"}, perms)
    assert json.loads(out)["error"] == "blocked: queued for the user"
    assert calls == [("bash", {"command": "git push"})]


def test_a_write_is_refused_and_the_file_unchanged(perms, tmp_path):
    perms.auto_veto = lambda name, args: (
        "blocked: measure" if name == "write_file" else None)
    out = ac._doc_executor.execute(
        "write_file", {"path": "a.txt", "content": "changed"}, perms)
    assert "blocked: measure" in out
    assert (tmp_path / "a.txt").read_text() == "hello\n"
    # A read passes the same veto.
    assert "hello" in ac._doc_executor.execute("read_file",
                                               {"path": "a.txt"}, perms)


def test_a_veto_that_raises_refuses(perms):
    def boom(name, args):
        raise RuntimeError("x")

    perms.auto_veto = boom
    out = ac._doc_executor.execute("read_file", {"path": "a.txt"}, perms)
    assert "error" in json.loads(out)


def test_mcp_calls_ask_it_too(perms):
    perms.auto_veto = lambda name, args: "blocked: sends"
    out = ac._mcp_call_refusal("mcp__slack__send_message", {"text": "x"},
                               perms)
    assert out is not None and "blocked: sends" in out


def test_a_subagent_copy_carries_the_veto(perms):
    veto = lambda name, args: None  # noqa: E731
    perms.auto_veto = veto
    child = dataclasses.replace(perms)
    assert child.auto_veto is veto
