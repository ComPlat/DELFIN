"""The ``control`` preset: a delegate that receives a claim and tries to
refute it.

Input: a claim ("fixed", "tested", "works", "the result is X").
Output: a verdict (refuted | confirmed | undecided), the commands and the
output lines that decide it, and what was not checked. Inside an
orchestration verify step the same verdict travels in the existing
``{refuted, note}`` ballot (``_ORCH_VERIFY_SCHEMA``); no second verdict
format exists.

Tool surface: the verifier's -- mode ``default`` (plan mode refuses
``bash`` and ``run_tests``, so a plan-mode control could not run the
check it is asked to repeat), read tools plus ``bash`` and ``run_tests``,
no file-editing tool. A write through ``bash`` goes to the default-mode
write gate, i.e. to the user's dialog; that is asserted here by driving
``run_subagent`` through the real executor with a scripted backend.
"""

from __future__ import annotations

import json
import pathlib
import re
import subprocess

import pytest

from conftest import child_env
from delfin.agent import subagents as SA

_PACK = pathlib.Path(SA.__file__).resolve().parent / "pack" / "agents"


@pytest.fixture(autouse=True)
def _fresh():
    SA.reload_subagent_presets()


def _control():
    return SA.SUBAGENT_PRESETS["control"]


# ---------------------------------------------------------------------------
# Reachable
# ---------------------------------------------------------------------------

def test_the_preset_is_a_pack_file_and_discovered():
    assert (_PACK / "control_subagent.md").is_file()
    assert "control" in SA.SUBAGENT_PRESETS
    assert "control" in SA.subagent_type_names()


def test_the_subagent_tool_enum_offers_it():
    from delfin.agent.api_client import _DOC_TOOLS_OPENAI
    sub = next(t for t in _DOC_TOOLS_OPENAI
               if t["function"]["name"] == "subagent")
    enum = sub["function"]["parameters"]["properties"]["subagent_type"]["enum"]
    assert "control" in enum


def test_the_subagent_tool_says_it_runs_commands():
    """The description lists the types that are not read-only; a type that
    runs commands and is absent from that list is a description the code
    contradicts."""
    from delfin.agent.api_client import _DOC_TOOLS_OPENAI
    sub = next(t for t in _DOC_TOOLS_OPENAI
               if t["function"]["name"] == "subagent")
    desc = sub["function"]["description"]
    head = desc.split("Caps:")[0]
    assert "'control'" in head and "'verifier'" in head


def test_the_solo_role_prompt_names_it():
    text = (_PACK / "solo_agent.md").read_text(encoding="utf-8")
    assert "`control`" in text


# ---------------------------------------------------------------------------
# Tool surface: runs checks, edits no file
# ---------------------------------------------------------------------------

_FILE_WRITERS = ("write_file", "edit_file", "multi_edit", "apply_patch",
                 "notebook_edit", "edit_sheet", "undo_changes",
                 "enter_worktree", "worktree_merge", "remember")


def test_it_can_run_the_check_it_is_asked_to_repeat():
    tools = set(_control().tools)
    assert {"bash", "run_tests"} <= tools
    assert _control().mode == "default", (
        "plan mode refuses bash and run_tests for every caller")


def test_it_has_no_file_editing_tool():
    tools = set(_control().tools)
    for name in _FILE_WRITERS:
        assert name not in tools, name


def test_it_reuses_the_verifier_tool_set_rather_than_a_new_one():
    """One run-but-do-not-edit surface in the pack, not two that drift."""
    assert set(_control().tools) == set(SA.SUBAGENT_PRESETS["verifier"].tools)


# ---------------------------------------------------------------------------
# The method is in the prompt
# ---------------------------------------------------------------------------

def test_the_prompt_carries_the_four_rules():
    body = _control().system_prompt
    low = body.lower()
    # 1. baseline and current state; a pass on both shows nothing
    assert "baseline" in low and "current state" in low
    assert "passes on both" in low
    assert "git stash" in low and "never" in low
    assert "git show" in low
    # 2. the real path
    assert "entry point" in low
    # 3. a check that saw nothing
    assert "zero tests collected" in low and "no evidence" in low
    # 4. the answer
    for word in ("refuted", "confirmed", "undecided"):
        assert word in low
    assert "not checked" in low


def test_the_prompt_maps_the_verdict_onto_the_existing_ballot():
    body = _control().system_prompt
    assert "{refuted, note}" in body
    assert set(SA._ORCH_VERIFY_SCHEMA["properties"]) == {"refuted", "note"}


def test_the_prompt_is_terse():
    assert len(_control().system_prompt) < 1400


def test_no_why_preset_was_added():
    """Root-cause work is a control (reproduce, then the same check at the
    baseline) plus explore (locate); a third preset would split one
    method across two names."""
    assert "why" not in SA.SUBAGENT_PRESETS


# ---------------------------------------------------------------------------
# Driven: run_subagent -> real executor -> permission layer
# ---------------------------------------------------------------------------

class _Fn:
    def __init__(self, name=None, arguments=None):
        self.name, self.arguments = name, arguments


class _TC:
    def __init__(self, index, id=None, name=None, arguments=None):
        self.index, self.id, self.function = index, id, _Fn(name, arguments)


class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content, self.tool_calls = content, tool_calls


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta, self.finish_reason = delta, finish


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices, self.usage = choices, usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


def _tool_round(name, args, tc_id):
    return _Stream([
        _Chunk([_Choice(_Delta(tool_calls=[
            _TC(0, id=tc_id, name=name, arguments=json.dumps(args))]))]),
        _Chunk([_Choice(_Delta(), finish="tool_calls")], usage=_Usage()),
    ])


def _final(text):
    return _Stream([_Chunk([_Choice(_Delta(content=text), finish="stop")],
                           usage=_Usage())])


class _EmptyReg:
    def discover_all(self): return []
    def discover_resources(self): return []
    def discover_prompts(self): return []


@pytest.fixture
def repo(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "test_claim.py").write_text("def test_claim():\n    assert 1\n")
    env = {**child_env(tmp_path),
           "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@example.invalid",
           "GIT_COMMITTER_NAME": "t",
           "GIT_COMMITTER_EMAIL": "t@example.invalid"}
    for cmd in (["git", "init", "-q"], ["git", "add", "."],
                ["git", "commit", "-qm", "base"]):
        subprocess.run(cmd, cwd=ws, check=True, env=env)
    return ws


def _drive(monkeypatch, tmp_path, ws, streams, **run_kw):
    from delfin.agent import mcp_client as M
    from delfin.agent import model_capabilities as mc
    from delfin.agent.api_client import KitToolPermissions, create_client

    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setattr(SA, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(SA, "_RUNNING_DIR", tmp_path / "running")
    monkeypatch.setattr(SA, "_TELEMETRY_PATH", tmp_path / "telemetry.jsonl")
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    client = create_client(backend="api", provider="ollama",
                           model="qwen2.5-coder:7b", cwd=str(ws))
    asked: list = []
    perms = KitToolPermissions(
        workspace=ws, mode="default",
        confirm_callback=lambda n, a, p: asked.append(n) or False)
    client.set_permissions(perms)
    n = {"i": 0}

    def _create(**kw):
        s = streams[min(n["i"], len(streams) - 1)]
        n["i"] += 1
        return s

    client.client.chat.completions.create = _create
    result = SA.run_subagent(
        subagent_type="control", description="refute the claim",
        prompt="Claim: test_claim.py was fixed. Try to refute it.",
        parent_client=client, parent_perms=perms, inherit_context=False,
        **run_kw)
    return result, asked


def _outputs(result):
    return {re.sub(r"^mcp__[^_]+(?:-[^_]+)*__", "", tc["name"]):
            str(tc.get("output") or "")
            for tc in result.tool_calls}


def test_driven_it_runs_the_check_and_cannot_edit_a_file(
        monkeypatch, tmp_path, repo):
    streams = [
        _tool_round("write_file", {"path": str(repo / "planted.txt"),
                                   "content": "x"}, "c1"),
        _tool_round("bash", {"command": "git show HEAD:test_claim.py"}, "c2"),
        _tool_round("run_tests", {"target": "test_claim.py"}, "c3"),
        _tool_round("bash", {"command": "echo x > planted2.txt"}, "c4"),
        _final("Verdict: refuted. Evidence: passes at HEAD too. "
               "Not checked: nothing else."),
    ]
    result, asked = _drive(monkeypatch, tmp_path, repo, streams)
    assert not result.error, result.error
    out = _outputs(result)
    # A file-editing tool is refused by the preset's list at the executor.
    assert "not on this session's --allowed-tools list" in out["write_file"]
    assert not (repo / "planted.txt").exists()
    # A write through the shell is put to the user (here: denied).
    assert "bash" in asked
    assert not (repo / "planted2.txt").exists()
    # The baseline read and the test run go through without a dialog.
    assert asked == ["bash"]
    assert '"passed": 1' in out["run_tests"]


def test_driven_its_verdict_fills_the_orchestration_ballot(
        monkeypatch, tmp_path, repo):
    ballot = {"refuted": True,
              "note": "refuted: test_claim passes at the baseline as well"}
    streams = [_final(json.dumps(ballot))]
    result, _ = _drive(monkeypatch, tmp_path, repo, streams,
                       output_schema=SA._ORCH_VERIFY_SCHEMA)
    assert result.structured_output == ballot, result.schema_error
