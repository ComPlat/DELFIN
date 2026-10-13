"""The ``verifier`` preset: a delegate that receives a claim about a change
and tries to refute it.

Input: a claim ("fixed", "tested", "works", "breaks nothing", "the result
is X").
Output: a verdict (refuted | confirmed | undecided), the commands and the
output lines that decide it, and what was not checked. Inside an
orchestration verify step the same verdict travels in the existing
``{refuted, note}`` ballot (``_ORCH_VERIFY_SCHEMA``); no second verdict
format exists.

One preset, not two. A separate ``control`` preset with the verifier's
tool list, mode and job ("gate work before it is reported") was the
alternative and is not built: the two would differ only in whether the
base run happens when the check passes, which is a rule of the method,
not a different surface. The rules (same check at base and now, the real
path, an empty check is no evidence) live in the verifier's prompt.

Tool surface: mode ``default`` (plan mode refuses ``bash`` and
``run_tests``, so a plan-mode verifier could not run the check it is
asked to repeat), read tools plus ``bash`` and ``run_tests``, no
file-editing tool. A write through ``bash`` goes to the default-mode
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


def _verifier():
    return SA.SUBAGENT_PRESETS["verifier"]


# ---------------------------------------------------------------------------
# One preset carries the method
# ---------------------------------------------------------------------------

def test_there_is_no_second_preset_with_the_same_job():
    assert "control" not in SA.SUBAGENT_PRESETS
    assert not (_PACK / "control_subagent.md").exists()


def test_the_solo_role_prompt_says_what_the_verifier_is_for():
    text = (_PACK / "solo_agent.md").read_text(encoding="utf-8")
    row = next(line for line in text.splitlines()
               if line.startswith("| `verifier`"))
    assert "claim" in row and "base" in row


# ---------------------------------------------------------------------------
# Tool surface: runs checks, edits no file
# ---------------------------------------------------------------------------

_FILE_WRITERS = ("write_file", "edit_file", "multi_edit", "apply_patch",
                 "notebook_edit", "edit_sheet", "undo_changes",
                 "enter_worktree", "worktree_merge", "remember")


def test_it_can_run_the_check_it_is_asked_to_repeat():
    tools = set(_verifier().tools)
    assert {"bash", "run_tests"} <= tools
    assert _verifier().mode == "default", (
        "plan mode refuses bash and run_tests for every caller")


def test_it_has_no_file_editing_tool():
    tools = set(_verifier().tools)
    for name in _FILE_WRITERS:
        assert name not in tools, name


# ---------------------------------------------------------------------------
# The method is in the prompt
# ---------------------------------------------------------------------------

def test_the_prompt_carries_the_rules():
    body = _verifier().system_prompt
    low = body.lower()
    # the covering checks, not a suite named after the change
    assert "cover" in low
    # base and current state; a pass on both shows nothing
    assert "base" in low and "current state" in low
    assert "passes on" in low and "both" in low
    assert "git stash" in low and "never" in low
    assert "git show" in low
    # the real path
    assert "entry point" in low
    # a check that saw nothing
    assert "zero tests collected" in low and "no evidence" in low
    # the failure text is quoted
    assert "message" in low
    # the answer
    for word in ("refuted", "confirmed", "undecided"):
        assert word in low
    assert "not checked" in low
    assert "never edit" in low


def test_the_prompt_maps_the_verdict_onto_the_existing_ballot():
    body = _verifier().system_prompt
    assert "{refuted, note}" in body
    assert set(SA._ORCH_VERIFY_SCHEMA["properties"]) == {"refuted", "note"}


def test_the_prompt_is_terse():
    assert len(_verifier().system_prompt) < 1800



def test_no_why_preset_was_added():
    """Root-cause work is a verifier run (reproduce, then the same check at the
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


def _drive(monkeypatch, tmp_path, ws, streams, *, can_ask=True,
           subagent_type="verifier", **run_kw):
    """``can_ask=False`` is the parent ``delfin agent run``, the benchmark
    and the scheduler hand a delegate: permissions with no confirm
    callback (cli.cmd_chat binds one only for a terminal or --supervised;
    the dashboard always binds one)."""
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
        confirm_callback=((lambda n, a, p: asked.append(n) or False)
                          if can_ask else None))
    client.set_permissions(perms)
    n = {"i": 0}

    def _create(**kw):
        s = streams[min(n["i"], len(streams) - 1)]
        n["i"] += 1
        return s

    client.client.chat.completions.create = _create
    result = SA.run_subagent(
        subagent_type=subagent_type, description="refute the claim",
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


# ---------------------------------------------------------------------------
# With nobody to ask, a shell write by a preset that edits no file is refused
# ---------------------------------------------------------------------------

_SHELL_WRITES = ("echo x > planted1.txt", "touch planted2.txt",
                 "cp test_claim.py planted3.txt")


def _bash_outputs(result):
    return [str(tc.get("output") or "") for tc in result.tool_calls
            if tc["name"].endswith("bash")]


def test_driven_without_a_callback_a_shell_write_does_not_land(
        monkeypatch, tmp_path, repo):
    """Default mode without a confirm callback lets a write inside the
    workspace through: right for an unattended parent that may edit, and
    the opposite of a preset whose tool list holds no file-editing tool.
    Measured on the base of this change: all three forms landed."""
    streams = [_tool_round("bash", {"command": c}, f"w{i}")
               for i, c in enumerate(_SHELL_WRITES)]
    streams += [
        _tool_round("bash", {"command": "git show HEAD:test_claim.py"}, "r1"),
        _tool_round("run_tests", {"target": "test_claim.py"}, "r2"),
        _final("Verdict: undecided."),
    ]
    result, asked = _drive(monkeypatch, tmp_path, repo, streams,
                           can_ask=False)
    assert not result.error, result.error
    for name in ("planted1.txt", "planted2.txt", "planted3.txt"):
        assert not (repo / name).exists(), name
    outs = _bash_outputs(result)
    for out in outs[:3]:
        assert "edits no file" in out, out
    # The check itself still runs.
    assert "def test_claim" in outs[3]
    assert '"passed": 1' in _outputs(result)["run_tests"]


def test_driven_without_a_callback_a_writing_preset_keeps_the_contract(
        monkeypatch, tmp_path, repo):
    """The refusal is keyed on the preset, not on the missing callback:
    an unattended general-purpose delegate still writes in its workspace."""
    streams = [_tool_round("bash", {"command": "echo x > planted1.txt"}, "w"),
               _final("done")]
    result, _ = _drive(monkeypatch, tmp_path, repo, streams, can_ask=False,
                       subagent_type="general-purpose")
    assert not result.error, result.error
    assert (repo / "planted1.txt").exists()
