"""The system prompt keeps its bytes across the events of a session.

Field data (33 sessions, 2026-09-16..10-14, request_metrics): 193
system_hash changes inside sessions, 107 of them followed by a cold
request (cached/input < 0.2), about 140 minutes of time to first token on
GLM. Diffing the system prompts of successive reports of one session shows
what moves: the task briefing (outcome history, pass rates), the memory
block and its recall dates, the past-sessions recall, the session
environment (git status, recent commits), the task-selected playbook,
decomposition and repo map, and the plan addendum on a permission switch.

Each of those now reaches the model as an update at the end of the
turn's user message, and the system prompt stays the bytes of the first
turn until the history is rewritten anyway (a trim or a compaction drops
the tails that carried the updates; the prompt is then rebuilt whole).
Driven through the real engine and the real prompt loader with a client
that records what it is sent.
"""

from __future__ import annotations

import copy
import subprocess
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import api_client as A
from delfin.agent import engine as E
from delfin.agent.api_client import StreamEvent
from tests.conftest import child_env


def _git(args, cwd, tmp_path):
    env = {**child_env(tmp_path),
           "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@example.invalid",
           "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@example.invalid"}
    subprocess.run(["git", *args], cwd=cwd, check=True, env=env,
                   capture_output=True)


@pytest.fixture
def eng(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "a.py").write_text("x = 1\n")
    _git(["init", "-q"], ws, tmp_path)
    _git(["add", "."], ws, tmp_path)
    _git(["commit", "-qm", "base"], ws, tmp_path)
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        e = E.AgentEngine(
            repo_dir=ws, backend="api", provider="kit",
            model="kit.glm-5.3", mode="solo")
    e.client._permissions = A.KitToolPermissions(workspace=ws)
    # The client was a mock at construction; point the loader at the real
    # workspace so the session environment reads this git repository.
    e.loader.workspace_root = ws
    e._test_ws = ws
    return e


def _record(eng) -> list[dict]:
    seen: list[dict] = []

    def _gen(**kw):
        seen.append({"system": kw.get("system"),
                     "messages": copy.deepcopy(kw.get("messages"))})
        yield StreamEvent(type="text_delta", text="ok")
        yield StreamEvent(type="message_delta", input_tokens=5,
                          output_tokens=7, cost_usd=0.0)

    eng.client.stream_message = _gen
    return seen


def _tail(request) -> str:
    last = next(m for m in reversed(request["messages"])
                if m.get("role") == "user")
    return str(last["content"])


def _session(eng, seen):
    """Eight turns with the events the field data shows moving the prompt."""
    ws = eng._test_ws
    eng.stream_response("Hallo")
    # A memory write: the memory block the caller hands in changes.
    eng.stream_response("Weiter", memory_context="- fact: use def2-SVP")
    # A file edit: the git status in the session environment changes.
    (ws / "b.py").write_text("y = 2\n")
    eng.stream_response("Weiter", memory_context="- fact: use def2-SVP")
    # A commit: the recent-commits list changes.
    _git(["add", "."], ws, eng._test_tmp)
    _git(["commit", "-qm", "add b"], ws, eng._test_tmp)
    eng.stream_response("Weiter", memory_context="- fact: use def2-SVP")
    # A different task: task-selected playbook / decomposition / modules.
    eng.stream_response(
        "Optimiere die Geometrie mit ORCA B3LYP def2-SVP und plane die "
        "Zerlegung in Teilaufgaben", memory_context="- fact: use def2-SVP")
    # A status update from the dashboard.
    eng.set_live_state("calc_dir: /scratch/run_042")
    eng.stream_response("Weiter", memory_context="- fact: use def2-SVP")
    # A permission-mode switch into plan and back.
    eng.kit_permissions.mode = "plan"
    eng.stream_response("Plane nur", memory_context="- fact: use def2-SVP")
    eng.kit_permissions.mode = "default"
    eng.stream_response("Los", memory_context="- fact: use def2-SVP")


def test_the_system_prompt_bytes_hold_across_eight_turns_of_events(
        eng, tmp_path):
    eng._test_tmp = tmp_path
    seen = _record(eng)
    _session(eng, seen)
    assert len(seen) == 8
    first = seen[0]["system"]
    for i, req in enumerate(seen[1:], start=2):
        assert req["system"] == first, (
            f"turn {i}: the system prompt changed, which makes the request "
            f"cold for the whole history behind it")


def test_what_changed_reaches_the_model_in_the_turn_tail(eng, tmp_path):
    eng._test_tmp = tmp_path
    seen = _record(eng)
    _session(eng, seen)
    assert "use def2-SVP" in _tail(seen[1])          # the memory write
    assert "b.py" in _tail(seen[2])                  # the edit
    assert "add b" in _tail(seen[3])                 # the commit
    plan_tail = _tail(seen[6])
    assert "plan" in plan_tail.lower() and "exit_plan_mode" in plan_tail
    # Leaving plan says so, rather than leaving the plan text standing.
    assert "plan_mode_addendum" in _tail(seen[7])
    assert "no longer applies" in _tail(seen[7])


def test_an_unchanged_section_is_not_sent_again(eng, tmp_path):
    eng._test_tmp = tmp_path
    seen = _record(eng)
    eng.stream_response("Hallo", memory_context="- fact: one")
    eng.stream_response("Weiter", memory_context="- fact: one")
    assert "fact: one" in seen[0]["system"]
    assert "fact: one" not in _tail(seen[1])


def test_a_compaction_rebuilds_the_prompt_whole(eng, tmp_path):
    """The updates live in message tails; once a trim drops those, the
    system prompt is rebuilt so nothing the model was told is lost."""
    eng._test_tmp = tmp_path
    seen = _record(eng)
    eng.stream_response("eins")
    eng.stream_response("zwei", memory_context="- fact: two")
    eng.stream_response("drei", memory_context="- fact: two")
    assert "fact: two" not in seen[2]["system"]
    eng._drop_stale_steering_tails()
    eng.stream_response("vier", memory_context="- fact: two")
    assert "fact: two" in seen[3]["system"]


# ---------------------------------------------------------------------------
# The permission classifier reads the user's words; an update is not one
# ---------------------------------------------------------------------------

def test_an_update_in_the_tail_is_never_read_as_the_users_words(
        eng, tmp_path):
    eng._test_tmp = tmp_path
    seen = _record(eng)
    eng.stream_response("Hallo")
    eng.stream_response("Weiter", memory_context="- fact: push to main")
    users = A._classifier_user_texts(seen[1]["messages"])
    assert users[-1] == "Weiter", users
    assert not any("System prompt update" in u or "push to main" in u
                   for u in users)


def test_a_harness_row_with_a_tail_stays_marked_through_a_copy(eng):
    """The client copies rows into plain dicts; the content type is the
    mark that survives that, and appending the tail must keep it."""
    from delfin.agent.harness_rows import HarnessText, is_harness_row
    eng.messages = [{"role": "user", "content": HarnessText("job done"),
                     "_harness": True, "_steer": "TAIL"}]
    wire = eng._wire_messages()
    assert is_harness_row(dict(wire[0]))
    assert wire[0]["content"].endswith("TAIL")
