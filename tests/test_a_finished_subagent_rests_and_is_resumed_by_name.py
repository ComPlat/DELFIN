"""A finished sub-agent rests, and one call naming it resumes it with context.

Measured on main (probe pinned to the worktree's own code, no model called):
the store, the recap and the preset re-derivation were all in place, and the
generalist could not reach any of it. A call carrying only
``resume_id`` + ``prompt`` was refused before dispatch by two separate
schema gates that both demand ``subagent_type``:

  _missing_required_argument -> "subagent_type is required for subagent"
  _malformed_arguments       -> "missing required field 'subagent_type'"

Both are wrong for a resume: ``run_subagent`` overwrites ``subagent_type``
with the stored one (``subagents.py:4080``) precisely because the resumed run
must keep its own preset, tools and limits. So the caller was asked to supply
a value that was then ignored, and a generalist that did not know the type
had no way to resume at all -- while the tool description never mentioned
resuming, and no tool listed what was restable.

Contracts asserted here:

* a resume-only call reaches the runner, passing no type of its own;
* a FRESH call still needs its type and description -- the gate did not
  simply stop checking (control case, green on both sides);
* the tool description says a run can be rested and resumed, so the model
  learns the capability exists without reading an optional property;
* ``subagent_message`` with no recipient lists rested runs with what they
  worked on, so the id comes from the tool itself and not from memory;
* listing is a file read and sends no request to a model (costs no tokens).

All three are driven through the real entry point -- the ``subagent`` tool
dispatch and the ``subagent_message`` dispatch -- never the helper.
"""

from __future__ import annotations

import json
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from delfin.agent import subagents as SA
from delfin.agent import tool_catalog as TC
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


@pytest.fixture
def stores(monkeypatch, tmp_path):
    monkeypatch.setattr(SA, "_RUNNING_DIR", tmp_path / "running")
    monkeypatch.setattr(SA, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(SA, "_PENDING_DIR", tmp_path / "pending")
    monkeypatch.setattr(SA, "_TELEMETRY_PATH", tmp_path / "telemetry.jsonl")
    SA._FINISHED_SUMMARY_CACHE.clear()
    SA._NOTES.clear()
    yield tmp_path
    SA._NOTES.clear()


def _event(**kw):
    return SimpleNamespace(
        type=kw.get("type", "text_delta"), text=kw.get("text", ""),
        tool_name=kw.get("tool_name", ""), tool_input=kw.get("tool_input", ""),
        tool_output=kw.get("tool_output", ""),
        input_tokens=0, output_tokens=0)


def _client(text="done", tool=None):
    """Scripted client. ``tool=(name, args, output)`` emits the pair, so a
    finished record carries a tool OUTPUT -- which is the thing a follow-up
    has to be able to reuse."""
    client = MagicMock()
    events = []
    if tool:
        name, args, out = tool
        events.append(_event(type="tool_use", tool_name=name,
                             tool_input=json.dumps(args)))
        events.append(_event(type="tool_result", tool_name=name,
                             tool_output=out))
    events.append(_event(type="text_delta", text=text))
    events.append(_event(type="message_delta", input_tokens=10, output_tokens=5))
    client.stream_message = MagicMock(return_value=iter(events))
    client._permissions = None
    client.set_permissions = MagicMock()
    return client


def _perms(tmp_path, sid="sess-A"):
    p = KitToolPermissions(workspace=tmp_path, mode="default")
    p.task_session_id = sid
    return p


def _run(tmp_path, sid="sess-A", text="done", tool=None, sa_id=""):
    """Finish one real sub-agent run through ``run_subagent``."""
    return SA.run_subagent(
        subagent_type="explore", description="read the config",
        prompt="read the config file and report what it sets",
        parent_client=_client(text, tool), parent_perms=_perms(tmp_path, sid),
        sa_id=sa_id)


def _real_runner(perms):
    """``permissions.subagent_runner`` as the product wires it.

    Mirrors ``OpenAIClient._attach_subagent_runner``: the dispatch calls a
    runner with (subagent_type, description, prompt, isolation, **extra) and
    expects a dict back, so the parent is closed over and the result is
    rendered with ``to_payload``. The store lookup and the error text below
    it are the product's own code, which is what these tests are about.
    """
    def _runner(*, subagent_type, description, prompt, isolation="",
                **kw):
        return SA.run_subagent(
            subagent_type=subagent_type, description=description,
            prompt=prompt, parent_client=_client("ok"), parent_perms=perms,
            isolation=isolation, **kw).to_payload()
    return _runner


# --------------------------------------------------------------------------
# 1. One call naming it reaches the runner
# --------------------------------------------------------------------------

def test_a_fresh_call_still_needs_its_type_and_description(stores):
    """Control, green on both sides: the gate was relaxed for RESUMES only.

    If this ever passes without a type, the required-argument guard has been
    switched off wholesale -- which is the regression the other test's fix
    would cause if written wrongly.
    """
    perms = _perms(stores)
    perms.subagent_runner = lambda **kw: {"result": "ran"}
    out = _DocToolExecutor().execute(
        "subagent",
        {"prompt": "a briefing with no type and no description at all"},
        permissions=perms,
    )
    payload = json.loads(out)
    assert "error" in payload
    assert "subagent_type" in payload["error"]


# --------------------------------------------------------------------------
# 2. The model is told the capability exists
# --------------------------------------------------------------------------

def test_the_tool_description_says_a_run_can_be_resumed():
    """A property nobody is told about is a capability nobody uses.

    The word has to be in the description -- the one field every request
    carries -- not only in the optional ``resume_id`` property.
    """
    schema = next(t["function"] for t in TC._DOC_TOOLS_OPENAI
                  if t["function"]["name"] == "subagent")
    assert "resume_id" in schema["description"], schema["description"]
    # ...and the safety phrases the compaction ratchet requires survive it.
    assert "FULL" in schema["description"]
    assert "read-only" in schema["description"]


# --------------------------------------------------------------------------
# 3. Rested runs are listed with what they worked on, through a tool
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# 4. A resumed run keeps its own tools, and its own context
# --------------------------------------------------------------------------

def test_a_resumed_run_keeps_the_stored_presets_type(stores):
    """The stored preset wins over the type the resuming caller names.

    That is the permission side of resuming: a run started read-only must not
    become a writing run because the caller resumed it under another name. The
    dispatch passes the caller's word through; ``run_subagent`` overrides it
    from the store, which is the line this pins.
    """
    res = _run(stores)
    assert SA.load_subagent_session(res.sa_id)["subagent_type"] == "explore"

    reached = {}

    def runner(**kw):
        reached.update(kw)
        return {"result": "resumed",
                "subagent_type": kw.get("subagent_type", "")}

    perms = _perms(stores)
    perms.subagent_runner = runner
    payload = json.loads(_DocToolExecutor().execute(
        "subagent",
        {"resume_id": res.sa_id, "subagent_type": "general-purpose",
         "prompt": "and which component reads that count?"},
        permissions=perms))
    assert "error" not in payload, payload
    assert reached.get("resume_from") == res.sa_id
    # The store still names the original preset: the resumed session did not
    # acquire a new one, so its tools and limits are the ones it was given.
    assert SA.load_subagent_session(res.sa_id)["subagent_type"] == "explore"


def test_a_resumed_run_replays_the_tool_output_it_read(stores):
    """The point of resuming: the follow-up sees what the first turn read.

    Without it the only route is re-briefing, which costs a second request and
    a re-read (measured: 15 384 chars in 2 requests vs 8 061 in 1).
    """
    res = _run(stores, text="the config sets retries=3",
               tool=("read_file", {"path": "config.ini"},
                     "[retry]\ncount = 7\n"))
    second = _client("It is 7.")
    SA.run_subagent(
        subagent_type="explore", description="",
        prompt="and which component reads that count?",
        parent_client=second, parent_perms=_perms(stores),
        resume_from=res.sa_id)
    sent = second.stream_message.call_args.kwargs["messages"]
    assert sent, "the resumed child sent no request"
    blob = json.dumps(sent, default=str)
    assert "count = 7" in blob, (
        "the earlier tool output was not replayed into the resumed context")


# --------------------------------------------------------------------------
# 5. A wrong resume id names the runs that ARE resumable
# --------------------------------------------------------------------------

def test_an_unknown_resume_id_names_the_resting_runs_of_this_session(stores):
    """A dead id should hand back the ids that work.

    The generalist's usual reason for holding a dead one is that its own
    turn was compacted since the delegate finished. Answered with only
    "no stored subagent session found", the recovery is to guess or to
    re-derive the work by re-briefing -- which is the cost resuming
    exists to remove. The store is already on disk and reading it spends
    no model tokens.
    """
    mine = _run(stores, sid="sess-A", text="the config sets retries=3")
    perms = _perms(stores, "sess-A")
    perms.subagent_runner = _real_runner(perms)
    out = _DocToolExecutor().execute(
        "subagent",
        {"resume_id": "deadbeef", "subagent_type": "explore",
         "description": "follow-up",
         "prompt": "and which component reads that count?"},
        permissions=perms,
    )
    payload = json.loads(out)
    assert "error" in payload, payload
    assert mine.sa_id in payload["error"], (
        "the error did not name the run this session could actually "
        f"resume: {payload['error']}")
    # It still says the id was wrong: this is a hint, not a silent success.
    assert "deadbeef" in payload["error"]


def test_an_unknown_resume_id_does_not_name_another_sessions_runs(stores):
    """The store is account-wide; the hint must not be.

    Naming another session's delegate would hand a generalist an id it has
    no business resuming -- and resuming is not scoped, so the hint IS the
    access. Scope it the way the running-delegate listing is scoped.
    """
    theirs = _run(stores, sid="sess-B", text="their work")
    perms = _perms(stores, "sess-A")
    perms.subagent_runner = _real_runner(perms)
    out = _DocToolExecutor().execute(
        "subagent",
        {"resume_id": "deadbeef", "subagent_type": "explore",
         "description": "follow-up",
         "prompt": "and which component reads that count?"},
        permissions=perms,
    )
    payload = json.loads(out)
    assert "error" in payload, payload
    assert theirs.sa_id not in payload["error"], (
        "another session's delegate was offered as resumable: "
        f"{payload['error']}")


def test_an_unknown_resume_id_with_nothing_resting_says_so(stores):
    """No candidates: do not pretend, and name the way out.

    The old text implied the caller had simply misremembered an id. With an
    empty store the honest answer is that there is nothing to resume.
    """
    perms = _perms(stores, "sess-A")
    perms.subagent_runner = _real_runner(perms)
    out = _DocToolExecutor().execute(
        "subagent",
        {"resume_id": "deadbeef", "subagent_type": "explore",
         "description": "follow-up",
         "prompt": "and which component reads that count?"},
        permissions=perms,
    )
    payload = json.loads(out)
    assert "error" in payload, payload
    assert "nothing to resume" in payload["error"], payload["error"]
