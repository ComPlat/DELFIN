"""A delegate could only report at the end, and nothing could reach it.

A session that delegates work had one channel: the report, after the run.
A delegate that hit a decision it could not take had to finish and be
re-spawned, and a session that learned something the delegate needed had
no way to say so. Sessions have had a mailbox between them
(``session_message``) for a while; delegates had none in either
direction.

``subagent_message`` is that channel, and the direction follows who asks:
a delegate's permissions carry ``subagent_id`` (stamped by
``run_subagent``) and for it the tool speaks upward; for a session it
speaks downward, or with no ``to`` lists the delegates running.

Both ends are drained at a ROUND boundary -- after a round's tool results,
before the next request -- which is where the plan-mode redirect and the
repeated-error abort already put text for the next round. So a message
arrives while the work is still going, and no new place in the loop was
invented for it.

Both directions are capped, and upward is the tighter number. A delegate's
job is to do the work and report; interrupting the session is for a
decision it cannot take. The precedent is measured: of 541 tool calls
across twelve sessions in one afternoon, 71 were ``session_message``, much
of it the same announcement sent twice.

In-process by design, like the cancel registry: a delegate runs inside the
session that spawned it, so that session is the only writer that can reach
it and the only reader of its replies. A file-backed queue would promise
delivery across processes, which nothing here can honour.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import subagents as SA
from delfin.agent.api_client import _DOC_TOOLS_OPENAI, _doc_executor


class _Perms:
    """A session has no subagent_id; a delegate has one."""

    def __init__(self, subagent_id=""):
        if subagent_id:
            self.subagent_id = subagent_id


@pytest.fixture(autouse=True)
def _clean():
    SA._NOTES.clear()
    SA._REPLIES.clear()
    yield
    SA._NOTES.clear()
    SA._REPLIES.clear()


@pytest.fixture
def running(monkeypatch):
    """Two delegates of this session are live."""
    entries = {
        "aa11": {"type": "explore", "description": "find the resolver"},
        "bb22": {"type": "plan", "description": "draft the migration"},
    }
    monkeypatch.setattr(SA, "read_running", lambda **kw: dict(entries))
    return entries


def _call(args, perms):
    return json.loads(_doc_executor._execute_subagent_message(args, perms))


# ---------------------------------------------------------------------------
# Session -> delegate
# ---------------------------------------------------------------------------

def test_a_session_can_write_to_a_running_delegate(running):
    out = _call({"to": "aa11", "message": "the schema changed"}, _Perms())
    assert out == {"status": "sent", "to": "aa11"}
    assert SA.notes_waiting("aa11") == 1


def test_the_delegate_reads_it_and_only_once(running):
    SA.post_note("aa11", "the schema changed")
    assert SA.take_notes("aa11") == ["the schema changed"]
    assert SA.take_notes("aa11") == []


def test_notes_are_delivered_in_the_order_they_were_left(running):
    for text in ("first", "second", "third"):
        SA.post_note("aa11", text)
    assert SA.take_notes("aa11") == ["first", "second", "third"]


def test_no_recipient_lists_the_delegates_and_what_they_owe(running):
    SA.post_note("bb22", "x")
    out = _call({}, _Perms())
    rows = {row["sa_id"]: row for row in out["delegates"]}
    assert set(rows) == {"aa11", "bb22"}
    assert rows["bb22"]["notes_queued"] == 1
    assert rows["aa11"]["type"] == "explore"


def test_writing_to_a_delegate_that_is_not_running_says_where_to_look(
        running):
    out = _call({"to": "zz99", "message": "x"}, _Perms())
    assert "error" in out
    assert "subagent_result" in out["error"], (
        "a finished delegate's answer is in its report; say so")
    assert out["delegates"] == ["aa11", "bb22"]


def test_a_note_needs_something_to_say(running):
    assert "error" in _call({"to": "aa11", "message": "  "}, _Perms())
    assert SA.notes_waiting("aa11") == 0


def test_a_session_cannot_bury_a_delegate_in_notes(running):
    for i in range(SA._MAX_NOTES):
        assert SA.post_note("aa11", f"n{i}") is True
    assert SA.post_note("aa11", "one too many") is False
    out = _call({"to": "aa11", "message": "one too many"}, _Perms())
    assert "maximum queued" in out["error"]
    assert SA.notes_waiting("aa11") == SA._MAX_NOTES


def test_a_long_note_is_cut_not_refused(running):
    SA.post_note("aa11", "x" * (SA._MAX_NOTE_CHARS + 500))
    assert len(SA.take_notes("aa11")[0]) == SA._MAX_NOTE_CHARS


# ---------------------------------------------------------------------------
# Delegate -> session
# ---------------------------------------------------------------------------

def test_a_delegate_can_answer_the_session_running_it(running):
    out = _call({"message": "the file does not exist"}, _Perms("aa11"))
    assert out["status"] == "sent"
    assert out["left"] == SA._MAX_REPLIES - 1
    assert SA.take_replies() == [("aa11", "the file does not exist")]


def test_the_reply_says_which_delegate_sent_it(running):
    _call({"message": "from a"}, _Perms("aa11"))
    _call({"message": "from b"}, _Perms("bb22"))
    assert dict(SA.take_replies()) == {"aa11": "from a", "bb22": "from b"}


def test_a_delegate_cannot_address_a_sibling(running):
    """It has exactly one correspondent. `to` is ignored rather than
    honoured -- a delegate-to-delegate channel is one nobody authorised."""
    out = _call({"to": "bb22", "message": "hello sibling"}, _Perms("aa11"))
    assert out["status"] == "sent"
    assert SA.take_replies() == [("aa11", "hello sibling")]
    assert SA.notes_waiting("bb22") == 0


def test_a_delegate_is_told_the_channel_is_finite(running):
    lefts = [_call({"message": f"m{i}"}, _Perms("aa11"))["left"]
             for i in range(SA._MAX_REPLIES)]
    assert lefts == list(range(SA._MAX_REPLIES - 1, -1, -1))


def test_past_the_cap_it_is_told_to_use_its_report(running):
    for i in range(SA._MAX_REPLIES):
        _call({"message": f"m{i}"}, _Perms("aa11"))
    out = _call({"message": "one more"}, _Perms("aa11"))
    assert "error" in out
    assert "report" in out["error"]
    assert SA.replies_waiting() == SA._MAX_REPLIES


def test_upward_is_the_tighter_cap():
    """A delegate interrupting its session is rarer than the reverse."""
    assert SA._MAX_REPLIES < SA._MAX_NOTES


def test_a_reply_needs_something_to_say(running):
    assert "error" in _call({"message": ""}, _Perms("aa11"))
    assert SA.replies_waiting() == 0


def test_a_blank_id_reads_as_a_session_not_as_a_nameless_delegate(running):
    """Whitespace is not an id. post_reply would refuse it anyway, but the
    direction would already have been chosen wrongly -- and a message
    meant for the parent would have been read as one going down."""
    out = _call({}, _Perms("   "))
    assert {row["sa_id"] for row in out["delegates"]} == {"aa11", "bb22"}


def test_a_reply_from_a_run_with_no_id_is_refused():
    assert "error" in SA.post_reply("", "x")
    assert "error" in SA.post_reply("   ", "x")
    assert SA.replies_waiting() == 0


# ---------------------------------------------------------------------------
# Nothing outlives the run it belonged to
# ---------------------------------------------------------------------------

def test_both_queues_are_dropped_when_a_run_ends(running):
    SA.post_note("aa11", "down")
    _call({"message": "up"}, _Perms("aa11"))
    SA.drop_notes("aa11")
    assert SA.notes_waiting("aa11") == 0
    assert SA.replies_waiting() == 0


def test_a_note_does_not_reach_a_later_delegate_with_the_same_id(running):
    """Ids are reused on a resume; a note from the previous run is not for
    the new one."""
    SA.post_note("aa11", "for the old run")
    SA.drop_notes("aa11")
    assert SA.take_notes("aa11") == []


# ---------------------------------------------------------------------------
# The wiring
# ---------------------------------------------------------------------------

def test_the_tool_is_advertised_and_dispatched():
    names = {t["function"]["name"] for t in _DOC_TOOLS_OPENAI}
    assert "subagent_message" in names
    import inspect
    src = inspect.getsource(_doc_executor.__class__)
    assert 'if name == "subagent_message":' in src


def test_the_schema_says_it_is_capped():
    desc = next(t["function"]["description"] for t in _DOC_TOOLS_OPENAI
                if t["function"]["name"] == "subagent_message")
    assert "Capped" in desc or "capped" in desc


def test_a_delegate_is_given_its_own_id():
    """Without it the tool cannot tell a delegate from a session, and a
    delegate's message would be read as a session writing downward."""
    import inspect
    src = inspect.getsource(SA)
    i_id = src.index("_sa_id = resume_from if prior")
    i_stamp = src.index("sub_perms.subagent_id")
    assert i_id < i_stamp, (
        "the stamp must come after the id exists, or it raises into a bare "
        "except and the field stays quietly unset")


def test_both_directions_are_drained_at_a_round_boundary():
    import inspect
    from delfin.agent import api_client as AC
    src = inspect.getsource(AC)
    i = src.index("_note_src = getattr(self, \"_subagent_notes\", None)")
    arm = src[i:i + 3000]
    assert "take_replies" in arm, "the upward half is not drained"
    assert "_wrap_untrusted" in arm, "text from elsewhere must be fenced"
    # The seam sits with the other end-of-round steering, not in a new place.
    assert "_plan_redirect_sent" in src[i:i + 4000]


def test_a_turn_without_delegates_pays_almost_nothing():
    """The guard is a getattr and an empty dict read; it must not reach
    for anything expensive before knowing there is mail."""
    import inspect
    from delfin.agent import api_client as AC
    src = inspect.getsource(AC)
    i = src.index("_note_src = getattr(self, \"_subagent_notes\", None)")
    arm = src[i:i + 900]
    assert "if callable(_note_src):" in arm


# ---------------------------------------------------------------------------
# The rule is a rule when it is in the prompt
# ---------------------------------------------------------------------------

def test_every_preset_can_reach_its_session():
    """_narrow_allowed_tools only narrows, so a preset that does not list
    the channel cannot use it at all."""
    for name, preset in SA.SUBAGENT_PRESETS.items():
        tools = getattr(preset, "tools", ()) or ()
        if not tools:
            continue                      # unrestricted: it has everything
        assert "subagent_message" in tools, name


def test_the_channel_is_listed_as_a_way_of_reporting():
    """Not as a working tool: it is for a decision the delegate cannot
    take, and grouping it with the report is what gets it into every
    read-only preset at once."""
    assert "subagent_message" in SA._REPORTING_TOOLS


def test_the_rule_reaches_a_delegate_built_from_any_preset(monkeypatch):
    """Built, not grepped. The rule is added to the assembled prompt and
    not to each preset, so a preset written by hand gets it too."""
    import inspect

    src = inspect.getsource(SA.run_subagent)
    i_assemble = src.index("system_prompt = (")
    arm = src[i_assemble:i_assemble + 400]
    assert "_CHANNEL_RULE" in arm, (
        "the rule is not in the assembled prompt, so only the presets that "
        "happen to mention it would carry it")


def test_the_rule_states_the_default_and_the_two_exceptions():
    rule = SA._CHANNEL_RULE
    low = rule.lower()
    assert "report" in low
    assert "decision you cannot take" in low
    assert "capped" in low
    # and that an incoming note carries no authority
    assert "not from the user" in low
