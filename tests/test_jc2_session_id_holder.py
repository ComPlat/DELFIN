"""JC2 target 2b: the session-id holder semantics (committed test).

The holder — session_store.current_agent_session_id + the setter — is
the documented per-process fact "which agent session does this process
serve".  Its CONSUMER (the mcp_client spawn injecting
DELFIN_AGENT_SESSION_ID) is a protected-file patch; this file pins the
holder's own contract:

- set / overwrite / clear-to-empty round trip;
- an empty value is the "no agent session" fact (falsy), never an
  invented id;
- the value survives unrelated session_store activity (the holder is
  process state, not per-session state).
"""

from __future__ import annotations

import delfin.agent.session_store as ss


def test_holder_set_overwrite_clear_round_trip():
    assert ss.current_agent_session_id == ""
    ss.note_agent_session_id("sess-a")
    assert ss.current_agent_session_id == "sess-a"
    ss.note_agent_session_id("sess-b")
    assert ss.current_agent_session_id == "sess-b", (
        "a later dispatch in the same process overwrites the holder")
    ss.note_agent_session_id("")
    assert ss.current_agent_session_id == "", (
        "clearing is the documented 'no agent session' fact")


def test_holder_rejects_non_string_input_cleanly():
    ss.note_agent_session_id(None)
    assert ss.current_agent_session_id == ""
    ss.note_agent_session_id(12345)
    assert ss.current_agent_session_id == "12345", (
        "the setter stringifies (the ids arrive as strings or ints)")
    ss.note_agent_session_id("")


def test_holder_is_process_state_not_per_session_state():
    ss.note_agent_session_id("sess-x")
    try:
        saved = ss.save_session({"session_id": "other", "messages": [],
                                 "title": "unrelated"})
        assert saved is not None
    finally:
        pass
    assert ss.current_agent_session_id == "sess-x", (
        "the holder is process state; an unrelated store write leaves it")
    ss.note_agent_session_id("")
