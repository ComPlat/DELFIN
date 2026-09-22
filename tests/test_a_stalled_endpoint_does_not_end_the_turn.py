"""A transient endpoint error must not end the turn for good.

Measured 2026-09-21, a run of six sessions: the largest share of the
operator's work was not granting things but FINDING standstill. A turn
died on::

    Error code: 400 - Open WebUI: Server Connection Error

at the KIT endpoint -- a temporary one. The turn ended, the session then
sat at the prompt for 24 minutes, and nobody noticed. Similarly "Stopped
while waiting for the endpoint".

The retry ladder exists (api_client._is_transient_api_error plus
_STREAM_RETRY_MAX in stream_message), but this 400 matched no rule: the
class name is BadRequestError (the openai SDK wraps any 400 in it),
the status 400 is not a transient status, and none of the
proxy-hiccup markers matched. So the one failure a shared gateway
produces when it cannot reach the model behind it -- which is
exactly the temporary kind -- was treated as a deterministic client
error.

The rule added here is narrow on purpose, in the spirit of the
"Extra data" and infra-exhaustion rules beside it: a chat request
cannot be malformed in a way that makes the gateway say its own
SERVER CONNECTION failed. Only retried when the message says so;
a real bad request (model not found, context length, bad params)
still fails on the first attempt.

And when the ladder does give up, somebody must hear it: a session
that stood still for 24 minutes with the operator watching six
termals is the finding. The engine records the turn as an error; the
terminal now leaves a note for the operator session when one is open,
so the standstill is announced on the channel the written
instructions prescribe.

What is judged here, and the instrument:

  the 400 is transient        _is_transient_api_error on the exact
                             message measured, yes; on a genuine
                             bad request (401-shaped, model not
                             found), no

  give-up reaches the operat- a fake client whose stream dies with
  or (engine)                the same error, retried to the end;
                             the engine emits an event the UI can
                             hook, and writes the operator inbox
                             when a session is open under that key
"""

from __future__ import annotations

import pytest

MEASURED = ("Error code: 400 - {'detail': 'Open WebUI: Server "
            "Connection Error'}")


def test_a_401_still_fails_at_once():
    """A real client error is never retried. The 400-marker rules
    themselves are red controls in
    test_a_gateway_that_lost_the_model_is_not_a_bad_request.py --
    api_client.py is owned by another session this run."""
    from delfin.agent.api_client import _is_transient_api_error

    class _Unauthorised(Exception):
        status_code = 401

    assert not _is_transient_api_error(
        _Unauthorised("Error code: 401 - invalid key"))


# -- the engine does not die silently -------------------------------------
#
# The harness is the crashed-turn test's: a real AgentEngine over a
# MagicMock client whose stream raises. The engine's contract today is
# that the exception propagates to the UI; what was missing on
# 2026-09-21 is that NOTHING was said to the operator when nobody sat
# at the terminal -- so the two assertions below are about the note,
# and about the turn log still carrying the classified error.


def _boom(system, messages, **kw):
    raise RuntimeError(MEASURED)
    yield  # pragma: no cover -- generator shape only


class _GatewayTimeout(Exception):
    """Already classified transient today (5xx) — the shape the engine
    note is judged with, independent of the 400-marker rule above."""


def _boom_transient(system, messages, **kw):
    raise _GatewayTimeout("Error code: 503 - Service Unavailable")
    yield  # pragma: no cover -- generator shape only


def _engine(agent_tree, home, stream=None):
    from unittest.mock import MagicMock, patch
    from delfin.agent.engine import AgentEngine

    client = MagicMock()
    client.model = "kit.test-model"
    client.stream_message = MagicMock(side_effect=stream or _boom)
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = AgentEngine(repo_dir=agent_tree, backend="api",
                          provider="kit", mode="quick",
                          pack_dir=agent_tree)
    eng.client = client
    return eng


@pytest.fixture
def home(monkeypatch, tmp_path):
    from pathlib import Path as _P
    monkeypatch.setattr(_P, "home", lambda: tmp_path)
    return tmp_path


@pytest.fixture
def agent_tree(tmp_path):
    import textwrap
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# solo mode")
    manifest = textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """)
    (lite_dir / "manifest.yaml").write_text(manifest)
    return tmp_path


def test_the_engine_records_and_classifies_the_given_up_turn(
        agent_tree, home):
    eng = _engine(agent_tree, home)
    with pytest.raises(RuntimeError):
        eng.stream_response("go on")
    from delfin.agent import turn_metrics as tm
    entries = tm.read(eng.trace_session())
    assert entries, "the turn should have been recorded"
    assert "Server Connection Error" in entries[-1]["error"]


def test_a_given_up_turn_leaves_a_note_for_the_operator(
        agent_tree, home, monkeypatch):
    from delfin.agent import session_messages as msgs
    from delfin.agent import session_presence as pres

    monkeypatch.setattr(pres, "_DIR", home / "presence")
    monkeypatch.setattr(msgs, "_DIR", home / "inbox")
    pres.announce("operator", session_id="op", title="operator",
                  workspace=str(home))
    eng = _engine(agent_tree, home, stream=_boom_transient)
    with pytest.raises(_GatewayTimeout):
        eng.stream_response("go on")
    # One message in the operator's inbox, saying the turn gave up.
    note = msgs.take("operator")
    assert len(note) == 1
    assert "503" in note[0]["text"]


def test_no_note_when_no_operator_is_open(agent_tree, home, monkeypatch):
    """Silence was the finding; but a note to nobody is not a report
    either. Without an open operator session nothing is written, and
    the turn still raises exactly as before."""
    from delfin.agent import session_messages as msgs
    from delfin.agent import session_presence as pres

    monkeypatch.setattr(pres, "_DIR", home / "presence")
    monkeypatch.setattr(msgs, "_DIR", home / "inbox")
    eng = _engine(agent_tree, home, stream=_boom_transient)
    with pytest.raises(_GatewayTimeout):
        eng.stream_response("go on")
    assert msgs.take("operator") == []


def test_the_note_names_the_session_that_stalled(
        agent_tree, home, monkeypatch):
    """A note with no sender answers "some session is stuck" -- the
    operator's next question is WHICH one, over six panes. The sender
    is the presence key on the permissions object, and the chain that
    sets it is: the engine mints the id, the id (or the -n name)
    decides the key, the key is written to the permissions LAST. The
    turn below reads it back through the real chain."""
    from delfin.agent import cli as C
    from delfin.agent import session_messages as msgs
    from delfin.agent import session_presence as pres

    monkeypatch.setattr(pres, "_DIR", home / "presence")
    monkeypatch.setattr(msgs, "_DIR", home / "inbox")
    pres.announce("operator", session_id="op", title="operator",
                  workspace=str(home))
    eng = _engine(agent_tree, home, stream=_boom_transient)
    # What cmd_chat's _open_session does to a fresh session, minus the
    # argument parser: mint/keep the id, then write the key.
    import types
    args = types.SimpleNamespace(session_name="runde2-s8",
                                 session="", new_session=True,
                                 fork_session=False)
    monkeypatch.setattr(
        C, "_claim_session", lambda _sid: True)
    assert C._open_session(eng, args, home) is True
    assert getattr(eng.kit_permissions, "presence_key", "") == "runde2-s8"
    with pytest.raises(_GatewayTimeout):
        eng.stream_response("go on")
    note = msgs.take("operator")
    assert len(note) == 1
    assert note[0]["from"] == "runde2-s8"
