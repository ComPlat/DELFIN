"""T3 phase-2 adversarial review (reviewer nacht-s17) of builder hash b359852b.

Two findings:

1. GREEN-ON-A-NOOP: the commit adds `session_messages.deliverable()` but wires
   nothing to it. The end-user `session_message` tool (api_client.py:
   _execute_session_message) still resolves the recipient only among
   `open_sessions()` (api_client.py:17620) and refuses every other address
   (api_client.py:17654-17668). A message to a known-but-closed session is
   still refused, and - crucially - a message to `operator` (a reserved
   mailbox, never a presence record in the current code) is still refused
   when no heartbeat has faked an open `operator` record. Repro through the
   real path: with only the sender present, `to="operator"` returns the
   refusal.

2. CASE MISMATCH: `deliverable("OPERATOR")` is True (that is the point of the
   case-insensitive reserved check, session_messages.py:66-67) but
   `_inbox(key)` (session_messages.py:27-29) preserves case, so
   `send("OPERATOR")` writes `OPERATOR.jsonl` while any reader using
   `take("operator")` reads `operator.jsonl` -- the message is queued in an
   inbox nobody reads. Deliverable says "yes" yet the delivery target file
   differs.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import session_messages as M
from delfin.agent import session_presence as P


def _inbox_path(key):
    """Mirror of session_messages._inbox for ageing files in tests."""
    return M._inbox(key)


@pytest.fixture(autouse=True)
def _dirs(tmp_path, monkeypatch):
    monkeypatch.setattr(M, "_DIR", tmp_path / "inbox")
    monkeypatch.setattr(P, "_DIR", tmp_path / "presence")
    P._last_written.clear()
    P._git_cache.clear()


# ---- Finding 2: case mismatch between the deliverable gate and the mailbox.

def test_operator_case_variant_queues_outside_the_read_mailbox(tmp_path):
    """deliverable('OPERATOR') is True, but the message lands in OPERATOR.jsonl
    which take('operator') never reads."""
    # The reserved check is case-insensitive (that is its stated intent).
    assert M.deliverable("OPERATOR") is True
    assert M.deliverable("operator") is True
    # A sender addresses the operator with an uppercase variant.
    M.send("OPERATOR", "hello operator", from_key="reviewer")
    # A reader using the canonical lower-case mailbox must still see it: the
    # reserved operator mailbox is ONE mailbox, not one per case.
    taken = M.take("operator")
    texts = [m.get("text") or "" for m in taken]
    assert any("hello operator" in t for t in texts), (
        "message addressed to OPERATOR must be readable via take('operator'); "
        "the reserved operator mailbox is ONE mailbox, not one per case. "
        "Got %r" % texts
    )


def test_known_session_case_variant_not_deliverable(tmp_path):
    """A session announced as `nacht-s12` is reachable; its case variant
    `NACHT-S12` names a DIFFERENT (empty, never-announced) mailbox. That is
    acceptable; what must not happen is the reserved case being treated as
    case-insensitive while the mailbox is not."""
    # Announce a known session.
    P.announce("nacht-s12", title="t", workspace=str(tmp_path))
    assert M.deliverable("nacht-s12") is True
    assert M.deliverable("NACHT-S12") is False


# ---- Finding 1: deliverable() has no caller; the end-user path still refuses.

def _call_session_message(to, text, me="reviewer", tmp_path=None):
    """Drive the real _execute_session_message path with a minimal executor."""
    from delfin.agent import api_client as AC

    perms = AC.KitToolPermissions(workspace=tmp_path or __import__(
        "pathlib").Path("."))
    perms.presence_key = me
    executor = AC._DocToolExecutor()
    return executor._execute_session_message(
        {"to": to, "message": text}, perms)


def test_operator_is_STILL_refused_by_the_end_user_path_after_phase2(tmp_path):
    """The wave-11/12 user-visible bug: session_message(to='operator') is
    refused unless an open record exists. Phase 2 was supposed to make the
    operator a reserved, always-deliverable mailbox. Neither the merged tree
    nor any patch wires deliverable() into _execute_session_message, so this
    repro must FAIL the phase-2 promise."""
    out = json.loads(_call_session_message("operator", "hi", tmp_path=tmp_path))
    assert "error" not in out, (
        "session_message(to='operator') must NOT refuse now that 'operator' "
        "is a reserved always-deliverable mailbox; got: %r" % out.get("error")
    )
    assert out.get("status") == "sent"


@pytest.mark.skip(reason=(
    "PARKED follow-up per the T3 wave-scope ruling (QS s26 / Operator): the "
    "withdraw-tombstone write side is NOT part of the three in-scope T3 phases. "
    "A gracefully withdrawn session is intentionally non-deliverable this wave; "
    "the tombstone design (session_presence.withdraw writes a bounded tombstone "
    "instead of unlinking the record, read by session_messages._known_key) is "
    "tracked as .gate/t3_withdraw.patch follow-up. Re-enable when that lands."))
def test_queued_to_known_but_closed_session_is_delivered_on_next_start(tmp_path):
    """A message to a known-but-closed session queues and is delivered on its
    next start (take), instead of being lost. Phase-2 promise.

    Parked: a session that calls withdraw() (clean close) has its presence
    record unlinked (session_presence.py:111), so _known_key sees neither
    record nor inbox and the end-user path still refuses it. Confirmed real by
    s16; ruled OUT of this wave's scope. See the skip reason."""
    P.announce("nacht-s12", title="t", workspace=str(tmp_path))
    P.withdraw("nacht-s12")  # closed: presence withdrawn, not open
    out = json.loads(_call_session_message("nacht-s12", "wake me", "reviewer",
                                           tmp_path))
    # The end-user path must not refuse a known-but-closed session.
    assert "error" not in out, out.get("error")
    # ... and the queued message must be readable on next start.
    M.take("nacht-s12")  # drop the (nonexistent) pre-queue so next asserts read
    _announce_again(P, "nacht-s12", str(tmp_path))
    msgs = M.take("nacht-s12")
    assert any("wake me" in (m.get("text") or "") for m in msgs), (
        "the queued message must be delivered on the session's next start"
    )


def _announce_again(pres, key, ws):
    pres.announce(key, title="t", workspace=ws)


# ---- QS-lead angles (nacht-s26): close-mid-write survival, cold operator,
#      the inbox-exists "known" latch.

def test_operator_never_refused_cold(tmp_path):
    """Even with NO presence records at all (cold), a message to the operator
    must queue, not be refused. This is the phase-2 headline promise."""
    out = json.loads(_call_session_message("operator", "hi", "reviewer",
                                           tmp_path))
    assert "error" not in out, out.get("error")
    assert out.get("status") == "sent", out
    # And it must be readable by the operator on its (future) next start.
    msgs = M.take("operator")
    assert any("hi" in (m.get("text") or "") for m in msgs)


def test_a_stale_inbox_stops_marking_a_typo_deliverable(tmp_path):
    """QS lead (nacht-s26) decided: BOUND the inbox-exists known latch with a
    TTL. `_known_key` (session_messages.py:46) treats an address with an inbox
    file as 'known' while the inbox is fresh; the builder aligned this with
    session_presence._STALE_S (the same window (15 min) a dormant presence
    record is reaped at), via _inbox_known (session_messages.py:80). So a live
    queue keeps its recipient known and an idle typo ages out to refused.
    This test ages beyond _STALE_S and expects the typo to revert to refused.
    Green = the TTL is correctly implemented; red = it is not."""
    from os import utime
    # A never-declared, never-messaged address is a typo and is refused.
    assert M.deliverable("totally-missing-typo") is False
    # One stray delivery creates the inbox file; fresh, it is still deliverable.
    M.send("once-fired", "hello", from_key="reviewer")
    assert M.deliverable("once-fired") is True
    # Age the inbox beyond the implemented TTL bound: the typo must revert to
    # refused (loud) rather than stay silently accept-and-drop forever.
    stale = time.time() - P._STALE_S - 60
    utime(_inbox_path("once-fired"), (stale, stale))
    assert M.deliverable("once-fired") is False, (
        "an inbox older than the TTL bound must not keep a never-real "
        "address deliverable; permanent accept-and-drop is the silent-loss "
        "hazard the QS lead's TTL decision removes"
    )
