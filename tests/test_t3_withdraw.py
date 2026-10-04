"""T3 phase-2 (withdraw tombstone) — RED control.

QS ruling (nacht-s26): a gracefully closed (withdrawn) session must be
deliverable during the close/restart gap, like a crashed one, otherwise the
end-user path refuses 'no other open session' for any message sent between
clean-exit withdraw() and re-announce (repl.py:2900 and
dashboard/agent_sessions.py:626 both call P.withdraw() on clean exit — the
gap is real).

Contract. session_presence.withdraw(key) writes a tombstone at
P._path(key) (the SAME path the presence record used) instead of unlinking
it:

    {"key": key, "tombstone": true, "withdrawn_at": <unix time.time()>,
     "pid": 0, "host": <hostname>}

pid 0 makes it inert to open_sessions()/_alive (a dead record): it never
counts as an open session, and _reap may remove it after _REAP_AFTER_S. The
read side (session_messages._known_key) treats this tombstone as marking an
existing session only while fresh: now - withdrawn_at < _STALE_S. Once it
ages past _STALE_S the address reverts to refuseable unless a real record or
inbox exists -- a typo is never made permanently deliverable, and the
operator stays always deliverable.

The test writes the tombstone exactly as the (operator-built) withdraw()
will, so it verifies the read side end-to-end and pins the write contract.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import session_messages as M
from delfin.agent import session_presence as P


@pytest.fixture(autouse=True)
def _dirs(tmp_path, monkeypatch):
    monkeypatch.setattr(M, "_DIR", tmp_path / "inbox")
    monkeypatch.setattr(P, "_DIR", tmp_path / "presence")
    P._last_written.clear()
    P._git_cache.clear()


def _write_tombstone(key: str, withdrawn_at: float | None = None) -> None:
    """Write the withdrawal tombstone the way operator-built withdraw() will."""
    from delfin.agent.state_paths import ensure_dir, write_text_atomic
    ensure_dir(P._DIR)
    rec = {
        "key": key,
        "tombstone": True,
        "withdrawn_at": withdrawn_at if withdrawn_at is not None else time.time(),
        "pid": 0,
        "host": "testhost",
    }
    write_text_atomic(P._path(key), json.dumps(rec, sort_keys=True))


def test_withdrawn_session_is_known_while_tombstone_fresh():
    """The core finding: a gracefully closed session is deliverable during the
    close/restart gap. Tombstone fresh -> _known_key True -> deliverable."""
    _write_tombstone("nacht-s12", withdrawn_at=time.time() - 1)
    assert M._known_key("nacht-s12") is True
    assert M.deliverable("nacht-s12") is True


def test_withdrawn_session_ages_back_to_refuseable():
    """Once the tombstone is older than _STALE_S the address is unknown again
    (no real record, no inbox) -- a typo is never permanently deliverable."""
    stale = time.time() - P._STALE_S - 5
    _write_tombstone("nacht-s13", withdrawn_at=stale)
    assert M._known_key("nacht-s13") is False
    assert M.deliverable("nacht-s13") is False


def test_never_announced_never_withdrawn_typo_stays_refused():
    """A gibberish address with no record and no tombstone is refused."""
    assert M._known_key("zz-typo") is False
    assert M.deliverable("zz-typo") is False


def test_operator_stays_deliverable_after_withdraw():
    """The operator is a reserved mailbox, always deliverable."""
    _write_tombstone("operator", withdrawn_at=time.time() - 1)
    assert M.deliverable("operator") is True


def test_live_announced_session_still_deliverable():
    """An actively announced session remains deliverable (guard the base case)."""
    P.announce("nacht-s14", session_id="s", title="t")
    assert M._known_key("nacht-s14") is True
    assert M.deliverable("nacht-s14") is True


def test_real_withdraw_writes_a_fresh_tombstone():
    """END-TO-END red control for the operator-built withdraw(): a clean exit
    must write a tombstone (not unlink), so the read side can bridge the
    close/restart gap. Fails on the current unlink-withdraw().

    This is the write side of the contract the read-side tests above assume.
    """
    P.announce("nacht-s15", session_id="s", title="t")
    P.withdraw("nacht-s15")
    rec = json.loads(P._path("nacht-s15").read_text(encoding="utf-8"))
    assert rec.get("tombstone") is True
    assert float(rec.get("withdrawn_at") or 0) <= time.time()
    assert M.deliverable("nacht-s15") is True
