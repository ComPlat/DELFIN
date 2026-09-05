"""The marker that says a sub-agent owes its parent a report survives.

The defect, measured before this existed: ``_note_pending_report`` wrote
the marker with ``Path.write_text``, which opens with ``"w"`` and so
truncates before it writes — and the very first thing it does is run the
reaper over that same directory. The reaper read an in-flight marker,
failed to parse it, substituted an empty record, and an empty record has
owner pid 0, which counts as ALIVE, with ``started_at`` 0.0, which makes
its age about 1.8e9 seconds — past every TTL. So it unlinked it. The
sibling reaper over the running registry skips an unreadable entry;
only this one deleted.

Twenty-four concurrent reservations lost two markers per trial, four
over five trials. A lost marker means the parent is never told the
background sub-agent finished: the report is not late, it is gone.

What is pinned here: the write is atomic, an unreadable record is skipped
rather than treated as expired, and the ordinary cases the reaper exists
for — a marker whose process is gone, one past its TTL — still go.
"""

from __future__ import annotations

import json
import os
import threading

import pytest

from delfin.agent import subagents as sa


@pytest.fixture(autouse=True)
def _iso(monkeypatch, tmp_path):
    monkeypatch.setattr(sa, "_PENDING_DIR", tmp_path / "pending")
    monkeypatch.setattr(sa, "_RUNNING_DIR", tmp_path / "running")
    monkeypatch.setattr(sa, "_SESSIONS_DIR", tmp_path / "sessions")


# ---------------------------------------------------------------------------
# The mechanism, named
# ---------------------------------------------------------------------------

def test_an_unreadable_record_is_skipped_rather_than_reaped():
    """A marker being written this instant is unreadable for as long as
    the write takes. Reading that as 'expired' deletes a live one."""
    sa._PENDING_DIR.mkdir(parents=True, exist_ok=True)
    half_written = sa._pending_path("sa-torn")
    half_written.write_text('{"sa_id": "sa-torn", "start', encoding="utf-8")

    assert sa.reap_pending_reports() == []
    assert half_written.exists()


def test_a_marker_is_never_visible_half_written():
    """The write goes through a temp file, so the reaper's own glob can
    only ever see a complete record."""
    sa._note_pending_report("sa-atomic", subagent_type="explore",
                            description="read the docs")
    path = sa._pending_path("sa-atomic")
    rec = json.loads(path.read_text(encoding="utf-8"))
    assert rec["sa_id"] == "sa-atomic"
    assert rec["owner_pid"] == os.getpid()
    assert list(sa._PENDING_DIR.glob("*.json")) == [path]


# ---------------------------------------------------------------------------
# The measurement
# ---------------------------------------------------------------------------

def test_concurrent_reservations_lose_no_markers():
    """The reproduction: 24 at once, five times over. Every reserved id
    must still be owed a report at the end of it."""
    lost_total = 0
    for trial in range(5):
        ids = [f"sa-{trial}-{i:03d}" for i in range(24)]
        threads = [
            threading.Thread(target=sa._note_pending_report, args=(sa_id,),
                             kwargs={"subagent_type": "explore",
                                     "description": "d"})
            for sa_id in ids
        ]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        present = {p.stem for p in sa._PENDING_DIR.glob("*.json")}
        lost = set(ids) - present
        lost_total += len(lost)
        for sa_id in ids:
            sa._claim_pending_report(sa_id)
    assert lost_total == 0


def test_a_reservation_still_owes_the_parent_a_report():
    """The contract the marker exists for, end to end: reserving records
    the debt and claiming it takes it exactly once."""
    sa.reserve_running("sa-owed", subagent_type="explore", description="d")
    assert sa._pending_path("sa-owed").exists()
    assert sa._claim_pending_report("sa-owed") is True
    assert sa._claim_pending_report("sa-owed") is False


# ---------------------------------------------------------------------------
# What the reaper is still for
# ---------------------------------------------------------------------------

def test_a_marker_whose_process_is_gone_is_still_reaped():
    sa._PENDING_DIR.mkdir(parents=True, exist_ok=True)
    sa._atomic_write_json(sa._pending_path("sa-dead"), {
        "sa_id": "sa-dead", "started_at": 0.0,
        "owner_pid": 2 ** 22 - 1, "owner_start_ticks": 12345,
    })
    assert sa.reap_pending_reports() == ["sa-dead"]
    assert not sa._pending_path("sa-dead").exists()


def test_a_marker_past_its_ttl_is_still_reaped():
    import time

    sa._PENDING_DIR.mkdir(parents=True, exist_ok=True)
    sa._atomic_write_json(sa._pending_path("sa-old"), {
        "sa_id": "sa-old",
        "started_at": time.time() - sa._PENDING_TTL_S - 60,
        **sa._owner_stamp(),
    })
    assert sa.reap_pending_reports() == ["sa-old"]


def test_a_fresh_marker_of_this_process_is_left_alone():
    sa._note_pending_report("sa-live")
    assert sa.reap_pending_reports() == []
    assert sa._pending_path("sa-live").exists()


def test_the_drain_does_not_claim_a_record_it_could_not_read():
    """Same substitution, one step further on: the drain would claim the
    marker — the claim IS the unlink — on a record it never read."""
    sa._PENDING_DIR.mkdir(parents=True, exist_ok=True)
    torn = sa._pending_path("sa-torn2")
    torn.write_text("{not json", encoding="utf-8")
    assert sa.drain_finished_subagents() == []
    assert torn.exists()


# ---------------------------------------------------------------------------
# ... and it goes to the session that is owed it
# ---------------------------------------------------------------------------

def _marker(sa_id: str, stamp: dict) -> None:
    sa._PENDING_DIR.mkdir(parents=True, exist_ok=True)
    (sa._PENDING_DIR / f"{sa_id}.json").write_text(json.dumps({
        "sa_id": sa_id, "type": "explore", "description": "d",
        "started_at": 1.0, **stamp,
    }))


@pytest.fixture
def _terminal(monkeypatch):
    """Every marker looks like a run that has ended, so the only thing
    deciding the outcome is ownership."""
    monkeypatch.setattr(sa, "read_running",
                        lambda include_dead=False: ({} if not include_dead
                                                    else {"any": {}}))
    monkeypatch.setattr(sa, "load_subagent_session", lambda sa_id: {"sa_id": sa_id})
    monkeypatch.setattr(sa, "get_subagent_result",
                        lambda sa_id: {"sa_id": sa_id, "status": "finished",
                                       "final_text": f"report for {sa_id}"})


def test_a_live_sessions_report_is_not_claimed_by_another_session(_terminal):
    """The docstring said the owner stamp keeps two concurrent sessions on
    one machine from draining each other's reports. It did not: the drain
    only asked ``_entry_owner_alive`` — is the writer still around — which
    a DIFFERENT live session passes.

    Measured before this: a marker stamped with another live pid was
    drained here and its report handed to a parent that never started the
    delegate, marker unlinked. The exactly-once contract held perfectly
    and delivered to the wrong session, while the right one waited for a
    report that could no longer arrive."""
    theirs = dict(sa._owner_stamp())
    theirs["owner_pid"] = 1               # pid 1 is alive on any POSIX host
    theirs.pop("owner_start_ticks", None)
    _marker("sa-theirs", theirs)

    assert sa.drain_finished_subagents(limit=10) == []
    # Still there, for the session that is actually owed it.
    assert (sa._PENDING_DIR / "sa-theirs.json").exists()


def test_our_own_report_still_arrives(_terminal):
    """The half that must not be broken by the half above."""
    _marker("sa-ours", dict(sa._owner_stamp()))
    got = sa.drain_finished_subagents(limit=10)
    assert [g["sa_id"] for g in got] == ["sa-ours"]
    assert not (sa._PENDING_DIR / "sa-ours.json").exists()


def test_a_recycled_pid_is_not_the_same_process():
    """Same pid number, different start time. After a reboot or a pid wrap
    the number is reissued, and matching on the number alone would hand a
    stranger's report to whoever inherited their pid.

    Asserted against ``_entry_owned_by_us`` DIRECTLY, not through the
    drain. Driving it through the drain proved nothing: the liveness check
    that runs first already rejects this record, so that version of the
    test stayed green with the ownership check deleted."""
    recycled = dict(sa._owner_stamp())
    recycled["owner_start_ticks"] = 999_999_999
    assert sa._entry_owned_by_us(recycled) is False
    # ... while the unmodified stamp of this very process is us.
    assert sa._entry_owned_by_us(dict(sa._owner_stamp())) is True


def test_a_marker_written_before_the_owner_field_is_still_delivered(_terminal):
    """Backward compatibility, on purpose: nobody can match a marker with
    no owner recorded, and leaving it to the TTL reaper would strand a
    report that used to arrive."""
    _marker("sa-legacy", {})
    got = sa.drain_finished_subagents(limit=10)
    assert [g["sa_id"] for g in got] == ["sa-legacy"]
