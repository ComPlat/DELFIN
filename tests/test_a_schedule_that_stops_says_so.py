"""The scheduler went quiet in exactly the case worth telling about.

THE SILENT DISABLE. An ordinary run failure emitted an attention event.
But the workspace gate — the check that a scheduled entry's directory
still exists — raised OUTSIDE the block that emitted, so a moved or
deleted workspace killed the schedule with no event and no notification,
leaving only a line in a detached daemon log nobody reads. Same for a
one-shot disabled for going stale. A schedule that stops without saying so
is indistinguishable from one that is simply not due yet, which is how a
recurring check on an eight-hour calculation silently stops happening.

THE NON-ATOMIC SAVE. This was the one agent state file written with a
plain ``write_text``: a reader catching it mid-write sees a truncated or
empty schedule, i.e. entries simply gone.

THE OVERWRITING SAVE. An in-process Scheduler loads once and never
reloads, so its next save flattened every entry another process (the
headless daemon, the CLI, a second dashboard) had created in the meantime.

No agent turns are run here — the fire callback is a stub.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pytest

from delfin.agent import scheduler as S
from delfin.agent import scheduler_daemon as SD


@pytest.fixture
def cron(tmp_path) -> Path:
    return tmp_path / "cron.json"


@pytest.fixture
def emitted(monkeypatch) -> list:
    events: list[tuple] = []
    monkeypatch.setattr(
        "delfin.agent.attention.emit_attention",
        lambda kind, **kw: events.append((kind, kw.get("title", ""),
                                          kw.get("detail", ""))))
    return events


# ---------------------------------------------------------------------------
# Every disable is announced
# ---------------------------------------------------------------------------

def test_a_vanished_workspace_disables_the_schedule_out_loud(cron, emitted):
    """The one case the silent path covered, and the one worth an alarm."""
    sch = S.Scheduler(path=cron)
    ent = sch.schedule_interval(every_seconds=3600, prompt="check the opt",
                                fire_immediately=True,
                                workspace="/nowhere/at/all")

    sch.tick(fire_callback=SD.make_fire_callback(log=lambda _m: None))

    assert sch._entries[ent.id].disabled is True
    assert "no longer exists" in sch._entries[ent.id].disabled_reason
    assert emitted, "a disabled schedule must reach the user"
    kind, title, detail = emitted[-1]
    assert kind == "run_failed"
    assert ent.id in title
    assert "no longer exists" in detail


def test_a_stale_one_shot_is_disabled_out_loud(cron, emitted):
    sch = S.Scheduler(path=cron)
    ent = sch.schedule_once(delay_seconds=60, prompt="look at it",
                            workspace=str(cron.parent))
    ent.next_fire_at = time.time() - S._STALE_ONCE_GRACE_S - 60

    sch.tick(fire_callback=lambda _e: None)

    assert ent.disabled is True
    assert emitted and ent.id in emitted[-1][1]


def test_repeated_failures_disable_out_loud(cron, emitted):
    sch = S.Scheduler(path=cron)
    ent = sch.schedule_interval(every_seconds=60, prompt="p",
                                fire_immediately=True,
                                workspace=str(cron.parent))

    def _boom(_entry):
        raise RuntimeError("backend down")

    for _ in range(S._MAX_CONSECUTIVE_FAILURES):
        sch.tick(fire_callback=_boom)

    assert ent.disabled is True
    disable_events = [e for e in emitted if "disabled" in e[1].lower()]
    assert disable_events, "the final, permanent stop must be announced"


def test_an_ordinary_failure_short_of_disabling_does_not_cry_wolf(
        cron, emitted):
    sch = S.Scheduler(path=cron)
    sch.schedule_interval(every_seconds=60, prompt="p", fire_immediately=True,
                          workspace=str(cron.parent))

    sch.tick(fire_callback=lambda _e: (_ for _ in ()).throw(OSError("blip")))

    assert [e for e in emitted if "disabled" in e[1].lower()] == []


def test_a_successful_fire_disables_nothing_and_says_nothing(cron, emitted):
    sch = S.Scheduler(path=cron)
    sch.schedule_interval(every_seconds=60, prompt="p", fire_immediately=True,
                          workspace=str(cron.parent))

    assert sch.tick(fire_callback=lambda _e: None) == 1
    assert emitted == []


def test_a_broken_notifier_never_stops_the_disable(cron, monkeypatch):
    """The notification is best-effort; the disable is not."""
    monkeypatch.setattr(
        "delfin.agent.attention.emit_attention",
        lambda *a, **kw: (_ for _ in ()).throw(OSError("inbox gone")))
    sch = S.Scheduler(path=cron)
    ent = sch.schedule_interval(every_seconds=60, prompt="p",
                                fire_immediately=True, workspace="/nowhere")

    sch.tick(fire_callback=SD.make_fire_callback(log=lambda _m: None))

    assert ent.disabled is True


# ---------------------------------------------------------------------------
# The save is atomic and does not flatten another writer
# ---------------------------------------------------------------------------

def test_the_schedule_is_written_atomically(cron):
    sch = S.Scheduler(path=cron)
    sch.schedule_once(delay_seconds=60, prompt="p", workspace=str(cron.parent))

    src = Path(S.__file__).read_text(encoding="utf-8")
    body = src[src.index("    def _save(self)"):]
    body = body[:body.index("    # --- disable notification")]
    assert "os.replace" in body, (
        "a reader must never catch a half-written schedule")
    assert "self.path.write_text" not in body
    assert json.loads(cron.read_text())["entries"]


def test_a_save_does_not_delete_an_entry_another_process_created(cron):
    """Scheduling into a file that another writer silently reverts."""
    first = S.Scheduler(path=cron)
    first.schedule_once(delay_seconds=600, prompt="mine",
                        workspace=str(cron.parent))

    second = S.Scheduler(path=cron)          # e.g. the headless daemon
    second.schedule_once(delay_seconds=600, prompt="theirs",
                         workspace=str(cron.parent))

    # The first instance, which has never seen "theirs", saves again.
    first.schedule_once(delay_seconds=600, prompt="mine again",
                        workspace=str(cron.parent))

    prompts = {e["prompt"] for e in json.loads(cron.read_text())["entries"]}
    assert prompts == {"mine", "theirs", "mine again"}


def test_a_deleted_entry_stays_deleted_across_the_merge(cron):
    sch = S.Scheduler(path=cron)
    ent = sch.schedule_once(delay_seconds=600, prompt="p",
                            workspace=str(cron.parent))

    assert sch.delete(ent.id) is True
    sch.schedule_once(delay_seconds=600, prompt="q", workspace=str(cron.parent))

    prompts = [e["prompt"] for e in json.loads(cron.read_text())["entries"]]
    assert prompts == ["q"]


def test_an_entry_another_process_created_can_be_deleted(cron):
    first = S.Scheduler(path=cron)
    other = S.Scheduler(path=cron)
    ent = other.schedule_once(delay_seconds=600, prompt="theirs",
                              workspace=str(cron.parent))

    assert first.delete(ent.id) is True
    assert json.loads(cron.read_text())["entries"] == []


def test_a_long_lived_scheduler_picks_up_a_new_entry(cron):
    """It loaded once at construction and never looked at the file again."""
    running = S.Scheduler(path=cron)
    S.Scheduler(path=cron).schedule_interval(
        every_seconds=60, prompt="added later", fire_immediately=True,
        workspace=str(cron.parent))

    fired: list = []

    assert running.tick(fire_callback=fired.append) == 1
    assert [e.prompt for e in fired] == ["added later"]


def test_a_fired_one_shot_does_not_come_back_from_the_file(cron):
    sch = S.Scheduler(path=cron)
    sch.schedule_once(delay_seconds=1, prompt="once",
                      workspace=str(cron.parent))
    time.sleep(1.1)

    assert sch.tick(fire_callback=lambda _e: None) == 1
    assert sch.tick(fire_callback=lambda _e: None) == 0
    assert json.loads(cron.read_text())["entries"] == []
    assert S.Scheduler(path=cron).list_entries() == []
