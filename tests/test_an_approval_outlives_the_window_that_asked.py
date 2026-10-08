"""An approval can be answered after the window that asked for it is gone.

The dashboard's confirm dialog was answerable in exactly one place: the
widget that rendered it. ``KitConfirmBroker`` parks the request on a
``threading.Event`` in the tab's own state, so a browser reload -- which
builds a new tab around a NEW broker while the worker thread is still
blocked on the old one's Event -- left the request unanswerable. Five
minutes later it died of the timeout and was recorded as absence: the
action blocked, and the user never saw why.

The room this uses is not new and neither are its checks. A terminal
session already publishes its whole question to
``~/.delfin/terminal_confirmations``, where ``file_confirm``'s single
checked reader answers it -- not a symlink, this user's, not
group-writable, naming this very question, not older than it -- and
``delfin-agent approvals`` lists and answers what is there. The
dashboard simply was not in that room. One publish puts it there, and
four surfaces start covering it.

Driven through the real ``callback`` on a worker thread, because the
property is about what the blocked thread returns, and a test that calls
the publisher and reads the room would pass while the thread stayed
stuck.
"""

from __future__ import annotations

import threading
import time

import pytest

from delfin.agent import kit_confirm as kc
from delfin.agent import terminal_confirm as tc


@pytest.fixture
def room(tmp_path, monkeypatch):
    """The approval room, in this test's own directory."""
    monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "pending")
    return tmp_path / "pending"


class _Asked:
    """A request made on a worker thread, as the engine makes it."""

    def __init__(self, broker, tool_name="bash", args=None, preview="rm -rf x",
                 timeout_s=8.0):
        self.broker = broker
        self.answer = None
        self.returned = threading.Event()
        broker._timeout_s = timeout_s
        self._thread = threading.Thread(
            target=self._run, args=(tool_name, dict(args or {"command": "ls"}),
                                    preview),
            daemon=True)

    def _run(self, tool_name, args, preview):
        try:
            self.answer = self.broker.callback(tool_name, args, preview)
        finally:
            self.returned.set()

    def start(self, *, published: bool = True):
        """Start the worker and wait until the request is really parked.

        Waited on rather than slept past, so the test does not depend on
        a machine's speed. `published` waits for the room entry too: the
        request reaches the pending queue BEFORE it is published, so a
        test that read the room after seeing the queue read it too early
        -- which is what the first run of this file did.
        """
        self._thread.start()
        for _ in range(400):
            if self.broker.pending_requests():
                if not published or tc.pending_for_supervisors():
                    return self
            time.sleep(0.01)
        raise AssertionError(
            "the request never reached "
            + ("the approval room" if published else "the broker"))

    def wait(self, timeout=10.0) -> bool:
        return self.returned.wait(timeout=timeout)


def _published(room) -> dict:
    rows = tc.pending_for_supervisors()
    assert rows, f"nothing published; room holds {list(room.glob('*'))}"
    assert len(rows) == 1, rows
    return rows[0]


class TestTheRequestReachesTheRoom:
    def test_it_is_published_and_says_which_surface_waits(self, room):
        asked = _Asked(kc.KitConfirmBroker()).start()
        try:
            row = _published(room)
            assert tc.surface_of(row) == tc.DASHBOARD
            assert row.get("tool") == "bash"
            assert row.get("preview") == "rm -rf x"
            assert row.get("command") == "ls"
            assert row.get("pid") and row.get("proc_start"), (
                "without these a supervisor cannot tell a waiting session "
                "from one that died waiting")
        finally:
            asked.broker.resolve_all(False)
            asked.wait()

    def test_a_dashboard_request_is_not_read_as_a_terminals(self, room):
        """`pending_at_terminals` feeds cli_resume and approval_answers,
        which answer questions about TERMINAL sessions. A dashboard row
        under that name would be a different answer to the same question.
        """
        asked = _Asked(kc.KitConfirmBroker()).start()
        try:
            assert tc.pending_at_terminals() == []
            assert len(tc.pending_in_dashboards()) == 1
            assert len(tc.pending_for_supervisors()) == 1
        finally:
            asked.broker.resolve_all(False)
            asked.wait()

    def test_a_record_without_the_field_still_reads_as_a_terminals(self):
        """One may be in flight, written before the field existed."""
        assert tc.surface_of({}) == tc.TERMINAL
        assert tc.surface_of({"surface": ""}) == tc.TERMINAL
        assert tc.surface_of({"surface": "DASHBOARD"}) == tc.DASHBOARD


class TestItCanBeAnsweredFromOutside:
    def test_an_approval_from_the_room_unblocks_the_worker(self, room):
        """The property the whole change is for."""
        asked = _Asked(kc.KitConfirmBroker()).start()
        row = _published(room)

        assert tc.answer_waiting(row["id"], True, by="supervisor")
        assert asked.wait(), "the worker is still blocked on its own Event"
        assert asked.answer is True
        assert not asked.broker.last_timed_out

    def test_a_refusal_from_the_room_is_a_refusal_not_an_absence(self, room):
        """A denial and a timeout are different facts: an expiry sets
        ``last_timed_out`` so the model is told "ask later", and a real
        refusal must not borrow that."""
        broker = kc.KitConfirmBroker()
        asked = _Asked(broker).start()
        row = _published(room)

        assert tc.answer_waiting(row["id"], False, by="supervisor",
                                 reason="not that directory")
        assert asked.wait()
        assert asked.answer is False
        settled = broker.history()[-1]
        assert settled.answered_from_outside is True
        assert settled.outside_reason == "not that directory", (
            "the reason reaches the model in the same turn as the refusal, "
            "so it does not have to guess the same path differently")
        assert not broker.last_timed_out

    def test_the_record_comes_out_of_the_room_once_it_is_settled(self, room):
        asked = _Asked(kc.KitConfirmBroker()).start()
        row = _published(room)
        assert tc.answer_waiting(row["id"], True)
        assert asked.wait()
        assert tc.pending_for_supervisors() == [], (
            "a settled question read as open is worse than none at all")

    def test_a_click_still_answers_and_also_clears_the_room(self, room):
        """The widget is not replaced by this, only joined."""
        broker = kc.KitConfirmBroker()
        asked = _Asked(broker).start()
        assert _published(room)
        broker.resolve_all(True)
        assert asked.wait()
        assert asked.answer is True
        assert tc.pending_for_supervisors() == []


class TestTheFirstAnswerDecides:
    def test_a_late_click_does_not_overturn_a_refusal_from_the_room(self, room):
        """Two channels, one decision.

        A window that was replaced after a reload can still be clicked;
        its click must not turn a supervisor's refusal into an approval.
        """
        broker = kc.KitConfirmBroker()
        asked = _Asked(broker).start()
        row = _published(room)

        assert tc.answer_waiting(row["id"], False, by="supervisor")
        assert asked.wait()
        assert asked.answer is False

        # The stale panel, clicking Approve on a request already decided.
        req = broker.history()[-1]
        broker.resolve(req.seq, True)
        assert req.decision is False, (
            "the first answer must take it; this is an approval path")


class TestTheDialogWorksWithoutTheRoom:
    def test_an_unwritable_room_does_not_break_the_prompt(self, tmp_path,
                                                          monkeypatch):
        """The same rule the attention event follows: the dialog has to
        work on a host where nothing can be published."""
        monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "pending")

        def _no(*a, **k):
            raise OSError("read-only")

        monkeypatch.setattr(tc, "publish_waiting", _no)
        broker = kc.KitConfirmBroker()
        asked = _Asked(broker).start(published=False)
        assert broker.pending_requests(), "the request still has to be asked"
        broker.resolve_all(True)
        assert asked.wait()
        assert asked.answer is True
        assert broker.history()[-1].published is None
