"""A user who is away costs the agent one approval window, not one each.

The broker waits five minutes for a decision and then treats the silence
as absence rather than refusal -- which is right. What it also did was
wait the full five minutes again for the next request, and the one after
that, although the same nobody was going to answer. Two sessions on
2026-09-17 spent a window each that way, one of them on a `git add`.

So the first expiry records the absence, and while it stands, further
requests are answered at once with the same "nobody is there" verdict.
The inbox entry is still written, so the user sees what was wanted when
they come back, and one real decision ends the absence.

  the first request waits         the window is what it always was
  the next one does not           same verdict, no second window
  it is absence, not refusal      last_timed_out stays true
  an answer ends it               the next request waits again
  the window is a window          after it passes, waiting resumes
"""

from __future__ import annotations

import threading
import time

import pytest

from delfin.agent.kit_confirm import KitConfirmBroker


@pytest.fixture
def broker(monkeypatch):
    # No inbox writing in a unit test: emit_attention reaches the file
    # system, and this is about the wait.
    import delfin.agent.attention as attention
    monkeypatch.setattr(attention, "emit_attention", lambda *a, **k: "")
    return KitConfirmBroker(default_timeout_s=0.3)


def _ask(broker) -> tuple[bool, float]:
    started = time.monotonic()
    answer = broker.callback("bash", {"command": "git push"}, "$ git push")
    return answer, time.monotonic() - started


def test_the_first_request_waits_its_window(broker):
    answer, waited = _ask(broker)
    assert answer is False
    assert broker.last_timed_out is True
    assert waited >= 0.3


def test_the_next_request_does_not_wait_again(broker):
    _ask(broker)
    answer, waited = _ask(broker)
    assert answer is False
    assert waited < 0.1, "a second window for the same absent user"
    # Still absence, not a refusal: the caller must keep telling the model
    # "the user is away", not "the user said no".
    assert broker.last_timed_out is True


def test_an_answer_ends_the_absence(broker):
    _ask(broker)

    def _answer_it():
        for _ in range(100):
            with broker._lock:
                pending = list(broker._pending)
            if pending:
                pending[0].decision = True
                pending[0].event.set()
                return
            time.sleep(0.01)

    # The absence has to pass first, or the next request would not wait at
    # all and there would be nothing to answer.
    time.sleep(0.31)
    threading.Thread(target=_answer_it, daemon=True).start()
    answer, _waited = _ask(broker)
    assert answer is True
    assert broker.last_timed_out is False

    # And with somebody there, the request after it waits again.
    _answer, waited = _ask(broker)
    assert waited >= 0.3


def test_the_absence_lasts_one_window_and_no_longer(broker):
    _ask(broker)
    time.sleep(0.31)
    _answer, waited = _ask(broker)
    assert waited >= 0.3, "after the window passes, ask properly again"


def test_the_user_still_learns_what_was_wanted(monkeypatch):
    """The wait is skipped; the inbox entry is not."""
    seen = []
    import delfin.agent.attention as attention
    monkeypatch.setattr(attention, "emit_attention",
                        lambda kind, **kw: seen.append(kind) or "")
    broker = KitConfirmBroker(default_timeout_s=0.2)
    broker.callback("bash", {"command": "git push"}, "$ git push")
    broker.callback("bash", {"command": "rm -rf build"}, "$ rm -rf build")
    assert seen == ["confirm_pending", "confirm_pending"]
