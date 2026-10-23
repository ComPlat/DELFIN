"""A broker question is published under the redirect, never in the home.

The measured incident (2026-09-21, 23:30, a compute node): a full suite
run left two ``$ ls`` confirmations "waiting at a terminal" for the
operator -- no session attached, the host column naming the node. A test
had asked through ``TerminalConfirmBroker``; the question was published
under the user's REAL ``~/.delfin/terminal_confirmations``, because the
room is a module constant resolved at import and stood in none of the
tables the suite's redirect walks.

This file measures the case itself. It checks where the room points
BEFORE asking anything, so a red run never asks -- and never writes into
the real home, not even while failing.
"""

from __future__ import annotations

import threading
import time
from pathlib import Path

from delfin.agent import terminal_confirm as tc


def _real_room() -> Path:
    return Path.home() / ".delfin" / "terminal_confirmations"


def test_a_broker_question_lands_under_the_redirect():
    real_room = _real_room()

    # Checked BEFORE the question: if the room still points into the real
    # home, asking would publish the question there -- the very leak this
    # file measures. A red run fails here, having written nothing.
    assert tc._PENDING_DIR != real_room and real_room not in tc._PENDING_DIR.parents, (
        f"the confirmation room still points into the real home "
        f"({tc._PENDING_DIR}); a question asked now would be published "
        "there -- the incident of 2026-09-21 again")

    broker = tc.TerminalConfirmBroker(
        session_id="s4-measured-case", session_key="s4",
        poll_s=0.01)

    answer: dict = {}

    def _ask():
        answer["r"] = broker.ask_user({"question": "nothing"})

    t = threading.Thread(target=_ask, daemon=True)
    t.start()

    req = None
    for _ in range(500):                 # up to ~5 s, asked within ms
        req = broker.take()
        if req is not None:
            break
        time.sleep(0.01)
    assert req is not None, "the broker never took the question"

    # _enqueue appends the request BEFORE it publishes (terminal_confirm.py:
    # queue first, then _publish_pending), so the take above can win the
    # race and see the request half-built. Publishing takes milliseconds;
    # wait for it rather than reading the attribute at the worst moment.
    published = None
    for _ in range(500):
        published = req.published
        if published is not None:
            break
        time.sleep(0.01)
    assert published is not None, (
        "the question was not published at all -- the supervisor-readable "
        "copy (see _publish_pending) is part of what a question is")
    assert real_room not in published.parents, (
        f"the question was published into the real home ({published})")
    assert published.exists(), (
        "the published copy does not exist where it says it does")

    broker.resolve(req, {"answers": []})
    t.join(timeout=5)
    assert not t.is_alive(), "the asking thread never came back"
