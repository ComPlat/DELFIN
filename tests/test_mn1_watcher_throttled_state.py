"""MN1 phase 3: an output watcher's state survives a throttled look.

watch_checks.check (output event, :301-302) answers THROTTLED with
``previous or ""``. Two lookups later that empty string replaces the
"never looked" sentinel None, and since the output event's state is the
SHA-256 of the output, "" is never equal to a real hash -- the next real
look then fires "output changed" on UNCHANGED output. A watcher that
fires on an event that did not happen starts a run (schedules.watch_seen
marks run_requested), so a throttled scheduler can start an agent turn
for nothing. Control: red before the fix on 3a6d7a66c.
"""
from __future__ import annotations

import pytest

from delfin.agent import schedules, watch_checks


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setattr(schedules, "_DIR", tmp_path / "user_schedules")
    monkeypatch.setattr(schedules, "_hostname", lambda: "here")
    folder = tmp_path / "ws"
    folder.mkdir()
    return folder


class _ThrottledClient:
    """A scheduler client that is always throttled, like the real one at
    peak: it answers :data:`THROTTLED` for every query."""

    def __init__(self):
        from delfin import scheduler_client
        self.throttled = scheduler_client.THROTTLED
        self.queries = 0

    def query(self, kind, argv, runner=None):
        self.queries += 1
        return self.throttled


def _output_spec():
    return {"kind": "watch", "what": "output", "target": "squeue -u max",
            "event": "changes", "every": 60, "text": "wenn ausgabe ...", "tz": None}


@pytest.fixture
def _squeue_on_path(monkeypatch, tmp_path):
    """``squeue`` must exist on PATH for _ask_scheduler to run at all; a
    tiny stand-in so the fake client's answer is reached."""
    stub = tmp_path / "squeue"
    stub.write_text("#!/bin/sh\necho nothing\n")
    stub.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path) + ":" +
                       _safe_orig_path())


def _safe_orig_path():
    import os
    return os.environ.get("PATH", "/usr/bin:/bin")


def test_throttled_look_keeps_the_no_state_sentinel(ws, monkeypatch, _squeue_on_path):
    from delfin import scheduler_client
    fake = _ThrottledClient()
    monkeypatch.setattr(scheduler_client, "default_client", fake)
    spec = _output_spec()
    # First look, throttled: the state must stay "never looked" (None),
    # not become "".
    state, fired, note = watch_checks.check(spec, str(ws), None)
    assert not fired
    assert state is None


def test_a_throttled_look_does_not_fake_a_change(ws, monkeypatch, _squeue_on_path):
    """The real sequence: look 1 learns state H, look 2 is throttled,
    look 3 sees the SAME output. With the bug, look 2 stored "" and look
    3 fires 'changed' on identical output."""
    from delfin import scheduler_client
    import hashlib
    out = "unsteady output\n"          # identical across all looks
    state_hash = hashlib.sha256(out.encode()).hexdigest()

    class _OneAnswerClient:
        throttled = scheduler_client.THROTTLED

        def query(self, kind, argv, runner=None):
            return out                 # always the same answer

    monkeypatch.setattr(scheduler_client, "default_client", _OneAnswerClient())
    spec = _output_spec()
    _, fired1, _ = watch_checks.check(spec, str(ws), None)
    assert not fired1
    _, fired2, _ = watch_checks.check(spec, str(ws), state_hash)
    assert not fired2, "a throttled look must not fabricate a state"
    _, fired3, _ = watch_checks.check(spec, str(ws), state_hash)
    assert not fired3, ("unchanged output after a throttled look must not "
                        "fire 'changed'")
