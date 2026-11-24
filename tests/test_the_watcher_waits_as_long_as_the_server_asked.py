"""The watcher retried before the server's own "try again in 1s".

MEASURED on 2026-11-23, wave 40: six agents shared one KIT key, whose
binding limit is 60 requests a minute, and the action classifier (the
model's second opinion on an action that leaves the workspace) spends
from the same allowance. Its failure came first in the log:

    permission_denied classifier_held
        "classifier error: the model gave no answer"

KIT had answered ``429 - Requests-per-minute limit exceeded (60
req/min). Try again at ... UTC (in 1s)``. The watcher retries exactly
once, after a FIXED ``_CLASSIFIER_RETRY_DELAY_S`` of 0.5 s -- half the
wait the server named -- so the retry met the same limit and the action
was held for want of a verdict.

The main stream already obeys the opposite rule ("The server's own
Retry-After outranks the ladder: the ladder is a guess made here, the
header is the one number the service knows"). The watcher now does too:
a wait the server STATES -- a Retry-After header, or KIT's "(in Ns)" in
the body -- is waited before the one retry, and a plain per-minute 429
that states none waits a short rate-limit step instead of half a second.

TWO LIMITS STAY, and they are tested here because they are what keeps
the watcher safe:

  * the watcher sits in front of an action, so its wait is capped
    (``_CLASSIFIER_RETRY_WAIT_MAX_S``): a stated wait beyond it is not
    sat out, the call fails as before, and the action is held -- an
    unjudged action that leaves the workspace is never let through;
  * any other model error (a dropped connection, a 502) keeps the short
    0.5 s retry it had.
"""
from __future__ import annotations

import pytest

from delfin.agent import api_client as A

ALLOW = "DECISION: allow\nREASON: part of the task."
TASK = [{"role": "user", "content": "rebuild the report in this workspace"}]


class RateLimitError(Exception):
    """Shaped like the client's own: a status code and KIT's body."""
    status_code = 429

    def __init__(self, seconds):
        super().__init__(
            "Error code: 429 - Requests-per-minute limit exceeded (60 "
            f"req/min). Try again at 2026-11-23 18:47 UTC (in {seconds}s).")


class Endpoint:
    """A rate limit with a memory: every call before ``free_at`` (virtual
    seconds) is refused with the remaining wait, as KIT does; after it the
    model answers. Time advances only through the watcher's sleep."""

    def __init__(self, free_at, answer=ALLOW, exc=None):
        self.now = 0.0
        self.free_at = free_at
        self.answer = answer
        self.exc = exc
        self.calls = 0

    def runner(self, prompt, model):
        self.calls += 1
        if self.exc is not None and self.calls == 1:
            raise self.exc
        if self.now < self.free_at:
            remaining = max(1, int(round(self.free_at - self.now)))
            raise RateLimitError(remaining)
        return self.answer

    def sleep(self, seconds):
        self.now += seconds


@pytest.fixture(autouse=True)
def _a_host_without_a_cage(monkeypatch):
    monkeypatch.setattr(A, "_bash_runs_caged", lambda perms: False)
    monkeypatch.setattr(
        A, "_classifier_settings",
        lambda: {"model": "m", "rules": [], "broad": False})


def _unattended(tmp_path, endpoint, monkeypatch):
    monkeypatch.setattr(A, "_classifier_sleep", endpoint.sleep, raising=False)
    ws = tmp_path / "ws"
    ws.mkdir(exist_ok=True)
    perms = A.KitToolPermissions(workspace=ws.resolve(), mode="default")
    perms.classifier_runner = endpoint.runner
    perms.classifier_messages = TASK
    perms.confirm_callback = None
    return perms


def _gate(perms):
    # An action that LEAVES the workspace: a watcher failure holds it when
    # nobody can be asked, which is what the measured run hit.
    return A._DocToolExecutor()._run_permission_gate(
        "web_fetch", {"url": "https://github.com/a"}, perms)


# --- the case that was measured ----------------------------------------

def test_the_server_stated_wait_is_waited_and_the_action_is_judged(
        tmp_path, monkeypatch):
    ep = Endpoint(free_at=1.0)
    out = _gate(_unattended(tmp_path, ep, monkeypatch))
    assert out is None, f"held although the server said 1s: {out}"
    assert ep.now >= 1.0, f"retried after {ep.now}s, the server said 1s"
    assert ep.calls == 2, "exactly one retry, as before"


def test_a_retry_after_header_is_obeyed_too(tmp_path, monkeypatch):
    class _Resp:
        headers = {"retry-after": "3"}

    class Limited(Exception):
        status_code = 429
        response = _Resp()

    ep = Endpoint(free_at=3.0, exc=Limited("429 Too Many Requests"))
    out = _gate(_unattended(tmp_path, ep, monkeypatch))
    assert out is None, out
    assert ep.now >= 3.0, ep.now


def test_a_per_minute_limit_with_no_stated_wait_waits_a_rate_limit_step(
        tmp_path, monkeypatch):
    class Limited(Exception):
        status_code = 429

    ep = Endpoint(free_at=1.5, exc=Limited("429 Too Many Requests"))
    out = _gate(_unattended(tmp_path, ep, monkeypatch))
    assert out is None, out
    assert ep.now > A._CLASSIFIER_RETRY_DELAY_S


# --- the limits that keep it safe --------------------------------------

def test_a_wait_beyond_the_cap_is_not_sat_out_and_the_action_is_held(
        tmp_path, monkeypatch):
    stated = getattr(A, "_CLASSIFIER_RETRY_WAIT_MAX_S", 10.0) * 10
    ep = Endpoint(free_at=stated)
    out = _gate(_unattended(tmp_path, ep, monkeypatch))
    assert out is not None and "watcher" in out.lower(), out
    assert ep.now <= getattr(A, "_CLASSIFIER_RETRY_WAIT_MAX_S", 10.0), (
        f"the watcher sat {ep.now}s in front of an action")


def test_any_other_model_error_keeps_the_short_retry(tmp_path, monkeypatch):
    ep = Endpoint(free_at=0.0, exc=ConnectionError("endpoint refused"))
    out = _gate(_unattended(tmp_path, ep, monkeypatch))
    assert out is None, out
    assert ep.now == pytest.approx(A._CLASSIFIER_RETRY_DELAY_S)


def test_on_the_success_path_nothing_waits(tmp_path, monkeypatch):
    ep = Endpoint(free_at=0.0)
    assert _gate(_unattended(tmp_path, ep, monkeypatch)) is None
    assert ep.now == 0.0 and ep.calls == 1
