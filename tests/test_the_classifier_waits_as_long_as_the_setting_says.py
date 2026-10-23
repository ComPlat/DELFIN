"""How long the action classifier waits for a verdict.

Input: ``agent.classifier.timeout_s`` in the settings (absent, a number,
out of bounds, a boolean, text, NaN) and a classifier runner that answers
after a delay, driven through the real second-opinion path
(``_DocToolExecutor._classifier_second_opinion``) for a message to
another agent -- an action that leaves the workspace, so a missing
verdict asks. No test contacts a model.

Output, pinned here:

* without the setting the wait is ``_CLASSIFIER_TIMEOUT_S`` (20 s);
* the setting is clamped to ``_CLASSIFIER_TIMEOUT_BOUNDS_S`` (5..120 s);
  a boolean, text or NaN falls back to the default;
* a runner slower than the default but faster than the setting gives
  its verdict (allow) and asks nobody; the same runner without the
  setting times out and the person is asked -- so a longer wait spares
  a question and never lets more through;
* the wait applies to the client's HTTP timeout too.

Measured in wave 16 (2026-10-16): with nine sessions on the KIT proxy
the classifier passed 20 s for most messages between sessions, and each
became a question for the person.
"""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import api_client as A


ALLOW = "DECISION: allow\nREASON: part of the task."
TASK = [{"role": "user", "content": "review the builder's branch"}]


class Runner:
    def __init__(self, answer, delay):
        self.answer = answer
        self.delay = delay
        self.calls = 0

    def __call__(self, prompt, model):
        self.calls += 1
        time.sleep(self.delay)
        return self.answer


class Dialog:
    last_timed_out = False

    def __init__(self):
        self.shown = []

    def ask(self, name, args, preview):
        self.shown.append((name, preview))
        return True


@pytest.fixture
def settings(tmp_path, monkeypatch):
    """A settings file of its own; write(value) sets timeout_s."""
    path = tmp_path / "settings.json"

    def write(value=None):
        cfg = {} if value is None else {"timeout_s": value}
        path.write_text(json.dumps({"agent": {"classifier": cfg}}))

    write()
    import delfin.user_settings as us
    real = us.load_settings
    monkeypatch.setattr(us, "load_settings",
                        lambda settings_path=None: real(path))
    return write


@pytest.mark.parametrize("value, expected", [
    (None, A._CLASSIFIER_TIMEOUT_S),
    (60, 60.0),
    (45.5, 45.5),
    (1, A._CLASSIFIER_TIMEOUT_BOUNDS_S[0]),
    (10_000, A._CLASSIFIER_TIMEOUT_BOUNDS_S[1]),
    (True, A._CLASSIFIER_TIMEOUT_S),
    ("soon", A._CLASSIFIER_TIMEOUT_S),
    (float("nan"), A._CLASSIFIER_TIMEOUT_S),
])
def test_the_setting_is_read_and_bounded(settings, value, expected):
    settings(value)
    assert A._classifier_timeout_s() == expected


def _second_opinion(tmp_path, monkeypatch, runner, dialog):
    monkeypatch.setattr(A, "_CLASSIFIER_RETRY_DELAY_S", 0.0, raising=False)
    monkeypatch.setattr(
        A, "_classifier_settings",
        lambda: {"model": "m", "rules": [], "broad": False})
    ws = tmp_path / "ws"
    ws.mkdir(exist_ok=True)
    perms = A.KitToolPermissions(workspace=ws.resolve(), mode="acceptEdits")
    perms.classifier_runner = runner
    perms.classifier_messages = TASK
    perms.confirm_callback = dialog.ask
    return A._DocToolExecutor()._classifier_second_opinion(
        perms, "session_message",
        {"to": "ec9c6393", "message": "phase 1 done, please review"},
        perms.mode)


def test_a_slow_verdict_within_the_setting_asks_nobody(
        tmp_path, monkeypatch, settings):
    # Scaled down: default 0.2 s, setting 0.2 s * 5 (the lower bound is
    # lifted for the test so it stays fast), runner 0.4 s.
    monkeypatch.setattr(A, "_CLASSIFIER_TIMEOUT_S", 0.2)
    monkeypatch.setattr(A, "_CLASSIFIER_TIMEOUT_BOUNDS_S", (0.1, 120.0))
    settings(1.0)
    dialog = Dialog()
    runner = Runner(ALLOW, delay=0.4)
    out = _second_opinion(tmp_path, monkeypatch, runner, dialog)
    assert out is None
    assert dialog.shown == []
    assert runner.calls == 1


def test_the_same_verdict_without_the_setting_times_out_and_asks(
        tmp_path, monkeypatch, settings):
    monkeypatch.setattr(A, "_CLASSIFIER_TIMEOUT_S", 0.2)
    settings(None)
    dialog = Dialog()
    runner = Runner(ALLOW, delay=0.4)
    _second_opinion(tmp_path, monkeypatch, runner, dialog)
    assert len(dialog.shown) == 1
    assert "watcher" in dialog.shown[0][1].lower()


def test_a_refusal_within_the_setting_is_still_a_refusal(
        tmp_path, monkeypatch, settings):
    monkeypatch.setattr(A, "_CLASSIFIER_TIMEOUT_S", 0.2)
    monkeypatch.setattr(A, "_CLASSIFIER_TIMEOUT_BOUNDS_S", (0.1, 120.0))
    settings(1.0)
    dialog = Dialog()
    runner = Runner("DECISION: ask\nREASON: not part of the task.", delay=0.4)
    _second_opinion(tmp_path, monkeypatch, runner, dialog)
    assert len(dialog.shown) == 1


def test_the_http_timeout_follows_the_setting(monkeypatch, settings):
    settings(60)
    seen = {}

    class Completions:
        def create(self, **kw):
            class R:
                choices = [type("C", (), {"message": type(
                    "M", (), {"content": "DECISION: allow"})()})()]
            return R()

    class Client:
        def with_options(self, **kw):
            seen.update(kw)
            return type("X", (), {"chat": type(
                "Y", (), {"completions": Completions()})()})()

    client = A.OpenAIClient.__new__(A.OpenAIClient)
    client.client = Client()
    client.model = "m"
    client._provider = "kit"
    client._classifier_complete("p", "m")
    assert seen["timeout"] == 60.0
