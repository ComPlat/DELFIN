"""A yes is not asked for twice, and a running job is not started twice.

Session 78c36f8c (2026-10-19): `python3 -m http.server 8038 ...` was asked
about in three turns, approved each time -- the second opinion's
``allowed`` set is per turn, so the next message forgot the yes. The same
session started eight file servers for one game, each beside or over the
last; nothing said one was already up.
"""

from __future__ import annotations

import json
import sys
import time

import pytest

from delfin.agent import api_client as A
from delfin.agent import bash_jobs as BJ

SERVER = "python3 -m http.server 8038 --bind 127.0.0.1 --directory site"


# -- the user's yes holds for the session ----------------------------------

class _Dialog:
    last_timed_out = False

    def __init__(self):
        self.shown = []

    def ask(self, name, args, preview):
        self.shown.append(name)
        return True


@pytest.fixture
def perms(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_classifier_settings",
                        lambda: {"model": "m", "rules": [], "broad": False})
    ws = tmp_path / "wt"
    ws.mkdir()
    p = A.KitToolPermissions(workspace=ws.resolve(), mode="bypassPermissions")

    def _down(prompt, model):
        raise RuntimeError("endpoint error")

    p.classifier_runner = _down
    p.classifier_messages = [{"role": "user", "content": "zeig mir das Spiel"}]
    p._dialog = _Dialog()
    p.confirm_callback = p._dialog.ask
    return p


def _second_opinion(perms, cmd):
    return A._DocToolExecutor()._classifier_second_opinion(
        perms, "bash", {"command": cmd}, "bypassPermissions")


def test_a_yes_is_not_asked_again_in_a_later_turn(perms):
    assert _second_opinion(perms, SERVER) is None
    assert perms._dialog.shown == ["bash"]
    # A new user message: a new turn, the per-turn cache is cleared.
    perms.classifier_messages = perms.classifier_messages + [
        {"role": "assistant", "content": "läuft"},
        {"role": "user", "content": "starte es neu"}]
    assert _second_opinion(perms, SERVER) is None
    assert perms._dialog.shown == ["bash"]


def test_a_different_command_is_still_asked(perms):
    _second_opinion(perms, SERVER)
    _second_opinion(perms, SERVER.replace("8038", "9000"))
    assert perms._dialog.shown == ["bash", "bash"]


def test_a_no_is_not_kept_as_a_yes(perms):
    perms.confirm_callback = lambda *a: False
    assert _second_opinion(perms, SERVER) is not None
    assert not perms.classifier_user_allowed


# -- a running job is handed back, not started twice -----------------------

@pytest.fixture
def ws(tmp_path):
    (tmp_path / "site").mkdir()
    (tmp_path / "site" / "index.html").write_text("<p>hi</p>")
    return tmp_path


def _start(ws, command, sid="s-twin"):
    return BJ.get_registry().start(
        command=command, cwd=str(ws), workspace=str(ws), session_id=sid,
        argv=[sys.executable, "-c", "import time; time.sleep(60)"])


def _stop(job):
    BJ.get_registry().kill(job.job_id)
    deadline = time.time() + 10
    while job.poll() is None and time.time() < deadline:
        time.sleep(0.05)


def _perms_for(ws, sid="s-twin"):
    p = A.KitToolPermissions(workspace=str(ws), mode="acceptEdits")
    p.task_session_id = sid
    return p


def test_the_same_server_is_handed_back(ws):
    job = _start(ws, SERVER)
    try:
        twin = A._bj_running_twin(SERVER, str(ws), _perms_for(ws))
        assert twin["status"] == "already_running"
        assert twin["job_id"] == job.job_id
        assert twin["url"] == "http://127.0.0.1:8038/"
        assert "bash_kill" in twin["hint"]
    finally:
        _stop(job)


def test_a_second_server_on_the_same_address_is_handed_back(ws):
    job = _start(ws, SERVER)
    try:
        other = "python3 -m http.server 8038 --bind 127.0.0.1 --directory ."
        assert A._bj_running_twin(other, str(ws), _perms_for(ws))["job_id"] \
            == job.job_id
    finally:
        _stop(job)


def test_another_sessions_job_or_another_command_is_not(ws):
    job = _start(ws, SERVER, sid="someone-else")
    try:
        assert A._bj_running_twin(SERVER, str(ws), _perms_for(ws)) is None
    finally:
        _stop(job)
    job = _start(ws, "python run.py")
    try:
        assert A._bj_running_twin("python other.py", str(ws),
                                  _perms_for(ws)) is None
    finally:
        _stop(job)


def test_a_job_that_ended_is_started_again(ws):
    job = _start(ws, SERVER)
    _stop(job)
    assert A._bj_running_twin(SERVER, str(ws), _perms_for(ws)) is None
