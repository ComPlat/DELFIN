"""The Pi's surfaces speak branches: /commands in the TUI, words on the phone.

Input: a fake git checkout (tests/test_the_pi_runs_branches_like_the_
dashboard_does.py's fixture), the REPL command table and a Rooms gateway
rig (tests/test_the_messenger_gateway_only_serves_the_paired_phone.py's
fake bridge). Output:

* ``/branches [view]``, ``/test <n|name>``, ``/pull``, ``/cleanup`` answer
  with lines; numbering matches the order ``/branches`` printed them, so
  ``/test 2`` tests the second branch shown;
* the messenger answers the short words ``branches``, ``test <n>``,
  ``pull``, ``cleanup`` the same way, numbered answers included;
* nothing here switches the shared checkout; a deletion is only made by
  /cleanup with its confirmation (the phone: ``cleanup`` then ``ok``).
"""

from __future__ import annotations

import asyncio
import os
import subprocess
from pathlib import Path

import pytest


def _run(cwd, *args, check=True):
    out = subprocess.run(["git", "-c", "user.email=t@t", "-c", "user.name=t",
                          *args], cwd=str(cwd), capture_output=True, text=True,
                         env={**os.environ, "GIT_CONFIG_GLOBAL": os.devnull,
                              "GIT_CONFIG_SYSTEM": os.devnull})
    if check and out.returncode != 0:
        raise AssertionError(f"git {' '.join(args)}: {out.stderr}")
    return out.stdout.strip()


def _commit(cwd, name):
    (cwd / name).write_text(name + "\n")
    _run(cwd, "add", name)
    _run(cwd, "commit", "-q", "-m", name)


@pytest.fixture
def repo(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    origin = tmp_path / "origin.git"
    _run(tmp_path, "init", "-q", "--bare", "-b", "main", str(origin))
    co = tmp_path / "checkout"
    _run(tmp_path, "clone", "-q", str(origin), str(co))
    _commit(co, "a.txt")
    _run(co, "push", "-q", "-u", "origin", "main")
    _run(co, "branch", "old/merged")
    _run(co, "checkout", "-q", "-b", "feat")
    _commit(co, "f.txt")
    _run(co, "checkout", "-q", "main")
    return co


# -- TUI slash commands ------------------------------------------------------

def _ctx(repo):
    from delfin.agent.repl_commands import ReplContext
    return ReplContext(workspace=Path(repo))


def test_slash_branches_lists_and_numbers(repo):
    from delfin.agent import repl_commands
    result = repl_commands.dispatch("/branches", _ctx(repo))
    assert result.handled and result.output
    assert "main" in result.output
    # the numbered list: /test N can name one of these
    assert "1." in result.output or "1)" in result.output
    for view in ("all", "hidden", "merged", "remote"):
        result = repl_commands.dispatch(f"/branches {view}", _ctx(repo))
        assert result.handled and result.output


def test_slash_test_by_number_and_name(repo, monkeypatch):
    from delfin.agent import repl_commands
    from delfin.agent import branch_ops
    calls = []
    monkeypatch.setattr(branch_ops, "run_tests",
                        lambda *a, **k: calls.append(a) or
                        {"passed": 1, "failed": 0, "first_failure": ""})
    listing = repl_commands.dispatch("/branches", _ctx(repo)).output
    # the second row shown is a real branch of this repo: /test 2 resolves it
    result = repl_commands.dispatch("/test 2", _ctx(repo))
    assert result.handled
    assert "feat" in result.output or "old/merged" in result.output, listing
    assert calls, "the test must actually run"
    calls.clear()
    result = repl_commands.dispatch("/test feat", _ctx(repo))
    assert result.handled and calls


def test_slash_test_refuses_out_of_range(repo):
    from delfin.agent import repl_commands
    result = repl_commands.dispatch("/test 99", _ctx(repo))
    assert result.handled
    assert "no branch" in result.output.lower() or "unknown" in \
        result.output.lower() or "not" in result.output.lower()


def test_slash_pull(repo, monkeypatch):
    from delfin.agent import repl_commands
    from delfin.agent import branch_ops
    seen = []
    real = branch_ops.pull_main

    def spy(repo_path):
        seen.append(Path(repo_path))
        return real(repo_path)

    monkeypatch.setattr(branch_ops, "pull_main", spy)
    result = repl_commands.dispatch("/pull", _ctx(repo))
    assert result.handled and result.output
    assert seen == [Path(repo)]
    # the checkout stayed on main: a pull never switches
    assert _run(repo, "branch", "--show-current") == "main"


def test_slash_cleanup_lists_then_needs_confirm(repo):
    from delfin.agent import repl_commands
    listing = repl_commands.dispatch("/cleanup", _ctx(repo))
    assert listing.handled and listing.output
    assert "old/merged" in listing.output
    # a bare yes is not a confirmation of anything: /cleanup ok does it
    result = repl_commands.dispatch("/cleanup ok", _ctx(repo))
    assert result.handled
    assert _run(repo, "branch", "--list", "old/merged") == ""


def test_slash_cleanup_unmerged_needs_second_confirm(repo):
    from delfin.agent import repl_commands
    _run(repo, "checkout", "-q", "-b", "wip")
    _commit(repo, "w.txt")
    _run(repo, "checkout", "-q", "main")
    listing = repl_commands.dispatch("/cleanup", _ctx(repo)).output
    # unmerged branches are offered only with include_unmerged; the plain
    # listing does not offer wip, and /cleanup wip is not a deletion route
    assert "wip" not in listing
    result = repl_commands.dispatch("/cleanup wip", _ctx(repo))
    assert result.handled
    assert _run(repo, "branch", "--list", "wip") != ""


# -- messenger commands ------------------------------------------------------

class Daemon:
    def __init__(self):
        self.sent = []

    async def send_chat(self, chat, text):
        self.sent.append((chat, text))
        return {"chatItems": [{"chatItem": {"meta": {"itemId":
                                                     len(self.sent)}}}]}

    async def create_group(self, title):
        return f"#{title}"


class Bridge:
    def __init__(self):
        self.calls = []

    def ask(self, command):
        self.calls.append(command)
        if command["cmd"] == "sessions":
            return {"sessions": [{"id": "s1", "name": "t", "state": "idle",
                                  "messengerEnabled": True}]}
        if command.get("action") == "status":
            return {"ok": True, "title": "t", "state": "idle"}
        if command.get("action") == "prompt":
            return {"ok": True, "accepted": True}
        return {"ok": True}


@pytest.fixture
def rig(tmp_path, monkeypatch, repo):
    from delfin.messenger import store
    monkeypatch.setenv("DELFIN_MESSENGER_DIR", str(tmp_path))
    store.save("7")
    store.put_room("#10", "s1", "t", owner="7", ready=True)
    daemon, bridge = Daemon(), Bridge()
    from delfin.messenger.gateway import Rooms
    rooms = Rooms(daemon, bridge)
    # The session's workspace is the fake checkout: the branch commands run
    # in the session's workspace, where the service record says it is.
    monkeypatch.setattr(rooms, "session", fake_session_with(repo))
    return rooms, daemon, bridge, repo


def fake_session_with(repo):
    async def session(sid, action, **kwargs):
        if action == "branch":
            return kwargs["reply"]
        if action == "status":
            return {"ok": True, "title": "t", "state": "idle"}
        return {"ok": True}

    return session


def message(chat="#10", text="branches", sender="7", mid="1"):
    return {"chat": chat, "text": text, "sender": sender, "id": mid}


def test_messenger_branches_lists(rig):
    rooms, daemon, bridge, repo = rig
    bridge.reply_map = {}
    asyncio.run(rooms.handle(message(text="branches")))
    text = daemon.sent[-1][1]
    assert "main" in text
    assert "1." in text or "1)" in text


def test_messenger_test_by_number(rig, monkeypatch):
    rooms, daemon, bridge, repo = rig
    # the numbered answer the service gives for `test 2`
    reply = {"ok": True, "text": "fix-it: 3 passed, 0 failed in 5.0s",
             "passed": 3, "failed": 0}
    calls = []
    real_session = fake_session_with(repo)

    async def session(sid, action, **kwargs):
        if action == "branch":
            calls.append(kwargs)
            kwargs["reply"] = reply
            return reply
        return await real_session(sid, action, **kwargs)

    monkeypatch.setattr(rooms, "session", session)
    asyncio.run(rooms.handle(message(text="test 2")))
    assert calls and calls[0].get("action") == "test"
    text = daemon.sent[-1][1]
    assert "passed" in text


def test_messenger_pull(rig):
    rooms, daemon, bridge, repo = rig
    asyncio.run(rooms.handle(message(text="pull")))
    assert daemon.sent[-1][1]


def test_messenger_cleanup_two_step(rig):
    rooms, daemon, bridge, repo = rig
    # step 1: the listing
    asyncio.run(rooms.handle(message(text="cleanup", mid="m1")))
    listing = daemon.sent[-1][1]
    assert "old/merged" in listing
    # step 2: ok <number> is the confirmation; without it nothing is deleted
    asyncio.run(rooms.handle(message(text="ok", mid="m2")))
    assert _run(repo, "branch", "--list", "old/merged") != ""
    asyncio.run(rooms.handle(message(text="ok 1", mid="m3")))
    assert _run(repo, "branch", "--list", "old/merged") == ""
