"""A TUI session on a backend that announces its own id starts all the same.

`delfin-agent tui` without --provider ran on the Claude CLI (no settings
file, no API key), whose session id arrives only with the first turn's
stream. The engine's start id was "", the worker keyed its event log by it
and ended at once ("not a usable session id: ''"), so the window said only
"the worker exited (1)".

* the Claude CLI is told the id up front (``--session-id``), so the
  session has one from the start and resumes by the same id later;
* a backend that cannot be told one still gets a worker with an address
  and a log, and keeps its own (empty) engine id for the stream to fill.
"""

from __future__ import annotations

import io
import json
import uuid
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent.service import protocol


# -- the worker --------------------------------------------------------------

def test_a_worker_starts_on_an_engine_with_no_id_yet(tmp_path, monkeypatch):
    from delfin.agent.service.fake_engine import FakeEngine
    from delfin.agent.service.worker import Worker

    monkeypatch.setenv("DELFIN_AGENT_SERVICE_DIR", str(tmp_path / "service"))
    engine = FakeEngine(str(tmp_path))
    engine.session_id = ""          # as the CLI backends start

    worker = Worker(engine, writer=protocol.LineWriter(io.BytesIO()),
                    reader=io.BytesIO(), workspace=tmp_path, banner=[],
                    notices=[], fake=True, save=False)

    sid = worker.info()["sessionId"]
    assert sid, "a session must have an address before its first turn"
    assert worker.log.path.name == f"{sid}.jsonl"
    assert engine.session_id == "", (
        "the engine's id is the backend's to announce; a minted one would "
        "be passed to the backend as a session to resume")


def test_a_worker_starts_on_a_backend_with_a_gate_and_no_permissions(tmp_path, monkeypatch):
    """The Claude CLI: its file gate takes the confirm callback (the bind
    answers True), and there is no kit_permissions object behind it."""
    from delfin.agent.service.fake_engine import FakeEngine
    from delfin.agent.service.worker import Worker

    monkeypatch.setenv("DELFIN_AGENT_SERVICE_DIR", str(tmp_path / "service"))
    engine = FakeEngine(str(tmp_path))
    bound: list = []
    engine.set_kit_confirm_callback = lambda cb: bound.append(cb) or True
    engine.kit_permissions = None

    worker = Worker(engine, writer=protocol.LineWriter(io.BytesIO()),
                    reader=io.BytesIO(), workspace=tmp_path, banner=[],
                    notices=[], fake=True, save=False)

    assert bound == [worker.broker.callback]
    assert worker.info()["approval"] is not None


# -- the Claude CLI client ---------------------------------------------------

class _Proc:
    def __init__(self, lines=()):
        self.stdout = iter(lines)

    def poll(self):
        return None

    def wait(self, timeout=None):
        return 0


def _client(monkeypatch, tmp_path) -> tuple[A.CLIClient, list[list[str]]]:
    monkeypatch.setattr(A.shutil, "which", lambda p: p)
    client = A.CLIClient(claude_path="claude", cwd=str(tmp_path))
    monkeypatch.setattr(client, "_write_gate_state", lambda: None)
    monkeypatch.setattr(client, "_start_gate_relay", lambda: None)
    started: list[list[str]] = []

    def popen(cmd, **_kw):
        started.append(list(cmd))
        return _Proc()

    monkeypatch.setattr(A.subprocess, "Popen", popen)
    return client, started


def _flag(cmd: list[str], name: str) -> str | None:
    return cmd[cmd.index(name) + 1] if name in cmd else None


def test_a_new_cli_conversation_begins_under_the_id_it_was_given(monkeypatch, tmp_path):
    client, started = _client(monkeypatch, tmp_path)
    sid = str(uuid.uuid4())

    assert client.preset_session_id(sid) is True
    client._ensure_proc("system", session_id=sid)

    assert _flag(started[0], "--session-id") == sid
    assert "--resume" not in started[0], "there is nothing to resume yet"


def test_once_announced_the_cli_conversation_is_resumed(monkeypatch, tmp_path):
    client, started = _client(monkeypatch, tmp_path)
    sid = str(uuid.uuid4())
    client.preset_session_id(sid)
    client._proc = _Proc([json.dumps({"type": "system", "subtype": "init",
                                      "session_id": sid}),
                          json.dumps({"type": "result", "session_id": sid,
                                      "usage": {}})])

    events = list(client._read_turn(client._proc))

    assert [e.text for e in events if e.type == "session_init"] == [sid]
    client._proc = None             # the process ended; the next turn respawns
    client._ensure_proc("system", session_id=sid)
    assert _flag(started[-1], "--resume") == sid
    assert "--session-id" not in started[-1]


def test_a_live_cli_conversation_is_not_given_another_id(monkeypatch, tmp_path):
    client, _ = _client(monkeypatch, tmp_path)
    client._session_id = str(uuid.uuid4())

    assert client.preset_session_id(str(uuid.uuid4())) is False


# -- the engine ---------------------------------------------------------------

def _engine(ws: Path):
    from delfin.agent.engine import AgentEngine
    return AgentEngine(repo_dir=ws, backend="api", provider="kit",
                       api_key="test-key-not-used",
                       model="kit.qwen3.5-397b-A17b", mode="solo")


class _Announcing:
    supplies_session_id = True

    def __init__(self, takes: bool):
        self.takes = takes
        self.given: list[str] = []

    def preset_session_id(self, sid: str) -> bool:
        self.given.append(sid)
        return self.takes


@pytest.mark.parametrize("takes", [True, False])
def test_the_engine_gives_an_announcing_backend_its_id_where_it_can(tmp_path, takes):
    eng = _engine(tmp_path)
    eng.client = _Announcing(takes)

    sid = eng._fresh_session_id()

    if takes:
        assert sid == eng.client.given[0]
        assert str(uuid.UUID(sid)) == sid, "the CLI takes a canonical UUID only"
    else:
        assert sid == "", "the stream announces it; a guess would split it"


def test_the_claude_cli_client_takes_an_id_and_codex_does_not():
    assert A.CLIClient.preset_session_id is not A._BaseClient.preset_session_id
    assert A.CodexCLIClient.preset_session_id is A._BaseClient.preset_session_id
