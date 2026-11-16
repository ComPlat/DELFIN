"""A chain run through DELFIN's tool server is indexed under its session.

Input: the builtin delfin-tools MCP server, started as the agent starts
it (mcp_client.MCPServer with the builtin command), and run_application
calls carrying the dispatching session's id.

Before: the session id had no way into the server process. A spawn-time
environment variable would not have worked either: the server is
started while the tool list is assembled, before any dispatch, and one
server is shared by every session of a process and workspace.
Now the id travels with each tools/call request in ``_meta`` under
"delfin/agent_session_id", only to the builtin server; the server reads
it from the SDK's request metadata (a pydantic model, not a dict) and
passes it down to pipeline.run, which records it in the chain index.
A model-supplied agent_session_id argument never counts.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import mcp_client

_APP_MODULE = '''
from delfin.tools import (Application, Pipeline, register_application,
                          register_callable)

@register_callable("jc2_record_fn")
def _record(results, last, wd):
    return {"ok": True}

def _pipeline():
    pipe = Pipeline("jc2-record")
    pipe.add_compute(_record, label="record")
    return pipe

register_application(
    Application.from_pipeline(
        _pipeline(), name="jc2-record", description="jc2 recorder")
)
'''


@pytest.fixture()
def home(tmp_path, monkeypatch):
    home = tmp_path / "home"
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr("pathlib.Path.home", lambda: home)
    apps = home / ".delfin" / "applications"
    apps.mkdir(parents=True, exist_ok=True)
    (apps / "jc2_record_app.py").write_text(_APP_MODULE, encoding="utf-8")
    return home


@pytest.fixture()
def server(home):
    cfg = mcp_client._BUILTIN_SERVERS["delfin-tools"]
    srv = mcp_client.MCPServer(
        name="delfin-tools", command=cfg["command"], args=list(cfg["args"]),
        env={"DELFIN_APPLICATIONS_DIR":
             str(home / ".delfin" / "applications")},
    )
    assert srv.initialize(), srv.last_error
    try:
        yield srv
    finally:
        srv.stop()


def _index(home: Path) -> list[dict]:
    idx = home / ".delfin" / "chain_index.jsonl"
    assert idx.is_file(), "the run did not index its chain"
    return [json.loads(line) for line in
            idx.read_text(encoding="utf-8").splitlines() if line.strip()]


def test_two_sessions_on_one_server_each_get_their_own_id(home, server,
                                                          tmp_path):
    # The server is already up (started before any dispatch, as in a
    # session); two sessions then dispatch through the same server.
    for sid, run in (("sess-1", "a"), ("sess-2", "b")):
        out = server.call_tool(
            "run_application",
            {"name": "jc2-record", "inputs": {"work_dir": str(tmp_path / run)}},
            agent_session_id=sid)
        assert json.loads(out).get("ok"), out
    by_dir = {Path(e["work_dir"]).name: e.get("agent_session_id")
              for e in _index(home)}
    assert by_dir == {"a": "sess-1", "b": "sess-2"}


def test_a_run_without_a_session_still_completes_and_records_none(
        home, server, tmp_path):
    out = server.call_tool(
        "run_application",
        {"name": "jc2-record", "inputs": {"work_dir": str(tmp_path / "plain")}})
    assert json.loads(out).get("ok"), out
    assert [e.get("agent_session_id") for e in _index(home)] == [None]


def test_a_model_supplied_id_does_not_count(home, server, tmp_path):
    out = server.call_tool(
        "run_application",
        {"name": "jc2-record",
         "inputs": {"work_dir": str(tmp_path / "forged"),
                    "agent_session_id": "forged"}},
        agent_session_id="sess-real")
    assert json.loads(out).get("ok"), out
    assert [e.get("agent_session_id") for e in _index(home)] == ["sess-real"]


def test_only_the_builtin_server_receives_the_id(monkeypatch):
    sent = []

    def fake_send(self, method, params):
        sent.append((self.name, params))
        return {"result": {"content": [{"type": "text", "text": "ok"}]}}

    monkeypatch.setattr(mcp_client.MCPServer, "initialize", lambda self: True)
    monkeypatch.setattr(mcp_client.MCPServer, "_send", fake_send)
    for name in ("delfin-tools", "someone-elses-server"):
        mcp_client.MCPServer(name=name, command="true", args=[]).call_tool(
            "run_application", {"name": "x"}, agent_session_id="sess-A")
    mcp_client.MCPServer(name="delfin-tools", command="true",
                         args=[]).call_tool("run_application", {"name": "x"})
    assert sent[0][1]["_meta"] == {"delfin/agent_session_id": "sess-A"}
    assert "_meta" not in sent[1][1]
    assert "_meta" not in sent[2][1]
