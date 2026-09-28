"""`delfin-agent session export` writes a replayable notebook.

The CLI is a thin shell over ``session_export``; these tests pin the
argument surface (`session export <id> --notebook out.ipynb`) and the
error paths, never touching the real ``~/.delfin``.
"""
from __future__ import annotations

import json

import pytest

from delfin.agent import cli as agent_cli, session_store as ss, tool_trace

nbformat = pytest.importorskip("nbformat")


@pytest.fixture
def fake_home(tmp_path, monkeypatch):
    """Sessions + tool traces under tmp_path, never the real ~/.delfin."""
    d = tmp_path / "agent_sessions"
    monkeypatch.setattr(ss, "_SESSIONS_DIR", d)
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    return tmp_path


def _store_session(fake_home, session_id: str) -> None:
    ss.save_session(session_id, mode="solo", model="kit", provider="glm",
                    chat_messages=[
                        {"role": "user", "content": "build benzene"},
                        {"role": "assistant", "content": "running xtb"},
                    ], title="xtb test")
    p = tool_trace.trace_path(
        session_id, root=tmp_of(fake_home))
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({
        "ts": 1.0, "tool": "mcp__delfin-ops__smiles_to_xyz",
        "input": json.dumps({"smiles": "C1=CC=CC=C1"}),
        "output": "benzene.xyz", "ok": True}) + "\n",
        encoding="utf-8")


def tmp_of(fake_home):
    return fake_home / "tool_traces"


def test_session_export_writes_a_notebook(fake_home, tmp_path):
    _store_session(fake_home, "cli-exp-1")
    out = tmp_path / "out.ipynb"
    rc = agent_cli.main(
        ["session", "export", "cli-exp-1", "--notebook", str(out)])
    assert rc == 0
    nb = nbformat.read(out, as_version=4)
    nbformat.validate(nb)
    code = [c for c in nb.cells if c.cell_type == "code"]
    assert any("smiles_to_xyz" in c.source for c in code)


def test_session_export_defaults_to_stdout_message(fake_home, tmp_path,
                                                   monkeypatch):
    _store_session(fake_home, "cli-exp-2")
    monkeypatch.chdir(tmp_path)  # default output lands next to the caller
    rc = agent_cli.main(["session", "export", "cli-exp-2"])
    assert rc == 0


def test_session_export_unknown_id_fails_cleanly(fake_home, capsys):
    rc = agent_cli.main(["session", "export", "does-not-exist",
                         "--notebook", "x.ipynb"])
    assert rc == 1
    assert "not found" in capsys.readouterr().err.lower()


def test_session_export_accepts_latest(fake_home, tmp_path):
    _store_session(fake_home, "cli-exp-3")
    out = tmp_path / "out.ipynb"
    rc = agent_cli.main(["session", "export", "latest",
                         "--notebook", str(out)])
    assert rc == 0
    assert out.exists()
