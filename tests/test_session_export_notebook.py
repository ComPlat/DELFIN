"""A session exported as a notebook keeps its chemistry steps replayable.

Each chemical step becomes a Markdown cell (what and why, taken from the
chat) and a code cell with the native DELFIN call reconstructed from the
tool trace. Internal agent tooling (tasks, memory) is left out, and every
text is scrubbed with ``memory_store._without_secrets``.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import session_export, tool_trace

nbformat = pytest.importorskip("nbformat")


def _write_trace(root: Path, session_id: str, entries: list[dict]) -> None:
    """Write a tool trace JSONL into ``root`` (never the real ~/.delfin)."""
    p = tool_trace.trace_path(session_id, root=root)
    p.parent.mkdir(parents=True, exist_ok=True)
    for e in entries:
        p.open("a", encoding="utf-8").write(json.dumps(e) + "\n")


def _fake_session(session_id: str) -> dict:
    return {
        "session_id": session_id,
        "title": "TADF screening",
        "workspace": "/tmp/proj",
        "chat_messages": [
            {"role": "user", "content": "Optimize this emitter at PBE0/def2-SVP"},
            {"role": "assistant", "content": "Building the geometry and running xtb."},
        ],
        "engine_messages": [
            {"role": "user", "content": "Optimize this emitter at PBE0/def2-SVP"},
            {"role": "assistant", "content": "Building the geometry and running xtb."},
        ],
        "tool_calls": 3,
    }


def _trace_entries() -> list[dict]:
    return [
        {  # chemistry: SMILES -> xyz
            "ts": 1.0, "tool": "mcp__delfin-ops__smiles_to_xyz",
            "input": json.dumps({"smiles": "C1=CC=CC=C1"}),
            "output": "wrote benzene.xyz", "ok": True,
        },
        {  # chemistry: xtb run
            "ts": 2.0, "tool": "mcp__delfin-ops__submit_calculation",
            "input": json.dumps({"folder": "calc/benzene", "engine": "xtb"}),
            "output": "job 12345", "ok": True,
        },
        {  # chemistry: parse
            "ts": 3.0, "tool": "mcp__delfin-ops__extract_energy_table",
            "input": json.dumps({"folders": ["calc/benzene"]}),
            "output": '{"benzene": -230.7}', "ok": True,
        },
        {  # internal agent tooling: excluded
            "ts": 4.0, "tool": "task_create",
            "input": json.dumps({"subject": "run opt"}),
            "output": "task 1", "ok": True,
        },
        {  # internal agent tooling: excluded
            "ts": 5.0, "tool": "remember",
            "input": json.dumps({"text": "prefers xtb"}),
            "output": "saved", "ok": True,
        },
        {  # generic bash: excluded (not chemistry)
            "ts": 6.0, "tool": "bash",
            "input": json.dumps({"command": "grep -rn foo delfin/"}),
            "output": "", "ok": True,
        },
    ]


def test_a_session_export_writes_an_nbformat_valid_notebook(tmp_path):
    session_id = "sess-export-1"
    _write_trace(tmp_path, session_id, _trace_entries())
    nb = session_export.export_session(
        _fake_session(session_id), trace_root=tmp_path)
    nbformat.validate(nb)
    code = [c for c in nb.cells if c.cell_type == "code"]
    md = [c for c in nb.cells if c.cell_type == "markdown"]
    # three chemistry steps -> three code cells, each paired with markdown
    assert len(code) == 3
    assert len(md) >= 4  # header + one per step
    for c in code:
        assert "delfin" in c.source or "smiles_to_xyz" in c.source


def test_each_chemistry_cell_is_playable_python(tmp_path):
    session_id = "sess-export-2"
    _write_trace(tmp_path, session_id, _trace_entries())
    nb = session_export.export_session(
        _fake_session(session_id), trace_root=tmp_path)
    for c in nb.cells:
        if c.cell_type == "code":
            compile(c.source, "<cell>", "exec")  # must be valid Python


def test_internal_agent_tools_are_left_out(tmp_path):
    session_id = "sess-export-3"
    _write_trace(tmp_path, session_id, _trace_entries())
    nb = session_export.export_session(
        _fake_session(session_id), trace_root=tmp_path)
    joined = "\n".join(c.source for c in nb.cells)
    assert "task_create" not in joined
    assert "remember" not in joined


def test_exported_text_is_secret_free(tmp_path):
    session_id = "sess-export-4"
    entries = _trace_entries() + [{
        "ts": 7.0, "tool": "mcp__delfin-ops__smiles_to_xyz",
        "input": json.dumps({"smiles": "C1=CC=CC=C1",
                             "api_key": "sk-super-secret-123456"}),
        "output": "wrote benzene.xyz", "ok": True,
    }]
    _write_trace(tmp_path, session_id, entries)
    nb = session_export.export_session(
        _fake_session(session_id), trace_root=tmp_path)
    joined = "\n".join(c.source for c in nb.cells)
    assert "sk-super-secret" not in joined


def test_missing_trace_still_exports_the_chat(tmp_path):
    session_id = "sess-export-5"
    nb = session_export.export_session(
        _fake_session(session_id), trace_root=tmp_path)
    nbformat.validate(nb)
    assert any("PBE0/def2-SVP" in c.source for c in nb.cells)


def test_write_notebook_creates_the_file(tmp_path):
    session_id = "sess-export-6"
    _write_trace(tmp_path, session_id, _trace_entries())
    out = tmp_path / "out.ipynb"
    session_export.export_session_to_file(
        _fake_session(session_id), out, trace_root=tmp_path)
    nb = nbformat.read(out, as_version=4)
    nbformat.validate(nb)


def test_failed_calls_are_marked_not_hidden(tmp_path):
    session_id = "sess-export-7"
    entries = _trace_entries() + [{
        "ts": 8.0, "tool": "mcp__delfin-ops__extract_energy_table",
        "input": json.dumps({"folders": ["calc/failed"]}),
        "output": "", "ok": False, "error": "no output file",
    }]
    _write_trace(tmp_path, session_id, entries)
    nb = session_export.export_session(
        _fake_session(session_id), trace_root=tmp_path)
    joined = "\n".join(c.source for c in nb.cells)
    assert "no output file" in joined
