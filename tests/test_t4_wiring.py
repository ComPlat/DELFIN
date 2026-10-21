"""A misspelled tool name comes back with a grounding hint, not a bare error.

Models called tools that do not exist and got "unknown tool" with nothing to
go on. The executor now asks action_grounding.check with the registered tool
surface before dispatch and attaches its hint to the result. Advisory only:
the call itself is unchanged, nothing is blocked by the hint.
"""
from __future__ import annotations

import json

from delfin.agent import api_client as A


def _run(tmp_path, name, arguments):
    perms = A.KitToolPermissions(workspace=tmp_path, mode="acceptEdits")
    return A._DocToolExecutor().execute(name, arguments, perms)


def test_a_misspelled_tool_gets_the_closest_real_one(tmp_path):
    out = json.loads(_run(tmp_path, "red_file", {"path": "a.txt"}))
    hint = out.get("grounding")
    assert hint, out
    assert hint.get("closest") == "read_file", hint


def test_a_real_tool_gets_no_hint(tmp_path):
    (tmp_path / "a.txt").write_text("x", encoding="utf-8")
    out = _run(tmp_path, "read_file", {"path": str(tmp_path / "a.txt")})
    assert "grounding" not in out


def test_every_registered_doc_tool_is_known(tmp_path, monkeypatch):
    """The hint's tool list reads the nested OpenAI shape, so no registered
    tool can be called 'unknown'. Captures the list; runs one read only."""
    from delfin.agent import action_grounding as G
    seen = {}
    real = G.check

    def spy(tool, args, workspace, known_tools=None):
        seen["known"] = set(known_tools or ())
        return real(tool, args, workspace, known_tools=known_tools)

    monkeypatch.setattr(G, "check", spy)
    (tmp_path / "a.txt").write_text("x", encoding="utf-8")
    _run(tmp_path, "read_file", {"path": str(tmp_path / "a.txt")})
    names = {str((t.get("function") or {}).get("name") or t.get("name") or "")
             for t in A._DOC_TOOLS_OPENAI} - {""}
    assert "read_file" in names
    assert names <= seen["known"], sorted(names - seen["known"])[:10]
