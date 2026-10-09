"""make_plot dispatch probe — protected api_client wiring (package V2, phase 3).

RED CONTROL for .gate/v2_make_plot.patch. make_plot is NOT yet registered in
api_client (protection: compiled by the operator). These dispatch cases fail
now ("unknown tool") and go green once the operator builds the patch which:
  - registers make_plot in _DOC_TOOLS_OPENAI,
  - adds _DocToolExecutor._execute_make_plot,
  - routes name == "make_plot" in _dispatch.
The executor delegates to delfin.agent.chat_plots.make_plot(spec, out_dir=ws)
and returns the shared DELFIN_CARD: result; errors are JSON {"error": ...}.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.api_client import KitToolPermissions


def _perms(ws):
    return KitToolPermissions(workspace=ws, mode="default")


def _spec():
    return {"kind": "line", "title": "E per method", "x": ["a", "b"],
            "y": [1.0, 2.0]}


def _dispatch(**kwargs):
    from delfin.agent.api_client import _DocToolExecutor
    return _DocToolExecutor()._dispatch("make_plot", kwargs,
                                        _perms(str(__import__("pathlib").Path(
                                            kwargs.pop("_ws", "/tmp")))))


def test_tool_is_registered_in_the_advertised_schema():
    from delfin.agent import api_client as A
    names = [t.get("function", {}).get("name") for t in A._DOC_TOOLS_OPENAI]
    assert "make_plot" in names


def test_dispatch_returns_a_delfin_card_result(tmp_path):
    from delfin.agent.api_client import _DocToolExecutor
    ex = _DocToolExecutor()
    out = ex._dispatch("make_plot", {"spec": _spec()}, _perms(str(tmp_path)))
    assert out.startswith("DELFIN_CARD:")
    obj = json.loads(out[len("DELFIN_CARD:"):])
    assert obj["escape"] == "html"
    assert obj["html"].startswith("<iframe")
    assert 'sandbox="allow-scripts"' in obj["html"]
    assert "allow-same-origin" not in obj["html"]
    assert "<" not in obj["text"]


def test_dispatch_writes_only_inside_the_workspace(tmp_path):
    from delfin.agent.api_client import _DocToolExecutor
    ws = tmp_path / "ws"
    ws.mkdir()
    ex = _DocToolExecutor()
    spec = _spec()
    spec["filename"] = "../escape.svg"
    out = ex._dispatch("make_plot", {"spec": spec}, _perms(str(ws)))
    assert not out.startswith("DELFIN_CARD:")  # error, not a card
    assert "error" in out.lower()
    assert not (tmp_path / "escape.svg").exists()


def test_bad_spec_is_refused_clearly(tmp_path):
    from delfin.agent.api_client import _DocToolExecutor
    ex = _DocToolExecutor()
    out = ex._dispatch("make_plot", {"spec": {"kind": "pie"}},
                       _perms(str(tmp_path)))
    assert not out.startswith("DELFIN_CARD:")
    assert "kind" in out.lower()


def test_requires_permissions():
    from delfin.agent.api_client import _DocToolExecutor
    ex = _DocToolExecutor()
    out = ex._dispatch("make_plot", {"spec": _spec()}, None)
    assert "workspace" in out.lower() or "permission" in out.lower()
