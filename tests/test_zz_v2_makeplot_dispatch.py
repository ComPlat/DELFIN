"""make_plot dispatch probe — protected api_client wiring (package V2, phase 3)
with the operator-required gate assertions.

RED CONTROL for .gate/v2_make_plot.patch. make_plot is NOT yet registered in
api_client (protection: compiled by the operator). These dispatch cases fail
now ("unknown tool") and go green once the operator builds the patch which:
  - registers make_plot in _DOC_TOOLS_OPENAI,
  - adds _DocToolExecutor._execute_make_plot,
  - routes name == "make_plot" in _dispatch,
  - adds make_plot to _GATED_TOOLS + _WRITE_TOOL_NAMES,
  - gates the write (plan mode + _gate_write_path), the read (secret-deny via
    _check_read_access) and wraps the model-facing text with untrusted.wrap.
Assertions are in BOTH directions: fires when it should (card + fenced text +
marker at byte 0) and does NOT fire when it should not (plan mode, .env data
file, glob-denied filename, bad spec, no permissions).
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.api_client import KitToolPermissions


def _perms(ws):
    return KitToolPermissions(workspace=ws, mode="default")


def _perms_plan(ws):
    return KitToolPermissions(workspace=ws, mode="plan")


def _perms_glob(ws, globs):
    return KitToolPermissions(workspace=ws, mode="default",
                              write_allow_globs=globs)


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
    # The untrusted fence itself contains '<src:'; assert instead that no
    # html / iframe / script reaches the model-facing text.
    assert "<iframe" not in obj["text"]
    assert "<script" not in obj["text"]
    # The model-facing text goes through untrusted.wrap (finding 3) while the
    # DELFIN_CARD: marker stays at byte 0 (operator ruling).
    assert "UNTRUSTED EXTERNAL CONTENT" in obj["text"]
    assert "END UNTRUSTED EXTERNAL CONTENT" in obj["text"]


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


def test_plan_mode_is_refused(tmp_path):
    from delfin.agent.api_client import _DocToolExecutor
    out = _DocToolExecutor()._dispatch("make_plot", {"spec": _spec()},
                                       _perms_plan(str(tmp_path)))
    assert not out.startswith("DELFIN_CARD:")
    assert "plan" in out.lower()
    assert not list(tmp_path.glob("*.svg")), "plan mode must not write"


def test_data_file_secret_path_is_refused_by_read_gate(tmp_path):
    from delfin.agent.api_client import _DocToolExecutor
    ws = tmp_path / "ws"
    ws.mkdir()
    dd = ws / ".ssh"
    dd.mkdir()
    (dd / "data.csv").write_text("x,y\n1,2", encoding="utf-8")
    spec = {"kind": "line", "data": {"file": ".ssh/data.csv"},
            "title": "t", "filename": "from_secret.svg"}
    out = _DocToolExecutor()._dispatch("make_plot", {"spec": spec},
                                       _perms(str(ws)))
    assert not out.startswith("DELFIN_CARD:")
    assert "read denied" in out.lower() or "denied" in out.lower()
    assert not (ws / "from_secret.svg").exists()


def test_filename_outside_write_globs_is_refused(tmp_path):
    from delfin.agent.api_client import _DocToolExecutor
    ws = tmp_path / "ws"
    ws.mkdir()
    spec = _spec()
    spec["filename"] = "x.svg"  # a .svg, but globs only allow allowed/*
    out = _DocToolExecutor()._dispatch("make_plot", {"spec": spec},
                                       _perms_glob(str(ws), ("allowed/*",)))
    assert not out.startswith("DELFIN_CARD:")
    assert not (ws / "x.svg").exists(), "refused by the write gate"


def test_existing_allowed_file_written_only_with_gate_ok(tmp_path):
    from delfin.agent.api_client import _DocToolExecutor
    ws = tmp_path / "ws"
    ws.mkdir()
    spec = _spec()
    spec["filename"] = "plot.svg"
    out = _DocToolExecutor()._dispatch("make_plot", {"spec": spec},
                                       _perms(str(ws)))
    assert out.startswith("DELFIN_CARD:")  # gate allowed it
    assert (ws / "plot.svg").exists()
