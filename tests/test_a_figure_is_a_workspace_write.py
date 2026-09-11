"""A figure the prompt offers is a figure the session may draw.

Asked for a figure of the energies per method, both models called
plot_energy_distribution and were refused -- in default mode every
MCP tool the gate does not know as read-only needs an approval nobody
could give -- and then drew the same PNG through bash and matplotlib,
which the gate let through. The refusal had only moved the write onto
the untyped path. Both named it as the one change (2026-09-11).

A plot tool writes one PNG into the agent workspace and reads
everything else. It is judged as that write: allowed exactly when
write_file to the same path would be, and otherwise handed to the
side-effect gate that asks, as before. Nothing that used to ask now
passes in silence.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from delfin import api
from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions, _doc_executor

_SERVER = Path(api.__file__).resolve().parent / "ops_server" / "server.py"
_PLOTS = sorted(A._MCP_ARTIFACT_TOOL_BASES)


def _gate(name, args, perms):
    return _doc_executor._gate_mcp_tool(name, args, perms)


def test_every_plot_tool_the_server_registers_is_an_artifact_writer():
    registered = set(re.findall(r'mcp\.tool\(name="(plot_[a-z_]+)"\)', _SERVER.read_text(encoding="utf-8")))
    assert registered, "no plot tools registered?"
    assert registered <= A._MCP_ARTIFACT_TOOL_BASES, registered - A._MCP_ARTIFACT_TOOL_BASES


@pytest.mark.parametrize("base", _PLOTS)
def test_a_plot_into_the_workspace_passes_in_default_mode(tmp_path, monkeypatch, base):
    ws = tmp_path / "ws"
    (ws / "agent_workspace").mkdir(parents=True)
    monkeypatch.setattr(api, "_default_plot_dir", lambda: str(ws / "agent_workspace"))
    perms = KitToolPermissions(workspace=ws, mode="default")
    assert _gate(f"mcp__delfin-ops__{base}", {"folders": str(ws)}, perms) is None


def test_an_explicit_output_inside_the_workspace_passes_too(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    perms = KitToolPermissions(workspace=ws, mode="default")
    verdict = _gate("mcp__delfin-ops__plot_energy_distribution",
                    {"folders": str(ws), "output_path": str(ws / "fig.png")}, perms)
    assert verdict is None


def test_a_session_nobody_can_ask_gets_its_figure_inside_its_workspace(tmp_path, monkeypatch):
    """The default plot directory is the dashboard's ~/agent_workspace. A
    CLI or benchmark session cannot write there and has no dialog: it
    used to be refused, and drew the same PNG through bash. Now the tool
    is told to write under the session's own workspace, on the args the
    call is dispatched with."""
    ws = tmp_path / "ws"
    ws.mkdir()
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.setattr(api, "_default_plot_dir", lambda: str(elsewhere))
    perms = KitToolPermissions(workspace=ws, mode="default")
    assert getattr(perms, "confirm_callback", None) is None
    args = {"folders": str(ws)}
    verdict = _gate("mcp__delfin-ops__plot_energy_distribution", args, perms)
    assert verdict is None
    assert Path(args["output_path"]).is_relative_to(ws / "agent_workspace")
    assert args["output_path"].endswith(".png")


def test_a_session_that_can_ask_is_still_asked(tmp_path, monkeypatch):
    ws = tmp_path / "ws"
    ws.mkdir()
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.setattr(api, "_default_plot_dir", lambda: str(elsewhere))
    asked = []
    perms = KitToolPermissions(workspace=ws, mode="default",
                               confirm_callback=lambda name, a, preview: asked.append(name) or False)
    args = {"folders": str(ws)}
    verdict = _gate("mcp__delfin-ops__plot_energy_distribution", args, perms)
    assert verdict is not None                      # the user said no
    assert asked, "a session with a dialog was not asked"
    assert "output_path" not in args                # nothing was redirected behind their back


def test_every_plot_wrapper_accepts_an_output_path():
    from delfin.ops_server import server as ops
    import inspect
    for base in _PLOTS:
        fn = getattr(ops, f"tool_{base}")
        assert "output_path" in inspect.signature(fn).parameters, base
        assert "output_path" in inspect.signature(getattr(api, base)).parameters, base


def test_a_mutating_tool_is_not_an_artifact(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    perms = KitToolPermissions(workspace=ws, mode="default")
    assert _gate("mcp__delfin-ops__delete_calc_folder", {"folder": str(ws / "x")}, perms) is not None


def test_the_target_is_the_tools_own_file_in_the_plot_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(api, "_default_plot_dir", lambda: str(tmp_path))
    assert A._artifact_target_path("plot_uvvis_spectrum", {}) == str(tmp_path / "plot_uvvis_spectrum.png")
    assert A._artifact_target_path("plot_uvvis_spectrum", {"output_path": "/x/y.png"}) == "/x/y.png"
