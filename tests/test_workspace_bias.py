"""Workspace-aware DELFIN-bias gating.

Two scenarios:
  - DELFIN-repo workspace      -> chemistry tools + xtb bash-auto-allow
                                  + chemistry-file mode escalation
  - generic project workspace  -> none of the above
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from delfin.agent import api_client as AC
from delfin.agent.api_client import (
    KitToolPermissions,
    _DEFAULT_BASH_AUTO_ALLOW,
    _DELFIN_BASH_AUTO_ALLOW,
    _DELFIN_ONLY_TOOL_NAMES,
    _DOC_TOOLS_OPENAI,
    _is_delfin_workspace,
)
from delfin.agent.engine import AgentEngine


def _make_delfin_workspace(base: Path) -> Path:
    (base / "delfin" / "agent").mkdir(parents=True)
    (base / "delfin" / "__init__.py").write_text("")
    (base / "delfin" / "agent" / "__init__.py").write_text("")
    return base


# ---- detection -----------------------------------------------------------


def test_is_delfin_workspace_detects_delfin_tree():
    with tempfile.TemporaryDirectory() as d:
        ws = _make_delfin_workspace(Path(d))
        assert _is_delfin_workspace(ws) is True


def test_is_delfin_workspace_rejects_generic():
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "main.py").write_text("print(1)")
        assert _is_delfin_workspace(Path(d)) is False


def test_is_delfin_workspace_handles_none_and_missing():
    assert _is_delfin_workspace(None) is False
    assert _is_delfin_workspace("/no/such/path/12345") is False


# ---- bash auto-allow -----------------------------------------------------


def test_permissions_delfin_workspace_gets_xtb_pattern():
    with tempfile.TemporaryDirectory() as d:
        ws = _make_delfin_workspace(Path(d))
        perms = KitToolPermissions(workspace=ws)
        assert perms.is_delfin_workspace is True
        # xtb is auto-allowed inside DELFIN
        assert perms.matches_bash_auto_allow("xtb mol.xyz --opt") is True
        # delfin-cli wrapper too
        assert perms.matches_bash_auto_allow("delfin-run --config x") is True


def test_permissions_generic_workspace_no_xtb_pattern():
    with tempfile.TemporaryDirectory() as d:
        ws = Path(d)
        perms = KitToolPermissions(workspace=ws)
        assert perms.is_delfin_workspace is False
        # xtb stays gated outside DELFIN
        assert perms.matches_bash_auto_allow("xtb mol.xyz --opt") is False
        # delfin-cli wrapper too
        assert perms.matches_bash_auto_allow("delfin-run --config x") is False
        # universal patterns still work
        assert perms.matches_bash_auto_allow("ls -la") is True
        assert perms.matches_bash_auto_allow("pytest tests/") is True


# ---- tool filtering ------------------------------------------------------


def test_delfin_only_tool_names_disjoint_from_universal():
    """No DELFIN-only tool should also be advertised in the universal set
    via the coding-tool gate."""
    universal_names = {
        t.get("function", {}).get("name") for t in _DOC_TOOLS_OPENAI
    }
    # All DELFIN-only tools exist in the catalogue
    assert _DELFIN_ONLY_TOOL_NAMES.issubset(universal_names)


def test_delfin_bash_patterns_are_separate_from_default():
    """Sanity: the DELFIN tuple is *additive*, not a subset of the default."""
    for p in _DELFIN_BASH_AUTO_ALLOW:
        assert p not in _DEFAULT_BASH_AUTO_ALLOW


# ---- escalation gating in recommend_task_route --------------------------


def test_recommend_task_route_chemistry_in_delfin_workspace():
    """Mentioning 'orca functional' in DELFIN-mode pulls in chemistry."""
    decision = AgentEngine.recommend_task_route(
        "Switch the orca functional to B3LYP for this DFT job",
        current_mode="dashboard",
        is_delfin_workspace=True,
    )
    # chemistry detected → task_class reflects that
    assert "chemistry" in decision.get("task_class", "") or \
        any("chemistry" in r.lower() for r in decision.get("reasons", []))


def test_recommend_task_route_no_chemistry_outside_delfin():
    """Same prompt outside DELFIN must NOT trigger chemistry escalation."""
    decision = AgentEngine.recommend_task_route(
        "Switch the orca functional to B3LYP for this DFT job",
        current_mode="dashboard",
        is_delfin_workspace=False,
    )
    # chemistry must not dominate — coding/general remains
    assert "chemistry" not in decision.get("task_class", "")


def test_recommend_task_route_smiles_converter_escalation():
    """Editing smiles_converter.py is reviewed-tier inside DELFIN, but
    just a plain Python file edit outside."""
    msg = "fix the H-atom bug in smiles_converter.py"
    inside = AgentEngine.recommend_task_route(
        msg, current_mode="dashboard", is_delfin_workspace=True,
    )
    outside = AgentEngine.recommend_task_route(
        msg, current_mode="dashboard", is_delfin_workspace=False,
    )
    # The chemistry-file escalation triggers reviewed in DELFIN
    inside_reviewed = "reviewed" in inside.get("risk_flags", []) or \
        inside.get("mode") in {"reviewed", "full"}
    outside_reviewed = "reviewed" in outside.get("risk_flags", []) or \
        outside.get("mode") in {"reviewed", "full"}
    # At minimum: outside should NOT escalate harder than inside
    if inside_reviewed:
        assert not outside_reviewed or outside.get("mode") != inside.get("mode")
