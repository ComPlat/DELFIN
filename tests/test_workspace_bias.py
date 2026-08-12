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


def test_every_delfin_only_name_exists_in_the_catalogue():
    """A name that matches nothing filters nothing -- this is the premise
    the disjointness test below rests on, and nothing more.

    It used to carry the name ``test_delfin_only_tool_names_disjoint_from
    _universal`` and the docstring "No DELFIN-only tool should also be
    advertised in the universal set", while asserting ``issubset`` -- that
    the catalogue contains its own entries, which is definitional. Deleting
    the filter it claimed to guard left it green.
    """
    universal_names = {
        t.get("function", {}).get("name") for t in _DOC_TOOLS_OPENAI
    }
    assert _DELFIN_ONLY_TOOL_NAMES.issubset(universal_names)


def _advertised_tool_names(workspace, monkeypatch) -> set:
    """The tool names the model is actually offered for this workspace.

    Read off the request the client builds, not off the catalogue: the
    filter under test runs while the request is assembled.
    """
    import json

    from delfin.agent import mcp_client as _mcp
    from delfin.agent import model_capabilities as _mc

    class _Delta:
        def __init__(self, content=None):
            self.content = content
            self.tool_calls = None
            self.reasoning_content = None

    class _Choice:
        def __init__(self, delta, finish):
            self.index, self.delta, self.finish_reason = 0, delta, finish

    class _Chunk:
        def __init__(self, choices):
            self.choices, self.usage = choices, None

    class _Stream:
        def __init__(self, chunks):
            self._chunks = chunks

        def __iter__(self):
            return iter(self._chunks)

        def close(self):
            pass

    class _NoRegistry:
        def discover_all(self):
            return []

        def discover_resources(self):
            return []

        def discover_prompts(self):
            return []

    caps = _mc.ModelCapabilities(model="m", provider="ollama",
                                 context_window=200_000, supports_tools=True)
    # The doc/calc indexes decide availability independently of the
    # workspace bias; pin them present so this measures the bias alone.
    monkeypatch.setattr(AC._doc_executor, "_ensure_loaded", lambda: True)
    monkeypatch.setattr(AC._doc_executor, "_ensure_calc_loaded", lambda: True)
    monkeypatch.setattr(_mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(_mcp, "get_registry", lambda *a, **k: _NoRegistry())
    monkeypatch.setattr(AC.time, "sleep", lambda *a, **k: None)

    # A model whose profile does NOT set core_tools_only: the weak-model
    # trim drops the chemistry tools for its own reasons and would make
    # the disjointness below true for the wrong one.
    from delfin.agent.model_profiles import get_profile as _get_profile
    model = "qwen2.5-coder:32b"
    assert not _get_profile(model, None).core_tools_only

    client = AC.create_client(backend="api", provider="ollama",
                              model=model, cwd=str(workspace))
    client.set_permissions(KitToolPermissions(workspace=workspace,
                                              mode="acceptEdits"))
    seen: dict = {}

    def _create(**kwargs):
        seen["tools"] = kwargs.get("tools") or []
        return _Stream([_Chunk([_Choice(_Delta(content="ok"), "stop")])])

    client.client.chat.completions.create = _create
    list(client.stream_message("sys", [{"role": "user", "content": "hallo"}],
                               max_tokens=32))
    assert "tools" in seen, "the client sent no request"
    return {t.get("function", {}).get("name") for t in seen["tools"]}


def test_a_generic_workspace_is_offered_no_chemistry_tool(tmp_path,
                                                          monkeypatch):
    """What the filter is for: search_calcs on a workspace with no calcs
    would return nothing, and costs a slot in every request forever."""
    ws = tmp_path / "generic"
    (ws / "src").mkdir(parents=True)
    (ws / "main.py").write_text("print(1)\n")
    assert not _is_delfin_workspace(ws)

    offered = _advertised_tool_names(ws, monkeypatch)
    assert offered, "no tools were advertised at all"
    assert offered & _DELFIN_ONLY_TOOL_NAMES == set(), (
        "a non-DELFIN workspace was offered DELFIN-only tools: "
        f"{sorted(offered & _DELFIN_ONLY_TOOL_NAMES)}")


def test_a_delfin_workspace_is_still_offered_them(tmp_path, monkeypatch):
    """The other half: a filter that drops them everywhere passes the test
    above and breaks the product."""
    ws = _make_delfin_workspace(tmp_path / "delfin_repo")
    assert _is_delfin_workspace(ws)

    offered = _advertised_tool_names(ws, monkeypatch)
    assert _DELFIN_ONLY_TOOL_NAMES <= offered, (
        "a DELFIN workspace lost its chemistry tools: "
        f"{sorted(_DELFIN_ONLY_TOOL_NAMES - offered)}")


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


def test_recommend_task_route_smiles_converter_escalation(monkeypatch):
    """Editing smiles_converter.py is reviewed-tier inside DELFIN, but
    just a plain Python file edit outside.

    The router ALSO escalates adaptively from the recorded task-success
    rate under ~/.delfin. On a machine whose history says coding tasks
    succeed 17% of the time, the outside case escalates too, and the test
    measures this installation's history rather than the workspace bias
    it is about. Neutralise the profile so only the bias is under test.
    """
    monkeypatch.setattr(
        "delfin.agent.provider_profile.load_provider_profile",
        lambda *a, **kw: {})
    msg = "fix the H-atom bug in smiles_converter.py"
    inside = AgentEngine.recommend_task_route(
        msg, current_mode="dashboard", is_delfin_workspace=True,
    )
    outside = AgentEngine.recommend_task_route(
        msg, current_mode="dashboard", is_delfin_workspace=False,
    )
    # Membership was tested against the risk_flags DICT, which iterates its
    # KEYS -- so both sides read as escalated and the check rested entirely
    # on the modes differing. Read the values.
    assert inside["risk_flags"]["reviewed"], "the bias no longer fires"
    assert not outside["risk_flags"]["reviewed"], "it fires outside DELFIN too"
    assert inside["task_class"] == "chemistry"
    assert outside["task_class"] == "coding"
