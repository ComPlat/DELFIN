"""Agent MCP tools for MLP evaluation (wave19 m4).

The handlers live in delfin.api, are wrapped in delfin.ops_server.server,
and are gated like every other tool: only the read-only one
(``mlp_list_backends``) is on _MCP_READONLY_TOOL_BASES; the run tools
(single point, optimise, rank, pre-optimise) start denied by default.

Every run test uses a FAKE calculator layer: the handlers call through
``delfin.mlp_tools.calculators``, which these tests monkeypatch.  No torch,
no GPU, no backend package is needed to test DELFIN's own contract.
"""

from __future__ import annotations

import pytest


# ---------------------------------------------------------------- fakes ---
class _FakeAtoms:
    """Minimal stand-in for ase.Atoms as the calculators layer sees it."""

    def __init__(self, symbols=("H", "O"), energy=-10.0):
        self.symbols = list(symbols)
        self.energy = energy
        self.calc = None

    def get_chemical_symbols(self):
        return self.symbols

    def __len__(self):
        return len(self.symbols)


@pytest.fixture
def fake_calcs(monkeypatch):
    """Replace the calculators layer with fakes; record the calls."""
    import delfin.mlp_tools.calculators as calcs

    calls = []

    def fake_single_point(atoms, backend="ani2x", *, device="cpu",
                          charge=0, mult=1):
        calls.append(("single_point", backend, device, charge, mult))
        return {"energy": -10.5, "forces": [[0.0, 0.0, 0.1], [0.0, 0.0, -0.1],
                                            [0.0, 0.1, 0.0]],
                "backend": backend, "n_atoms": len(atoms)}

    def fake_optimize(atoms, backend="ani2x", *, device="cpu", charge=0,
                      mult=1, fmax=0.05, steps=200, optimizer="LBFGS"):
        calls.append(("optimize", backend, fmax, steps, optimizer))
        atoms.positions_moved = True
        return {"converged": True, "energy": -10.9, "n_steps": 7,
                "backend": backend}

    def fake_rank(structures, backend="ani2x", *, device="cpu", charge=0,
                  mult=1):
        calls.append(("rank", backend, len(structures)))
        return [(1, -10.2), (0, -9.8)]

    monkeypatch.setattr(calcs, "single_point", fake_single_point)
    monkeypatch.setattr(calcs, "optimize_geometry", fake_optimize)
    monkeypatch.setattr(calcs, "rank_structures", fake_rank)
    return calls


XYZ = "3\nwater\nO 0.0 0.0 0.0\nH 0.9 0.0 0.0\nH -0.3 0.85 0.0\n"

ENSEMBLE = (
    "3\nconf a\nO 0.0 0.0 0.0\nH 0.9 0.0 0.0\nH -0.3 0.85 0.0\n"
    "3\nconf b\nO 0.0 0.0 0.0\nH 0.85 0.1 0.0\nH -0.3 0.8 0.1\n"
)


# ------------------------------------------------------- mlp_list_backends ---
def test_list_backends_reports_every_registered_backend():
    from delfin.api import mlp_list_backends

    out = mlp_list_backends()
    names = [b["name"].lower().replace("-", "").replace("_", "") for b in out["backends"]]
    # every backend the calculator factory registers is named, even when
    # its package is missing on this machine
    for expected in ("ani2x", "aimnet2", "maceoff", "chgnet", "m3gnet",
                     "schnetpack", "nequip", "alignn"):
        assert expected in names, (expected, names)
    for entry in out["backends"]:
        assert isinstance(entry["installed"], bool)
        assert "elements" in entry
        if not entry["installed"]:
            assert entry.get("install_hint"), entry


def test_list_backends_is_read_only_and_on_the_gate_allowlist():
    from delfin.agent.api_client import _MCP_READONLY_TOOL_BASES

    assert "mlp_list_backends" in _MCP_READONLY_TOOL_BASES


def test_run_tools_are_not_on_the_gate_allowlist():
    from delfin.agent.api_client import _MCP_READONLY_TOOL_BASES

    for tool in ("mlp_single_point", "mlp_optimize", "mlp_rank",
                 "mlp_preoptimize"):
        assert tool not in _MCP_READONLY_TOOL_BASES, tool


# ------------------------------------------------------ mlp_single_point ---
def test_single_point_runs_and_returns_energy_and_forces(fake_calcs):
    from delfin.api import mlp_single_point

    out = mlp_single_point(XYZ, backend="mace")
    assert out["ok"] is True
    assert out["energy_eV"] == -10.5
    assert out["backend"] == "mace"          # the requested spelling
    assert out["n_atoms"] == 3
    assert len(out["forces"]) == 3
    assert fake_calcs == [("single_point", "mace", "cpu", 0, 1)]


def test_single_point_passes_charge_and_mult(fake_calcs):
    from delfin.api import mlp_single_point

    out = mlp_single_point(XYZ, backend="aimnet2", charge=1, mult=2)
    assert out["ok"] is True
    assert fake_calcs == [("single_point", "aimnet2", "cpu", 1, 2)]


def test_single_point_unknown_backend_is_an_error_not_a_crash(fake_calcs):
    from delfin.api import mlp_single_point

    out = mlp_single_point(XYZ, backend="not_a_backend")
    assert out["ok"] is False
    assert "ani2x" in out["error"]           # names the accepted backends
    assert fake_calcs == []                  # nothing ran


def test_single_point_missing_package_reports_clearly(monkeypatch):
    from delfin.api import mlp_single_point
    import delfin.mlp_tools.calculators as calcs

    def boom(*a, **k):
        raise ModuleNotFoundError("No module named 'torchani'")

    monkeypatch.setattr(calcs, "create_calculator", boom)
    out = mlp_single_point(XYZ, backend="ani2x")
    assert out["ok"] is False
    assert out["available"] is False
    assert "torchani" in out["reason"]
    assert out["install_hint"]               # never a bare traceback


# ---------------------------------------------------------- mlp_optimize ---
def test_optimize_returns_geometry_and_convergence(fake_calcs):
    from delfin.api import mlp_optimize

    out = mlp_optimize(XYZ, backend="ani2x", fmax=0.02, steps=100)
    assert out["ok"] is True
    assert out["converged"] is True
    assert out["energy_eV"] == -10.9
    assert out["n_steps"] == 7
    assert out["xyz"].startswith("3\n")      # optimised geometry as xyz text
    assert fake_calcs == [("optimize", "ani2x", 0.02, 100, "LBFGS")]


# -------------------------------------------------------------- mlp_rank ---
def test_rank_orders_structures_lowest_first(fake_calcs):
    from delfin.api import mlp_rank

    out = mlp_rank(ENSEMBLE, backend="ani2x")
    assert out["ok"] is True
    assert out["n_structures"] == 2
    energies = [e["energy_eV"] for e in out["ranking"]]
    assert energies == sorted(energies)
    assert out["ranking"][0]["index"] == 1   # the lower energy first
    assert fake_calcs == [("rank", "ani2x", 2)]


def test_rank_rejects_a_single_structure_ensemble(fake_calcs):
    from delfin.api import mlp_rank

    out = mlp_rank(XYZ, backend="ani2x")
    assert out["ok"] is False
    assert "ensemble" in out["error"].lower()
    assert fake_calcs == []


# -------------------------------------------------------- mlp_preoptimize ---
def test_preoptimize_is_optimize_with_preopt_defaults_and_a_dft_note(fake_calcs):
    from delfin.api import mlp_preoptimize

    out = mlp_preoptimize(XYZ, backend="mace")
    assert out["ok"] is True
    assert out["converged"] is True
    assert out["xyz"].startswith("3\n")
    # the contract: an MLP geometry is a STARTING GUESS, never a final answer
    assert "DFT" in out["note"]
    assert fake_calcs == [("optimize", "mace", out["fmax"], out["steps"],
                           "LBFGS")]


def test_preoptimize_missing_package_reports_clearly(monkeypatch):
    from delfin.api import mlp_preoptimize
    import delfin.mlp_tools.calculators as calcs

    def boom(*a, **k):
        raise ModuleNotFoundError("No module named 'mace'")

    monkeypatch.setattr(calcs, "create_calculator", boom)
    out = mlp_preoptimize(XYZ, backend="mace")
    assert out["ok"] is False
    assert out["available"] is False
    assert out["install_hint"]


# ---------------------------------------------------- server registration ---
def test_every_mlp_tool_is_wrapped_in_the_ops_server():
    import re
    from pathlib import Path

    src = (Path(__file__).resolve().parents[1] / "delfin" / "ops_server"
           / "server.py").read_text(encoding="utf-8")
    for name in ("mlp_list_backends", "mlp_single_point", "mlp_optimize",
                 "mlp_rank", "mlp_preoptimize"):
        assert f'mcp.tool(name="{name}")' in src, name
        assert f"def tool_{name}(" in src, name
        # handlers never leak the _safe wrapper's exceptions raw: the run
        # tools return ok=False dicts instead
    assert not re.search(r"mcp\.tool\(name=\"mlp_.*\"\)\(.*raise", src)


def test_every_mlp_tool_is_in_the_tool_catalog():
    from delfin.api import _TOOL_CATALOG

    names = {t["name"] for t in _TOOL_CATALOG}
    for name in ("mlp_list_backends", "mlp_single_point", "mlp_optimize",
                 "mlp_rank", "mlp_preoptimize"):
        assert name in names, name


def test_api_handlers_exist_and_are_documented():
    import delfin.api as api

    for name in ("mlp_list_backends", "mlp_single_point", "mlp_optimize",
                 "mlp_rank", "mlp_preoptimize"):
        fn = getattr(api, name)
        assert fn.__doc__, name
