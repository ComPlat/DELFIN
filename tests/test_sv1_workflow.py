"""SV1 phase 2: the contributed explicit-solvation workflow.

Red tests for the workflow under delfin/workflows/contrib/
w21_solvation_workflow.py. Before the workflow existed:
- the workflows registry had no 'explicit_solvation' entry;
- nothing chained packmol_solvate -> xtb_opt -> xtb_solvator;
- a physically wrong argument set (bad solvator method, box too small
  for the solute) had no workflow-level refusal.

The workflow is exercised through the REAL registry discovery
(agent/s20-tm1b auto-discovery, merged into this branch) and the REAL
adapters — no shadow imports.
"""

from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Registry discovery (through the real registry, TM1 auto-discovery)
# ---------------------------------------------------------------------------

def test_explicit_solvation_is_registry_discovered():
    from delfin.workflows.registry import list_workflows

    workflows = list_workflows()
    assert "explicit_solvation" in workflows, (
        "the contrib workflow w21_solvation_workflow must be discovered by "
        "the registry (TM1 auto-discovery on agent/s20-tm1b)")


def test_explicit_solvation_refuses_missing_geometry():
    wf = _workflow()
    with pytest.raises(ValueError, match="geometry"):
        wf.run(config={}, geometry=None, solvent_xyz="water.xyz")


def test_explicit_solvation_refuses_missing_solvent_xyz():
    wf = _workflow()
    with pytest.raises(ValueError, match="solvent_xyz"):
        wf.run(config={}, geometry="solute.xyz", solvent_xyz=None)


def _workflow():
    from delfin.workflows.registry import get as get_workflow

    wf = get_workflow("explicit_solvation")
    assert wf is not None
    return wf


# ---------------------------------------------------------------------------
# The chain runs through the real adapters with stubbed drivers
# ---------------------------------------------------------------------------

def _write_solute(path: Path) -> Path:
    path.write_text("1\nx\nO 0.0 0.0 0.0\n", encoding="utf-8")
    return path


def _write_solvent(path: Path) -> Path:
    path.write_text("3\nw\nO 0.0 0.0 0.0\nH 0.8 0.0 0.0\nH -0.4 0.7 0.0\n",
                    encoding="utf-8")
    return path


def _stub_run_orca(monkeypatch, capture: dict):
    """A fake ORCA that 'produces' the SOLVATOR output file and records
    the input text — the real deliverable of the solvator step."""

    def fake(inp, out, working_dir=None):
        capture.setdefault("inp", []).append(
            Path(inp).read_text(encoding="utf-8"))
        # The real deliverable lands next to the input, in the run dir
        Path(inp).with_name("XTB_SOLVATOR.solvator.xyz").write_text(
            "4\nsolvated\nO 0 0 0\nH 1 0 0\nH -0.4 0.7 0\nH 0 1 0\n",
            encoding="utf-8")

    monkeypatch.setattr("delfin.orca.run_orca", fake)


def _stub_solvate(monkeypatch, out_xyz: Path):
    """A fake packmol that 'produces' the solvated box (the binary is
    not installed in this cage; the adapter, not packmol, is the test
    subject)."""
    from delfin.analysis_tools import packmol_wrapper

    def fake(**kwargs):
        out_xyz.write_text(
            "4\nbox\nO 0 0 0\nH 1 0 0\nH -0.4 0.7 0\nO 3 3 3\n",
            encoding="utf-8")
        return out_xyz

    monkeypatch.setattr(packmol_wrapper, "solvate", fake)


def _stub_xtb_opt(monkeypatch, out_xyz: Path):
    """A fake xtb that 'produces' the relaxed box."""
    from delfin.tools.adapters import xtb

    def fake(self, work_dir, *, geometry=None, cores=1, **kwargs):
        from delfin.tools._types import StepResult, StepStatus
        out_xyz.write_text("4\nrelaxed\nO 0 0 0\nH 1 0 0\nH -0.4 0.7 0\n"
                           "O 3 3 3\n", encoding="utf-8")
        return StepResult(step_name="xtb_opt", status=StepStatus.SUCCESS,
                          geometry=out_xyz, work_dir=str(work_dir))

    monkeypatch.setattr(xtb.XtbOptAdapter, "execute", fake)


def test_workflow_chain_runs_and_refuses_bad_solvator_method(
        tmp_path, monkeypatch):
    """A DFT method for SOLVATOR must fail at the WORKFLOW level with the
    adapter's physics message — the ORCA input writer must never run."""
    wf = _workflow()
    solute = _write_solute(tmp_path / "solute.xyz")
    solvent = _write_solvent(tmp_path / "water.xyz")
    capture: dict = {}
    _stub_run_orca(monkeypatch, capture)
    _stub_solvate(monkeypatch, tmp_path / "box.xyz")
    _stub_xtb_opt(monkeypatch, tmp_path / "relaxed.xyz")

    result = wf.run(
        config={}, geometry=str(solute), solvent_xyz=str(solvent),
        solvator_method="B3LYP", work_dir=str(tmp_path / "run"))
    assert result["ok"] is False
    assert "xTB Hamiltonian" in (result["errors"][0]), result
    # The solvator's ORCA input must never have been written
    assert all("B3LYP" not in text for text in capture.get("inp", []))


def test_workflow_chain_succeeds_and_reports_levels(
        tmp_path, monkeypatch):
    """The does-not-fire direction: with valid parameters the three-step
    chain runs through the REAL adapters (stubs only at the binary
    boundary) and the result states which level computed what — an
    xTB-relaxed structure must never be mistaken for the shell level."""
    wf = _workflow()
    solute = _write_solute(tmp_path / "solute.xyz")
    solvent = _write_solvent(tmp_path / "water.xyz")
    capture: dict = {}
    _stub_run_orca(monkeypatch, capture)
    _stub_solvate(monkeypatch, tmp_path / "box.xyz")
    _stub_xtb_opt(monkeypatch, tmp_path / "relaxed.xyz")

    result = wf.run(
        config={}, geometry=str(solute), solvent_xyz=str(solvent),
        nsolv=12, solvator_method="GFN2-xTB",
        work_dir=str(tmp_path / "run"))
    assert result["ok"] is True, result
    assert set(result["steps"]) == {
        "packmol_box", "xtb_relax", "xtb_solvator"}, result["steps"]
    # Charge and mult must be inherited by the solvator step (spin is a
    # property of the solute, not reset between steps)
    assert "*xyz 0 1" in capture["inp"][-1]
    assert "%SOLVATOR NSOLV 12" in capture["inp"][-1]
    assert result["levels"]["box_relaxation"].startswith("GFN-xTB")
    assert result["levels"]["explicit_shell"] == "GFN2-xTB"
    assert result["solvated_geometry"].endswith("XTB_SOLVATOR.solvator.xyz")


def test_run_cli_prints_why_it_failed(tmp_path, monkeypatch, capsys):
    """F1 (nacht-s19): exit 1 with NO message hides the cause (missing
    binary, wrong path) — run_cli must print the reason to stderr.
    Cross-review finding (nacht-s3): without a chdir the default
    work_dir lands in the checkout and trips the conftest leak guard."""
    monkeypatch.chdir(tmp_path)
    wf = _workflow()
    rc = wf.run_cli(["missing.xyz", "missing_solvent.xyz"])
    out = capsys.readouterr()
    assert rc == 1
    assert out.err.strip(), "run_cli failed silently"
    assert "explicit_solvation" in out.err or "geometry" in out.err


def test_run_cli_success_path(tmp_path, monkeypatch, capsys):
    """The CLI exit 0 on a run whose steps succeed (stubs at the
    binary boundary, as in the chain test)."""
    capture: dict = {}
    _stub_run_orca(monkeypatch, capture)
    _stub_solvate(monkeypatch, tmp_path / "box.xyz")
    _stub_xtb_opt(monkeypatch, tmp_path / "relaxed.xyz")
    solute = _write_solute(tmp_path / "solute.xyz")
    solvent = _write_solvent(tmp_path / "water.xyz")
    wf = _workflow()
    rc = wf.run_cli([
        str(solute), str(solvent), "--nsolv", "12",
        "--work-dir", str(tmp_path / "run")])
    assert rc == 0, capsys.readouterr().err
