"""Adversarial tests: solvation workflow grading surface (RUN1, nacht-s25).

These drive the REAL public path (SolvationWorkflow.run -> registry
get("packmol_solvate"/"xtb_opt"/"xtb_solvator") -> execute) with
contract-faithful fake adapters (same StepResult contract the real ones
use), so they probe what a real run would report.

Attack surface:
- T1: nsolv/n_solvent routing — 50 to packmol, 12 to solvator; an
  inversion is physically wrong (box size vs shell size).
- T2: the xTB-relax fallback must leave the "unrelaxed box used" warning
  on the FINAL returned dict (report, not hide).
- T3: cores passthrough — every step gets the cores the caller asked for
  (wave-22 rule: xTB-heavy steps at pal 1; a default other than 1 or a
  silent raise would violate it).
"""
from pathlib import Path

import pytest

from delfin.tools import _registry
from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus


class _RecordingAdapter(StepAdapter):
    """Contract-faithful adapter that records the params it received."""

    name = "unset"

    def __init__(self, status=StepStatus.SUCCESS, geometry="geom.xyz"):
        self.name = self.__class__.NAME
        self.calls = []
        self._status = status
        self._geometry = geometry

    def validate_params(self, **params):
        pass

    def execute(self, work_dir, **params):
        self.calls.append(dict(params))
        return StepResult(step_name=self.name, status=self._status,
                          geometry=Path(self._geometry) if self._geometry else None)


class _Packmol(_RecordingAdapter):
    NAME = "packmol_solvate"


class _XtbOpt(_RecordingAdapter):
    NAME = "xtb_opt"


class _XtbSolvator(_RecordingAdapter):
    NAME = "xtb_solvator"


@pytest.fixture()
def wired(tmp_path):
    """Force discovery, then install recording fakes; restore afterwards."""
    _registry.get("packmol_solvate")
    _registry.get("xtb_opt")
    _registry.get("xtb_solvator")
    pm = _Packmol(geometry="box.xyz")
    opt = _XtbOpt(geometry="relaxed.xyz")
    solv = _XtbSolvator(geometry="solvated.xyz")
    old = {k: _registry._REGISTRY.get(k)
           for k in ("packmol_solvate", "xtb_opt", "xtb_solvator")}
    _registry._REGISTRY.update({"packmol_solvate": pm, "xtb_opt": opt,
                                "xtb_solvator": solv})
    _registry._DISCOVERED = True
    yield pm, opt, solv, tmp_path
    for k, v in old.items():
        _registry._REGISTRY[k] = v


def _workflow():
    from delfin.workflows.contrib.w21_solvation_workflow import SolvationWorkflow
    return SolvationWorkflow()


def test_nsolv_vs_n_solvent_routing(wired):
    """T1: packmol gets n_solvent (box), solvator gets nsolv (shell)."""
    pm, opt, solv, tmp = wired
    wf = _workflow()
    out = wf.run(config={}, geometry="solute.xyz", solvent_xyz="water.xyz",
                 n_solvent=50, nsolv=12, work_dir=str(tmp))
    assert out["ok"], out
    assert pm.calls and pm.calls[0].get("n_solvent") == 50, pm.calls
    assert solv.calls and solv.calls[0].get("n_solvent") == 12, solv.calls


def test_unrelaxed_box_warning_survives_to_result(wired):
    """T2: failed xTB relax must leave the warning on the returned dict."""
    pm, opt, solv, tmp = wired
    opt._status = StepStatus.FAILED
    opt._geometry = None
    wf = _workflow()
    out = wf.run(config={}, geometry="solute.xyz", solvent_xyz="water.xyz",
                 work_dir=str(tmp))
    assert out.get("warnings"), out
    assert any("unrelaxed" in w for w in out["warnings"]), out["warnings"]
    # The solvator must have received the PACKMOL box, not None.
    assert solv.calls and solv.calls[0].get("geometry") == Path("box.xyz")


def test_cores_passthrough_to_all_steps(wired):
    """T3: cores reach every step unchanged (xTB-heavy steps at pal 1)."""
    pm, opt, solv, tmp = wired
    wf = _workflow()
    out = wf.run(config={}, geometry="solute.xyz", solvent_xyz="water.xyz",
                 cores=1, work_dir=str(tmp))
    assert out["ok"], out
    for step in (pm, opt, solv):
        assert step.calls and step.calls[0].get("cores") == 1, step.calls
