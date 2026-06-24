"""Tests for the flagship application library (Track 1).

Covers the registered applications, the input-default mechanism that lets a
caller pass only the required inputs (no template-placeholder leaks), the
CONTROL.txt-class composite, and a real end-to-end run of the open-source
xTB thermochemistry workflow.
"""

from __future__ import annotations

import importlib.util
import shutil

import pytest

from delfin.tools import platform

_HAS_XTB = shutil.which("xtb") is not None
_HAS_RDKIT = importlib.util.find_spec("rdkit") is not None

FLAGSHIP = [
    "opt_freq_energy", "redox_potential", "multi_level_energy",
    "xtb_thermochemistry", "conformer_energy", "full_workflow",
]


def test_flagship_applications_registered():
    apps = platform.list_applications()
    for name in FLAGSHIP:
        assert name in apps, f"{name} not registered"


@pytest.mark.parametrize("name", FLAGSHIP)
def test_application_validates_with_only_required_inputs(name):
    # only smiles + charge — every optional input must fill from its default so
    # no template placeholder ({method}/{basis}/{solvent}/…) leaks.
    rep = platform.validate_application(name, geometry=False, smiles="CCO", charge=0)
    assert rep.ok, rep.summary()
    leaks = [m for d in rep.diagnostics for m in d.messages if "placeholder" in m]
    assert leaks == [], leaks


def test_application_fills_input_defaults():
    app = platform.describe_application("xtb_thermochemistry")
    merged = app.with_input_defaults({"smiles": "CO", "charge": 0})
    assert merged["mult"] == 1
    assert merged["method"] == "gfn2"
    assert merged["solvent"] == ""           # default fills → resolves {solvent}


def test_full_workflow_is_control_txt_class():
    app = platform.describe_application("full_workflow")
    steps = [s.get("step") for s in app.spec.get("steps", [])]
    # composes the derived blocks of classic-engine complexity:
    # pre-opt → DFT opt → DFT freq → imaginary cleanup → refine
    for block in ("xtb_opt", "orca_opt", "orca_freq", "imag_fix", "orca_sp"):
        assert block in steps, f"{block} missing from full_workflow"
    assert {o.name for o in app.outputs} == {"energy_Eh", "gibbs_Eh"}


def test_native_applications_are_categorised():
    assert platform.describe_application("xtb_thermochemistry").category == "semiempirical"
    assert platform.describe_application("full_workflow").category == "dft"


@pytest.mark.skipif(not (_HAS_XTB and _HAS_RDKIT), reason="needs xtb + rdkit")
def test_xtb_thermochemistry_runs_end_to_end():
    # fully open-source path: SMILES → native xTB opt → xTB Hessian, no licensed engine
    res = platform.run_application("xtb_thermochemistry", smiles="CO", charge=0, cores=2)
    assert res.ok, res.error
    assert res.outputs["energy_Eh"] < 0
    assert res.outputs["free_energy_Eh"] < 0
    assert res.outputs["n_imaginary"] == 0
    assert res.outputs["zpve_Eh"] > 0
