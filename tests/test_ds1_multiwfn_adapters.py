"""DS1 P2 red test: the unreachable Multiwfn analyses get adapters.

F1: delfin/analysis_tools/multiwfn_wrapper.py has had population_analysis
(:139) and orbital_composition (:181) since the wrapper was written, but no
adapter exposes them — a workflow builder cannot use Multiwfn charges or
orbital compositions because the step registry only knows bond orders.
The step registry is the only discovery path (delfin/tools/_registry.py),
so "reachable" means "registered with a contract".

On the previous commit both cases FAIL: get() raises KeyError.
"""
from __future__ import annotations

import pytest

from delfin.tools._registry import get, list_steps


NEW_STEPS = ("multiwfn_population", "multiwfn_orbital_composition")


@pytest.mark.parametrize("name", NEW_STEPS)
def test_new_multiwfn_step_is_registered(name):
    assert name in list_steps(), (
        f"{name} is not reachable through the step registry; the wrapper "
        "function exists but no workflow can call it"
    )


@pytest.mark.parametrize("name", NEW_STEPS)
def test_new_multiwfn_step_has_a_contract(name):
    adapter = get(name)
    contract = adapter.contract()
    assert contract.required_params, f"{name} declares no required params"
    assert "Multiwfn" in contract.requires_binaries, (
        f"{name} must declare the Multiwfn binary requirement"
    )
    assert adapter.validity_limits, (
        f"{name} must declare validity_limits (wave-21 realism rule)"
    )


@pytest.mark.parametrize("name", NEW_STEPS)
def test_new_multiwfn_step_requires_its_input_file(name):
    with pytest.raises(ValueError):
        get(name).validate_params()


def test_population_execution_reaches_the_real_wrapper(monkeypatch, tmp_path):
    """The adapter calls population_analysis with the given method."""
    import delfin.analysis_tools.multiwfn_wrapper as mfnw

    seen = {}

    def fake(input_file, *, method="hirshfeld"):
        seen["input_file"] = input_file
        seen["method"] = method
        return {"output": "raw", "charges": [], "method": method}

    monkeypatch.setattr(mfnw, "population_analysis", fake)
    result = get("multiwfn_population").execute(
        tmp_path, input_file=str(tmp_path / "x.molden"), method="lowdin"
    )
    assert seen == {"input_file": str(tmp_path / "x.molden"), "method": "lowdin"}
    assert result.data["method"] == "lowdin"


def test_orbital_execution_reaches_the_real_wrapper(monkeypatch, tmp_path):
    import delfin.analysis_tools.multiwfn_wrapper as mfnw

    seen = {}

    def fake(input_file, *, orbital_indices=None):
        seen["input_file"] = input_file
        seen["orbital_indices"] = orbital_indices
        return {"output": "raw"}

    monkeypatch.setattr(mfnw, "orbital_composition", fake)
    result = get("multiwfn_orbital_composition").execute(
        tmp_path, input_file=str(tmp_path / "x.molden")
    )
    assert seen["input_file"] == str(tmp_path / "x.molden")
    assert "output" in result.data
