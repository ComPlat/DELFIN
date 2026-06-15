"""Contract consistency guards for the self-describing Baustein layer.

These tests pin the invariant that an adapter's *declared* contract
(:mod:`delfin.tools._spec`) cannot lie about its runtime requirements — the
validator and discovery layers rely on it.  They are deliberately generic so
they cover every adapter that opts into a contract, present and future.
"""

from __future__ import annotations

import pytest

from delfin.tools._registry import list_steps
from delfin.tools._spec import StepContract


_DUMMY_BY_TYPE = {
    "str": "x",
    "int": 0,
    "float": 0.0,
    "bool": False,
    "list": [],
    "dict": {},
    "path": "x.xyz",
}


def _dummy(pspec) -> object:
    return _DUMMY_BY_TYPE.get(pspec.type, "x")


def test_contract_available_on_every_adapter():
    """Every registered adapter yields a valid StepContract, even undeclared ones."""
    for name, adapter in list_steps().items():
        c = adapter.contract()
        assert isinstance(c, StepContract)
        assert c.name == name


def test_produces_geometry_is_folded_into_produces():
    """The primitive produces_geometry flag must appear as a 'geometry' port."""
    for name, adapter in list_steps().items():
        c = adapter.contract()
        if c.produces_geometry:
            assert "geometry" in c.produces, name
        else:
            assert "geometry" not in c.produces, name


def test_declared_required_params_are_enforced():
    """Each declared required param must actually be enforced by validate_params.

    Build a full kwargs dict from all declared params, drop one required param,
    and assert validate_params rejects it.  This stops a contract from claiming
    a param is required when the adapter would silently accept its absence.
    """
    for name, adapter in list_steps().items():
        contract = adapter.contract()
        required = contract.required_params
        if not required:
            continue
        full = {p.name: _dummy(p) for p in contract.params}
        for req in required:
            kwargs = dict(full)
            kwargs.pop(req, None)
            with pytest.raises(ValueError):
                adapter.validate_params(**kwargs)


def test_full_declared_kwargs_pass_validation():
    """Providing every declared param must satisfy validate_params.

    Guards the opposite direction: the contract must declare *enough* params
    that a fully-specified call is accepted (no hidden required key missing
    from the contract).
    """
    for name, adapter in list_steps().items():
        contract = adapter.contract()
        if not contract.params:
            continue
        full = {p.name: _dummy(p) for p in contract.params}
        adapter.validate_params(**full)  # must not raise


def test_enum_defaults_are_within_enum():
    """A param that declares both a default and an enum must be consistent."""
    for name, adapter in list_steps().items():
        for p in adapter.contract().params:
            if p.enum is not None and p.default is not None:
                assert p.default in p.enum, f"{name}.{p.name}"


# --- Discovery layer ------------------------------------------------------


def test_catalog_groups_known_categories():
    from delfin.tools import catalog

    cat = catalog(by="category")
    assert "dft" in cat
    dft_names = {c.name for c in cat["dft"]}
    assert {"orca_sp", "orca_opt", "orca_freq"} <= dft_names


def test_compatible_successors_respect_ports():
    from delfin.tools import compatible_successors

    # orca_opt_freq produces geometry + hessian → imag_fix (needs geometry) fits.
    assert "imag_fix" in compatible_successors("orca_opt_freq")
    # orca_freq does NOT produce geometry → imag_fix is not a pure-capability successor.
    assert "imag_fix" not in compatible_successors("orca_freq")


def test_describe_unknown_step_returns_none():
    from delfin.tools import describe

    assert describe("does_not_exist") is None


# --- Capability wiring ----------------------------------------------------


def test_autowire_gbw_to_moread_for_orca():
    from delfin.tools._wiring import autowire_kwargs

    kw: dict = {}
    autowire_kwargs("orca_sp", kw, {"gbw": "/x/calc.gbw"})
    assert kw["moread"] == "/x/calc.gbw"


def test_autowire_skips_non_orca_without_consumes():
    from delfin.tools._wiring import autowire_kwargs

    kw: dict = {}
    autowire_kwargs("xtb_opt", kw, {"gbw": "/x/calc.gbw"})
    assert "moread" not in kw


def test_autowire_explicit_kwarg_wins():
    from delfin.tools._wiring import autowire_kwargs

    kw = {"moread": "explicit.gbw"}
    autowire_kwargs("orca_sp", kw, {"gbw": "/x/calc.gbw"})
    assert kw["moread"] == "explicit.gbw"


def test_autowire_fires_via_consumes_tag():
    from delfin.tools._wiring import WiringRule, autowire_kwargs

    rules = (WiringRule("hessian", "hess_file", lambda n: False, "hess"),)
    kw: dict = {}
    autowire_kwargs(
        "imag_fix", kw, {"hess": "/x/calc.hess"},
        consumes=frozenset({"hessian"}), rules=rules,
    )
    assert kw["hess_file"] == "/x/calc.hess"


def test_can_autowire_moread_from_gbw():
    from delfin.tools._wiring import can_autowire

    assert can_autowire("moread", frozenset({"gbw"}))
    assert not can_autowire("moread", frozenset())


# --- Build-time validation ------------------------------------------------


def test_prebuilt_templates_validate_clean():
    """The shipped templates must pass static validation given real params.

    Only ERROR diagnostics fail .ok, so this is independent of whether ORCA /
    RDKit happen to be installed on the test host (those are warnings).
    """
    from delfin.tools.templates import classic_opt_freq, redox_potential

    assert redox_potential.validate(
        smiles="CCO", charge=0, method="B3LYP", basis="def2-SVP",
        mult_ox=2, mult_red=2,
    ).ok
    assert classic_opt_freq.validate(
        smiles="CCO", charge=0, method="B3LYP", basis="def2-SVP",
    ).ok


def test_validate_flags_missing_param_and_input():
    from delfin.tools import Pipeline

    p = Pipeline("bad")
    p.add("orca_opt")  # missing required charge + needs an upstream geometry
    rep = p.validate(geometry=False)
    assert not rep.ok
    err = rep.errors[0]
    assert "charge" in err.missing_params
    assert "geometry" in err.missing_inputs


def test_validate_unknown_step_is_error():
    from delfin.tools import Pipeline

    p = Pipeline("x")
    p.add("nope_not_real", charge=0)
    assert not p.validate(geometry=True).ok


def test_validate_geometry_propagates_through_trunk():
    from delfin.tools import Pipeline

    p = Pipeline("chain")
    p.add("smiles_to_xyz", smiles="CCO")   # produces geometry
    p.add("orca_opt", charge=0)            # consumes geometry — satisfied upstream
    rep = p.validate(geometry=False)
    assert rep.ok, rep.summary()


def test_validate_dynamic_compute_relaxes_downstream():
    from delfin.tools import Pipeline

    p = Pipeline("dyn")
    p.add_compute(lambda r, l, w: {"x": 1}, label="calc")  # opaque to the validator
    p.add("orca_sp", charge=0)             # geometry would be missing, but relaxed
    rep = p.validate(geometry=False)
    assert rep.ok, rep.summary()


# --- Serialization round-trip ---------------------------------------------


def _trunk_sig(p):
    return [(s.step_name, s.label, s.kwargs) for s in p._trunk]


def test_roundtrip_prebuilt_templates():
    from delfin.tools import PipelineTemplate
    from delfin.tools.templates import (
        classic_opt_freq, conformer_screening, multi_level_opt, redox_potential,
    )

    for tpl in (classic_opt_freq, redox_potential, conformer_screening, multi_level_opt):
        rebuilt = PipelineTemplate.from_dict(tpl.to_dict())
        assert _trunk_sig(rebuilt) == _trunk_sig(tpl), tpl.name
        assert set(rebuilt._branches) == set(tpl._branches), tpl.name
        for b in tpl._branches:
            assert _trunk_sig(rebuilt._branches[b]) == _trunk_sig(tpl._branches[b]), b


def test_roundtrip_json_concrete_pipeline():
    from delfin.tools import Pipeline

    p = Pipeline("rt", defaults={"method": "B3LYP"})
    p.add("smiles_to_xyz", smiles="CCO")
    p.add("xtb_opt", charge=0)
    p.add("orca_opt", charge=0)
    rebuilt = Pipeline.from_json(p.to_json())
    assert _trunk_sig(rebuilt) == _trunk_sig(p)


def test_roundtrip_if_loop_dsl():
    from delfin.tools import Pipeline
    from delfin.tools._serialize import build_condition, build_until

    p = Pipeline("flow")
    p.add("orca_freq", charge=0)
    p.add_if(build_condition("last.data.n_imaginary > 0"), "imag_fix",
             charge=0, mult=1, solvent="water", metals=[],
             main_basisset="def2-SVP", metal_basisset="def2-TZVP", hess_file="x.hess")
    p.add_loop("orca_opt", until=build_until("result.data.energy_change < 0.0001"),
               max_iter=5, charge=0)
    d = p.to_dict()
    assert d["steps"][1]["type"] == "if"
    assert d["steps"][1]["condition"] == "last.data.n_imaginary > 0"
    assert d["steps"][2]["type"] == "loop"
    rebuilt = Pipeline.from_dict(d)
    assert _trunk_sig(rebuilt) == _trunk_sig(p)


def test_strict_serialize_raises_on_raw_lambda():
    from delfin.tools import Pipeline, PipelineSerializationError

    p = Pipeline("lam")
    p.add_compute(lambda r, l, w: {"x": 1}, label="raw")
    with pytest.raises(PipelineSerializationError):
        p.to_dict(strict=True)
    d = p.to_dict(strict=False)
    assert d.get("_non_serializable_steps")


def test_register_callable_roundtrip():
    from delfin.tools import Pipeline, register_callable

    @register_callable("mytest_compute_xyz")
    def f(results, last, work_dir):
        return {"k": 1}

    p = Pipeline("named")
    p.add_compute(f, label="c")
    rebuilt = Pipeline.from_dict(p.to_dict())  # strict ok — has a ref
    assert rebuilt._trunk[0].compute_fn is f


def test_from_dict_accepts_legacy_flat_kwargs():
    """The old flat YAML form (kwargs as siblings of 'step') must still load."""
    from delfin.tools import Pipeline

    data = {
        "name": "legacy",
        "steps": [
            {"step": "smiles_to_xyz", "smiles": "CCO"},
            {"step": "orca_opt", "charge": 0, "method": "B3LYP"},
        ],
    }
    p = Pipeline.from_dict(data)
    assert p._trunk[0].kwargs == {"smiles": "CCO"}
    assert p._trunk[1].kwargs == {"charge": 0, "method": "B3LYP"}


# --- Result vocabulary ----------------------------------------------------


def test_energy_accessors_normalize_units():
    from delfin.tools import StepResult, StepStatus

    r1 = StepResult("x", StepStatus.SUCCESS, data={"energy_Eh": -1.0})
    assert r1.energy_eh() == -1.0
    assert abs(r1.energy_ev() - (-27.211386245988)) < 1e-9

    r2 = StepResult("y", StepStatus.SUCCESS, data={"energy_eV": -27.211386245988})
    assert abs(r2.energy_eh() - (-1.0)) < 1e-9

    r3 = StepResult("z", StepStatus.SUCCESS, data={})
    assert r3.energy_eh() is None and r3.energy_ev() is None


def test_default_error_kind_is_none():
    from delfin.tools import ErrorKind, StepResult, StepStatus

    # The historical two-arg construction stays valid; error_kind defaults NONE.
    r = StepResult("x", StepStatus.SUCCESS)
    assert r.error_kind is ErrorKind.NONE


# --- OCCUPIER building block (thin wrapper over the legacy engine) ---------


def test_occupier_registered_and_described():
    from delfin.tools import describe, list_steps

    assert "occupier" in list_steps()
    c = describe("occupier")
    assert c.category == "dft"
    assert "orca" in c.requires_binaries
    assert c.produces_geometry is False


def test_occupier_fails_cleanly_without_control(tmp_path):
    """Without a CONTROL.txt the adapter must fail fast (no engine subprocess)."""
    from delfin.tools import run_step
    from delfin.tools._types import ErrorKind, StepStatus

    r = run_step("occupier", work_dir=str(tmp_path))
    assert r.status == StepStatus.FAILED
    assert r.error_kind == ErrorKind.MISSING_INPUT


# --- ESD building block (thin wrapper over the legacy engine) -------------


def test_esd_registered_and_described():
    from delfin.tools import describe, list_steps

    assert "esd" in list_steps()
    c = describe("esd")
    assert c.category == "dft"
    assert "orca" in c.requires_binaries
    assert c.produces_geometry is False
    assert {"completed", "failed", "skipped"} <= {d.name for d in c.data_keys}


def test_esd_fails_cleanly_without_control(tmp_path):
    """Without a CONTROL.txt the adapter must fail fast (no engine subprocess)."""
    from delfin.tools import run_step
    from delfin.tools._types import ErrorKind, StepStatus

    r = run_step("esd", work_dir=str(tmp_path))
    assert r.status == StepStatus.FAILED
    assert r.error_kind == ErrorKind.MISSING_INPUT
