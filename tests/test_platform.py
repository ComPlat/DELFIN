"""Tests for the platform facade: capabilities, applications, environment."""

from __future__ import annotations

import time

import pytest

from delfin.tools import platform
from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._types import StepResult, StepStatus


# --- Capabilities ---------------------------------------------------------


def test_list_and_describe_capabilities():
    caps = platform.list_capabilities()
    assert "orca_sp" in caps and "occupier" in caps
    c = platform.describe_capability("orca_sp")
    assert c is not None and c.name == "orca_sp"
    assert platform.describe_capability("nope_xyz") is None


# --- Applications ---------------------------------------------------------


def test_builtin_application_registered_and_described():
    assert "opt_freq_energy" in platform.list_applications()
    app = platform.describe_application("opt_freq_energy")
    assert app is not None
    assert app.required_inputs == ("smiles", "charge")
    assert {o.name for o in app.outputs} == {"energy_Eh", "gibbs_Eh"}


def test_application_to_from_dict_roundtrip():
    from delfin.tools import Application

    app = platform.describe_application("opt_freq_energy")
    rebuilt = Application.from_dict(app.to_dict())
    assert rebuilt.name == app.name
    assert rebuilt.required_inputs == app.required_inputs
    assert [o.name for o in rebuilt.outputs] == [o.name for o in app.outputs]
    assert rebuilt.spec == app.spec


def test_run_application_unknown():
    res = platform.run_application("does_not_exist")
    assert res.ok is False and "unknown application" in res.error


def test_run_application_missing_required_input():
    # opt_freq_energy requires smiles + charge; omit them → fail fast, no run.
    res = platform.run_application("opt_freq_energy", method="B3LYP")
    assert res.ok is False
    assert "missing required input" in res.error
    assert res.pipeline_result is None


# --- End-to-end run with a fake adapter (no real QM) ----------------------


class _FakeEnergyAdapter(StepAdapter):
    name = "fake_energy_step"
    description = "fake energy producer for platform tests"
    produces_geometry = False

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, time.monotonic(),
            data={"energy_Eh": -1.5},
        )


register(_FakeEnergyAdapter())


def test_run_application_extracts_named_outputs(tmp_path):
    from delfin.tools import Application, OutputSpec, Pipeline

    pipe = Pipeline("fake_app_pipe")
    pipe.add("fake_energy_step", label="calc")
    app = Application.from_pipeline(
        pipe,
        name="fake_app",
        outputs=(OutputSpec("energy_Eh", step="calc", key="energy_Eh", unit="Eh"),),
    )
    platform.register_application(app)

    res = platform.run_application("fake_app", work_dir=tmp_path)
    assert res.ok is True
    assert res.outputs == {"energy_Eh": -1.5}


# --- Environment (probe + install policy) ---------------------------------


def test_probe_covers_every_capability():
    readiness = platform.probe()
    names = {r.step_name for r in readiness}
    assert {"orca_sp", "smiles_to_xyz"} <= names


def test_install_policy_orca_and_turbomole_are_manual():
    from delfin.tools._environment import tool_info

    assert tool_info("orca", "binary").policy == "manual"
    assert tool_info("turbomole", "binary").policy == "manual"
    assert tool_info("define", "binary").policy == "manual"
    # open-source QM binaries are auto-installable
    assert tool_info("xtb", "binary").policy == "auto"


def test_orca_is_never_scheduled_for_auto_install():
    """Legal guarantee: license-restricted tools never land in auto_binaries."""
    plan = platform.install_plan()
    assert "orca" not in plan["auto_binaries"]
    for tm in ("turbomole", "define", "x2t", "ridft"):
        assert tm not in plan["auto_binaries"]
    # every manual entry really is policy=manual and carries a hint
    from delfin.tools._environment import tool_info
    for entry in plan["manual"]:
        assert tool_info(entry["tool"], "binary").policy == "manual"
        assert entry["hint"]


def test_missing_orca_requirement_is_annotated_manual():
    """Wherever ORCA shows up as missing, it is flagged manual with a hint."""
    for cr in platform.probe():
        for r in cr.missing:
            if r.name.lower() == "orca":
                assert r.policy == "manual"
                assert "orca" in r.install_hint.lower()


def test_install_tools_plans_without_executing():
    out = platform.install_tools(run=False)
    assert out["executed"] is False
    assert "plan" in out
    assert set(out["plan"]).issuperset({"auto_binaries", "manual", "python", "installer"})


# --- Key vocabulary -------------------------------------------------------


def test_central_keys_and_allowed_values():
    keys = platform.list_keys()
    assert {"functional", "basis", "solvent", "relativity", "dispersion"} <= set(keys)
    func = platform.describe_key("functional")
    assert func is not None and func.enum and "PBE0" in func.enum
    assert "ZORA" in platform.describe_key("relativity").enum
    assert platform.describe_key("does_not_exist") is None


def test_key_helper_builds_paramspec_with_enum():
    from delfin.tools import key

    p = key("functional", required=True)
    assert p.name == "functional" and p.required is True
    assert p.enum and "B3LYP" in p.enum
    # unknown key still yields a usable ParamSpec
    q = key("totally_unknown")
    assert q.name == "totally_unknown" and q.enum is None


# --- Manifest -------------------------------------------------------------


def test_manifest_structure():
    m = platform.manifest()
    assert m["delfin_tools_manifest"] == "1"
    assert {"capabilities", "applications", "keys", "schemas"} <= set(m)
    cap_names = {c["name"] for c in m["capabilities"]}
    assert {"orca_sp", "occupier", "esd"} <= cap_names
    key_names = {k["name"] for k in m["keys"]}
    assert "functional" in key_names
    func = next(k for k in m["keys"] if k["name"] == "functional")
    assert "PBE0" in func["enum"] and func["enum_source"]
    assert m["schemas"]["pipeline_spec"]["title"] == "DelfinPipelineSpec"
    assert m["schemas"]["application"]["title"] == "DelfinApplication"


def test_manifest_json_is_valid_json():
    import json

    data = json.loads(platform.manifest_json())
    assert data["delfin_tools_manifest"] == "1"


# --- MCP request handlers (no mcp dependency needed) ----------------------


def test_mcp_handlers_return_valid_json():
    import json

    from delfin.tools import mcp_server as srv

    # read-only handlers
    assert "orca_sp" in json.loads(srv.h_list_capabilities())
    assert json.loads(srv.h_describe_capability("orca_sp"))["name"] == "orca_sp"
    assert "error" in json.loads(srv.h_describe_capability("nope_xyz"))
    assert any(k["name"] == "functional" for k in json.loads(srv.h_list_keys()))
    assert json.loads(srv.h_describe_key("functional"))["name"] == "functional"
    assert "opt_freq_energy" in json.loads(srv.h_list_applications())
    assert json.loads(srv.h_describe_application("opt_freq_energy"))["name"] == "opt_freq_energy"
    assert "capabilities" in json.loads(srv.h_get_manifest())

    # validate handler
    v = json.loads(srv.h_validate_application("opt_freq_energy",
                                              {"smiles": "CCO", "charge": 0}))
    assert "ok" in v and "diagnostics" in v

    # run handler on an unknown app fails fast with structured error
    r = json.loads(srv.h_run_application("does_not_exist"))
    assert r["ok"] is False and r["error"]

    # environment handlers
    assert isinstance(json.loads(srv.h_probe()), list)
    plan = json.loads(srv.h_install_plan())
    assert "orca" not in plan["auto_binaries"]
