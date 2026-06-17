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


def test_flagship_applications_registered_and_validate():
    """The flagship template-based apps are registered and validate statically."""
    apps = platform.list_applications()
    assert {"redox_potential", "multi_level_energy", "opt_freq_energy"} <= set(apps)

    redox = platform.describe_application("redox_potential")
    assert redox.required_inputs == ("smiles", "charge")
    out = {o.name: (o.branch, o.step, o.key) for o in redox.outputs}
    assert out["e_oxidation_Eh"] == ("oxidation", "ox_freq", "gibbs_Eh")
    assert out["e_reduction_Eh"] == ("reduction", "red_freq", "gibbs_Eh")
    assert platform.validate_application(
        "redox_potential", smiles="CCO", charge=0, mult_ox=2, mult_red=2).ok
    assert platform.validate_application(
        "multi_level_energy", smiles="CCO", charge=0).ok


def test_application_to_from_dict_roundtrip():
    from delfin.tools import Application

    app = platform.describe_application("opt_freq_energy")
    rebuilt = Application.from_dict(app.to_dict())
    assert rebuilt.name == app.name
    assert rebuilt.required_inputs == app.required_inputs
    assert [o.name for o in rebuilt.outputs] == [o.name for o in app.outputs]
    assert rebuilt.spec == app.spec


def test_user_application_json_is_discovered(tmp_path, monkeypatch):
    """A serialized application dropped in the user dir is auto-loaded."""
    import json

    import delfin.tools._application as A
    from delfin.tools import Application, OutputSpec, Pipeline

    pipe = Pipeline("uapp_pipe")
    pipe.add("fake_energy_step", label="calc")   # registered fake adapter (below)
    app = Application.from_pipeline(
        pipe, name="user_demo_app",
        outputs=(OutputSpec("energy_Eh", step="calc", key="energy_Eh"),),
    )
    (tmp_path / "user_demo_app.json").write_text(json.dumps(app.to_dict()))

    monkeypatch.setenv("DELFIN_APPLICATIONS_DIR", str(tmp_path))
    A._load_user_applications()
    assert "user_demo_app" in platform.list_applications()


def test_run_application_unknown():
    res = platform.run_application("does_not_exist")
    assert res.ok is False and "unknown application" in res.error


# --- agent build-loop toolkit (validate / save / diagnostics / new module) ---


def test_validate_spec_draft_and_application():
    from delfin.tools import Application, OutputSpec, Pipeline

    pipe = Pipeline("vs_pipe")
    pipe.add("fake_energy_step", label="calc")
    assert platform.validate_spec(pipe.to_dict()).ok                    # raw pipeline spec
    app = Application.from_pipeline(
        pipe, name="vs_app",
        outputs=(OutputSpec("energy_Eh", step="calc", key="energy_Eh"),))
    assert platform.validate_spec(app.to_dict()).ok                     # application dict
    # a broken draft is reported, not raised
    assert not platform.validate_spec({"name": "bad", "steps": [{"step": "orca_opt"}]}).ok


def test_save_application_persists_and_registers(tmp_path, monkeypatch):
    from delfin.tools import Application, OutputSpec, Pipeline

    monkeypatch.setenv("DELFIN_APPLICATIONS_DIR", str(tmp_path))
    pipe = Pipeline("sa_pipe")
    pipe.add("fake_energy_step", label="calc")
    app = Application.from_pipeline(
        pipe, name="saved_demo_app",
        outputs=(OutputSpec("energy_Eh", step="calc", key="energy_Eh"),))
    platform.save_application(app.to_dict())
    assert (tmp_path / "saved_demo_app.json").is_file()
    assert "saved_demo_app" in platform.list_applications()


def test_run_diagnostics_unknown_run():
    assert "error" in platform.run_diagnostics("does_not_exist")


def test_new_capability_template_and_user_adapter_discovery(tmp_path, monkeypatch):
    import delfin.tools._registry as R
    from delfin.tools import list_steps

    code = platform.new_capability_template("my_custom_step", category="meta")
    assert "class MyCustomStepAdapter(StepAdapter)" in code
    assert "register(" in code

    monkeypatch.setenv("DELFIN_ADAPTERS_DIR", str(tmp_path))
    (tmp_path / "my_custom_step.py").write_text(code)
    R._load_user_adapters()
    assert "my_custom_step" in list_steps()


def test_register_module_builds_and_integrates(tmp_path, monkeypatch):
    from delfin.tools import list_steps

    monkeypatch.setenv("DELFIN_ADAPTERS_DIR", str(tmp_path))
    code = platform.new_capability_template("built_block", category="analysis")
    res = platform.register_module("built_block", code)
    assert res["ok"] is True
    assert res["provides"] == ["built_block"]
    assert res["capability"]["category"] == "analysis"      # contract is introspectable
    assert "built_block" in list_steps()                    # usable in pipelines now
    assert (tmp_path / "built_block.py").is_file()          # persisted → survives restart


def test_register_module_rejects_broken_source(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_ADAPTERS_DIR", str(tmp_path))
    res = platform.register_module("broken_block", "this is not valid python @#$")
    assert res["ok"] is False
    assert res["diagnostics"] and "failed to load" in res["diagnostics"][0]
    assert not (tmp_path / "broken_block.py").exists()      # poison file cleaned up


def test_register_module_requires_a_registration(tmp_path, monkeypatch):
    monkeypatch.setenv("DELFIN_ADAPTERS_DIR", str(tmp_path))
    res = platform.register_module("noop_block", "VALUE = 1\n")
    assert res["ok"] is False
    assert "registered no new capability" in res["diagnostics"][0]
    assert not (tmp_path / "noop_block.py").exists()


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


def test_tool_inventory_reports_per_tool_status():
    inv = platform.tool_inventory()
    by = {t["name"]: t for t in inv}
    # ORCA is tracked with an installed-status (it's not 'missing' just because it
    # isn't license-installable — that confused a user; the inventory shows it).
    assert "orca" in by
    assert isinstance(by["orca"]["available"], bool)
    assert "orca_sp" in by["orca"]["used_by"]
    assert "morfeus" in by and by["morfeus"]["kind"] == "python"


def test_install_tools_plans_without_executing():
    out = platform.install_tools(run=False)
    assert out["executed"] is False
    assert "plan" in out
    assert set(out["plan"]).issuperset({"auto_binaries", "manual", "python", "installer"})


def test_install_tools_selective_empty_installs_nothing():
    # An empty / non-installable selection executes nothing (deterministic).
    out = platform.install_tools(run=True, select=[])
    assert out["executed"] is False
    assert out["actions"] == []
    # selecting a license-restricted tool also installs nothing
    out2 = platform.install_tools(run=True, select=["orca", "define"])
    assert out2["executed"] is False


def test_open_source_tools_have_specific_install_guidance():
    """Multiwfn/packmol/genarris get a concrete source + hint, not the generic default."""
    from delfin.tools._environment import tool_info

    for tool in ("multiwfn", "packmol", "genarris"):
        info = tool_info(tool, "binary")
        assert info.source, f"{tool} should carry a source URL"
        assert info.install_hint


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
