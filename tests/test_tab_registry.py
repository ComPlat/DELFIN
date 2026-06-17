"""Tests for the dashboard tab registry + auto-generated application form."""

from __future__ import annotations

import pytest

from delfin.dashboard import tab_application as TA
from delfin.dashboard import tab_registry as R


# --- registry -------------------------------------------------------------


def test_register_and_build_specs_defensively():
    R.register_tab("t_ok", "OK", lambda ctx: ("W", {"p": 1}), order=120)

    def _boom(ctx):
        raise RuntimeError("kaboom")

    R.register_tab("t_boom", "Boom", _boom, order=130)

    specs = {s["id"]: s for s in R.registered_tab_specs(ctx=None)}
    assert specs["t_ok"]["widget"] == "W" and specs["t_ok"]["available"] is True
    # a failing factory is included but marked unavailable — never raises
    assert specs["t_boom"]["widget"] is None
    assert specs["t_boom"]["available"] is False
    assert "kaboom" in specs["t_boom"]["reason"]


def test_registered_tabs_are_ordered():
    R.register_tab("z_late", "Z", lambda ctx: "w", order=9999)
    R.register_tab("a_early", "A", lambda ctx: "w", order=1)
    orders = [t.order for t in R.registered_tabs()]
    assert orders == sorted(orders)


def test_pipelines_tab_is_auto_discovered():
    # importing nothing extra: discovery imports tab_application which self-registers
    ids = {t.id for t in R.registered_tabs()}
    titles = {t.title for t in R.registered_tabs()}
    assert "pipelines" in ids and "Pipelines" in titles


# --- application form logic (no ipywidgets needed) ------------------------


def test_form_field_specs_maps_paramspecs():
    from delfin.tools import ParamSpec

    inputs = (
        ParamSpec("smiles", "str", required=True),
        ParamSpec("charge", "int", required=True),
        ParamSpec("method", "str", default="B3LYP", enum=("B3LYP", "PBE0")),
        ParamSpec("flag", "bool", default=False),
        ParamSpec("fmax", "float", default=0.05),
    )
    fields = {f["name"]: f for f in TA.form_field_specs(inputs)}
    assert fields["smiles"]["kind"] == "text" and fields["smiles"]["required"]
    assert fields["charge"]["kind"] == "int"
    assert fields["method"]["kind"] == "enum" and "PBE0" in fields["method"]["options"]
    assert fields["flag"]["kind"] == "bool"
    assert fields["fmax"]["kind"] == "float"


def test_collect_inputs_required_and_optional():
    fields = [
        {"name": "smiles", "kind": "text", "required": True},
        {"name": "charge", "kind": "int", "required": True},
        {"name": "method", "kind": "text", "required": False},
    ]
    inputs, missing = TA.collect_inputs(
        fields, {"smiles": "CCO", "charge": 0, "method": ""})
    assert inputs == {"smiles": "CCO", "charge": 0}   # empty optional dropped
    assert missing == []

    _, missing2 = TA.collect_inputs(fields, {"smiles": "", "charge": 0})
    assert missing2 == ["smiles"]


# --- widget build ---------------------------------------------------------


def test_application_create_tab_builds():
    if not TA.HAS_WIDGETS:
        pytest.skip("ipywidgets not installed")
    widget, refs = TA.create_tab(ctx=None)
    assert widget is not None
    assert "application_form" in refs
