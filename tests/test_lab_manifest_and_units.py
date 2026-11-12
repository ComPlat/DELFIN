"""The manifest, the unit table and the interlock grammar refuse what is unclear.

A manifest error found at load time is a write that never happens: these
pin the closed schema (unknown keys, inline secrets, unknown units,
missing safe values, interlocks naming channels that do not exist, drivers
that are only planned) and the three-valued logic that makes an unknown
reading fail closed.
"""

from __future__ import annotations

import copy

import pytest
import yaml

from delfin.lab import expr, manifest, units


def _doc():
    return yaml.safe_load(manifest.demo_text())


def _device(doc, dev_id):
    for m in doc["modules"]:
        for d in m["devices"]:
            if d["id"] == dev_id:
                return d
    raise KeyError(dev_id)


def test_the_demo_manifest_loads_and_has_a_stable_checksum():
    a = manifest.parse_document(_doc())
    b = manifest.parse_document(copy.deepcopy(_doc()))
    assert a.checksum == b.checksum and len(a.checksum) == 64
    assert a.devices["hotplate"].writes["power"].states == ("on", "off")
    assert a.level(a.devices["mains_bypass"],
                   a.devices["mains_bypass"].writes["power"]) == "forbidden"
    doc = _doc()
    doc["version"] = 2
    assert manifest.parse_document(doc).checksum != a.checksum


@pytest.mark.parametrize("address", [
    {"host": "192.0.2.10", "password": "hunter2"},
    {"host": "192.0.2.10", "api_token": "abc"},
    {"host": "192.0.2.10", "secret_ref": "hunter2"},
])
def test_a_manifest_with_an_inline_secret_is_refused(address):
    doc = _doc()
    _device(doc, "lamp")["address"] = address
    with pytest.raises(manifest.ManifestError, match="secret"):
        manifest.parse_document(doc)
    _device(doc, "lamp")["address"] = {"host": "192.0.2.10",
                                       "secret_ref": "cred:lamp_plug"}
    manifest.parse_document(doc)                   # control


@pytest.mark.parametrize("mutate,needle", [
    (lambda d: _device(d, "hotplate")["writes"]["setpoint"].update(
        unit="C"), "ambiguous"),
    (lambda d: _device(d, "hotplate")["writes"]["setpoint"].pop("unit"),
     "no unit"),
    (lambda d: _device(d, "hotplate")["safe_state"].pop("setpoint"),
     "safe_state has no value"),
    (lambda d: _device(d, "hotplate")["safe_state"].update(setpoint=500),
     "outside its range"),
    (lambda d: _device(d, "hotplate")["interlocks"][0].update(
        requires="stirer.power == on"), "unknown channel stirer.power"),
    # shelly_rpc, tasmota and mqtt are implemented since their drivers landed;
    # home_assistant stays the planned-driver example
    (lambda d: _device(d, "lamp").update(driver="home_assistant"), "planned"),
    (lambda d: _device(d, "lamp").update(colour="red"), "unknown key"),
    (lambda d: _device(d, "pump")["writes"]["flow"].update(level="maybe"),
     "level must be"),
])
def test_an_unclear_manifest_is_refused(mutate, needle):
    doc = _doc()
    mutate(doc)
    with pytest.raises(manifest.ManifestError) as exc:
        manifest.parse_document(doc)
    assert needle in str(exc.value)


def test_units_convert_within_a_dimension_only():
    assert units.convert(313.15, "K", "degC") == pytest.approx(40.0)
    assert units.convert(212, "degF", "degC") == pytest.approx(100.0)
    assert units.convert(60, "mL/h", "mL/min") == pytest.approx(1.0)
    assert units.canonical("°C") == "degC"
    for bad in (None, "", "C", "F", "m", "parsec"):
        with pytest.raises(units.UnitError):
            units.canonical(bad)
    with pytest.raises(units.UnitError):
        units.convert(1, "rpm", "degC")


def test_interlocks_are_three_valued_and_unknown_fails_closed():
    tree = expr.parse("a.x > 40 degC and (b.on == on or c.t < 300 K)")
    vals = {"a.x": (50.0, "degC"), "b.on": ("off", None),
            "c.t": (20.0, "degC")}

    def look(ref):
        return vals.get(ref.name, expr.UNKNOWN)

    assert expr.evaluate(tree, look) is True           # 20 degC < 300 K
    vals["c.t"] = expr.UNKNOWN
    assert expr.evaluate(tree, look) is expr.UNKNOWN
    vals["a.x"] = (30.0, "degC")
    assert expr.evaluate(tree, look) is False          # False and UNKNOWN
    with pytest.raises(TypeError):
        bool(expr.UNKNOWN)
    for bad in ("", "a.x >", "__import__('os')", "a.x > 1 and", "x > 1",
                "a.x < on"):
        with pytest.raises((expr.ExprError, units.UnitError)):
            expr.parse(bad)
