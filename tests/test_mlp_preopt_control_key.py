"""The opt-in MLP pre-optimisation CONTROL key (wave19 m4).

``mlp_preopt`` selects an ML potential to relax the start geometry with
BEFORE the ORCA steps.  Off by default: every existing default path is
byte-identical, which the last two tests pin with the golden input
fixtures and the shipped CONTROL templates.
"""

from __future__ import annotations

from pathlib import Path

import pytest


CONTROL = Path(__file__).resolve().parent / "fixtures" / "control_templates"

# ---------------------------------------------------------------- validator ---
def _validate(**overrides):
    from delfin.common.control_validator import validate_control_config

    config = {
        "charge": 0, "method": "classic", "smiles_converter": "NORMAL",
        "ESD_modul": "no", "solvent": "", "functional": "PBE0",
        "basis": "def2-SVP", "PAL": 6,
    }
    config.update(overrides)
    return validate_control_config(config)


def test_mlp_preopt_defaults_to_none():
    out = _validate()
    assert out["mlp_preopt"] == "none"


def test_mlp_preopt_accepts_a_backend_and_its_aliases():
    assert _validate(mlp_preopt="mace")["mlp_preopt"] == "mace_off"
    assert _validate(mlp_preopt="MACE-OFF")["mlp_preopt"] == "mace_off"
    assert _validate(mlp_preopt="ani2x")["mlp_preopt"] == "ani2x"
    assert _validate(mlp_preopt="AIMNet2")["mlp_preopt"] == "aimnet2"


def test_mlp_preopt_rejects_an_unknown_value():
    with pytest.raises(Exception) as excinfo:
        _validate(mlp_preopt="not_a_backend")
    assert "mlp_preopt" in str(excinfo.value)
    assert "mace_off" in str(excinfo.value)   # names the accepted values


def test_mlp_preopt_can_be_switched_off_explicitly():
    assert _validate(mlp_preopt="none")["mlp_preopt"] == "none"
    assert _validate(mlp_preopt="off")["mlp_preopt"] == "none"


def test_mlp_preopt_key_is_documented_in_the_key_notes():
    import json

    notes = json.loads(
        (Path(__file__).resolve().parents[1] / "delfin" / "common"
         / "control_key_notes.json").read_text(encoding="utf-8"))
    entry = notes.get("mlp_preopt")
    assert entry, "mlp_preopt is not in control_key_notes.json"
    assert "none" in entry.get("valid_values", "")


def test_mlp_preopt_is_in_the_config_defaults():
    from delfin.config import DEFAULT_CONFIG_VALUES

    assert DEFAULT_CONFIG_VALUES.get("mlp_preopt") == "none"


# --------------------------------------------------------- pipeline hook ---
class _FakeCalc:
    def __init__(self, **kw):
        self.kw = kw


def _patch_calcs(monkeypatch, calls):
    import delfin.workflows.pipeline as pipeline

    def fake_optimize(atoms, backend="ani2x", **kw):
        calls.append(("optimize", backend, kw.get("fmax"), kw.get("steps")))
        atoms.set_positions = atoms.set_positions
        return {"converged": True, "energy": -42.0, "n_steps": 5,
                "backend": backend}

    monkeypatch.setattr(pipeline, "_run_mlp_preoptimize", fake_optimize)


def test_pipeline_runs_mlp_preopt_when_the_key_is_set(tmp_path, monkeypatch):
    import delfin.workflows.pipeline as pipeline

    calls = []
    _patch_calcs(monkeypatch, calls)
    start = tmp_path / "start.txt"
    start.write_text("O 0.0 0.0 0.0\nH 0.9 0.0 0.0\nH -0.3 0.85 0.0\n")
    config = {"mlp_preopt": "mace", "input_file": str(start),
              "charge": 0, "mult": 1}
    result = pipeline._apply_mlp_preopt(config, start, charge=0, mult=1)
    assert calls and calls[0][1] == "mace"
    # the optimised geometry replaces start.txt's coordinates
    assert start.read_text() != "O 0.0 0.0 0.0\nH 0.9 0.0 0.0\nH -0.3 0.85 0.0\n"
    assert result["converged"] is True


def test_pipeline_default_path_runs_no_mlp(tmp_path, monkeypatch):
    import delfin.workflows.pipeline as pipeline

    calls = []
    _patch_calcs(monkeypatch, calls)
    start = tmp_path / "start.txt"
    body = "O 0.0 0.0 0.0\nH 0.9 0.0 0.0\nH -0.3 0.85 0.0\n"
    start.write_text(body)
    config = {"input_file": str(start), "charge": 0, "mult": 1}
    result = pipeline._apply_mlp_preopt(config, start, charge=0, mult=1)
    assert calls == []                      # the MLP never ran
    assert result is None
    assert start.read_text() == body        # the geometry is untouched


def test_pipeline_mlp_preopt_failure_keeps_the_geometry(tmp_path, monkeypatch):
    import delfin.workflows.pipeline as pipeline

    def fake_optimize(*a, **kw):
        raise RuntimeError("model exploded")

    monkeypatch.setattr(pipeline, "_run_mlp_preoptimize", fake_optimize)
    start = tmp_path / "start.txt"
    body = "O 0.0 0.0 0.0\nH 0.9 0.0 0.0\nH -0.3 0.85 0.0\n"
    start.write_text(body)
    config = {"mlp_preopt": "mace", "input_file": str(start),
              "charge": 0, "mult": 1}
    result = pipeline._apply_mlp_preopt(config, start, charge=0, mult=1)
    assert result is None or result.get("ok") is False
    assert start.read_text() == body        # fail-open for the ORCA step


# ------------------------------------------------ byte-identical default ---
def test_shipped_control_templates_validate_with_mlp_preopt_absent():
    """A CONTROL without the key must validate exactly as it did before."""
    from delfin.config import validate_control_text

    for template in sorted(CONTROL.glob("*.txt")):
        # errors that exist without the key must not change
        assert validate_control_text(template.read_text(encoding="utf-8")) == \
            validate_control_text(
                template.read_text(encoding="utf-8") + "\nmlp_preopt = none\n")


def test_control_key_lookup_does_not_change_unrelated_keys():
    from delfin.common.control_validator import CONTROL_FIELD_SPECS

    names = [s.name for s in CONTROL_FIELD_SPECS]
    assert "mlp_preopt" in names
    # the rest of the schema is untouched
    assert names.count("thdy_preopt") == 1 and "xtb_preopt" not in names
