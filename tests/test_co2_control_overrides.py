"""Main CONTROL.txt → CO2 Coordinator CONTROL.txt override tests.

Users should be able to tweak any CO2-Coordinator setting from the
main DELFIN CONTROL.txt instead of editing the separate CO2-Coordinator
CONTROL.txt manually.

Two ways to override, both equivalent:
  - Plain key:   ``place_axis=x``
  - CO2 prefix:  ``co2_place_axis=x``  (wins if both are set — explicit)

Fields users have defaults for (from the CO2 template):
  xyz, out, co2, orientation_distance, rot_step_deg, rot_range_deg,
  orientation_job, scan_job, scan_end, scan_steps, metal, metal_index,
  align_bond_index, neighbors, place_axis, mode, perp_axis,
  place_optimize, place_samples, place_clearance_scale, no_place_co2,
  parallel_orientation_scan, max_workers

Plus anything the pipeline already owns (functional, main_basisset,
metal_basisset, solvent, PAL, maxcore, ...).
"""

from __future__ import annotations

import pytest

from delfin.co2.chain_setup import (
    _build_co2_control,
    _co2_override_value,
    _CO2_DEFAULTS,
    _TRANSFER_KEYS,
)


# ---------------------------------------------------------------------------
# _co2_override_value helper
# ---------------------------------------------------------------------------


def test_plain_key_override():
    ctrl = {"place_axis": "x"}
    assert _co2_override_value(ctrl, "place_axis") == "x"


def test_prefixed_key_override():
    ctrl = {"co2_place_axis": "-z"}
    assert _co2_override_value(ctrl, "place_axis") == "-z"


def test_prefixed_key_wins_over_plain():
    """co2_<key> is the explicit form and takes priority."""
    ctrl = {"place_axis": "x", "co2_place_axis": "-y"}
    assert _co2_override_value(ctrl, "place_axis") == "-y"


def test_missing_key_returns_empty():
    assert _co2_override_value({}, "place_axis") == ""
    assert _co2_override_value({"place_axis": ""}, "place_axis") == ""
    assert _co2_override_value({"co2_place_axis": ""}, "place_axis") == ""


def test_empty_prefixed_falls_back_to_plain():
    """Empty co2_<key> should not shadow a plain key with a real value."""
    ctrl = {"place_axis": "y", "co2_place_axis": ""}
    assert _co2_override_value(ctrl, "place_axis") == "y"


# ---------------------------------------------------------------------------
# _TRANSFER_KEYS coverage
# ---------------------------------------------------------------------------


def test_every_co2_default_is_transferable():
    """Every field that has a CO2 default MUST also be overridable from
    the main CONTROL.txt. If _CO2_DEFAULTS gets a new key, _TRANSFER_KEYS
    must include it — otherwise users can't override the default."""
    missing = set(_CO2_DEFAULTS.keys()) - set(_TRANSFER_KEYS)
    # charge/multiplicity/broken_sym are set dynamically, not from CONTROL
    missing -= {"charge", "multiplicity", "broken_sym"}
    assert not missing, (
        f"CO2 defaults {sorted(missing)} are not in _TRANSFER_KEYS — "
        f"users would have no way to override them from the main CONTROL.txt."
    )


def test_transfer_keys_include_pipeline_level_of_theory():
    """LoT fields that the main pipeline already owns must also reach
    the CO2 coordinator so the two stay on the same level of theory."""
    required_lot = {
        "functional", "main_basisset", "metal_basisset",
        "first_coordination_sphere_metal_basisset",
        "solvent", "implicit_solvation_model",
    }
    missing = required_lot - set(_TRANSFER_KEYS)
    assert not missing, (
        f"Level-of-theory keys {sorted(missing)} missing from _TRANSFER_KEYS"
    )


def test_charge_multiplicity_broken_sym_not_transferred():
    """These are set per-species, not from CONTROL."""
    for dynamic_key in ("charge", "multiplicity", "broken_sym"):
        assert dynamic_key not in _TRANSFER_KEYS, (
            f"{dynamic_key} is dynamic per species and must not be "
            f"transferred from the main CONTROL.txt."
        )


# ---------------------------------------------------------------------------
# End-to-end: _build_co2_control
# ---------------------------------------------------------------------------


def test_build_co2_control_uses_main_control_overrides():
    """Overrides from the main CONTROL.txt must appear in the generated
    CO2 Coordinator CONTROL.txt — both plain and prefixed forms."""
    delfin_ctrl = {
        "place_axis": "-z",                  # plain override
        "co2_mode": "end-on",                # prefixed override
        "scan_end": "2.5",                   # plain, CO2-specific
        "co2_scan_steps": "50",              # prefixed, CO2-specific
        "functional": "BP86",                # pipeline LoT
        "main_basisset": "def2-TZVP",        # pipeline LoT
    }
    txt = _build_co2_control(
        charge=0, multiplicity=1, broken_sym="", delfin_ctrl=delfin_ctrl,
    )

    # Plain overrides survive
    assert "place_axis=-z" in txt, "place_axis override lost in CO2 CONTROL"
    assert "scan_end=2.5" in txt, "scan_end override lost in CO2 CONTROL"

    # Prefixed overrides survive
    assert "mode=end-on" in txt, "co2_mode override lost in CO2 CONTROL"
    assert "scan_steps=50" in txt, "co2_scan_steps override lost in CO2 CONTROL"

    # Pipeline LoT forwarded
    assert "functional=BP86" in txt
    assert "main_basisset=def2-TZVP" in txt


def test_build_co2_control_prefix_wins_over_plain():
    """If both ``key`` and ``co2_key`` are set in the main CONTROL.txt,
    the prefixed form wins."""
    delfin_ctrl = {
        "place_axis": "x",        # plain
        "co2_place_axis": "-y",   # prefix
    }
    txt = _build_co2_control(
        charge=0, multiplicity=1, broken_sym="", delfin_ctrl=delfin_ctrl,
    )
    assert "place_axis=-y" in txt
    assert "place_axis=x" not in txt


def test_build_co2_control_falls_back_to_defaults():
    """If nothing is overridden, the CO2 defaults must appear verbatim."""
    txt = _build_co2_control(
        charge=0, multiplicity=1, broken_sym="", delfin_ctrl={},
    )
    for key, default in _CO2_DEFAULTS.items():
        assert f"{key}={default}" in txt, (
            f"Default value for {key!r} missing from generated CO2 CONTROL"
        )


def test_build_co2_control_dynamic_fields_always_fresh():
    """charge/multiplicity/broken_sym come from the args, not from
    whatever the main CONTROL.txt happened to have."""
    txt = _build_co2_control(
        charge=-2, multiplicity=3, broken_sym="2,1",
        delfin_ctrl={"charge": "999", "multiplicity": "7"},
    )
    assert "charge=-2" in txt
    assert "multiplicity=3" in txt
    assert "charge=999" not in txt
    assert "broken_sym=2,1" in txt


def test_co2_multiplicity_override_wins_over_detected_spin(tmp_path, monkeypatch, capsys):
    """An explicit ``co2_multiplicity`` in the main CONTROL.txt overrides
    the detected per-species multiplicity when the CO2 coordinator
    CONTROL.txt is generated (fixes: override multiplicities were not
    carried over to the coordinator)."""
    from delfin.co2 import chain_setup

    # Job dir with a CONTROL.txt carrying the override + a species xyz
    (tmp_path / "CONTROL.txt").write_text(
        "charge=0\nco2_multiplicity=3\n", encoding="utf-8"
    )
    (tmp_path / "initial.xyz").write_text("1\n\nH 0 0 0\n", encoding="utf-8")

    # Pretend spin detection found mult=1 / no BS
    monkeypatch.setattr(
        chain_setup, "_spin_from_state_json", lambda jd, d: (1, None)
    )

    chain_setup.setup_co2_from_delfin(tmp_path, 0)

    co2_ctrl = (tmp_path / "CO2_coordination" / "CONTROL.txt").read_text()
    assert "multiplicity=3" in co2_ctrl, (
        "co2_multiplicity override from main CONTROL.txt did not reach "
        f"the CO2 coordinator CONTROL.txt:\n{co2_ctrl}"
    )
    assert "multiplicity=1" not in co2_ctrl


def test_co2_multiplicity_override_invalid_keeps_detected(tmp_path, monkeypatch, capsys):
    """A non-int ``co2_multiplicity`` is ignored with a warning; the
    detected multiplicity survives."""
    from delfin.co2 import chain_setup

    (tmp_path / "CONTROL.txt").write_text(
        "charge=0\nco2_multiplicity=abc\n", encoding="utf-8"
    )
    (tmp_path / "initial.xyz").write_text("1\n\nH 0 0 0\n", encoding="utf-8")

    monkeypatch.setattr(
        chain_setup, "_spin_from_state_json", lambda jd, d: (2, "1,0")
    )

    chain_setup.setup_co2_from_delfin(tmp_path, 0)

    co2_ctrl = (tmp_path / "CO2_coordination" / "CONTROL.txt").read_text()
    assert "multiplicity=2" in co2_ctrl
    assert "WARN" in capsys.readouterr().out


def test_plain_multiplicity_still_not_forwarded(tmp_path, monkeypatch):
    """A plain ``multiplicity`` key (base-species info) must NOT leak into
    the coordinator CONTROL — only the explicit co2_ prefix overrides."""
    from delfin.co2 import chain_setup

    (tmp_path / "CONTROL.txt").write_text(
        "charge=0\nmultiplicity=7\n", encoding="utf-8"
    )
    (tmp_path / "initial.xyz").write_text("1\n\nH 0 0 0\n", encoding="utf-8")

    monkeypatch.setattr(
        chain_setup, "_spin_from_state_json", lambda jd, d: (2, None)
    )

    chain_setup.setup_co2_from_delfin(tmp_path, 0)

    co2_ctrl = (tmp_path / "CO2_coordination" / "CONTROL.txt").read_text()
    assert "multiplicity=2" in co2_ctrl
    assert "multiplicity=7" not in co2_ctrl
