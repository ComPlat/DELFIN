"""OPI-based tests for the orca-mlips ExtOpt builder (delfin/mlp_tools/orca_mlips.py).

Runs in the py3.12 OPI venv (CI sets DELFIN_OPI_REQUIRED=1); the input is
built through OPI and compared against the ORCA 6.1.1 manual, chapter 4.1.10
(structurereactivity/optimizations.txt lines 1537-1701) and the
orca-mlips README/OPTIONS.md
(https://github.com/t-0hmura/orca-mlips). No ORCA run, no orca-mlips
package needed: the builder only writes an input file.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.mlp_tools.orca_mlips import (
    build_orca_mlips_input,
    detect_orca_mlips,
    ext_params_for,
    mlips_progext_names,
)

WATER = (["O", "H", "H"], [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]])


def _read(work_dir, base_name="extopt"):
    return (Path(work_dir) / f"{base_name}.inp").read_text()


# ---------------------------------------------------------------- detection


def test_detect_reports_missing_package_without_raising(monkeypatch):
    """No orca-mlips installed -> available=False + install hint, never a crash."""
    from delfin.mlp_tools import orca_mlips as om

    monkeypatch.setattr(om, "_find_spec", lambda name: None)
    info = detect_orca_mlips()
    assert info["available"] is False
    assert "orca-mlips" in info["install_hint"]
    assert info["progext_names"] == []


def test_detect_lists_progext_names_from_entry_points(monkeypatch):
    """Entry points discovered from the installed package's metadata."""
    import importlib.metadata as im
    from types import SimpleNamespace
    from delfin.mlp_tools import orca_mlips as om

    monkeypatch.setattr(om, "_find_spec", lambda name: object())

    class _Dist:
        version = "1.3.0"

    monkeypatch.setattr(im, "distribution", lambda name: _Dist())
    eps = [SimpleNamespace(name="orca-mlips-uma", value="orca_mlips.uma:main")]
    monkeypatch.setattr(im, "entry_points", lambda **kw: eps)
    info = detect_orca_mlips()
    assert info["available"] is True
    assert "orca-mlips-uma" in info["progext_names"]
    assert info["version"] == "1.3.0"


# ---------------------------------------------------------------- Ext_Params


def test_ext_params_defaults():
    params = ext_params_for("uma")
    assert "--backend uma" in params
    assert "--device" in params


def test_ext_params_solvent():
    params = ext_params_for("mace", solvent="water", solvent_model="cpcmx")
    assert "--solvent water" in params
    assert "--solvent-model cpcmx" in params


def test_ext_params_rejects_bad_solvent_model():
    with pytest.raises(ValueError, match="solvent-model"):
        ext_params_for("mace", solvent="water", solvent_model="magic")


def test_ext_params_device_and_model():
    params = ext_params_for("orb", model="orb-v3-conservative-omol", device="cuda")
    assert "--model orb-v3-conservative-omol" in params
    assert "--device cuda" in params


def test_ext_params_rejects_unknown_backend():
    with pytest.raises(ValueError, match="backend"):
        ext_params_for("ani2x")


# ---------------------------------------------------------------- ProgExt names


def test_progext_names_prefixed_alias_documented():
    """README troubleshooting: 'ProgExt "uma" runs the wrong plugin -- use orca-mlips-uma'."""
    names = mlips_progext_names(["uma", "orb"])
    assert names == ("orca-mlips-uma", "orca-mlips-orb")


# ---------------------------------------------------------------- input build


def test_opt_input_produces_extopt_and_progext(tmp_path):
    elements, coords = WATER
    build_orca_mlips_input(
        "opt", elements=elements, coordinates=coords, work_dir=tmp_path
    )
    text = _read(tmp_path)
    assert "extopt" in text.lower()
    assert 'progext "orca-mlips-uma"' in text.lower()
    assert "ext_params" in text.lower()
    assert "*xyzfile" in text.lower()
    assert "extopt.xyz" in text


def test_optts_two_step_hessian_block(tmp_path):
    """Two-step Hessian: --dump-hessian writes the .hess; %geom inhess Read +
    InHessName feed it back (manual optimizations_TS.txt lines 394-400)."""
    elements, coords = WATER
    build_orca_mlips_input(
        "optts", elements=elements, coordinates=coords, work_dir=tmp_path,
        hessian="two_step",
    )
    text = _read(tmp_path)
    assert "inhess read" in text.lower()
    assert 'inhessname "extopt.hess"' in text.lower()
    # Ext_Params carries the dump option
    assert "--dump-hessian extopt.hess" in text
    assert "freq" in text.lower()  # the Freq run reads the Hessian


def test_freq_input_uses_hessian_too(tmp_path):
    elements, coords = WATER
    build_orca_mlips_input(
        "freq", elements=elements, coordinates=coords, work_dir=tmp_path,
        hessian="two_step",
    )
    text = _read(tmp_path)
    assert "inhess read" in text.lower()


def test_unknown_run_type_rejected(tmp_path):
    elements, coords = WATER
    with pytest.raises(ValueError, match="run_type"):
        build_orca_mlips_input("neb", elements=elements, coordinates=coords,
                               work_dir=tmp_path)


def test_unknown_backend_rejected(tmp_path):
    elements, coords = WATER
    with pytest.raises(ValueError, match="backend"):
        build_orca_mlips_input(
            "opt", elements=elements, coordinates=coords, backend="ani2x",
            work_dir=tmp_path,
        )


def test_geometry_and_charge(tmp_path):
    elements, coords = WATER
    build_orca_mlips_input(
        "opt", elements=elements, coordinates=coords, charge=1, mult=2,
        work_dir=tmp_path,
    )
    text = _read(tmp_path)
    assert "*xyzfile 1 2 extopt.xyz" in text
