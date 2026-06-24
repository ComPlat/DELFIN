"""Tests for the enriched native CREST/CENSO building blocks (Track 3)."""

from __future__ import annotations

from delfin.tools import platform
from delfin.tools._registry import get
from delfin.tools.adapters.crest import (
    _count_xyz_frames,
    _parse_crest_outputs,
    build_crest_cmd,
)


# --- CREST -----------------------------------------------------------------


def test_crest_cmd_builder():
    cmd = build_crest_cmd("m.xyz", charge=-1, uhf=1, cores=4, ewin=4.0,
                          method="gfnff", solvent="water", extra_args=["--nci"])
    assert cmd[:2] == ["crest", "m.xyz"]
    assert cmd[cmd.index("--chrg") + 1] == "-1"
    assert cmd[cmd.index("--uhf") + 1] == "1"
    assert cmd[cmd.index("-T") + 1] == "4"
    assert "--gfnff" in cmd
    assert cmd[cmd.index("--alpb") + 1] == "water"
    assert cmd[-1] == "--nci"                       # raw mode flag passes through
    # gfn2 is the default level → no extra flag; gfn1 adds its flag
    assert "--gfnff" not in build_crest_cmd("m.xyz", method="gfn2")
    assert "--gfn1" in build_crest_cmd("m.xyz", method="gfn1")


def test_crest_contract_enriched():
    c = get("crest_conformers").contract()
    assert {d.name for d in c.data_keys} == {"n_conformers", "lowest_energy_Eh"}
    assert set(c.param("method").enum) >= {"gfn2", "gfn1", "gfnff"}
    assert "ensemble" in c.produces


def test_crest_output_parsing(tmp_path):
    ens = tmp_path / "crest_conformers.xyz"
    ens.write_text("1\nc1\nH 0 0 0\n1\nc2\nH 0 0 1\n1\nc3\nH 0 0 2\n")
    assert _count_xyz_frames(ens) == 3
    (tmp_path / "crest_best.xyz").write_text("1\n  -5.123456\nH 0 0 0\n")
    out = _parse_crest_outputs(tmp_path)
    assert out["n_conformers"] == 3
    assert out["lowest_energy_Eh"] == -5.123456


# --- CENSO + the crest→censo chain ----------------------------------------


def test_censo_wires_ensemble():
    assert get("censo_sort").contract().wires == (("ensemble", "ensemble"),)
    assert {d.name for d in get("censo_sort").contract().data_keys} >= {"output_dir"}


def test_crest_to_censo_auto_wires_ensemble():
    # CREST produces the ensemble → CENSO's required ensemble path auto-wires
    rep = platform.validate_spec({"name": "t", "steps": [
        {"step": "smiles_to_xyz", "kwargs": {"smiles": "CCO"}},
        {"step": "xtb_optimize", "kwargs": {"charge": 0}},
        {"step": "crest_conformers", "kwargs": {"charge": 0}},
        {"step": "censo_sort", "kwargs": {"charge": 0}}]}, geometry=False)
    assert rep.ok, rep.summary()
    # standalone CENSO still needs the ensemble
    bad = platform.validate_spec(
        {"name": "t", "steps": [{"step": "censo_sort", "kwargs": {"charge": 0}}]},
        geometry=False)
    assert any("ensemble" in d.missing_params for d in bad.diagnostics)


def test_conformer_dG_application():
    assert "conformer_dG" in platform.list_applications()
    rep = platform.validate_application("conformer_dG", geometry=False,
                                        smiles="CCO", charge=0)
    assert rep.ok, rep.summary()
    leaks = [m for d in rep.diagnostics for m in d.messages if "placeholder" in m]
    assert leaks == []
    app = platform.describe_application("conformer_dG")
    assert {o.name for o in app.outputs} == {"n_conformers", "ranked_dir"}
