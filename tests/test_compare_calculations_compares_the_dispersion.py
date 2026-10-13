"""compare_calculations compares the dispersion correction and the solvent.

Input: two synthetic ORCA run folders (``.inp`` + ``.out``, optionally
``DELFIN_Data.json``). Output asserted: ``dispersion_match`` /
``solvent_match`` and a note naming the differing field, through
``delfin.api.compare_calculations`` and the ops-server tool wrapper.

The defect: a PBE0/def2-SVP pair one D3BJ apart came back with
``method_match`` and ``basis_match`` True, no note, and an energy delta
(-3.52 kcal/mol on a butane pair) that is mostly the dispersion term of
one run. Functionals whose name carries the dispersion (wB97X-V, *-3c,
GFN-xTB) must not be reported as a mismatch because one folder also
records a separate dispersion keyword.
"""
from __future__ import annotations

import json

import pytest

from delfin.api import compare_calculations


def _run(folder, *, functional="PBE0", disp="", solv="CPCM(DMF)",
         energy=-157.963184019, write_inp=True, data_disp=None):
    folder.mkdir(parents=True)
    kw = " ".join(t for t in ("!", functional, disp, "def2-SVP", solv,
                              "TightSCF", "Opt") if t)
    if write_inp:
        (folder / "mol.inp").write_text(
            kw + "\n* xyz 0 1\n  C 0.0 0.0 0.0\n*\n", encoding="utf-8")
    (folder / "mol.out").write_text("\n".join([
        "                                 * O   R   C   A *",
        f"|  1> {kw}",
        "|  2> * xyz 0 1",
        "Number of atoms                             ...     1",
        "  *** SCF CONVERGED ***",
        f"FINAL SINGLE POINT ENERGY     {energy:.9f}",
        "                             ****ORCA TERMINATED NORMALLY****",
    ]) + "\n", encoding="utf-8")
    if data_disp is not None:
        (folder / "DELFIN_Data.json").write_text(json.dumps({"metadata": {
            "functional": functional, "basis_set": "def2-SVP",
            "dispersion_correction": data_disp,
        }}), encoding="utf-8")
    return str(folder)


def test_a_pair_one_d3bj_apart_is_flagged(tmp_path):
    a = _run(tmp_path / "anti", disp="")
    b = _run(tmp_path / "gauche", disp="D3BJ", energy=-157.968792351)
    d = compare_calculations(a, b)
    assert d.method_match and d.basis_match       # functional/basis do match
    assert d.dispersion_match is False
    # "none": the keyword line was read and carries no correction.
    assert d.a_dispersion == "none" and d.b_dispersion == "D3BJ"
    notes = " ".join(d.notes)
    assert "dispersion" in notes and "D3BJ" in notes
    # The delta is still returned (additive change) but is called out.
    assert d.delta_spe_kcal == pytest.approx(-3.519, abs=1e-2)
    assert any("delta" in n and "dispersion" in n for n in d.notes)


def test_an_identical_pair_matches_with_no_note(tmp_path):
    a = _run(tmp_path / "a", disp="D3BJ")
    b = _run(tmp_path / "b", disp="D3BJ", energy=-157.9694)
    d = compare_calculations(a, b)
    assert d.method_match and d.basis_match
    assert d.dispersion_match is True and d.solvent_match is True
    assert d.notes == []


@pytest.mark.parametrize("functional", ["wB97X-V", "wB97M-V", "r2SCAN-3c",
                                        "wB97X-D3"])
def test_built_in_dispersion_is_no_mismatch(tmp_path, functional):
    # One folder also records a dispersion keyword (a CONTROL default the
    # composite method ignores); the method name carries it on both.
    a = _run(tmp_path / "a", functional=functional, data_disp="D4")
    b = _run(tmp_path / "b", functional=functional, energy=-157.97)
    d = compare_calculations(a, b)
    assert d.method_match
    assert d.dispersion_match is True
    assert not any("dispersion" in n for n in d.notes)


def test_built_in_against_an_added_correction_is_a_mismatch(tmp_path):
    a = _run(tmp_path / "a", functional="wB97X-V")
    b = _run(tmp_path / "b", functional="PBE0", disp="D3BJ")
    d = compare_calculations(a, b)
    assert d.dispersion_match is False


def test_the_output_echo_alone_is_read(tmp_path):
    # No .inp, no CONTROL.txt: the keyword line echoed in the .out is
    # what ran.
    a = _run(tmp_path / "a", write_inp=False)
    b = _run(tmp_path / "b", disp="D3BJ", write_inp=False, energy=-157.9688)
    d = compare_calculations(a, b)
    assert d.dispersion_match is False
    assert any("D3BJ" in n for n in d.notes)


def test_a_solvent_difference_is_flagged(tmp_path):
    a = _run(tmp_path / "a", disp="D3BJ", solv="")
    b = _run(tmp_path / "b", disp="D3BJ")
    d = compare_calculations(a, b)
    assert d.dispersion_match is True
    assert d.solvent_match is False
    assert any("solvent" in n and "DMF" in n for n in d.notes)


def test_the_ops_tool_carries_the_fields(tmp_path):
    from delfin.ops_server import server as ops_server
    a = _run(tmp_path / "anti")
    b = _run(tmp_path / "gauche", disp="D3BJ", energy=-157.968792351)
    data = json.loads(ops_server.tool_compare_calculations(a, b))
    # Additive: every earlier key is still there.
    for key in ("folder_a", "folder_b", "method_match", "basis_match",
                "a_functional", "b_functional", "a_basis", "b_basis",
                "a_gibbs", "b_gibbs", "delta_gibbs_kcal", "a_spe", "b_spe",
                "delta_spe_kcal", "a_n_imag", "b_n_imag", "notes"):
        assert key in data, key
    assert data["dispersion_match"] is False
    assert data["b_dispersion"] == "D3BJ"
    assert any("dispersion" in n for n in data["notes"])


# Built-in dispersion is keyed by the correction the name carries, not
# by "has one": a single BUILT-IN key made each of these pairs equal.
# Several of the names are not in the indexer's functional list, so the
# folder's functional is None and only the keyword line names them.
@pytest.mark.parametrize("fa,da,fb,db", [
    ("B3LYP-D3", "", "B3LYP-D3BJ", ""),      # zero vs Becke-Johnson damping
    ("wB97X-V", "", "wB97X-D3BJ", ""),       # VV10 vs D3BJ
    ("B97M-V", "", "B97M-D3BJ", ""),         # VV10 vs D3BJ
    ("wB97X-D", "", "wB97X-D3", ""),         # Chai-Head-Gordon D2 vs D3
    ("r2SCAN-3c", "", "B97-3c", ""),         # each 3c has its own fit
    ("r2SCAN-3c", "", "r2SCAN", "D4"),       # 3c-fitted D4 + gCP vs plain D4
    ("GFN2-xTB", "", "GFN1-xTB", ""),
    ("PBE0", "D4", "PBE0", "D3BJ"),
])
def test_different_dispersion_physics_is_a_mismatch(tmp_path, fa, da, fb, db):
    a = _run(tmp_path / "a", functional=fa, disp=da)
    b = _run(tmp_path / "b", functional=fb, disp=db, energy=-157.97)
    d = compare_calculations(a, b)
    assert d.dispersion_match is False, (d.a_dispersion, d.b_dispersion)
    assert any("dispersion" in n for n in d.notes)


@pytest.mark.parametrize("fa,da,fb,db", [
    ("B3LYP-D3", "", "B3LYP", "D3"),          # suffix vs separate keyword
    ("B3LYP", "D3", "B3LYP", "D3ZERO"),        # ORCA's D3 is D3zero
    ("B3LYP", "D3BJ", "B3LYP", "D3(BJ)"),
    ("B3LYP-D3BJ", "", "B3LYP", "D3BJ"),
])
def test_one_correction_spelled_two_ways_matches(tmp_path, fa, da, fb, db):
    a = _run(tmp_path / "a", functional=fa, disp=da)
    b = _run(tmp_path / "b", functional=fb, disp=db, energy=-157.97)
    d = compare_calculations(a, b)
    assert d.dispersion_match is True, (d.a_dispersion, d.b_dispersion)
    assert not any("dispersion" in n for n in d.notes)


def _bare_out(folder, energy=-157.97):
    # An output with no echoed keyword line, no .inp, no CONTROL.txt:
    # nothing in the folder says which correction ran.
    folder.mkdir(parents=True)
    (folder / "mol.out").write_text("\n".join([
        "  *** SCF CONVERGED ***",
        f"FINAL SINGLE POINT ENERGY     {energy:.9f}",
        "                             ****ORCA TERMINATED NORMALLY****",
    ]) + "\n", encoding="utf-8")
    return str(folder)


def test_an_unreadable_dispersion_never_claims_a_match(tmp_path):
    a = _bare_out(tmp_path / "a")
    b = _bare_out(tmp_path / "b", energy=-157.98)
    d = compare_calculations(a, b)
    assert d.dispersion_match is None
    assert d.a_dispersion is None and d.b_dispersion is None
    assert any("dispersion" in n and "not" in n for n in d.notes)


def test_unknown_against_none_is_not_a_match(tmp_path):
    a = _bare_out(tmp_path / "a")
    b = _run(tmp_path / "b")
    d = compare_calculations(a, b)
    assert d.dispersion_match is None
    assert d.b_dispersion == "none"


def test_the_solvation_model_is_compared(tmp_path):
    a = _run(tmp_path / "a", disp="D3BJ", solv="CPCM(water)")
    b = _run(tmp_path / "b", disp="D3BJ", solv="SMD(water)")
    d = compare_calculations(a, b)
    assert d.solvent_match is False
    assert d.a_solvation == "CPCM" and d.b_solvation == "SMD"
    assert any("CPCM" in n and "SMD" in n for n in d.notes)


def test_smd_switched_on_in_the_cpcm_block_is_smd(tmp_path):
    # ``! CPCM(water)`` with ``%cpcm smd true end`` runs SMD.
    a = _run(tmp_path / "a", disp="D3BJ", solv="SMD(water)")
    b = _run(tmp_path / "b", disp="D3BJ", solv="CPCM(water)", energy=-157.97)
    out = tmp_path / "b" / "mol.out"
    out.write_text(out.read_text(encoding="utf-8").replace(
        "|  2> * xyz 0 1",
        "|  2> %cpcm\n|  3>   smd true\n|  4> end\n|  5> * xyz 0 1"),
        encoding="utf-8")
    d = compare_calculations(a, b)
    assert d.b_solvation == "SMD"
    assert d.solvent_match is True, d.notes
