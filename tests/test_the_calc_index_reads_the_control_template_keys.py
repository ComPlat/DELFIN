"""The calc index reads the level-of-theory keys that CONTROL.txt carries.

Input: a calculation folder whose only description is a CONTROL.txt
written from ``delfin.define.TEMPLATE`` (no DELFIN_Data.json, no .inp).
Output: an index record with functional, basis set, metal basis,
dispersion, relativity, solvation, solvent, charge and multiplicity.

The template writes ``functional``, ``main_basisset``,
``metal_basisset``, ``disp_corr`` and ``relativity``; the indexer read
``basis_set`` only, a key the template never writes (the benchmark
fixtures write it, which hid the gap). Such a run was indexed with no
method at all, so ``search_calcs(functional=...)`` never returned it.

Basis selection follows ``delfin.utils.set_main_basisset``: the ``*_rel``
values replace the plain ones only when the structure holds a 4d/5d/4f/5f
metal, independent of ``relativity``. With no structure to decide from,
the basis is left empty for the .inp header to fill.
"""

from __future__ import annotations

import json
import re

import pytest

import delfin.agent.api_client as A
from delfin import define
from delfin.doc_server import calc_indexer as ci
from delfin.doc_server.calc_indexer import build_calc_index
from delfin.doc_server.calc_search import CalcSearchEngine


def _template_kv() -> dict:
    kv = {}
    for line in define.TEMPLATE.splitlines():
        line = line.strip()
        if "=" in line and not line.startswith(("#", "-")):
            k, _, v = line.partition("=")
            kv[k.strip()] = v.strip()
    return kv


def _template_control(**overrides) -> str:
    """define.TEMPLATE with its placeholders filled, keys overridden."""
    text = define.TEMPLATE.replace("[CHARGE]", "0").replace("[SOLVENT]", "DMF")
    text = re.sub(r"=\[[^\]\n]*\]", "=", text)
    for key, value in overrides.items():
        text, n = re.subn(rf"(?m)^{re.escape(key)}=.*$", f"{key}={value}", text)
        assert n == 1, key
    return text


_XYZ_ORGANIC = "3\n\nO 0 0 0\nH 0 0 0.96\nH 0.93 0 -0.24\n"
_XYZ_RU = "2\n\nRu 0 0 0\nCl 0 0 2.3\n"
_XYZ_FE = "2\n\nFe 0 0 0\nCl 0 0 2.2\n"


def _calc(root, name, control, geometry=None, **extra):
    d = root / name
    d.mkdir(parents=True)
    (d / "CONTROL.txt").write_text(control, encoding="utf-8")
    if geometry is not None:
        (d / "input.txt").write_text(geometry, encoding="utf-8")
    for fname, body in extra.items():
        (d / fname.replace("__", ".")).write_text(body, encoding="utf-8")
    return d


def _record(tmp_path, name, control, geometry=None, **extra):
    calc = tmp_path / "calc"
    _calc(calc, name, control, geometry, **extra)
    idx = build_calc_index(calc_dir=calc, quiet=True)
    (rec,) = [r for r in idx["records"] if r["calc_id"] == name]
    return rec, idx


@pytest.fixture(autouse=True)
def _private_home(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))


def test_the_keys_read_are_the_keys_the_template_writes():
    """Pinned to define.TEMPLATE so a renamed key fails here, not in a
    search that silently returns nothing."""
    template = _template_kv()
    missing = [k for k in ci.CONTROL_TEMPLATE_KEYS if k not in template]
    assert not missing, missing
    for key in ("functional", "main_basisset", "main_basisset_rel",
                "metal_basisset", "metal_basisset_rel", "disp_corr",
                "relativity", "implicit_solvation_model", "solvent",
                "charge", "multiplicity_0", "input_file"):
        assert key in ci.CONTROL_TEMPLATE_KEYS, key


def test_a_template_control_alone_files_the_run_under_its_method(tmp_path):
    rec, _ = _record(tmp_path, "organic", _template_control(),
                     geometry=_XYZ_ORGANIC)
    assert rec["index_source"] == "CONTROL.txt"
    assert rec["functional"] == "PBE0"
    assert rec["basis_set"] == "def2-SVP"
    assert rec["dispersion"] == "D4"
    assert rec["solvation"] == "CPCM"
    assert rec["solvent"] == "DMF"
    assert rec["charge"] == "0"
    assert rec.get("relativity", "") == ""
    assert "pbe0" in rec["_search_text"] and "d4" in rec["_search_text"]


def test_a_heavy_metal_takes_the_rel_basis_and_the_hamiltonian(tmp_path):
    rec, _ = _record(tmp_path, "ru", _template_control(multiplicity_0="1"),
                     geometry=_XYZ_RU)
    assert rec["basis_set"] == "ZORA-def2-SVP"
    assert rec["metal_basis"] == "SARC-ZORA-TZVP"
    assert rec["relativity"] == "ZORA"
    assert rec["multiplicity"] == "1"


def test_a_3d_metal_keeps_the_plain_basis_although_relativity_is_set(tmp_path):
    rec, _ = _record(tmp_path, "fe", _template_control(), geometry=_XYZ_FE)
    assert rec["basis_set"] == "def2-SVP"
    assert rec["metal_basis"] == "def2-TZVP"
    assert rec.get("relativity", "") == ""


# input.txt as DELFIN writes it: coordinate lines, no atom count, no
# comment (define.convert_xyz_to_input_txt drops the first two lines).
_TXT_RU = "Ru 0.0 0.0 0.0\nCl 0.0 0.0 2.3\n"
_TXT_FE = "Fe 0.0 0.0 0.0\nCl 0.0 0.0 2.2\n"
_TXT_DY = "Dy 0.0 0.0 0.0\nCl 0.0 0.0 2.6\n"
_TXT_ORGANIC = "O 0.0 0.0 0.0\nH 0.0 0.0 0.96\nH 0.93 0.0 -0.24\n"


def _delfin_writes(tmp_path, control, geometry):
    """(main basis, metal basis, relativity token) of the ORCA input that
    DELFIN's own writer builds from this CONTROL and geometry, with the
    relativity reset that cli.py applies to a structure without a heavy
    metal."""
    from delfin.config import read_control_file
    from delfin.utils import _should_use_rel, search_transition_metals
    from delfin.xyz_io import read_and_modify_file
    work = tmp_path / "delfin_writes"
    work.mkdir(exist_ok=True)
    ctl = work / "CONTROL.txt"
    ctl.write_text(control.replace("method=\n", "method=classic\n").replace(
        "smiles_converter=\n", "smiles_converter=MANTA\n"), encoding="utf-8")
    geom = work / "input.txt"
    geom.write_text(geometry, encoding="utf-8")
    cfg = read_control_file(str(ctl))
    metals = search_transition_metals(str(geom))
    if not _should_use_rel(metals):
        cfg["relativity"] = "none"
    out = work / "out.inp"
    read_and_modify_file(str(geom), str(out), 0, 1, "DMF", metals, None, None, cfg, "")
    text = out.read_text(encoding="utf-8")
    toks = next(ln for ln in text.splitlines() if ln.startswith("!")).split()
    rel = next((t for t in toks if t in ("ZORA", "X2C", "DKH")), "")
    m = re.search(r'NewGTO "([^"]+)"', text)
    return toks[3 if rel else 2], (m.group(1) if m else ""), rel


@pytest.mark.parametrize("name,geometry,overrides", [
    ("organic", _TXT_ORGANIC, {}),
    ("fe", _TXT_FE, {}),
    ("ru", _TXT_RU, {}),
    ("dy", _TXT_DY, {}),
    ("ru_relativity_none", _TXT_RU, {"relativity": "none"}),
    ("ru_empty_rel_basis", _TXT_RU, {"main_basisset_rel": ""}),
    ("organic_empty_basis", _TXT_ORGANIC, {"main_basisset": ""}),
], ids=lambda v: v if isinstance(v, str) and "\n" not in v else "")
def test_the_index_records_the_basis_delfin_writes(tmp_path, name, geometry, overrides):
    """The index and DELFIN's input writer agree on main basis, metal
    basis and relativity, for the input.txt format DELFIN itself writes
    and for a CONTROL key left empty (DELFIN falls back to def2-SVP /
    ZORA-def2-SVP there)."""
    control = _template_control(**overrides)
    rec, _ = _record(tmp_path, name, control, geometry=geometry)
    got = (rec["basis_set"], rec.get("metal_basis", ""), rec.get("relativity", ""))
    assert got == _delfin_writes(tmp_path, control, geometry)


def test_an_oversized_geometry_is_not_read(tmp_path, monkeypatch):
    """A geometry above the size bound is not opened; the basis is then
    left to the .inp header, as with no geometry at all."""
    d = _calc(tmp_path, "big", _template_control(), geometry=_XYZ_RU)
    assert ci._metals_in(d, {"input_file": "input.txt"}) == ["Ru"]
    monkeypatch.setattr(ci, "_GEOMETRY_MAX_BYTES", 8, raising=False)
    assert ci._metals_in(d, {"input_file": "input.txt"}) is None


def test_a_heavy_metal_in_the_smiles_counts(tmp_path):
    rec, _ = _record(tmp_path, "ru_smiles",
                     _template_control(SMILES="Cl[Ru](Cl)(N)N"))
    assert rec["basis_set"] == "ZORA-def2-SVP"


def test_without_a_structure_the_basis_comes_from_the_inp(tmp_path):
    """The *_rel choice cannot be made without the atoms; the ORCA input
    states what ran."""
    rec, _ = _record(tmp_path, "no_geom", _template_control(),
                     run__inp="! PBE0 ZORA-def2-SVP D4\n* xyz 0 1\n*\n")
    assert rec["functional"] == "PBE0"
    assert rec["basis_set"] == "ZORA-def2-SVP"


def test_the_inp_header_states_what_ran_over_control(tmp_path):
    """DELFIN regenerates every .inp from CONTROL.txt at run time, so the
    two differ only when CONTROL was edited after the run (or the .inp by
    hand); the .inp is then what produced the outputs in the folder.
    Its functional, basis and dispersion win over CONTROL; CONTROL fills
    what the header does not state; DELFIN_Data.json wins over both."""
    inp = "! B3LYP def2-TZVP D3BJ RIJCOSX def2/J\n* xyz 0 1\n*\n"
    rec, _ = _record(tmp_path, "edited", _template_control(),
                     geometry=_TXT_ORGANIC, initial__inp=inp)
    assert (rec["functional"], rec["basis_set"], rec["dispersion"]) == (
        "B3LYP", "def2-TZVP", "D3BJ")
    parts = ci.method_parts_of_folder(tmp_path / "calc" / "edited")
    assert (parts["functional"], parts["basis"], parts["dispersion"]) == (
        "B3LYP", "def2-TZVP", "D3BJ")
    assert parts["solvent"] == "DMF"


def test_an_inp_without_a_dft_method_leaves_control_in_place(tmp_path):
    rec, _ = _record(tmp_path, "xtb_first", _template_control(),
                     geometry=_TXT_ORGANIC,
                     XTB__inp="! XTB2 OPT ALPB(DMF)\n* xyz 0 1\n*\n")
    assert (rec["functional"], rec["basis_set"], rec["dispersion"]) == (
        "PBE0", "def2-SVP", "D4")
    parts = ci.method_parts_of_folder(tmp_path / "calc" / "xtb_first")
    assert (parts["functional"], parts["basis"], parts["dispersion"]) == (
        "PBE0", "def2-SVP", "D4")


def test_an_xtb_inp_sorting_first_does_not_hide_the_dft_inp(tmp_path):
    """XTB.inp sorts before initial.inp; the index reads the header that
    names a functional, as method_parts_of_folder does."""
    name = "xtb_and_dft"
    rec, _ = _record(tmp_path, name, _template_control(),
                     geometry=_TXT_ORGANIC,
                     XTB__inp="! XTB2 OPT ALPB(DMF)\n* xyz 0 1\n*\n",
                     initial__inp="! B3LYP def2-TZVP D3BJ\n* xyz 0 1\n*\n")
    parts = ci.method_parts_of_folder(tmp_path / "calc" / name)
    assert (rec["functional"], rec["basis_set"], rec["dispersion"]) == (
        parts["functional"], parts["basis"], parts["dispersion"]) == (
        "B3LYP", "def2-TZVP", "D3BJ")


@pytest.mark.parametrize("method", [
    "B3LYP-D3", "B3LYP-D3BJ", "PBE0-D4", "PBE-D3BJ", "TPSSh-D3BJ",
    "B97M-V", "wB97X-D", "wB97M-V", "r2SCAN-3c", "B97-3c", "PBEh-3c",
    "HF-3c", "GFN1-xTB", "GFN2-xTB", "XTB", "XTB1", "XTB2",
])
def test_the_inp_header_names_the_method(method):
    rec = ci._extract_from_inp_header(f"! {method} OPT\n* xyz 0 1\n*\n")
    assert rec.get("inp_functional") == method


def test_an_xtb_header_does_not_override_the_control_functional(tmp_path):
    """XTB.inp is DELFIN's pre-optimisation; it names a method now, but
    the run's DFT functional stays the one CONTROL (or a DFT .inp)
    states."""
    rec, _ = _record(tmp_path, "xtb_only", _template_control(),
                     geometry=_TXT_ORGANIC,
                     XTB__inp="! XTB2 OPT ALPB(DMF)\n* xyz 0 1\n*\n")
    parts = ci.method_parts_of_folder(tmp_path / "calc" / "xtb_only")
    assert rec["functional"] == parts["functional"] == "PBE0"
    rec, _ = _record(tmp_path, "xtb_bare", "NAME = x\n",
                     XTB__inp="! XTB2 OPT ALPB(DMF)\n* xyz 0 1\n*\n")
    assert rec["functional"] == "XTB2"


def test_delfin_data_json_wins_over_the_inp_header(tmp_path):
    rec, _ = _record(
        tmp_path, "json_inp", _template_control(), geometry=_TXT_ORGANIC,
        initial__inp="! B3LYP def2-TZVP D3BJ\n* xyz 0 1\n*\n",
        DELFIN_Data__json=json.dumps({"metadata": {
            "functional": "PBE0", "basis_set": "def2-SVP",
            "dispersion_correction": "D4"}}))
    assert (rec["functional"], rec["basis_set"], rec["dispersion"]) == (
        "PBE0", "def2-SVP", "D4")


def test_a_geometry_outside_the_folder_is_not_read(tmp_path):
    """input_file naming a path out of the calc folder (or a symlink out
    of it) is not opened; the basis then comes from nowhere."""
    outside = tmp_path / "outside.xyz"
    outside.write_text(_XYZ_RU, encoding="utf-8")
    rec, _ = _record(tmp_path, "escape",
                     _template_control(input_file="../../outside.xyz"))
    assert rec["basis_set"] == ""
    d = tmp_path / "calc" / "escape"
    (d / "input.txt").symlink_to(outside)
    assert ci._metals_in(d, {"input_file": "input.txt"}) is None
    (d / "input.txt").unlink()
    (d / "input.txt").write_text(_XYZ_RU, encoding="utf-8")
    assert ci._metals_in(d, {"input_file": "input.txt"}) == ["Ru"]


def test_the_legacy_basis_set_key_is_still_read(tmp_path):
    rec, _ = _record(tmp_path, "legacy",
                     "NAME = x\nfunctional = B3LYP\nbasis_set = def2-TZVP\n")
    assert rec["functional"] == "B3LYP"
    assert rec["basis_set"] == "def2-TZVP"


def test_delfin_data_json_still_wins_over_control(tmp_path):
    rec, _ = _record(
        tmp_path, "both", _template_control(), geometry=_XYZ_ORGANIC,
        DELFIN_Data__json=json.dumps({"metadata": {
            "functional": "B3LYP", "basis_set": "def2-TZVP",
            "dispersion_correction": "D3BJ"}}))
    assert rec["index_source"] == "DELFIN_Data.json"
    assert rec["functional"] == "B3LYP"
    assert rec["basis_set"] == "def2-TZVP"
    assert rec["dispersion"] == "D3BJ"


def test_method_parts_of_folder_reads_the_template(tmp_path):
    d = _calc(tmp_path, "parts", _template_control(), geometry=_XYZ_ORGANIC)
    parts = ci.method_parts_of_folder(d)
    assert parts["functional"] == "PBE0"
    assert parts["basis"] == "def2-SVP"
    assert parts["dispersion"] == "D4"
    assert parts["solvent"] == "DMF"


def test_method_parts_keeps_the_old_dispersion_key_as_fallback(tmp_path):
    d = _calc(tmp_path, "old_disp", "functional = PBE0\ndispersion_correction = D3BJ\n")
    assert ci.method_parts_of_folder(d)["dispersion"] == "D3BJ"
    d = _calc(tmp_path, "both_disp",
              "disp_corr = D4\ndispersion_correction = D3BJ\n")
    assert ci.method_parts_of_folder(d)["dispersion"] == "D4"


def test_the_engine_finds_and_summarises_the_run(tmp_path):
    _, idx = _record(tmp_path, "organic", _template_control(),
                     geometry=_XYZ_ORGANIC)
    eng = CalcSearchEngine(idx)
    hits = eng.search(functional="PBE0", basis_set="def2-SVP")["results"]
    assert [h["calc_id"] for h in hits] == ["organic"]
    assert ("PBE0", 1) in eng.summary()["functionals"]


def test_the_agent_tools_return_the_method(tmp_path, monkeypatch):
    """Through the executor the agent calls, not the engine alone."""
    calc = tmp_path / "calc"
    _calc(calc, "ru", _template_control(), geometry=_XYZ_RU)
    monkeypatch.setattr(A._DocToolExecutor, "_configured_calc_dirs",
                        staticmethod(dict))
    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    ex._calc_engine = None
    ex._calc_dirs = {"calc": str(calc), "archive": str(tmp_path / "none")}
    ex._calc_roots = {}
    ex._calc_evidence = []

    found = json.loads(ex._execute_calc(
        "search_calcs", {"functional": "PBE0", "basis_set": "ZORA-def2-SVP"}))
    assert [r["calc_id"] for r in found["results"]] == ["ru"]

    info = json.loads(ex._execute_calc("get_calc_info", {"calc_id": "ru"}))
    assert info["functional"] == "PBE0"
    assert info["basis_set"] == "ZORA-def2-SVP"
    assert info["dispersion"] == "D4"
    assert info["metal_basis"] == "SARC-ZORA-TZVP"
    assert info["relativity"] == "ZORA"
