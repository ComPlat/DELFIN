"""Tests for the proton coordinator (HER competition study).

All tests are mock-based: no real ORCA, xTB or SLURM job is ever launched.
They protect the chemistry contract (donor templates, charge/spin
bookkeeping, per-donor folders) and the CONTROL wiring.
"""
import os
import textwrap
from pathlib import Path

import pytest

from delfin.proton.proton_coordinator import (
    CANONICAL_DONOR_DIRS,
    DONOR_TEMPLATES,
    _run_for_donor,
    build_donor_atoms,
    build_donor_opt_input,
    combine_complex_donor,
    detect_metal_index,
    orient_donor_toward_metal,
    read_xyz_robust,
    resolve_donor_name,
    rotate_donor_about_axis,
    write_orca_scan_input,
    write_orca_sp_input,
)
from delfin.proton.chain_setup import (
    bookkeeping,
    species_delta_to_name,
    write_default_control,
)

import numpy as np


# ---------------------------------------------------------------------------
# Donor templates
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("donor_key", ["tfe", "bih", "phenol", "hoac", "h3o"])
def test_donor_template_atom_counts_match_xyz(donor_key):
    """Every template's declared n_atoms matches its xyz block."""
    tpl = DONOR_TEMPLATES[donor_key]
    atoms = build_donor_atoms(donor_key)
    assert len(atoms) == tpl["n_atoms"]


@pytest.mark.parametrize("donor_key", ["tfe", "bih", "phenol", "hoac", "h3o"])
def test_donor_h_is_last_atom_and_bonded(donor_key):
    """donor_h_index points at an H whose neighbour list contains donor_heavy_index."""
    tpl = DONOR_TEMPLATES[donor_key]
    atoms = build_donor_atoms(donor_key)
    symbols = atoms.get_chemical_symbols()
    assert symbols[tpl["donor_h_index"]] == "H"
    dist = np.linalg.norm(
        atoms.positions[tpl["donor_h_index"]]
        - atoms.positions[tpl["donor_heavy_index"]])
    # X-H bond lengths: C-H ~1.1, O-H ~1.0; allow a sketch-geometry margin
    assert 0.7 < dist < 1.4


def test_donor_aliases_resolve():
    assert resolve_donor_name("H3O+") == "h3o"
    assert resolve_donor_name("trifluorethanol") == "tfe"
    assert resolve_donor_name("acetic_acid") == "hoac"
    assert resolve_donor_name("something_else") is None


def test_donor_folders_are_per_donor():
    """The user requested one folder per donor, not a shared one."""
    dirs = set(CANONICAL_DONOR_DIRS.values())
    assert dirs == {"TFE_coordinator", "BIH_coordinator", "phenol_coordinator",
                    "HOAc_coordinator", "H3O+_coordinator"}


# ---------------------------------------------------------------------------
# Bookkeeping: charge/spin conservation incl. radical donors
# ---------------------------------------------------------------------------

def test_bookkeeping_protic_tfe_conserves_charge():
    bk = bookkeeping("tfe", complex_charge=-2, complex_multiplicity=1)
    assert bk["combined_charge"] == -2          # neutral donor added
    assert bk["donor_product"]["charge"] == -1  # TFE-O(-) conjugate base
    # the complex GAINS the proton: [M-H] charge = complex_charge + 1
    assert bk["complex_product"]["charge"] == -1
    # conservation: fragments sum == combined
    assert (bk["complex_product"]["charge"]
            + bk["donor_product"]["charge"]) == bk["combined_charge"]


def test_bookkeeping_h3o_protic_products():
    # H3O+ carries its own +1: combined = complex + 1, and the proton that
    # leaves goes to the complex (+1), while the residual water cluster is 0.
    bk = bookkeeping("h3o", complex_charge=1, complex_multiplicity=2)
    assert bk["combined_charge"] == 2             # complex(+1) + H3O+(+1)
    assert bk["donor_product"]["charge"] == 0     # residual water cluster
    assert bk["donor_product"]["multiplicity"] == 1
    assert bk["complex_product"]["charge"] == 2   # complex + H+
    assert (bk["complex_product"]["charge"]
            + bk["donor_product"]["charge"]) == bk["combined_charge"]


def test_bookkeeping_bih_hat_gives_radical_doublet():
    bk = bookkeeping("bih", complex_charge=0, complex_multiplicity=1)
    assert bk["scheme"] == "hat"
    assert bk["donor_product"]["name"] == "bi_radical"
    assert bk["donor_product"]["multiplicity"] == 2  # open-shell BI radical
    assert bk["donor_product"]["charge"] == 0


def test_bookkeeping_phenol_hat_gives_phenoxyl_radical():
    bk = bookkeeping("phenol", complex_charge=0, complex_multiplicity=1)
    assert bk["scheme"] == "hat"
    assert bk["donor_product"]["name"] == "phenoxyl_radical"
    assert bk["donor_product"]["multiplicity"] == 2


def test_bookkeeping_unknown_donor_raises():
    with pytest.raises(KeyError):
        bookkeeping("not_a_donor")


# ---------------------------------------------------------------------------
# Geometry: placement and rotation
# ---------------------------------------------------------------------------

def _tiny_complex():
    """Fake metal complex: one metal atom as the whole 'complex'."""
    return read_xyz_robust(str(_tiny_complex_xyz()))


def _tiny_complex_xyz(tmp_path=None):
    """Write a 1-atom 'complex' (a metal ion) to a temp file."""
    import tempfile
    d = tmp_path or Path(tempfile.mkdtemp())
    p = d / "complex.xyz"
    p.write_text("1\n\nCu 0.0 0.0 0.0\n")
    return p


def test_orient_donor_places_h_at_distance_toward_metal():
    donor = build_donor_atoms("tfe")
    metal_pos = np.array([0.0, 0.0, 0.0])
    placed = orient_donor_toward_metal(donor, donor_h_index=6,
                                       donor_heavy_index=5,
                                       metal_pos=metal_pos, distance=3.0)
    dist = np.linalg.norm(placed.positions[6] - metal_pos)
    assert abs(dist - 3.0) < 1e-6


def test_combine_complex_donor_offset():
    complex_atoms = _tiny_complex()
    donor = build_donor_atoms("tfe")
    combined, h_idx = combine_complex_donor(complex_atoms, donor)
    assert h_idx == len(complex_atoms) + len(donor) - 1
    assert len(combined) == len(complex_atoms) + len(donor)
    assert combined.get_chemical_symbols()[h_idx] == "H"


def test_rotate_donor_moves_only_donor_atoms():
    complex_atoms = _tiny_complex()
    donor = build_donor_atoms("tfe")
    combined, _ = combine_complex_donor(complex_atoms, donor)
    donor_start = len(complex_atoms)
    before = combined.positions.copy()
    rot = rotate_donor_about_axis(combined, donor_start, 90.0,
                                  axis_pos=combined.positions[0],
                                  axis_dir=np.array([0.0, 0.0, 1.0]))
    # complex atoms untouched
    assert np.allclose(rot.positions[:donor_start],
                       before[:donor_start])
    # donor atoms moved (axis through origin, donor offset from origin)
    assert not np.allclose(rot.positions[donor_start:],
                           before[donor_start:])


def test_detect_metal_index_finds_cu():
    atoms = _tiny_complex()
    assert detect_metal_index(atoms) == 0


# ---------------------------------------------------------------------------
# ORCA input writers (content only, never executed)
# ---------------------------------------------------------------------------

# Donor pre-optimisation (sketch geometry -> OPT before scan)
# ---------------------------------------------------------------------------

def test_donor_opt_input_charge_and_multiplicity(tmp_path):
    # H3O+ carries its own +1; all other donors are neutral ground states.
    for donor_key, expect_charge, expect_mult in [
        ("tfe", 0, 1), ("bih", 0, 1), ("phenol", 0, 1),
        ("hoac", 0, 1), ("h3o", 1, 1),
    ]:
        inp = build_donor_opt_input(donor_key, str(tmp_path),
                                    "OPT PBE0 def2-SVP",
                                    charge=DONOR_TEMPLATES[donor_key]
                                    .get("donor_charge", 0),
                                    multiplicity=DONOR_TEMPLATES[donor_key]
                                    .get("donor_multiplicity", 1),
                                    PAL=4, maxcore=1000)
        text = Path(inp).read_text()
        assert text.startswith("! OPT PBE0 def2-SVP\n"), donor_key
        assert f"* xyz {expect_charge} {expect_mult}" in text, donor_key


def test_run_for_donor_reuses_preopt_xyz(tmp_path, monkeypatch, capsys):
    # A completed donor OPT must be picked up; without it the scan runs from
    # the sketch with a warning (not validated chemistry).
    monkeypatch.chdir(tmp_path)
    p = _tiny_complex_xyz(tmp_path)
    complex_xyz = str(p)
    donor_tpl = DONOR_TEMPLATES["tfe"]
    donor_atoms = build_donor_atoms("tfe")
    complex_atoms = read_xyz_robust(complex_xyz)
    metal_idx = detect_metal_index(complex_atoms)

    calls = {}

    def fake_scan(*args, **kwargs):
        from ase.io import write as ase_write
        ase_write("scan_best.xyz", args[0])
        calls["atoms"] = args[0]
        return "scan_best.xyz", 0.0, -1.0

    monkeypatch.setattr("delfin.proton.proton_coordinator.orientation_scan",
                        fake_scan)
    monkeypatch.setattr("delfin.proton.proton_coordinator.relaxed_distance_scan",
                        lambda *a, **k: None)

    args = {"xyz": "complex.xyz", "donor": "tfe", "out": "aligned.xyz",
            "orientation_distance": 3.0, "rot_step_deg": 30,
            "rot_range_deg": 30, "scan_end": 2.0, "scan_steps": 3,
            "keywords": "SP PBE0 def2-SVP", "PAL": 4, "maxcore": 1000,
            "charge": -1, "multiplicity": 2}

    # Run 1: no preopt.xyz -> OPT input written, sketch used with warning
    _run_for_donor(args, "tfe", xyz_in_abs=complex_xyz,
                   charge=-1, multiplicity=2, PAL=4, maxcore=1000,
                   keywords="SP PBE0 def2-SVP",
                   orientation_distance=3.0, rot_step_deg=30,
                   rot_range_deg=30, scan_end=2.0, scan_steps=3,
                   control_path=str(tmp_path / "CONTROL.txt"))
    out1 = capsys.readouterr().out
    assert (tmp_path / "tfe_preopt.inp").exists()
    assert "[WARN]" in out1 and "SKETCH" in out1

    # Run 2: preopt.xyz present -> used without warning
    # (write a rotated donor copy as the fake "optimised" geometry)
    preopt = donor_atoms.copy()
    preopt.positions[0] += [0.1, 0.0, 0.0]
    from ase.io import write as ase_write
    ase_write(str(tmp_path / "tfe_preopt.xyz"), preopt)
    _run_for_donor(args, "tfe", xyz_in_abs=complex_xyz,
                   charge=-1, multiplicity=2, PAL=4, maxcore=1000,
                   keywords="SP PBE0 def2-SVP",
                   orientation_distance=3.0, rot_step_deg=30,
                   rot_range_deg=30, scan_end=2.0, scan_steps=3,
                   control_path=str(tmp_path / "CONTROL.txt"))
    out2 = capsys.readouterr().out
    assert "[INFO]" in out2 and "pre-optimised" in out2
    assert "SKETCH" not in out2   # only the unrelated plot warning may remain


def test_write_orca_sp_input_content(tmp_path):
    atoms = _tiny_complex()
    donor = build_donor_atoms("tfe")
    combined, _ = combine_complex_donor(atoms, donor)
    inp = write_orca_sp_input(combined, str(tmp_path), "angle_000",
                              "SP PBE0 def2-SVP", charge=-1, multiplicity=2,
                              PAL=4, maxcore=1000,
                              broken_sym="%scf rot 1 1 90 end\n")
    text = Path(inp).read_text()
    assert text.startswith("! SP PBE0 def2-SVP\n")
    assert "%pal nprocs 4 end" in text
    assert "* xyz -1 2" in text          # combined charge/mult line
    assert text.rstrip().endswith("*")
    assert "Cu" in text


def test_write_orca_scan_input_scan_block(tmp_path):
    atoms = _tiny_complex()
    donor = build_donor_atoms("tfe")
    combined, h_idx = combine_complex_donor(atoms, donor)
    inp = write_orca_scan_input(combined, str(tmp_path), metal_index=0,
                                donor_h_index=h_idx,
                                start_distance=3.0, end_distance=1.5,
                                steps=10, keywords="OPT PBE0 def2-SVP",
                                charge=-1, multiplicity=2)
    text = Path(inp).read_text()
    assert "B  0  %d = 3.00, 1.5, 10" % h_idx in text
    assert "%geom" in text and "Scan" in text


# ---------------------------------------------------------------------------
# CONTROL wiring
# ---------------------------------------------------------------------------

def test_write_default_control_lists_donors(tmp_path):
    target = tmp_path / "CONTROL.txt"
    write_default_control(str(target))
    text = target.read_text()
    assert "proton_donor=tfe,phenol,bih,h3o" in text
    assert "implicit_solvation_model=CPCM" in text
    # never overwrite without --force
    write_default_control(str(target))
    assert text == target.read_text()


def test_species_delta_naming_matches_co2_convention():
    assert species_delta_to_name(0) == "initial"
    assert species_delta_to_name(-2) == "red_step_2"
    assert species_delta_to_name(+1) == "ox_step_1"


def test_build_control_keys_roundtrip_through_minimal_parser(tmp_path):
    # Regression: build_control once wrote proton_-prefixed keys while the
    # coordinator's minimal parser reads them un-prefixed — charge/mult
    # silently fell back to defaults (charge=0, mult=1) in the recalc path.
    # The generated CONTROL must round-trip through _minimal_read_control_file.
    from delfin.proton.chain_setup import build_control
    from delfin.proton.proton_coordinator import _minimal_read_control_file

    parent = {"charge": "-2", "multiplicity": "3", "functional": "PBE0",
              "PAL": "48", "maxcore": "3800", "solvent": "MeCN",
              "implicit_solvation_model": "CPCM"}
    target = tmp_path / "TFE_coordinator" / "CONTROL.txt"
    build_control(parent, "tfe", control_path=target)

    parsed = _minimal_read_control_file(str(target))
    assert parsed["charge"] == "-2"       # NOT lost to a proton_ prefix
    assert parsed["multiplicity"] == "3"
    assert parsed["PAL"] == "48"
    assert parsed["functional"] == "PBE0"
    assert parsed["proton_donor"] == "tfe"


def test_control_validator_accepts_proton_keys():
    from delfin.common.control_validator import validate_control_config
    result = validate_control_config({
        "charge": -1,
        "method": "classic",
        "proton_coordination": "on",
        "proton_donor": "tfe,phenol,bih,h3o",
        "proton_species_delta": 0,
    })
    assert result.get("proton_coordination") == "on"
    assert result.get("proton_donor") == "tfe,phenol,bih,h3o"


def test_control_validator_defaults_off():
    from delfin.common.control_validator import validate_control_config
    result = validate_control_config({"charge": 0, "method": "classic"})
    assert result.get("proton_coordination", "off") == "off"


def test_config_alias_maps_proton_keys():
    from delfin.config import _CONTROL_KEY_ALIASES
    assert _CONTROL_KEY_ALIASES["protoncoordination"] == "proton_coordination"
    assert _CONTROL_KEY_ALIASES["protondonor"] == "proton_donor"
