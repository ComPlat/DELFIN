"""Tests for the fragment reference (reaction energy evaluation).

Mock-only: no ORCA runs. The splitting logic and the reaction-energy
arithmetic are checked against fake scan/fragment energies.
"""
import json
from pathlib import Path

import pytest
from ase import Atoms
from ase.io import write as ase_write

from delfin.proton.chain_setup import bookkeeping
from delfin.proton.fragments import run_fragment_reference


@pytest.fixture
def scan_end_frame(tmp_path):
    # Complex = 1 Cu atom; donor = 3 atoms with H as the LAST one
    # (donor_h_global = 3 in the combined frame 0..3).
    atoms = Atoms("CuCOH", positions=[[0., 0., 0.], [2., 0., 0.],
                                      [1., 2., 0.], [1.5, 0.5, 0.]])
    p = tmp_path / "complex_with_donor.xyz"
    ase_write(str(p), atoms)
    return str(p)


def test_fragment_split_and_inputs(tmp_path, scan_end_frame, monkeypatch):
    monkeypatch.chdir(tmp_path)
    bk = bookkeeping("tfe", complex_charge=-1, complex_multiplicity=2)
    reaction = run_fragment_reference(
        combined_path=scan_end_frame, complex_n_atoms=1,
        donor_h_global=3, bookkeeping=bk, keywords="OPT PBE0 def2-SVP",
        PAL=4, maxcore=1000)
    assert reaction["status"] == "incomplete"
    # [M-H] = Cu + transferred H (2 atoms); donor product = C,O (2 atoms)
    hyd = Atoms()
    from ase.io import read as ase_read
    hyd = ase_read(str(tmp_path / "fragment_reference" / "metal_hydride.xyz"))
    don = ase_read(str(tmp_path / "fragment_reference" / "tfe_alkoxide.xyz"))
    assert sorted(hyd.get_chemical_symbols()) == ["Cu", "H"]
    assert sorted(don.get_chemical_symbols()) == ["C", "O"]
    # charges/multiplicities land in the ORCA inputs
    assert Path(str(tmp_path / "fragment_reference" / "metal_hydride.inp"))\
        .read_text().startswith("! OPT PBE0 def2-SVP")
    # the metal hydride inherits the combined spin state (mult 2 here)
    assert "* xyz 0 2" in Path(
        str(tmp_path / "fragment_reference" / "metal_hydride.inp")).read_text()
    # reactions.json written
    data = json.loads((tmp_path / "fragment_reference" / "reactions.json")
                      .read_text())
    assert data["status"] == "incomplete"


def test_fragment_deltae_complete(tmp_path, scan_end_frame, monkeypatch):
    monkeypatch.chdir(tmp_path)
    bk = bookkeeping("tfe", complex_charge=-1, complex_multiplicity=2)
    # Fake ORCA outputs: scan end + the two fragments
    Path("relaxed_surface_scan").mkdir()
    Path("relaxed_surface_scan/scan.relaxscanact.out").write_text(
        "Final Single Point Energy   -100.000000\n")
    Path("fragment_reference").mkdir()
    Path("fragment_reference/metal_hydride.out").write_text(
        "Final Single Point Energy   -49.700000\n")
    Path("fragment_reference/tfe_alkoxide.out").write_text(
        "Final Single Point Energy   -49.900000\n")
    reaction = run_fragment_reference(
        combined_path=scan_end_frame, complex_n_atoms=1,
        donor_h_global=3, bookkeeping=bk, keywords="OPT PBE0 def2-SVP",
        PAL=4, maxcore=1000)
    assert reaction["status"] == "complete"
    # (-49.7 + -49.9) - (-100.0) = 0.4 Eh = 251.0 kcal/mol
    assert reaction["deltaE_products_minus_scanend_Eh"] == pytest.approx(0.4)
    assert reaction["deltaE_products_minus_scanend_kcalmol"] == \
        pytest.approx(0.4 * 627.509474)
