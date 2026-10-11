"""Tests for delfin/co2/reduction_surface.py — the CO2 reduction surface chain.

All tests are dry runs: rs_run_xtb stays off and ORCA is never executed
(or only via monkeypatched stubs). No real ORCA/xTB/SLURM job is started.
"""
import json
import os

import numpy as np
import pytest

from delfin.co2 import reduction_surface as rs
from delfin.co2 import CO2_Coordinator6 as _coord


def _write_adduct_xyz(path, m_c_dist=2.0, with_h=False):
    """M···C(x)O2 placement geometry: metal first, substrate tail last.

    Atom order: M, 2 ligand atoms (to give the H donor candidates context),
    then the substrate O, C, O block (the same O,C,O sequence as co2.xyz).
    """
    d = m_c_dist
    lines = [
        "8",
        "complex + CO2 test geometry",
        "Ni  0.0 0.0 0.0",
        "C   0.0 1.9 0.0",
        "C   0.0 2.5 0.0",
        "O   0.0 0.0 {:.4f}".format(d + 1.16),   # substrate O (anchored)
        "C   0.0 0.0 {:.4f}".format(d + 2.16),
        "O   0.0 0.0 {:.4f}".format(d + 2.16 + 1.16),
        "H   1.4 0.0 3.0",                        # complex ligand H (donor 1)
        "H  -1.4 0.0 3.0",                        # complex ligand H (donor 2)
    ]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _make_coordinator_dir(tmp_path, extra=""):
    d = tmp_path / "CO2_coordination"
    d.mkdir()
    _write_adduct_xyz(d / "complex_aligned_with_CO2.xyz")
    _coord.write_default_files(str(d / "CONTROL.txt"), str(d / "co2.xyz"))
    text = (d / "CONTROL.txt").read_text()
    text += (
        "\nreduction_surface=true\nadduct_start_xyz=complex_aligned_with_CO2.xyz\n"
        "rs_run_xtb=false\nreduction_surface_ts=false\nrs_run_dft=off\n"
    )
    if extra:
        text += extra
    (d / "CONTROL.txt").write_text(text)
    return d


class TestStationDeltas:
    def test_default_deltas_are_derived(self):
        assert rs._station_deltas({}) == [0, 1, 2, 2]

    def test_rs_charges_override(self):
        ctrl = {"rs_charges": "0,1,1,1"}
        assert rs._station_deltas(ctrl) == [0, 1, 1, 1]

    def test_rs_charges_wrong_length_falls_back(self):
        ctrl = {"rs_charges": "0,1"}
        assert rs._station_deltas(ctrl) == [0, 1, 2, 2]

    def test_rs_charges_nonint_falls_back(self):
        ctrl = {"rs_charges": "a,b,c,d"}
        assert rs._station_deltas(ctrl) == [0, 1, 2, 2]


class TestSurgery:
    SUBSTRATE_START = 3  # Ni, C, C | O, C, O, H(donor)

    def _atoms(self, with_h=True, n_h=2):
        lines = ["Ni  0.0 0.0 0.0", "C  0.0 1.9 0.0", "C   0.0 2.5 0.0",
                 "O   0.0 0.0 3.16", "C   0.0 0.0 4.16", "O   0.0 0.0 5.32"]
        if with_h:
            lines.append("H   1.4 0.0 3.0")
            if n_h > 1:
                lines.append("H  -1.4 0.0 3.0")
        from ase import Atoms
        return Atoms([l.split()[0] for l in lines],
                     positions=[[float(x) for x in l.split()[1:4]] for l in lines])

    def test_co2h_adds_one_h_at_distal_substrate_oxygen(self):
        atoms = self._atoms()
        anchor = rs._pick_anchor_oxygen(atoms, self.SUBSTRATE_START, 0)
        assert anchor == 5  # the distal O (farthest from the metal at idx 0)
        out = rs._apply_station_surgery(atoms, 1, self.SUBSTRATE_START,
                                        anchor_idx=anchor)
        assert len(out) == len(atoms)  # re-bound, not appended
        symbols = out.get_chemical_symbols()
        # tail keeps the O,C,O block, the donor H follows AFTER it
        assert symbols[3:6] == ["O", "C", "O"]
        assert symbols[6] == "H"
        # and sits at ~1.0 A from the DISTAL O (idx 5, farthest from the
        # metal) — the one that later leaves as H2O. Protonating the
        # proximal O (idx 3) would leave a hydroxamate, not a CO ligand.
        d = np.linalg.norm(out.positions[6] - out.positions[5])
        assert d == pytest.approx(1.00, abs=0.05)
        # the proximal O (idx 3) stays un-protonated
        d_prox = np.linalg.norm(out.positions[6] - out.positions[3])
        assert d_prox > 1.5

    def test_co2h2_adds_second_h_at_other_oxygen(self):
        atoms = self._atoms()
        anchor = rs._pick_anchor_oxygen(atoms, self.SUBSTRATE_START, 0)
        step1 = rs._apply_station_surgery(atoms, 1, self.SUBSTRATE_START,
                                          anchor_idx=anchor)
        # substrate start re-detected like the chain does (tail shifts)
        out = rs._apply_station_surgery(step1, 2, self.SUBSTRATE_START - 1,
                                        anchor_idx=anchor)
        symbols = out.get_chemical_symbols()
        # tail is O, C, O + two substrate H at the end
        assert symbols[3:6] == ["O", "C", "O"]
        assert symbols[6] == "H" and symbols[7] == "H"
        # user's mechanism: ONE H on EACH O — first H on the distal/anchor
        # O (idx 5), second H on the OTHER (proximal, metal-bound) O (idx 3)
        d_distal = np.linalg.norm(out.positions[6] - out.positions[5])
        d_prox = np.linalg.norm(out.positions[7] - out.positions[3])
        assert d_distal == pytest.approx(1.00, abs=0.05)
        assert d_prox == pytest.approx(1.00, abs=0.05)
        # and the second H is NOT on the anchor O
        d_wrong = np.linalg.norm(out.positions[7] - out.positions[5])
        assert d_wrong > 1.5

    def test_surgery_without_hydrogen_raises(self):
        atoms = self._atoms(with_h=False)
        anchor = rs._pick_anchor_oxygen(atoms, self.SUBSTRATE_START, 0)
        with pytest.raises(ValueError, match="no (unbound hydrogen donor|hydrogen)"):
            rs._apply_station_surgery(atoms, 1, self.SUBSTRATE_START,
                                      anchor_idx=anchor)


class TestSplitWater:
    def test_split_yields_body_and_water(self):
        # station_co2h2-like geometry: Ni, ligand C | O, C, O + 2 substrate H
        lines = [
            "Ni  0.0 0.0 0.0",
            "C  0.0 1.9 0.0",
            "O   0.0 0.0 3.16", "C   0.0 0.0 4.16", "O   0.0 0.0 5.32",
            "H   0.3 0.0 5.32",
            "H  -0.3 0.0 5.32",
        ]
        from ase import Atoms
        atoms = Atoms([l.split()[0] for l in lines],
                      positions=[[float(x) for x in l.split()[1:4]] for l in lines])
        body, water = rs._split_water(atoms, substrate_start=2)
        assert water is not None
        assert water.get_chemical_formula() == "H2O"
        # body: Ni + ligand C + substrate O + C (the bound CO)
        assert body.get_chemical_formula() == "C2NiO"
        assert len(body) == 4


class TestRunChain:
    def test_chain_runs_dry_and_writes_result(self, tmp_path):
        d = _make_coordinator_dir(tmp_path)
        result = rs.run_reduction_surface(str(d), workdir=str(d))

        assert result["status"] == "ok"
        names = [st["station"] for st in result["stations"]]
        assert names == ["co2", "co2h", "co2h2", "co"]
        # charges are cumulative
        assert [st["charge"] for st in result["stations"]] == [0, 1, 2, 2]
        # result JSON written
        assert (d / "reduction_surface_result.json").exists()
        # station folders with geometries
        for name in names:
            assert (d / "reduction_surface" / name / f"station_{name}.xyz").exists()
        # H2O-release station wrote body + water
        assert (d / "reduction_surface" / "co" / "station_co_body.xyz").exists()
        assert (d / "reduction_surface" / "co" / "station_h2o.xyz").exists()
        # the leaving water is exactly H2O, the body keeps the bound CO
        water = (d / "reduction_surface" / "co" / "station_h2o.xyz").read_text()
        assert "O" in water and water.count("H") >= 2
        # profile plot was written automatically
        assert (d / "reduction_surface" / "reduction_surface_profile.png").exists()
        assert (d / "reduction_surface" / "reduction_surface_profile.pdf").exists()
        assert result["profile_png"].endswith("reduction_surface_profile.png")

    def test_station_geometries_are_chained(self, tmp_path):
        d = _make_coordinator_dir(tmp_path)
        result = rs.run_reduction_surface(str(d), workdir=str(d))
        # co2h has one more H than co2; co2h2 two more; co keeps the same
        # heavy skeleton with the water split off
        def n_atoms(name):
            xyz = d / "reduction_surface" / name / f"station_{name}.xyz"
            return int(xyz.read_text().split()[0])
        assert n_atoms("co2h") == n_atoms("co2")
        assert n_atoms("co2h2") == n_atoms("co2")
        # the co station keeps the full co2h2 geometry (water still bound
        # for the NEB endpoint); the split body loses the H2O (3 atoms)
        body = d / "reduction_surface" / "co" / "station_co_body.xyz"
        assert int(body.read_text().split()[0]) == n_atoms("co2h2") - 3

    def test_missing_start_geometry_raises(self, tmp_path):
        d = _make_coordinator_dir(tmp_path)
        os.remove(d / "complex_aligned_with_CO2.xyz")
        with pytest.raises(FileNotFoundError):
            rs.run_reduction_surface(str(d), workdir=str(d))

    def test_station_failure_marks_incomplete(self, tmp_path, monkeypatch):
        d = _make_coordinator_dir(tmp_path)
        # break the second station's surgery
        orig = rs._apply_station_surgery

        def boom(atoms, idx, start, anchor_idx=None):
            if idx == 2:
                raise ValueError("no hydrogen found")
            return orig(atoms, idx, start, anchor_idx=anchor_idx)
        monkeypatch.setattr(rs, "_apply_station_surgery", boom)
        result = rs.run_reduction_surface(str(d), workdir=str(d))
        assert result["status"] == "incomplete"
        assert result["stations"][2]["status"] == "error"

    def test_custom_rs_charges_reach_the_stations(self, tmp_path):
        d = _make_coordinator_dir(tmp_path, extra="rs_charges=0,1,1,1\n")
        result = rs.run_reduction_surface(str(d), workdir=str(d))
        assert [st["charge"] for st in result["stations"]] == [0, 1, 1, 1]


class TestNebTs:
    def test_input_text_matches_orca_syntax(self):
        text = rs.build_neb_ts_input("reactant.xyz", "product.xyz", -2, 3,
                                     "PBE0 D4 RIJCOSX def2/J def2-SVP",
                                     nimages=8, maxcore=2000, pal=12)
        assert text.startswith("! NEB-TS PBE0 D4 RIJCOSX def2/J def2-SVP\n")
        assert '%neb\nNEB_End_XYZFile "product.xyz"\nNImages 8\nend\n' in text
        assert "*xyzfile -2 3 reactant.xyz\n" in text
        assert "%maxcore 2000\n" in text
        assert "%pal nprocs 12 end\n" in text

    def test_ts_runs_when_requested(self, tmp_path, monkeypatch):
        d = _make_coordinator_dir(tmp_path, extra="reduction_surface_ts=true\n")
        calls = []

        def fake_pair(reactant_xyz, product_xyz, charge, multiplicity,
                      control, workdir, step_name):
            calls.append(step_name)
            return {"step": step_name, "status": "ok",
                    "ts_energy_eh": -1.0,
                    "ts_energy_kcal_rel_reactant": None}

        monkeypatch.setattr(rs, "run_neb_ts_station_pair", fake_pair)
        # energies so the barrier arithmetic has something to act on
        result = rs.run_reduction_surface(str(d), workdir=str(d))
        assert calls == ["tsco2_co2h", "tsco2h_co2h2", "tsco2h2_co"]
        assert len(result["transition_states"]) == 3
        assert all(ts["status"] == "ok" for ts in result["transition_states"])

    def test_ts_off_by_default(self, tmp_path):
        d = _make_coordinator_dir(tmp_path)
        result = rs.run_reduction_surface(str(d), workdir=str(d))
        assert result["transition_states"] == []


class TestEnergyParsing:
    def test_final_single_point_energy_parsed(self, tmp_path):
        out = tmp_path / "x.out"
        out.write_text("blah\nFINAL SINGLE POINT ENERGY     -123.456789\n")
        assert rs._orca_output_energy(str(out)) == pytest.approx(-123.456789)

    def test_last_energy_wins(self, tmp_path):
        out = tmp_path / "x.out"
        out.write_text("FINAL SINGLE POINT ENERGY  -1.0\n"
                       "FINAL SINGLE POINT ENERGY  -2.5\n")
        assert rs._orca_output_energy(str(out)) == pytest.approx(-2.5)

    def test_missing_energy_is_none(self, tmp_path):
        out = tmp_path / "x.out"
        out.write_text("nothing here\n")
        assert rs._orca_output_energy(str(out)) is None


class TestPlot:
    def test_profile_plot_renders(self, tmp_path):
        import matplotlib
        matplotlib.use("Agg")
        from delfin.co2.reduction_surface_plot import write_profile_plot

        result = {
            "stations": [
                {"station": "co2", "label": "M···CO2", "status": "ok",
                 "energy_eh": -100.0, "energy_eh_xtb": None},
                {"station": "co2h", "label": "M···CO2H", "status": "ok",
                 "energy_eh": -100.01, "energy_eh_xtb": None},
                {"station": "co2h2", "label": "M···CO2H2", "status": "error",
                 "energy_eh": None, "energy_eh_xtb": None},
                {"station": "co", "label": "M···CO + H2O", "status": "ok",
                 "energy_eh": -100.02, "energy_eh_xtb": None},
            ],
            "transition_states": [
                {"step": "tsco2_co2h", "status": "ok", "ts_energy_eh": -99.99,
                 "ts_energy_kcal_rel_reactant": 6.275},
            ],
        }
        png, pdf = write_profile_plot(result, str(tmp_path))
        assert os.path.exists(png)
        assert os.path.exists(pdf)
        assert os.path.getsize(png) > 0

    def test_profile_plot_xtb_panel(self, tmp_path):
        import matplotlib
        matplotlib.use("Agg")
        from delfin.co2.reduction_surface_plot import write_profile_plot

        result = {
            "stations": [
                {"station": "co2", "label": "M···CO2", "status": "ok",
                 "energy_eh": -100.0, "energy_eh_xtb": -10.0},
                {"station": "co2h", "label": "M···CO2H", "status": "ok",
                 "energy_eh": -100.01, "energy_eh_xtb": -10.01},
            ],
            "transition_states": [],
        }
        png, pdf = write_profile_plot(result, str(tmp_path))
        assert os.path.exists(png) and os.path.getsize(png) > 0
