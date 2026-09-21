"""Tests for delfin/co2/adduct_flow.py — the automatic coordination chain.

All tests are dry runs: run_xtb stays off, so no ORCA/xTB is ever executed.
"""
import json
import os

import numpy as np
import pytest

from delfin.co2 import adduct_flow
from delfin.co2.CO2_Coordinator6 import write_default_files


def _write_adduct_xyz(path, metal="Ni", m_c_dist=2.0, substrate_anchor="C"):
    """Simple linear M...C(x) geometry with the given M-C distance."""
    d = m_c_dist
    lines = [
        "3",
        f"adduct {metal}-CO2 test geometry",
        f"{metal}  0.0 0.0 0.0",
        f"{substrate_anchor}  0.0 0.0 {d}",
        "O  0.0 0.0 {:.4f}".format(d + 1.16),
    ]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _make_coordinator_dir(tmp_path, m_c_dist=2.0, coord_max_dist="3.0",
                          substrate_atom_index="1", extra=None):
    d = tmp_path / "CO2_coordination"
    d.mkdir()
    _write_adduct_xyz(d / "complex_aligned_with_CO2.xyz", m_c_dist=m_c_dist)
    write_default_files(str(d / "CONTROL.txt"), str(d / "co2.xyz"))
    # enable the flow keys on top of the default template
    text = (d / "CONTROL.txt").read_text()
    text += (
        f"\nadduct_flow=true\nadduct_start_xyz=complex_aligned_with_CO2.xyz\n"
        f"substrate_atom_index={substrate_atom_index}\n"
        f"coord_max_dist={coord_max_dist}\nrun_xtb=false\n"
    )
    if extra:
        text += extra
    (d / "CONTROL.txt").write_text(text)
    return d


class TestRunAdductFlow:
    def test_coordinated_adduct_prepares_occupier_job(self, tmp_path):
        d = _make_coordinator_dir(tmp_path, m_c_dist=2.0)
        result = adduct_flow.run_adduct_flow(str(d), workdir=str(d))

        assert result["status"] == "coordinated"
        # result JSON written
        assert (d / "adduct_flow_result.json").exists()
        # OCCUPIER job dir prepared but NOT executed
        job = d / "occupier_job"
        assert job.is_dir()
        assert (job / "input.xyz").exists()
        assert (job / "input.txt").exists()
        ctrl = (job / "CONTROL.txt").read_text()
        assert "method=OCCUPIER" in ctrl
        assert "enable_auto_recovery=yes" in ctrl
        assert "max_recovery_attempts=3" in ctrl

    def test_distant_substrate_reports_no_coordination(self, tmp_path):
        d = _make_coordinator_dir(tmp_path, m_c_dist=4.0)
        result = adduct_flow.run_adduct_flow(str(d), workdir=str(d))

        assert result["status"] == "no_coordination"
        # no OCCUPIER job for a failed coordination test (key stays None)
        assert result.get("occupier_job_path") is None
        assert not (d / "occupier_job").exists()

    def test_missing_placement_geometry_raises(self, tmp_path):
        d = _make_coordinator_dir(tmp_path)
        os.remove(d / "complex_aligned_with_CO2.xyz")
        with pytest.raises(FileNotFoundError):
            adduct_flow.run_adduct_flow(str(d), workdir=str(d))

    def test_missing_substrate_index_raises(self, tmp_path):
        d = _make_coordinator_dir(tmp_path, substrate_atom_index="")
        with pytest.raises(ValueError, match="substrate_atom_index"):
            adduct_flow.run_adduct_flow(str(d), workdir=str(d))

    def test_coordination_distance_is_reported(self, tmp_path):
        d = _make_coordinator_dir(tmp_path, m_c_dist=2.5)
        result = adduct_flow.run_adduct_flow(str(d), workdir=str(d))

        assert result["status"] == "coordinated"
        assert result["metal_substrate_distance_A"] == pytest.approx(2.5, abs=1e-6)


class TestOccupierJobPrep:
    def test_control_forwards_level_of_theory(self, tmp_path):
        d = _make_coordinator_dir(tmp_path, extra="functional=PBE0\nbasisset=def2-SVP\n")
        result = adduct_flow.run_adduct_flow(str(d), workdir=str(d))
        ctrl = (d / result["occupier_job_path"] / "CONTROL.txt").read_text()
        assert "functional=PBE0" in ctrl
        assert "basis=def2-SVP" in ctrl
