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
        with pytest.raises(ValueError, match="substrate_atom"):
            adduct_flow.run_adduct_flow(str(d), workdir=str(d))

    def test_general_substrate_co_end_on_flow(self, tmp_path):
        d = tmp_path / "CO_coordination"
        d.mkdir()
        # CO molecule: C and O
        co_xyz = d / "co.xyz"
        co_xyz.write_text("2\nCO\nC 0 0 0\nO 0 0 1.13\n")
        # Bare metal complex
        comp_xyz = d / "complex.xyz"
        comp_xyz.write_text("1\nFe\nFe 0 0 0\n")
        # Placed geometry
        _write_adduct_xyz(d / "complex_aligned_with_CO.xyz", m_c_dist=1.9)
        write_default_files(str(d / "CONTROL.txt"), str(co_xyz))
        text = (d / "CONTROL.txt").read_text()
        text += (
            "\nadduct_flow=true\nadduct_start_xyz=complex_aligned_with_CO.xyz\n"
            "substrate_atom=C\ncoord_max_dist=3.0\nrun_xtb=false\nmode=atom-on:C\n"
        )
        (d / "CONTROL.txt").write_text(text)
        result = adduct_flow.run_adduct_flow(str(d), workdir=str(d))
        assert result["status"] == "coordinated"
        assert result["metal_substrate_distance_A"] == pytest.approx(1.9, abs=0.01)

    def test_substrate_atom_symbol_resolves_anchor(self, tmp_path):
        d = tmp_path / "CO2_coordination_sym"
        d.mkdir()
        _write_adduct_xyz(d / "complex_aligned_with_CO2.xyz", m_c_dist=2.1)
        write_default_files(str(d / "CONTROL.txt"), str(d / "co2.xyz"))
        text = (d / "CONTROL.txt").read_text()
        text += (
            "\nadduct_flow=true\nadduct_start_xyz=complex_aligned_with_CO2.xyz\n"
            "substrate_atom=C\ncoord_max_dist=3.0\nrun_xtb=false\n"
        )
        (d / "CONTROL.txt").write_text(text)
        result = adduct_flow.run_adduct_flow(str(d), workdir=str(d))
        assert result["status"] == "coordinated"
        assert result["metal_substrate_distance_A"] == pytest.approx(2.1, abs=1e-4)

    def test_automatic_inverse_rss_triggered_after_occupier(self, tmp_path, monkeypatch):
        d = _make_coordinator_dir(
            tmp_path,
            extra="run_occupier=true\nscan_dissoc=true\ndissoc_distance=4.0\nscan_steps=10\n"
        )
        api_calls = []
        rss_calls = []

        def fake_api_run(control_file="CONTROL.txt", **kwargs):
            api_calls.append(control_file)
            job_dir = os.path.dirname(control_file)
            occ = os.path.join(job_dir, "adduct_OCCUPIER", "opt")
            os.makedirs(occ, exist_ok=True)
            with open(os.path.join(occ, "optimized.xyz"), "w") as f:
                f.write("3\nopt\nNi 0 0 0\nC 0 0 1.95\nO 0 0 3.11\n")
            return 0

        def fake_write_orca_input_and_run(atoms, xyz_path, metal_index, co2_c_index,
                                         start_distance, end_distance, steps, **kwargs):
            rss_calls.append({
                "xyz_path": xyz_path,
                "metal_index": metal_index,
                "co2_c_index": co2_c_index,
                "start_distance": start_distance,
                "end_distance": end_distance,
                "steps": steps,
            })
            scan_dir = os.path.join(os.path.dirname(xyz_path) or ".", "relaxed_surface_scan")
            os.makedirs(scan_dir, exist_ok=True)
            with open(os.path.join(scan_dir, "scan.relaxscanact.dat"), "w") as f:
                f.write("# dummy scan data\n1.95 -100.0\n4.00 -99.8\n")

        import delfin.api as api_mod
        from delfin.co2 import CO2_Coordinator6 as coord_mod
        monkeypatch.setattr(api_mod, "run", fake_api_run, raising=False)
        monkeypatch.setattr(coord_mod, "write_orca_input_and_run", fake_write_orca_input_and_run)

        result = adduct_flow.run_adduct_flow(str(d), workdir=str(d))

        assert len(api_calls) == 1
        assert len(rss_calls) == 1
        assert rss_calls[0]["end_distance"] == pytest.approx(4.0)
        assert rss_calls[0]["start_distance"] == pytest.approx(1.95, abs=1e-3)
        assert rss_calls[0]["steps"] == 10
        assert result["inverse_rss"]["status"] == "ok"
        assert result["inverse_rss"]["start_distance"] == pytest.approx(1.95, abs=1e-3)
        assert result["inverse_rss"]["end_distance"] == pytest.approx(4.0)

    def test_coordination_distance_is_reported(self, tmp_path):
        d = _make_coordinator_dir(tmp_path, m_c_dist=2.5)
        result = adduct_flow.run_adduct_flow(str(d), workdir=str(d))

        assert result["status"] == "coordinated"
        assert result["metal_substrate_distance_A"] == pytest.approx(2.5, abs=1e-6)

    def test_run_occupier_true_triggers_pipeline(self, tmp_path, monkeypatch):
        d = _make_coordinator_dir(tmp_path, extra="run_occupier=true\n")
        calls = []

        def fake_api_run(control_file="CONTROL.txt", **kwargs):
            calls.append(control_file)
            # simulate a successful OCCUPIER pass writing a result folder
            job_dir = os.path.dirname(control_file)
            occ = os.path.join(job_dir, "adduct_OCCUPIER", "opt")
            os.makedirs(occ, exist_ok=True)
            with open(os.path.join(occ, "optimized.xyz"), "w") as f:
                f.write("3\nopt\nNi 0 0 0\nC 0 0 2.0\nO 0 0 3.2\n")
            return 0

        import delfin.api as api_mod
        monkeypatch.setattr(api_mod, "run", fake_api_run, raising=False)
        # adduct_flow imports api lazily inside _run_occupier_pass, so
        # patching the module attribute is sufficient.
        result = adduct_flow.run_adduct_flow(str(d), workdir=str(d))

        assert len(calls) == 1
        assert calls[0].endswith("CONTROL.txt")
        assert result["occupier_run"]["status"] == "ok"
        assert result["occupier_opt_xyz"].endswith("optimized.xyz")

    def test_run_occupier_failure_is_reported(self, tmp_path, monkeypatch):
        d = _make_coordinator_dir(tmp_path, extra="run_occupier=true\n")

        def fake_api_run(control_file="CONTROL.txt", **kwargs):
            return 2

        import delfin.api as api_mod
        monkeypatch.setattr(api_mod, "run", fake_api_run, raising=False)
        result = adduct_flow.run_adduct_flow(str(d), workdir=str(d))

        assert result["status"] == "coordinated"  # chain itself succeeded
        assert result["occupier_run"]["status"] == "failed"
        assert "exit code 2" in result["occupier_run"]["error"]


class TestOccupierJobPrep:
    def test_control_forwards_level_of_theory(self, tmp_path):
        d = _make_coordinator_dir(tmp_path, extra="functional=PBE0\nbasisset=def2-SVP\n")
        result = adduct_flow.run_adduct_flow(str(d), workdir=str(d))
        ctrl = (d / result["occupier_job_path"] / "CONTROL.txt").read_text()
        assert "functional=PBE0" in ctrl
        assert "basis=def2-SVP" in ctrl
