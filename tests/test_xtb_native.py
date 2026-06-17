"""Tests for the native xTB building blocks (the standalone ``xtb`` binary)."""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import pytest

from delfin.tools import list_steps, platform
from delfin.tools._registry import get
from delfin.tools.adapters import xtb_native as X

_HAS_XTB = shutil.which("xtb") is not None
_H2O = "3\n\nO 0.0 0.0 0.0\nH 0.0 0.0 0.9572\nH 0.9266 0.0 -0.2400\n"


# --- registration + contracts --------------------------------------------


def test_native_xtb_blocks_registered():
    caps = list_steps()
    for name in ("xtb_sp", "xtb_optimize", "xtb_hessian", "xtb_md", "xtb_fukui"):
        assert name in caps
    # the ORCA-driven blocks are untouched and still present
    assert "xtb_opt" in caps and "xtb_goat" in caps


def test_native_xtb_contracts_declare_full_surface():
    c = platform.describe_capability("xtb_hessian")
    assert c.requires_binaries == frozenset({"xtb"})
    assert {"hessian", "qc_output"} <= set(c.produces)
    assert {"energy_Eh", "free_energy_Eh", "zpve_Eh", "n_imaginary"} <= {
        d.name for d in c.data_keys}
    # method enum spans all xTB levels incl. the force field
    method = platform.describe_capability("xtb_sp").param("method")
    assert set(method.enum) == {"gfn0", "gfn1", "gfn2", "gfnff"}


# --- command construction (pure, no binary) -------------------------------


def test_cmd_method_flags_and_run_types():
    assert X.build_xtb_cmd("m.xyz", method="gfn0")[2:4] == ["--gfn", "0"]
    assert "--gfnff" in X.build_xtb_cmd("m.xyz", method="gfnff")
    opt_cmd = X.build_xtb_cmd("m.xyz", run="opt", opt_level="tight")
    assert "--opt" in opt_cmd and opt_cmd[opt_cmd.index("--opt") + 1] == "tight"
    assert "--ohess" in X.build_xtb_cmd("m.xyz", run="ohess")
    assert "--vfukui" in X.build_xtb_cmd("m.xyz", run="fukui")
    assert "--omd" in X.build_xtb_cmd("m.xyz", run="md")


def test_cmd_solvent_charge_and_extras():
    cmd = X.build_xtb_cmd("m.xyz", charge=-2, uhf=1, solvent="water",
                          solvent_model="gbsa", accuracy=0.5, etemp=400,
                          cores=4, extra_args=["--ceasefiles"])
    assert "--chrg" in cmd and cmd[cmd.index("--chrg") + 1] == "-2"
    assert "--uhf" in cmd and cmd[cmd.index("--uhf") + 1] == "1"
    assert cmd[cmd.index("--gbsa") + 1] == "water"
    assert "--acc" in cmd and "--etemp" in cmd
    assert cmd[cmd.index("-P") + 1] == "4"
    assert cmd[-1] == "--ceasefiles"          # raw escape hatch preserved


# --- output parsing (no binary) ------------------------------------------


def test_parsers_pull_energy_gap_thermo():
    sp = "| TOTAL ENERGY  -5.084832242975 Eh |\n| HOMO-LUMO GAP  15.166 eV |"
    assert X._parse_energy(sp) == pytest.approx(-5.084832242975)
    assert X._parse_gap(sp) == pytest.approx(15.166)
    thermo = (":: total free energy  -5.0680 Eh ::\n"
              ":: zero point energy  0.0201 Eh ::\n"
              "#  imaginary freq.   1 :")
    out = X._parse_thermo(thermo)
    assert out["free_energy_Eh"] == pytest.approx(-5.0680)
    assert out["zpve_Eh"] == pytest.approx(0.0201)
    assert out["n_imaginary"] == 1


# --- execution against the real binary (skipped if xtb absent) ------------


@pytest.mark.skipif(not _HAS_XTB, reason="xtb binary not installed")
def test_xtb_sp_optimize_hessian_execute():
    d = Path(tempfile.mkdtemp())
    geo = d / "h2o.xyz"
    geo.write_text(_H2O)

    def _wd(name: str) -> Path:
        wd = d / name
        wd.mkdir()
        return wd

    sp = get("xtb_sp").execute(_wd("sp"), geometry=geo, cores=2, charge=0, method="gfn2")
    assert sp.status.value == "success"
    assert sp.data["energy_Eh"] < 0 and "homo_lumo_gap_eV" in sp.data

    opt = get("xtb_optimize").execute(_wd("opt"), geometry=geo, cores=2, charge=0)
    assert opt.status.value == "success" and opt.geometry.name == "xtbopt.xyz"

    hess = get("xtb_hessian").execute(_wd("hess"), geometry=geo, cores=2, charge=0)
    assert hess.status.value == "success"
    assert hess.data["n_imaginary"] == 0 and "free_energy_Eh" in hess.data
    assert "hess" in hess.artifacts            # hessian artifact for downstream wiring


@pytest.mark.skipif(not _HAS_XTB, reason="xtb binary not installed")
def test_xtb_optimize_then_hessian_pipeline():
    # opt → hessian flow through a real pipeline run
    from delfin.tools import Pipeline

    d = Path(tempfile.mkdtemp())
    geo = d / "h2o.xyz"
    geo.write_text(_H2O)
    p = Pipeline("xtb_chain")
    p.add("xtb_optimize", charge=0, label="opt")
    p.add("xtb_hessian", charge=0, label="hess")
    res = p.run(geometry=geo, work_dir=d, cores=2)
    assert res.ok
