"""Fixture test for tests/fixtures/orca611/ru_nh3_6_cpcm/.

The fixture is a real ORCA 6.1.1 single point of [Ru(NH3)6]2+ (def2-TZVP with
its def2-ECP, CPCM(water)) built by delfin/tools/orca_basis_solv_opi.py and
run through delfin.orca.run_orca in scratch/.  It carries the ECP lines the
task asked for.  Two things are locked here:
  1. fixture hygiene: no absolute /home path, no user name, no host name
     (the orca611 fixture convention: Working dir. -> the fixture's own
     relative path, Host name -> fixture-host),
  2. the run it records really terminated normally with the def2-ECP for Ru
     replacing 28 core electrons, and the adapter can read it back through
     the fake ORCA that a real run copies fixture outputs from.
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("opi")

FIXTURE = Path("tests/fixtures/orca611/ru_nh3_6_cpcm")


def test_fixture_exists_and_is_scrubbed():
    for name in ("ru_nh3_6_cpcm.inp", "ru_nh3_6_cpcm.out",
                 "ru_nh3_6_cpcm.property.txt"):
        path = FIXTURE / name
        assert path.is_file(), f"fixture file missing: {name}"
        txt = path.read_text(errors="replace")
        assert "/home/" not in txt, f"{name} carries an absolute /home path"
        assert "qmchem" not in txt, f"{name} carries user or host identity"


def test_fixture_records_terminated_run_with_ru_ecp():
    txt = (FIXTURE / "ru_nh3_6_cpcm.out").read_text(errors="replace")
    assert "ORCA TERMINATED NORMALLY" in txt
    # def2-ECP for Ru: replaces 28 core electrons, lmax 3 (Rb-Rn rule; Ru=Z 44)
    assert "Type Ru ECP Def2-ECP (replacing 28 core electrons" in txt
    assert "* core charge reduced due to ECP" in txt
    # the solvation block the builder produced for the run
    assert "CPCM" in txt


def test_adapter_reads_the_fixture_back_through_fake_orca(monkeypatch, tmp_path):
    import delfin.orca as orca_mod
    from delfin.tools.adapters.orca_basis_solv_opi import (
        OrcaBasisSolvOpiAdapter,
    )

    script = Path("tests/scripts/fake_orca.sh").resolve()
    monkeypatch.setattr(orca_mod, "find_orca_executable", lambda: str(script))
    monkeypatch.setenv(
        "FAKE_ORCA_FIXTURE", str(FIXTURE.resolve())
    )
    adapter = OrcaBasisSolvOpiAdapter()
    res = adapter.execute(
        tmp_path,
        elements=["Ru"] + ["N"] * 6 + ["H"] * 18,
        coordinates=[
            (0.0, 0.0, 0.0),
            (2.1, 0.0, 0.0), (-2.1, 0.0, 0.0),
            (0.0, 2.1, 0.0), (0.0, -2.1, 0.0),
            (0.0, 0.0, 2.1), (0.0, 0.0, -2.1),
            *[(2.43667, 0.95224, 0.0), (2.43667, -0.95224, 0.0),
              (2.43667, 0.0, 0.95224), (-2.43667, 0.95224, 0.0),
              (-2.43667, -0.95224, 0.0), (-2.43667, 0.0, 0.95224),
              (0.95224, 2.43667, 0.0), (-0.95224, 2.43667, 0.0),
              (0.0, 2.43667, 0.95224), (0.95224, -2.43667, 0.0),
              (-0.95224, -2.43667, 0.0), (0.0, -2.43667, 0.95224),
              (0.95224, 0.0, 2.43667), (0.0, 0.95224, 2.43667),
              (-0.95224, 0.0, 2.43667), (0.95224, 0.0, -2.43667),
              (0.0, 0.95224, -2.43667), (-0.95224, 0.0, -2.43667)],
        ],
        charge=2,
        mult=1,
        functional="PBE",
        basis="def2-TZVP",
        solvent_model="CPCM",
        solvent="water",
        maxcore=2000,
        cores=2,
    )
    assert res.status.value == "success", res.error
    assert res.data["meaningful"] is True
    assert res.data["energy_Eh"] is not None and res.data["energy_Eh"] < 0
