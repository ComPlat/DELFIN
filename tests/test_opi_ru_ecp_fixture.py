"""Fixture tests for tests/fixtures/orca611/ru_nh3_6_cpcm/.

The real ORCA 6.1.1 single point of [Ru(NH3)6]2+ (PBE/def2-TZVP with its
auto-loaded def2-ECP, CPCM(water)) lives in the private fixture store
(LICENSE_GUARD.md rule family 2: real licensed-program output is not
committed); tests read it through tests/test_lc1_fixtures_helper.py and SKIP
without DELFIN_PRIVATE_FIXTURES, so CI (no ORCA binary) is unaffected.  The
committed directory keeps the input DELFIN's builder wrote plus a
first-line-marked SYNTHETIC output so the reader path has CI coverage without
the ORCA program.

Locked here:
  1. committed-fixture hygiene and the synthetic marker as FIRST line,
  2. the synthetic .out quotes the measured ECP evidence (def2-ECP for Ru
     replacing 28 core electrons, core-charge-reduced note, normal
     termination),
  3. against the real output (private store) the same ECP evidence holds,
  4. the adapter drives the synthetic run back through the fake ORCA (no
     ORCA binary needed).
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("opi")

from test_lc1_fixtures_helper import private_fixtures_or_skip  # noqa: E402

FIXTURE = Path("tests/fixtures/orca611/ru_nh3_6_cpcm")
SYNTHETIC_MARKER = "# SYNTHETIC test fixture written for DELFIN -- not ORCA output"
ECP_GROUP_LINE = "Group 1, Type Ru ECP Def2-ECP (replacing 28 core electrons, lmax=3)"


def _text(p: Path) -> str:
    return p.read_text(encoding="utf-8", errors="replace")


def test_committed_fixture_hygiene_and_synthetic_marker():
    inp = _text(FIXTURE / "ru_nh3_6_cpcm.inp")
    assert "def2-TZVP" in inp and "CPCM" in inp
    for path in sorted(FIXTURE.iterdir()):
        body = _text(path)
        assert "/home/" not in body and "/Users/" not in body, path
    out = _text(FIXTURE / "ru_nh3_6_cpcm.out")
    assert out.startswith(SYNTHETIC_MARKER), "the marker must be the FIRST line"
    assert "Host name:     fixture-host" in out
    assert "Working dir.:  <worktree>" in out


def test_synthetic_output_quotes_the_measured_ecp_evidence():
    out = _text(FIXTURE / "ru_nh3_6_cpcm.out")
    assert ECP_GROUP_LINE in out
    assert "* core charge reduced due to ECP" in out
    assert "****ORCA TERMINATED NORMALLY****" in out


def test_real_output_in_private_store_shows_the_ecp_lines():
    # The committed copy is synthetic; the real run is only in the store.
    real = private_fixtures_or_skip(str(FIXTURE)) / "ru_nh3_6_cpcm.out"
    txt = _text(real)
    assert "ORCA TERMINATED NORMALLY" in txt
    assert ECP_GROUP_LINE in txt
    assert "* core charge reduced due to ECP" in txt


def test_adapter_drives_the_synthetic_run_through_fake_orca(monkeypatch, tmp_path):
    script = Path("tests/scripts/fake_orca.sh").resolve()
    if not script.is_file():
        pytest.skip("fake_orca.sh helper not present")
    monkeypatch.setenv("FAKE_ORCA_FIXTURE", str(FIXTURE.resolve()))
    import delfin.orca as orca_mod
    from delfin.tools.adapters.orca_basis_solv_opi import (
        OrcaBasisSolvOpiAdapter,
    )

    monkeypatch.setattr(orca_mod, "find_orca_executable", lambda: str(script))
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
