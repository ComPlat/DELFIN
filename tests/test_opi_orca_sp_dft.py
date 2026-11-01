"""Tests for the OPI-backed single-point DFT building block.

The building block (``delfin.tools.adapters.orca_sp_opi``) is DELFIN's first
ORCA step built on OPI (orca-pi 2.0.0): every keyword and block option is
checked against DELFIN's ORCA knowledge store before the input is written
(``delfin.orca_knowledge``), the input is assembled with OPI's input classes
rather than string templates, the run is launched through DELFIN's existing
``delfin.orca.run_orca`` (never a new runner), and the result is read back
with OPI while honouring the rule that no value is reported until the
convergence flags say it is meaningful.

These fixtures were measured once against ORCA 6.1.1 runs kept under
``tests/fixtures/orca611/`` (H2O opt+freq, O2 triplet, SCF-not-converged).
A fake ORCA script copies those runs into the working directory, so the tests
never invoke the real ORCA binary.
"""

from __future__ import annotations

import pytest

# The main py3.11 CI jobs run ALL tests without OPI installed; skip there.
# The OPI job (DELFIN_OPI_REQUIRED=1, tests/conftest.py) turns that skip
# into a failure, so the OPI-required job still gates on these tests.
pytest.importorskip("opi")

from delfin.tools.adapters.orca_sp_opi import (  # noqa: E402
    validate_method,
)


@pytest.fixture(scope="module")
def valid_method() -> dict:
    """A method the knowledge store accepts (B3LYP/def2-SVP/D4/RIJCOSX+def2/J)."""
    return {
        "functional": "B3LYP",
        "basis": "def2-SVP",
        "dispersion": "D4",
        "ri": "RIJCOSX",
        "aux_basis": "def2/J",
    }


class TestValidateMethod:
    """Every keyword is checked against delfin.orca_knowledge before use.
    An unknown keyword or block option is an error that stops the step; the
    knowledge-store cards' findings come back as warnings."""

    def test_known_method_produces_no_errors(self, valid_method):
        errors, warnings = validate_method(**valid_method)
        assert errors == []
        assert warnings == []

    def test_unknown_functional_is_an_error(self):
        errors, warnings = validate_method(functional="NOTAFUNCTIONAL", basis="def2-SVP")
        assert any("NOTAFUNCTIONAL" in e for e in errors)
        assert warnings == []

    def test_unknown_basis_is_an_error(self):
        errors, _ = validate_method(functional="B3LYP", basis="not-a-basis")
        assert any("not-a-basis" in e for e in errors)

    def test_unknown_dispersion_is_an_error(self):
        errors, _ = validate_method(functional="B3LYP", basis="def2-SVP", dispersion="D999")
        assert any("D999" in e for e in errors)

    def test_unknown_ri_is_an_error(self):
        errors, _ = validate_method(functional="B3LYP", basis="def2-SVP", ri="RI-NOTREAL")
        assert any("RI" in e for e in errors)

    def test_unknown_scf_option_is_an_error(self):
        errors, _ = validate_method(
            functional="B3LYP",
            basis="def2-SVP",
            scf_options={"maxiter_wrong": 10},
        )
        assert any("scf" in e for e in errors)

    def test_knowledge_card_warning_is_reported_not_silenced(self):
        """A finding from a verified knowledge card comes back as a warning,
        never silently dropped."""
        card = {
            "id": "test.card-warns",
            "topic": "scf",
            "kind": "recommendation",
            "statement": "A TIGHTSCF calc must set %scf maxiter (test card, not shipped).",
            "orca": ">=6.0",
            "status": "verified",
            "severity": "warning",
            "applies": {"keyword": "TIGHTSCF"},
            "check": {"requires_keyword": ["MAXITER"]},
            "sources": [{"type": "manual", "version": "6.1", "page": "x"}],
        }
        errors, warnings = validate_method(
            functional="B3LYP",
            basis="def2-SVP",
            extra_keywords=["TIGHTSCF"],
            cards=[card],
        )
        assert errors == []
        assert any(w.get("id") == "test.card-warns" for w in warnings)


# ---------------------------------------------------------------------------
# Run + read path: a fake ORCA stands in for the binary, copying a fixture
# run's outputs into the working directory.  The fixtures were measured once
# against ORCA 6.1.1 (H2O opt+freq, O2 triplet, SCF-not-converged).
# ---------------------------------------------------------------------------

from pathlib import Path  # noqa: E402

import pytest  # noqa: E402

import delfin.orca as orca_mod  # noqa: E402
from delfin.tools.adapters.orca_sp_opi import run_sp_dft  # noqa: E402

# Real ORCA fixtures live in the private store outside the repo
# (LICENSE_GUARD.md rule family 2); skip without DELFIN_PRIVATE_FIXTURES.
from test_lc1_fixtures_helper import lc1_fixture_dir  # noqa: E402

FIXTURE_DIR = Path("tests/fixtures/orca611")


@pytest.fixture()
def fake_orca(monkeypatch, tmp_path):
    """Point DELFIN's launcher at tests/scripts/fake_orca.sh so a run copies a
    fixture's outputs into the working directory instead of launching ORCA."""
    script = Path("tests/scripts/fake_orca.sh").resolve()
    assert script.is_file()

    def patch(fixture_name: str) -> Path:
        fix_dir = lc1_fixture_dir(f"tests/fixtures/orca611/{fixture_name}").resolve()
        assert fix_dir.is_dir(), f"fixture {fixture_name!r} missing"
        monkeypatch.setattr(orca_mod, "find_orca_executable", lambda: str(script))
        import delfin.qm_runtime as qmr

        qmr.clear_resolver_cache()
        monkeypatch.setenv("FAKE_ORCA_FIXTURE", str(fix_dir))
        return fix_dir

    return patch


class TestRunAndRead:
    """run_sp_dft drives DELFIN's own launcher with a fake ORCA and reads the
    result back with OPI, gating every number on the convergence flags."""

    def test_converged_o2_reports_meaningful_values(self, fake_orca, tmp_path):
        fake_orca("o2_triplet")
        res = run_sp_dft(
            elements=["O", "O"],
            coordinates=[(0.0, 0.0, 0.0), (0.0, 0.0, 1.207)],
            charge=0,
            mult=3,
            functional="B3LYP",
            basis="def2-SVP",
            ri="RIJCOSX",
            aux_basis="def2/J",
            work_dir=tmp_path / "o2",
        )
        assert res.ok is True
        assert res.meaningful is True
        assert res.terminated_normally is True
        assert res.scf_converged is True
        assert res.energy_Eh is not None and res.energy_Eh < 0
        assert res.hl_gap_eV is not None and res.hl_gap_eV > 0
        assert res.s2 is not None and res.s2[0] > 2.0  # open-shell <S^2> ~ 2.007
        assert res.failure_lines == ()
        # the input file really was written by OPI into the working dir
        assert (res.work_dir / "calc.inp").is_file()

    def test_converged_h2o_reports_meaningful_values(self, fake_orca, tmp_path):
        fake_orca("h2o_optfreq")
        res = run_sp_dft(
            elements=["O", "H", "H"],
            coordinates=[(0.0, 0.0, 0.1173), (0.0, 0.7572, -0.4692),
                         (0.0, -0.7572, -0.4692)],
            charge=0,
            mult=1,
            functional="B3LYP",
            basis="def2-SVP",
            work_dir=tmp_path / "h2o",
        )
        assert res.ok is True and res.meaningful is True
        assert res.energy_Eh is not None and res.energy_Eh < 0
        assert res.hl_gap_eV is not None
        assert res.s2 is None or res.s2[0] == 0.0  # closed shell

    def test_failed_run_never_reports_numbers(self, fake_orca, tmp_path):
        """A run OPI marks failed must not produce a HOMO-LUMO gap from
        unconverged orbitals: every number stays None, the failure is
        explained by .out lines."""
        fake_orca("scf_not_converged")
        res = run_sp_dft(
            elements=["O", "O"],
            coordinates=[(0.0, 0.0, 0.0), (0.0, 0.0, 1.207)],
            charge=0,
            mult=3,
            functional="B3LYP",
            basis="def2-SVP",
            work_dir=tmp_path / "fail",
        )
        assert res.ok is False
        assert res.meaningful is False
        assert res.terminated_normally is False
        assert res.scf_converged is False
        assert res.energy_Eh is None
        assert res.hl_gap_eV is None
        assert res.s2 is None
        assert any("SCF NOT CONVERGED" in ln for ln in res.failure_lines)

    def test_unknown_keyword_never_launches(self, fake_orca, tmp_path):
        """Validation errors stop the step before any ORCA launch."""
        fake_orca("o2_triplet")
        res = run_sp_dft(
            elements=["O", "O"],
            coordinates=[(0.0, 0.0, 0.0), (0.0, 0.0, 1.207)],
            charge=0,
            mult=3,
            functional="NOTAFUNCTIONAL",
            basis="def2-SVP",
            work_dir=tmp_path / "nolaunch",
        )
        assert res.ok is False
        assert "NOTAFUNCTIONAL" in (res.error or "")
        assert not (res.work_dir / "calc.out").exists()


