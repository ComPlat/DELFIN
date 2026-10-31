"""Tests for the OPI-backed DLPNO-CCSD(T) single-point building block.

``delfin.tools.adapters.orca_dlpno_opi`` is DELFIN's local-correlation
single-point step on OPI (orca-pi 2.0.0).  Manual ground truth (ORCA 6.1
manual, delfin-knowledge/orca/manual-6.1-text/contents/modelchemistries/
mdci.txt):

* Section 3.10.10.1 "Basic Usage": the example
  ``! cc-pVTZ cc-pVTZ/C DLPNO-CCSD(T) TightSCF`` -- the auxiliary basis is
  the ``/C`` fitting basis of the same family as the orbital basis.
* "Sensible defaults have been determined from extensive benchmark
  calculations and are accessible via LoosePNO, NormalPNO and TightPNO
  keywords in the simple input line.  These keywords represent the
  recommended way to control the accuracy of DLPNO calculations"; the
  TightPNO example is annotated "# Tight settings for increased accuracy,
  e.g. when investigating weak interactions or conformational equilibria".
* CBS: "Using the keyword ``! Extrapolate(X/Y,basis)`` ... ORCA will
  calculate the SCF and optionally the MP2 or MDCI energies with two basis
  sets and separately extrapolate"; the printed exponents are
  "Alpha(2/3): 4.420, Beta(2/3): 2.460" (cc family) and the def2-family
  table values alpha_23 = 10.39, beta_23 = 2.40.  "The extrapolation
  method should only be used with very tight SCF convergence criteria."
  "Auxiliary basis sets are automatically chosen and can not be changed."
* DLPNO-CCSD(T1): "the improved iterative perturbative triples correction"
  (vs the semi-canonical DLPNO-CCSD(T)).
* T1 diagnostic: "A rule of thumb says, that for a value of the diagnostic
  of larger than 0.02 the results are not to be trusted."

The measured-run tests use a fake ORCA that copies the ORCA 6.1.1 fixture
``tests/fixtures/orca611/butane_dlpno_ccsdt_pvtz`` (anti-butane
DLPNO-CCSD(T)/cc-pVTZ single point) into the working directory.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("opi")

from delfin.tools.adapters.orca_dlpno_opi import (  # noqa: E402
    build_dlpno_input,
    run_dlpno,
    validate_dlpno,
)

import delfin.orca as orca_mod  # noqa: E402

# Real ORCA fixtures live in the private store outside the repo
# (LICENSE_GUARD.md rule family 2); skip without DELFIN_PRIVATE_FIXTURES.
from test_lc1_fixtures_helper import lc1_fixture_dir  # noqa: E402

#: anti-butane from ideal geometry (C-C 1.53 A, C-H 1.09 A, anti dihedral).
BUTANE_ELEMENTS = ("C", "C", "C", "C", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H")
BUTANE_COORDS = (
    (-1.930, 0.000, 0.770),
    (-0.770, 0.000, 0.385),
    (0.770, 0.000, -0.385),
    (1.930, 0.000, -0.770),
    (-1.930, 0.886, 1.385),
    (-1.930, -0.886, 1.385),
    (-2.780, 0.000, 0.085),
    (-0.770, 0.886, -0.230),
    (-0.770, -0.886, -0.230),
    (0.770, 0.886, 0.230),
    (0.770, -0.886, 0.230),
    (1.930, 0.886, -1.385),
    (1.930, -0.886, -1.385),
    (2.780, 0.000, -0.085),
)


@pytest.fixture()
def fake_orca(monkeypatch):
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


class TestValidateDlpno:
    """Keyword and option validation against the knowledge store."""

    def test_default_method_is_dlno_ccsdt_tightpno(self):
        """The module's suggested default: DLPNO-CCSD(T)/TightPNO -- the
        manual's TightPNO example names conformational equilibria."""
        errors, _ = validate_dlpno(basis="cc-pVTZ")
        assert errors == []
        res = suggested_default()
        assert res["method"] == "DLPNO-CCSD(T)"
        assert res["pno"] == "TightPNO"
        assert "conformational" in res["reason"].lower()

    def test_unknown_method_is_an_error(self):
        errors, _ = validate_dlpno(method="CCSD-NOPE", basis="cc-pVTZ")
        assert any("CCSD-NOPE" in e for e in errors)

    def test_unknown_pno_is_an_error(self):
        errors, _ = validate_dlpno(pno="TIGHTPNO-WRONG", basis="cc-pVTZ")
        assert any("TIGHTPNO-WRONG" in e for e in errors)

    def test_aux_basis_must_match_orbital_family(self):
        """The /C auxiliary basis belongs to the orbital basis's family:
        cc-pVnZ -> cc-pVnZ/C, def2-XZVPP -> def2-XZVPP/C.  A mismatched
        pairing is refused before any input is written."""
        errors, _ = validate_dlpno(basis="cc-pVTZ", aux_basis="def2-TZVPP/C")
        assert any("def2-TZVPP/C" in e for e in errors)

    def test_matching_aux_basis_accepted(self):
        errors, _ = validate_dlpno(basis="cc-pVTZ", aux_basis="cc-pVTZ/C")
        assert errors == []

    def test_rijcosx_or_rijk_accepted(self):
        errors, _ = validate_dlpno(basis="cc-pVTZ", ri="RIJCOSX")
        assert errors == []
        errors, _ = validate_dlpno(basis="cc-pVTZ", ri="RIJK")
        assert errors == []

    def test_cbs_extrapolation_pairs(self):
        """Two-point CBS: (cc-pVTZ, cc-pVQZ) and (def2-TZVPP, def2-QZVPP)
        map onto the manual's Extrapolate(X/Y,basis) keyword."""
        errors, _ = validate_dlpno(cbs=("cc-pVTZ", "cc-pVQZ"))
        assert errors == []
        errors, _ = validate_dlpno(cbs=("def2-TZVPP", "def2-QZVPP"))
        assert errors == []
        errors, _ = validate_dlpno(cbs=("cc-pVTZ", "def2-QZVPP"))
        assert any("family" in e.lower() for e in errors)

    def test_cbs_pair_must_be_successive_cardinals(self):
        errors, _ = validate_dlpno(cbs=("cc-pVDZ", "cc-pVQZ"))
        assert errors != []


def suggested_default():
    from delfin.tools.adapters.orca_dlpno_opi import suggested_default as sd

    return sd()


class TestBuildDlpnoInput:
    def test_written_with_dlpno_keywords(self, tmp_path):
        inp = build_dlpno_input(
            tmp_path,
            elements=BUTANE_ELEMENTS,
            coordinates=BUTANE_COORDS,
            basis="cc-pVTZ",
            base_name="dlpno",
        )
        assert inp.is_file()
        text = inp.read_text()
        assert "dlpno-ccsd(t)" in text.lower()
        assert "tightpno" in text.lower()
        assert "cc-pvtz/c" in text.lower()

    def test_cbs_keyword_written(self, tmp_path):
        inp = build_dlpno_input(
            tmp_path,
            elements=BUTANE_ELEMENTS,
            coordinates=BUTANE_COORDS,
            cbs=("cc-pVTZ", "cc-pVQZ"),
            base_name="dlpno",
        )
        text = inp.read_text()
        assert "extrapolate" in text.lower()

    def test_unknown_keyword_aborts_before_any_file(self, tmp_path):
        with pytest.raises(Exception):
            build_dlpno_input(
                tmp_path,
                elements=BUTANE_ELEMENTS,
                coordinates=BUTANE_COORDS,
                basis="cc-pVTZ",
                extra_keywords=["NOTAWORD"],
                base_name="dlpno",
            )
        assert not (tmp_path / "dlpno.inp").exists()


class TestRunDlpno:
    def test_converged_run_reports_ccsd_t_energy_and_t1(self, fake_orca, tmp_path):
        fake_orca("butane_dlpno_ccsdt_pvtz")
        res = run_dlpno(
            elements=BUTANE_ELEMENTS,
            coordinates=BUTANE_COORDS,
            charge=0,
            mult=1,
            basis="cc-pVTZ",
            work_dir=tmp_path / "dlpno",
        )
        assert res.ok is True
        assert res.meaningful is True
        assert res.terminated_normally is True
        assert res.energy_Eh is not None and res.energy_Eh < 0
        assert res.ccsd_t_energy_Eh is not None and res.ccsd_t_energy_Eh < 0
        assert res.t1_diagnostic is not None
        assert isinstance(res.t1_diagnostic, float)
        assert res.triples_correction == "(T)"
        assert res.aux_basis_used == "cc-pVTZ/C"

    def test_failed_run_reports_no_numbers(self, fake_orca, tmp_path):
        fake_orca("scf_not_converged")
        res = run_dlpno(
            elements=BUTANE_ELEMENTS,
            coordinates=BUTANE_COORDS,
            charge=0,
            mult=1,
            basis="cc-pVTZ",
            work_dir=tmp_path / "fail",
        )
        assert res.ok is False
        assert res.meaningful is False
        assert res.energy_Eh is None
        assert res.ccsd_t_energy_Eh is None
        assert res.t1_diagnostic is None
        assert res.failure_lines != ()

    def test_unknown_keyword_never_launches(self, fake_orca, tmp_path):
        fake_orca("butane_dlpno_ccsdt_pvtz")
        res = run_dlpno(
            elements=BUTANE_ELEMENTS,
            coordinates=BUTANE_COORDS,
            charge=0,
            mult=1,
            basis="cc-pVTZ",
            extra_keywords=["NOTAWORD"],
            work_dir=tmp_path / "nolaunch",
        )
        assert res.ok is False
        assert "NOTAWORD" in (res.error or "")
        assert not (res.work_dir / "dlpno.out").exists()
