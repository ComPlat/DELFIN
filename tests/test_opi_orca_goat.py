"""Tests for the OPI-backed GOAT conformer-search building block.

``delfin.tools.adapters.orca_goat_opi`` is DELFIN's GOAT step on OPI
(orca-pi 2.0.0): the ``!GOAT`` simple keyword and the ``%goat`` block come
from OPI's typed input classes, every option is checked against DELFIN's
ORCA knowledge store before the input is written, the run is launched
through DELFIN's existing ``delfin.orca.run_orca``, and the final ensemble
is read back from the run directory (``<base>.finalensemble.xyz`` plus the
``# Final ensemble info #`` table in the ``.out``).

Manual ground truth (ORCA 6.1 manual, text extraction at
delfin-knowledge/orca/manual-6.1-text/contents/structurereactivity/goat.txt,
section 4.10 "GOAT: global geometry optimization and ensemble generator"):

* "The command to call the global optimizer is simply ``!GOAT``, like you
  would with ``!OPT``, and its options can be given under the ``%GOAT``
  block as usual."  "Here we will just use GFN2 (or ``!XTB``)."
* Table 4.9 ``%goat`` block keywords: MAXITER ("defines an arbitrary number
  of max GOAT geom opt iters per worker"), NWORKERS ("default auto ... any
  number multiple of 4 (number of temperatures)"), MAXEN ("the maximum
  relative energy of a conformer to be taken, in kcal/mol. 6 kcal/mol by
  default"), ENDIFF ("minimum energy difference needed to differentiate
  conformers, in kcal/mol"), CONFTEMP ("temperature used to compute the
  free energy, in Kelvin", default 298.15), RMSD (default 0.125 Angs),
  TEMPLIST ("a list of temperatures, defines the number of basic workers,
  in Kelvin").
* Ensemble output: "The final relative energies of the ensemble are
  printed afterwards, together with a ``BaseName.finalensemble.xyz``
  file" -- the ``# Final ensemble info #`` table has the columns
  "Conformer  Energy (kcal/mol)  Degen.  % total  % cumul." and the
  printed lines include "Lowest energy conformer    : <E> Eh".

These tests run against a fake ORCA that copies a measured ORCA 6.1.1
GOAT fixture run (``tests/fixtures/orca611/butane_goat_gfn2``) into the
working directory, so no test invokes the real ORCA binary.
"""

from __future__ import annotations

from pathlib import Path

import pytest

# The main py3.11 CI jobs run ALL tests without OPI installed; skip there.
# The OPI job (DELFIN_OPI_REQUIRED=1, tests/conftest.py) turns that skip
# into a failure, so the OPI-required job still gates on these tests.
pytest.importorskip("opi")

from delfin.tools.adapters.orca_goat_opi import (  # noqa: E402
    build_goat_input,
    run_goat,
    validate_goat,
)

import delfin.orca as orca_mod  # noqa: E402

# Real ORCA fixtures live in the private store outside the repo
# (LICENSE_GUARD.md rule family 2); skip without DELFIN_PRIVATE_FIXTURES.
from test_lc1_fixtures_helper import lc1_fixture_dir  # noqa: E402

#: n-butane built from ideal geometry (all-trans, C-C 1.53 A, C-H 1.09 A).
BUTANE_ELEMENTS = ("C", "C", "C", "C", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H")
BUTANE_COORDS = (
    (-1.930, 0.000, 0.770),   # C1
    (-0.770, 0.000, 0.385),   # C2
    (0.770, 0.000, -0.385),   # C3
    (1.930, 0.000, -0.770),   # C4
    (-1.930, 0.886, 1.385),   # H on C1
    (-1.930, -0.886, 1.385),
    (-2.780, 0.000, 0.085),
    (-0.770, 0.886, -0.230),  # H on C2
    (-0.770, -0.886, -0.230),
    (0.770, 0.886, 0.230),    # H on C3
    (0.770, -0.886, 0.230),
    (1.930, 0.886, -1.385),   # H on C4
    (1.930, -0.886, -1.385),
    (2.780, 0.000, -0.085),
)


@pytest.fixture()
def fake_orca(monkeypatch):
    """Point DELFIN's launcher at tests/scripts/fake_orca.sh so a run copies
    a fixture's outputs into the working directory instead of launching
    ORCA."""
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


class TestValidateGoat:
    """Every keyword and %goat block option is checked against the knowledge
    store; an unknown one is an error that stops the step, never dropped."""

    def test_default_gfn2_xtb_produces_no_errors(self):
        errors, warnings = validate_goat(method="XTB")
        assert errors == []
        assert warnings == []

    def test_unknown_method_is_an_error(self):
        errors, _ = validate_goat(method="NOT-A-METHOD")
        assert any("NOT-A-METHOD" in e for e in errors)

    def test_unknown_goat_option_is_an_error(self):
        errors, _ = validate_goat(goat_options={"maxiter_wrong": 10})
        assert any("goat" in e and "maxiter_wrong" in e for e in errors)

    def test_known_goat_options_produce_no_errors(self):
        errors, _ = validate_goat(
            goat_options={"maxiter": 128, "maxen": 6.0, "nworkers": 4, "conftemp": 298.15}
        )
        assert errors == []

    def test_unknown_simple_keyword_is_an_error(self):
        errors, _ = validate_goat(method="XTB", extra_keywords=["NOTAWORD"])
        assert any("NOTAWORD" in e for e in errors)


class TestBuildGoatInput:
    """The .inp is assembled with OPI's typed classes and written by OPI's
    Calculator -- the file's text is asserted only through the markers OPI
    itself writes, never by re-parsing a hand-written template."""

    def test_written_with_goat_keyword_and_block(self, tmp_path):
        inp = build_goat_input(
            tmp_path,
            elements=BUTANE_ELEMENTS,
            coordinates=BUTANE_COORDS,
            method="XTB",
            goat_options={"maxiter": 128, "nworkers": 4},
            base_name="goat",
        )
        assert inp.is_file()
        text = inp.read_text()
        assert "%goat" in text.lower()
        assert "maxiter" in text.lower()
        assert "nworkers" in text.lower()

    def test_unknown_option_aborts_before_any_file(self, tmp_path):
        with pytest.raises(Exception):
            build_goat_input(
                tmp_path,
                elements=BUTANE_ELEMENTS,
                coordinates=BUTANE_COORDS,
                method="XTB",
                goat_options={"definitelynotanoption": 1},
                base_name="goat",
            )
        assert not (tmp_path / "goat.inp").exists()


class TestRunGoat:
    """run_goat drives DELFIN's own launcher with a fake ORCA and reads the
    ensemble back, gated on terminated_normally."""

    def test_converged_run_reports_ensemble(self, fake_orca, tmp_path):
        fake_orca("butane_goat_gfn2")
        res = run_goat(
            elements=BUTANE_ELEMENTS,
            coordinates=BUTANE_COORDS,
            charge=0,
            mult=1,
            method="XTB",
            work_dir=tmp_path / "goat",
        )
        assert res.ok is True
        assert res.meaningful is True
        assert res.terminated_normally is True
        assert res.ensemble_file is not None
        assert len(res.energies_kcal_mol) >= 1
        # lowest conformer defines zero
        assert min(res.energies_kcal_mol) == pytest.approx(0.0, abs=1e-6)
        assert res.boltzmann_weights is not None
        assert len(res.boltzmann_weights) == len(res.energies_kcal_mol)
        assert res.boltzmann_weights[0] == pytest.approx(
            max(res.boltzmann_weights), rel=1e-6
        )
        assert res.structures is not None
        assert len(res.structures) == len(res.energies_kcal_mol)
        assert res.temperature_k == 298.15
        assert res.ensemble_file is not None

    def test_failed_run_reports_no_ensemble(self, fake_orca, tmp_path):
        fake_orca("scf_not_converged")
        res = run_goat(
            elements=BUTANE_ELEMENTS,
            coordinates=BUTANE_COORDS,
            charge=0,
            mult=1,
            method="XTB",
            work_dir=tmp_path / "fail",
        )
        assert res.ok is False
        assert res.meaningful is False
        assert res.ensemble_file is None
        assert res.energies_kcal_mol == ()
        assert res.failure_lines != ()

    def test_unknown_method_never_launches(self, fake_orca, tmp_path):
        fake_orca("butane_goat_gfn2")
        res = run_goat(
            elements=BUTANE_ELEMENTS,
            coordinates=BUTANE_COORDS,
            charge=0,
            mult=1,
            method="NOT-A-METHOD",
            work_dir=tmp_path / "nolaunch",
        )
        assert res.ok is False
        assert "NOT-A-METHOD" in (res.error or "")
        assert not (res.work_dir / "goat.out").exists()
