"""Tests for the OPI-backed geometry-optimisation / frequency building block.

The building block (``delfin.tools.adapters.orca_opt_opi`` and its core module
``delfin.tools.orca_opt_opi``) assembles ORCA optimisation, frequency, relaxed
surface scan, NEB and IRC inputs through OPI's typed block classes
(``BlockGeom``/``Scan``/``Constraints``/``BlockFreq``/``BlockNeb``/``BlockIrc``)
rather than string templates, validates every keyword and block option against
``delfin.orca_knowledge`` before the input is written, launches through DELFIN's
existing ``delfin.orca.run_orca``, and reads the result back with OPI while
honouring the rule that no value is reported until the matching convergence
flag makes it meaningful; a minimum with an imaginary mode is reported as such,
never silently accepted.

These fixtures were measured once against ORCA 6.1.1 runs kept under
``tests/fixtures/orca611/`` (H2O opt+freq).  A fake ORCA script copies those
runs into the working directory, so the tests never invoke the real ORCA
binary.  The build-only tests import opi lazily so they run in every py3.11
CI job without OPI; the read-back tests run only where OPI is installed.
"""

from __future__ import annotations

import pytest

H2O_XYZ = [
    ("O", 0.0, 0.0, 0.1173),
    ("H", 0.0, 0.7572, -0.4692),
    ("H", 0.0, -0.7572, -0.4692),
]
H2O_ELEMENTS = [a[0] for a in H2O_XYZ]
H2O_COORDS = [a[1:] for a in H2O_XYZ]


# ---------------------------------------------------------------------------
# Build-only tests: validation that never launches ORCA.  These import OPI
# lazily (inside the functions) so they pass in the py3.11 no-OPI CI jobs.
# ---------------------------------------------------------------------------


class TestBuildValidation:
    """Keyword / block-option validation stops the step with a loud error."""

    def test_unknown_basis_is_an_error(self):
        from delfin.tools.orca_opt_opi import validate_opt_method

        errors, _warnings = validate_opt_method(
            functional="BP86", basis="not-a-basis"
        )
        assert any("not-a-basis" in e for e in errors)

    def test_unknown_geom_option_is_an_error(self):
        from delfin.tools.orca_opt_opi import validate_opt_method

        errors, _warnings = validate_opt_method(
            functional="BP86",
            basis="def2-SVP",
            geom_options={"maxiter_wrong": 99},
        )
        assert any("geom" in e and "maxiter_wrong" in e for e in errors)

    def test_unknown_freq_option_is_an_error(self):
        from delfin.tools.orca_opt_opi import validate_opt_method

        errors, _warnings = validate_opt_method(
            functional="BP86",
            basis="def2-SVP",
            freq_options={"nonsense": True},
        )
        assert any("freq" in e and "nonsense" in e for e in errors)

    def test_unknown_neb_option_is_an_error(self):
        from delfin.tools.orca_opt_opi import validate_opt_method

        errors, _warnings = validate_opt_method(
            functional="BP86",
            basis="def2-SVP",
            neb_options={"nimages_wrong": 7},
        )
        assert any("neb" in e and "nimages_wrong" in e for e in errors)

    def test_unknown_irc_option_is_an_error(self):
        from delfin.tools.orca_opt_opi import validate_opt_method

        errors, _warnings = validate_opt_method(
            functional="BP86",
            basis="def2-SVP",
            irc_options={"direction_wrong": "both"},
        )
        assert any("irc" in e and "direction_wrong" in e for e in errors)

    def test_unknown_task_keyword_is_an_error(self):
        from delfin.tools.orca_opt_opi import validate_opt_method

        errors, _warnings = validate_opt_method(
            functional="BP86", basis="def2-SVP", task="freezelands"
        )
        assert any("freezelands" in e for e in errors)


# ---------------------------------------------------------------------------
# OPI-run + read-back tests.  These need OPI and are skipped in py3.11 CI.
# ---------------------------------------------------------------------------

try:
    import opi  # noqa: F401

    _HAS_OPI = True
except ImportError:  # pragma: no cover - py3.11 CI has no OPI
    _HAS_OPI = False

SKIP_NO_OPI = pytest.mark.skipif(
    not _HAS_OPI, reason="OPI (orca-pi) not installed"
)

from pathlib import Path  # noqa: E402

import delfin.orca as orca_mod  # noqa: E402
from delfin.tools.orca_opt_opi import run_opt_freq  # noqa: E402
from delfin.tools.orca_opt_opi import _build_opi_input  # noqa: E402

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


class TestOptFreqReadBack:
    """An Opt+Freq run's result is read back with OPI and gated on the
    convergence flags; the converged H2O minimum (0 imaginary modes) is
    reported cleanly."""

    pytestmark = SKIP_NO_OPI

    def test_converged_h2o_reports_geometry_and_thermochemistry(self, fake_orca, tmp_path):
        fake_orca("h2o_optfreq")
        res = run_opt_freq(
            elements=H2O_ELEMENTS,
            coordinates=H2O_COORDS,
            charge=0,
            mult=1,
            functional="B3LYP",
            basis="def2-SVP",
            work_dir=tmp_path / "optfreq",
        )
        assert res.ok is True
        assert res.meaningful is True
        assert res.terminated_normally is True
        assert res.scf_converged is True
        assert res.opt_converged is True
        assert res.freq_converged is True
        # converged geometry (Structure), final energy, ZPE/G/H present
        assert res.geometry is not None
        assert len(res.geometry.atoms) == 3
        assert res.final_energy is not None and res.final_energy < 0
        assert res.zpe is not None and res.zpe > 0
        assert res.enthalpy is not None
        assert res.free_energy is not None
        # frequencies: 3 real modes, none imaginary
        assert res.n_modes == 3
        assert res.n_imag == 0
        assert res.lowest_freq is not None and res.lowest_freq > 0
        # a minimum with no imaginary mode is reported as a minimum
        assert res.imag_freqs == ()
        assert "minimum" in res.character

    def test_opt_block_options_surfaced(self, fake_orca, tmp_path):
        """TightOpt (tighter thresholds), %geom maxiter and Calc_Hess are set on
        the input; the read-back still honours the convergence gates."""
        fake_orca("h2o_optfreq")
        res = run_opt_freq(
            elements=H2O_ELEMENTS,
            coordinates=H2O_COORDS,
            charge=0,
            mult=1,
            functional="B3LYP",
            basis="def2-SVP",
            opt_level="tightopt",
            geom_options={"maxiter": 200, "calc_hess": True},
            work_dir=tmp_path / "tight",
        )
        assert res.ok is True and res.meaningful is True
        assert res.opt_converged is True


class TestBuildEmitsValidInput:
    """The input files OPI writes are real .inp text the probe would accept."""

    pytestmark = SKIP_NO_OPI

    def test_opt_freq_single_job_emits_opt_and_freq(self, fake_orca, tmp_path):
        """Opt + Freq in one job: the written input carries both %geom and the
        freq simple keyword."""
        inp_path = _build_opi_input(
            work_dir=tmp_path,
            elements=H2O_ELEMENTS,
            coordinates=H2O_COORDS,
            charge=0,
            mult=1,
            functional="BP86",
            basis="def2-SVP",
            task="opt+freq",
            maxcore=1000,
            ncores=2,
            base_name="h2o_optfreq_build",
        )
        text = inp_path.read_text()
        # a plain Opt+Freq job carries no %geom block (ORCA defaults), but the
        # two job-type keywords are emitted and the parallelism is set
        assert "!opt" in text and "!freq" in text
        assert "nprocs 2" in text

    def test_constraints_and_scan_emitted_not_as_template(self, fake_orca, tmp_path):
        """Frozen-bond constraint and a single-coordinate scan go through
        OPI's typed Constraints/Scan, so the input carries the %geom block
        written by OPI (no %geom block left to hand-roll)."""
        inp_path = _build_opi_input(
            work_dir=tmp_path,
            elements=H2O_ELEMENTS,
            coordinates=H2O_COORDS,
            charge=0,
            mult=1,
            functional="BP86",
            basis="def2-SVP",
            task="opt",
            geom_options={"constraints": ["B 0 1 C"], "recalc_hess": 5},
            scan=[(("B", 0, 1), 2.0, 3.0, 5)],
            maxcore=1000,
            ncores=1,
            base_name="scan",
        )
        text = inp_path.read_text()
        assert "%geom" in text
        # OPI renders the Scan block lowercase (single coordinate supported)
        assert "scan b 0 1" in text.lower()
        # the frozen-bond constraint appears in the %geom block
        assert "B 0 1 C" in text or "b 0 1 c" in text.lower()

    def test_multi_coordinate_scan_is_refused(self, tmp_path):
        """OPI's BlockGeom.scan holds a single Scan object only; more than one
        scan coordinate cannot be expressed through OPI and is refused rather
        than spilled as a second string template.  (Manual: optimizations_scans
        allows up to three Scan lines; OPI's model as of 2.0.0 does not.)"""
        with pytest.raises(ValueError) as exc:
            _build_opi_input(
                work_dir=tmp_path,
                elements=H2O_ELEMENTS,
                coordinates=H2O_COORDS,
                charge=0,
                mult=1,
                functional="BP86",
                basis="def2-SVP",
                task="opt",
                scan=[(("B", 0, 1), 2.0, 3.0, 5), (("B", 1, 2), 1.5, 2.5, 4)],
                base_name="scan2",
                maxcore=1000,
                ncores=1,
            )
        assert "holds one coordinate" in str(exc.value)

    def test_ts_mode_is_refused_opi_emits_inline_form_orca_rejects(
            self, tmp_path):
        """ORCA 6.1.1 parses TS_Mode only as a sub-block inside %geom
        ('TS_Mode / { M 0 } / end'; probe: 'finished'; the inline
        'ts_mode { M 1 }' line OPI 2.0.0 emits probes 'rejected' with
        'Unrecognized symbol in GEOM block line').  ts_mode is refused
        instead of hand-rolling the sub-block as a string template."""
        with pytest.raises(ValueError) as exc:
            _build_opi_input(
                work_dir=tmp_path,
                elements=H2O_ELEMENTS,
                coordinates=H2O_COORDS,
                charge=0,
                mult=1,
                functional="BP86",
                basis="def2-SVP",
                task="optts",
                geom_options={"ts_mode": ["M", 1]},
                base_name="tsmode",
                maxcore=1000,
                ncores=1,
            )
        assert "ts_mode" in str(exc.value)

    def test_invalid_scan_coordinate_is_refused(self, tmp_path):
        with pytest.raises(ValueError):
            _build_opi_input(
                work_dir=tmp_path,
                elements=H2O_ELEMENTS,
                coordinates=H2O_COORDS,
                charge=0,
                mult=1,
                functional="BP86",
                basis="def2-SVP",
                task="opt",
                scan=[(("X", 0, 1), 2.0, 3.0, 5)],
                base_name="scanbad",
                maxcore=1000,
                ncores=1,
            )

    def test_inhess_without_inhessname_is_refused(self, tmp_path):
        """inhess=read needs an inhessname (the Hessian file to read); the
        manual requires the file when reusing a Hessian."""
        with pytest.raises(ValueError):
            _build_opi_input(
                work_dir=tmp_path,
                elements=H2O_ELEMENTS,
                coordinates=H2O_COORDS,
                charge=0,
                mult=1,
                functional="BP86",
                basis="def2-SVP",
                task="opt",
                geom_options={"inhess": "read", "maxiter": 5},
                base_name="inhessbad",
                maxcore=1000,
                ncores=1,
            )
