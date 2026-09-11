"""ORCA ending normally is not always a result, and the recovery treats it so.

"ORCA TERMINATED NORMALLY" also ends an optimisation that ran out of cycles,
a TD-DFT root that collapsed, and an ESD rate that is negative or cut off
before its correlation function decayed.  The recovery counted every such
run as done, so its strategies for them -- continue the optimisation, TDA,
a longer window -- never ran.  Now those runs are retried; when the retries
are used up the last result is kept, as before, and the report says what is
wrong with it.  An output that recalc kept is not judged again here: an
archived ISC of a TADF emitter took two days.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin import orca
from delfin.common import orca_input as oi
from delfin.orca_recovery import OrcaErrorDetector, OrcaErrorType as E, OrcaInputModifier, RecoveryStrategy

_NORMAL = "\n                             ****ORCA TERMINATED NORMALLY****\nTOTAL RUN TIME: 0 days 0 hours 0 minutes 12 seconds\n"
_GAVE_UP = """
                    The optimization did not converge but reached the maximum
                    number of optimization cycles.
                    Please check your results very carefully.
""" + _NORMAL
_CONVERGED = "\n                    ***********************HURRAY********************\n                    ***        THE OPTIMIZATION HAS CONVERGED     ***\n" + _NORMAL


def _esd(max_time_fs: float, rate: str, warn: bool = False) -> str:
    return ("                           ORCA EXCITED STATE DYNAMICS\n"
            "Homogeneous linewidth is:\t\t\t50.00 cm-1\n"
            "Number of points:\t\t\t\t131072\n"
            f"Maximum time:\t\t\t\t\t{max_time_fs:.2f} fs\n"
            f"The calculated ISC rate constant is\t\t{rate} s-1\n"
            + ("WARNING: negative rates are unphysical! It means something went wrong with the CorrFunc integration.\n" if warn else "")
            + _NORMAL)


def _banner(n: int) -> str:
    return f"\n                 $$$$$$$$$$$$$$$$  JOB NUMBER  {n} $$$$$$$$$$$$$$\n"


@pytest.mark.parametrize("text, expected", [
    (_GAVE_UP, (E.GEOMETRY_NOT_CONVERGED, 0)),
    (_CONVERGED, None),
    (_esd(290.27, "5.698207e+07"), (E.ESD_WINDOW_TRUNCATED, 0)),
    (_esd(2933.77, "-8.06e-03", warn=True), (E.ESD_RATE_UNPHYSICAL, 0)),
    (_esd(2933.77, "2.459783e+07"), None),
    (_banner(1) + _CONVERGED.replace(_NORMAL, "") + _banner(2) + _esd(290.27, "1.0e+06"), (E.ESD_WINDOW_TRUNCATED, 1)),
    (_banner(1) + _GAVE_UP.replace(_NORMAL, "") + _banner(2) + "single point\n" + _NORMAL, (E.GEOMETRY_NOT_CONVERGED, 0)),
])
def test_what_ended_normally_but_is_not_a_result_is_named(tmp_path, text, expected):
    out = tmp_path / "job.out"
    out.write_text(text)
    inp = tmp_path / "job.inp"
    inp.write_text("! PBE0 def2-SVP OPT ESD(ISC)\n")

    found = OrcaErrorDetector.unusable_result(out, inp)

    assert (found[:2] if found else None) == expected


def test_a_single_point_is_not_scanned(tmp_path):
    out = tmp_path / "job.out"
    out.write_text(_GAVE_UP)
    inp = tmp_path / "job.inp"
    inp.write_text("! PBE0 def2-SVP\n")

    assert OrcaErrorDetector.unusable_result(out, inp) is None


class _Orca:
    """run_orca stand-in: writes the next scripted output, returns ORCA's own verdict."""

    def __init__(self, outputs):
        self.outputs = list(outputs)
        self.inputs = []

    def __call__(self, inp, out, timeout=None, **kwargs):
        self.inputs.append(Path(inp).read_text())
        text = self.outputs.pop(0)
        Path(out).write_text(text)
        return "ORCA TERMINATED NORMALLY" in text


OPT_INPUT = "! PBE0 def2-SVP OPT\n%geom maxiter 5 end\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n"
ESD_INPUT = ("! PBE0 def2-SVP ESD(ISC)\n%base \"S1_T1_ISC_ms0\"\n%ESD\n  NPOINTS 131072\n  MAXTIME 12000\nEND\n"
             "* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n")


def _recover(tmp_path, monkeypatch, name, text, outputs):
    inp = tmp_path / name
    inp.write_text(text)
    fake = _Orca(outputs)
    monkeypatch.setattr(orca, "run_orca", fake)
    ok = orca.run_orca_with_intelligent_recovery(str(inp), str(tmp_path / (inp.stem + ".out")),
                                                 working_dir=tmp_path,
                                                 config={"enable_auto_recovery": "yes", "max_recovery_attempts": 2})
    return ok, fake


def test_an_optimisation_that_ran_out_of_cycles_is_continued(tmp_path, monkeypatch):
    ok, fake = _recover(tmp_path, monkeypatch, "initial.inp", OPT_INPUT, [_GAVE_UP, _CONVERGED])

    assert ok
    assert len(fake.inputs) == 2
    retry = oi.parse_job(fake.inputs[1])
    assert oi.block_value(retry, "geom", "maxiter") == "250"
    assert oi.setting_value(retry, "base") == '"initial"'


def test_a_window_cut_short_is_handed_back_to_orca(tmp_path, monkeypatch):
    ok, fake = _recover(tmp_path, monkeypatch, "S1_T1_ISC_ms0.inp", ESD_INPUT,
                        [_esd(290.27, "5.698207e+07"), _esd(2933.77, "2.459783e+07")])

    assert ok
    retry = fake.inputs[1]
    assert "%ESD" in retry.upper()
    assert "NPOINTS" not in retry and "MAXTIME" not in retry


def test_a_negative_rate_is_rerun_with_more_points(tmp_path, monkeypatch):
    ok, fake = _recover(tmp_path, monkeypatch, "S1_S0_IC.inp", ESD_INPUT.replace("ESD(ISC)", "ESD(IC)"),
                        [_esd(2933.77, "-8.06e-03", warn=True), _esd(2933.77, "3.79e-03")])

    assert ok
    retry = oi.parse_job(fake.inputs[1])
    assert oi.block_value(retry, "esd", "npoints") == str(4 * 131072)
    assert oi.block_value(retry, "esd", "maxtime") is None


def test_when_the_retries_are_used_up_the_last_result_is_kept(tmp_path, monkeypatch):
    ok, fake = _recover(tmp_path, monkeypatch, "initial.inp", OPT_INPUT, [_GAVE_UP, _GAVE_UP, _GAVE_UP])

    assert ok, "a run that ended normally keeps its result, as before"
    assert len(fake.inputs) == 3


def test_a_failed_run_still_fails(tmp_path, monkeypatch):
    ok, _ = _recover(tmp_path, monkeypatch, "initial.inp", OPT_INPUT, ["something unknown went wrong\n"])

    assert not ok


def test_an_output_recalc_kept_is_not_judged_again(tmp_path, monkeypatch):
    inp = tmp_path / "S1_T1_ISC_ms0.inp"
    inp.write_text(ESD_INPUT)
    out = tmp_path / "S1_T1_ISC_ms0.out"
    out.write_text(_esd(290.27, "5.698207e+07"))
    calls = []
    monkeypatch.setattr(orca, "run_orca", lambda *a, **k: calls.append(a) or True)   # recalc: kept, untouched

    ok = orca.run_orca_with_intelligent_recovery(str(inp), str(out), working_dir=tmp_path,
                                                 config={"enable_auto_recovery": "yes"})

    assert ok and len(calls) == 1
    assert not list(tmp_path.glob("*.retry*.inp"))


def test_the_esd_fix_reaches_every_job_with_an_esd_block(tmp_path):
    two = ESD_INPUT + "\n$new_job\n" + ESD_INPUT.replace("S1_T1_ISC_ms0", "S1_T1_ISC_ms0_b")
    inp = tmp_path / "S1_T1_ISC_ms0.inp"
    inp.write_text(two)
    out = tmp_path / "S1_T1_ISC_ms0.out"
    out.write_text(_esd(290.27, "1.0e+06"))
    strategy = RecoveryStrategy(E.ESD_WINDOW_TRUNCATED, 1, {})
    strategy.output_file = out
    strategy.job_index = 0

    new = OrcaInputModifier(inp, {}).apply_recovery(strategy).read_text()

    assert new.count("%ESD") + new.count("%esd") == 2
    assert "MAXTIME" not in new and "NPOINTS" not in new


# ------------------------------------------------------------ LEANSCF = the SCF

_MAIN_SCF_FAILED = """
                         GEOMETRY OPTIMIZATION CYCLE   1
--------------
SCF ITERATIONS
--------------
               *        SCF NOT CONVERGED AFTER 125 CYCLES         *
[file orca_leanscf/orca_leanscf.cpp, line 306]: Error (ORCA_LEANSCF): unfortunately, the SCF has not converged.
ORCA finished by error termination in LEANSCF
"""
_HESSIAN_SCF_FAILED = """
--------------
SCF ITERATIONS
--------------
                            ORCA SCF RESPONSE CALCULATION
[file orca_leanscf/orca_leanscf.cpp, line 306]: Error (ORCA_LEANSCF): unfortunately, the SCF has not converged.
ORCA finished by error termination in LEANSCF
"""


@pytest.mark.parametrize("attempt", [1, 2, 3, 4])
def test_an_scf_that_leanscf_could_not_converge_keeps_its_frequencies(tmp_path, attempt):
    """In ORCA 6 every SCF runs in orca_leanscf: all 472 archived LEANSCF
    failures were the main SCF (426 in the first optimisation cycle).  Their
    retries got VeryTightSCF, then had FREQ taken away -- neither helps an
    SCF converge, and the second leaves IMAG without a Hessian."""
    out = tmp_path / "job.out"
    out.write_text(_MAIN_SCF_FAILED)
    leanscf = RecoveryStrategy(E.LEANSCF_NOT_CONVERGED, attempt, {})
    leanscf.output_file = out
    scf = RecoveryStrategy(E.SCF_NO_CONVERGENCE, attempt, {})

    mods = leanscf.get_modifications()

    assert mods == scf.get_modifications()
    assert not mods.get("skip_freq")
    assert "VeryTightSCF" not in mods.get("keywords_add", [])


def test_an_scf_that_failed_inside_the_hessian_keeps_the_frequency_path(tmp_path):
    out = tmp_path / "job.out"
    out.write_text(_HESSIAN_SCF_FAILED)
    strategy = RecoveryStrategy(E.LEANSCF_NOT_CONVERGED, 3, {})
    strategy.output_file = out

    assert strategy.get_modifications().get("skip_freq") is True



def test_each_failed_run_keeps_its_own_backup(tmp_path, monkeypatch):
    # measured: an optimisation that ran out of cycles, then an SCF that did
    # not converge -- both "attempt 1" of their error type, one output.old1.out
    ok, _ = _recover(tmp_path, monkeypatch, "input5.inp", OPT_INPUT,
                     [_GAVE_UP, _MAIN_SCF_FAILED, _CONVERGED])

    assert ok
    assert "did not converge" in (tmp_path / "input5.old1.out").read_text()
    assert "LEANSCF" in (tmp_path / "input5.old2.out").read_text()

# ----------------------------------------------------------- refused input

_PARITY = """[file orca_main/main_util_tools.cpp, line 715]: Error : multiplicity (4) is even and number of electrons (242) is even -> impossible
[file orca_main/run.cpp, line 2]: ORCA finished with error return - aborting the run
"""


def test_an_input_orca_refuses_is_named_and_not_retried(tmp_path, monkeypatch, caplog):
    ok, fake = _recover(tmp_path, monkeypatch, "input9.inp", OPT_INPUT, [_PARITY])

    assert not ok
    assert len(fake.inputs) == 1, "a retry cannot make 242 electrons fit a quartet"
    assert "multiplicity (4) is even" in caplog.text
    assert OrcaErrorDetector.analyze_output(tmp_path / "input9.out") is E.INPUT_ERROR


# ------------------------------------------------------ rate jobs' SCF

def test_rate_jobs_converge_the_scf_like_the_states_they_start_from():
    from delfin.esd_input_generator import _rate_scf

    assert _rate_scf("! PBE0 def2-SVP ESD(ISC)").endswith(" TightSCF")
    assert _rate_scf("! PBE0 def2-SVP VeryTightSCF ESD(ISC)") == "! PBE0 def2-SVP VeryTightSCF ESD(ISC)"
