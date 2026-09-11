"""A recovery retry is an input ORCA reads, for the job that failed, with nothing else lost.

Measured before this was rebuilt: of 14 350 retries written for 543 archived
inputs (every strategy, every attempt), 7 519 were wrong, and ORCA 6.1.1
refused the ones tried at once:

* the OCCUPIER inputs' one-line ``%scf maxiter 125 end`` came back nested in
  a new ``%scf`` ("Unrecognized symbol in SCF block"); a second one-line
  ``%scf BrokenSym 1,1 end`` was lost or written twice;
* ``SOSCF`` alone opens a sub-block, so the settings after it were read
  inside it ("Unknown identifier in SCF::SOSCF block");
* ``Convergence TightSCF`` / ``VeryTightConv`` are not ORCA values ("Invalid
  assignment in SCF block"); ``Tight`` / ``VeryTight`` are;
* a metal's ``NewGTO ... end`` on its coordinate line was dropped when the
  coordinates were updated, so the retry ran the metal in another basis;
* in the archive (JRD1_test/ESD/T1), a T1 optimisation plus its TD-DFT check
  came back as one job with the second job's ``%base`` and ``%tddft``; ORCA
  stopped at once and the 191 MB output of the first attempt was replaced.

And a job without ``%base`` wrote its retry's gbw/xyz as ``<job>.retry1.*``,
while every later step read the failed run's ``<job>.*``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from delfin import orca
from delfin.common import orca_input as oi
from delfin.orca_recovery import OrcaErrorType as E
from delfin.orca_recovery import OrcaInputModifier, RecoveryStrategy

OCCUPIER_BS = """! PBE0 def2-SVP D4 RIJCOSX def2/J CPCM(dmf) OPT FREQ PModel MOREAD TightSCF
%moinp "input_red_step_1_OCCUPIER.gbw"
%maxcore 6000
%pal nprocs 18 end
%scf maxiter 125 end
%scf BrokenSym 1,1 end
%freq
  Temp 298.15
end
* xyz -1 1
C          -1.11736974265122      3.57845588429044      0.53176128123210
Co         -0.02912657244695     -0.34055506228113     -0.00413799049569   NewGTO "def2-TZVP" end
N           1.18159758310338     -1.73604082305673      0.43982649453832
*
"""

T1_WITH_CHECK = """! CAM-B3LYP UKS def2-TZVPP D4 RIJCOSX def2/J CPCM(chcl3) OPT MOREAD
%base "T1"
%moinp "S0.gbw"
%pal nprocs 20 end
%maxcore 6000
* xyz 0 3
C  0.0 0.0 0.0
O  0.0 0.0 1.2
*

$new_job
! CAM-B3LYP RKS def2-TZVPP D4 RIJCOSX def2/J CPCM(chcl3)
%base "T1_TDDFT"
%pal nprocs 20 end
%maxcore 6000

%tddft
  nroots 15
  tda FALSE
  triplets true
end

* xyzfile 0 1 T1.xyz
"""

DELTA_SCF = """! CAM-B3LYP UKS 6-31+G* D4 RIJCOSX def2/J CPCM(chcl3) OPT deltaSCF MOREAD FreezeAndRelease
%base "S1_second_deltaSCF"
%moinp "S1_first_TDDFT.gbw"
%pal nprocs 20 end
%maxcore 6000
%scf
  DOMOM true
  alphaconf 0,1
  betaconf 0
  SOSCFHESSUP LSR1
  maxiter 300
end

* xyzfile 0 1 S1_first_TDDFT.xyz
"""

PLAN = [(E.SCF_NO_CONVERGENCE, a) for a in (1, 2, 3)] + [(E.LEANSCF_NOT_CONVERGED, a) for a in (1, 2, 3, 4)] + \
       [(E.TRAH_SEGFAULT, a) for a in (1, 2)] + [(E.DIIS_ERROR, a) for a in (1, 2, 3)] + \
       [(E.GEOMETRY_NOT_CONVERGED, a) for a in (1, 2, 3)] + [(E.MPI_CRASH, a) for a in (1, 2, 3)] + \
       [(E.CIS_FAILURE, a) for a in (1, 2, 3)] + [(E.FREQUENCY_FAILURE, a) for a in (1, 2)] + \
       [(E.MEMORY_ERROR, a) for a in (1, 2)]


def _retry(tmp_path: Path, name: str, text: str, error, attempt: int, *, failed_job: int = 1,
           own_files: dict | None = None, memory_note: str = "") -> tuple[str, dict]:
    inp = tmp_path / name
    inp.write_text(text)
    for fname, content in (own_files or {}).items():
        path = tmp_path / fname
        path.write_bytes(content) if isinstance(content, bytes) else path.write_text(content)
    out = tmp_path / (inp.stem + ".out")
    out.write_text("".join(f"\n $$$$$$$$$$$$$$$$  JOB NUMBER  {j} $$$$$$$$$$$$$$\n"
                           for j in range(1, failed_job + 1)) + memory_note)
    strategy = RecoveryStrategy(error, attempt, {})
    strategy.output_file = out
    new = OrcaInputModifier(inp, {}).apply_recovery(strategy)
    assert new != inp, "no retry was written"
    return new.read_text(), strategy.get_modifications()


def _scf(job_text):
    return {st.key: " ".join(st.lines[0].split()[1:]) for it in oi.parse_job(job_text)
            if it.kind == "block" and it.name == "scf" for st in it.block.statements if st.key}


@pytest.mark.parametrize("name, text", [("red_step_1.inp", OCCUPIER_BS), ("T1.inp", T1_WITH_CHECK),
                                        ("S1_second.inp", DELTA_SCF)])
@pytest.mark.parametrize("error, attempt", PLAN)
def test_every_retry_reads_as_orca_input_and_keeps_every_job(tmp_path, name, text, error, attempt):
    retry, _ = _retry(tmp_path, name, text, error, attempt)
    jobs = oi.split_jobs(retry)

    assert len(jobs) == len(oi.split_jobs(text))
    for job, original in zip(jobs, oi.split_jobs(text)):
        items = oi.parse_job(job)                      # raises on anything ORCA could not read
        base = oi.setting_value(items, "base")
        assert base is not None
        own = oi.setting_value(oi.parse_job(original), "base")
        assert base == (own or f'"{Path(name).stem}"')
        assert oi.has_block(items, "tddft") == oi.has_block(oi.parse_job(original), "tddft")
        for key, value in _scf(job).items():
            assert value.split()[0] not in ("TightSCF", "VeryTightConv")
            assert not key.startswith("sub "), "a bare word in %scf opens a sub-block"


@pytest.mark.parametrize("error, attempt", PLAN)
def test_broken_symmetry_and_the_metal_basis_survive_every_retry(tmp_path, error, attempt):
    retry, _ = _retry(tmp_path, "red_step_1.inp", OCCUPIER_BS, error, attempt,
                      own_files={"red_step_1.gbw": b"\0" * 150_000,
                                 "red_step_1.xyz": "3\n\nC 1.0 2.0 3.0\nCo 0.1 0.2 0.3\nN -1 -2 -3\n"})

    assert _scf(retry).get("brokensym") == "1,1"
    assert "maxiter" in _scf(retry)
    co = next(line for line in retry.splitlines() if line.split()[:1] == ["Co"])
    assert co.endswith('NewGTO "def2-TZVP" end')


def test_the_scf_escalation_writes_orcas_own_switches(tmp_path):
    retry, _ = _retry(tmp_path, "red_step_1.inp", OCCUPIER_BS, E.SCF_NO_CONVERGENCE, 3)
    scf = _scf(retry)

    assert scf["cnvsoscf"] == "true" and scf["maxiter"] == "800" and scf["dampfac"] == "0.95"
    tddft_retry, _ = _retry(tmp_path, "T1.inp", T1_WITH_CHECK, E.CIS_FAILURE, 2, failed_job=2)
    assert _scf(oi.split_jobs(tddft_retry)[1])["convergence"] == "Tight"


def test_a_failed_second_job_is_the_one_changed(tmp_path):
    retry, _ = _retry(tmp_path, "T1.inp", T1_WITH_CHECK, E.CIS_FAILURE, 1, failed_job=2)
    first, second = oi.split_jobs(retry)

    assert "tda true" in second
    assert "%tddft" not in first, "%tddft in an optimisation makes ORCA optimise an excited state"
    assert '%base "T1"' in first and '%base "T1_TDDFT"' in second
    assert "* xyzfile 0 1 T1.xyz" in second


def test_an_scf_failure_in_the_first_job_leaves_the_second_alone(tmp_path):
    retry, _ = _retry(tmp_path, "T1.inp", T1_WITH_CHECK, E.SCF_NO_CONVERGENCE, 2, failed_job=1)
    first, second = oi.split_jobs(retry)

    assert "KDIIS" in first.splitlines()[0]
    assert second == oi.split_jobs(T1_WITH_CHECK)[1]


def test_a_retry_writes_its_files_under_the_jobs_own_name(tmp_path):
    retry, _ = _retry(tmp_path, "input4.inp", OCCUPIER_BS.replace("red_step_1", "input4"), E.GEOMETRY_NOT_CONVERGED, 1)

    assert oi.setting_value(oi.parse_job(retry), "base") == '"input4"'


def test_a_retry_of_a_retry_still_writes_under_the_original_name(tmp_path):
    (tmp_path / "input4.retry1.inp").write_text(OCCUPIER_BS)
    strategy = RecoveryStrategy(E.SCF_NO_CONVERGENCE, 2, {})
    new = OrcaInputModifier(tmp_path / "input4.retry1.inp", {}).apply_recovery(strategy)

    assert oi.setting_value(oi.parse_job(new.read_text()), "base") == '"input4"'


def test_a_job_continues_from_its_own_orbitals_and_never_from_another_jobs(tmp_path):
    own, _ = _retry(tmp_path, "S1_second.inp", DELTA_SCF, E.SCF_NO_CONVERGENCE, 1,
                    own_files={"S1_second_deltaSCF.gbw": b"\0" * 150_000})
    assert '%moinp "S1_second_deltaSCF_old.gbw"' in own

    other = tmp_path / "other"
    other.mkdir()
    kept, _ = _retry(other, "S1_second.inp", DELTA_SCF, E.SCF_NO_CONVERGENCE, 1,
                     own_files={"S1_first_TDDFT.gbw": b"\0" * 150_000, "T1.gbw": b"\0" * 150_000})
    assert '%moinp "S1_first_TDDFT.gbw"' in kept and "T1" not in kept

    gone = tmp_path / "gone"
    gone.mkdir()
    fresh, _ = _retry(gone, "S1_second.inp", DELTA_SCF, E.SCF_NO_CONVERGENCE, 1)
    assert "%moinp" not in fresh and "MOREAD" not in fresh.upper().split("\n")[0].split()


def test_an_optimisation_continues_from_its_own_last_geometry_only(tmp_path):
    retry, _ = _retry(tmp_path, "red_step_1.inp", OCCUPIER_BS, E.GEOMETRY_NOT_CONVERGED, 1,
                      own_files={"red_step_1.xyz": "3\n\nC 1.0 2.0 3.0\nCo 0.1 0.2 0.3\nN -1 -2 -3\n"})
    assert "0.10000000" in retry

    other = tmp_path / "other"
    other.mkdir()
    kept, _ = _retry(other, "red_step_1.inp", OCCUPIER_BS, E.GEOMETRY_NOT_CONVERGED, 1,
                     own_files={"red_step_1.xyz": "3\n\nC 1.0 2.0 3.0\nNi 0.1 0.2 0.3\nN -1 -2 -3\n"})
    assert "-0.02912657244695" in kept, "an xyz of other atoms is not this job's geometry"


def test_frequencies_are_retried_numerically(tmp_path):
    retry, mods = _retry(tmp_path, "red_step_1.inp", OCCUPIER_BS, E.FREQUENCY_FAILURE, 1)
    words = retry.splitlines()[0].split()

    assert mods["freq_method"] == "NumFreq"
    assert "NumFreq" in words and "FREQ" not in words


def test_maxcore_is_raised_to_what_orca_asked_for_but_never_lowered(tmp_path):
    note = "       ====>        Please increase MaxCore to more than:         29.9 MB\n"
    kept, _ = _retry(tmp_path, "red_step_1.inp", OCCUPIER_BS, E.MEMORY_ERROR, 1, memory_note=note)
    assert "%maxcore 6000" in kept

    big = OCCUPIER_BS.replace("%maxcore 6000", "%maxcore 10")
    raised, _ = _retry(tmp_path, "red_step_1.inp", big, E.MEMORY_ERROR, 1, memory_note=note)
    assert "%maxcore 45" in raised


def test_the_mpi_settings_of_a_crash_fix_reach_orca(tmp_path, monkeypatch):
    inp = tmp_path / "job.inp"
    inp.write_text(OCCUPIER_BS)
    seen = []

    def fake_run(path, output_log, timeout=None, **kwargs):
        seen.append(kwargs.get("extra_env"))
        Path(output_log).write_text("mpirun noticed that process rank 3 exited on signal 9 (Killed).\n"
                                    "ORCA finished by error termination in SCF\n" if len(seen) == 1 else
                                    "****ORCA TERMINATED NORMALLY****\n")
        return len(seen) > 1

    monkeypatch.setattr(orca, "run_orca", fake_run)
    monkeypatch.setattr(orca.OrcaErrorDetector, "analyze_output", classmethod(
        lambda cls, out: E.MPI_CRASH if "error termination" in Path(out).read_text() else None))
    ok = orca.run_orca_with_intelligent_recovery(str(inp), str(tmp_path / "job.out"), working_dir=tmp_path,
                                                 config={"enable_auto_recovery": "yes"})

    assert ok
    assert seen[0] is None
    assert seen[1]["OMPI_MCA_btl_vader_single_copy_mechanism"] == "none"


def test_input_that_cannot_be_taken_apart_is_not_retried(tmp_path, caplog):
    inp = tmp_path / "c.inp"
    inp.write_text("%compound\n  New_Step\n  ! B3LYP\n  Step_End\nend\n")
    with caplog.at_level(logging.ERROR, logger="delfin.orca_recovery"):
        new = OrcaInputModifier(inp, {}).apply_recovery(RecoveryStrategy(E.SCF_NO_CONVERGENCE, 1, {}))

    assert new == inp
    assert "cannot rewrite" in caplog.text


# ---------------------------------------------------------------- coverage

def test_every_step_of_a_workflow_runs_through_the_recovery():
    """enable_auto_recovery=yes (4 265 archived CONTROL files) reached only the
    ESD jobs and the OCCUPIER runs; initial/ox/red of classic and after
    OCCUPIER, and IMAG's jobs, ran without it -- a failed SCF ended the step."""
    import ast
    import inspect

    from delfin.workflows.engine import classic, occupier

    for module in (classic, occupier):
        tree = ast.parse(inspect.getsource(module))
        calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)]
        assert not [c for c in calls if c.func.id == "run_orca"], module.__name__
        recovered = [c for c in calls if c.func.id == "run_orca_with_intelligent_recovery"]
        assert recovered, module.__name__
        assert all(any(k.arg == "config" for k in c.keywords) for c in recovered), module.__name__


def test_imag_jobs_run_through_the_recovery_with_the_runs_settings(tmp_path, monkeypatch):
    from delfin import imag

    seen = {}

    def fake(inp, out, timeout=None, **kwargs):
        seen.update(kwargs)
        return True

    monkeypatch.setattr(orca, "run_orca_with_intelligent_recovery", fake)
    config = {"enable_auto_recovery": "yes"}
    assert imag._pipeline_run_orca(tmp_path / "a.inp", tmp_path / "a.out", working_dir=tmp_path, config=config)
    assert seen["config"] is config and seen["isolate"] is True
