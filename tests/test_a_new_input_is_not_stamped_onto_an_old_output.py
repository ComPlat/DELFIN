"""A recalc does not stamp a new input onto an output computed from another.

A finished job without a fingerprint (every archived run before the
fingerprints, and some after) was given one on recalc whenever its output
was complete.  The recalc wrapper did that without asking whether the output
belonged to the input beside it -- ORCA copies its input into the head of
the output, so it can be asked.  Measured in a dry run over Jerome's
archived jobs: after functional=B3LYP, the FoBs of complexes_1_200/7-26_23
were written anew with B3LYP and then skipped, their PBE0 results kept under
a B3LYP fingerprint.  orca.run_orca already asked; the wrapper, which decides
first in a recalc, now does too after an edit.  An output a recovery retry of
the job wrote still counts, and without an edit a finished output stays as it
always did.
"""

from __future__ import annotations

from delfin import cli_recalc, orca, smart_recalc

INPUT = "! {functional} def2-SVP\n* xyz 0 1\nH 0.0 0.0 0.0\n*\n"


def _output(functional: str) -> str:
    echo = "".join(f"|{i:3d}> {line}\n" for i, line in enumerate(INPUT.format(functional=functional).splitlines(), 1))
    return "INPUT FILE\n" + echo + "|  5> ****END OF INPUT****\n...\n****ORCA TERMINATED NORMALLY****\n"


def _wrapper(tmp_path, monkeypatch, edited=True):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("DELFIN_RECALC", "1")
    monkeypatch.setenv("DELFIN_SMART_RECALC", "1")
    monkeypatch.setenv("DELFIN_RECALC_CONTROL_EDITED", "1" if edited else "0")
    ran = []
    monkeypatch.setattr(orca, "run_orca", lambda *a, **k: ran.append(a) or True)
    wrappers, _ = cli_recalc.setup_recalc_mode()
    return wrappers["run_orca"], ran


def test_an_output_from_another_functional_is_computed_again(tmp_path, monkeypatch):
    run, ran = _wrapper(tmp_path, monkeypatch)
    (tmp_path / "job.inp").write_text(INPUT.format(functional="B3LYP"))
    (tmp_path / "job.out").write_text(_output("PBE0"))

    assert run("job.inp", "job.out") is True
    assert len(ran) == 1
    assert not (tmp_path / "job.inp.fprint").exists()


def test_an_output_of_this_input_is_kept_and_stamped(tmp_path, monkeypatch):
    run, ran = _wrapper(tmp_path, monkeypatch)
    (tmp_path / "job.inp").write_text(INPUT.format(functional="PBE0"))
    (tmp_path / "job.out").write_text(_output("PBE0"))

    assert run("job.inp", "job.out") is True
    assert ran == []
    assert (tmp_path / "job.inp.fprint").exists()


def test_an_output_a_retry_of_the_job_wrote_is_kept(tmp_path, monkeypatch):
    run, ran = _wrapper(tmp_path, monkeypatch)
    (tmp_path / "job.inp").write_text(INPUT.format(functional="PBE0"))
    (tmp_path / "job.retry1.inp").write_text(INPUT.format(functional="PBE0 SlowConv"))
    (tmp_path / "job.out").write_text(_output("PBE0 SlowConv"))

    assert run("job.inp", "job.out") is True
    assert ran == []
    assert smart_recalc.output_belongs_to_job(tmp_path / "job.inp", tmp_path / "job.out")


def test_without_an_edit_a_finished_output_stays_as_it_always_did(tmp_path, monkeypatch):
    # DELFIN writing an input a little differently since is no reason to
    # compute a finished job again
    run, ran = _wrapper(tmp_path, monkeypatch, edited=False)
    (tmp_path / "job.inp").write_text(INPUT.format(functional="PBE0") + "%scf maxiter 150 end\n")
    (tmp_path / "job.out").write_text(_output("PBE0"))

    assert run("job.inp", "job.out") is True
    assert ran == []
