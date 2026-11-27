from pathlib import Path

from delfin import cli_recalc, orca, smart_recalc


def test_run_orca_wrapper_bootstraps_missing_fingerprint(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("DELFIN_RECALC", "1")
    monkeypatch.setenv("DELFIN_SMART_RECALC", "1")

    inp = tmp_path / "job.inp"
    out = tmp_path / "job.out"
    inp.write_text("! SP\n* xyz 0 1\nH 0.0 0.0 0.0\n*\n", encoding="utf-8")
    out.write_text("...\n****ORCA TERMINATED NORMALLY****\n", encoding="utf-8")

    def fail_run_orca(*args, **kwargs):
        raise AssertionError("run_orca should not be called for complete legacy output")

    monkeypatch.setattr(orca, "run_orca", fail_run_orca)
    wrappers, _ = cli_recalc.setup_recalc_mode()

    assert wrappers["run_orca"]("job.inp", "job.out") is True
    assert (tmp_path / "job.inp.fprint").exists()


def test_run_orca_wrapper_reruns_when_fingerprint_differs(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("DELFIN_RECALC", "1")
    monkeypatch.setenv("DELFIN_SMART_RECALC", "1")

    inp = tmp_path / "job.inp"
    out = tmp_path / "job.out"
    fprint = tmp_path / "job.inp.fprint"
    inp.write_text("! SP\n* xyz 0 1\nH 0.0 0.0 0.0\n*\n", encoding="utf-8")
    out.write_text("...\n****ORCA TERMINATED NORMALLY****\n", encoding="utf-8")
    fprint.write_text("stale-fingerprint\n", encoding="utf-8")

    calls = []

    def fake_run_orca(*args, **kwargs):
        calls.append((args, kwargs))
        return True

    monkeypatch.setattr(orca, "run_orca", fake_run_orca)
    wrappers, _ = cli_recalc.setup_recalc_mode()

    assert wrappers["run_orca"]("job.inp", "job.out") is True
    assert len(calls) == 1


def test_required_orca_outputs_does_not_require_gbw_for_generic_orca_job(tmp_path):
    inp = tmp_path / "job.inp"
    out = tmp_path / "job.out"
    inp.write_text("! SP\n* xyz 0 1\nH 0.0 0.0 0.0\n*\n", encoding="utf-8")
    out.write_text("...\n****ORCA TERMINATED NORMALLY****\n", encoding="utf-8")

    assert smart_recalc.required_orca_outputs(inp_path=inp, out_path=out) == []
    assert smart_recalc.outputs_complete(inp, out) is True


# A forced rerun (--recalc --occupier-override) recomputes a job that already
# finished.  The output it replaces is the record of that job: deleting it up
# front means a rerun that fails, or writes junk, leaves the junk in its place
# and the finished result is gone.  Measured through delfin.cli.main: the good
# output disappeared and the run still reported success.

GOOD_OUTPUT = ("* O R C A *\nProgram Version 6.1.1\n" + "SCF done\n" * 50
               + "                             ****ORCA TERMINATED NORMALLY****\n")


def _forced_wrapper(tmp_path, monkeypatch, rerun):
    """The wrapper cli.main installs for a forced output, around a fake rerun."""
    inp = tmp_path / "job.inp"
    out = tmp_path / "job.out"
    inp.write_text("! PBE0 def2-SVP\n* xyz 0 1\nH 0 0 0\n*\n", encoding="utf-8")
    out.write_text(GOOD_OUTPUT, encoding="utf-8")

    # setup_recalc_mode binds its "real" function at call time, which is how
    # cli.main composes it; the fake stands in for ORCA itself.
    monkeypatch.setattr(orca, "run_orca", rerun)
    wrappers, _ = cli_recalc.setup_recalc_mode(force_outputs={out.resolve()})
    return wrappers["run_orca"], inp, out


def test_a_failed_forced_rerun_leaves_the_previous_output_in_place(tmp_path, monkeypatch):
    """The rerun fails and its junk must not become the job's output."""
    def failing_rerun(inp_file, out_file, *args, **kwargs):
        Path(out_file).write_text("orca: segfault-ish junk, no terminator\n", encoding="utf-8")
        return False

    wrapper, inp, out = _forced_wrapper(tmp_path, monkeypatch, failing_rerun)

    assert wrapper(str(inp), str(out)) is False, "a failed rerun is not a success"
    assert out.is_file(), "the finished output must survive a failed rerun"
    assert out.read_text(encoding="utf-8") == GOOD_OUTPUT, \
        "the previous output is restored byte for byte, not left as the rerun's junk"


def test_a_forced_rerun_that_exits_zero_but_writes_junk_also_restores(tmp_path, monkeypatch):
    """ORCA can exit 0 and still leave nothing usable; the output is judged, not the code."""
    def junk_rerun(inp_file, out_file, *args, **kwargs):
        Path(out_file).write_text("orca: no terminator here\n", encoding="utf-8")
        return True

    wrapper, inp, out = _forced_wrapper(tmp_path, monkeypatch, junk_rerun)

    assert wrapper(str(inp), str(out)) is False
    assert out.read_text(encoding="utf-8") == GOOD_OUTPUT


def test_a_successful_forced_rerun_keeps_the_new_output_and_no_copy(tmp_path, monkeypatch):
    """Control case: a good rerun replaces the output and leaves no kept copy."""
    new_output = GOOD_OUTPUT.replace("SCF done", "SCF recomputed")

    def good_rerun(inp_file, out_file, *args, **kwargs):
        Path(out_file).write_text(new_output, encoding="utf-8")
        return True

    wrapper, inp, out = _forced_wrapper(tmp_path, monkeypatch, good_rerun)

    assert wrapper(str(inp), str(out)) is True
    assert out.read_text(encoding="utf-8") == new_output, "the rerun's output is the result"
    assert not list(tmp_path.glob("job.out.before-rerun*")), \
        "a complete rerun spends the kept copy instead of leaving a stale duplicate"

