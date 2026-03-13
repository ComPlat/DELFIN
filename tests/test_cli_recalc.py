from delfin import cli_recalc, orca


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
