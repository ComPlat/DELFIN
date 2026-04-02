from delfin import orca


def test_run_orca_triggers_auto_nmr_report_after_success(monkeypatch, tmp_path):
    inp = tmp_path / "job.inp"
    out = tmp_path / "job.out"
    inp.write_text("! NMR\n%EPRNMR\nEND\n", encoding="utf-8")

    monkeypatch.setattr(orca, "find_orca_executable", lambda: "/opt/orca/orca")
    monkeypatch.setattr(orca, "_apply_control_overrides_to_input", lambda *args, **kwargs: None)
    monkeypatch.setattr(orca.smart_recalc, "should_skip", lambda *args, **kwargs: False)
    monkeypatch.setattr(orca.smart_recalc, "recalc_enabled", lambda: False)
    monkeypatch.setattr(orca.smart_recalc, "smart_mode_enabled", lambda: False)
    monkeypatch.setattr(orca, "_run_orca_subprocess", lambda *args, **kwargs: True)

    stored = []
    post = []
    monkeypatch.setattr(orca.smart_recalc, "store_fingerprint", lambda path, extra_deps=None: stored.append(path))
    monkeypatch.setattr(orca, "_auto_generate_nmr_report", lambda input_path, output_path: post.append((input_path, output_path)) or True)

    assert orca.run_orca(str(inp), str(out)) is True
    assert stored == [inp.resolve()]
    assert post == [(inp.resolve(), out.resolve())]


def test_run_orca_skip_path_still_generates_missing_nmr_report(monkeypatch, tmp_path):
    inp = tmp_path / "job.inp"
    out = tmp_path / "job.out"
    inp.write_text("! NMR\n%EPRNMR\nEND\n", encoding="utf-8")
    out.write_text("...\n****ORCA TERMINATED NORMALLY****\n", encoding="utf-8")

    monkeypatch.setattr(orca, "find_orca_executable", lambda: "/opt/orca/orca")
    monkeypatch.setattr(orca, "_apply_control_overrides_to_input", lambda *args, **kwargs: None)
    monkeypatch.setattr(orca.smart_recalc, "should_skip", lambda *args, **kwargs: True)

    post = []
    monkeypatch.setattr(orca, "_auto_generate_nmr_report", lambda input_path, output_path: post.append((input_path, output_path)) or True)

    assert orca.run_orca(str(inp), str(out)) is True
    assert post == [(inp.resolve(), out.resolve())]
