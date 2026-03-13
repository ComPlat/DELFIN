from delfin import xtb_crest


def test_goat_bootstraps_missing_fingerprint_from_complete_output(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("DELFIN_RECALC", "1")
    monkeypatch.setenv("DELFIN_SMART_RECALC", "1")

    input_path = tmp_path / "input.txt"
    input_path.write_text("H 0.0 0.0 0.0\n", encoding="utf-8")

    goat_dir = tmp_path / "XTB2_GOAT"
    goat_dir.mkdir()
    (goat_dir / "output_XTB_GOAT.out").write_text(
        "header\n****ORCA TERMINATED NORMALLY****\n",
        encoding="utf-8",
    )
    (goat_dir / "XTB_GOAT.globalminimum.xyz").write_text(
        "1\ncomment\nH 0.0 0.0 0.0\n",
        encoding="utf-8",
    )

    def fail_run_orca(*args, **kwargs):
        raise AssertionError("GOAT should not rerun for a complete legacy output")

    monkeypatch.setattr(xtb_crest, "run_orca", fail_run_orca)

    xtb_crest.XTB_GOAT(
        multiplicity=1,
        charge=0,
        config={"xTB_method": "XTB2", "PAL": 4, "input_file": "input.txt"},
    )

    assert (goat_dir / "XTB_GOAT.inp.fprint").exists()
    assert input_path.read_text(encoding="utf-8").strip() == "H 0.0 0.0 0.0"
