"""PV1 phase 4 red controls: CONTROL-driven inputs, no generated outputs, --verify.

Synthetic run layouts written by hand (no ORCA output in the repo).  Controls
the operator's real-run findings (PR #625):
 1. input_files must carry the file CONTROL.txt names (input_file=input.txt),
    not only what matches a suffix list;
 2. ORCA-generated outputs (initial_trj.xyz etc.) must never appear as inputs;
 3. `main(["--verify", folder])` recomputes every stored sha256 and the banner
    version, reports each mismatch and exits non-zero when one is found.
"""
import json

import pytest

from delfin.provenance_record import (
    PROVENANCE_FILE,
    SCHEMA_VERSION,
    build_record,
    main,
    read_record,
    write_record,
)


def make_run(tmp_path, banner=True):
    """A hand-written run layout: CONTROL naming input.txt + a generated xyz."""
    run = tmp_path / "run"
    run.mkdir()
    (run / "input.txt").write_text("H 0 0 0\nO 0 0 0.96\nH 0 0.92 1.51\n")
    (run / "CONTROL.txt").write_text("input_file=input.txt\nNAME=rr\n")
    # DELFIN/ORCA-generated artefacts that must never be listed as inputs.
    (run / "initial.xyz").write_text("H 0 0 0\nO 0 0 0.96\nH 0 0.92 1.51\n")
    (run / "initial_trj.xyz").write_text("trajectory frames\n")
    out = run / "initial.out"
    out.write_text(
        "* O   R   C   A *\nProgram Version 6.1.0\n\n"
        if banner
        else "not an orca output\n"
    )
    return run


def test_named_input_file_is_recorded(tmp_path):
    run = make_run(tmp_path)
    record = build_record(run)
    paths = [f["path"] for f in record["input_files"]]
    assert "CONTROL.txt" in paths
    # (1) the file CONTROL names must be an input even though .txt is not in
    # any suffix list.
    assert "input.txt" in paths


def test_generated_outputs_are_not_inputs(tmp_path):
    run = make_run(tmp_path)
    record = build_record(run)
    paths = [f["path"] for f in record["input_files"]]
    # (2) initial.xyz / initial_trj.xyz are ORCA/DELFIN-generated, not inputs.
    assert "initial.xyz" not in paths
    assert "initial_trj.xyz" not in paths


def test_structure_origin_points_to_named_input(tmp_path):
    run = make_run(tmp_path)
    record = build_record(run)
    origin = record["structure_origin"]
    assert origin["kind"] == "file"
    assert origin["path"] == "input.txt"


def test_verify_passes_on_untouched_run(tmp_path, capsys):
    run = make_run(tmp_path)
    write_record(build_record(run), run)
    rc = main(["--verify", str(run)])
    assert rc == 0


def test_verify_reports_tampered_input(tmp_path, capsys):
    run = make_run(tmp_path)
    write_record(build_record(run), run)
    (run / "input.txt").write_text("tampered\n")  # CONTROL.txt edited after
    rc = main(["--verify", str(run)])
    out = capsys.readouterr().out + capsys.readouterr().err
    assert rc != 0
    assert "input.txt" in out
    assert "sha256" in out.lower()


def test_verify_reports_missing_banner_version(tmp_path, capsys):
    run = make_run(tmp_path, banner=False)
    write_record(build_record(run), run)
    rc = main(["--verify", str(run)])
    assert rc == 0  # recorded as provenance_missing, nothing to contradict
