"""PV1 Phase 1 — controls: what cannot be determined is recorded, never guessed.

Target 4 of the brief: a missing banner / a dirty tree / an unknown origin is
recorded as such (``provenance_missing`` / ``unknown``), never fabricated; and
no user name, host name or absolute home path may appear inside the record
(PR-public rule).  All run folders are synthetic, built under tmp_path — never
ORCA/TURBOMOLE-generated data.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

from delfin.provenance_record import MISSING, UNKNOWN_ORIGIN, build_record

BANNER_NO_VERSION = "                                 * O   R   C   A *\n\n"
GARBAGE_OUT = "this output has\nno banner at all\n"


def _run(tmp_path, out_text: str, with_control: bool = True,
         with_data_json: bool = True) -> Path:
    run = tmp_path / "run"
    run.mkdir(exist_ok=True)
    control = "method=orca\nsmiles=O\n" if with_control else "method=orca\n"
    (run / "CONTROL.txt").write_text(control, encoding="utf-8")
    (run / "input.inp").write_text("! B3LYP\n", encoding="utf-8")
    (run / "output.out").write_text(out_text, encoding="utf-8")
    if with_data_json:
        (run / "DELFIN_Data.json").write_text(
            json.dumps({"input": {"smiles": "O"}}), encoding="utf-8")
    return run


class TestMissingBanner:
    def test_orca_banner_without_version_is_recorded_missing(self, tmp_path):
        run = _run(tmp_path, BANNER_NO_VERSION)
        rec = build_record(str(run))
        assert rec["program"] == "orca"          # the banner names ORCA...
        assert rec["program_version"] == MISSING  # ...but never a made-up version

    def test_no_output_banner_is_recorded_missing(self, tmp_path):
        run = _run(tmp_path, GARBAGE_OUT)
        rec = build_record(str(run))
        assert rec["program"] == MISSING
        assert rec["program_version"] == MISSING


class TestUnknownOrigin:
    def test_no_smiles_and_no_structure_file_is_unknown(self, tmp_path):
        run = _run(tmp_path, GARBAGE_OUT, with_data_json=False)
        # remove the only structure input so nothing names the origin
        (run / "input.inp").unlink()
        rec = build_record(str(run))
        assert rec["structure_origin"]["kind"] == "unknown"
        assert rec["structure_origin"]["value"] == UNKNOWN_ORIGIN

    def test_smiles_origin_from_delfin_data(self, tmp_path):
        run = _run(tmp_path, GARBAGE_OUT)
        rec = build_record(str(run))
        assert rec["structure_origin"]["kind"] == "smiles"
        assert rec["structure_origin"]["smiles"] == "O"


class TestSubmitCommand:
    def test_no_submit_script_is_not_found_not_guessed(self, tmp_path):
        run = _run(tmp_path, GARBAGE_OUT)
        rec = build_record(str(run))
        assert rec["submit_command"]["script"] is None
        assert rec["submit_command"]["command"] == MISSING


class TestNoPrivateDataInRecord:
    def test_serialized_record_has_no_home_path(self, tmp_path):
        run = _run(tmp_path, GARBAGE_OUT)
        rec = build_record(str(run))
        text = json.dumps(rec)
        assert "~" not in text
        home = os.path.expanduser("~")
        if home and home != "~":
            # the record must not name the current user's home at all
            assert home.rstrip("/") not in text

    def test_record_has_no_absolute_host_appearance(self, tmp_path):
        run = _run(tmp_path, GARBAGE_OUT)
        rec = build_record(str(run))
        assert rec["host_class"] == "host_class_redacted"
        # every input file path is relative to the run
        for f in rec["input_files"]:
            assert not f["path"].startswith("/")
