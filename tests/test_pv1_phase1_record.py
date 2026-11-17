"""PV1 Phase 1 — build, write, read a complete provenance record.

The module ``delfin.provenance_record`` does not exist yet; every test here
starts from the public functions ``build_record`` /
``write_record`` / ``read_record`` and asserts the record a finished run
gets is complete.  The run folder is SYNTHETIC (built under tmp_path): no
ORCA/TURBOMOLE-generated data is allowed in the repo, so the "program"
banner is a stub line shaped like an ORCA banner.

These tests are the red control: on the current commit they fail at import
of ``delfin.provenance_record``.
"""
from __future__ import annotations

import json
import re

from delfin.provenance_record import (
    PROVENANCE_FILE,
    SCHEMA_VERSION,
    build_record,
    read_record,
    write_record,
)

BANNER = (
    "                                 * O   R   C   A *\n"
    "\n"
    "           --- An Ab Initio, DFT and Semiempirical electronic structure package ---\n"
    "\n"
    "\n"
    "Program Version 6.1.1, authors: F. Neese, ...\n"
)
INP = (
    "! B3LYP def2-SVP\n"
    "* xyzfile 0 1\n"
    "mol.xyz\n"
    "*\n"
)
XYZ = (
    "3\n"
    "water\n"
    "O  0.0000  0.0000  0.1173\n"
    "H  0.0000  0.7572 -0.4692\n"
    "H  0.0000 -0.7572 -0.4692\n"
)


def _make_run_folder(tmp_path) -> str:
    """Build a synthetic finished run folder and return its path."""
    run = tmp_path / "run_h2o"
    run.mkdir()
    (run / "CONTROL.txt").write_text("method=orca\nsmiles=O\n", encoding="utf-8")
    (run / "input.inp").write_text(INP, encoding="utf-8")
    (run / "mol.xyz").write_text(XYZ, encoding="utf-8")
    (run / "input.out").write_text(BANNER, encoding="utf-8")
    # A DELFIN-style pipeline state that names the structure origin.
    (run / "DELFIN_Data.json").write_text(
        json.dumps({"input": {"smiles": "O"}}), encoding="utf-8"
    )
    (run / "final_energy.csv").write_text("E,hartree\n-76.5\n", encoding="utf-8")
    return str(run)


class TestBuildCompleteRecord:
    def test_public_functions_exist(self):
        assert callable(build_record)
        assert callable(write_record)
        assert callable(read_record)

    def test_schema_version_is_set(self):
        assert isinstance(SCHEMA_VERSION, str) and SCHEMA_VERSION

    def test_record_is_complete_for_synthetic_finished_run(self, tmp_path):
        run = _make_run_folder(tmp_path)
        record = build_record(run)

        # Required top-level keys — every target-1 field present.
        assert record["schema_version"] == SCHEMA_VERSION
        assert isinstance(record["input_files"], list) and record["input_files"]
        for f in record["input_files"]:
            assert set(f) >= {"path", "sha256"}
            # path is relative to the run, never absolute, never a home path.
            assert not f["path"].startswith("/")
            assert "~" not in f["path"]
        assert record["program"] in ("orca", "provenance_missing")
        assert record["program_version"] or record["program_version"] == "provenance_missing"
        assert "submit_command" in record
        assert "commit" in record["delfin"]
        assert "dirty" in record["delfin"]
        assert isinstance(record["python_packages"], dict)
        assert "structure_origin" in record
        assert "host_class" in record and record["host_class"] == "host_class_redacted"
        assert isinstance(record["timestamps"], dict)
        assert "created" in record["timestamps"]

    def test_program_and_version_read_from_banner(self, tmp_path):
        run = _make_run_folder(tmp_path)
        record = build_record(run)
        assert record["program"] == "orca"
        assert record["program_version"] == "6.1.1"

    def test_input_files_carry_real_hashes(self, tmp_path):
        run = _make_run_folder(tmp_path)
        record = build_record(run)
        by_path = {f["path"]: f["sha256"] for f in record["input_files"]}
        assert "input.inp" in by_path
        assert "mol.xyz" in by_path
        assert re.fullmatch(r"[0-9a-f]{64}", by_path["input.inp"])


class TestReadWrite:
    def test_write_then_read_roundtrip(self, tmp_path):
        run = _make_run_folder(tmp_path)
        record = build_record(run)
        written = write_record(record, run)
        assert written.name == PROVENANCE_FILE
        assert written.parent == tmp_path / "run_h2o"

        loaded = read_record(str(written))
        assert loaded["schema_version"] == SCHEMA_VERSION
        assert loaded["program_version"] == "6.1.1"
        assert loaded["structure_origin"] == record["structure_origin"]
