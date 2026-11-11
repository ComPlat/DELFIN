"""PV1 phase 1 CLI: provenance_record.main() prints a record for ``delfin provenance``.

These tests exercise the public ``delfin provenance <result>`` path (the
``main`` function of delfin.provenance_record) against a synthetic finished run
folder -- no ORCA data. The red control for the package is the module-absent
collection error in test_pv1_phase1_record.py; here we assert the readable
entry point behaves.
"""

import tempfile
from pathlib import Path

from delfin.provenance_record import (
    build_record,
    main,
    write_record,
)


class TestMainPrintsRecord:
    def _make_run(self):
        run = Path(tempfile.mkdtemp()) / "run_h2o"
        run.mkdir(parents=True)
        (run / "CONTROL.txt").write_text("method=orca\nsmiles=O\n")
        (run / "input.inp").write_text("! B3LYP def2-SVP\n* xyz 0 1\nO 0 0 0\n*\n")
        (run / "mol.xyz").write_text("1\nwater\nO 0 0 0\n")
        (run / "output.out").write_text(
            "* O R C A *\nProgram Version 6.1.1, authors: F. Neese\n"
        )
        return run

    def test_main_prints_schema_and_program(self, capsys):
        run = self._make_run()
        record = build_record(run)
        write_record(record, run)
        rc = main([str(run)])
        out = capsys.readouterr().out
        assert rc == 0
        assert "schema_version" in out
        assert "provenance_missing" not in out  # banner present -> program is orca

    def test_main_json_flag_prints_raw_json(self, capsys):
        run = self._make_run()
        record = build_record(run)
        write_record(record, run)
        rc = main(["--json", str(run)])
        out = capsys.readouterr().out
        assert rc == 0
        assert '"schema_version"' in out
        assert '"host_class_redacted"' in out

    def test_main_missing_result_returns_two(self, capsys):
        missing = Path(tempfile.mkdtemp()) / "does_not_exist"
        rc = main([str(missing)])
        err = capsys.readouterr().err
        assert rc == 2
        assert "provenance unavailable" in err

    def test_main_accepts_file_path_directly(self, capsys):
        run = self._make_run()
        record = build_record(run)
        written = write_record(record, run)
        rc = main([str(written)])
        out = capsys.readouterr().out
        assert rc == 0
        assert "host_class_redacted" in out
