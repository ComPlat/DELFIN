"""PV1 protected-hook red controls for the two target-2 hooks (delfin/cli.py).

The two hooks live in the protected CLI dispatch and are prepared for the
operator as ONE .gate/px orig/new pair (orig_cli.py -> new_cli.py):
  1. ``delfin provenance <result>`` CLI dispatch -> delfin.provenance_record.main
  2. run-finalize auto-write of PROVENANCE.json at the end-of-run completion block
     (right after the DELFIN_Data.json write).

These tests drive the top-level ``delfin.cli.main`` console dispatch -- NOT
``delfin.provenance_record.main`` directly (the on-demand path is already covered
by test_pv1_phase1_cli.py). They are RED on the unpatched delfin/cli.py (the
``provenance`` subcommand is not routed) and turn GREEN once the operator builds
the pair.
"""

import tempfile
from pathlib import Path

from delfin import cli as delfin_cli
from delfin.provenance_record import build_record, write_record


class TestProvenanceCliDispatch:
    def _make_run(self) -> Path:
        """Synthetic finished run folder (banner -> program orca), no ORCA data."""
        run = Path(tempfile.mkdtemp()) / "run_h2o"
        run.mkdir(parents=True)
        (run / "CONTROL.txt").write_text("method=orca\nsmiles=O\n")
        (run / "input.inp").write_text("! B3LYP def2-SVP\n* xyz 0 1\nO 0 0 0\n*\n")
        (run / "mol.xyz").write_text("1\nwater\nO 0 0 0\n")
        (run / "output.out").write_text(
            "* O R C A *\nProgram Version 6.1.1, authors: F. Neese\n"
        )
        # On-demand record is pre-written so the dispatch's job is only to route
        # to provenance_record.main and print the existing record.
        write_record(build_record(run), run)
        return run

    def test_provenance_subcommand_routes_to_record(self, capsys):
        """RED: delfin.cli.main does not yet dispatch ``provenance``."""
        run = self._make_run()
        rc = delfin_cli.main(["provenance", str(run)])
        out = capsys.readouterr().out
        assert rc == 0
        assert "schema_version" in out
        assert "structure_origin" in out
