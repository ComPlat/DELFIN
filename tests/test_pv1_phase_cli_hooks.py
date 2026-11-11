"""PV1 protected-hook red control for target-2 hook #1 (delfin/cli.py).

The two target-2 hooks live in the protected CLI dispatch and are prepared for
the operator as ONE .gate/px orig/new pair (orig_cli.py -> new_cli.py, see
.gate/px/new_cli.py):
  1. ``delfin provenance <result>`` CLI dispatch -> delfin.provenance_record.main
     (new_cli.py:1590-1593) -- covered by the single test method below.
  2. run-finalize auto-write of PROVENANCE.json at the end-of-run completion
     block (new_cli.py:2284-2291, inside the ``_finalize`` closure). NOT driven
     from here: ``_finalize`` is only reachable through the full ``delfin.cli.main``
     orchestration (global job-manager init + signal handlers), so hook #2 is
     verified by the operator's build of the pair (build_record + write_record
     themselves are green-tested in test_pv1_phase1_record.py /
     test_pv1_adversarial_public_safety.py).

This test drives the top-level ``delfin.cli.main`` console dispatch -- NOT
``delfin.provenance_record.main`` directly (the on-demand path is already covered
by test_pv1_phase1_cli.py). It is RED on the unpatched delfin/cli.py (the
``provenance`` subcommand is not routed, so rc==2 "Unknown argument(s)") and
turns GREEN once the operator builds the pair.
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
