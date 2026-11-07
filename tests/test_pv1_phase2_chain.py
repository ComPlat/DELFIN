"""PV1 phase 2: a result links back through its chain of ancestor step ids.

Target 3: links a result back through a chain (step ids), so a plot can name
every calculation it came from.  These are the red controls for the
implementation in the next commit -- on the current module compute_step_id
does not exist, so collection fails here.
"""

import json
from pathlib import Path

from delfin.provenance_record import (
    MISSING,
    build_record,
    compute_step_id,
    read_record,
    write_record,
)


def _make_run(root: Path, name: str = "run_calc") -> Path:
    run = Path(root) / name
    run.mkdir(parents=True, exist_ok=True)
    (run / "CONTROL.txt").write_text("method=orca\nsmiles=O\n")
    (run / "input.inp").write_text("! B3LYP def2-SVP\n* xyz 0 1\nO 0 0 0\n*\n")
    (run / "mol.xyz").write_text("1\nwater\nO 0 0 0\n")
    (run / "output.out").write_text("* O R C A *\nProgram Version 6.1.1\n")
    return run


class TestStepId:
    def test_step_id_is_stable_across_rebuilds(self, tmp_path):
        run = _make_run(tmp_path, "a")
        rec1 = build_record(str(run))
        rec2 = build_record(str(run))
        assert rec1["step_id"] == rec2["step_id"]
        assert isinstance(rec1["step_id"], str) and len(rec1["step_id"]) >= 16

    def test_step_id_differs_when_inputs_differ(self, tmp_path):
        run_a = _make_run(tmp_path, "a")
        run_b = _make_run(tmp_path, "b")
        (run_b / "input.inp").write_text("! PBE0 def2-TZVP\n* xyz 0 1\nO 0 0 0\n*\n")
        a = compute_step_id(build_record(str(run_a)))
        b = compute_step_id(build_record(str(run_b)))
        assert a != b

    def test_result_type_is_part_of_the_id(self, tmp_path):
        run = _make_run(tmp_path, "a")
        calc = compute_step_id(build_record(str(run), result_type="calculation"))
        plot = compute_step_id(build_record(str(run), result_type="plot"))
        assert calc != plot

    def test_unknown_missing_sentinel_is_exported(self):
        assert isinstance(MISSING, str) and MISSING


class TestChain:
    def test_build_with_parents_records_their_step_ids(self, tmp_path):
        run_a = _make_run(tmp_path, "a")
        rec_a = build_record(str(run_a))
        run_b = _make_run(tmp_path, "b")
        rec_b = build_record(str(run_b))
        write_record(rec_a, run_a)
        write_record(rec_b, run_b)

        run_plot = _make_run(tmp_path, "plot")
        rec_plot = build_record(str(run_plot), parents=[rec_a, rec_b], result_type="plot")
        assert rec_plot["result"]["type"] == "plot"
        assert rec_plot["chain"] == [rec_a["step_id"], rec_b["step_id"]]

    def test_chain_is_transitive_and_deduped(self, tmp_path):
        run_x = _make_run(tmp_path, "x")
        rec_x = build_record(str(run_x))
        run_y = _make_run(tmp_path, "y")
        rec_y = build_record(str(run_y), parents=[rec_x])
        run_z = _make_run(tmp_path, "z")
        rec_z = build_record(str(run_z), parents=[rec_y, rec_x])
        chain = rec_z["chain"]
        assert rec_x["step_id"] in chain
        assert rec_y["step_id"] in chain
        assert len(chain) == len(set(chain))

    def test_plot_roundtrip_names_every_source_calc(self, tmp_path):
        run_a = _make_run(tmp_path, "a")
        rec_a = build_record(str(run_a))
        run_plot = _make_run(tmp_path, "plot")
        rec_plot = build_record(str(run_plot), parents=[rec_a], result_type="plot")
        written = write_record(rec_plot, run_plot)
        loaded = read_record(str(written))
        assert loaded["chain"] == [rec_a["step_id"]]
        assert loaded["step_id"] == rec_plot["step_id"]


class TestNoChain:
    def test_single_calculation_has_empty_chain(self, tmp_path):
        run = _make_run(tmp_path, "a")
        rec = build_record(str(run))
        assert rec["chain"] == []

    def test_no_absolute_paths_in_input_files(self, tmp_path):
        run = _make_run(tmp_path, "a")
        rec = build_record(str(run))
        for f in rec["input_files"]:
            assert not f["path"].startswith("/")
            assert "~" not in f["path"]


class TestSchema:
    def test_schema_version_still_present(self, tmp_path):
        run = _make_run(tmp_path, "a")
        rec = build_record(str(run))
        assert rec["schema_version"] is not None


def test_full_chain_serialises_without_losing_step_ids(tmp_path):
    run_a = _make_run(tmp_path, "a")
    rec_a = build_record(str(run_a))
    run_plot = _make_run(tmp_path, "plot")
    rec_plot = build_record(str(run_plot), parents=[rec_a], result_type="plot")
    text = json.dumps(rec_plot)
    assert rec_a["step_id"] in text
    assert rec_plot["step_id"] in text
