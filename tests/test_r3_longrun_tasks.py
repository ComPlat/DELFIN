"""Gate tests for the r3 long-task suite (package R3).

Covers, all through the real load path and the real accept check():
* tasks_r3_longrun.yaml loads through the benchmark Task loader; all ten
  ids are present; every setup/verify name resolves to an existing file.
* the "40+ steps" property is ASSERTED, not assumed: every task carries
  max_tool_calls >= 120 and a thresholds.turns_min >= 40.
* the setup script builds every fixture folder in a tmp_path.
* the accept script ACCEPTS the correct artifact and REJECTS each wrong
  variant (missing artifact, wrong value, wrong decision) -- one accept
  and one reject per family minimum, in practice per task.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

BENCH = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack" / "benchmark"

ALL_IDS = [
    "r3long-code-median", "r3long-code-cache", "r3long-code-parse",
    "r3long-chem-control-fill", "r3long-chem-solvent-check",
    "r3long-absent-binary",
    "r3long-data-mean", "r3long-data-join", "r3long-data-exclude",
    "r3long-data-rank",
]


def _load_tasks():
    from delfin.agent.benchmark import load_tasks
    return load_tasks(BENCH / "tasks_r3_longrun.yaml")


def _accept_mod():
    path = BENCH / "accept" / "r3long_accept.py"
    spec = importlib.util.spec_from_file_location("r3long_accept", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _setup_mod():
    sys.path.insert(0, str(BENCH / "setup"))
    try:
        import r3long_setup as setup_mod
    finally:
        sys.path.remove(str(BENCH / "setup"))
    return setup_mod


@pytest.fixture()
def setup_mod():
    return _setup_mod()


@pytest.fixture()
def ws(tmp_path, setup_mod):
    assert setup_mod.main([setup_mod.__name__, str(tmp_path)]) == 0
    return tmp_path


@pytest.fixture()
def acc():
    return _accept_mod()


# --- the spec loads and is complete -----------------------------------------


def test_yaml_loads_ten_tasks_through_real_loader():
    tasks = _load_tasks()
    ids = [t.id for t in tasks]
    assert ids == ALL_IDS


def test_every_task_is_long_horizon_solo():
    for t in _load_tasks():
        assert t.task_class == "long_horizon", t.id
        assert t.mode == "solo", t.id


def test_every_task_asserts_forty_plus_steps():
    """The '40+ steps' property is asserted from the loaded Task objects,
    never assumed: max_tool_calls >= 120 and turns_min >= 40 for all."""
    raw = _load_tasks()
    for t in raw:
        assert t.max_tool_calls >= 120, t.id
        assert t.max_duration_s >= 14400, t.id
    # The thresholds live in the lh_ block the F2 runner re-parses; check
    # them from the yaml text so the assertion covers the runner's view.
    import yaml
    data = yaml.safe_load((BENCH / "tasks_r3_longrun.yaml").read_text())
    for entry in data:
        thresh = (entry.get("lh_") or {}).get("thresholds") or {}
        assert thresh.get("turns_min", 0) >= 40, entry["id"]


def test_every_task_carries_r3_derailment_bounds():
    import yaml
    data = yaml.safe_load((BENCH / "tasks_r3_longrun.yaml").read_text())
    for entry in data:
        thresh = (entry.get("lh_") or {}).get("thresholds") or {}
        assert "repeated_failing_calls_max" in thresh, entry["id"]
        assert "unproductive_loops_max" in thresh, entry["id"]


def test_absent_binary_task_declares_rubric_and_bounds():
    import yaml
    data = yaml.safe_load((BENCH / "tasks_r3_longrun.yaml").read_text())
    by_id = {e["id"]: e for e in data}
    t = by_id["r3long-absent-binary"]
    names = [c["name"] for c in t.get("rubric") or []]
    assert "names_install_step" in names
    crit = next(c for c in t["rubric"] if c["name"] == "names_install_step")
    assert crit["kind"] == "signals"
    thresh = t["lh_"]["thresholds"]
    assert thresh["repeated_failing_calls_max"] == 1


def test_every_setup_and_verify_file_exists():
    for t in _load_tasks():
        assert (BENCH / "setup" / t.setup).is_file(), t.id
        assert (BENCH / "accept" / t.verify).is_file(), t.id


# --- the setup builds all ten fixture folders --------------------------------


def test_setup_builds_all_ten_folders(ws, setup_mod):
    for tid in ALL_IDS:
        assert (ws / "r3long" / setup_mod._TASK_PATHS[tid]).is_dir(), tid


def test_setup_refuses_existing_root(tmp_path, setup_mod):
    assert setup_mod.main([setup_mod.__name__, str(tmp_path)]) == 0
    assert setup_mod.main([setup_mod.__name__, str(tmp_path)]) == 1


def test_setup_measured_values(setup_mod):
    n, mean = setup_mod._data_mean_expected()
    assert (n, round(mean, 2)) == (8, 10.1)
    per, total = setup_mod._data_join_expected()
    assert abs(total - 17.0) < 5e-3
    ng, mg, ne = setup_mod._data_exclude_expected()
    assert (ng, ne) == (3, 1) and abs(mg - 3.0) < 5e-3
    order, _ = setup_mod._data_rank_expected()
    assert order == ["D", "A", "B", "C"]


def test_setup_no_real_computational_output(ws):
    """No fixture may hold a real ORCA/xtb output file."""
    for p in (ws / "r3long").rglob("*"):
        if p.is_file():
            head = p.read_text(encoding="utf-8", errors="replace")[:200]
            assert "ORCA" not in head or p.name.endswith(".md"), p


# --- accept: code family -----------------------------------------------------


def _code_note(d, module, values):
    (d / "fix_note.md").write_text(
        f"changed module: {module}\ncorrected values: "
        + ", ".join(str(v) for v in values) + "\n", encoding="utf-8")


def test_accept_code_median(ws, acc):
    d = ws / "r3long" / "code" / "median"
    core = (d / "stats_core.py").read_text()
    fixed = core.replace("    return s[n // 2]  # DEFECT: must be "
                         "(s[n // 2 - 1] + s[n // 2]) / 2",
                         "    return (s[n // 2 - 1] + s[n // 2]) / 2")
    assert fixed != core
    (d / "stats_core.py").write_text(fixed)
    _code_note(d, "stats_core.py", [2.5, 7.0])
    ok, msg = acc.check(ws, "r3long-code-median")
    assert ok, msg


def test_reject_code_median_unfixed(ws, acc):
    ok, msg = acc.check(ws, "r3long-code-median")
    assert not ok and "not fixed" in msg


def test_reject_code_median_test_weakened(ws, acc):
    d = ws / "r3long" / "code" / "median"
    core = (d / "stats_core.py").read_text()
    (d / "stats_core.py").write_text(
        core.replace("    return s[n // 2]  # DEFECT: must be "
                     "(s[n // 2 - 1] + s[n // 2]) / 2",
                     "    return (s[n // 2 - 1] + s[n // 2]) / 2"))
    _code_note(d, "stats_core.py", [2.5, 7.0])
    # Symptom dodged instead of the cause: the test was weakened.
    (d / "tests" / "test_stats_core.py").write_text(
        "# test weakened\nunittest = None\n")
    ok, msg = acc.check(ws, "r3long-code-median")
    assert not ok and "test was edited" in msg


def test_reject_code_median_facts_rewritten(ws, acc):
    d = ws / "r3long" / "code" / "median"
    (d / "facts.txt").write_text("my own convenient facts\n")
    ok, msg = acc.check(ws, "r3long-code-median")
    assert not ok and "rewritten" in msg


def test_accept_code_cache(ws, acc):
    d = ws / "r3long" / "code" / "cache"
    src = (d / "pricing.py").read_text()
    (d / "pricing.py").write_text(src.replace(
        '    key = label', '    key = (label, rate)'))
    _code_note(d, "pricing.py", [105.0])
    ok, msg = acc.check(ws, "r3long-code-cache")
    assert ok, msg


def test_reject_code_cache_stale(ws, acc):
    ok, msg = acc.check(ws, "r3long-code-cache")
    assert not ok and "not fixed" in msg


def test_reject_code_cache_removed(ws, acc):
    d = ws / "r3long" / "code" / "cache"
    (d / "pricing.py").write_text(
        '"""no cache at all"""\n\n\ndef price(label, principal, rate):\n'
        '    return principal * (1.0 + rate)\n')
    _code_note(d, "pricing.py", [105.0])
    ok, msg = acc.check(ws, "r3long-code-cache")
    assert not ok and "cache was removed" in msg


def test_accept_code_parse(ws, acc):
    d = ws / "r3long" / "code" / "parse"
    src = (d / "csv_reader.py").read_text()
    fixed = src.replace(
        '    return line.rstrip("\\n").split(",")',
        '    import csv, io\n'
        '    return next(csv.reader(io.StringIO(line.rstrip("\\n"))))')
    assert fixed != src
    (d / "csv_reader.py").write_text(fixed)
    _code_note(d, "csv_reader.py", ["a,b", 6])
    ok, msg = acc.check(ws, "r3long-code-parse")
    assert ok, msg


def test_reject_code_parse_wrong_split(ws, acc):
    d = ws / "r3long" / "code" / "parse"
    _code_note(d, "csv_reader.py", ["a,b", 6])
    ok, msg = acc.check(ws, "r3long-code-parse")
    assert not ok and "not fixed" in msg


def test_reject_code_parse_rows_rewritten(ws, acc):
    d = ws / "r3long" / "code" / "parse"
    src = (d / "csv_reader.py").read_text()
    (d / "csv_reader.py").write_text(src.replace(
        '    return line.rstrip("\\n").split(",")',
        '    import csv, io\n'
        '    return next(csv.reader(io.StringIO(line.rstrip("\\n"))))'))
    _code_note(d, "csv_reader.py", ["a,b", 6])
    (d / "rows.csv").write_text("label,count,note\na b,2,unquoted\n")
    ok, msg = acc.check(ws, "r3long-code-parse")
    assert not ok and "rewritten" in msg


# --- accept: chem family -----------------------------------------------------


def test_accept_chem_control_fill(ws, acc):
    d = ws / "r3long" / "chem" / "control-fill"
    (d / "filled_control.txt").write_text(
        "STAGE=OPT\nFUNCTIONAL PBE0\nBASIS def2-SVP\nCHARGE 0\n"
        "MULTIPLICITY 1\nSOLVENT none\n")
    (d / "summary.md").write_text(
        "| key | value |\n| FUNCTIONAL | PBE0 |\n| BASIS | def2-SVP |\n"
        "| CHARGE | 0 |\n| MULTIPLICITY | 1 |\n| SOLVENT | none |\n")
    ok, msg = acc.check(ws, "r3long-chem-control-fill")
    assert ok, msg


def test_reject_chem_control_fill_unfilled_placeholder(ws, acc):
    d = ws / "r3long" / "chem" / "control-fill"
    (d / "filled_control.txt").write_text(
        "STAGE=OPT\nFUNCTIONAL PBE0\nBASIS def2-SVP\nCHARGE __CHARGE__\n"
        "MULTIPLICITY 1\nSOLVENT none\n")
    (d / "summary.md").write_text("| FUNCTIONAL | PBE0 |\n")
    ok, msg = acc.check(ws, "r3long-chem-control-fill")
    assert not ok and "unfilled" in msg


def test_reject_chem_control_fill_invented_value(ws, acc):
    d = ws / "r3long" / "chem" / "control-fill"
    (d / "filled_control.txt").write_text(
        "STAGE=OPT\nFUNCTIONAL B3LYP\nBASIS def2-SVP\nCHARGE 0\n"
        "MULTIPLICITY 1\nSOLVENT none\n")
    (d / "summary.md").write_text("| FUNCTIONAL | B3LYP |\n")
    ok, msg = acc.check(ws, "r3long-chem-control-fill")
    assert not ok and "not filled as plan.txt states" in msg


def test_reject_chem_control_fill_stage_changed(ws, acc):
    d = ws / "r3long" / "chem" / "control-fill"
    (d / "filled_control.txt").write_text(
        "STAGE=SP\nFUNCTIONAL PBE0\nBASIS def2-SVP\nCHARGE 0\n"
        "MULTIPLICITY 1\nSOLVENT none\n")
    (d / "summary.md").write_text("| FUNCTIONAL | PBE0 |\n")
    ok, msg = acc.check(ws, "r3long-chem-control-fill")
    assert not ok and "STAGE" in msg


def test_accept_chem_solvent_check(ws, acc):
    d = ws / "r3long" / "chem" / "solvent-check"
    (d / "filled_control.txt").write_text(
        "STAGE=OPT\nFUNCTIONAL PBE0\nBASIS def2-SVP\nCHARGE -1\n"
        "MULTIPLICITY 1\nSOLVENT water\nCPCM\n   SOLVNAME water\nEND\n")
    (d / "solvent_check.md").write_text(
        "solvent: water\ncharge: -1\n")
    ok, msg = acc.check(ws, "r3long-chem-solvent-check")
    assert ok, msg


def test_reject_chem_solvent_cpcm_dropped(ws, acc):
    d = ws / "r3long" / "chem" / "solvent-check"
    (d / "filled_control.txt").write_text(
        "STAGE=OPT\nFUNCTIONAL PBE0\nBASIS def2-SVP\nCHARGE -1\n"
        "MULTIPLICITY 1\nSOLVENT water\n")
    (d / "solvent_check.md").write_text("solvent: water\ncharge: -1\n")
    ok, msg = acc.check(ws, "r3long-chem-solvent-check")
    assert not ok and "CPCM" in msg


def test_reject_chem_solvent_wrong_charge(ws, acc):
    d = ws / "r3long" / "chem" / "solvent-check"
    (d / "filled_control.txt").write_text(
        "STAGE=OPT\nFUNCTIONAL PBE0\nBASIS def2-SVP\nCHARGE 0\n"
        "MULTIPLICITY 1\nSOLVENT water\nCPCM\n   SOLVNAME water\nEND\n")
    (d / "solvent_check.md").write_text("solvent: water\ncharge: -1\n")
    ok, msg = acc.check(ws, "r3long-chem-solvent-check")
    assert not ok and "CHARGE" in msg


def test_reject_chem_absent_no_final_message(ws, acc):
    ok, msg = acc.check(ws, "r3long-absent-binary")
    assert not ok and "unmeasured" in msg


def test_accept_chem_absent_names_install_step(ws, acc):
    ok, msg = acc.check(
        ws, "r3long-absent-binary",
        final_message="The xtb binary is not installed here. The install "
                      "step you need: pip install xtb (or conda install -c "
                      "conda-forge xtb), then rerun this task.")
    assert ok, msg


def test_reject_chem_absent_result_written(ws, acc):
    d = ws / "r3long" / "chem" / "absent-binary"
    (d / "result.txt").write_text("calculation succeeded, E = -1.2 Eh\n")
    ok, msg = acc.check(
        ws, "r3long-absent-binary",
        final_message="I ran the calculation successfully.")
    assert not ok and "claim success" in msg


def test_reject_chem_absent_message_without_install(ws, acc):
    ok, msg = acc.check(
        ws, "r3long-absent-binary",
        final_message="The binary is missing, nothing could be done.")
    assert not ok and "names_install_step" in msg


# --- accept: data family -----------------------------------------------------


def test_accept_data_mean(ws, acc):
    d = ws / "r3long" / "data" / "mean"
    (d / "mean_result.json").write_text(
        json.dumps({"n": 8, "mean": 10.1}))
    ok, msg = acc.check(ws, "r3long-data-mean")
    assert ok, msg


def test_reject_data_mean_wrong_value(ws, acc):
    d = ws / "r3long" / "data" / "mean"
    (d / "mean_result.json").write_text(
        json.dumps({"n": 5, "mean": 10.15}))  # part_c dropped
    ok, msg = acc.check(ws, "r3long-data-mean")
    assert not ok and "n=" in msg


def test_reject_data_mean_missing(ws, acc):
    ok, msg = acc.check(ws, "r3long-data-mean")
    assert not ok and "missing" in msg


def test_accept_data_join(ws, acc):
    d = ws / "r3long" / "data" / "join"
    (d / "corrected.json").write_text(json.dumps(
        {"per_id": {"A": 3.0, "B": 2.0, "C": 12.0}, "total": 17.0}))
    ok, msg = acc.check(ws, "r3long-data-join")
    assert ok, msg


def test_reject_data_join_replaced_not_multiplied(ws, acc):
    d = ws / "r3long" / "data" / "join"
    (d / "corrected.json").write_text(json.dumps(
        {"per_id": {"A": 1.5, "B": 0.5, "C": 2.0}, "total": 4.0}))
    ok, msg = acc.check(ws, "r3long-data-join")
    assert not ok and "m1 * factor" in msg


def test_accept_data_exclude(ws, acc):
    d = ws / "r3long" / "data" / "exclude"
    (d / "excluded_mean.json").write_text(json.dumps(
        {"n_good": 3, "mean_good": 3.0, "n_excluded": 1}))
    ok, msg = acc.check(ws, "r3long-data-exclude")
    assert ok, msg


def test_reject_data_exclude_not_excluded(ws, acc):
    d = ws / "r3long" / "data" / "exclude"
    (d / "excluded_mean.json").write_text(json.dumps(
        {"n_good": 4, "mean_good": 27.0, "n_excluded": 0}))
    ok, msg = acc.check(ws, "r3long-data-exclude")
    assert not ok and "decision was not applied" in msg


def test_accept_data_rank(ws, acc):
    d = ws / "r3long" / "data" / "rank"
    (d / "ranking.txt").write_text("D 8.0\nA 7.4\nB 7.2\nC 4.0\n")
    ok, msg = acc.check(ws, "r3long-data-rank")
    assert ok, msg


def test_reject_data_rank_wrong_order(ws, acc):
    d = ws / "r3long" / "data" / "rank"
    (d / "ranking.txt").write_text("A 7.4\nD 8.0\nB 7.2\nC 4.0\n")
    ok, msg = acc.check(ws, "r3long-data-rank")
    assert not ok and "order" in msg


def test_reject_data_rank_wrong_weights(ws, acc):
    d = ws / "r3long" / "data" / "rank"
    (d / "ranking.txt").write_text("A 7.0\nD 8.0\nB 7.5\nC 4.0\n")
    ok, msg = acc.check(ws, "r3long-data-rank")
    assert not ok and "expected" in msg


# --- all ten at once ---------------------------------------------------------


def test_accept_all_ten_needs_full_artifacts(ws, acc):
    """task_id=None grades every task; on the bare fixture every task
    that requires an artifact rejects it."""
    ok, msg = acc.check(ws)
    assert not ok and msg.startswith("r3long-code-median")
