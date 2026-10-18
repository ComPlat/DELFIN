"""Gate tests for the w17_longrun package (KG4 phase 3).

Covers, all through the real load path and the real accept check()
functions:
* tasks_w17_longrun.yaml loads through the benchmark Task loader; each
  task's setup/verify names resolve to existing fixture files.
* the setup script builds all four fixture folders and the graded
  values are the measured ones (mean re-check).
* each accept script: accepts the right artifact and rejects each wrong
  variant (missing artifact, wrong value, unattested rubric fact).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

BENCH = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack" / "benchmark"


def _load_tasks():
    from delfin.agent.benchmark import load_tasks
    return load_tasks(BENCH / "tasks_w17_longrun.yaml")


def _accept(name):
    import importlib.util
    path = BENCH / "accept" / name
    spec = importlib.util.spec_from_file_location(name[:-3], path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --- the spec loads and is complete ----------------------------------------


def test_yaml_loads_four_tasks_through_real_loader():
    tasks = _load_tasks()
    ids = [t.id for t in tasks]
    assert ids == ["lr_wait_for_job", "lr_steered_goal",
                   "lr_forced_compaction", "lr_restart_midway"]


def test_every_setup_and_verify_file_exists():
    tasks = _load_tasks()
    for t in tasks:
        assert (BENCH / "setup" / t.setup).is_file(), t.id
        assert (BENCH / "accept" / t.verify).is_file(), t.id


# --- the setup builds all four folders with measured values ----------------


@pytest.fixture()
def fixture_ws(tmp_path):
    sys.path.insert(0, str(BENCH / "setup"))
    try:
        import w17_longrun_setup as setup_mod
    finally:
        sys.path.remove(str(BENCH / "setup"))
    rc = setup_mod.main([str(setup_mod.__name__), str(tmp_path)])
    assert rc == 0
    return tmp_path / "w17_longrun"


def test_setup_builds_all_four_folders_with_measured_mean(fixture_ws):
    for d in ("lr_wait_for_job", "lr_steered_goal", "lr_forced_compaction",
              "lr_restart_midway"):
        assert (fixture_ws / d).is_dir(), d
    rows = [ln.split(",") for ln in
            (fixture_ws / "lr_wait_for_job" / "samples.csv")
            .read_text().strip().splitlines()[1:]]
    mean = sum(float(r[1]) for r in rows) / len(rows)
    assert abs(mean - 42.19) <= 5e-3


def test_setup_refuses_existing_root(tmp_path):
    sys.path.insert(0, str(BENCH / "setup"))
    try:
        import w17_longrun_setup as setup_mod
    finally:
        sys.path.remove(str(BENCH / "setup"))
    assert setup_mod.main([setup_mod.__name__, str(tmp_path)]) == 0
    assert setup_mod.main([setup_mod.__name__, str(tmp_path)]) == 1


# --- accept: wait ----------------------------------------------------------


def _wait_ws(tmp_path):
    d = tmp_path / "w17_longrun" / "lr_wait_for_job"
    d.mkdir(parents=True)
    return d


def test_wait_accepts_correct_result_with_attested_waiting(tmp_path):
    d = _wait_ws(tmp_path)
    (d / "result.json").write_text(
        json.dumps({"n": 6, "mean_kj_per_mol": 42.19}))
    ok, msg = _accept("w17_longrun_wait_verify.py").check(
        tmp_path, {"idle_with_mail_s": 18.0})
    assert ok, msg


def test_wait_rejects_wrong_mean(tmp_path):
    d = _wait_ws(tmp_path)
    (d / "result.json").write_text(
        json.dumps({"n": 6, "mean_kj_per_mol": 99.0}))
    ok, msg = _accept("w17_longrun_wait_verify.py").check(
        tmp_path, {"idle_with_mail_s": 18.0})
    assert not ok and "expected" in msg


def test_wait_rejects_unattested_waiting(tmp_path):
    d = _wait_ws(tmp_path)
    (d / "result.json").write_text(
        json.dumps({"n": 6, "mean_kj_per_mol": 42.19}))
    ok, msg = _accept("w17_longrun_wait_verify.py").check(tmp_path, None)
    assert not ok and "unmeasured" in msg


def test_wait_rejects_missing_result(tmp_path):
    ok, msg = _accept("w17_longrun_wait_verify.py").check(
        tmp_path, {"idle_with_mail_s": 18.0})
    assert not ok and "missing" in msg


# --- accept: steer ---------------------------------------------------------


def _steer_ws(tmp_path):
    d = tmp_path / "w17_longrun" / "lr_steered_goal"
    d.mkdir(parents=True)
    (d / "STEERING_MESSAGE.txt").write_text(
        "PERSON (mid-task steering, new goal): now do X.\n")
    return d


def test_steer_accepts_new_goal_only(tmp_path):
    d = _steer_ws(tmp_path)
    (d / "answer.txt").write_text(
        "a: 3.4 nm\nb: absent\nc: 12.0 kJ/mol\nd: absent\n")
    ok, msg = _accept("w17_longrun_steer_verify.py").check(
        tmp_path, steer_followed=True)
    assert ok, msg


def test_steer_rejects_old_goal_still_done(tmp_path):
    d = _steer_ws(tmp_path)
    (d / "answer.txt").write_text(
        "a: 3.4 nm\nb: absent\nc: 12.0 kJ/mol\nd: absent\n")
    (d / "summary.txt").write_text("still did the old goal")
    ok, msg = _accept("w17_longrun_steer_verify.py").check(
        tmp_path, steer_followed=True)
    assert not ok and "OLD goal" in msg


def test_steer_rejects_unattested_follow(tmp_path):
    d = _steer_ws(tmp_path)
    (d / "answer.txt").write_text(
        "a: 3.4 nm\nb: absent\nc: 12.0 kJ/mol\nd: absent\n")
    ok, msg = _accept("w17_longrun_steer_verify.py").check(
        tmp_path, steer_followed=None)
    assert not ok and "not attested" in msg


def test_steer_rejects_wrong_values(tmp_path):
    d = _steer_ws(tmp_path)
    (d / "answer.txt").write_text("a: 9.9 nm\nb: absent\n")
    ok, msg = _accept("w17_longrun_steer_verify.py").check(
        tmp_path, steer_followed=True)
    assert not ok and "wrong entries" in msg


# --- accept: compaction ----------------------------------------------------


_NOTES_LINE = ("## Note {i}\n\nDecision {i}: the cutoff for batch {i} is "
               "{cut:.2f} eV. Rationale: derived from run {i} of the "
               "screening series (see batch_{i}.dat).\n")


def _compaction_ws(tmp_path):
    d = tmp_path / "w17_longrun" / "lr_forced_compaction"
    d.mkdir(parents=True)
    notes = "\n".join(_NOTES_LINE.format(i=i, cut=0.5 + 0.05 * (i % 7))
                      for i in range(1, 41))
    (d / "notes.md").write_text(notes)
    digest = "\n".join(
        f"Decision {i}: cutoff {0.5 + 0.05 * (i % 7):.2f} eV"
        for i in range(1, 41))
    (d / "digest.md").write_text(digest + "\n")
    return d


def test_compaction_accepts_complete_digest_with_attested_event(tmp_path):
    d = _compaction_ws(tmp_path)
    ok, msg = _accept("w17_longrun_compaction_verify.py").check(
        tmp_path, n_compactions=1)
    assert ok, msg


def test_compaction_rejects_missing_decisions(tmp_path):
    d = _compaction_ws(tmp_path)
    digest = (d / "digest.md").read_text().splitlines()
    (d / "digest.md").write_text("\n".join(digest[:-3]) + "\n")
    ok, msg = _accept("w17_longrun_compaction_verify.py").check(
        tmp_path, n_compactions=1)
    assert not ok and "37 decisions" in msg


def test_compaction_rejects_no_attested_event(tmp_path):
    _compaction_ws(tmp_path)
    ok, msg = _accept("w17_longrun_compaction_verify.py").check(
        tmp_path, n_compactions=0)
    assert not ok and "no compaction" in msg


def test_compaction_rejects_unmeasured_count(tmp_path):
    _compaction_ws(tmp_path)
    ok, msg = _accept("w17_longrun_compaction_verify.py").check(
        tmp_path, n_compactions=None)
    assert not ok and "not attested" in msg


# --- accept: restart -------------------------------------------------------


def _restart_ws(tmp_path):
    d = tmp_path / "w17_longrun" / "lr_restart_midway"
    d.mkdir(parents=True)
    return d


def test_restart_accepts_continuation_without_redo(tmp_path):
    d = _restart_ws(tmp_path)
    (d / "energy.txt").write_text(
        "alpha: 2.10 kJ/mol\nbeta: 5.80 kJ/mol\n"
        "gamma: 1.25 kJ/mol\ndelta: 3.10 kJ/mol\n")
    ok, msg = _accept("w17_longrun_restart_verify.py").check(
        tmp_path, redone_after_restart=False, restart_attested=True)
    assert ok, msg


def test_restart_rejects_redone_checkpoint(tmp_path):
    d = _restart_ws(tmp_path)
    (d / "energy.txt").write_text(
        "alpha: 2.10 kJ/mol\nbeta: 5.80 kJ/mol\n"
        "gamma: 1.25 kJ/mol\ndelta: 3.10 kJ/mol\n")
    ok, msg = _accept("w17_longrun_restart_verify.py").check(
        tmp_path, redone_after_restart=True, restart_attested=True)
    assert not ok and "redone" in msg


def test_restart_rejects_unattested_restart(tmp_path):
    d = _restart_ws(tmp_path)
    (d / "energy.txt").write_text(
        "alpha: 2.10 kJ/mol\nbeta: 5.80 kJ/mol\n"
        "gamma: 1.25 kJ/mol\ndelta: 3.10 kJ/mol\n")
    ok, msg = _accept("w17_longrun_restart_verify.py").check(
        tmp_path, redone_after_restart=False, restart_attested=None)
    assert not ok and "restart" in msg


def test_restart_rejects_wrong_value(tmp_path):
    d = _restart_ws(tmp_path)
    (d / "energy.txt").write_text(
        "alpha: 9.99 kJ/mol\nbeta: 5.80 kJ/mol\n"
        "gamma: 1.25 kJ/mol\ndelta: 3.10 kJ/mol\n")
    ok, msg = _accept("w17_longrun_restart_verify.py").check(
        tmp_path, redone_after_restart=False, restart_attested=True)
    assert not ok and "checkpoint not continued" in msg
