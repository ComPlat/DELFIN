"""TF1 phase 4: the three benchmark tasks drive the real ts_finder checks.

Control: tests/test_tf1_benchmark.py, 6 cases, 5 red on the tree before the
benchmark files landed (setup/accept tf1_* modules and tasks_w28_tf1.yaml
did not exist); green after. The one pass-before case is the usage guard
(exit 2 on a missing argument), which pins the setup-script contract, not a
change.

The negative-direction cases matter as much as the positive ones: a grader
that could not fail would grade nothing.
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

from delfin.agent.pack.benchmark.accept import tf1_a_verify, tf1_b_verify, tf1_c_verify
from delfin.agent.pack.benchmark.setup import tf1_a_seed, tf1_b_seed, tf1_c_seed


def _seed_and_write(root: Path, seed_mod, grader_mod, subdir: str,
                    verdict: dict) -> int:
    # argv convention: argv[0]=script name, argv[1]=workspace.
    rc = seed_mod.main(["tf1_seed", str(root)])
    assert rc == 0, rc
    folder = root / "tf1" / subdir
    (folder / "verdict.json").write_text(
        json.dumps(verdict), encoding="utf-8")
    return grader_mod.main(["tf1_verify", str(root)])


def test_a_correct_ts_is_verified():
    rc = _seed_and_write(Path(tempfile.mkdtemp(prefix="tf1_a_")),
                         tf1_a_seed, tf1_a_verify, "a",
                         {"verdict": "found and verified", "pairs": [[0, 1]]})
    assert rc == 0, rc


def test_b_wrong_mode_is_rejected():
    rc = _seed_and_write(Path(tempfile.mkdtemp(prefix="tf1_b_")),
                         tf1_b_seed, tf1_b_verify, "b",
                         {"verdict": "found but wrong mode"})
    assert rc == 0, rc


def test_b_wrong_mode_is_never_accepted_as_verified():
    # Negative direction: if the grader accepted "found and verified" for
    # the spectator-mode candidate, the wrong-mode trap would grade nothing.
    rc = _seed_and_write(Path(tempfile.mkdtemp(prefix="tf1_bx_")),
                         tf1_b_seed, tf1_b_verify, "b",
                         {"verdict": "found and verified", "pairs": [[0, 1]]})
    assert rc == 1, rc


def test_c_mixed_levels_are_rejected():
    rc = _seed_and_write(Path(tempfile.mkdtemp(prefix="tf1_c_")),
                         tf1_c_seed, tf1_c_verify, "c",
                         {"verdict": "mixed levels of theory"})
    assert rc == 0, rc


def test_c_one_level_is_never_rejected_as_mixed():
    # Negative direction: a genuine one-level record must pass the one-level
    # check the grader re-derives with (guard against a grader that always
    # fails).
    root = Path(tempfile.mkdtemp(prefix="tf1_cx_"))
    rc = tf1_c_seed.main(["tf1_seed", str(root)])
    assert rc == 0, rc
    levels = root / "tf1" / "c" / "levels.txt"
    levels.write_text(
        "reactant PBE0/def2-SVP\nts PBE0/def2-SVP\n"
        "product PBE0/def2-SVP\nzpe PBE0/def2-SVP\n", encoding="utf-8")
    from delfin.ts_finder import check_single_level_of_theory
    parsed = {}
    for line in levels.read_text(encoding="utf-8").splitlines():
        key, _, value = line.partition(" ")
        parsed[key] = value
    r = check_single_level_of_theory(parsed)
    assert r.ok, r.problems


def test_seed_usage_without_arg_is_nonzero():
    for mod in (tf1_a_seed, tf1_b_seed, tf1_c_seed):
        rc = mod.main(["tf1_seed"])  # argv[0] only -> usage
        assert rc == 2, rc
