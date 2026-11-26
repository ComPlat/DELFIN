"""A round that ran returns its record, even when it could not be kept.

``bench_nightly.run`` has one documented job: run a round and return the
:class:`RunRecord` it produced. ``main()`` is the layer that turns a lost
record into a non-zero exit code, and it says why in a comment: "exiting 0
would let a supervisor log a successful nightly that left nothing behind".

The history store added a second failure mode to ``run()`` itself -- a round
whose history root resolves inside a git work tree raises out of ``run()``
before returning its record. The CLI is loud about it (exit 4), which is
right. A CALLER of ``run()`` is not: it paid for the round, it holds no
record, and it cannot tell a lost record from a crash in the round. The
distinction the record is for -- "the round ran, its bytes are gone" -- is
exactly what raising destroys.

So the failure is recorded ON the record and reported by ``main()``. The
round's cost, its instances and its reason survive; the caller decides.
"""

from __future__ import annotations

import pytest


class _Sampler:
    def sample(self, seed, n):
        from delfin.agent.benchmark import Task

        return [Task(id=f"t{i}", task_class="generalist", mode="solo",
                     prompt="do the thing") for i in range(n)]


class _Grader:
    def grade(self, task, result):
        return True, 5


class _Accounting:
    def estimate(self, task):
        return 0.1

    def charge(self, core_hours):
        pass


class _Result:
    success = True
    n_samples = 1
    per_run_success = [True]
    error = ""
    unmeasured = False
    input_tokens = 0
    output_tokens = 0


def _runner(task, **kw):
    return _Result()


def _config(**extra):
    cfg = {"model": "M", "sampler": _Sampler(), "grader": _Grader(),
           "accounting": _Accounting(), "runner": _runner,
           "instances": 2, "mode": "practice", "public_seed": 7}
    cfg.update(extra)
    return cfg


def test_a_round_whose_record_cannot_be_written_still_returns_it(
        tmp_path, monkeypatch):
    """The defect: run() raised and the round it had just run was lost."""
    from delfin.agent import bench_history as bh
    from delfin.agent import bench_nightly as bn

    (tmp_path / ".git").mkdir()
    monkeypatch.setattr(bh, "default_root",
                        lambda: tmp_path / "checkout" / "rounds")

    record = bn.run(_config())

    assert record.instances, (
        "the round ran and its instances vanished with the exception")
    assert record.history_error, (
        "the record does not say why it could not be kept, so a caller "
        "cannot tell a lost record from a healthy round")
    assert "repository" in record.history_error


def test_a_healthy_round_records_no_history_error():
    """The control: the new field must stay empty when the write works."""
    from delfin.agent import bench_nightly as bn

    record = bn.run(_config(persist=False))
    assert not record.history_error
    assert record.instances


def test_the_cli_still_exits_non_zero_when_the_record_is_lost(
        tmp_path, monkeypatch, capsys):
    """``main()`` keeps its contract: a lost record is not a success.

    Without this, moving the failure onto the record would have made a lost
    round exit 0 -- the silent loss the store exists to prevent.
    """
    from delfin.agent import bench_history as bh
    from delfin.agent import bench_nightly as bn

    (tmp_path / ".git").mkdir()
    monkeypatch.setattr(bh, "default_root",
                        lambda: tmp_path / "checkout" / "rounds")

    from delfin.agent import bench_holdout_glue

    cfg = _config()
    monkeypatch.setattr(bench_holdout_glue, "build_default_config",
                        lambda **kw: dict(cfg))
    rc = bn.main(["--model", "M", "--instances", "2", "--practice"])
    assert rc == 4, rc
    err = capsys.readouterr().err
    assert "NOT kept" in err, err
