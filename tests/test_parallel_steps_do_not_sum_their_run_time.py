"""A check that fails on correct work gets switched off.

The contract's wall-clock half asks whether the engine's stated run time
could have happened in this workspace. It used to SUM the stated time
over every ORCA step and compare the sum against the span the
workspace's files occupy.

DELFIN runs engine jobs in parallel by design. Four one-hour jobs started
together state four hours of run time and occupy a one-hour span, so the
sum made an honest run report ``wallclock_plausible=False``. A referee
running an ordinary redox workflow meets this on the first attempt, and a
flag that fires on correct work is a flag somebody switches off.

The fix takes the LONGEST step instead. Summing would require knowing
the steps ran one after another, and the files carry no schedule: the
only parallelism-safe statement is about one step at a time. That is
weaker on purpose -- four one-hour jobs run serially inside a one-hour
span now pass, and that is the price of never failing an honest run.

What still contradicts a claim, and is asserted below: ONE step whose
stated run time does not fit in the whole workspace's lifetime. No
schedule can rescue that.
"""
from __future__ import annotations

import os
import time

import pytest

from delfin.agent import run_evidence


def _out(path, hours, minutes=0, seconds=0.0):
    """An ORCA output stating a total run time."""
    path.write_text(
        "                             ****ORCA TERMINATED NORMALLY****\n"
        f"TOTAL RUN TIME: 0 days {hours} hours {minutes} minutes "
        f"{seconds} seconds\n",
        encoding="utf-8")
    return path


def _workspace(tmp_path, jobs, span_seconds):
    """A workspace whose files span *span_seconds*, with *jobs* outputs.

    ``jobs`` is a list of (name, hours). Every input is made older than
    its output so the first check (an output older than its input)
    cannot fire and confuse the measurement.
    """
    now = time.time()
    steps = []
    for i, (name, hours) in enumerate(jobs):
        inp = tmp_path / f"{name}.inp"
        inp.write_text("! PBE0 def2-SVP\n", encoding="utf-8")
        out = _out(tmp_path / f"{name}.out", hours)
        # inputs at the start of the span, outputs at its end
        os.utime(inp, (now - span_seconds, now - span_seconds))
        os.utime(out, (now, now))
        steps.append({"program": "orca", "input": inp.name,
                      "output": out.name, "step": name})
    return steps


def _flag(tmp_path, steps):
    ev = run_evidence.collect(tmp_path)
    flag, problems, _notes = run_evidence._wallclock(ev, steps)
    return flag, problems


# --- the defect this file exists for ----------------------------------

def test_four_parallel_hour_jobs_in_one_hour_are_plausible(tmp_path):
    """The honest run that the sum rejected. Four jobs state one hour
    each; they ran together, so the workspace spans one hour."""
    steps = _workspace(tmp_path,
                       [("j1", 1), ("j2", 1), ("j3", 1), ("j4", 1)],
                       span_seconds=3700)
    flag, problems = _flag(tmp_path, steps)
    assert flag is not False, f"an honest parallel run was rejected: {problems}"
    assert not problems


def test_one_step_longer_than_the_whole_span_still_contradicts(tmp_path):
    """What no schedule can rescue: a single job claiming more time than
    the workspace has existed."""
    steps = _workspace(tmp_path, [("slow", 5)], span_seconds=600)
    flag, problems = _flag(tmp_path, steps)
    assert flag is False
    assert problems and "slow.out" in problems[0]
    assert "one step's" in problems[0]


def test_the_problem_names_the_step_rather_than_a_total(tmp_path):
    """A message about a sum tells nobody which file to look at."""
    steps = _workspace(tmp_path, [("fast", 0), ("slow", 9)],
                       span_seconds=300)
    _flag_, problems = _flag(tmp_path, steps)
    assert problems
    assert "slow.out" in problems[0] and "fast.out" not in problems[0]


def test_many_short_steps_never_add_up_to_a_contradiction(tmp_path):
    """Ten six-minute jobs in a ten-minute span: the sum would be an
    hour and would fail. Each one fits, so the run passes."""
    steps = _workspace(tmp_path, [(f"j{i}", 0) for i in range(10)],
                       span_seconds=600)
    for i in range(10):
        _out(tmp_path / f"j{i}.out", 0, minutes=6)
    now = time.time()
    for i in range(10):
        os.utime(tmp_path / f"j{i}.out", (now, now))
    flag, problems = _flag(tmp_path, steps)
    assert flag is not False, problems


# --- the price of the fix, stated rather than hidden -------------------

def test_a_serial_run_that_overflows_is_no_longer_caught(tmp_path):
    """The deliberate weakening, asserted so it is documented.

    Four one-hour jobs run one after another cannot fit in one hour, and
    the maximum rule does not see it: catching that needs the schedule,
    which the workspace's files do not carry. Recorded here so the limit
    is known rather than discovered.
    """
    steps = _workspace(tmp_path,
                       [("s1", 1), ("s2", 1), ("s3", 1), ("s4", 1)],
                       span_seconds=3700)
    flag, problems = _flag(tmp_path, steps)
    assert flag is not False
    assert not problems


# --- the other half must keep working ---------------------------------

def test_an_output_older_than_its_input_is_still_a_contradiction(tmp_path):
    inp = tmp_path / "a.inp"
    inp.write_text("! PBE0\n", encoding="utf-8")
    out = _out(tmp_path / "a.out", 0, minutes=1)
    now = time.time()
    os.utime(out, (now - 600, now - 600))   # output BEFORE its input
    os.utime(inp, (now, now))
    flag, problems = _flag(
        tmp_path, [{"program": "orca", "input": "a.inp", "output": "a.out",
                    "step": "a"}])
    assert flag is False
    assert any("did not produce" in p for p in problems)


def test_nothing_to_judge_stays_none(tmp_path):
    """No ORCA output, no stated time: the flag is None, never True.
    A check that cannot run does not pass."""
    flag, _problems = _flag(tmp_path, [])
    assert flag is None
