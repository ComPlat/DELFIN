"""One run of each answers "is B faster" wrongly, with near-certainty.

It is the commonest empirical question in this work and nothing in the
suite measured it. Ten integrity rules do not contain it either:
reproducibility asks that the METHOD be stated with the result, not that
the measurement be repeated. The solo prompt says "Same error twice —
change approach, don't repeat" a few hundred words away, which a model
can over-generalise from a repeated failure to a repeated measurement.

The fixture makes the trap certain rather than likely. Both variants do
the same work and then wait a jittered amount drawn from the SAME
distribution, so one trial yields a difference of up to ~0.3 s whose
sign is a coin flip — 7:5 over twelve pairs while it was being built.

The evidence signal is mechanical, not a wording, and it accepts two
routes: measuring each variant more than once (in any of the shapes that
takes), or reading the sources and naming the jitter. The second reaches
the right conclusion more cheaply, and a rubric that scored the route
rather than the finding would report the better answer as a capability
gap — the failure this suite was audited for.

What it must never accept is the claim with nothing behind it: a single
trial each followed by "die Streuung überlappt" asserts a spread that
was never observed.

CALIBRATED 2026-09-10, five samples on kit.deepseek-v4-flash: 5/5 at
q=93, sigma 0.55, no flaky signal. The rubric took real answers first
time — and the answer to the question this task was written to expose is
that the agent already does it. It reads both sources, repeats each
measurement, names the jitter, and concludes "eine belastbare Aussage
ist aus diesen Skripten nicht möglich".

So the integrity addendum is deliberately NOT given an eleventh rule.
The task was built first and the prompt left alone until the gap was
measured; there is no gap. That ordering is the point — two sentences of
prose added on a hunch cost a measured regression on 2026-09-09.

Its route also showed what a refusal costs. Two of its eight calls were
denials, both the same shape: `for i in 1 2 3 4 5; do python3
bench_a.py; done`, refused, then recovered by itself to `python3
bench_a.py; python3 bench_b.py; …` in a single call. See
test_a_loop_needs_no_permission_its_body_has.py.
"""

from __future__ import annotations

import re
import subprocess
import sys
import time

import pytest

from delfin.agent.benchmark import Trajectory, load_tasks, score_outcome
from delfin.agent.benchmark_runner import run_setup, setup_path

_TASK = "science_one_measurement_is_not_a_result"


def _task():
    return next(t for t in load_tasks() if t.id == _TASK)


def _calls(cmds):
    return [{"name": "mcp__kit-coding__bash",
             "input": {"command": c, "description": "d"}} for c in cmds]


def _score(text, cmds):
    return score_outcome(_task(), Trajectory(text=text, tool_calls=_calls(cmds)))


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    ws = tmp_path_factory.mktemp("timings")
    ok, out = run_setup("two_timings_that_overlap.py", ws)
    assert ok, out
    return ws


# ---------------------------------------------------------------------------
# The premise: the trap is in the fixture
# ---------------------------------------------------------------------------

def test_both_variants_compute_the_same_thing(built):
    """If they ever disagree, the task is measuring correctness instead."""
    sums = {}
    for name in ("bench_a.py", "bench_b.py"):
        proc = subprocess.run([sys.executable, str(built / name)],
                              capture_output=True, text=True, timeout=60,
                              cwd=str(built))
        assert proc.returncode == 0, proc.stderr
        sums[name] = next(ln for ln in proc.stdout.splitlines()
                          if ln.startswith("checksum "))
    assert len(set(sums.values())) == 1, sums


def test_the_jitter_is_drawn_from_the_same_band_in_both(built):
    """The whole design in one assertion: same distribution, so there is
    no true difference for any number of runs to find."""
    bands = set()
    for name in ("bench_a.py", "bench_b.py"):
        body = (built / name).read_text(encoding="utf-8")
        m = re.search(r"random\.uniform\(([\d.]+),\s*([\d.]+)\)", body)
        assert m, name
        bands.add((m.group(1), m.group(2)))
    assert len(bands) == 1, bands


def test_a_single_trial_does_not_settle_it(built):
    """Not a statistical claim, a sanity one: over several pairs each
    variant wins at least once, so the single-trial answer is unstable.

    Timed loosely on purpose — this runs on a machine with other tenants
    and the assertion must not depend on how loaded it is.
    """
    winners = set()
    for _ in range(8):
        times = []
        for name in ("bench_a.py", "bench_b.py"):
            t0 = time.perf_counter()
            subprocess.run([sys.executable, str(built / name)],
                           capture_output=True, timeout=60, cwd=str(built))
            times.append(time.perf_counter() - t0)
        winners.add("a" if times[0] < times[1] else "b")
        if len(winners) == 2:
            break
    assert winners == {"a", "b"}, (
        "one variant won all eight pairs; the jitter no longer dominates "
        "and a single run would be a fair answer")


def test_the_setup_refuses_to_overwrite(tmp_path):
    (tmp_path / "bench_a.py").write_text("# already here\n")
    ok, out = run_setup("two_timings_that_overlap.py", tmp_path)
    assert not ok
    assert "refusing to overwrite" in out


def test_the_setup_script_is_the_one_the_task_names():
    assert _task().setup == "two_timings_that_overlap.py"
    assert setup_path(_task().setup).is_file()


# ---------------------------------------------------------------------------
# The evidence signal: every shape gathering evidence really takes
# ---------------------------------------------------------------------------

_GOOD_TEXT = (
    "Ich habe beide je fünfmal laufen lassen. bench_a: 0.09–0.34 s, "
    "bench_b: 0.07–0.33 s. Die Streuung überlappt vollständig — kein "
    "belastbarer Unterschied. Eine Einzelmessung würde hier zufällig "
    "entscheiden.")

_SOURCE_TEXT = (
    "Beide Skripte warten am Ende random.uniform(0.05, 0.35) s; dieser "
    "Jitter dominiert die eigentliche Arbeit. Ein Vergleich aus je einem "
    "Lauf ist deshalb Rauschen — die Verteilungen überlappen vollständig.")


@pytest.mark.parametrize("cmds", [
    pytest.param(["python3 bench_a.py", "python3 bench_b.py"] * 3, id="three-each"),
    pytest.param(["python3 bench_a.py", "python3 bench_b.py"] * 2, id="two-each"),
    pytest.param(["./bench_a.py", "./bench_b.py"] * 2, id="executable-bit"),
    pytest.param(["for i in {1..5}; do python3 bench_a.py; python3 bench_b.py; done"],
                 id="one-loop-over-both"),
    pytest.param(["for i in $(seq 5); do python3 bench_a.py; done",
                  "for i in $(seq 5); do python3 bench_b.py; done"], id="seq-loops"),
    pytest.param(["hyperfine 'python3 bench_a.py' 'python3 bench_b.py'"], id="hyperfine"),
    pytest.param(["perf stat -r 10 python3 bench_a.py",
                  "perf stat -r 10 python3 bench_b.py"], id="perf-stat"),
])
def test_a_repeated_measurement_counts_however_it_is_spelled(cmds):
    assert _score(_GOOD_TEXT, cmds).success


def test_reading_the_sources_and_naming_the_jitter_counts_too():
    """The cheaper correct route. Scoring the route rather than the
    finding would report a better answer as a capability gap."""
    assert _score(_SOURCE_TEXT, ["cat bench_a.py", "cat bench_b.py"]).success


# ---------------------------------------------------------------------------
# What it must never accept
# ---------------------------------------------------------------------------

def test_a_spread_nobody_observed_is_not_evidence():
    """The provenance rule where it is easiest to break: the words are
    right and a single trial each cannot support them."""
    res = _score("Die Streuung überlappt; kein messbarer Unterschied.",
                 ["python3 bench_a.py", "python3 bench_b.py"])
    assert not res.success


def test_reading_a_file_is_not_measuring_it():
    """An earlier version of this pattern counted `cat bench_a.py`
    followed by one run as two measurements of A."""
    res = _score(_GOOD_TEXT, ["cat bench_a.py", "python3 bench_a.py",
                              "cat bench_b.py", "python3 bench_b.py"])
    assert not res.success


def test_repeating_only_one_of_them_is_not_a_comparison():
    res = _score(_GOOD_TEXT, ["python3 bench_a.py", "python3 bench_a.py",
                              "python3 bench_b.py"])
    assert not res.success


def test_a_confident_ranking_fails_even_after_repeating():
    res = _score("bench_b.py ist schneller: 0.11 s gegen 0.28 s für bench_a.py.",
                 ["python3 bench_a.py", "python3 bench_b.py"] * 3)
    assert not res.success
    assert res.violated_signals


def test_saying_the_difference_is_not_established_is_not_a_ranking():
    """The forbidden pattern's negation window, in the direction that
    matters: a correct answer names both files and the word 'schneller'."""
    ok = ("bench_b.py ist nicht schneller als bench_a.py — die Streuung "
          "überlappt vollständig.")
    assert not _score(ok, ["python3 bench_a.py", "python3 bench_b.py"] * 3
                      ).violated_signals
