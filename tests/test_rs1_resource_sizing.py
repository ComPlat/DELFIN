"""RS1 phase-3 tests: docs/resource_sizing.md cites measured numbers.

Every number in the sizing advice must trace to a recorded run or a
labeled secondary source (wave-22 rule 3). The SV1 primary artifacts
were NOT in git (.gate/ excluded from every wave tree) — the operator
supplied the exact lines of both SV1 job summaries plus three
corrections to the earlier prose (run 1 FAILED exit 138, not a clean
walltime kill; core-walltime counts LOGICAL CPUs incl. hyperthreads;
run 2 NSOLV 3 vs run 1 NSOLV 12 — different tasks, so "pal 1 finished
the same task in 38 s" was WRONG). The exact lines are fixtures here
so the citations in docs/resource_sizing.md are checkable, not prose.

Red on the previous commit: docs/resource_sizing.md does not exist.
"""
import re
from pathlib import Path

import pytest

from delfin.agent import cluster_efficiency as ce

DOCS = Path(__file__).resolve().parents[1] / "docs" / "resource_sizing.md"

# --- Fixtures: exact lines from the SV1 run summaries -----------------
# Source: operator message (wave 22), md5 of run 1 .out
# 2f195b869915ee805aa9602ec31be629. Primary artifacts live under
# .gate/practice/ (not committed); these lines are the committed
# citation record until a fixture directory exists.

SV1_RUN1_LINES = (
    # .gate/practice/sv1_tiny_solvator/delfin_7922181.out (NSOLV 12)
    "CPUs: 2",
    "Time Limit: 00:30:00",
    "State: FAILED (exit code 138)",
    "CPU Utilized: 00:01:31",
    "CPU Efficiency: 1.53% of 01:39:08 core-walltime",
    "Job Wall-clock time: 00:24:47",
)

SV1_RUN2_LINES = (
    # .gate/practice/sv1_tiny_solvator_r2/delfin_7922985.out (NSOLV 3)
    "CPUs: 1",
    "Time Limit: 00:45:00",
    "CPU Utilized: 00:02:33",
    "CPU Efficiency: 35.75% of 00:07:08 core-walltime",
    "Job Wall-clock time: 00:03:34",
)

# --- Red checks on the doc ---------------------------------------------


def test_resource_sizing_doc_exists_and_cites_the_sv1_lines():
    text = DOCS.read_text(encoding="utf-8")
    for line in ("1.53%", "35.75%", "exit code 138", "NSOLV 12",
                 "NSOLV 3", "hyperthread"):
        assert line in text, f"doc must cite {line!r}"


def test_doc_states_the_core_walltime_logical_cpu_convention():
    text = DOCS.read_text(encoding="utf-8")
    # The operator correction: 24:47 x 2 logical CPUs x 2 threads =
    # 01:39:08 core-walltime; the percentage denominator is logical
    # CPUs, so the doc must say so, not imply physical cores.
    assert "logical" in text.lower()
    assert "01:39:08" in text and "00:24:47" in text


def test_doc_makes_no_same_task_claim_between_run1_and_run2():
    text = DOCS.read_text(encoding="utf-8")
    # The retracted claim "pal 1 finished the same task in 38 s" must
    # not appear; the two runs were different tasks (NSOLV 12 vs 3).
    assert not re.search(r"same task in 38\s*s", text)
    assert "NSOLV 12" in text and "NSOLV 3" in text


def test_threshold_justification_matches_the_cited_populations():
    # LOW_PARALLEL_EFFICIENCY = 0.10 must sit between the cited
    # populations (1.53 % pal 2 FAILED run, 35.75 % pal 1 run) and be
    # justified in the doc, far from either.
    assert 0.0153 < ce.LOW_PARALLEL_EFFICIENCY < 0.3575
    text = DOCS.read_text(encoding="utf-8")
    assert "LOW_PARALLEL_EFFICIENCY" in text
    assert "0.10" in text
