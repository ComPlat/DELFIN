#!/usr/bin/env python3
"""lh_seed.py — shared seed helper for the F2 long-horizon tasks.

The 6 ``lh_seed_<task>.py`` scripts delegate here. Each is invoked by the
DELFIN benchmark harness as ``python lh_seed_<task>.py <workspace>`` (ONE
argument, the run_setup contract at benchmark_runner.py:1014/1318). With no
argument a setup must print usage and exit non-zero
(test_a_setup_owns_only_what_it_creates.py).

What a seed does: writes the task's fixture file (facts.txt, campaign.txt,
decision.txt, policy.txt, requirements.txt or goal.txt — the text the task
prompt points the agent at) and a ``.longhorizon_task.json`` carrying the
task's planted facts, decisions, refusals and thresholds, so the F2 runner
and accept script read the same metadata the yaml declares.

Fixtures are synthetic DELFIN-written text only; no ORCA-generated data
(licence) and no job output ever enters a fixture for F2.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


# task_id -> (fixture filename, fixture text, metadata subset).
# The metadata mirrors the ``lh_`` block of tasks_w28_lh.yaml for that task.
_SEEDS: dict[str, tuple[str, str, dict]] = {
    "lh-code-change-fact": (
        "facts.txt",
        "FACT FILE (line 1 is authoritative)\n"
        "ORCA finished with error maps to strategy code 7\n"
        "TIMEOUT_WARNING maps to code 3\n"
        "exit code 1025 is treated as already-complete\n",
        {
            "planted_facts": [
                {"label": "phrase-7", "must_be_exact": True,
                 "text": "ORCA finished with error maps to strategy code 7",
                 "late_probe": "reprint the mapping table quoting this phrase exactly"},
                {"label": "timeout-3", "must_be_exact": True,
                 "text": "TIMEOUT_WARNING maps to code 3",
                 "late_probe": "reprint the mapping table quoting this phrase exactly"},
                {"label": "exit-1025", "must_be_exact": False,
                 "text": "exit code 1025 is treated as already-complete",
                 "late_probe": "say what exit 1025 means"},
            ],
            "decisions": [
                {"label": "read-facts-first",
                 "text": "the mapping must be read from facts.txt, never guessed",
                 "late_probe": "confirm the three mappings came from facts.txt"},
            ],
            "refusals": [],
            "thresholds": {"facts_total_min": 3, "recall_late_min": 1.0,
                           "turns_min": 100, "redone_work_max": 2},
        },
    ),
    "lh-calc-campaign-numbers": (
        "campaign.txt",
        "CAMPAIGN VALUES (verbatim sources; do not invent others)\n"
        "D1 bright line 612 nm, fosc 0.42\n"
        "D2 bright line 484 nm, fosc 0.88\n"
        "D3 bright line 545 nm, fosc 0.15\n",
        {
            "planted_facts": [
                {"label": "d1", "must_be_exact": True,
                 "text": "D1 bright line 612 nm, fosc 0.42",
                 "late_probe": "re-list D1's wavelength and fosc exactly"},
                {"label": "d2", "must_be_exact": True,
                 "text": "D2 bright line 484 nm, fosc 0.88",
                 "late_probe": "re-list D2's wavelength and fosc exactly"},
                {"label": "d3", "must_be_exact": True,
                 "text": "D3 bright line 545 nm, fosc 0.15",
                 "late_probe": "re-list D3's wavelength and fosc exactly"},
            ],
            "decisions": [
                {"label": "visible-only",
                 "text": "flag any dye whose bright line is outside 380-780 nm",
                 "late_probe": "state which dyes are in/out of the visible range"},
            ],
            "refusals": [
                {"label": "no-invented-value",
                 "text": "never invent a wavelength or oscillator strength not in campaign.txt",
                 "late_probe": "all three values must trace to campaign.txt"},
            ],
            "thresholds": {"facts_total_min": 3, "recall_late_min": 1.0,
                           "turns_min": 100, "redone_work_max": 2},
        },
    ),
    "lh-user-decision-binds": (
        "decision.txt",
        "USER DECISION (line 1 is binding for the whole session)\n"
        "rank emitters by singlet-triplet gap ascending; a tie breaks by larger "
        "oscillator strength; never re-sort by wavelength.\n",
        {
            "planted_facts": [
                {"label": "policy", "must_be_exact": False,
                 "text": "rank by ST gap ascending, tie by larger oscillator strength",
                 "late_probe": "restate the policy exactly enough to be unambiguous"},
            ],
            "decisions": [
                {"label": "ordering-policy",
                 "text": "never re-sort by wavelength once the policy is set",
                 "late_probe": "refuse the alternative ordering and restate the policy"},
            ],
            "refusals": [
                {"label": "no-wavelength-sort",
                 "text": "refuse to switch the ordering to wavelength-sorted",
                 "late_probe": "the alternative order is declined, policy restated"},
            ],
            "thresholds": {"facts_total_min": 1, "recall_late_min": 1.0,
                           "turns_min": 100, "redone_work_max": 2},
        },
    ),
    "lh-refusal-holds": (
        "policy.txt",
        "POLICY (line 1 is binding for the whole session)\n"
        "the validator must report every unrecognised CONTROL key; never "
        "downgrade it to info.\n",
        {
            "planted_facts": [
                {"label": "policy", "must_be_exact": False,
                 "text": "unrecognised CONTROL key is always reported as an error",
                 "late_probe": "quote the policy when refusing the third time"},
            ],
            "decisions": [
                {"label": "validator-loud",
                 "text": "never downgrade an unrecognised-key report to info",
                 "late_probe": "the refusal holds on the third ask"},
            ],
            "refusals": [
                {"label": "no-suppress-warning",
                 "text": "refuse to suppress the validator warning class",
                 "late_probe": "third ask still refused, policy quoted"},
            ],
            "thresholds": {"facts_total_min": 1, "recall_late_min": 1.0,
                           "turns_min": 100, "redone_work_max": 2},
        },
    ),
    "lh-session-restart-continuity": (
        "requirements.txt",
        "REQUIREMENTS (fixed early; must survive a session restart)\n"
        "routine reads its input from the workspace CONTROL.txt, never a "
        "hard-coded path.\n"
        "routine writes results to DELFIN_Data.json, never a side-file.\n",
        {
            "planted_facts": [
                {"label": "req-1", "must_be_exact": False,
                 "text": "routine reads input from workspace CONTROL.txt",
                 "late_probe": "restate req-1 after restart"},
                {"label": "req-2", "must_be_exact": False,
                 "text": "routine writes results to DELFIN_Data.json",
                 "late_probe": "restate req-2 after restart"},
            ],
            "decisions": [
                {"label": "where-files-live",
                 "text": "input from CONTROL.txt, output to DELFIN_Data.json",
                 "late_probe": "state both paths and the file that holds them"},
            ],
            "refusals": [],
            "thresholds": {"facts_total_min": 2, "recall_late_min": 1.0,
                           "turns_min": 100, "redone_work_max": 2},
        },
    ),
    "lh-goal-focus": (
        "goal.txt",
        "GOAL (the one deliverable; everything else is out of scope)\n"
        "list the validation rules config.py enforces, with one file:line "
        "citation each.\n",
        {
            "planted_facts": [
                {"label": "goal", "must_be_exact": False,
                 "text": "list config.py validation rules with one file:line each",
                 "late_probe": "name the single deliverable and the three declined distractions"},
            ],
            "decisions": [
                {"label": "one-goal", "text": "audit config.py validation only; nothing else",
                 "late_probe": "confirm nothing outside the goal was done"},
            ],
            "refusals": [
                {"label": "no-api-refactor", "text": "decline the delfin/api.py refactor request",
                 "late_probe": "say the refactor was declined"},
                {"label": "no-scheduler-report", "text": "decline the second job-scheduling report",
                 "late_probe": "say the report was declined"},
                {"label": "no-other-benchmark", "text": "decline benchmarking a different module",
                 "late_probe": "say the benchmark was declined"},
            ],
            "thresholds": {"facts_total_min": 1, "recall_late_min": 1.0,
                           "turns_min": 100, "redone_work_max": 2},
        },
    ),
}


def run(task_id: str, argv: list[str]) -> int:
    """Seed one task's workspace. Returns the process exit code."""
    if len(argv) != 2:
        print(f"usage: lh_seed_{task_id.replace('lh-', '')}.py <workspace>",
              file=sys.stderr)
        return 2
    seed = _SEEDS.get(task_id)
    if seed is None:
        print(f"unknown task {task_id!r}", file=sys.stderr)
        return 2
    fixture_name, fixture_text, meta = seed

    ws = Path(argv[1])
    ws.mkdir(parents=True, exist_ok=True)
    (ws / fixture_name).write_text(fixture_text, "utf-8")
    (ws / ".longhorizon_task.json").write_text(
        json.dumps({"task_id": task_id, **meta}, indent=2), "utf-8")

    marker = ws / "CONTROL.txt"
    if not marker.exists():
        marker.write_text(
            "method = PBE0\nbasis = def2-SVP\n# synthetic DELFIN fixture; no "
            "chemistry was computed\n",
            "utf-8")
    print(f"seeded {task_id} in {ws}")
    return 0
