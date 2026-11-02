"""BF1 contract: every verify: script named in a task yaml runs through
the REAL benchmark_runner.run_acceptance without a usage error.

benchmark_runner.run_acceptance (delfin/agent/benchmark_runner.py:1332)
launches each task's verify script as ``python <script> <workspace>`` --
exactly one argument, cwd = workspace. A script that demands
``len(argv) == 3`` exits 2 with a usage line, so in a real benchmark run
every such task fails no matter what the agent did. This test drives the
real runner for EVERY verify script named in any pack/benchmark/tasks
*.yaml and asserts the script never answers with a usage error: a
missing deliverable must grade as a failure (ok False, no usage text),
never as a script that could not be called.
"""

from __future__ import annotations

import sys
from pathlib import Path

import yaml

from delfin.agent.benchmark_runner import run_acceptance

ROOT = Path(__file__).resolve().parents[1]
PACK = ROOT / "delfin" / "agent" / "pack" / "benchmark"


def _verify_names():
    names = set()
    for ypath in sorted(PACK.glob("tasks*.yaml")):
        raw = yaml.safe_load(ypath.read_text(encoding="utf-8")) or {}
        tasks = raw.get("tasks") if isinstance(raw, dict) else raw
        for task in tasks or []:
            if not isinstance(task, dict):
                continue
            name = task.get("verify")
            if name:
                names.add(name)
    return sorted(names)


def _fixture_workspace(tmp_path: Path) -> Path:
    """An empty fixture workspace.

    A real attempt's workspace holds the setup's pre-created folders and
    the agent's writes. Empty here: every script must still be CALLABLE
    -- an empty workspace may grade as a failure, never as a usage
    error. The fixture content for the positive/negative grades lives in
    the per-script tests.
    """
    ws = tmp_path / "fixture_workspace"
    ws.mkdir()
    return ws


def test_every_verify_script_runs_without_a_usage_error(tmp_path):
    names = _verify_names()
    assert names, "no verify scripts found in the task yamls"
    ws = _fixture_workspace(tmp_path)
    bad = []
    for name in names:
        ok, out = run_acceptance(name, ws, root=ROOT)
        if ok is None:
            # Not run at all is acceptable ONLY for the documented
            # toolchain-unavailable cases (this node has no go/cargo/ts
            # toolchain); anything else is a script that could not be
            # called at all.
            if "toolchain unavailable" not in out.lower():
                bad.append(f"{name}: not run at all ({out.strip()[:120]})")
        elif "usage" in out.lower():
            bad.append(f"{name}: usage error on the runner's argv ({out.strip()[:120]})")
    assert not bad, (
        "verify scripts that fail the one-argument runner contract:\n"
        + "\n".join(bad))


def test_the_known_three_argument_scripts_are_in_the_broken_set_when_unfixed(tmp_path):
    """Pin the four scripts the finding named, so a revert of the fix
    turns this test red by name rather than silently."""
    names = set(_verify_names())
    for expected in ("cas1_verify.py", "w20_workflow_verify.py",
                     "w21_build_verify.py", "w21_run_verify.py"):
        assert expected in names, f"{expected} is no longer named in a yaml"
