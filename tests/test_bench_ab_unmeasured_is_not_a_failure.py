"""The worker that scores one side must not read "nothing ran" as "it failed".

THE DEFECT, measured 2026-11-25 through `cli.main` on this branch with no
provider key in the environment:

    bench ab: origin/main (41fbf18e280e) -> origin/main (41fbf18e280e)
    task          base  head  rate   interval            verdict
    w27-sw1-01    0/2   0/2   +0.00  [-0.66, +0.66]     no measurable difference
    ...1 no measurable difference, 0 not measured
    EXIT 0

Every attempt died at engine init: no attempt put the question, and the report
filed two failures per side and called the pair indistinguishable. That reads
as a result about the code and is a result about the environment.

ROOT CAUSE, in `bench_ab._WORKER_SOURCE`:

    flags = list(result.per_run_success or [])
    ...
    if not flags:
        payload["unmeasured"] = max(1, int(result.n_samples or 1))

`benchmark.aggregate_replicates` marks a result `unmeasured=True` only when
EVERY replicate is unmeasured, and in that case KEEPS them with
`success=False` (benchmark.py:2399-2401). So `per_run_success` is [False, False]
-- non-empty -- the `if not flags` branch never fires, and the side is scored
0/2. The worker reads the wrong field: it asks "were there flags?" when the
question is "did any replicate measure anything?"

INPUT: the worker source is executed exactly as the subprocess runs it --
`python -P -c <source> task model provider backend repeats out` -- with
`run_task` replaced by a stub that returns an all-unmeasured BenchmarkResult.
No model is called, and nothing outside a tmpdir is written.

The control below is the same harness with a result that DID measure: it must
keep its flags and report unmeasured=0, so "always mark unmeasured" fails here.
"""
from __future__ import annotations

import json
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from delfin.agent import bench_ab as AB
from delfin.agent.benchmark import BenchmarkResult


def _run_worker(tmp_path: Path, *, all_unmeasured: bool) -> dict:
    """Execute bench_ab's worker with `run_task` stubbed; return its payload.

    The stub is installed by a sitecustomize-style prelude: the worker imports
    `delfin.agent.benchmark_runner` itself, so the stub must be in place when it
    does. `sys.modules` is pre-seeded by a wrapper source that runs BEFORE the
    worker's own text, which is exactly the seam the worker itself uses.
    """
    out = tmp_path / "attempt.json"
    if all_unmeasured:
        # What aggregate_replicates returns when every replicate was unmeasured:
        # the replicates are KEPT, with success False, and unmeasured True.
        result = {
            "task_id": "t-alpha", "task_class": "generalist", "model": "m",
            "success": False, "n_samples": 2, "unmeasured": True,
            "per_run_success": [False, False], "error":
                "engine init failed: Command 'pip install anthropic' returned 1",
        }
    else:
        result = {
            "task_id": "t-alpha", "task_class": "generalist", "model": "m",
            "success": True, "n_samples": 2, "unmeasured": False,
            "per_run_success": [True, True], "error": "",
        }

    prelude = textwrap.dedent(f"""
        import sys, types
        _stub = types.ModuleType("delfin.agent.benchmark_runner")

        class _R:
            pass

        def run_task(task, **kw):
            r = _R()
            for k, v in {result!r}.items():
                setattr(r, k, v)
            return r

        _stub.run_task = run_task
        sys.modules["delfin.agent.benchmark_runner"] = _stub
    """)

    # The worker source is used verbatim -- the same string the module runs.
    source = prelude + "\n" + AB._WORKER_SOURCE
    proc = subprocess.run(
        [sys.executable, "-P", "-c", source,
         "dash_nav_calc_typo", "m", "", "api", "2", str(out)],
        cwd=str(tmp_path), env={"PATH": "/usr/bin:/bin",
                                 "PYTHONPATH": str(Path(AB.__file__).resolve().parents[2])},
        capture_output=True, text=True, timeout=120)
    assert out.exists(), (
        "the worker wrote no result file; stderr was:\n"
        + (proc.stderr or proc.stdout or "")[-1500:]
    )
    return json.loads(out.read_text(encoding="utf-8"))


def test_an_all_unmeasured_run_is_reported_unmeasured_not_as_two_failures(tmp_path):
    """The regression this file exists for, at the exact line that had it.

    An all-unmeasured result must reach the payload as `unmeasured` and must
    NOT contribute two failing passes. What would make this fail: the payload
    saying passes=[False, False] with unmeasured=0 -- which is what the branch
    did before the fix, and what the real run printed as "0/2 ... no measurable
    difference".
    """
    payload = _run_worker(tmp_path, all_unmeasured=True)

    assert payload["passes"] == [], (
        "an attempt that never reached the model must not be counted as a "
        f"failed attempt; passes was {payload['passes']!r}"
    )
    assert payload["unmeasured"] == 2, (
        "both attempts must be reported unmeasured, got "
        f"{payload['unmeasured']!r} with passes {payload['passes']!r}"
    )


def test_a_run_that_measured_keeps_its_passes_and_is_not_unmeasured(tmp_path):
    """The control: a measuring run must keep its flags and not be unmeasured.

    Without this, "mark everything unmeasured" would pass the test above.
    """
    payload = _run_worker(tmp_path, all_unmeasured=False)

    assert payload["passes"] == [True, True], (
        f"a measuring run lost its passes: {payload['passes']!r}"
    )
    assert payload["unmeasured"] == 0, (
        f"a measuring run was marked unmeasured {payload['unmeasured']!r}"
    )
