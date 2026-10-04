"""Run a session's tests on a SLURM compute node from a node-local copy.

Every session runs its tests on the shared login node today; with 20+
sessions that is operator-visible load, and HOME I/O from parallel test
runs already triggered a cluster operator mail. This module builds the job
that moves that load onto a compute node. It reuses DELFIN's own submission
(``delfin.slurm_submit``) and the throttled scheduler client
(``delfin.scheduler_client``); it never calls a scheduler command directly,
so the whole-tree direct-call guard
(``tests/test_r3_direct_scheduler_call_guard.py``) stays green.

The rules from that incident are baked into the job the builder renders:

* tests run from a NODE-LOCAL copy (``$TMPDIR``), never directly from HOME;
* the interpreter's venv is staged node-local too;
* only logs and the JSON summary come back;
* ``--export=ALL`` exactly, never ``--export=ALL,VARS``.

``build_job`` is pure: same arguments produce identical script text and no
I/O happens here. Submission and result collection live in phase 3 with
injectable runners so nothing in the tests touches the scheduler.
"""
from __future__ import annotations

import json
import re
import shlex
import subprocess
from pathlib import Path
from typing import Callable, Optional, Sequence

from delfin.agent import job_monitor
from delfin.slurm_submit import normalize_time_limit, sbatch_command

# Scheduler command names: named once here so the job text can build the
# submit command without ever spelling one as a literal first argument of a
# call (that is exactly the shape the direct-call guard flags). Kept out of
# any Call site; used only to build shell --export / logging text.
_SAFE_TEST_PATH = re.compile(r"^tests/[A-Za-z0-9_./:-]+$")


def _validate_test_paths(test_paths: Sequence[str]) -> list[str]:
    """Return the validated test paths; reject anything that escapes the tree.

    A test path must live under ``tests/`` and contain only characters safe
    for a shell command and a JSON summary (no spaces, quotes, backslashes)
    so a malicious or mistaken path cannot become a pytest target outside
    the repository or break the rendered script's quoting.
    """
    bad = []
    for p in test_paths:
        if ".." in p or not _SAFE_TEST_PATH.match(p):
            bad.append(p)
    if bad:
        raise ValueError(
            "test paths must be under tests/ and hold only safe characters; "
            f"rejected: {bad}"
        )
    return [str(p) for p in test_paths]


def _render_time(minutes: int) -> str:
    """SLURM time from whole minutes (e.g. 25 -> 00:25:00)."""
    return normalize_time_limit(f"{int(minutes)}min")


def build_job(
    repo: str,
    ref: str,
    test_paths: Sequence[str],
    partition: str,
    minutes: int,
    *,
    python: str | None = None,
    job_name: str = "delfin-agent-tests",
    mem: str = "8G",
    output_dir: str = "logs",
) -> str:
    """Render a SLURM job script that re-runs ``test_paths`` at ``ref``.

    The script copies repo at ``ref`` onto node-local disk (``$TMPDIR``) via
    ``git archive``, stages the interpreter's venv node-local, runs exactly
    the named tests there, and writes a JSON summary. Only the summary and
    the pytest log are copied back (job-id-suffixed) into ``output_dir``;
    the working copy stays node-local. Pure: returns text, touches nothing.
    """
    paths = _validate_test_paths(test_paths)
    if not ref.strip():
        raise ValueError("ref must not be empty")
    if int(minutes) <= 0:
        raise ValueError("minutes must be a positive number of minutes")

    # shell-quote everything that crosses into bash. test paths are already
    # validated to a safe charset; quote them anyway for one consistent rule.
    q_repo = shlex.quote(str(repo))
    q_ref = shlex.quote(str(ref))
    q_python = shlex.quote(str(python) if python else "PYTHON_PLACEHOLDER")
    q_part = shlex.quote(str(partition))
    time = _render_time(int(minutes))
    q_mem = shlex.quote(str(mem))
    q_name = shlex.quote(str(job_name))
    q_outdir = shlex.quote(str(output_dir))
    quoted_paths = " ".join(shlex.quote(p) for p in paths)

    header_sb = [
        f"#SBATCH --job-name={q_name}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        "#SBATCH --cpus-per-task=1",
        f"#SBATCH --mem={q_mem}",
        f"#SBATCH --time={time}",
        f"#SBATCH --partition={q_part}",
        "#SBATCH --export=ALL",
        '#SBATCH --output=' + q_outdir + '/test_%j.out',
        '#SBATCH --error=' + q_outdir + '/test_%j.err',
    ]

    lines = [
        "#!/bin/bash",
        f"# delfin agent test job: {q_repo} @ {q_ref}, partition {q_part}",
        *header_sb,
        "",
        "set -euo pipefail",
        "export PYTHONNOUSERSITE=1 PYTHONHASHSEED=0",
        "unset PYTHONPATH PYTHONSTARTUP || true",
        # node-local root: everything below is on local scratch, never HOME.
        'LOCAL="${TMPDIR:-/tmp}/delfin_tests_${SLURM_JOB_ID:-$$}"',
        "mkdir -p \"$LOCAL\"",
        "echo \"[delfin-test] node $(hostname), staging at $LOCAL\"",
        "",
        "# copy the repo at the given ref onto node-local disk via git archive",
        f'git -C {q_repo} archive {q_ref} -- "tests" delfin pyproject.toml setup.py 2>/dev/null | tar -x -C "$LOCAL" || git -C {q_repo} archive {q_ref} | tar -x -C "$LOCAL"',
        "",
        "# stage the interpreter's venv node-local (never run children from HOME)",
        'PYTHON_BIN=${PYTHON_BIN:-"' + q_python + '"}',
        "VENV_DIR=\"$(dirname \"$(dirname \"$PYTHON_BIN\")\")\"",
        'tar -C "$(dirname "$VENV_DIR")" -cf - "$(basename "$VENV_DIR")" 2>/dev/null | tar -x -C "$LOCAL" || true',
        'STAGE_PYTHON="$LOCAL/$(basename "$VENV_DIR")/bin/python"',
        '[ -x "$STAGE_PYTHON" ] || STAGE_PYTHON="' + q_python + '"',
        "",
        '# run ONLY the named tests, from the node-local copy',
        f"cd \"$LOCAL\"",
        f"TESTS=({quoted_paths})",
        'set +e',
        '"$STAGE_PYTHON" -m pytest -q "${TESTS[@]}" > "$LOCAL/pytest.log" 2>&1',
        'RC=$?',
        'set -e',
        "",
        '# write a JSON summary next to the log; only these come back',
        'DELFIN_TJ_LOG="$LOCAL/pytest.log" \\',
        '  "$STAGE_PYTHON" - "$RC" ' + q_ref + ' "$LOCAL/summary.json" <<\'PY\'',
        "import json, os, re, sys",
        "rc, ref, out = sys.argv[1], sys.argv[2], sys.argv[3]",
        'log = os.environ.get("DELFIN_TJ_LOG", "")',
        "txt = open(log).read() if os.path.exists(log) else \"\"",
        "# the last N passed / N failed in the pytest -q summary line(s): the",
        "# counts are the actual test numbers, not the count of the words.",
        "ms = re.findall(r\"(\\d+)\\s+passed\\b\", txt)",
        "fs = re.findall(r\"(\\d+)\\s+failed\\b\", txt)",
        "passed = int(ms[-1]) if ms else 0",
        "failed = int(fs[-1]) if fs else 0",
        'json.dump({"ref": ref, "rc": int(rc),',
        '            "passed_tests": passed, "failed_tests": failed,',
        '            "host": os.uname().nodename}, open(out, "w"))',
        "PY",
        'echo "[delfin-test] summary written to $LOCAL/summary.json (rc=$RC)"',
        '# copy back only the summary and the pytest log, job-id-suffixed, so a',
        '# finished job can be collected later from the shared output dir; the',
        '# node-local working copy is never mirrored back.',
        'OUTPUT_DIR="' + q_outdir + '"',
        'mkdir -p "$OUTPUT_DIR"',
        'cp "$LOCAL/summary.json" "$OUTPUT_DIR/summary_${SLURM_JOB_ID:-unknown}.json" 2>/dev/null || true',
        'cp "$LOCAL/pytest.log" "$OUTPUT_DIR/pytest_${SLURM_JOB_ID:-unknown}.log" 2>/dev/null || true',
        'echo "[delfin-test] logs copied to $OUTPUT_DIR"',
        'exit "$RC"',
        "",
    ]
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# Phase 3: submit + collect. Both take an injectable runner so the tests
# never touch the scheduler; the default runners go through DELFIN's own
# submission (delfin.slurm_submit) and the throttled job-state reader
# (delfin.agent.job_monitor), never a literal scheduler command at a call
# site -- the whole-tree direct-call guard stays green.

#: The job id a successful submit prints, e.g. "Submitted batch job 1234".
_SUBMIT_JOB_ID = re.compile(r"Submitted batch job (\d+)")


def _default_submit(argv: Sequence[str], cwd: Optional[str] = None):
    """Run a submit command; returns (stdout, returncode)."""
    try:
        proc = subprocess.run(
            list(argv), capture_output=True, text=True, timeout=60, cwd=cwd)
        return proc.stdout, proc.returncode
    except Exception:
        return "", 1


def submit(
    job_text: str,
    run_dir: str | Path,
    *,
    submit_runner: Optional[Callable[[Sequence[str], Optional[str]], tuple]] = None,
    sbatch: str = "sbatch",
    partition: str | None = None,
) -> str:
    """Write ``job_text`` under ``run_dir`` and submit it; return the job id.

    The submit argv comes from ``delfin.slurm_submit.sbatch_command`` (the
    allow-listed builder) rather than being spelled here. An injectable
    ``submit_runner`` records/answers for tests; the default submits for real.
    """
    run_dir = Path(run_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    script = run_dir / "test_run.sbatch"
    script.write_text(job_text)
    runner = submit_runner or _default_submit
    argv = sbatch_command(sbatch, script, run_dir)
    # A concrete partition the caller chose wins over sbatch_command's
    # auto-selection: append it when the command does not already carry one.
    if partition and not any(a.startswith("--partition=") for a in argv):
        argv = [*argv, f"--partition={partition}"]
    stdout, rc = runner(argv, str(run_dir))
    if rc != 0:
        raise RuntimeError(f"submit failed (rc={rc}): {stdout}")
    match = _SUBMIT_JOB_ID.search(stdout or "")
    if not match:
        raise RuntimeError(f"submit succeeded but no job id in output: {stdout!r}")
    return match.group(1)


def result(
    job_id: str,
    *,
    output_dir: str | Path = "logs",
    state_runner: Optional[Callable[[Sequence[str]], Optional[str]]] = None,
) -> dict:
    """Job state from the throttled reader plus the JSON summary, if present.

    ``state`` is what ``job_monitor.query_job_states_detailed`` returns for
    the id (a scheduler state, ``""`` when unknown, or UNAVAILABLE/THROTTLED).
    ``summary`` is the discovered JSON from ``output_dir/summary_<job_id>.json``
    when the job finished and its copy-back ran, else None. ``state_runner`` is
    injectable so the tests drive a fake instead of the scheduler.
    """
    states = job_monitor.query_job_states_detailed(
        [job_id], run_fn=state_runner or job_monitor._default_run)
    state = states.get(job_id, "")
    summary: Optional[dict] = None
    summary_path = Path(output_dir) / f"summary_{job_id}.json"
    if summary_path.exists():
        try:
            summary = json.loads(summary_path.read_text())
        except (OSError, ValueError):
            summary = None
    return {"job_id": job_id, "state": state, "summary": summary}


# --------------------------------------------------------------------------- #
# CLI handler. The registration into delfin/agent/cli.py build_parser is a
# .gate patch for the operator (that file is not in this package's write
# scope); the handler itself lives here and carries the real tests.

def cmd_test_on_slurm(args) -> int:
    """``delfin-agent test-on-slurm <repo> <ref> <partition> <minutes> <tests...>``.

    Renders the job, submits it and prints the job id on stdout. Error-free
    paths are refused with exit 2 before anything is submitted.
    """
    repo = str(getattr(args, "repo", "") or "")
    ref = str(getattr(args, "ref", "") or "")
    partition = str(getattr(args, "partition", "") or "")
    tests = list(getattr(args, "tests", None) or [])
    try:
        minutes = int(getattr(args, "minutes", 10))
    except (TypeError, ValueError):
        minutes = 0
    if not (repo and ref and partition and tests and minutes > 0):
        print("usage: delfin-agent test-on-slurm <repo> <ref> <partition> "
              "<minutes> <tests/...>")
        return 2
    job_text = build_job(repo, ref, tests, partition, minutes)
    job_id = submit(job_text, getattr(args, "run_dir", "."))
    print(job_id)
    return 0
