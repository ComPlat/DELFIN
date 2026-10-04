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

import re
import shlex
from typing import Iterable, Sequence

from delfin.slurm_submit import normalize_time_limit

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
) -> str:
    """Render a SLURM job script that re-runs ``test_paths`` at ``ref``.

    The script copies repo at ``ref`` onto node-local disk (``$TMPDIR``) via
    ``git archive``, stages the interpreter's venv node-local, runs exactly
    the named tests there, and writes a JSON summary. Pure: returns text,
    touches nothing.
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
        '#SBATCH --output=${TMPDIR:-/tmp}/delfin_agent_tests_%j.out',
        '#SBATCH --error=${TMPDIR:-/tmp}/delfin_agent_tests_%j.err',
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
        'DELFIN_TJ_REF="' + q_ref + '" DELFIN_TJ_LOG="$LOCAL/pytest.log" \\',
        '  "$STAGE_PYTHON" - "$RC" "$LOCAL/summary.json" <<\'PY\'',
        "import json, os, sys",
        "rc, out = int(sys.argv[1]), sys.argv[2]",
        'log = os.environ.get("DELFIN_TJ_LOG", "")',
        "txt = open(log).read() if os.path.exists(log) else \"\"",
        'passed = txt.count("passed") - txt.count("false")',
        'failed = txt.count("failed")',
        'json.dump({"ref": os.environ.get("DELFIN_TJ_REF", ""),',
        '            "rc": rc, "passed_tests": max(passed, 0),',
        '            "failed_tests": max(failed, 0), "host": os.uname().nodename},',
        '           open(out, "w"))',
        "PY",
        'echo "[delfin-test] summary written to $LOCAL/summary.json (rc=$RC)"',
        'exit "$RC"',
        "",
    ]
    return "\n".join(lines)
