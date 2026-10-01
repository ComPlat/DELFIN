"""A test's child processes must not write into the real user state.

Measured after the suite runs of 2026-09-22: nine confirmation requests
in the real ``~/.delfin/terminal_confirmations``, none with a
``session_id``, the host column naming the compute node. The suite
redirects its own process's sinks in memory (``_isolate_user_state``,
from the one product table), but that redirect exists only in the test
process. A child process a test starts -- the voila servers, the agent
scripts, any ``sys.executable -c`` importing DELFIN -- resolves its
sinks at import time in the REAL home.

The fix is ``conftest.child_env``: a private HOME (what every
``Path.home()`` sink follows, the confirmation room among them) plus
``DELFIN_SCRATCH_STATE`` (the product's own environment route, honoured
by the agent CLI). No product code changes: the child is simply pointed
at a home of its own.

The guards here are read-only on purpose: the probe children resolve
and print their sink paths, and nobody writes into the real home, not
even while a guard is red.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

from conftest import child_env

_REPO = Path(__file__).resolve().parents[1]

# The measured sink first, then one more from the product's sink table
# to pin that the mechanism -- not one attribute -- is what a child
# inherits.
_PROBE = (
    "from delfin.agent import terminal_confirm as tc\n"
    "from delfin.agent import outcome_tracker as ot\n"
    "print(tc._PENDING_DIR)\n"
    "print(ot._DEFAULT_PATH)\n"
)


def _child_answer(env: dict) -> list[Path]:
    """Start a child the way the suite's own server tests do. Read-only:
    nothing is asked, nothing is written."""
    env = dict(env)
    env.setdefault("PYTHONPATH", str(_REPO))
    out = subprocess.run([sys.executable, "-c", _PROBE], env=env,
                         capture_output=True, text=True, timeout=120,
                         cwd=str(_REPO))
    assert out.returncode == 0, out.stderr[-800:]
    return [Path(line) for line in out.stdout.strip().splitlines()]


def test_a_bare_child_still_points_at_the_real_home(tmp_path):
    """The incident itself, as a control on the unredirected form: a
    child started with the raw ``os.environ`` (HOME among it) resolves
    its sinks into the real home. This stays red-free by construction --
    it asserts the LEAK's direction, so a future reader sees what
    ``child_env`` is for; the guard that must hold is the next one."""
    real_home = Path.home()
    leaked = [p for p in _child_answer(dict(os.environ))
              if real_home in p.parents]
    assert leaked, (
        "a bare child no longer resolves any sink into the real home; "
        "the control premise is gone and this file should be revisited")


def test_a_child_with_child_env_resolves_outside_the_real_home(tmp_path):
    """The guard. A child started through ``child_env`` -- the one way
    the suite hands a test's children their own state -- must not
    resolve a single sink into the real home, the confirmation room of
    2026-09-22 included. On the unredirected state of the suite this
    failed: the child answered with the real home's paths."""
    real_home = Path.home()
    for resolved in _child_answer(child_env(tmp_path)):
        assert resolved != real_home and real_home not in resolved.parents, (
            f"a child process of a test resolves {resolved} into the real "
            "home; the suite's in-memory redirect stops at the process "
            "boundary, so a child publishing a confirmation lands it in "
            "the operator's room (measured 2026-09-22)")


# ---------------------------------------------------------------------------
# The scan: every future child start is covered too
# ---------------------------------------------------------------------------

# A call that starts a child process of the test's own interpreter
# (quoted or bare -- the suite spells both), or an installed delfin
# entry point. Plain shell commands (sleep, git in a tmp repo) are not
# DELFIN children and resolve no sinks.
_CHILD_START = re.compile(
    r"(subprocess\.(run|Popen|check_output)|C\.run)\s*\(\s*\[?\s*[\"']?\s*"
    r"(sys\.executable)"
)
# How far AFTER the start call the child's own code can reach. Only
# after: a window before the call catches the test file's own imports
# (measured: it flagged test_orphaned_orca_processes... and others
# whose children are bare sleeps), which say nothing about the child.
_CHILD_WINDOW = 800

#: Files allowed to start DELFIN children without child_env, each with
#: the reason its children cannot reach the real user state.
_CHILD_ENV_EXEMPT: dict[str, str] = {
    # This file's own control starts a bare child on purpose, with the
    # raw environment: it is the measured leak itself, read-only.
    "test_a_child_process_writes_into_no_real_home.py":
        "starts its bare-child control on purpose, read-only",
    # The flagged child runs `python -m pytest` INSIDE a benchmark
    # fixture repo (ensemble_tools) whose tests import no delfin; the
    # window matched the test function's own delfin import instead.
    "test_the_control_task_scores_what_it_means_to.py":
        "its child runs pytest in a fixture repo that imports no delfin",
}


def _delfin_children(path: Path, text: str) -> int:
    """The line of the first child start whose started code is DELFIN's.

    A bare ``sys.executable`` child that never imports delfin (a sleep,
    a probe, a subprocess to re-read a file) resolves no ``~/.delfin``
    sink and is not this guard's business -- the incident's children
    were the ones running DELFIN's own modules.

    Where the started code sits: right after the call (an argv list
    with ``-m delfin...``, or an inline ``-c`` script), or just above
    it (a script variable the call passes). The backward search is
    bounded by the enclosing ``def``/``class`` so the file's own
    module-level imports -- which say nothing about the child -- are
    never in the window (measured: without the bound, files whose
    children are bare sleeps were flagged by their own import lines).
    The forward window deliberately crosses function boundaries: the
    suite spells child scripts as module-level constants too, and a
    window that stops at the next def missed one of those.
    """
    import re as _re
    for match in _CHILD_START.finditer(text):
        head = text[:match.start()]
        starts = [m.start() for m in _re.finditer(r"\n(?=def |class )", head)]
        func_start = starts[-1] + 1 if starts else 0
        window = text[func_start:match.start() + _CHILD_WINDOW]
        if _re.search(r"(-m[\"'],\s*[\"']delfin|import\s+delfin|from\s+delfin)",
                      window):
            return head.count("\n") + 1
    return 0


def test_every_delfin_child_start_passes_child_env():
    """The future half of the guard: scan the suite's sources for child
    starts whose code is DELFIN's own, and require the environment they
    pass to come from ``child_env`` (or name why not, in the table
    above). A test written next week that starts a DELFIN child with
    the raw environment is caught here, in CI, before it leaves
    confirmations in the operator's room."""
    offenders = []
    for path in sorted(_REPO.joinpath("tests").glob("*.py")):
        text = path.read_text(encoding="utf-8", errors="replace")
        if path.name in _CHILD_ENV_EXEMPT:
            continue
        if "child_env" not in text:
            line = _delfin_children(path, text)
            if line:
                offenders.append(f"{path.name}:{line}")
    assert not offenders, (
        "test files start DELFIN child processes without child_env "
        "(the child resolves its ~/.delfin sinks against the real "
        "home):\n  " + "\n  ".join(offenders))
