"""A test's child processes must not write into the real user state.

Measured after the suite runs of 2026-09-22: nine confirmation requests
in the real ``~/.delfin/terminal_confirmations``, none with a
``session_id``, the host column naming the compute node. The suite
redirects its own process's sinks in memory (``_isolate_user_state``,
from the one product table), but that redirect exists only in the test
process. A child process a test starts -- the voila servers, the agent
scripts, any ``sys.executable -c`` importing DELFIN -- resolves its
sinks at import time in the REAL home.

The probe here is read-only on purpose: a child resolves the module
constant and prints it, and nobody writes into the real home, not even
while the test is red. Where the child's answer points is the whole
question.

DELFIN already has the environment route for exactly this
(``DELFIN_SCRATCH_STATE``, honoured by the agent CLI before any sink is
imported); the suite's children simply never receive it, and the entry
points they actually start do not all honour it.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

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


def _child_answer() -> list[Path]:
    """Start a child the way the suite's own server tests do: full
    ``os.environ`` (the real HOME among it), the checkout on
    PYTHONPATH. Read-only: nothing is asked, nothing is written."""
    env = dict(os.environ)
    env["PYTHONPATH"] = str(_REPO)
    out = subprocess.run([sys.executable, "-c", _PROBE], env=env,
                         capture_output=True, text=True, timeout=120,
                         cwd=str(_REPO))
    assert out.returncode == 0, out.stderr[-800:]
    return [Path(line) for line in out.stdout.strip().splitlines()]


def test_a_child_process_resolves_its_sinks_outside_the_real_home():
    """The control. A child a test starts inherits no in-memory redirect:
    whatever redirection reaches it must travel through the environment.
    Today nothing does, so the child answers with the real home's paths
    -- the nine left-over confirmations of 2026-09-22 came from exactly
    such children."""
    real_home = Path.home()
    for resolved in _child_answer():
        assert resolved != real_home and real_home not in resolved.parents, (
            f"a child process of a test resolves {resolved} into the real "
            "home; the suite's in-memory redirect stops at the process "
            "boundary, so a child publishing a confirmation lands it in "
            "the operator's room (measured 2026-09-22)")
