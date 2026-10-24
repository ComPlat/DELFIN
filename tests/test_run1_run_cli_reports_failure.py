"""RUN1 control: every w21 run_cli reports WHY it failed (QS via s25/s29).

w21_mlp_screen.py run_cli and w21_turbomole_workflow.py run_cli caught
Exception and returned 1 with NO message, while
w21_solvation_workflow.py:155-158 prints WHY it failed (wave-17 lesson
5: a bare failure with no message hides the cause).  A user whose run
dies sees nothing.  Expected: every w21 run_cli prints the failure
reason to stderr before returning 1.  Red through the gate on the
UNFIXED commit (rc=1, empty stderr, observed for both).
"""

import contextlib
import io
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from delfin.workflows.contrib.w21_mlp_screen import (  # noqa: E402
    MlpScreenWorkflow,
)
from delfin.workflows.contrib.w21_turbomole_workflow import (  # noqa: E402
    TurbomoleWorkflow,
)


def _run_and_capture_stderr(workflow, argv):
    """Call run_cli on argv that raises inside run(); return its stderr."""
    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        rc = workflow.run_cli(argv)
    return rc, err.getvalue()


def test_mlp_screen_run_cli_reports_failure_reason(tmp_path):
    # A nonexistent ensemble file makes run() raise inside the CLI.
    rc, err = _run_and_capture_stderr(
        MlpScreenWorkflow(),
        [str(tmp_path / "missing.xyz"), "--backend", "ani2x",
         "--work-dir", str(tmp_path / "wd")])
    assert rc == 1, (rc, err)
    assert "mlp_screen" in err and "failed" in err, err


def test_turbomole_run_cli_reports_failure_reason(tmp_path):
    # A nonexistent control input makes run() raise inside the CLI.
    rc, err = _run_and_capture_stderr(
        TurbomoleWorkflow(),
        [str(tmp_path / "missing.xyz"),
         "--work-dir", str(tmp_path / "wd")])
    assert rc == 1, (rc, err)
    assert "turbomole" in err and "failed" in err, err
