"""`delfin-agent tui` does not check the principles before the window opens
-- the check reads the digest's copy in api_client, 80 ms of every launch
-- because it is made where it protects: the agent service at its start,
each session's worker at its own, every agent where it is born. A session
over tampered principles does not start, and says why.
"""

from __future__ import annotations

import os
import subprocess
import sys

from conftest import child_env
from tui_worker_client import REPO, WorkerClient


def test_the_window_s_launch_does_not_load_the_principles_check(tmp_path):
    env = child_env(tmp_path)
    env["PYTHONPATH"] = str(REPO)
    out = subprocess.run(
        [sys.executable, "-c",
         "import sys\n"
         "from delfin.agent.service import cli_service\n"
         "cli_service._security_preamble()\n"
         "print('api_client' if 'delfin.agent.api_client' in sys.modules else 'light')"],
        env=env, capture_output=True, text=True, timeout=60)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip().splitlines()[-1] == "light"


def test_a_session_over_tampered_principles_does_not_start_and_says_why(tmp_path):
    hooks = tmp_path / "hooks"
    hooks.mkdir()
    # The principles as tampered: the guard's own check says so.
    (hooks / "sitecustomize.py").write_text(
        "from delfin.agent import principles_guard as g\n"
        "def _tampered(*a, **k):\n"
        "    raise g.PrinciplesTampered('the principles addendum was edited (a test)')\n"
        "g.enforce = _tampered\n")
    env = {"PYTHONPATH": os.pathsep.join([str(hooks), str(REPO)])}
    w = WorkerClient(tmp_path, env=env, initialize=False)
    try:
        ended, _ = w.wait_for("session.ended", timeout=60)
    finally:
        w.close()
    assert "refusing to start" in ended["reason"] and "a test" in ended["reason"], ended
    assert w.proc.returncode == 2
