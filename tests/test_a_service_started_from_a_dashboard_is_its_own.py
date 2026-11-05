"""A service a dashboard kernel starts does not carry the kernel along.

The dashboard starts the person's agent service through
``launcher.start_service`` (via ``ensure_service``). The kernel's
environment holds the server's token and port, its parent server's pid,
the editor's and Jupyter's variables, the lifeline: none of it is the
service's. Taken along, ``stop_all.status()`` read the service as a kernel
whose dashboard had died ("its dashboard server ... is gone") the moment
the page closed -- while running on is what the service is for.

Real path: the launcher starts the supervisor in a private home with such
an environment; /proc/<pid>/environ and a status() run beside it are read.
Control: a plain child with the same environment is reported by status()
as outlived -- the check would fire on that environment.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pytest

from conftest import child_env

REPO = Path(__file__).resolve().parents[1]
KERNEL_ONLY = ("JUPYTER_TOKEN", "JPY_PARENT_PID", "JPY_SESSION_NAME", "DELFIN_VOILA_PORT",
               "BROWSER", "VSCODE_IPC_HOOK_CLI", "DELFIN_LIFELINE_PID", "KERNEL_ID",
               "DELFIN_PREVIEW_PROXY")


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="reads /proc")
def test_the_service_drops_the_kernel_s_environment_and_status_says_what_it_is(tmp_path, gone_pid):
    runtime = Path(tempfile.mkdtemp(prefix="dsk", dir="/tmp"))
    env = child_env(tmp_path)
    env.update(DELFIN_RUNTIME_DIR=str(runtime), PYTHONPATH=str(REPO),
               JUPYTER_TOKEN="t0k", JPY_PARENT_PID=str(gone_pid), JPY_SESSION_NAME="x",
               DELFIN_VOILA_PORT="8866", BROWSER="none", VSCODE_IPC_HOOK_CLI="/x",
               KERNEL_ID="k", DELFIN_PREVIEW_PROXY="/p")
    env.pop("DELFIN_LIFELINE_PID", None)
    start = ("import sys, time\n"
             "from delfin.agent.service import launcher, rundir\n"
             "p = launcher.start_service(idle_exit_s=60)\n"
             "for _ in range(200):\n"
             "    if rundir.service_pid(): break\n"
             "    time.sleep(0.1)\n"
             "print(rundir.service_pid())\n")
    out = subprocess.run([sys.executable, "-c", start], env=env, capture_output=True,
                         text=True, timeout=60, cwd=str(tmp_path))
    pid = int(out.stdout.strip().splitlines()[-1])
    control = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)"], env=env)
    try:
        assert pid > 0, out.stderr
        environ = Path(f"/proc/{pid}/environ").read_bytes().split(b"\0")
        names = {e.split(b"=", 1)[0].decode() for e in environ if b"=" in e}
        assert not names & set(KERNEL_ONLY), names & set(KERNEL_ONLY)
        assert "PYTHONPATH" in names and "HOME" in names          # what it needs stays
        status = subprocess.run(
            [sys.executable, "-c", "import json; from delfin.agent import stop_all; "
                                   "print(json.dumps(stop_all.status()))"],
            env=env, capture_output=True, text=True, timeout=60, cwd=str(tmp_path))
        report = json.loads(status.stdout.strip().splitlines()[-1])
        outlived = {r["pid"]: r["why"] for r in report["outlived_their_start_here"]}
        # Control: the plain child with the kernel's environment is read as
        # one whose dashboard is gone.
        assert "dashboard server" in outlived.get(control.pid, ""), outlived
        assert pid not in outlived, outlived
        assert report["agent_service"]["pid"] == pid
        assert report["agent_service"]["runs_on_its_own"] is True
    finally:
        control.kill()
        control.wait(timeout=10)
        stop = ("from delfin.agent.service.client import ServiceClient\n"
                "with ServiceClient(timeout=5) as c: c.request('service.shutdown', timeout=10)\n")
        subprocess.run([sys.executable, "-c", stop], env=env, capture_output=True, timeout=60,
                       cwd=str(tmp_path))
        deadline = time.monotonic() + 30
        while time.monotonic() < deadline and os.path.exists(f"/proc/{pid}"):
            time.sleep(0.1)
        import shutil
        shutil.rmtree(runtime, ignore_errors=True)


def test_the_dropped_names_are_the_preview_server_s_less_what_the_service_needs():
    from delfin.agent.service import launcher
    from delfin.dashboard import branch_preview as bp

    # The same judgement as for a preview server; the service keeps which
    # DELFIN it imports.
    assert launcher._SERVICE_DROP_EXACT == bp._DROP_EXACT - {"PYTHONPATH", "PYTHONHOME"}
    assert launcher._SERVICE_DROP_PREFIX == bp._DROP_PREFIX
