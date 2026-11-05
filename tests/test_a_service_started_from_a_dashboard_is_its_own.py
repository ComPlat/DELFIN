"""A service a dashboard kernel starts does not carry the kernel along.

The dashboard starts the person's agent service through
``launcher.start_service`` (via ``ensure_service``). The kernel's
environment holds the server's token and port, its parent server's pid,
the editor's and Jupyter's variables, the lifeline: none of it is the
service's. Taken along, ``stop_all.status()`` read the service as a kernel
whose dashboard had died ("its dashboard server ... is gone") the moment
the page closed -- while running on is what the service is for.

Every machine: the launcher's own start is driven with the process start
recorded (subprocess.Popen stood in for), and stop_all's own judgement
(_why_outlived) is asked about the environment the service gets. Control:
the same judgement on the kernel's environment says the dashboard is gone.
"""

from __future__ import annotations

import os
import subprocess

from conftest import child_env

KERNEL_ONLY = ("JUPYTER_TOKEN", "JPY_PARENT_PID", "JPY_SESSION_NAME", "DELFIN_VOILA_PORT",
               "BROWSER", "VSCODE_IPC_HOOK_CLI", "DELFIN_LIFELINE_PID", "KERNEL_ID",
               "DELFIN_PREVIEW_PROXY")


def _kernel_env(tmp_path, gone_pid) -> dict:
    env = child_env(tmp_path)
    (tmp_path / "rt").mkdir(exist_ok=True)
    env.update(DELFIN_RUNTIME_DIR=str(tmp_path / "rt"), PYTHONPATH="/somewhere/delfin",
               JUPYTER_TOKEN="t0k", JPY_PARENT_PID=str(gone_pid), JPY_SESSION_NAME="x",
               DELFIN_VOILA_PORT="8866", BROWSER="none", VSCODE_IPC_HOOK_CLI="/x",
               DELFIN_LIFELINE_PID=str(gone_pid), KERNEL_ID="k", DELFIN_PREVIEW_PROXY="/p")
    return env


def test_the_service_drops_the_kernel_s_environment_and_is_not_read_as_orphaned(
        tmp_path, gone_pid, monkeypatch):
    from delfin.agent import stop_all
    from delfin.agent.service import launcher

    kernel = _kernel_env(tmp_path, gone_pid)
    for name in list(os.environ):
        monkeypatch.delenv(name, raising=False)
    for name, value in kernel.items():
        monkeypatch.setenv(name, value)
    started = []

    class _Popen:
        def __init__(self, cmd, **kw):
            started.append({"cmd": cmd, "env": dict(kw.get("env") or {})})
            self.pid = 4242

    monkeypatch.setattr(subprocess, "Popen", _Popen)
    launcher.start_service(idle_exit_s=60)
    (run,) = started
    assert "delfin.agent.service.supervisor" in " ".join(run["cmd"])
    env = run["env"]
    assert not set(env) & set(KERNEL_ONLY), set(env) & set(KERNEL_ONLY)
    assert env["PYTHONPATH"] == "/somewhere/delfin" and env["HOME"] == kernel["HOME"]
    # Control: the kernel's own environment is judged orphaned ...
    assert stop_all._why_outlived(4243, kernel)
    # ... the service's is not.
    assert stop_all._why_outlived(4242, env) == ""


def test_status_names_the_service_as_running_on_its_own(tmp_path, monkeypatch):
    import json

    from delfin.agent import stop_all
    from delfin.agent.service import rundir

    roster = tmp_path / "h.roster"
    roster.write_text(json.dumps({"sessions": [{"sessionId": "a"}, {"sessionId": "b"}],
                                  "dormant": [{"summary": {"sessionId": "c"}}]}))
    monkeypatch.setattr(rundir, "service_pid", lambda: 4242)
    monkeypatch.setattr(rundir, "read_run_file", lambda host=None: {
        "pid": 4242, "socket": "/run/x/agent.sock", "startedAt": 1.0})
    monkeypatch.setattr(rundir, "roster_file", lambda host=None: roster)
    got = stop_all._agent_service_here()
    assert got == {"pid": 4242, "socket": "/run/x/agent.sock", "started_at": 1.0,
                   "sessions": 2, "asleep": 1, "runs_on_its_own": True}
    # Control: no service, nothing named.
    monkeypatch.setattr(rundir, "service_pid", lambda: 0)
    assert stop_all._agent_service_here() is None


def test_the_dropped_names_are_the_preview_server_s_less_what_the_service_needs():
    from delfin.agent.service import launcher
    from delfin.dashboard import branch_preview as bp

    # The same judgement as for a preview server; the service keeps which
    # DELFIN it imports.
    assert launcher._SERVICE_DROP_EXACT == bp._DROP_EXACT - {"PYTHONPATH", "PYTHONHOME"}
    assert launcher._SERVICE_DROP_PREFIX == bp._DROP_PREFIX
