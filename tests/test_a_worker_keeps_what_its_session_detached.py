"""A worker is a child subreaper (Linux): a process an agent command
detaches is re-parented to the worker, not to pid 1, so the service's
ancestry check (supervisor._from_an_agent) still finds the worker above
it. The worker reaps those adopted orphans, and never the exit status of
a child it started itself.

Each arm runs in a fresh interpreter: the subreaper flag is per process
and must not leak into the test process. The control arm runs the same
detach without the flag and expects the orphan elsewhere.
"""

from __future__ import annotations

import subprocess
import sys
import textwrap
from pathlib import Path


from delfin.agent.service import subreaper

REPO = Path(__file__).resolve().parents[1]

def _linux() -> bool:
    """PR_SET_CHILD_SUBREAPER is Linux's. Elsewhere each probe test asserts
    the other half of the contract instead of skipping: nothing is
    claimed and no reaper runs."""
    if sys.platform.startswith("linux"):
        return True
    assert subreaper.become_subreaper() is False
    assert subreaper.keep_descendants() is False
    return False

PROBE = textwrap.dedent("""
    import os, subprocess, sys, time
    from delfin.agent.service import subreaper
    flag = sys.argv[1] == "on"
    if flag:
        assert subreaper.become_subreaper()
    reaper = subreaper.Reaper()
    # A shell that starts a background sleep, says its pid, and waits for
    # a line before it exits: the sleep is then an orphan.
    sh = subprocess.Popen(["sh", "-c", "sleep 0.6 & echo $!; read x"],
                          stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    orphan = int(sh.stdout.readline())
    reaper.reap_once()                       # sees it below the shell
    sh.stdin.write("go\\n"); sh.stdin.flush(); sh.wait()
    time.sleep(0.1)
    with open(f"/proc/{orphan}/stat") as fh:
        ppid = int(fh.read().rsplit(")", 1)[1].split()[1])
    print("PARENT", "me" if ppid == os.getpid() else "other")
    # A child of our own whose status its Popen waits for: never taken.
    own = subprocess.Popen(["sh", "-c", "sleep 0.3; exit 3"])
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline and os.path.exists(f"/proc/{orphan}"):
        reaper.reap_once()
        time.sleep(0.05)
    print("ORPHAN", "gone" if not os.path.exists(f"/proc/{orphan}") else "zombie")
    print("OWN", own.wait())
""")


def _probe(tmp_path, flag: str) -> dict:
    script = tmp_path / "probe.py"
    script.write_text(PROBE)
    out = subprocess.run([sys.executable, str(script), flag], capture_output=True,
                         text=True, timeout=60, cwd=str(tmp_path),
                         env={"PYTHONPATH": str(REPO), "PATH": "/usr/bin:/bin"})
    assert out.returncode == 0, out.stderr
    return dict(line.split(" ", 1) for line in out.stdout.split("\n") if line)


def test_a_detached_process_stays_below_the_worker_and_is_reaped(tmp_path):
    if not _linux():
        return
    seen = _probe(tmp_path, "on")
    assert seen["PARENT"] == "me"
    assert seen["ORPHAN"] == "gone"
    assert seen["OWN"] == "3"


def test_control_without_the_flag_the_orphan_leaves_the_tree(tmp_path):
    if not _linux():
        return
    seen = _probe(tmp_path, "off")
    assert seen["PARENT"] == "other"
    assert seen["OWN"] == "3"


def test_elsewhere_nothing_is_claimed(monkeypatch):
    monkeypatch.setattr(subreaper.sys, "platform", "darwin")
    assert subreaper.become_subreaper() is False
    assert subreaper.keep_descendants() is False


def test_the_worker_becomes_a_subreaper_before_its_session():
    from delfin.agent.service import worker
    import inspect
    source = inspect.getsource(worker.main)
    assert source.index("keep_descendants()") < source.index("open_real_session(")
