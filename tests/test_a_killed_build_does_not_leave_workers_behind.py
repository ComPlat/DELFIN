"""A build killed on the clock must not leave its workers running.

``subprocess.run(timeout=)`` kills the direct child and nothing else.  The
build child is not a leaf: ``smiles_converter`` opens a ProcessPoolExecutor for
batch UFF with up to ``DELFIN_MAX_PROCESS_WORKERS`` (64) workers.  Killing only
the parent reparents those to init, where they hold RAM until somebody notices
-- measured elsewhere in this tree as 128 orphans alive for three to five hours
after one such kill.

The second thing pinned here is the timeout itself.  It was hard-coded at
1800 s and nothing ever set it, and a cut build returns *nothing* rather than a
smaller answer.  Measured over 2000 systems at extreme with max_isomers=0,
1800 s cut 9.6 % of 51-80 atom complexes and 39.6 % of those above 80 atoms --
the band that holds porphyrins, terpyridines and cyclams, which is what this
pipeline exists to build.
"""

from __future__ import annotations

import os
import signal
import subprocess
import sys
import time

import pytest

from delfin.dashboard import input_processing


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except (OSError, ProcessLookupError):
        return False
    return True


def test_killing_a_build_takes_its_children_with_it(tmp_path):
    # A parent that spawns a long-lived child and then hangs, standing in for
    # the builder and its UFF pool.
    marker = tmp_path / "child.pid"
    script = (
        "import os, subprocess, sys, time\n"
        f"kid = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(300)'])\n"
        f"open({str(marker)!r}, 'w').write(str(kid.pid))\n"
        "time.sleep(300)\n"
    )
    proc = subprocess.Popen([sys.executable, "-c", script],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            start_new_session=True)
    deadline = time.monotonic() + 15
    while time.monotonic() < deadline and not marker.exists():
        time.sleep(0.05)
    assert marker.exists(), "the stand-in never spawned its child"
    child_pid = int(marker.read_text().strip())
    assert _alive(child_pid)

    input_processing._kill_process_group(proc)

    deadline = time.monotonic() + 10
    while time.monotonic() < deadline and (_alive(child_pid) or proc.poll() is None):
        time.sleep(0.05)
    assert proc.poll() is not None, "the build itself survived"
    assert not _alive(child_pid), (
        f"pid {child_pid} outlived the build it belonged to; that is the "
        "orphan this guards against")


def test_killing_something_already_gone_is_not_an_error():
    proc = subprocess.Popen([sys.executable, "-c", "pass"], start_new_session=True)
    proc.wait()
    input_processing._kill_process_group(proc)      # must not raise


def test_the_build_timeout_is_not_the_one_that_cut_the_porphyrins():
    # 1800 s returned nothing for 39.6 % of complexes above 80 atoms.
    assert input_processing._UI_ISOLATE_TIMEOUT >= 3600


def test_the_timeout_can_still_be_switched_off_entirely():
    # A complete manifold on a heavily substituted macrocycle is a long
    # deterministic construction, not a hang, and has to be allowed to finish.
    assert "DELFIN_UI_ISOLATE_TIMEOUT" in open(
        input_processing.__file__, encoding="utf-8").read()


def test_the_control_key_reaches_the_enforcement_point():
    """MANTA_TIME_BUDGET was parsed, defaulted, documented -- and never read.

    It was removed on the grounds that the builder has no interruption point,
    which was wrong: the build already runs in an isolated subprocess, and that
    subprocess is exactly the interruption point. The key now sets the timeout
    that subprocess is killed on.
    """
    from delfin.common.manta_settings import apply_construction_env

    applied = {}
    apply_construction_env({"MANTA_TIME_BUDGET": "7200"}, applied)
    assert applied["DELFIN_UI_ISOLATE_TIMEOUT"] == "7200"

    # 0 means no limit, and has to survive as 0 rather than becoming a default
    applied = {}
    apply_construction_env({"MANTA_TIME_BUDGET": "0"}, applied)
    assert applied["DELFIN_UI_ISOLATE_TIMEOUT"] == "0"

    # unset leaves the module default alone rather than pinning it
    applied = {}
    apply_construction_env({}, applied)
    assert "DELFIN_UI_ISOLATE_TIMEOUT" not in applied


def test_the_template_ships_the_budget_it_documents():
    from delfin import define

    assert "MANTA_TIME_BUDGET=3600" in define.TEMPLATE
