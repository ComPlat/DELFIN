"""R3 phase 4: every job template that runs Python stages its venv node-local.

Finding 2: the cluster operators reported 41.5 M open calls on the shared
workspace file system from a benchmark campaign that ran Python from a venv
on the shared file system per system. A DELFIN job template that runs Python
must therefore get that venv off the shared FS and onto node-local disk (or
state why it does not), and must keep per-item scratch off the shared file
system too.

These are *checking* tests (the templates are R3-owned for checks/comments
only — no behaviour change): they pin the staging/scratch contract the
templates already implement, so a future edit that silently drops node-local
staging fails.

* submit_delfin.sh runs Python (DELFIN pipeline), so it must stage the venv:
  DELFIN_STAGE_VENV defaults ON, the venv is unpacked to
  $STAGE_BASE/delfin_venv_${SLURM_JOB_ID}, and per-item scratch and caches
  are redirected off the shared FS (XDG_CACHE_HOME/MPLCONFIGDIR/PIP_CACHE_DIR
  -> $STAGE_BASE; DELFIN_SCRATCH/ORCA_TMPDIR -> BeeOND/TMPDIR node-local,
  network /scratch only as a warned last resort).
* submit_turbomole.sh runs no Python (it execs the TURBOMOLE binary, default
  ridft), so it documents why no venv staging is needed, and still keeps its
  scratch on node-local $TMPDIR.
"""
from __future__ import annotations

from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def _template(name: str) -> str:
    path = _repo_root() / "delfin" / "submit_templates" / name
    return path.read_text()


def test_delfin_template_stages_venv_to_node_local_disk():
    t = _template("submit_delfin.sh")
    # The staging gate defaults ON (DELFIN_STAGE_VENV=1)...
    assert "DELFIN_STAGE_VENV:-1" in t, "staging must default on"
    # ...and the venv is unpacked to node-local $STAGE_BASE, not run from HOME.
    assert "delfin_venv_" in t
    assert "VENV_LOCAL=\"$STAGE_BASE/" in t, "venv must unpack under $STAGE_BASE"
    # The staged copy is preferred as the Python interpreter.
    assert "VENV_LOCAL/bin/python" in t, "staged python must be used first"


def test_delfin_template_keeps_scratch_off_shared_fs():
    t = _template("submit_delfin.sh")
    # Per-item scratch goes to node-local BeeOND or TMPDIR...
    assert "DELFIN_SCRATCH=\"$BEEOND_MOUNTPOINT/" in t
    assert "DELFIN_SCRATCH=\"$TMPDIR/" in t
    assert "ORCA_TMPDIR=\"$TMPDIR/" in t
    # ...and only falls back to the network filesystem with an explicit warning.
    assert "/scratch/${USER:-$LOGNAME}/delfin_" in t
    assert "WARNING: Using /scratch (network filesystem)" in t


def test_delfin_template_redirects_caches_off_shared_fs():
    t = _template("submit_delfin.sh")
    assert "XDG_CACHE_HOME=\"${XDG_CACHE_HOME:-$STAGE_BASE/" in t
    assert "MPLCONFIGDIR=\"${MPLCONFIGDIR:-$STAGE_BASE/" in t
    assert "PIP_CACHE_DIR=\"${PIP_CACHE_DIR:-$STAGE_BASE/" in t


def test_turbomole_template_documents_no_python_then_keeps_scratch_local():
    t = _template("submit_turbomole.sh")
    # It runs no Python (execs the TURBOMOLE binary) — so no venv staging is
    # needed. It must execute the binary, not source a venv.
    assert "exec \"$TM_COMMAND\"" in t
    # Per-item scratch still stays off the shared FS: node-local $TMPDIR.
    assert "RUN_DIR=\"${TMPDIR:-$JOB_DIR}/turbomole_" in t
