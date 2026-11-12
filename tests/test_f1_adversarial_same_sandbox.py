"""Adversarial contract tests for Target 4 (same sandbox on the node) —
reviewer s2.

The brief: the durable (SLURM-run) supervisor must start the SAME
launcher/supervisor with the SAME settings as on the login node — never
a wider permission. On main 04bf9703b no durable job renderer exists, so
the module-existence case is red; the rest arms for the builder hash.

Contract points, derived from the launcher surface (launcher.py:65-93:
start_service drops _SERVICE_DROP_* env keys and runs the supervisor in
rundir.ensure_runtime_dir()):
- the rendered job calls `python -m delfin.agent.service.supervisor`
  (the same command launcher.start_service uses);
- it does NOT strip or widen environment beyond what the login path
  drops (_SERVICE_DROP_EXACT/_SERVICE_DROP_PREFIX semantics);
- it does not run as root / with sudo, and contains no `sudo`;
- the venv interpreter in the job is the same one that renders it.
"""

from __future__ import annotations

import sys

import pytest


def _module():
    try:
        from delfin.agent.service import durable_host
        return durable_host
    except ImportError:
        return None


def _rendered_job():
    mod = _module()
    if mod is None:
        pytest.skip("no durable_host module on this branch (red state)")
    render = getattr(mod, "render_durable_job", None)
    if render is None:
        pytest.skip("durable_host has no render_durable_job yet")
    return render()


class TestSameSandboxOnTheNode:
    def test_durable_host_module_exists(self):
        assert _module() is not None, (
            "Target 4: no durable_host module yet — the same-sandbox "
            "contract cannot hold")

    def test_job_starts_the_same_supervisor_command(self):
        job = _rendered_job()
        assert "delfin.agent.service.supervisor" in job, (
            "the durable job must start the same supervisor process the "
            "launcher starts (launcher.py:81), never a second runner")

    def test_job_contains_no_sudo_or_root_escalation(self):
        job = _rendered_job()
        lowered = job.lower()
        assert "sudo" not in lowered
        assert "#SBATCH --uid" not in job and "--uid" not in job

    def test_job_does_not_export_an_explicit_var_list(self):
        job = _rendered_job()
        assert "--export=ALL," not in job
        assert "--export=ALL\n" in job or "--export=ALL " in job or \
            "--export=ALL\"" in job, "the job must carry --export=ALL exactly"

    def test_job_runs_nothing_without_the_set_guard(self):
        job = _rendered_job()
        assert "set -euo pipefail" in job

    def test_job_interpreter_matches_the_rendering_interpreter(self):
        job = _rendered_job()
        assert sys.executable in job or "$PYTHON_BIN" in job, (
            "the durable job must use the same interpreter (same venv, "
            "same settings), not a system python")

    def test_job_carries_no_extra_env_planting(self):
        """Never wider: the job script must not PLANT privilege-related
        variables (an env line granting write scope the session did not
        have)."""
        job = _rendered_job()
        for bad in ("DELFIN_WRITE_ALLOW", "DELFIN_EXTRA_DIR",
                    "DELFIN_BYPASS", "DELFIN_UNSAFE"):
            assert bad not in job
