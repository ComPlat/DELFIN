"""On a cluster login node the agent's tests and builds go to SLURM.

Issue #640: the account was suspended for the load of test suites and
builds that agent sessions ran on a login node. One main session may stay
there (the person's decision, 11.11.); what it runs heavy goes through
`delfin-agent test-on-slurm` or sbatch. The bash gate refuses a test run, a
build or an install there in every mode -- the cluster's rule is not one a
dialog or a grant waves through -- and names the way that works.
"""
from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent import calc_policy as C


def _host(monkeypatch, kind):
    monkeypatch.setattr(C, "host_kind", lambda env=None: kind)


def _bash(ws, cmd, mode):
    perms = A.KitToolPermissions(workspace=ws, mode=mode,
                                 agent_role="solo_agent")
    perms.confirm_callback = lambda *a, **k: True
    return A._doc_executor._gate_bash_call("bash", {"command": cmd}, perms, [])


@pytest.mark.parametrize("mode", ["default", "acceptEdits", "bypassPermissions"])
@pytest.mark.parametrize("cmd", ["pytest -q tests", "python -m pytest -x",
                                 "npx tsc --noEmit", "make -j8",
                                 "pip install foo"])
def test_a_login_node_refuses_tests_and_builds_in_every_mode(
        tmp_path, monkeypatch, mode, cmd):
    _host(monkeypatch, "login_node")
    out = _bash(tmp_path, cmd, mode)
    assert out and "test-on-slurm" in out and "#640" in out


def test_elsewhere_the_same_command_is_not_this_refusal(tmp_path, monkeypatch):
    _host(monkeypatch, "slurm_job")
    out = _bash(tmp_path, "pytest -q tests", "bypassPermissions")
    assert not out or "#640" not in out


def test_reading_and_git_stay_free_on_a_login_node(tmp_path, monkeypatch):
    _host(monkeypatch, "login_node")
    out = _bash(tmp_path, "grep -n pytest setup.cfg", "bypassPermissions")
    assert not out or "#640" not in out
