"""The agent child environment names the session scratch.

Every process an agent session starts is built from _scrubbed_bash_env
(api_client.py:29856). Inside a cage the session scratch dir is exported
as TMPDIR/TMP/TEMP and bound read-write into the cage; code inside the
cage (e.g. the ORCA scratch choice, delfin/orca.py) must find it without
guessing it from TMPDIR, because a cage layer may remap TMPDIR. So the
environment carries the explicit name DELFIN_SESSION_SCRATCH, together
with the cage marker DELFIN_AGENT_SESSION (approval_guard.ENV_MARKER).
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client
from delfin.agent import approval_guard


def test_agent_child_env_names_the_session_scratch(monkeypatch, tmp_path):
    monkeypatch.delenv("DELFIN_SESSION_SCRATCH", raising=False)
    monkeypatch.delenv("TMPDIR", raising=False)
    env = api_client._scrubbed_bash_env()
    scratch = env["TMPDIR"]
    assert scratch and scratch != tmp_path
    assert env["DELFIN_SESSION_SCRATCH"] == scratch
    assert env[approval_guard.ENV_MARKER]


def test_when_the_scratch_cannot_be_made_no_variable_is_exported(
        monkeypatch):
    """The block sits in one try; if _session_scratch_dir fails, no
    TMPDIR/DELFIN_SESSION_SCRATCH may appear -- a half-export would point
    code at a directory that does not exist."""
    monkeypatch.setattr(api_client, "_session_scratch_dir",
                        lambda: (_ for _ in ()).throw(RuntimeError("no scratch")))
    monkeypatch.delenv("TMPDIR", raising=False)
    env = api_client._scrubbed_bash_env()
    assert "DELFIN_SESSION_SCRATCH" not in env
    assert "TMPDIR" not in env
