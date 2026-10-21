"""``gate explain`` gives the executor's reason, whatever ran before it.

``_bash_isolation_argv`` records once per process what holds the network
of the commands (``_announce_child_network``: kind ``network`` or
``network_missing``, ``blocked=False``). It describes the host, not the
command. ``boundaries.classify`` took the first recorded kind as the
reason code of an allowed command, so the FIRST bash command classified
in a process read ``network`` and every later one read ``allowed``. Run
alone, ``test_explain_agrees_with_the_executor[default-terminal]`` failed
on main: explain (first) said ``echo hi`` -> ``allow (network)``, the
executor (second) ``allow (allowed)``. In a long run another test had
spent the notice already, and the two agreed.

The dry run also spent the notice for real: a session started in the
same process after an explain never announced its network.

Input: ``_CHILD_NETWORK_ANNOUNCED`` reset to its process-start value
before the first of the two steps. Output: explain and the executor give the same verdict
and reason code for a plain command in either order; explain leaves the
flag as it found it.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from delfin.agent import api_client as ac
from delfin.agent import boundaries as B
from delfin.agent import security_events
from delfin.agent.api_client import _DocToolExecutor


def _ws(root: Path) -> Path:
    ws = root / "ws"
    ws.mkdir(parents=True)
    subprocess.run(["git", "init", "-q", str(ws)], check=True,
                   capture_output=True)
    return ws


def _executed(ws: Path, command: str, monkeypatch) -> dict:
    events: list = []
    real = security_events.record

    def _capture(kind, tool, detail, *, blocked=True):
        events.append({"kind": kind, "blocked": bool(blocked)})

    monkeypatch.setattr(security_events, "record", _capture)
    perms, _ = B.session_permissions(ws, permission_mode="default",
                                     asker="terminal",
                                     confirm_callback=lambda *a: False)
    try:
        result = _DocToolExecutor().execute("bash", {"command": command}, perms)
    finally:
        monkeypatch.setattr(security_events, "record", real)
    return B.classify("bash", result, [], events)


@pytest.mark.parametrize("command", ["echo hi", "pwd"])
@pytest.mark.parametrize("explain_first", [True, False])
def test_explain_and_executor_agree_with_a_fresh_notice(
        tmp_path, monkeypatch, command, explain_first):
    ws = _ws(tmp_path)
    monkeypatch.setattr(ac, "_CHILD_NETWORK_ANNOUNCED", False)
    # Fresh for the first step only: the order of the failing test, where
    # the first command classified in the process meets the notice.
    if explain_first:
        said = B.explain("bash", {"command": command}, cwd=ws,
                         permission_mode="default", asker="terminal")
        done = _executed(ws, command, monkeypatch)
    else:
        done = _executed(ws, command, monkeypatch)
        said = B.explain("bash", {"command": command}, cwd=ws,
                         permission_mode="default", asker="terminal")
    assert (said["verdict"], said["reason_code"]) == (
        done["verdict"], done["reason_code"])
    assert said["reason_code"] not in ("network", "network_missing"), said


def test_classify_does_not_read_a_host_notice_as_the_reason():
    for kind in ("network", "network_missing"):
        out = B.classify("bash", '{"exit_code": 0, "stdout": "hi"}', [],
                         [{"kind": kind, "blocked": False}])
        assert (out["verdict"], out["reason_code"]) == ("allow", "allowed")


def test_explain_leaves_the_notice_for_the_session(tmp_path, monkeypatch):
    ws = _ws(tmp_path)
    monkeypatch.setattr(ac, "_CHILD_NETWORK_ANNOUNCED", False)
    B.explain("bash", {"command": "echo hi"}, cwd=ws,
              permission_mode="default", asker="terminal")
    assert ac._CHILD_NETWORK_ANNOUNCED is False
