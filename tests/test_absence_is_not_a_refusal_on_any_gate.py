"""A window that nobody answered means the user is away, on every gate.

The distinction was drawn for bash commands and nowhere else. The other
gates checked it — to decide whether to remember the denial — and then
told the model "user denied" all the same. In one afternoon two
approvals for a protected path expired unanswered, five minutes each,
and the model was told it had been refused (2026-09-17).

  a protected path           timeout says away, not denied
  an ordinary write          the same
  reaching the network       the same
  an MCP side effect         the same
  a real denial              still says denied, and is remembered
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor


class _Broker:
    """A confirm callback with the broker's own shape.

    The gate reads ``last_timed_out`` off the bound method's owner, which
    is how it tells an expired window from a click on Deny.
    """

    def __init__(self, timed_out: bool) -> None:
        self.last_timed_out = timed_out

    def callback(self, name, args, preview):
        return False


@pytest.fixture
def ws(tmp_path):
    w = tmp_path / "ws"
    w.mkdir()
    (w / "app.py").write_text("x = 1\n")
    return w


def _perms(ws, timed_out: bool):
    return KitToolPermissions(workspace=ws, mode="default",
                              confirm_callback=_Broker(timed_out).callback)


def _write_gate(ws, timed_out: bool, path="app.py"):
    perms = _perms(ws, timed_out)
    return _doc_executor._run_permission_gate(
        "write_file", {"path": path, "content": "y = 2\n"}, perms)


def test_a_write_nobody_answered_is_not_a_denial(ws):
    out = _write_gate(ws, True)
    assert out is not None
    assert "TIMED OUT" in out and "NOT a denial" in out
    assert "user denied" not in out


def test_a_write_the_user_refused_still_says_so(ws):
    out = _write_gate(ws, False)
    assert out is not None and "user denied" in out


def test_a_refusal_is_remembered_and_an_absence_is_not(ws):
    perms = _perms(ws, True)
    _doc_executor._run_permission_gate(
        "write_file", {"path": "app.py", "content": "y\n"}, perms)
    assert not perms.denied_actions

    perms = _perms(ws, False)
    _doc_executor._run_permission_gate(
        "write_file", {"path": "app.py", "content": "y\n"}, perms)
    assert perms.denied_actions, "a real denial is what must not be retried"


def test_every_gate_that_asks_carries_the_absence_wording():
    """No gate may report absence as refusal: a model treats a denial as
    final and a timeout as "come back later"."""
    import inspect

    from delfin.agent import api_client

    source = inspect.getsource(api_client)
    assert source.count("this is NOT a denial") >= 4, (
        "the write, protected-path, egress and MCP gates each say it")
