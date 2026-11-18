"""The messenger bridge rebuilds its service connection instead of dying.

The bridge holds one connection to the agent service. Losing it used to
raise out of ``serve()`` and end the process. Measured on the Pi,
2026-11-17: a deploy restarts the agent service, the bridge exited with a
failure 0.33 s later, and systemd's ``Restart=on-failure`` put it back
5 s afterwards. So every deploy cost a traceback and a gap -- and five
such deaths inside the start-limit window would latch the unit in
``failed``, where neither a deploy nor anything else starts it again
(``pi_install.sh`` deliberately never starts a stopped unit). The phone
would then go quiet with nobody told, which is the failure that matters
once the Pi serves more than one person.

Driven through the real ``serve()`` with a stand-in service client whose
pump ends, which is exactly what a restart of the service looks like from
here.
"""

from __future__ import annotations

import asyncio
import os

import pytest

from delfin.agent.service.bridge import Bridge


class _DeadClient:
    """A service client that is closed from the start, so the pump thread
    returns at once -- what the bridge sees when the service went away."""

    closed = True
    events: list = []

    def __init__(self):
        self._cond = asyncio.Condition  # never used: closed short-circuits

    def close(self):
        pass


def test_a_lost_service_connection_is_rebuilt(tmp_path, monkeypatch):
    bridge = Bridge(allowed_uid=os.getuid() + 1)
    monkeypatch.setattr(Bridge, "SERVICE_RECONNECT_WAIT", 0.01)
    tries = []

    def _connect():
        tries.append(1)
        bridge._service = _DeadClient()
        return bridge._service

    monkeypatch.setattr(bridge, "_connect_service", _connect)
    path = tmp_path / "messenger.sock"
    state = {}

    async def _drive():
        task = asyncio.create_task(bridge.serve(path))
        await asyncio.sleep(0.3)
        state["done"] = task.done()
        state["bound"] = path.exists()
        task.cancel()
        await asyncio.gather(task, return_exceptions=True)

    asyncio.run(_drive())
    assert len(tries) > 1, (
        f"the bridge connected {len(tries)} time(s): a lost connection "
        "must be rebuilt, not fatal")
    assert not state["done"], "serve() ended instead of reconnecting"
    assert state["bound"], "the gateway's socket was dropped while reconnecting"


def test_it_gives_up_after_the_cap_so_systemd_sees_a_failure(
        tmp_path, monkeypatch):
    """The cap is the point where papering over stops: a service that
    cannot come back must end as a failure rather than loop for ever."""
    bridge = Bridge(allowed_uid=os.getuid() + 1)
    monkeypatch.setattr(Bridge, "SERVICE_RECONNECT_WAIT", 0.0)
    monkeypatch.setattr(Bridge, "SERVICE_RECONNECTS", 3)

    def _connect():
        bridge._service = _DeadClient()
        return bridge._service

    monkeypatch.setattr(bridge, "_connect_service", _connect)
    async def _drive():
        await asyncio.wait_for(bridge.serve(tmp_path / "messenger.sock"), 5)

    with pytest.raises(ConnectionError) as err:
        asyncio.run(_drive())
    assert "giving up" in str(err.value)


def test_the_cap_is_high_enough_for_a_deploy():
    """A deploy costs one reconnect. The cap must not be so tight that a
    handful of restarts in a session ends the bridge."""
    assert Bridge.SERVICE_RECONNECTS >= 10
    assert 0 < Bridge.SERVICE_RECONNECT_WAIT <= 5
