"""A detached dashboard server ends itself once nothing holds it.

Nothing holds it when it has no kernel and has had no request for the
idle grace: every open window and every kept dashboard is a kernel, so
either holds the server. The rule is jupyter-server's own
(``shutdown_no_activity_timeout``: no kernel and no API activity for
that long); the grace is a setting, ``DELFIN_DASHBOARD_IDLE_EXIT_MINUTES``
(default 10, ``0`` keeps the server up).

A foreground server keeps today's rule (tests/test_the_server_ends_what_
nobody_watches.py): it stops as soon as its last kernel ends and nothing
is kept, because its terminal is where it was started and Ctrl+C is no
longer there to press. A detached one waits out the grace instead, so a
window closed by mistake finds the same address and login again, and the
next ``delfin-voila`` within minutes opens it rather than a new server.
"""

from __future__ import annotations

import asyncio
import logging
from types import SimpleNamespace

import pytest

from delfin.dashboard import resume_server as R
from delfin.dashboard import server_record as SR
from delfin.dashboard import session as S


@pytest.mark.parametrize("raw, seconds", [
    ("", 600.0), ("5", 300.0), ("0.5", 30.0), ("0", 0.0),
    ("nonsense", 600.0), ("-3", 600.0),
])
def test_the_idle_grace_is_a_setting(monkeypatch, raw, seconds):
    monkeypatch.setenv(SR.IDLE_EXIT_ENV, raw)
    assert SR.idle_exit_seconds() == seconds


def test_the_server_flags_carry_the_grace(monkeypatch):
    monkeypatch.delenv(R.STAY_UP_ENV, raising=False)
    monkeypatch.setenv(SR.IDLE_EXIT_ENV, "5")
    assert "--ServerApp.shutdown_no_activity_timeout=300" in R.cull_config_args()


def test_zero_keeps_the_server_up(monkeypatch):
    monkeypatch.delenv(R.STAY_UP_ENV, raising=False)
    monkeypatch.setenv(SR.IDLE_EXIT_ENV, "0")
    assert not any("shutdown_no_activity_timeout" in a
                   for a in R.cull_config_args())


def test_stay_up_still_wins(monkeypatch):
    monkeypatch.setenv(R.STAY_UP_ENV, "1")
    monkeypatch.setenv(SR.IDLE_EXIT_ENV, "5")
    assert not any("shutdown_no_activity_timeout" in a
                   for a in R.cull_config_args())


class _Base:
    def __init__(self):
        self._kernel_connections = {}
        self._kernels = {}
        self.log = logging.getLogger("test")
        self.ended = []

    def __contains__(self, kid):
        return kid in self._kernels

    async def start_kernel(self, **kwargs):
        kid = f"k{len(self._kernels) + 1}"
        self._kernels[kid] = object()
        self._kernel_connections[kid] = 0
        return kid

    def notify_connect(self, kid):
        self._kernel_connections[kid] += 1

    def notify_disconnect(self, kid):
        self._kernel_connections[kid] -= 1

    async def shutdown_kernel(self, kid, *a, **k):
        self.ended.append(kid)
        self._kernels.pop(kid, None)

    async def shutdown_all(self, *a, **k):
        for kid in list(self._kernels):
            await self.shutdown_kernel(kid)


@pytest.fixture
def manager(tmp_path, monkeypatch):
    monkeypatch.setattr(S, "RECORD_DIR", str(tmp_path / "kept"))
    monkeypatch.delenv(R.STAY_UP_ENV, raising=False)
    m = R.resume_kernel_manager_class(_Base)()
    m.stops = []
    m.parent = SimpleNamespace(stop=lambda: m.stops.append("stop"))
    return m


def test_a_detached_server_does_not_stop_when_its_last_kernel_ends(
        manager, monkeypatch):
    monkeypatch.setenv(SR.DETACHED_ENV, "1")
    kid = asyncio.run(manager.start_kernel())
    asyncio.run(manager.shutdown_kernel(kid))
    assert manager.ended == [kid]
    assert manager.stops == [], "the idle grace ends it, not the last kernel"


def test_a_foreground_server_still_stops_at_once(manager, monkeypatch):
    monkeypatch.delenv(SR.DETACHED_ENV, raising=False)
    kid = asyncio.run(manager.start_kernel())
    asyncio.run(manager.shutdown_kernel(kid))
    assert manager.stops == ["stop"]
