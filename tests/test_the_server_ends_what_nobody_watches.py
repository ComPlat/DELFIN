"""A kernel whose last window closed is ended by the server, unless kept.

This replaced a heartbeat in the page and a watchdog thread in the
kernel. Four of the seven things a browser found in that design were
its own mechanics: the beat had to be displayed, the first beat raced
the field it wrote to, exiting on our own invited the restarter, and
disarming judged the page by a stale beat. The server already counts
websocket connections per kernel, and a count of zero for longer than
the grace is the same fact with none of those parts.

What is pinned:

  a window keeps it            connections > 0 -> never ended
  no window ends it            zero connections past the grace -> ended
  kept is the exemption        a record naming the kernel -> not ended
  a reload is not a departure  zero connections INSIDE the grace -> kept
  never connected waits longer a page still loading is not a departure
  ending drops the record      whatever ends it, the landing page must
                               not offer it afterwards
  activity says nothing        a busy kernel nobody watches is ended;
                               an idle kernel somebody watches is not
"""

from __future__ import annotations

import asyncio
import logging

import pytest

from delfin.dashboard import resume_server as R
from delfin.dashboard import session as S


class _Base:
    """The parts of a MappingKernelManager the wrapper touches."""

    def __init__(self):
        self._kernel_connections = {}
        self._kernels = {}
        self.log = logging.getLogger("test")
        self.ended = []

    def __contains__(self, kid):
        return kid in self._kernels

    async def start_kernel(self, **kwargs):
        kid = kwargs.get("kernel_id") or f"k{len(self._kernels) + 1}"
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
        self._kernel_connections.pop(kid, None)

    async def shutdown_all(self, *a, **k):
        for kid in list(self._kernels):
            await self.shutdown_kernel(kid)


def _run(coro):
    # asyncio.run rather than get_event_loop().run_until_complete: in a
    # full suite some earlier test has called asyncio.run, which closes
    # and unsets the loop, and get_event_loop then raises. Ten of these
    # passed alone and failed in the run that matters.
    return asyncio.run(coro)


@pytest.fixture
def manager(tmp_path, monkeypatch):
    monkeypatch.setattr(S, "RECORD_DIR", str(tmp_path / "kept"))
    monkeypatch.delenv(R.GRACE_ENV, raising=False)
    cls = R.resume_kernel_manager_class(_Base)
    return cls()


def _age(manager, kid, seconds):
    """Move a kernel's unwatched-since back in time."""
    since, had = manager._delfin_unwatched()[kid]
    manager._delfin_unwatched()[kid] = (since - seconds, had)


def _keep(kid, name="kept-one", root=None):
    S.write_record(name, kid=kid, root=root)


# ---------------------------------------------------------------------------

def test_a_window_keeps_the_kernel(manager):
    kid = _run(manager.start_kernel())
    manager.notify_connect(kid)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == []


def test_no_window_past_the_grace_ends_it(manager):
    kid = _run(manager.start_kernel())
    manager.notify_connect(kid)
    manager.notify_disconnect(kid)
    _age(manager, kid, R.GRACE_SECONDS + 1)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == [kid]


def test_no_window_inside_the_grace_is_a_reload_not_a_departure(manager):
    kid = _run(manager.start_kernel())
    manager.notify_connect(kid)
    manager.notify_disconnect(kid)
    _age(manager, kid, R.GRACE_SECONDS - 5)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == []


def test_a_kept_session_is_not_ended(manager, tmp_path):
    kid = _run(manager.start_kernel())
    manager.notify_connect(kid)
    manager.notify_disconnect(kid)
    _age(manager, kid, R.GRACE_SECONDS * 100)
    _keep(kid)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == []


def test_disarming_ends_it_after_the_grace(manager):
    """Switch the option off on a page, close it: the next check after
    the grace ends the kernel. No beat to be stale here."""
    kid = _run(manager.start_kernel())
    manager.notify_connect(kid)
    _keep(kid)
    manager.notify_disconnect(kid)
    _age(manager, kid, R.GRACE_SECONDS + 1)
    S.drop_record("kept-one")
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == [kid]


def test_one_window_leaving_while_another_stays(manager):
    """Two tabs on one kernel: the first closing is not a departure."""
    kid = _run(manager.start_kernel())
    manager.notify_connect(kid)
    manager.notify_connect(kid)
    manager.notify_disconnect(kid)
    assert manager._delfin_seconds_unwatched(kid) is None
    manager.notify_disconnect(kid)
    assert manager._delfin_seconds_unwatched(kid) is not None


def test_a_page_still_loading_is_given_longer(manager):
    kid = _run(manager.start_kernel())
    _age(manager, kid, R.GRACE_SECONDS + 1)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == [], "never connected, inside the allowance"
    _age(manager, kid, R.NEVER_CONNECTED_SECONDS)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == [kid], "never connected, past the allowance"


def test_a_window_arriving_resets_the_clock(manager):
    kid = _run(manager.start_kernel())
    _age(manager, kid, R.NEVER_CONNECTED_SECONDS * 2)
    manager.notify_connect(kid)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == []


def test_ending_a_kernel_drops_its_record(manager, tmp_path):
    """When the kernel really ends -- here: the server stopping -- its
    record goes with it. A plain shutdown request on a kept kernel is the
    page's goodbye and is ignored (see the tests at the end)."""
    kid = _run(manager.start_kernel())
    _keep(kid, "gone-soon")
    assert S.list_records()
    _run(manager.shutdown_all())
    assert kid in manager.ended
    assert S.list_records() == []


def test_a_resumed_kernel_is_not_re_registered_as_newborn(manager, tmp_path, monkeypatch):
    """start_kernel handing back an existing kernel must not restart its
    unwatched clock from zero as if it had just been born."""
    kid = _run(manager.start_kernel())
    manager.notify_connect(kid)
    _keep(kid, "back")
    monkeypatch.setattr(R, "requested_session", lambda url: "back")
    again = _run(manager.start_kernel(env={"VOILA_REQUEST_URL": "http://h/?session=back"}))
    assert again == kid
    assert manager._delfin_seconds_unwatched(kid) is None, "a window is still there"


def test_the_grace_can_be_overridden_for_a_suite(monkeypatch):
    monkeypatch.setenv(R.GRACE_ENV, "7")
    assert R.grace_seconds() == 7.0
    monkeypatch.setenv(R.GRACE_ENV, "0")
    assert R.grace_seconds() == R.GRACE_SECONDS
    monkeypatch.setenv(R.GRACE_ENV, "nonsense")
    assert R.grace_seconds() == R.GRACE_SECONDS


def test_the_launcher_flags_switch_the_culler_on(monkeypatch):
    monkeypatch.delenv(R.GRACE_ENV, raising=False)
    args = R.cull_config_args()
    assert f"--MappingKernelManager.cull_idle_timeout={int(R.GRACE_SECONDS)}" in args
    assert f"--MappingKernelManager.cull_interval={R.POLL_SECONDS}" in args
    assert "--MappingKernelManager.cull_connected=False" in args


def test_the_launcher_passes_the_flags_with_the_manager():
    import inspect

    from delfin import cli_voila

    src = inspect.getsource(cli_voila.main)
    assert "_resume.cull_config_args()" in src


# ---------------------------------------------------------------------------
# The page's goodbye is not the owner's
# ---------------------------------------------------------------------------
#
# Voila's frontend sends a shutdown for its kernel when the page unloads.
# That is the very event a kept session is kept through. Seen on a
# cluster on 2026-09-11: "kernel ... is shutting down; dropping the
# record" seven seconds before the return request.

def test_a_kept_kernel_ignores_the_pages_goodbye(manager, tmp_path):
    kid = _run(manager.start_kernel())
    S.write_record("kept-1", root=str(tmp_path / "kept"), kid=kid)
    _run(manager.shutdown_kernel(kid))               # the beacon on unload
    assert kid not in manager.ended, "the page's goodbye ended a kept kernel"
    assert any(r["session_name"] == "kept-1" for r in S.list_records(root=str(tmp_path / "kept")))


def test_an_unkept_kernel_still_obeys_the_page(manager):
    kid = _run(manager.start_kernel())
    _run(manager.shutdown_kernel(kid))
    assert kid in manager.ended


def test_a_restart_goes_through_for_a_kept_kernel(manager, tmp_path):
    kid = _run(manager.start_kernel())
    S.write_record("kept-2", root=str(tmp_path / "kept"), kid=kid)
    _run(manager.shutdown_kernel(kid, False, True))    # (now=False, restart=True)
    assert kid in manager.ended


def test_the_server_stopping_ends_kept_kernels_and_drops_their_records(manager, tmp_path):
    kid = _run(manager.start_kernel())
    S.write_record("kept-3", root=str(tmp_path / "kept"), kid=kid)
    _run(manager.shutdown_all())
    assert kid in manager.ended
    assert not any(r["session_name"] == "kept-3" for r in S.list_records(root=str(tmp_path / "kept")))
