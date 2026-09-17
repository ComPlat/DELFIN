"""A page that closed said so; a connection that dropped did not.

The server sees one fact -- no websocket -- where the person at the
browser lived two: they closed the page, or their link went quiet. The
short grace was written for the first (it is really about a reload) and
was being applied to the second, so a sleeping laptop or a ping timeout
ended a session nobody had closed.

A page that is unloaded sends a beacon; silence is a drop. What is
pinned:

  a closed window            the beacon arrived -> the short grace
  a dropped connection       nothing arrived -> far longer
  the long one still ends    an unwatched session is not immortal
  a beacon is about ITS window  an older one does not shorten a later drop
  coming back clears it      a reconnect forgets the last closing
  the beacon can only shorten no route to a longer stay
  the page carries the script pagehide + sendBeacon, naming the kernel
  the script survives the widget  an HTML widget would strip it
"""

from __future__ import annotations

import asyncio
import logging

import pytest

from delfin.dashboard import resume_server as R
from delfin.dashboard import session as S
from delfin.dashboard import turn_record as T
from delfin.dashboard import window_close as W


class _Base:
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
    return asyncio.run(coro)


@pytest.fixture
def manager(tmp_path, monkeypatch):
    monkeypatch.setattr(S, "RECORD_DIR", str(tmp_path / "kept"))
    monkeypatch.setattr(T, "RECORD_DIR", str(tmp_path / "turns"))
    monkeypatch.delenv(R.GRACE_ENV, raising=False)
    monkeypatch.delenv(R.DROPPED_ENV, raising=False)
    cls = R.resume_kernel_manager_class(_Base)
    return cls()


def _leave(manager, kid):
    manager.notify_connect(kid)
    manager.notify_disconnect(kid)


def _age(manager, kid, seconds):
    since, had = manager._delfin_unwatched()[kid]
    manager._delfin_unwatched()[kid] = (since - seconds, had)
    for key, when in list(manager._delfin_closed().items()):
        if key == kid:
            manager._delfin_closed()[key] = when - seconds


# -- the two kinds of quiet -------------------------------------------------

def test_a_window_that_closed_gets_the_short_grace(manager):
    kid = _run(manager.start_kernel())
    _leave(manager, kid)
    manager.delfin_window_closed(kid)
    _age(manager, kid, R.GRACE_SECONDS + 1)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == [kid]


def test_a_connection_that_dropped_is_given_far_longer(manager):
    kid = _run(manager.start_kernel())
    _leave(manager, kid)
    # No beacon: the page is still open somewhere.
    _age(manager, kid, R.GRACE_SECONDS + 1)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == []


def test_even_a_dropped_connection_does_not_last_for_ever(manager):
    kid = _run(manager.start_kernel())
    _leave(manager, kid)
    _age(manager, kid, R.DROPPED_SECONDS + 1)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == [kid]


def test_a_beacon_from_an_earlier_window_does_not_shorten_a_later_drop(manager):
    """Two tabs. One is closed; much later the other one's link dies."""
    kid = _run(manager.start_kernel())
    manager.delfin_window_closed(kid)
    manager._delfin_closed()[kid] -= R.DROPPED_SECONDS  # an hour ago, say
    _leave(manager, kid)
    _age(manager, kid, R.GRACE_SECONDS + 1)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == []


def test_coming_back_forgets_the_last_closing(manager):
    kid = _run(manager.start_kernel())
    _leave(manager, kid)
    manager.delfin_window_closed(kid)
    manager.notify_connect(kid)
    assert kid not in manager._delfin_closed()

    manager.notify_disconnect(kid)
    _age(manager, kid, R.GRACE_SECONDS + 1)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == []


def test_a_page_still_loading_is_not_a_closed_window(manager):
    """A kernel that never had a window keeps its own allowance."""
    kid = _run(manager.start_kernel())
    manager._delfin_unwatched()[kid] = (0.0, False)
    manager.delfin_window_closed(kid)
    assert manager._delfin_allowance(kid, False) == R.NEVER_CONNECTED_SECONDS


def test_a_running_turn_outranks_a_closed_window(manager):
    """Closing the tab on a working agent does not throw the run away."""
    kid = _run(manager.start_kernel())
    T.mark(True, kid=kid)
    try:
        _leave(manager, kid)
        manager.delfin_window_closed(kid)
        _age(manager, kid, R.GRACE_SECONDS + 1)
        _run(manager.cull_kernel_if_idle(kid))
        assert manager.ended == []
    finally:
        T.mark(False, kid=kid)


def test_the_override_is_read_from_the_environment(manager, monkeypatch):
    monkeypatch.setenv(R.DROPPED_ENV, "30")
    assert R.dropped_seconds() == 30
    kid = _run(manager.start_kernel())
    _leave(manager, kid)
    _age(manager, kid, 31)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == [kid]


# -- the beacon itself ------------------------------------------------------

def test_the_script_names_the_kernel_and_sends_on_pagehide():
    script = W.beacon_script("kernel-42")
    assert "pagehide" in script
    assert "sendBeacon" in script
    assert '"kernel-42"' in script
    assert W.ROUTE in script


def test_the_script_is_not_put_where_a_widget_would_strip_it():
    """An HTML widget sanitises its value and would strip the script,
    leaving every close looking like a drop. An Output does not."""
    widgets = pytest.importorskip("ipywidgets")
    shown = []
    import IPython.display as _d

    real = _d.display
    _d.display = lambda *args, **kw: shown.extend(args)
    try:
        out = W.beacon_widget("kernel-42")
    finally:
        _d.display = real
    assert isinstance(out, widgets.Output)
    assert not isinstance(out, widgets.HTML)
    # Displayed INTO the Output, script intact. (An Output captures
    # nothing outside a kernel, so what was displayed is read here
    # rather than what the widget holds.)
    assert shown and "sendBeacon" in str(getattr(shown[0], "data", ""))
    assert "kernel-42" in str(getattr(shown[0], "data", ""))


def test_outside_a_kernel_there_is_no_beacon(monkeypatch):
    monkeypatch.setattr(S, "kernel_id", lambda: "")
    assert W.beacon_widget() is None


def test_the_route_is_added_to_the_server():
    added = []

    class _App:
        settings = {"base_url": "/prefix/"}

        def add_handlers(self, host, handlers):
            added.append((host, handlers))

    class _Server:
        web_app = _App()

    route = W.register(_Server())
    assert route == "/prefix/" + W.ROUTE
    assert added and added[0][1][0][0] == route


def test_a_server_that_cannot_take_the_route_is_not_a_failure():
    class _Server:
        pass

    assert W.register(_Server()) is None


def test_the_handler_tells_the_manager_which_window_closed():
    """The one thing the route does, and the only thing it can do."""
    from tornado.web import Application

    told = []

    class _Manager:
        def delfin_window_closed(self, kid):
            told.append(kid)

    handler_cls = W.handler_class()
    app = Application(kernel_manager=_Manager())
    handler = handler_cls.__new__(handler_cls)
    handler.application = app
    handler.get_argument = lambda name, default=None: (
        "kernel-42" if name == W.KERNEL_ARG else default)
    handler.set_status = lambda code: told.append(code)

    async def _finish():
        return None

    handler.finish = _finish
    # The decorated post is the authenticated wrapper; the body under it
    # is what runs once a request is authenticated.
    _run(handler_cls.post.__wrapped__(handler))
    assert told == ["kernel-42", 204]
