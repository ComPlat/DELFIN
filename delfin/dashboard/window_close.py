"""Telling a closed window apart from a dropped connection.

The server ends a kernel whose last window has been gone for the grace.
From the server that is one fact -- the websocket count fell to zero --
but from the person at the browser it is two very different ones:

  they closed the page      the session is over; ending it frees the
                            port, which is the point of the short grace
  the connection dropped    a laptop slept, a VPN turned over, a ping
                            timed out. The page is still open and the
                            person is coming back to it

Only the first can say so. A page that is being unloaded gets one last
moment to send a beacon; a connection that is simply gone sends
nothing. So the short grace now belongs to the windows that announced
their own closing, and everything else -- every silence -- is read as a
drop and given far longer.

The beacon can only ever shorten a kernel's stay, never lengthen it,
and a lost beacon leaves the session running. Both directions of
failure end on the safe side.
"""

from __future__ import annotations

from typing import Any, Optional

#: The route the page's beacon calls, under the server's base url.
ROUTE = "delfin-api/window-closing"

#: The query key naming the kernel whose window closed.
KERNEL_ARG = "kernel"


def beacon_script(kernel_id: str) -> str:
    """The page's side of it: one beacon as the window goes.

    ``pagehide`` rather than ``unload``: it is the event a browser still
    fires for a page it puts in its back/forward cache, and the one
    ``sendBeacon`` is specified to survive. The url is derived from the
    page's own path because the server may be mounted under a prefix,
    and the beacon carries nothing but the kernel id.
    """
    import json

    kid = json.dumps(str(kernel_id))
    route = json.dumps(ROUTE)
    arg = json.dumps(KERNEL_ARG)
    return (
        "<script>(function(){\n"
        f"  var kid = {kid}, route = {route}, arg = {arg}, sent = false;\n"
        "  function base(){\n"
        "    var p = window.location.pathname || '/';\n"
        "    var i = p.indexOf('/voila/');\n"
        "    return i >= 0 ? p.slice(0, i + 1) : '/';\n"
        "  }\n"
        "  function bye(){\n"
        "    if (sent) { return; }\n"
        "    sent = true;\n"
        "    try {\n"
        "      navigator.sendBeacon(\n"
        "        base() + route + '?' + arg + '=' + encodeURIComponent(kid));\n"
        "    } catch (e) {}\n"
        "  }\n"
        "  window.addEventListener('pagehide', bye);\n"
        "})();</script>"
    )


def beacon_widget(kernel_id: str = ""):
    """A widget that carries the beacon into the page, or None.

    An Output rather than an HTML widget: ``widgets.HTML`` sanitises its
    value and strips the script, which would leave every close looking
    like a drop.
    """
    kid = kernel_id
    if not kid:
        from delfin.dashboard import session as _session
        try:
            kid = _session.kernel_id()
        except Exception:
            kid = ""
    if not kid:
        return None
    try:
        import ipywidgets as widgets
        from IPython.display import HTML, display
    except Exception:
        return None
    out = widgets.Output(layout=widgets.Layout(width="0px", height="0px"))
    with out:
        display(HTML(beacon_script(kid)))
    return out


def handler_class():
    """The server's side of it, built late so importing this module
    costs nothing outside a server."""
    from jupyter_server.base.handlers import JupyterHandler
    from tornado import web

    class WindowClosingHandler(JupyterHandler):
        """One page saying it is going away.

        Authenticated like every other api call. A beacon that cannot
        authenticate is simply lost, and a lost beacon reads as a drop,
        which keeps the session -- the safe way round.
        """

        def check_xsrf_cookie(self) -> None:
            # sendBeacon cannot carry the token header, and this accepts
            # no content: it names a kernel and shortens that kernel's
            # own grace. Nothing here can extend a stay or read a thing.
            return

        @web.authenticated
        async def post(self) -> None:
            kid = self.get_argument(KERNEL_ARG, "")
            manager = self.kernel_manager
            noted = getattr(manager, "delfin_window_closed", None)
            if kid and callable(noted):
                noted(kid)
            self.set_status(204)
            await self.finish()

    return WindowClosingHandler


def register(serverapp: Any) -> Optional[str]:
    """Add the route to a running server. Returns it, or None."""
    try:
        from jupyter_server.utils import url_path_join
    except Exception:
        return None
    try:
        web_app = serverapp.web_app
        base_url = web_app.settings.get("base_url", "/")
        route = url_path_join(base_url, ROUTE)
        web_app.add_handlers(".*$", [(route, handler_class())])
    except Exception:
        return None
    return route
