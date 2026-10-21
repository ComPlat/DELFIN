"""The dashboard server process cannot be read by the user's other processes.

The server is started by ``delfin-voila`` with the user's environment, which
may hold a provider key the user exported, and every kernel it starts runs
agent commands. A protection set before ``exec`` does not survive it, so
this server extension applies ``process_guard.protect`` inside the server
once it runs.

It also installs the HTTP-level protections of ``web_security`` (the
authenticated ``/voila/files/`` route, security headers, the cookie
policy), the browser-close beacon route and the preview proxy.
"""


def _jupyter_server_extension_points():
    return [{"module": "delfin.dashboard.server_guard"}]


def _load_jupyter_server_extension(serverapp) -> None:
    try:
        from delfin.agent import process_guard
        process_guard.protect("dashboard server")
    except Exception:
        pass
    # Authentication on Voila's file route, security headers on every
    # response, SameSite and a bounded lifetime for the login cookie.
    # Each part that fails is named in the log; the file route fails
    # closed (web_security.install_files_route).
    try:
        from delfin.dashboard import web_security
        web_security.install(serverapp)
    except Exception as exc:
        serverapp.log.warning(
            "[delfin] the dashboard's web protections were not installed: %s",
            exc)
    # The route a closing page sends its beacon to. Without it every
    # close and every dropped connection look the same, and a session
    # whose link went quiet for a minute would be ended as if the
    # person had closed it.
    try:
        from delfin.dashboard import window_close
        bridged = window_close.bridge_voila_shutdown()
        route = window_close.register(serverapp)
        if not bridged or not route:
            serverapp.log.warning(
                "[delfin] secure browser-close handling could not be fully "
                "installed; an unkept dashboard may take longer to stop.")
    except Exception as exc:
        # The default Voila route still has Jupyter's XSRF protection.  This
        # failure therefore loses cleanup convenience, never authentication;
        # make the degradation visible instead of silently weakening either.
        serverapp.log.warning(
            "[delfin] secure browser-close handling was not installed: %s", exc)


    # Previews served through this server's own port: the one port VS Code
    # forwards, or the Windows launcher tunnels, reaches them too.
    try:
        from delfin.dashboard import preview_proxy
        if preview_proxy.register(serverapp) is None:
            serverapp.log.warning(
                "[delfin] the preview proxy could not be installed; previews "
                "are linked by their own port")
    except Exception as exc:
        serverapp.log.warning("[delfin] the preview proxy was not installed: %s", exc)


async def _start_jupyter_server_extension(serverapp) -> None:
    """Runs once every extension is loaded: Voila's file route exists by now
    whichever order the extensions were loaded in."""
    try:
        from delfin.dashboard import web_security
        web_security.install_files_route(serverapp.web_app, serverapp.log)
    except Exception as exc:
        serverapp.log.warning(
            "[delfin] /voila/files/ could not be put behind authentication: %s",
            exc)


# Older jupyter_server spelling.
load_jupyter_server_extension = _load_jupyter_server_extension
