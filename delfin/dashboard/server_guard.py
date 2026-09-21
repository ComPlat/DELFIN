"""The dashboard server process cannot be read by the user's other processes.

The server is started by ``delfin-voila`` with the user's environment, which
may hold a provider key the user exported, and every kernel it starts runs
agent commands. A protection set before ``exec`` does not survive it, so
this server extension applies ``process_guard.protect`` inside the server
once it runs.
"""


def _jupyter_server_extension_points():
    return [{"module": "delfin.dashboard.server_guard"}]


def _load_jupyter_server_extension(serverapp) -> None:
    try:
        from delfin.agent import process_guard
        process_guard.protect("dashboard server")
    except Exception:
        pass
    # The route a closing page sends its beacon to. Without it every
    # close and every dropped connection look the same, and a session
    # whose link went quiet for a minute would be ended as if the
    # person had closed it.
    try:
        from delfin.dashboard import window_close
        window_close.register(serverapp)
    except Exception:
        pass


# Older jupyter_server spelling.
load_jupyter_server_extension = _load_jupyter_server_extension
