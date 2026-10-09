"""A dashboard terminal that runs the Claude Code CLI and nothing else.

Opt-in (``delfin-voila --claude-terminal``, loopback only). The stock
Jupyter terminal API cannot be restricted by configuration: ``POST
/api/terminals`` passes its JSON body to ``new_terminal`` as keyword
arguments, and ``new_terminal`` lets them override ``shell_command``,
``extra_env`` and ``cwd`` -- a client asking for ``{"shell_command":
["bash"]}`` gets a shell whatever the server was told to run. This
manager discards every client-supplied option except the window size, so
the process is always the configured CLI, in the configured directory,
with the server's environment minus the dashboard token.

Inputs, from the launcher's environment (``cli_voila``):

- ``DELFIN_CLAUDE_TERMINAL_ARGV``: JSON list, the command (absolute path
  of ``claude`` first).
- ``DELFIN_CLAUDE_TERMINAL_CWD``: the directory it starts in.

Rejected alternative: ``ServerApp.terminado_settings.shell_command`` alone
-- the default the body overrides (measured: see
tests/test_a_claude_terminal_runs_only_claude.py).
"""

from __future__ import annotations

import json
import os
from typing import Any

from jupyter_server_terminals.terminalmanager import TerminalManager

ARGV_ENV = "DELFIN_CLAUDE_TERMINAL_ARGV"
CWD_ENV = "DELFIN_CLAUDE_TERMINAL_CWD"
ENABLED_ENV = "DELFIN_CLAUDE_TERMINAL"
MAX_TERMINALS = 2

# Size fields the browser legitimately sends; everything else is dropped.
_SIZE_KEYS = ("height", "width", "winheight", "winwidth")
# Never handed to the CLI: the dashboard's own access token. Same user,
# but nothing the CLI or the commands it runs needs.
_DROP_ENV = ("JUPYTER_TOKEN", "JPY_API_TOKEN")


def configured_argv() -> list[str]:
    """The command from the launcher, or [] when absent or malformed."""
    try:
        argv = json.loads(os.environ.get(ARGV_ENV, "") or "[]")
    except ValueError:
        return []
    if (not isinstance(argv, list) or not argv
            or not all(isinstance(a, str) and a for a in argv)
            or not os.path.isabs(argv[0])):
        return []
    return argv


def configured_cwd() -> str:
    cwd = os.environ.get(CWD_ENV, "") or ""
    return cwd if cwd and os.path.isdir(cwd) else os.path.expanduser("~")


class ClaudeOnlyTerminalManager(TerminalManager):
    """Terminal manager whose terminals are always the configured CLI."""

    def __init__(self, **kwargs: Any) -> None:
        kwargs.pop("shell_command", None)
        kwargs["max_terminals"] = MAX_TERMINALS
        argv = configured_argv()
        if not argv:
            raise RuntimeError(
                f"{ARGV_ENV} is not set to an absolute command; the Claude "
                "terminal stays off rather than fall back to a shell.")
        super().__init__(shell_command=argv, **kwargs)
        self._delfin_argv = argv
        self._delfin_cwd = configured_cwd()

    def new_terminal(self, **kwargs: Any):
        safe: dict[str, Any] = {}
        for key in _SIZE_KEYS:
            value = kwargs.get(key)
            if isinstance(value, int) and 0 <= value <= 10_000:
                safe[key] = value
        return super().new_terminal(
            shell_command=list(self._delfin_argv),
            cwd=self._delfin_cwd, **safe)

    def new_named_terminal(self, **kwargs: Any):
        # terminado checks max_terminals only on the websocket route
        # (get_terminal); POST /api/terminals reaches this method directly
        # and was unbounded.
        if len(self.terminals) >= MAX_TERMINALS:
            from terminado.management import MaxTerminalsReached
            raise MaxTerminalsReached(MAX_TERMINALS)
        # The name is the only other field kept: a URL segment, never part
        # of the process. A taken name is not reused -- that would replace
        # a running terminal under the same route.
        name = kwargs.get("name")
        keep = {k: v for k, v in kwargs.items() if k in _SIZE_KEYS}
        if (isinstance(name, str) and name.isalnum() and len(name) <= 32
                and name not in self.terminals):
            keep["name"] = name
        return super().new_named_terminal(**keep)

    def make_term_env(self, *args: Any, **kwargs: Any) -> dict[str, str]:
        kwargs.pop("extra_env", None)
        env = super().make_term_env(*args, **kwargs)
        for key in _DROP_ENV:
            env.pop(key, None)
        return env
