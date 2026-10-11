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

try:
    from jupyter_server_terminals.terminalmanager import TerminalManager
except ImportError:  # the tab imports this module; it must load without it
    TerminalManager = object  # type: ignore[assignment,misc]

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
        if TerminalManager is object:
            raise RuntimeError("jupyter_server_terminals is not installed")
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


# -- the browser side -----------------------------------------------------------
# xterm.js, pinned to exact versions and to their bytes (SRI): an unpinned
# CDN script in a page that holds a terminal is a supply-chain surface, and
# the browser refuses a file whose hash does not match.
XTERM_VERSION = "5.5.0"
XTERM_FIT_VERSION = "0.10.0"
XTERM_JS = ("https://cdn.jsdelivr.net/npm/@xterm/xterm@" + XTERM_VERSION
            + "/lib/xterm.js")
XTERM_JS_SRI = ("sha384-M169f14mRZOXm3hD/v2Ti0ThIT/RnAQagXA9nlE15yHAtrW19gdeP"
                "Jh/HaTzUOe/")
XTERM_CSS = ("https://cdn.jsdelivr.net/npm/@xterm/xterm@" + XTERM_VERSION
             + "/css/xterm.css")
XTERM_CSS_SRI = ("sha384-8Xk9wy/gzEDUKrXtrmCFa2bBuK3BpjpDuL/p0SeKQX19Khl/M+lH"
                 "OgD/CyYf7efP")
XTERM_FIT_JS = ("https://cdn.jsdelivr.net/npm/@xterm/addon-fit@"
                + XTERM_FIT_VERSION + "/lib/addon-fit.js")
XTERM_FIT_SRI = ("sha384-iF+jqbuti4XlB64clWgFWYEscb+UnSRv3VgVikGYZu+otNFnSHr7"
                 "y7NcKfBnGizn")

# Attaches every visible .delfin-claude-term to the server's terminal: the
# running one if there is one (a reload or a second window continues the
# same CLI session), else a new one. The server decides what runs; the page
# sends nothing but keystrokes and the window size.
_INIT_JS = r"""
(function () {
    if (window.__delfinClaudeTerm) return;
    window.__delfinClaudeTerm = true;
    var CFG = __DELFIN_CLAUDE_TERM_CFG__;
    var loading = null;

    function addScript(url, sri) {
        return new Promise(function (resolve, reject) {
            var s = document.createElement('script');
            s.src = url; s.integrity = sri; s.crossOrigin = 'anonymous';
            /* The UMD builds register with an AMD loader when one is
               present and then set no global; hide it for the load. */
            var amd = window.define;
            window.define = undefined;
            s.onload = function () { window.define = amd; resolve(); };
            s.onerror = function () { window.define = amd; reject(); };
            document.head.appendChild(s);
        });
    }

    function load() {
        if (loading) return loading;
        var css = document.createElement('link');
        css.rel = 'stylesheet'; css.href = CFG.css;
        css.integrity = CFG.cssSri; css.crossOrigin = 'anonymous';
        document.head.appendChild(css);
        loading = addScript(CFG.js, CFG.jsSri).then(function () {
            return addScript(CFG.fit, CFG.fitSri);
        });
        return loading;
    }

    function base() {
        try {
            var el = document.getElementById('jupyter-config-data');
            return (JSON.parse(el.textContent).baseUrl || '/');
        } catch (e) { return '/'; }
    }

    function xsrf() {
        var m = document.cookie.match('\\b_xsrf=([^;]*)\\b');
        return m ? decodeURIComponent(m[1]) : '';
    }

    function terminalName() {
        var b = base();
        return fetch(b + 'api/terminals', {credentials: 'same-origin'})
            .then(function (r) { return r.ok ? r.json() : []; })
            .then(function (list) {
                if (list && list.length) return list[0].name;
                return fetch(b + 'api/terminals', {
                    method: 'POST', credentials: 'same-origin',
                    headers: {'Content-Type': 'application/json',
                              'X-XSRFToken': xsrf()},
                    body: '{}'
                }).then(function (r) {
                    if (!r.ok) throw new Error('terminal ' + r.status);
                    return r.json();
                }).then(function (m) { return m.name; });
            });
    }

    function attach(el) {
        el.setAttribute('data-ready', '1');
        load().then(function () {
            var term = el.__term;
            if (!term) {
                term = new window.Terminal({
                    cursorBlink: true, fontSize: 13, scrollback: 5000,
                    fontFamily: 'ui-monospace, SFMono-Regular, Menlo, monospace',
                    theme: {background: '#1e1e1e'}
                });
                var fit = new window.FitAddon.FitAddon();
                term.loadAddon(fit);
                term.open(el);
                el.__term = term; el.__fit = fit;
                try {
                    new ResizeObserver(function () {
                        if (el.offsetParent === null) return;
                        try { fit.fit(); } catch (e) {}
                        send(el, ['set_size', term.rows, term.cols,
                                  el.clientHeight, el.clientWidth]);
                    }).observe(el);
                } catch (e) {}
                term.onData(function (d) { send(el, ['stdin', d]); });
            }
            try { el.__fit.fit(); } catch (e) {}
            return terminalName().then(function (name) { connect(el, name); });
        }).catch(function (e) {
            el.removeAttribute('data-ready');
            el.textContent = 'The Claude terminal could not start: '
                + (e && e.message ? e.message : 'xterm.js did not load');
        });
    }

    function send(el, msg) {
        var ws = el.__ws;
        if (ws && ws.readyState === 1) ws.send(JSON.stringify(msg));
    }

    function connect(el, name) {
        var proto = location.protocol === 'https:' ? 'wss://' : 'ws://';
        var ws = new WebSocket(proto + location.host + base()
                               + 'terminals/websocket/' + encodeURIComponent(name));
        el.__ws = ws;
        var term = el.__term;
        ws.onopen = function () {
            send(el, ['set_size', term.rows, term.cols,
                      el.clientHeight, el.clientWidth]);
            term.focus();
        };
        ws.onmessage = function (ev) {
            var msg;
            try { msg = JSON.parse(ev.data); } catch (e) { return; }
            if (msg[0] === 'stdout') term.write(msg[1]);
            else if (msg[0] === 'disconnect') ws.close();
        };
        ws.onclose = function () {
            if (el.__ws !== ws) return;
            el.__ws = null;
            term.write('\r\n\x1b[2m[Claude Code ended -- press Enter to '
                       + 'start a new session]\x1b[0m\r\n');
            var sub = term.onData(function (d) {
                if (d !== '\r') return;
                sub.dispose();
                term.clear();
                terminalName().then(function (n) { connect(el, n); });
            });
        };
    }

    function sweep() {
        var els = document.querySelectorAll('.delfin-claude-term:not([data-ready])');
        for (var i = 0; i < els.length; i++) {
            if (els[i].offsetParent !== null) attach(els[i]);
        }
    }
    sweep();
    try {
        new MutationObserver(sweep).observe(document.body, {
            childList: true, subtree: true, attributes: true,
            attributeFilter: ['style', 'class']});
    } catch (e) {}
})();
"""


def init_js() -> str:
    """The page script, with the pinned URLs and hashes filled in."""
    cfg = {"js": XTERM_JS, "jsSri": XTERM_JS_SRI, "css": XTERM_CSS,
           "cssSri": XTERM_CSS_SRI, "fit": XTERM_FIT_JS,
           "fitSri": XTERM_FIT_SRI}
    return _INIT_JS.replace("__DELFIN_CLAUDE_TERM_CFG__", json.dumps(cfg))


def enabled() -> bool:
    """True in a dashboard launched with --claude-terminal."""
    return os.environ.get(ENABLED_ENV) == "1"
