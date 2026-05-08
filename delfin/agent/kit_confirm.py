"""Self-Modification Guard broker for KIT-Toolbox coding-agent tools.

Since the per-action confirm dialog has been retired in favour of mode-
based control (mode chip in the dashboard) plus sandbox + denylist
enforcement, this broker now only handles the one remaining mandatory
case: the Self-Modification Guard. The guard fires whenever the agent
tries to write/edit a path that matches ``path_protected_globs`` (the
agent's own safety layer: api_client.py, kit_confirm.py, engine.py,
tab_agent.py, plus the remember_permission tool).

The agent stream runs in a worker thread; tool execution happens inline
on that worker. When the guard fires, the worker thread blocks while
the UI thread (Jupyter/ipywidgets) renders Approve/Deny buttons. This
module provides the threading bridge.

Usage in a dashboard tab::

    broker = KitConfirmBroker()
    panel  = broker.build_widget()        # ipywidgets container (slim)
    engine.set_kit_confirm_callback(broker.callback)

The mode dropdown / session-allowlist / counters / history panel from
earlier revisions has been removed — those were redundant with the KIT
mode chip and ``~/.delfin/settings.json`` and just added noise.
"""

from __future__ import annotations

import json
import threading
import time
from collections import deque
from dataclasses import dataclass, field
from typing import Any, Optional


@dataclass
class _ConfirmRequest:
    seq: int
    tool_name: str
    args: dict
    preview: str
    event: threading.Event = field(default_factory=threading.Event)
    decision: Optional[bool] = None
    remember_pattern: Optional[str] = None
    persist_pattern: Optional[str] = None  # regex form for ~/.delfin/settings.json


class KitConfirmBroker:
    """Thread-safe broker between the agent worker and the UI thread.

    The ``callback`` method is the entry point passed to
    ``engine.set_kit_confirm_callback``. UI rendering is done lazily via
    ``build_widget`` (returns an ipywidgets VBox).
    """

    def __init__(self, default_timeout_s: float = 300.0,
                 persist_callback: Optional[Any] = None) -> None:
        self._lock = threading.Lock()
        self._pending: deque[_ConfirmRequest] = deque()
        self._history: list[_ConfirmRequest] = []
        self._seq = 0
        self._timeout_s = default_timeout_s
        self._on_request: Optional[Any] = None  # UI callback to refresh the panel
        # Optional persistence hook: callable(kind, pattern) -> (ok, msg).
        # When a remember_permission request goes through the broker and the
        # user ticks "Dauerhaft speichern", the broker forwards the
        # suggested pattern here. The dashboard wires this to
        # engine.persist_kit_pattern so it survives across sessions.
        self._persist_callback: Optional[Any] = persist_callback
        # ipywidgets handles (lazy)
        self._panel: Any = None
        self._requests_box: Any = None
        self._toast: Any = None

    # -- public API --------------------------------------------------------

    def callback(self, tool_name: str, args: dict, preview: str) -> bool:
        """Called from the agent worker thread. Blocks until decided.

        With the per-action confirm dialog retired, this is now invoked
        only by the Self-Modification Guard (and ``remember_permission``,
        which uses the same threading bridge). No session-level allowlist
        is consulted — every call blocks for an explicit user decision.
        """
        with self._lock:
            self._seq += 1
            req = _ConfirmRequest(
                seq=self._seq,
                tool_name=tool_name,
                args=dict(args),
                preview=preview,
            )
            self._pending.append(req)

        # Notify UI (best-effort — UI may or may not be live).
        try:
            if self._on_request is not None:
                self._on_request()
        except Exception:
            pass
        self._refresh_panel()

        # Block worker thread until decided or timeout.
        decided = req.event.wait(timeout=self._timeout_s)

        with self._lock:
            try:
                self._pending.remove(req)
            except ValueError:
                pass
            self._history.append(req)
            if not decided and req.decision is None:
                req.decision = False  # timeout = deny

            persist_cb = self._persist_callback
            persist_pat = req.persist_pattern

        # remember_permission piggy-backs on this broker; if the user ticked
        # "Dauerhaft speichern" the persist callback writes to settings.json.
        # Self-Modification Guard requests never set persist_pattern (a
        # "always allow rewriting api_client.py" rule would defeat the guard).
        if persist_cb and persist_pat:
            kind = "allow" if req.decision else "deny"
            try:
                ok, msg = persist_cb(kind, persist_pat)
            except Exception as exc:
                ok, msg = False, f"persist failed: {exc}"
            self._set_toast(
                f"{'OK' if ok else 'FAIL'} dauerhaft: {msg}"
            )

        self._refresh_panel()
        return bool(req.decision)

    def set_persist_callback(self, callback) -> None:
        """Bind/unbind the persistence hook (callable(kind, pattern) -> (ok, msg))."""
        self._persist_callback = callback

    def snapshot(self) -> dict:
        """Inspect current state (for tests/CLI)."""
        with self._lock:
            return {
                "pending": [
                    {"seq": r.seq, "tool": r.tool_name, "args": r.args}
                    for r in self._pending
                ],
                "history": [
                    {"seq": r.seq, "tool": r.tool_name,
                     "decision": r.decision}
                    for r in self._history
                ],
            }

    # -- ipywidgets UI -----------------------------------------------------

    def build_widget(self):
        """Construct the Self-Modification-Guard panel (lazy import).

        The panel is collapsed (display:none) when no requests are pending
        so it doesn't clutter the chat — it pops into view only when the
        agent tries to edit a protected path and asks for approval.
        """
        try:
            import ipywidgets as widgets
        except Exception as exc:
            raise RuntimeError(f"ipywidgets is required for build_widget: {exc}")

        if self._panel is not None:
            return self._panel

        self._status_label = widgets.HTML(value="")
        self._requests_box = widgets.VBox([])
        self._toast = widgets.HTML(value="")

        self._panel = widgets.VBox(
            [
                widgets.HTML(
                    "<b>Self-Modification Guard</b> "
                    "<small>(forciert Confirm wenn der Agent seine eigene "
                    "Sicherheitsschicht ändern will — egal in welchem Mode)</small>"
                ),
                self._toast,
                self._requests_box,
            ],
            layout=widgets.Layout(
                border="1px solid #c5221f",
                padding="6px",
                margin="6px 0",
                display="none",
            ),
        )
        self._refresh_panel()
        return self._panel

    def _refresh_panel(self) -> None:
        if self._requests_box is None:
            return
        try:
            import ipywidgets as widgets  # noqa: F401
        except Exception:
            return

        rows = []
        with self._lock:
            pending = list(self._pending)

        for req in pending:
            rows.append(self._build_request_row(req, widgets))

        try:
            self._requests_box.children = tuple(rows)
            # Auto-collapse when nothing is pending — the panel only matters
            # while the agent is waiting on a Self-Mod-Guard or a
            # remember_permission decision.
            if self._panel is not None:
                self._panel.layout.display = "" if pending else "none"
        except Exception:
            pass

    def _build_request_row(self, req: _ConfirmRequest, widgets):
        title = f"#{req.seq}  {req.tool_name}"
        persist_pat = ""
        if req.tool_name == "bash":
            cmd = req.args.get("command", "")
            subtitle = f"<code>{_html_escape(cmd[:200])}</code>"
            suggested_pattern = _suggest_bash_pattern(cmd)
            persist_pat = _suggest_persist_pattern(cmd)
        elif req.tool_name in ("write_file", "edit_file", "multi_edit"):
            subtitle = f"<code>{_html_escape(req.args.get('path', ''))}</code>"
            suggested_pattern = req.args.get("path", "")
        else:
            subtitle = ""
            suggested_pattern = ""

        preview = widgets.HTML(
            value=(
                f"<b>{title}</b><br>{subtitle}"
                f"<details><summary>Vorschau</summary>"
                f"<pre style='max-height:240px; overflow:auto'>"
                f"{_html_escape(req.preview[:4000])}</pre></details>"
            )
        )
        approve = widgets.Button(description="Erlauben", button_style="success")
        deny = widgets.Button(description="Ablehnen", button_style="danger")
        remember = widgets.Checkbox(
            value=False,
            description=(f"Diese Session merken ({suggested_pattern[:40]})"
                         if suggested_pattern else "Diese Session merken"),
            indent=False,
        )
        persist_box = widgets.Checkbox(
            value=False,
            description=(f"Dauerhaft speichern ({persist_pat[:40]})"
                         if persist_pat else "Dauerhaft speichern"),
            indent=False,
            disabled=(not persist_pat or self._persist_callback is None),
            tooltip=(
                "Schreibt das Pattern nach ~/.delfin/settings.json — "
                "der Agent fragt für passende Befehle in zukünftigen Sessions "
                "nicht mehr nach."
                if persist_pat else
                "Persistenz nur für Bash-Befehle verfügbar."
            ),
        )

        def _decide(ok: bool):
            with self._lock:
                req.decision = ok
                if remember.value and suggested_pattern:
                    req.remember_pattern = suggested_pattern
                if persist_box.value and persist_pat:
                    req.persist_pattern = persist_pat
            req.event.set()
            self._refresh_panel()

        approve.on_click(lambda _b: _decide(True))
        deny.on_click(lambda _b: _decide(False))

        return widgets.VBox(
            [preview, widgets.HBox([approve, deny, remember, persist_box])],
            layout=widgets.Layout(
                border="1px solid #ccc",
                padding="4px",
                margin="2px 0",
            ),
        )

    def _set_toast(self, msg: str) -> None:
        """Show a one-line status message in the panel (best-effort)."""
        if self._toast is None:
            return
        try:
            self._toast.value = f"<small style='color:#444'>{_html_escape(msg)}</small>"
        except Exception:
            pass


# -- helpers ---------------------------------------------------------------

def _html_escape(s: str) -> str:
    return (
        s.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&#39;")
    )


def _suggest_bash_pattern(cmd: str) -> str:
    """Pick a sensible 'remember' substring (session allowlist).

    Examples:
      'pytest -xvs tests/foo.py'  -> 'pytest'
      'git commit -m "x"'         -> 'git commit'
      'python3 -m delfin.cli x'   -> 'python3 -m delfin.cli'
    """
    parts = cmd.strip().split()
    if not parts:
        return ""
    if parts[0] in ("git", "uv", "pip", "conda", "npm", "yarn", "cargo") and len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    if parts[0] in ("python", "python3") and len(parts) >= 3 and parts[1] == "-m":
        return f"{parts[0]} -m {parts[2]}"
    return parts[0]


def _suggest_persist_pattern(cmd: str) -> str:
    """Regex form for ``~/.delfin/settings.json`` allow_patterns.

    Delegates to ``kit_settings.suggest_pattern_for_command`` so the same
    rules apply that the engine consumes on load.
    """
    try:
        from .kit_settings import suggest_pattern_for_command
        return suggest_pattern_for_command(cmd)
    except Exception:
        return ""
