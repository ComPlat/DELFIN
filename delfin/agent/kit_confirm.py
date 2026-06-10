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
    persist_pattern: Optional[str] = None  # value for ~/.delfin/settings.json
    persist_kind: Optional[str] = None     # 'allow_pattern' or 'extra_dir'


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
            persist_kind = req.persist_kind

        # If the user clicked "Erlauben + Dauerhaft", forward the value
        # to the engine's persist hook so it lands in settings.json AND
        # applies to the live perms (so the next call doesn't re-prompt).
        # Persist only fires when the user APPROVED — a deny + persist
        # combination is meaningless for our supported kinds.
        if persist_cb and persist_pat and req.decision:
            if persist_kind == "extra_dir":
                cb_kind = "extra_dir"
            else:
                cb_kind = "allow"
            try:
                ok, msg = persist_cb(cb_kind, persist_pat)
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
        # ---- 1. Classify the request and compute what (if anything) can
        # be persisted. The user wants the dauerhaft-button to be
        # actionable when it makes sense, and CLEARLY explain why it's
        # disabled when it doesn't. -----------------------------------
        tool = req.tool_name
        cmd = req.args.get("command", "") if tool == "bash" else ""
        path_arg = req.args.get("path", "") if tool != "bash" else ""
        rationale = (req.args.get("rationale", "")
                     or req.args.get("description", "")).strip()

        persist_pat = ""
        persist_kind = ""        # 'allow_pattern' / 'extra_dir' / ''
        persist_disabled_reason = ""
        if tool == "bash":
            persist_pat = _suggest_persist_pattern(cmd)
            persist_kind = "allow_pattern"
            if not persist_pat:
                persist_disabled_reason = (
                    "Kein generalisierbares Pattern aus dem Bash-Befehl "
                    "ableitbar — kannst du im Chat mit "
                    "remember_permission(...) selbst formulieren."
                )
        elif tool in ("write_file", "edit_file", "multi_edit"):
            # Two distinct cases — show the RIGHT reason for each, so the
            # user isn't told "protected file" for a perfectly normal edit.
            if _is_protected_path(path_arg):
                # Self-Mod-Guard: protected core files must be approved every
                # time; persisting that would defeat the guard's purpose.
                persist_disabled_reason = (
                    "Self-Modification Guard: writes to protected core files "
                    "(api_client.py, kit_confirm.py, engine.py, tab_agent.py) "
                    "must be approved explicitly every time. Click 'Erlauben' "
                    "for THIS action."
                )
            else:
                # Normal file: there is no per-file persist rule for edits.
                # The way to stop being asked is the permission MODE — point
                # the user there instead of a misleading 'protected' message.
                persist_disabled_reason = (
                    "There is no per-file persist rule for writes/edits. To "
                    "stop being asked before every edit: switch the Perms "
                    "dropdown to 'acceptEdits' — write/edit actions then run "
                    "without prompting (sandbox + self-modification guard "
                    "stay active)."
                )
        elif tool == "read_file":
            # outside-workspace read: persisting = adding the parent
            # directory as extra_workspace_dir.
            from pathlib import Path
            try:
                parent = str(Path(path_arg).expanduser().resolve().parent)
                persist_pat = parent
                persist_kind = "extra_dir"
            except Exception:
                persist_disabled_reason = "Pfad nicht auflösbar."
        elif tool in ("remember_permission", "remember_permission_bundle"):
            # The click IS the persistence — no separate dauerhaft option.
            persist_disabled_reason = (
                "Diese Tool-Aktion IST die Persistenz: ein Klick auf "
                "'Erlauben' schreibt die Regel direkt in settings.json. "
                "Ablehnen verwirft sie."
            )

        # ---- 2. Header: tool name + path/cmd snippet + rationale ----
        if tool == "bash":
            subtitle = f'$&nbsp;<code>{_html_escape(cmd[:200])}</code>'
        elif path_arg:
            subtitle = f'<code>{_html_escape(path_arg)}</code>'
        else:
            subtitle = ""
        rationale_html = (
            f'<div style="font-size:11px; color:#9aa5b1; margin:2px 0;">'
            f'<i>Begründung Agent:</i> {_html_escape(rationale)}</div>'
            if rationale else ""
        )

        # ---- 3. Vorschau IMMER sichtbar (nicht in <details>) --------
        preview_pre = (
            f'<pre style="max-height:280px; overflow:auto; '
            f'font-size:11px; line-height:1.4; padding:6px; '
            f'background:rgba(127,127,127,0.06); '
            f'border-left:2px solid #94a3b8; margin:4px 0;">'
            f'{_html_escape(req.preview[:6000])}'
            f'{"…" if len(req.preview) > 6000 else ""}'
            f'</pre>'
        )

        header = widgets.HTML(
            value=(
                f'<b>#{req.seq} &middot; {tool}</b>'
                f'{("&nbsp;&middot;&nbsp;" + subtitle) if subtitle else ""}'
                f'{rationale_html}'
                f'{preview_pre}'
            )
        )

        # ---- 4. Dauerhaft-Status row: explicit target file --------
        target_path = ""
        try:
            from . import kit_settings as _ks
            target_path = str(_ks.USER_SETTINGS_PATH)
        except Exception:
            target_path = "~/.delfin/settings.json"
        persist_status: Any
        if persist_pat and self._persist_callback is not None:
            kind_label = (
                "Bash-Allow-Pattern" if persist_kind == "allow_pattern"
                else "extra_workspace_dir" if persist_kind == "extra_dir"
                else persist_kind
            )
            persist_status = widgets.HTML(value=(
                f'<div style="font-size:10px; color:#6b7280; margin:2px 0;">'
                f'Bei <b>Dauerhaft</b>: <code>{_html_escape(persist_pat)}</code> '
                f'({kind_label}) → <code>{_html_escape(target_path)}</code>'
                f'</div>'
            ))
        elif persist_disabled_reason:
            persist_status = widgets.HTML(value=(
                f'<div style="font-size:10px; color:#a16207; margin:2px 0;">'
                f'<b>Dauerhaft</b> nicht verfügbar: '
                f'{_html_escape(persist_disabled_reason)}'
                f'</div>'
            ))
        else:
            persist_status = widgets.HTML(value="")

        # ---- 5. Buttons --------------------------------------------
        approve = widgets.Button(
            description="Erlauben (1×)", button_style="success",
            tooltip="Diese eine Aktion durchführen, in zukünftigen "
                    "Sessions wieder fragen.",
        )
        approve_persist = widgets.Button(
            description="Erlauben + Dauerhaft",
            button_style="info",
            tooltip=(
                f"Aktion erlauben UND Regel nach {target_path} schreiben "
                "(gilt in zukünftigen Sessions ohne erneutes Fragen)."
            ),
            disabled=(not persist_pat or self._persist_callback is None),
        )
        deny = widgets.Button(description="Ablehnen", button_style="danger")

        def _decide(ok: bool, persist: bool = False):
            with self._lock:
                req.decision = ok
                if persist and persist_pat:
                    req.persist_pattern = persist_pat
                    req.persist_kind = persist_kind
            req.event.set()
            self._refresh_panel()

        approve.on_click(lambda _b: _decide(True, persist=False))
        approve_persist.on_click(lambda _b: _decide(True, persist=True))
        deny.on_click(lambda _b: _decide(False, persist=False))

        return widgets.VBox(
            [header,
             persist_status,
             widgets.HBox([approve, approve_persist, deny])],
            layout=widgets.Layout(
                border="1px solid #ccc",
                padding="6px",
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


def _is_protected_path(path_arg: str) -> bool:
    """True if ``path_arg`` is one of the self-mod-guarded core files.

    Used so the confirm dialog only shows the (intentional) 'cannot persist
    writes to protected files' reason for ACTUAL protected paths — normal
    files get the real guidance (switch to acceptEdits) instead.
    """
    if not path_arg:
        return False
    try:
        from .api_client import _DEFAULT_PATH_PROTECTED_GLOBS as _globs
    except Exception:
        _globs = ()
    import fnmatch
    p = str(path_arg).replace("\\", "/")
    for g in _globs:
        if p == g or p.endswith("/" + g) or fnmatch.fnmatch(p, "*/" + g) \
                or fnmatch.fnmatch(p, g):
            return True
    return False


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
