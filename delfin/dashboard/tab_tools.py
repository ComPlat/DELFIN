"""Tools & Platform tab — surface the DELFIN tools platform in the dashboard.

Lets a human browse the same platform an agent binds to: capabilities (building
blocks), applications (callable workflows), the central key vocabulary, run
history, and — importantly — tool installation status.  Open-source tools can be
installed with one click (the bundled installer); license-restricted tools
(ORCA, Turbomole) are never installed — their official source and where to place
the files are shown instead.

The HTML renderers are plain functions (no ipywidgets), so they are testable
without a display; ``create_tab`` wires them into an ipywidgets panel exactly
like the other dashboard tabs.
"""

from __future__ import annotations

from typing import Any

try:
    import ipywidgets as widgets
    HAS_WIDGETS = True
except ImportError:  # pragma: no cover - exercised only without ipywidgets
    HAS_WIDGETS = False

_TABLE = ('<table style="width:100%;border-collapse:collapse;font-family:monospace;'
          'font-size:0.9em">')
_TH = '<tr style="background:#f0f0f0;border-bottom:2px solid #ccc">'


def _safe(fn):
    """Run a renderer, turning any error into a visible message instead of a crash."""
    try:
        return fn()
    except Exception as exc:  # noqa: BLE001
        return f"<p style='color:#dc3545'>unavailable: {exc}</p>"


# --- pure HTML renderers --------------------------------------------------


def render_environment_html() -> str:
    """Tool readiness + install plan, with manual-install guidance for ORCA etc."""
    from delfin.tools import platform

    readiness = platform.probe()
    plan = platform.install_plan()
    ready = sum(1 for r in readiness if r.ready)
    total = len(readiness)

    parts = [f"<h4>Tool environment — {ready}/{total} capabilities ready</h4>"]

    auto = plan.get("auto_binaries", [])
    if auto:
        parts.append(
            "<p><b>Auto-installable (open source):</b> "
            f"{', '.join(auto)} — use <i>Install open-source tools</i> above.</p>"
        )

    manual = plan.get("manual", [])
    if manual:
        items = "".join(
            f"<li><b>{m['tool']}</b>: {m['hint']}"
            + (f" &middot; <a href='{m['source']}' target='_blank'>{m['source']}</a>"
               if m.get("source") else "")
            + "</li>"
            for m in manual
        )
        parts.append(
            "<p><b>Manual / license-restricted "
            "(DELFIN does <u>not</u> install these):</b></p>"
            f"<ul>{items}</ul>"
        )

    python = plan.get("python", [])
    if python:
        items = "".join(f"<li>{p['tool']}: <code>{p['hint']}</code></li>" for p in python)
        parts.append(f"<p><b>Python packages:</b></p><ul>{items}</ul>")

    if ready == total:
        parts.append("<p style='color:#28a745'>All capability requirements satisfied.</p>")
    return "".join(parts)


def render_tool_status_html() -> str:
    """Per-tool install status — every required tool with ✓ installed / ✗ missing."""
    from delfin.tools import platform

    inv = platform.tool_inventory()
    ready = sum(1 for t in inv if t["available"])
    parts = [f"<h4>Tool status — {ready}/{len(inv)} installed</h4>", _TABLE,
             _TH + "<th>tool</th><th>kind</th><th>status</th>"
                   "<th>how to get it</th><th>used by</th></tr>"]
    for t in inv:
        if t["available"]:
            status = "<span style='color:#28a745'>✓ installed</span>"
            how = f"<span style='font-size:0.8em;color:#666'>{t['detail']}</span>"
        else:
            status = "<span style='color:#dc3545'>✗ missing</span>"
            if t["policy"] == "auto":
                how = "open source — tick it in 'Selective install' above"
            else:
                src = (f" &middot; <a href='{t['source']}' target='_blank'>{t['source']}</a>"
                       if t["source"] else "")
                how = f"{t['install_hint']}{src}"
        used = ", ".join(t["used_by"][:4]) + (" …" if len(t["used_by"]) > 4 else "")
        parts.append(
            f"<tr><td><b>{t['name']}</b></td><td>{t['kind']}</td><td>{status}</td>"
            f"<td style='font-size:0.85em'>{how}</td>"
            f"<td style='font-size:0.8em;color:#666'>{used}</td></tr>"
        )
    parts.append("</table>")
    return "".join(parts)


def render_capabilities_html() -> str:
    from delfin.tools import platform

    cat = platform.catalog(by="category")
    parts = ["<h4>Capabilities (building blocks)</h4>"]
    for category in sorted(cat):
        names = ", ".join(sorted(c.name for c in cat[category]))
        parts.append(f"<p><b>{category}</b>: {names}</p>")
    return "".join(parts)


def render_applications_html() -> str:
    from delfin.tools import platform

    apps = platform.list_applications()
    if not apps:
        return "<p><em>No applications registered.</em></p>"
    parts = ["<h4>Applications (callable workflows)</h4>"]
    for name in apps:
        app = platform.describe_application(name)
        ins = ", ".join(p.name + ("*" if p.required else "") for p in app.inputs)
        outs = ", ".join(o.name for o in app.outputs)
        parts.append(
            f"<p><b>{name}</b> — {app.description}<br>"
            f"<small>inputs: {ins or '—'} | outputs: {outs or '—'}</small></p>"
        )
    return "".join(parts)


def render_keys_html() -> str:
    from delfin.tools import platform

    parts = ["<h4>Key vocabulary (define once, reuse everywhere)</h4>", _TABLE,
             _TH + "<th>key</th><th>type</th><th>default</th><th>#allowed</th></tr>"]
    for name in platform.list_keys():
        ks = platform.describe_key(name)
        n_allowed = len(ks.enum) if ks.enum else ""
        parts.append(
            f"<tr><td>{ks.name}</td><td>{ks.type}</td>"
            f"<td>{ks.default if ks.default is not None else ''}</td>"
            f"<td>{n_allowed}</td></tr>"
        )
    parts.append("</table>")
    return "".join(parts)


def _fmt_outputs(data: dict) -> str:
    if not data:
        return ""
    items = list(data.items())[:4]
    text = ", ".join(f"{k}={v}" for k, v in items)
    return text + (" …" if len(data) > 4 else "")


def render_runs_html(limit: int = 20) -> str:
    from delfin.tools import platform

    runs = platform.list_runs()
    metrics = platform.run_metrics()
    parts = [
        f"<h4>Runs — {metrics['total_runs']} total</h4>",
        f"<p>by status: {metrics.get('by_status', {})}</p>",
        "<p><small>Named outputs are stored per run in "
        "<code>~/.delfin/runs/&lt;id&gt;.json</code>; raw files (ORCA .out, "
        "geometries, logs) land in <code>~/calc/&lt;app&gt;_&lt;id&gt;/</code> "
        "(the standard calc dir → Calculations browser).</small></p>",
    ]
    if runs:
        rows = "".join(
            f"<tr><td>{r.id}</td><td>{r.name}</td><td>{r.status}</td>"
            f"<td>{_fmt_outputs(r.outputs)}</td>"
            f"<td style='font-size:0.8em;color:#666'>{r.work_dir or ''}</td></tr>"
            for r in runs[:limit]
        )
        parts.append(_TABLE + _TH
                     + "<th>id</th><th>application</th><th>status</th>"
                       "<th>outputs</th><th>work dir (files)</th></tr>"
                     + rows + "</table>")
    else:
        parts.append("<p><em>No runs yet.</em></p>")
    return "".join(parts)


# --- run detail (the calc browser) ----------------------------------------


def _human_size(n: int) -> str:
    size = float(n)
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024:
            return f"{size:.0f}{unit}"
        size /= 1024
    return f"{size:.0f}TB"


def list_run_files(work_dir, limit: int = 200) -> list:
    """All files under a run's work dir as ``{path, size, bytes}`` (path relative)."""
    from pathlib import Path

    out: list = []
    if not work_dir:
        return out
    wd = Path(work_dir)
    if not wd.is_dir():
        return out
    for f in sorted(wd.rglob("*")):
        if len(out) >= limit:
            break
        if f.is_file():
            try:
                size = f.stat().st_size
            except OSError:
                size = 0
            out.append({"path": str(f.relative_to(wd)), "size": _human_size(size),
                        "bytes": size})
    return out


def tail_text(path, n: int = 40) -> str:
    """Last *n* lines of a text file (best-effort)."""
    from pathlib import Path

    try:
        lines = Path(path).read_text(errors="replace").splitlines()
    except OSError as exc:
        return f"(cannot read {path}: {exc})"
    return "\n".join(lines[-n:])


def render_run_detail_html(run_id: str, *, tail_lines: int = 40) -> str:
    """Detailed view of one run — status, error, outputs, events, the files in its
    ``~/calc`` work dir, and a tail of its main log. The visual face of
    ``platform.run_diagnostics``."""
    import html as _html

    from delfin.tools import platform

    if not run_id:
        return "<p><em>Select a run to see its details.</em></p>"
    diag = platform.run_diagnostics(run_id)
    if "status" not in diag:
        return f"<p style='color:#dc3545'>{diag.get('error', 'unknown run')}</p>"

    colour = {"success": "#28a745", "failed": "#dc3545", "cancelled": "#6c757d"}.get(
        diag["status"], "#0366d6")
    parts = [f"<h4>Run {diag['id']} — {diag['name']} "
             f"<span style='color:{colour}'>[{diag['status']}]</span></h4>"]
    if diag.get("error"):
        parts.append(f"<p style='color:#dc3545'>error: {_html.escape(str(diag['error']))}</p>")
    if diag.get("outputs"):
        parts.append(f"<p>outputs: <code>{_html.escape(_fmt_outputs(diag['outputs']))}</code></p>")
    parts.append(f"<p><small>work dir: <code>{diag.get('work_dir') or '—'}</code></small></p>")

    events = diag.get("events") or []
    if events:
        body = _html.escape("\n".join(str(e) for e in events))
        parts.append("<b>events</b><pre style='background:#f6f8fa;padding:6px;"
                     f"max-height:160px;overflow:auto'>{body}</pre>")

    files = list_run_files(diag.get("work_dir"))
    if files:
        rows = "".join(f"<tr><td>{_html.escape(f['path'])}</td>"
                       f"<td style='text-align:right'>{f['size']}</td></tr>" for f in files)
        parts.append("<b>files (~/calc)</b>" + _TABLE + _TH
                     + "<th>file</th><th>size</th></tr>" + rows + "</table>")
    else:
        parts.append("<p><em>No files in the work dir yet.</em></p>")

    logs = diag.get("log_files") or []
    if logs:
        main = max(logs, key=_file_size)
        body = _html.escape(tail_text(main, n=tail_lines))
        parts.append(f"<b>log tail — <code>{_html.escape(main)}</code></b>"
                     "<pre style='background:#0d1117;color:#c9d1d9;padding:6px;"
                     f"max-height:240px;overflow:auto'>{body}</pre>")
    return "".join(parts)


def _file_size(path) -> int:
    from pathlib import Path
    try:
        return Path(path).stat().st_size
    except OSError:
        return 0


def render_all_html() -> str:
    """One combined HTML view of the whole platform (used by the widget)."""
    sections = [render_environment_html, render_tool_status_html,
                render_capabilities_html, render_applications_html,
                render_keys_html, render_runs_html]
    return "<hr>".join(_safe(s) for s in sections)


# --- ipywidgets panel -----------------------------------------------------


class ToolsPanel:
    """Dashboard panel for the tools platform (status, install, browse, runs)."""

    def __init__(self):
        if not HAS_WIDGETS:
            raise ImportError("ipywidgets is required for the Tools panel.")

        self._title = widgets.HTML("<h3>Tools &amp; Platform</h3>")
        self._refresh_btn = widgets.Button(description="Refresh", icon="refresh")
        self._install_btn = widgets.Button(
            description="Install all open-source", icon="download",
            tooltip="Install all auto-installable open-source tools "
                    "(ORCA / Turbomole / Multiwfn are never installed).",
        )
        self._install_sel_btn = widgets.Button(
            description="Install selected", icon="download", button_style="primary",
        )
        self._install_box = widgets.VBox([])
        self._install_checks: dict = {}
        self._status = widgets.HTML("")
        self._body = widgets.HTML(render_all_html())

        # run detail (the calc browser): pick a run → status/files/log tail
        self._run_select = widgets.Dropdown(description="Run:", options=[], value=None)
        self._run_detail = widgets.HTML(render_run_detail_html(""))
        self._run_select.observe(self._on_run_select, names="value")

        self._refresh_btn.on_click(self._on_refresh)
        self._install_btn.on_click(lambda _b: self._run_install(None))
        self._install_sel_btn.on_click(self._on_install_selected)

        self._build_install_checklist()
        self._refresh_run_options()

        self.widget = widgets.VBox([
            self._title,
            widgets.HBox([self._refresh_btn, self._install_btn]),
            widgets.HTML("<b>Selective install</b> — tick what to install:"),
            self._install_box,
            self._status,
            self._body,
            widgets.HTML("<hr><b>Run detail</b> — inspect a run's files in "
                         "<code>~/calc</code> and its log:"),
            self._run_select,
            self._run_detail,
        ])

    def _refresh_run_options(self) -> None:
        from delfin.tools import platform
        try:
            runs = platform.list_runs()
        except Exception:  # noqa: BLE001
            runs = []
        self._run_select.options = [(f"{r.id} — {r.name} [{r.status}]", r.id) for r in runs]

    def _on_run_select(self, _change) -> None:
        self._run_detail.value = _safe(
            lambda: render_run_detail_html(self._run_select.value or ""))

    def _build_install_checklist(self) -> None:
        from delfin.tools import platform
        try:
            plan = platform.install_plan()
        except Exception:  # noqa: BLE001
            self._install_box.children = ()
            return
        self._install_checks = {}
        rows = []
        auto = plan.get("auto_binaries", [])
        if auto:
            cb = widgets.Checkbox(value=True, indent=False,
                                  description=f"QM tools bundle ({', '.join(auto)})")
            self._install_checks["__qm__"] = (cb, list(auto))
            rows.append(cb)
        for p in plan.get("python", []):
            cb = widgets.Checkbox(value=False, indent=False,
                                  description=f"pip: {p['tool']}")
            self._install_checks[p["tool"]] = (cb, [p["tool"]])
            rows.append(cb)
        if rows:
            rows.append(self._install_sel_btn)
            self._install_box.children = tuple(rows)
        else:
            self._install_box.children = (
                widgets.HTML("<em>Nothing auto-installable is missing.</em>"),)

    def refresh(self) -> None:
        self._body.value = render_all_html()
        self._build_install_checklist()
        self._refresh_run_options()
        self._run_detail.value = _safe(
            lambda: render_run_detail_html(self._run_select.value or ""))

    def _on_refresh(self, _btn) -> None:
        self.refresh()

    def _on_install_selected(self, _btn) -> None:
        select = []
        for cb, names in self._install_checks.values():
            if cb.value:
                select.extend(names)
        if not select:
            self._status.value = "<p>Nothing selected.</p>"
            return
        self._run_install(select)

    def _run_install(self, select) -> None:
        import threading

        self._install_btn.disabled = True
        self._install_sel_btn.disabled = True
        what = "selected tools" if select else "all open-source tools"
        self._status.value = (f"<p>Installing {what}… "
                              "(ORCA/Turbomole/Multiwfn are skipped)</p>")

        def _work():
            from delfin.tools import platform
            try:
                result = platform.install_tools(run=True, select=select)
                if not result.get("executed"):
                    self._status.value = ("<p>Nothing was installed (nothing "
                                          "selected/missing, or installer not found).</p>")
                elif result.get("ok"):
                    self._status.value = "<p style='color:#28a745'>Installed.</p>"
                else:
                    self._status.value = ("<p style='color:#dc3545'>Some installs "
                                          "failed — check the run logs.</p>")
            except Exception as exc:  # noqa: BLE001
                self._status.value = f"<p style='color:#dc3545'>install failed: {exc}</p>"
            finally:
                self._install_btn.disabled = False
                self._install_sel_btn.disabled = False
                self.refresh()

        threading.Thread(target=_work, daemon=True).start()


def create_tab(ctx: Any):
    """Create the Tools & Platform tab. Returns ``(widget, refs)``."""
    panel = ToolsPanel()
    return panel.widget, {"tools_panel": panel}


__all__ = [
    "render_environment_html",
    "render_tool_status_html",
    "render_capabilities_html",
    "render_applications_html",
    "render_keys_html",
    "render_runs_html",
    "render_run_detail_html",
    "list_run_files",
    "tail_text",
    "render_all_html",
    "ToolsPanel",
    "create_tab",
]
