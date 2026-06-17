"""Auto-generated "Run Application" tab.

Builds a form from any registered :class:`Application`'s input contract (the
``ParamSpec`` inputs), submits the run through the platform runtime, and shows
its status and named outputs.  This is the bridge between the tools platform and
the dashboard: a human gets a UI for the same applications an agent calls.

Registered via :func:`delfin.dashboard.tab_registry.register_tab`, so it is added
to the dashboard additively without touching ``create_dashboard``.

The form logic (input → field descriptors → collected values) is plain functions
with no ipywidgets dependency, so it is testable without a display.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

try:
    import ipywidgets as widgets
    HAS_WIDGETS = True
except ImportError:  # pragma: no cover
    HAS_WIDGETS = False


# --- pure form logic (no ipywidgets) --------------------------------------


def form_field_specs(inputs) -> List[Dict[str, Any]]:
    """Map an application's ``ParamSpec`` inputs to UI field descriptors."""
    fields: List[Dict[str, Any]] = []
    for p in inputs:
        if p.enum:
            kind = "enum"
        elif p.type == "bool":
            kind = "bool"
        elif p.type == "int":
            kind = "int"
        elif p.type == "float":
            kind = "float"
        else:
            kind = "text"
        fields.append({
            "name": p.name,
            "kind": kind,
            "required": p.required,
            "default": p.default,
            "options": list(p.enum) if p.enum else None,
            "description": p.description,
            "unit": p.unit,
        })
    return fields


def collect_inputs(
    fields: List[Dict[str, Any]], raw_values: Dict[str, Any],
) -> Tuple[Dict[str, Any], List[str]]:
    """Turn raw widget values into an inputs dict + a list of missing required.

    Empty optional text fields are dropped so the application/pipeline default
    applies; required empty text fields are reported as missing.
    """
    inputs: Dict[str, Any] = {}
    missing: List[str] = []
    for f in fields:
        value = raw_values.get(f["name"])
        if f["kind"] == "text":
            value = ("" if value is None else str(value)).strip()
            if value == "":
                if f["required"]:
                    missing.append(f["name"])
                continue
        elif value is None and f["required"]:
            missing.append(f["name"])
            continue
        inputs[f["name"]] = value
    return inputs, missing


# --- ipywidgets panel -----------------------------------------------------


def _make_field_widget(field: Dict[str, Any]):
    label = field["name"] + ("*" if field["required"] else "")
    default = field.get("default")
    tip = field.get("description", "") or ""
    layout = widgets.Layout(width="60%")
    style = {"description_width": "140px"}
    if field["kind"] == "enum":
        opts = field["options"] or []
        value = default if default in opts else (opts[0] if opts else None)
        return widgets.Dropdown(options=opts, value=value, description=label,
                                layout=layout, style=style)
    if field["kind"] == "bool":
        return widgets.Checkbox(value=bool(default), description=label, style=style)
    if field["kind"] == "int":
        return widgets.IntText(value=int(default) if default is not None else 0,
                               description=label, layout=layout, style=style)
    if field["kind"] == "float":
        return widgets.FloatText(value=float(default) if default is not None else 0.0,
                                 description=label, layout=layout, style=style)
    return widgets.Text(value="" if default is None else str(default),
                        placeholder=tip, description=label, layout=layout, style=style)


class ApplicationFormPanel:
    """Select an application, fill an auto-generated form, run it, see outputs."""

    def __init__(self):
        if not HAS_WIDGETS:
            raise ImportError("ipywidgets is required for the Run-Application panel.")

        from delfin.tools import platform

        self._platform = platform
        apps = platform.list_applications()

        self._title = widgets.HTML("<h3>Run Application</h3>")
        self._app_dd = widgets.Dropdown(
            options=apps or [], description="Application",
            style={"description_width": "140px"}, layout=widgets.Layout(width="60%"),
        )
        self._desc = widgets.HTML("")
        self._form_box = widgets.VBox([])
        self._run_btn = widgets.Button(description="Run", icon="play",
                                       button_style="primary")
        self._status = widgets.HTML("")

        self._fields: List[Dict[str, Any]] = []
        self._field_widgets: Dict[str, Any] = {}

        self._app_dd.observe(self._on_select, names="value")
        self._run_btn.on_click(self._on_run)

        if apps:
            self._rebuild_form(apps[0])
        else:
            self._desc.value = "<em>No applications registered.</em>"
            self._run_btn.disabled = True

        self.widget = widgets.VBox([
            self._title, self._app_dd, self._desc, self._form_box,
            self._run_btn, self._status,
        ])

    # -- form building -------------------------------------------------

    def _rebuild_form(self, app_name: str) -> None:
        app = self._platform.describe_application(app_name)
        if app is None:
            self._form_box.children = ()
            return
        outs = ", ".join(o.name for o in app.outputs) or "—"
        self._desc.value = f"<small>{app.description}<br>outputs: {outs}</small>"
        self._fields = form_field_specs(app.inputs)
        self._field_widgets = {f["name"]: _make_field_widget(f) for f in self._fields}
        self._form_box.children = tuple(self._field_widgets.values())

    def _on_select(self, change) -> None:
        if change.get("new"):
            self._rebuild_form(change["new"])

    # -- running -------------------------------------------------------

    def _on_run(self, _btn) -> None:
        import threading

        name = self._app_dd.value
        if not name:
            return
        raw = {n: w.value for n, w in self._field_widgets.items()}
        inputs, missing = collect_inputs(self._fields, raw)
        if missing:
            self._status.value = (
                f"<p style='color:#dc3545'>missing required input(s): "
                f"{', '.join(missing)}</p>"
            )
            return

        self._run_btn.disabled = True
        self._status.value = "<p>Submitting…</p>"

        def _work():
            try:
                run_id = self._platform.submit_application(name, **inputs)
                self._status.value = f"<p>Running… (run {run_id})</p>"
                rec = self._platform.wait_run(run_id, timeout=None)
                if rec is None:
                    self._status.value = "<p style='color:#dc3545'>run vanished.</p>"
                elif rec.status == "success":
                    outs = "".join(f"<li><b>{k}</b>: {v}</li>"
                                   for k, v in (rec.outputs or {}).items())
                    self._status.value = (
                        f"<p style='color:#28a745'>Done (run {rec.id}).</p>"
                        f"<ul>{outs or '<li>(no named outputs)</li>'}</ul>"
                    )
                else:
                    self._status.value = (
                        f"<p style='color:#dc3545'>{rec.status}: "
                        f"{rec.error or ''}</p>"
                    )
            except Exception as exc:  # noqa: BLE001
                self._status.value = f"<p style='color:#dc3545'>run failed: {exc}</p>"
            finally:
                self._run_btn.disabled = False

        threading.Thread(target=_work, daemon=True).start()


def create_tab(ctx: Any):
    """Create the Run-Application tab. Returns ``(widget, refs)``."""
    panel = ApplicationFormPanel()
    return panel.widget, {"application_form": panel}


# Register additively with the dashboard tab registry.
try:
    from delfin.dashboard.tab_registry import register_tab
    register_tab("run_application", "Run", create_tab, order=8500)
except Exception:  # pragma: no cover
    pass


__all__ = ["form_field_specs", "collect_inputs", "ApplicationFormPanel", "create_tab"]
