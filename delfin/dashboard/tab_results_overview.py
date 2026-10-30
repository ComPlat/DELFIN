"""Results Overview tab: pick calculations, pick quantities, get a table.

Pure logic (folder scan, value formatting, table building) lives in
:mod:`delfin.analysis.results_overview` and is testable without widgets.
This module only wires ipywidgets around it and registers itself
additively with :func:`delfin.dashboard.tab_registry.register_tab`, so
``create_dashboard`` stays untouched.
"""

from __future__ import annotations

import html as _html
import os
from typing import Any, List

try:
    import ipywidgets as widgets
    HAS_WIDGETS = True
except ImportError:  # pragma: no cover
    HAS_WIDGETS = False

from delfin.analysis.results_overview import (
    METRICS,
    build_table,
    format_value,
    list_calculations,
    list_rss_pngs,
    sort_table,
)


def rows_to_markdownish(columns: List[dict], rows: List[dict],
                        sorted_by: Any = None) -> str:
    """Build an HTML table for the status widget from parsed rows.

    ``sorted_by`` is ``(column_id, descending)``: the active sort column
    gets an arrow in its header and all headers are wrapped in clickable
    sort links (handled by :meth:`ResultsOverviewPanel._on_sort`).
    """
    sort_col, sort_desc = sorted_by if sorted_by else (None, False)
    head = []
    head.append("<th style='padding:4px 10px;text-align:left;border-bottom:"
                "2px solid #ccc'>Rechnung</th>")
    for col in columns:
        arrow = ""
        if col["id"] == sort_col:
            arrow = " &#9660;" if sort_desc else " &#9650;"
        head.append(
            f"<th style='padding:4px 10px;text-align:left;border-bottom:2px "
            f"solid #ccc'><a href='#' class='ro-sort' data-col="
            f"'{_html.escape(col['id'])}' style='color:inherit;text-"
            f"decoration:none'>{_html.escape(col['label'])}{arrow}</a></th>")
    head_html = "".join(head)
    body_rows = []
    for row in rows:
        cells = "".join(
            f"<td style='padding:4px 10px;border-bottom:1px solid #eee'>"
            f"{_html.escape(format_value(row['values'].get(col['id'])))}</td>"
            for col in columns
        )
        name = _html.escape(row["label"])
        body_rows.append(
            f"<tr><td style='padding:4px 10px;border-bottom:1px solid #eee'>"
            f"{name}</td>{cells}</tr>"
        )
    return (
        "<table style='border-collapse:collapse'>"
        f"<thead><tr><th style='padding:4px 10px;text-align:left;"
        f"border-bottom:2px solid #ccc'>Rechnung</th>{head_html}"
        f"</tr></thead>"
        f"<tbody>{''.join(body_rows)}</tbody></table>"
    )


class ResultsOverviewPanel:
    """Multi-select calculations x multi-select quantities -> table + PNGs."""

    def __init__(self, ctx: Any):
        self._ctx = ctx
        self._roots = self._resolve_roots(ctx)
        self._calcs = list_calculations(self._roots)

        labels = [c["label"] for c in self._calcs]
        paths = [c["path"] for c in self._calcs]
        self._calc_select = widgets.SelectMultiple(
            options=list(zip(labels, paths)),
            rows=12,
            layout=widgets.Layout(width="420px", height="300px"),
        )
        self._metric_select = widgets.SelectMultiple(
            options=[(m["label"], m["id"]) for m in METRICS],
            rows=len(METRICS),
            layout=widgets.Layout(width="420px"),
        )
        self._refresh_btn = widgets.Button(
            description="Tabelle erstellen",
            button_style="primary",
        )
        self._rescan_btn = widgets.Button(description="Rechnungen neu einlesen")
        self._show_pngs = widgets.Checkbox(
            value=False,
            description="PNGs inline anzeigen",
            indent=False,
        )
        self._status = widgets.HTML()
        self._pngs = widgets.HTML()
        self._sortbar = widgets.HBox()
        # Current table state for click-to-sort on the column headers.
        self._table: Any = None
        self._sort: Any = None  # (column_id, descending)
        self._paths: Any = None
        self._png_mode = "paths"
        self._refresh_btn.on_click(self._on_build)
        self._rescan_btn.on_click(self._on_rescan)
        self._show_pngs.observe(self._on_png_toggle, names="value")

        left = widgets.VBox([
            widgets.HTML("<b>Rechnungen</b> (Mehrfachauswahl)"),
            self._calc_select,
            self._rescan_btn,
        ])
        right = widgets.VBox([
            widgets.HTML("<b>Größen</b> (Mehrfachauswahl)"),
            self._metric_select,
            self._refresh_btn,
        ])
        self.widget = widgets.VBox([
            widgets.HBox([left, right]),
            self._show_pngs,
            self._status,
            self._sortbar,
            self._pngs,
        ])

    # -- helpers ------------------------------------------------------------

    @staticmethod
    def _resolve_roots(ctx: Any) -> List[str]:
        roots: List[str] = []
        for attr in ("archive_dir", "calc_dir"):
            path = getattr(ctx, attr, None)
            if path is None:
                continue
            path = str(path)
            if path and os.path.isdir(path) and path not in roots:
                roots.append(path)
        return roots or [os.path.expanduser("~/calc")]

    def _on_rescan(self, _btn=None) -> None:
        self._calcs = list_calculations(self._roots)
        labels = [c["label"] for c in self._calcs]
        paths = [c["path"] for c in self._calcs]
        self._calc_select.options = list(zip(labels, paths))
        self._status.value = (
            f"<p>{len(self._calcs)} Rechnungen gefunden unter "
            f"{_html.escape(', '.join(self._roots))}.</p>")

    def _on_build(self, _btn=None) -> None:
        paths = list(self._calc_select.value)
        mids = list(self._metric_select.value)
        if not paths or not mids:
            self._status.value = (
                "<p style='color:#dc3545'>Bitte Rechnungen und Größen "
                "auswählen.</p>")
            return
        try:
            table = build_table(paths, mids)
        except Exception as exc:  # noqa: BLE001
            self._status.value = (
                f"<p style='color:#dc3545'>Extraktion fehlgeschlagen: "
                f"{_html.escape(str(exc))}</p>")
            return
        self._table = table
        self._sort = None  # (column_id, descending)
        self._paths = paths
        self._render()

    def _on_sort(self, _btn) -> None:
        """Column-header click: sort by that column, toggle direction."""
        col = getattr(_btn, "ro_sort_col", None)
        if not col or self._table is None:
            return
        cur_col, cur_desc = self._sort or (None, False)
        if cur_col == col:
            desc = not cur_desc
        else:
            desc = False
        self._sort = (col, desc)
        self._table = sort_table(self._table, col, desc)
        self._render()

    def _on_png_toggle(self, change) -> None:
        self._png_mode = "inline" if change["new"] else "paths"
        if self._table is not None:
            self._render()

    def _render(self) -> None:
        """Re-render table + PNG section from current state."""
        if self._table is None:
            return
        self._status.value = rows_to_markdownish(
            self._table["columns"], self._table["rows"], self._sort)
        png_html = []
        for row in self._table["rows"]:
            rels = list_rss_pngs(row["path"])
            if not rels:
                continue
            if self._png_mode == "inline":
                # Only show the two standard scans inline, and only the
                # file that exists: large galleries make the tab unusable.
                preferred = [
                    "CO2_coordination/relaxed_surface_scan/scan_relative.png",
                    "CO2_coordination/relaxed_surface_scan/scan_absolute.png",
                ]
                shown = [p for p in rels if p["rel"] in preferred] or rels[:1]
                imgs = "".join(
                    f"<img src='file://{_html.escape(p['abs'], safe='')}'"
                    f" alt='' style='max-width:640px;display:block;"
                    f"margin:4px 0'>"
                    for p in shown
                )
                png_html.append(
                    f"<details open><summary>RSS-PNGs: "
                    f"{_html.escape(row['label'])}"
                    f" ({len(rels)})</summary>{imgs}</details>")
            else:
                items = "".join(
                    f"<li><code>{_html.escape(p['rel'])}</code></li>"
                    for p in rels
                )
                png_html.append(
                    f"<details><summary>RSS-PNGs: "
                    f"{_html.escape(row['label'])}"
                    f" ({len(rels)})</summary><ul>{items}</ul></details>")
        self._pngs.value = "<h4>Relaxed-Surface-Scan PNGs</h4>" + "".join(
            png_html) if png_html else ""
        # Click-to-sort wiring: one hidden button per sortable column.
        # HTML links in an HTML widget cannot call Python, so the panel
        # renders a button row under the table instead.
        sort_btns = [
            widgets.Button(description=f"⇅ {c['label'][:28]}")
            for c in self._table["columns"]
        ]
        for btn, col in zip(sort_btns, self._table["columns"]):
            btn.ro_sort_col = col["id"]
            btn.on_click(self._on_sort)
        self._sortbar.children = tuple(sort_btns)


def create_tab(ctx: Any):
    """Create the Results Overview tab. Returns ``(widget, refs)``."""
    panel = ResultsOverviewPanel(ctx)
    return panel.widget, {"results_overview": panel}


try:  # Register additively; create_dashboard stays untouched.
    from delfin.dashboard.tab_registry import register_tab
    register_tab("results_overview", "Results Overview", create_tab,
                 order=8100)
except Exception:  # pragma: no cover
    pass
