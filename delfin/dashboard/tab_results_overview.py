"""Results Overview tab: pick calculations, pick quantities, get a table.

Pure logic (folder scan, value formatting, table building) lives in
:mod:`delfin.analysis.results_overview` and is testable without widgets.
This module only wires ipywidgets around it and registers itself
additively with :func:`delfin.dashboard.tab_registry.register_tab`, so
``create_dashboard`` stays untouched.
"""

from __future__ import annotations

import base64 as _base64
import csv as _csv
import html as _html
import io as _io
import json as _json
import os
import zipfile as _zipfile
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
        # Gallery controls: relative scan only, size in percent of the
        # 640px default, images per row before wrapping.
        self._png_size = widgets.IntSlider(
            value=70, min=20, max=150, step=10,
            description="Bildgroesse %", continuous_update=False,
            layout=widgets.Layout(width="260px"))
        self._png_cols = widgets.IntSlider(
            value=2, min=1, max=10, step=1,
            description="Bilder/Zeile", continuous_update=False,
            layout=widgets.Layout(width="260px"))
        self._png_opts = widgets.HBox([
            self._png_size, self._png_cols])
        self._csv_btn = widgets.Button(description="Export CSV")
        self._xlsx_btn = widgets.Button(description="Export XLSX")
        self._zip_btn = widgets.Button(description="Export PNG-ZIP")
        self._export_status = widgets.HTML()
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
        self._png_size.observe(self._on_png_opts_change, names="value")
        self._png_cols.observe(self._on_png_opts_change, names="value")
        self._csv_btn.on_click(self._on_export_csv)
        self._xlsx_btn.on_click(self._on_export_xlsx)
        self._zip_btn.on_click(self._on_export_zip)

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
            self._png_opts,
            self._status,
            self._sortbar,
            widgets.HBox([
                self._csv_btn, self._xlsx_btn, self._zip_btn,
                self._export_status]),
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

    def _on_png_opts_change(self, _change) -> None:
        if self._table is not None and self._png_mode == "inline":
            self._render()

    # -- export -------------------------------------------------------------

    def _table_bytes_csv(self) -> bytes:
        """CSV of the current table (format_value strings, Excel-friendly)."""
        assert self._table is not None
        buf = _io.StringIO()
        w = _csv.writer(buf, delimiter=";")
        w.writerow(["Rechnung"]
                   + [c["label"] for c in self._table["columns"]])
        for row in self._table["rows"]:
            w.writerow([row["label"]]
                       + [format_value(row["values"].get(c["id"]))
                          for c in self._table["columns"]])
        return buf.getvalue().encode("utf-8-sig")  # BOM: Excel opens UTF-8

    def _table_bytes_xlsx(self) -> bytes:
        """XLSX of the current table; missing values as empty cells."""
        import openpyxl  # already a DELFIN dependency
        assert self._table is not None
        wb = openpyxl.Workbook()
        ws = wb.active
        ws.title = "Results Overview"
        ws.append(["Rechnung"]
                  + [c["label"] for c in self._table["columns"]])
        for row in self._table["rows"]:
            cells = []
            for c in self._table["columns"]:
                v = row["values"].get(c["id"])
                cells.append(v if v is not None else None)
            ws.append([row["label"]] + cells)
        out = _io.BytesIO()
        wb.save(out)
        return out.getvalue()

    def _pngs_bytes_zip(self) -> bytes:
        """ZIP of every selected calc's relative scan PNG."""
        assert self._paths is not None
        out = _io.BytesIO()
        with _zipfile.ZipFile(out, "w", _zipfile.ZIP_DEFLATED) as zf:
            for path in self._paths:
                label = os.path.basename(path.rstrip("/"))
                for p in list_rss_pngs(path):
                    try:
                        with open(p["abs"], "rb") as fh:
                            data = fh.read(8 * 1024 * 1024)
                    except OSError:
                        continue
                    zf.writestr(f"{label}/{p['rel']}", data)
        return out.getvalue()

    def _export_bytes(self, filename: str, payload: bytes,
                      mime: str) -> None:
        """Download via the same base64-blob route the browser tab uses."""
        run_js = getattr(self._ctx, "run_js", None)
        if not callable(run_js):
            self._export_status.value = (
                "<span style='color:#dc3545'>Kein JS-Kontext: Export "
                "im Dashboard nicht moeglich.</span>")
            return
        b64 = _base64.b64encode(payload).decode("ascii")
        js = (
            "(function(){"
            f"const fileName={_json.dumps(filename)};"
            f"const mimeType={_json.dumps(mime)};"
            f"const b64={_json.dumps(b64)};"
            "try{"
            "const bin=atob(b64);"
            "const len=bin.length;"
            "const bytes=new Uint8Array(len);"
            "for(let i=0;i<len;i++){bytes[i]=bin.charCodeAt(i);}"
            "const blob=new Blob([bytes],{type:mimeType});"
            "const url=URL.createObjectURL(blob);"
            "const a=document.createElement('a');"
            "a.href=url;a.download=fileName;"
            "document.body.appendChild(a);a.click();"
            "setTimeout(function(){URL.revokeObjectURL(url);"
            "if(a.parentNode){a.parentNode.removeChild(a);}},1500);"
            "}catch(err){console.error('Download failed:', err);}"
            "})();"
        )
        run_js(js)
        self._export_status.value = (
            f"<span style='color:#28a745'>{_html.escape(filename)} "
            f"heruntergeladen ({len(payload) // 1024} KB).</span>")

    def _on_export_csv(self, _btn=None) -> None:
        if self._table is None:
            self._export_status.value = (
                "<span style='color:#dc3545'>Erst eine Tabelle "
                "erstellen.</span>")
            return
        self._export_bytes("results_overview.csv", self._table_bytes_csv(),
                           "text/csv;charset=utf-8")

    def _on_export_xlsx(self, _btn=None) -> None:
        if self._table is None:
            self._export_status.value = (
                "<span style='color:#dc3545'>Erst eine Tabelle "
                "erstellen.</span>")
            return
        self._export_bytes("results_overview.xlsx", self._table_bytes_xlsx(),
                           "application/vnd.openxmlformats-"
                           "officedocument.spreadsheetml.sheet")

    def _on_export_zip(self, _btn=None) -> None:
        if not self._paths:
            self._export_status.value = (
                "<span style='color:#dc3545'>Erst Rechnungen "
                "auswaehlen.</span>")
            return
        payload = self._pngs_bytes_zip()
        if not payload:
            self._export_status.value = (
                "<span style='color:#dc3545'>Keine RSS-PNGs "
                "gefunden.</span>")
            return
        self._export_bytes("rss_pngs.zip", payload, "application/zip")

    def _render(self) -> None:
        """Re-render table + PNG section from current state."""
        if self._table is None:
            return
        self._status.value = rows_to_markdownish(
            self._table["columns"], self._table["rows"], self._sort)
        # Gallery: ONE cross-calculation grid. Each calc contributes its
        # relative scan (per-calc chunking can never fill a row, since
        # one calc has only one relative PNG).
        png_html = []
        cells = []
        for row in self._table["rows"]:
            rels = list_rss_pngs(row["path"])
            if self._png_mode == "paths":
                if rels:
                    items = "".join(
                        f"<li><code>{_html.escape(p['rel'])}</code></li>"
                        for p in rels)
                    png_html.append(
                        f"<details><summary>RSS-PNGs: "
                        f"{_html.escape(row['label'])}"
                        f" ({len(rels)})</summary><ul>{items}</ul></details>")
                continue
            rel_only = [
                p for p in rels
                if p["rel"].endswith("scan_relative.png")
            ]
            shown = rel_only or rels[:1]
            if not shown:
                continue
            width = int(640 * self._png_size.value / 100)
            for p in shown:
                try:
                    with open(p["abs"], "rb") as fh:
                        data = fh.read(8 * 1024 * 1024)
                    b64 = _base64.b64encode(data).decode("ascii")
                    cells.append(
                        f"<td style='padding:6px;text-align:center;"
                        f"vertical-align:top'>"
                        f"<img src='data:image/png;base64,{b64}'"
                        f" alt='' style='width:{width}px;height:auto'>"
                        f"<br><small>{_html.escape(row['label'])}</small></td>")
                except OSError:
                    cells.append(
                        f"<td style='color:#dc3545'>{_html.escape(p['rel'])}"
                        f"</td>")
        if self._png_mode == "inline" and cells:
            cols = self._png_cols.value
            rows_png = []
            # Explicit rows: the dashboard panel can be narrower than
            # cols*width, so a max-width container would be clamped and
            # show fewer images per row than asked.
            for i in range(0, len(cells), cols):
                rows_png.append("<tr>" + "".join(cells[i:i + cols]) + "</tr>")
            png_html.append(
                "<table style='border-collapse:collapse'>"
                + "".join(rows_png) + "</table>")
        self._pngs.value = (
            "<h4>Relaxed-Surface-Scan PNGs</h4>" + "".join(png_html)
            if png_html else "")
        # Click-to-sort wiring: the column headers carry an onclick
        # attribute that calls window.__roSort(idx); that function (sent
        # once via ctx.add_init_js) clicks the matching hidden button.
        # HTML widgets cannot call Python directly, so the buttons stay
        # as the Python-side targets but are invisible (no extra row).
        sort_btns = [
            widgets.Button(description=f"⇅ {c['label'][:28]}", layout=widgets.Layout(
                display="none"))
            for c in self._table["columns"]
        ]
        for btn, col in zip(sort_btns, self._table["columns"]):
            btn.ro_sort_col = col["id"]
            # CSS class is the JS->Python bridge: the header link carries
            # data-col=<id>, the matching hidden button this class.
            btn.add_class(f"ro-sort-btn-{col['id']}")
            btn.on_click(self._on_sort)
        self._sortbar.children = tuple(sort_btns)
        self._sortbar.add_class("ro-sortbar")

    def init_js(self) -> str:
        """Startup JS: header clicks click the matching hidden button.

        A MutationObserver re-binds after every re-render of the HTML
        widget (each render replaces the table DOM).
        """
        return """
(function () {
  if (window.__roSortBind) { return; }
  window.__roSortBind = function (root) {
    (root || document).querySelectorAll('a.ro-sort').forEach(function (a) {
      if (a.__roBound) { return; }
      a.__roBound = true;
      a.addEventListener('click', function (ev) {
        ev.preventDefault();
        var col = a.getAttribute('data-col');
        var btn = document.querySelector('.ro-sortbar .ro-sort-btn-' + col);
        if (btn) { btn.click(); }
      });
    });
  };
  new MutationObserver(function (muts) {
    muts.forEach(function (m) {
      m.addedNodes.forEach(function (n) {
        if (n.nodeType === 1) { window.__roSortBind(n); }
      });
    });
  }).observe(document.body, {childList: true, subtree: true});
  window.__roSortBind();
})();
"""


def create_tab(ctx: Any):
    """Create the Results Overview tab. Returns ``(widget, refs)``."""
    panel = ResultsOverviewPanel(ctx)
    init = getattr(ctx, "add_init_js", None)
    if callable(init):
        init(panel.init_js())
    return panel.widget, {"results_overview": panel}


try:  # Register additively; create_dashboard stays untouched.
    from delfin.dashboard.tab_registry import register_tab
    register_tab("results_overview", "Results Overview", create_tab,
                 order=8100)
except Exception:  # pragma: no cover
    pass
