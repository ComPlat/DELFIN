"""Pipeline Status tab: live view of pipeline execution progress.

Shows a real-time table of pipeline steps, their status, and elapsed time.
Can be used standalone or embedded in the DELFIN Voila dashboard.

Standalone usage::

    from delfin.dashboard.tab_pipeline_status import PipelineStatusWidget

    widget = PipelineStatusWidget()
    display(widget)

    # Connect to a pipeline
    pipe = Pipeline("my_workflow")
    pipe.add("smiles_to_xyz", smiles="CCO")
    pipe.add("xtb_opt", charge=0)
    pipe.on_step(widget.on_step)  # live updates
    result = pipe.run(cores=4)
    widget.finalize(result)

Dashboard tab integration::

    def create_tab(ctx):
        widget = PipelineStatusWidget()
        return widget.widget, {"pipeline_status": widget}
"""

from __future__ import annotations

import time
from typing import Any, Dict, Optional

try:
    import ipywidgets as widgets
    HAS_WIDGETS = True
except ImportError:
    HAS_WIDGETS = False


STATUS_COLORS = {
    "success": "#28a745",
    "failed": "#dc3545",
    "skipped": "#6c757d",
    "running": "#007bff",
    "pending": "#ffc107",
}


class PipelineStatusWidget:
    """Live pipeline status widget for Jupyter/Voila dashboards.

    Usage::

        widget = PipelineStatusWidget()
        display(widget.widget)

        pipe.on_step(widget.on_step)
        result = pipe.run(cores=4)
        widget.finalize(result)
    """

    def __init__(self):
        if not HAS_WIDGETS:
            raise ImportError(
                "ipywidgets is required for the pipeline dashboard widget. "
                "Install with: pip install ipywidgets"
            )

        self._steps: list[dict] = []
        self._start_time = time.monotonic()

        # Header
        self._title = widgets.HTML(
            value="<h3>Pipeline Status</h3>",
        )

        # Progress bar
        self._progress = widgets.IntProgress(
            value=0, min=0, max=1,
            description="Progress:",
            bar_style="info",
            layout=widgets.Layout(width="100%"),
        )

        # Step table
        self._table = widgets.HTML(
            value=self._render_table(),
        )

        # Summary
        self._summary = widgets.HTML(value="")

        # Layout
        self.widget = widgets.VBox([
            self._title,
            self._progress,
            self._table,
            self._summary,
        ])

    def set_pipeline(self, name: str, n_steps: int) -> None:
        """Initialize the widget for a pipeline with known step count."""
        self._title.value = f"<h3>Pipeline: {name}</h3>"
        self._progress.max = n_steps
        self._progress.value = 0
        self._steps = []
        self._start_time = time.monotonic()
        self._summary.value = ""
        self._table.value = self._render_table()

    def on_step(self, result) -> None:
        """Callback for ``pipe.on_step()`` — updates widget after each step."""
        self._steps.append({
            "step_name": result.step_name,
            "status": result.status.value,
            "elapsed": result.elapsed_seconds,
            "error": result.error,
            "data_keys": list(result.data.keys()) if result.data else [],
        })
        self._progress.value = len(self._steps)
        if self._progress.max < len(self._steps):
            self._progress.max = len(self._steps)

        # Update bar style
        has_failure = any(s["status"] == "failed" for s in self._steps)
        self._progress.bar_style = "danger" if has_failure else "info"

        self._table.value = self._render_table()

    def finalize(self, pipeline_result) -> None:
        """Show final summary after pipeline completes."""
        self._progress.value = self._progress.max
        elapsed = time.monotonic() - self._start_time

        if pipeline_result.ok:
            self._progress.bar_style = "success"
            status_html = '<span style="color: #28a745; font-weight: bold;">OK</span>'
        else:
            self._progress.bar_style = "danger"
            status_html = '<span style="color: #dc3545; font-weight: bold;">FAILED</span>'

        self._summary.value = (
            f"<p><strong>Result:</strong> {status_html} | "
            f"<strong>Total time:</strong> {elapsed:.1f}s | "
            f"<strong>Steps:</strong> {len(pipeline_result.results)} | "
            f"<strong>Branches:</strong> {len(pipeline_result.branch_results)}</p>"
        )

    def _render_table(self) -> str:
        """Render the step status table as HTML."""
        if not self._steps:
            return "<em>No steps executed yet.</em>"

        rows = []
        for i, s in enumerate(self._steps):
            color = STATUS_COLORS.get(s["status"], "#333")
            status_badge = (
                f'<span style="color: {color}; font-weight: bold;">'
                f'{s["status"].upper()}</span>'
            )
            elapsed = f'{s["elapsed"]:.1f}s' if s["elapsed"] else ""
            error = f' <span style="color: #dc3545; font-size: 0.85em;">{s["error"]}</span>' if s["error"] else ""
            data_info = ", ".join(s["data_keys"][:5]) if s["data_keys"] else ""

            rows.append(
                f"<tr>"
                f"<td>{i + 1}</td>"
                f"<td>{s['step_name']}</td>"
                f"<td>{status_badge}{error}</td>"
                f"<td>{elapsed}</td>"
                f"<td style='font-size: 0.85em; color: #666;'>{data_info}</td>"
                f"</tr>"
            )

        return (
            '<table style="width: 100%; border-collapse: collapse; font-family: monospace;">'
            '<tr style="background: #f0f0f0; border-bottom: 2px solid #ccc;">'
            "<th>#</th><th>Step</th><th>Status</th><th>Time</th><th>Data</th></tr>"
            + "".join(rows)
            + "</table>"
        )


def create_tab(ctx):
    """Create the Pipeline Status tab for the DELFIN dashboard.

    Returns ``(tab_widget, refs_dict)``.
    """
    widget = PipelineStatusWidget()
    return widget.widget, {"pipeline_status": widget}
