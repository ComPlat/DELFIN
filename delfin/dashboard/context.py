"""DashboardContext: shared state passed to every tab module."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, List, Optional

import ipywidgets as widgets
from IPython.display import clear_output, Javascript

from .backend_base import JobBackend


@dataclass
class DashboardContext:
    """Central shared state for the dashboard.

    Every tab's ``create_tab(ctx)`` function receives this object so that
    cross-tab communication and shared resources are available.
    """
    # Paths
    calc_dir: Path = field(default_factory=lambda: Path.home() / 'calc')
    archiv_dir: Path = field(default_factory=lambda: Path.home() / 'archiv')
    notebook_dir: Path = field(default_factory=Path.cwd)
    repo_dir: Optional[Path] = None

    # Job backend
    backend: Optional[JobBackend] = None

    # ORCA
    orca_base: Optional[str] = None
    orca_candidates: List[Path] = field(default_factory=list)

    # Only-GOAT template path (BwUni loads from file, local uses inline)
    only_goat_template_path: Optional[Path] = None

    # SLURM submit templates directory (BwUni only)
    submit_templates_dir: Optional[Path] = None

    # Run script for local backend
    run_script: Optional[Path] = None

    # Shared widgets (created by create_dashboard)
    js_output: widgets.Output = field(default_factory=widgets.Output)
    busy_indicator: widgets.HTML = field(default_factory=lambda: widgets.HTML(value=''))
    busy_css: Optional[widgets.HTML] = None

    # Cross-tab widget references (set by ORCA Builder, read by Calc Browser)
    orca_pal_widget: Any = None
    orca_maxcore_widget: Any = None

    # Templates
    default_control: str = ''
    only_goat_template: str = ''

    def run_js(self, script):
        """Execute JavaScript in a way that works in both Jupyter and Voila."""
        if not script:
            return
        with self.js_output:
            clear_output(wait=True)
            from IPython.display import display
            display(Javascript(script))

    def set_busy(self, is_busy):
        """Show or hide the busy spinner."""
        if is_busy:
            self.busy_indicator.value = (
                "<span class='delfin-busy' title='Working'></span>"
            )
        else:
            self.busy_indicator.value = ''
