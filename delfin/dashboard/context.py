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
    archive_dir: Path = field(default_factory=lambda: Path.home() / 'archive')
    agent_dir: Path = field(default_factory=lambda: Path.home() / 'agent_workspace')
    primary_calc_dir: Optional[Path] = None
    default_calc_dir: Path = field(default_factory=lambda: Path.home() / 'calc')
    default_archive_dir: Path = field(default_factory=lambda: Path.home() / 'archive')
    runtime_settings: dict = field(default_factory=dict)
    runtime_backend: str = 'auto'
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

    # Shared widgets (created by create_dashboard)
    js_output: widgets.Output = field(default_factory=widgets.Output)
    busy_indicator: widgets.HTML = field(default_factory=lambda: widgets.HTML(value=''))
    busy_css: Optional[widgets.HTML] = None

    # Cross-tab widget references (set by ORCA Builder, read by Calc Browser)
    orca_pal_widget: Any = None
    orca_maxcore_widget: Any = None

    # Shared clipboard for cross-tab cut/copy/paste
    clipboard_paths: list = field(default_factory=list)
    clipboard_mode: str = ''  # 'cut' or 'copy'
    shared_clipboard: dict = field(
        default_factory=lambda: {'paths': [], 'mode': ''}
    )
    shared_explorer_state: dict = field(
        default_factory=lambda: {'refresh_hooks': {}, 'refresh_running': False}
    )
    tabs_widget: Any = None
    tab_indices: dict = field(default_factory=dict)
    tab_specs: list = field(default_factory=list)
    rebuild_dashboard_tabs: Any = None

    # Agent engine (set by tab_agent, used cross-tab)
    agent_engine: Any = None
    # Agent status indicator (shown in top header bar, updated by tab_agent)
    agent_status_html: widgets.HTML = field(
        default_factory=lambda: widgets.HTML(value='')
    )

    # Cross-tab refs for agent dashboard control
    submit_refs: dict = field(default_factory=dict)
    orca_builder_refs: dict = field(default_factory=dict)
    job_status_refs: dict = field(default_factory=dict)
    recalc_refs: dict = field(default_factory=dict)
    calc_browser_refs: dict = field(default_factory=dict)

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

    def select_tab(self, title):
        """Switch to a dashboard tab by title when available."""
        if not self.tabs_widget:
            return False
        try:
            index = int(self.tab_indices.get(title))
        except Exception:
            return False
        self.tabs_widget.selected_index = index
        return True
