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
    office_dir: Path = field(default_factory=lambda: Path.home() / 'office')
    agent_dir: Path = field(default_factory=lambda: Path.home() / 'agent_workspace')
    primary_calc_dir: Optional[Path] = None
    default_calc_dir: Path = field(default_factory=lambda: Path.home() / 'calc')
    default_archive_dir: Path = field(default_factory=lambda: Path.home() / 'archive')
    default_office_dir: Path = field(default_factory=lambda: Path.home() / 'office')
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

    # Startup JavaScript, collected from the tabs that need it. A tab appends
    # its own script here; create_dashboard sends them as ONE run_js call,
    # because run_js clears the output first and a second call would wipe the
    # first script before the browser ran it. Collected rather than listed by
    # the caller: when the caller kept the list, a tab that was added later was
    # left out of it and its script silently never ran.
    init_js_parts: list = field(default_factory=list)

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
    remote_archive_refs: dict = field(default_factory=dict)

    # Templates
    default_control: str = ''
    only_goat_template: str = ''

    def add_init_js(self, script):
        """Register startup JavaScript for a tab.

        Call this from the tab that owns the script. ``create_dashboard``
        sends everything registered here in one go, so a tab cannot be
        forgotten by whoever assembles the dashboard.
        """
        if script and str(script).strip():
            self.init_js_parts.append(str(script))

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
