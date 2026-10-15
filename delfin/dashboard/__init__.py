"""DELFIN Dashboard: modular, backend-agnostic dashboard for Jupyter / Voila.

Usage::

    from delfin.dashboard import create_dashboard
    ctx = create_dashboard(backend='auto')
"""

import base64
import contextlib
import fcntl
import importlib
import importlib.resources
import html
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from datetime import datetime
import warnings
from pathlib import Path

# Suppress noisy but harmless third-party warnings before they are imported
os.environ.setdefault("TORCHANI_NO_WARN_EXTENSIONS", "1")
warnings.filterwarnings("ignore", message="PySisiphus is not installed", module="aimnet2calc")

# These imports are safe without ipywidgets and needed by module-level helpers
from delfin.runtime_setup import (
    apply_runtime_environment,
    describe_installation_source,
    detect_local_runtime_limits,
    discover_orca_installations,
    get_packaged_submit_templates_dir,
    resolve_backend_choice,
    resolve_orca_base,
    resolve_submit_templates_dir,
)
from delfin.dashboard import session as _session
from delfin.quota import home_usage
from delfin.user_settings import load_remote_archive_enabled, load_settings


def create_dashboard(backend='auto', calc_dir=None, orca_base=None):
    """Create and display the full DELFIN Dashboard.

    Parameters
    ----------
    backend : str
        ``'auto'`` (default), ``'local'``, or ``'slurm'``.
    calc_dir : str or Path, optional
        Calculations directory.  Defaults to ``~/calc``.
    orca_base : str, optional
        Path to the ORCA installation directory.

    Returns
    -------
    DashboardContext
        The context object (useful for programmatic access / debugging).
    """
    # -- lazy imports (keep dashboard importable without ipywidgets) --------
    import ipywidgets as widgets
    from IPython.display import display

    from .constants import DEFAULT_CONTROL, ONLY_GOAT_TEMPLATE
    from .context import DashboardContext
    from .helpers import (
        _append_js,
        apply_branding,
        create_busy_css,
        create_page_css,
        debounce_input,
        disable_spellcheck_global,
        hand_over_tabs,
        keep_scroll_position,
    )
    from .molecule_viewer import RIGHT_MOUSE_TRANSLATE_PATCH_JS, vendored_3dmol_js

    from . import (
        tab_agent,
        tab_agent_activity,
        tab_archive_statistics,
        tab_home,
        tab_office,
        tab_calculations_browser,
        tab_literature,
        tab_remote_archive,
        tab_settings,
        tab_job_status,
        tab_orca_builder,
        tab_recalc,
        tab_submit,
        tab_tools,
    )

    # -- ensure delfin is importable & up-to-date --------------------------
    _ensure_delfin_importable()
    _force_reload_delfin()

    # -- resolve paths -----------------------------------------------------
    home = Path.home()
    notebook_dir = _get_notebook_dir()
    settings_payload = {}
    configured_paths = {}
    runtime_settings = {}
    try:
        settings_payload = load_settings()
        configured_paths = (settings_payload.get('paths', {}) or {})
        runtime_settings = (settings_payload.get('runtime', {}) or {})
    except Exception:
        settings_payload = {}
        configured_paths = {}
        runtime_settings = {}

    # SLURM: look for a root_dir/calc pattern; local: ~/calc
    root_dir = _find_root_dir(notebook_dir)
    if root_dir and (root_dir / 'calc').exists():
        default_calc_dir = root_dir / 'calc'
    else:
        default_calc_dir = home / 'calc'
    default_archive_dir = default_calc_dir.parent / 'archive'
    # Deliberately NOT derived from the launch directory like calc and
    # archive: those belong to a project checkout, office documents do not.
    # The agent keys its memory and its saved sessions by this folder, so a
    # launch-dependent default would split both in two the moment the
    # dashboard is started somewhere else — which reads as lost memory
    # rather than as a second folder. An explicit Office path in the
    # settings still wins.
    default_office_dir = home / 'office'

    if calc_dir is None:
        calc_dir = configured_paths.get('calculations_dir') or default_calc_dir
    calc_dir = Path(calc_dir).expanduser()
    calc_dir.mkdir(parents=True, exist_ok=True)
    archive_dir = configured_paths.get('archive_dir') or default_archive_dir
    archive_dir = Path(archive_dir).expanduser()
    archive_dir.mkdir(parents=True, exist_ok=True)
    office_dir = configured_paths.get('office_dir') or default_office_dir
    office_dir = Path(office_dir).expanduser()
    office_dir.mkdir(parents=True, exist_ok=True)
    agent_dir = Path(configured_paths.get('agent_dir') or Path.home() / 'agent_workspace')
    agent_dir.mkdir(parents=True, exist_ok=True)

    repo_dir = _find_delfin_root()

    # -- auto-detect backend -----------------------------------------------
    backend = resolve_backend_choice(
        backend,
        runtime_settings.get('backend', 'auto'),
        slurm_available=shutil.which('sbatch') is not None,
    )

    # -- build backend object ----------------------------------------------
    if backend == 'slurm':
        from .backend_slurm import SlurmJobBackend

        submit_templates_dir = resolve_submit_templates_dir(
            runtime_settings,
            get_packaged_submit_templates_dir(),
        )
        orca_candidates = _find_orca_candidates(notebook_dir)
        auto_orca_candidates = [str(path) for path in orca_candidates]
        auto_orca_candidates.sort(
            key=lambda value: ('6_1_1' not in Path(value).name, value)
        )
        orca_base = resolve_orca_base(
            orca_base,
            runtime_settings,
            'slurm',
            auto_candidates=auto_orca_candidates,
            local_default='',
        )
        apply_runtime_environment(
            qm_tools_root=runtime_settings.get('qm_tools_root', ''),
            orca_base=orca_base,
            csp_tools_root=runtime_settings.get('csp_tools_root', ''),
            mlp_tools_root=runtime_settings.get('mlp_tools_root', ''),
            tool_binaries=runtime_settings.get('tool_binaries', {}) or {},
        )

        slurm_profile = (runtime_settings.get('slurm', {}) or {}).get('profile', '')
        backend_obj = SlurmJobBackend(
            submit_templates_dir=submit_templates_dir,
            orca_base=orca_base,
            tool_binaries=runtime_settings.get('tool_binaries', {}) or {},
            slurm_profile=slurm_profile,
            partitions=(runtime_settings.get('slurm', {}) or {}).get('partitions', ''),
            gpu_partitions=(runtime_settings.get('slurm', {}) or {}).get('gpu_partitions', ''),
        )
        only_goat_template_path = _find_only_goat_template(notebook_dir)
    else:
        from .backend_local import LocalJobBackend

        detected_local_cores, detected_local_ram_mb = detect_local_runtime_limits()
        orca_base = resolve_orca_base(
            orca_base,
            runtime_settings,
            'local',
            auto_candidates=[],
            local_default='',
        )
        orca_candidates = []
        apply_runtime_environment(
            qm_tools_root=runtime_settings.get('qm_tools_root', ''),
            orca_base=orca_base,
            csp_tools_root=runtime_settings.get('csp_tools_root', ''),
            mlp_tools_root=runtime_settings.get('mlp_tools_root', ''),
            tool_binaries=runtime_settings.get('tool_binaries', {}) or {},
        )
        backend_obj = LocalJobBackend(
            orca_base=orca_base,
            tool_binaries=runtime_settings.get('tool_binaries', {}) or {},
            max_cores=int((runtime_settings.get('local', {}) or {}).get('max_cores', detected_local_cores)),
            max_ram_mb=int((runtime_settings.get('local', {}) or {}).get('max_ram_mb', detected_local_ram_mb)),
            allow_oversubscribe=bool((runtime_settings.get('local', {}) or {}).get('allow_oversubscribe', False)),
            oversubscribe_factor=float((runtime_settings.get('local', {}) or {}).get('oversubscribe_factor', 1.0)),
            allow_live_load_bypass=bool((runtime_settings.get('local', {}) or {}).get('allow_live_load_bypass', False)),
            live_cpu_target_factor=float((runtime_settings.get('local', {}) or {}).get('live_cpu_target_factor', 0.95)),
            live_min_free_ram_mb=int((runtime_settings.get('local', {}) or {}).get('live_min_free_ram_mb', 64_000)),
        )
        only_goat_template_path = None

    # -- shared widgets ----------------------------------------------------
    js_output = widgets.Output()
    busy_indicator = widgets.HTML(value='')
    busy_css = create_busy_css()

    # -- build context -----------------------------------------------------
    ctx = DashboardContext(
        calc_dir=calc_dir,
        archive_dir=archive_dir,
        office_dir=office_dir,
        agent_dir=agent_dir,
        primary_calc_dir=calc_dir,
        default_calc_dir=default_calc_dir,
        default_archive_dir=default_archive_dir,
        default_office_dir=default_office_dir,
        runtime_settings=runtime_settings,
        runtime_backend=backend,
        notebook_dir=notebook_dir,
        repo_dir=repo_dir,
        backend=backend_obj,
        orca_base=orca_base,
        orca_candidates=[str(p) for p in orca_candidates] if orca_candidates else [],
        only_goat_template_path=only_goat_template_path,
        submit_templates_dir=getattr(backend_obj, 'submit_templates_dir', None),
        js_output=js_output,
        busy_indicator=busy_indicator,
        busy_css=busy_css,
        default_control=DEFAULT_CONTROL,
        only_goat_template=ONLY_GOAT_TEMPLATE,
    )

    # -- create tabs -------------------------------------------------------
    tab1, refs1 = tab_submit.create_tab(ctx)
    tab2, refs2 = tab_recalc.create_tab(ctx)
    tab3, refs3 = tab_job_status.create_tab(ctx)
    tab4, refs4 = tab_orca_builder.create_tab(ctx)

    # Wire cross-tab widget references for agent dashboard control
    ctx.orca_pal_widget = refs4.get('orca_pal')
    ctx.orca_maxcore_widget = refs4.get('orca_maxcore')
    ctx.submit_refs = refs1
    ctx.orca_builder_refs = refs4
    ctx.job_status_refs = refs3
    ctx.recalc_refs = refs2

    # TURBOMOLE tab (only for SLURM backends)
    tab_tm = None
    if ctx.backend.supports_turbomole:
        from . import tab_turbomole_builder
        tab_tm, _ = tab_turbomole_builder.create_tab(ctx)

    try:
        remote_archive_enabled = load_remote_archive_enabled()
    except Exception:
        remote_archive_enabled = False

    tab5, refs5 = tab_calculations_browser.create_tab(ctx)
    ctx.calc_browser_refs = refs5
    # A tab hidden in Settings is built anyway today: the preference drops it
    # from the bar and nothing else. Archive and Office are two more builds of
    # the calculations browser -- the heaviest thing a cold start does -- so
    # hiding one and paying for it is the worst of both. They are built when
    # they are wanted, which for a visible tab is still at startup.
    try:
        _hidden_at_start = {
            str(item) for item in
            (((load_settings().get('ui', {}) or {}).get('tabs', {}) or {})
             .get('hidden', []) or [])
        }
    except Exception:
        _hidden_at_start = set()

    refs6, refs_off, refs_home = {}, {}, {}

    def _clone_of_browser_now(build, refs):
        """Build one of the browser clones and fill in the refs Settings holds."""
        widget, made = build(ctx)
        refs.update(made or {})
        return widget

    def _clone_of_browser(tab_id, build, refs):
        """The tab, or None if Settings has it hidden.

        *refs* is the dict the Settings tab is handed further down. It is
        filled in rather than replaced, because Settings captures it by value
        and applying a new Archive or Office folder calls back through it.
        """
        if tab_id in _hidden_at_start:
            return None
        return _clone_of_browser_now(build, refs)

    tab6 = _clone_of_browser('archive', tab_archive_statistics.create_tab, refs6)
    tab_off = _clone_of_browser('office', tab_office.create_tab, refs_off)
    tab_hm = _clone_of_browser('home', tab_home.create_tab, refs_home)
    tab7, refs7 = (tab_remote_archive.create_tab(ctx) if remote_archive_enabled else (None, {}))
    ctx.remote_archive_refs = refs7
    tab_lit, _ = tab_literature.create_tab(ctx)
    # The agent tab is a list of sessions, each a whole agent tab. Should the
    # list itself fail to build, the single tab it wraps still works.
    try:
        from . import agent_sessions
        tab_ag, _ = agent_sessions.create_tab(ctx)
    except Exception:
        tab_ag, _ = tab_agent.create_tab(ctx)
    tab_ag_act = tab_agent_activity.create_tab(ctx)
    # Tools & Platform tab (defensive: never let it break dashboard startup)
    try:
        tab_tools_widget, refs_tools = tab_tools.create_tab(ctx)
    except Exception:
        tab_tools_widget, refs_tools = None, {}
    ctx.tools_refs = refs_tools
    try:
        from delfin.agent.credentials import load_credential as _lc_anth
        _anth_key = _lc_anth("ANTHROPIC_API_KEY")
    except Exception:
        _anth_key = os.environ.get("ANTHROPIC_API_KEY", "")
    _agent_backend_available = bool(shutil.which("claude") or _anth_key)

    tab_specs = [
        {
            'id': 'submit',
            'title': 'Submit Job',
            'widget': tab1,
            'default_order': 10,
            'default_visible': True,
            'available': True,
            'fixed': False,
            'reason': '',
        },
        {
            'id': 'recalc',
            'title': 'Recalc',
            'widget': tab2,
            'default_order': 20,
            'default_visible': True,
            'available': True,
            'fixed': False,
            'reason': '',
        },
        {
            'id': 'orca_builder',
            'title': 'ORCA Builder',
            'widget': tab4,
            'default_order': 30,
            'default_visible': True,
            'available': True,
            'fixed': False,
            'reason': '',
        },
        {
            'id': 'turbomole_builder',
            'title': 'TURBOMOLE Builder',
            'widget': tab_tm,
            'default_order': 40,
            'default_visible': tab_tm is not None,
            'available': tab_tm is not None,
            'fixed': False,
            'reason': '' if tab_tm is not None else 'Only available for supported SLURM backends.',
        },
        {
            'id': 'chemdarwin',
            'title': 'ChemDarwin',
            'widget': None,
            'default_order': 50,
            'default_visible': False,
            'available': False,
            'fixed': False,
            'reason': 'Only available when the ChemDarwin tab module is installed.',
        },
        {
            'id': 'agent',
            'title': 'DELFIN Agent',
            'widget': tab_ag,
            'default_order': 55,
            'default_visible': _agent_backend_available,
            'available': True,
            'fixed': False,
            'reason': '' if _agent_backend_available else 'Install the Anthropic CLI or set ANTHROPIC_API_KEY.',
        },
        {
            'id': 'agent_activity',
            'title': 'Agent Activity',
            'widget': tab_ag_act,
            'default_order': 56,
            'default_visible': _agent_backend_available,
            'available': True,
            'fixed': False,
            'reason': '' if _agent_backend_available else 'Install the Anthropic CLI or set ANTHROPIC_API_KEY.',
        },
        {
            'id': 'job_status',
            'title': 'Job Status',
            'widget': tab3,
            'default_order': 60,
            'default_visible': True,
            'available': True,
            'fixed': False,
            'reason': '',
        },
        {
            'id': 'calculations',
            'title': 'Calculations',
            'widget': tab5,
            'default_order': 70,
            'default_visible': True,
            'available': True,
            'fixed': False,
            'reason': '',
        },
        {
            'id': 'archive',
            'title': 'Archive',
            'widget': tab6,
            'default_order': 80,
            'default_visible': True,
            'available': True,
            'fixed': False,
            # Hidden in Settings means not built; showing it again builds it.
            'build': lambda: _clone_of_browser_now(tab_archive_statistics.create_tab, refs6),
            'reason': '',
        },
        {
            'id': 'home',
            'title': 'Home',
            'widget': tab_hm,
            'default_order': 81,
            'default_visible': False,
            'available': True,
            'fixed': False,
            'build': lambda: _clone_of_browser_now(tab_home.create_tab, refs_home),
            'reason': 'Off by default: it browses the whole home directory.',
        },
        {
            'id': 'office',
            'title': 'Office',
            'widget': tab_off,
            'default_order': 82,
            'default_visible': True,
            'available': True,
            'fixed': False,
            # Hidden in Settings means not built; showing it again builds it.
            'build': lambda: _clone_of_browser_now(tab_office.create_tab, refs_off),
            'reason': '',
        },
        {
            'id': 'literature',
            'title': 'Literature',
            'widget': tab_lit,
            'default_order': 85,
            'default_visible': True,
            'available': True,
            'fixed': False,
            'reason': '',
        },
        {
            'id': 'remote_archive',
            'title': 'Remote Archive',
            'widget': tab7,
            'default_order': 90,
            'default_visible': False,
            'available': tab7 is not None,
            'fixed': False,
            'reason': '' if tab7 is not None else 'Enable Remote Archive in Settings first.',
        },
        {
            'id': 'tools',
            'title': 'Tools',
            'widget': tab_tools_widget,
            'default_order': 9_000,
            'default_visible': True,
            'available': tab_tools_widget is not None,
            'fixed': False,
            'reason': '' if tab_tools_widget is not None else 'Tools panel unavailable.',
        },
        {
            'id': 'settings',
            'title': 'Settings',
            'widget': None,
            'default_order': 10_000,
            'default_visible': True,
            'available': True,
            'fixed': True,
            'reason': 'Always visible.',
        },
    ]

    try:
        from . import tab_chemdarwin
        tab_cd, _ = tab_chemdarwin.create_tab(ctx)
    except Exception:
        tab_cd = None

    for spec in tab_specs:
        if spec['id'] == 'chemdarwin':
            spec['widget'] = tab_cd
            spec['available'] = tab_cd is not None
            spec['reason'] = '' if tab_cd is not None else spec['reason']
            if tab_cd is not None:
                spec['default_visible'] = True
            break

    # Additively append dynamically-registered tabs (tab_registry). No-op when
    # nothing is registered; defensive so a registered tab can never break
    # dashboard startup. Built-in tabs above are untouched.
    try:
        from delfin.dashboard.tab_registry import registered_tab_specs
        tab_specs.extend(registered_tab_specs(ctx))
    except Exception:
        pass

    ctx.tab_specs = tab_specs
    tab8, _refs8 = tab_settings.create_tab(ctx, calc_refs=refs5, archive_refs=refs6, office_refs=refs_off)
    for spec in ctx.tab_specs:
        if spec['id'] == 'settings':
            spec['widget'] = tab8
            break

    # Every tab's startup script goes out in ONE ctx.run_js() call: run_js
    # clears its output first, so a second call would wipe the first tab's
    # script before the browser ever executed it.
    #
    # The scripts are collected from ctx (each tab registered its own via
    # ctx.add_init_js) rather than listed here. When they were listed here,
    # the Office tab was added without being added to the list, and its
    # splitter, upload target and viewer resize silently never bound.
    # The bundled 3Dmol goes first: every viewer's loader checks for the
    # library before fetching it, so having it already there is what makes
    # the dashboard work without outbound network.
    # Supersampling factor for every viewer on the page, read once here so a
    # quality change takes effect on the next render rather than per call site.
    try:
        from .molecule_viewer import get_viewer_profile
        _pixel_ratio = float(get_viewer_profile().get('pixel_ratio') or 2.0)
    except Exception:
        _pixel_ratio = 2.0
    _calc_init = '\n'.join([
        vendored_3dmol_js(),
        f'window.__delfinViewerPixelRatio = {_pixel_ratio};',
        RIGHT_MOUSE_TRANSLATE_PATCH_JS,
        *ctx.init_js_parts,
    ])
    if _calc_init.strip():
        ctx.keep_js(_calc_init)

    def _resolve_tab_preferences():
        specs = ctx.tab_specs or []
        spec_map = {spec.get('id'): spec for spec in specs}
        try:
            prefs = (load_settings().get('ui', {}) or {}).get('tabs', {}) or {}
        except Exception:
            prefs = {}
        raw_order = [str(item) for item in (prefs.get('order', []) or []) if str(item).strip()]
        raw_hidden = {str(item) for item in (prefs.get('hidden', []) or []) if str(item).strip()}
        order = []
        seen = set()
        for tab_id in raw_order:
            spec = spec_map.get(tab_id)
            if not spec or spec.get('fixed') or tab_id in seen:
                continue
            seen.add(tab_id)
            order.append(tab_id)
        for spec in sorted(specs, key=lambda item: item.get('default_order', 10_000)):
            tab_id = spec.get('id')
            if not tab_id or spec.get('fixed') or tab_id in seen:
                continue
            seen.add(tab_id)
            order.append(tab_id)
        hidden = set()
        saved_ids = set(raw_order)
        for tab_id in order:
            spec = spec_map.get(tab_id) or {}
            if tab_id in saved_ids:
                if tab_id in raw_hidden:
                    hidden.add(tab_id)
            else:
                if not bool(spec.get('default_visible', True)):
                    hidden.add(tab_id)
        return order, hidden

    def _build_tab_now(spec):
        """Build a tab that was left out because Settings had it hidden.

        Its startup script cannot go out through ctx.run_js: that clears the
        output before it writes, which would wipe the script every tab already
        on the page is running from.
        """
        build = spec.pop('build', None)
        if not callable(build):
            return
        already = len(ctx.init_js_parts)
        try:
            spec['widget'] = build()
        except Exception:
            spec['widget'] = None
            return
        fresh = '\n'.join(ctx.init_js_parts[already:])
        if fresh.strip():
            _append_js(ctx, fresh)

    def _sorted_visible_specs():
        order, hidden = _resolve_tab_preferences()
        order_index = {tab_id: idx for idx, tab_id in enumerate(order)}
        visible_specs = []
        fixed_specs = []
        for spec in ctx.tab_specs:
            if not spec.get('available'):
                continue
            if spec.get('fixed'):
                fixed_specs.append(spec)
                continue
            if spec['id'] in hidden:
                continue
            visible_specs.append(spec)
        visible_specs.sort(
            key=lambda spec: (
                0 if spec['id'] in order_index else 1,
                order_index.get(spec['id'], spec.get('default_order', 10_000)),
                spec.get('default_order', 10_000),
            )
        )
        fixed_specs.sort(key=lambda spec: spec.get('default_order', 10_000))
        return visible_specs + fixed_specs

    refresh_job_list = refs3.get('refresh_job_list')
    _tab_change_state = {'handler': None}

    def _install_tab_observer():
        if not ctx.tabs_widget:
            return
        previous = _tab_change_state.get('handler')
        if previous is not None:
            try:
                ctx.tabs_widget.unobserve(previous, names='selected_index')
            except Exception:
                pass

        job_status_idx = ctx.tab_indices.get('Job Status')

        def _on_tab_change(change):
            if change.get('new') == job_status_idx and refresh_job_list:
                refresh_job_list()

        ctx.tabs_widget.observe(_on_tab_change, names='selected_index')
        _tab_change_state['handler'] = _on_tab_change

    def _rebuild_dashboard_tabs(selected_title=None):
        specs = _sorted_visible_specs()
        for spec in specs:
            if spec.get('widget') is None and callable(spec.get('build')):
                _build_tab_now(spec)
        children = [spec['widget'] for spec in specs if spec.get('widget') is not None]
        titles = [spec['title'] for spec in specs if spec.get('widget') is not None]
        selected_index = titles.index(selected_title) if selected_title in titles else None
        if ctx.tabs_widget is None:
            ctx.tabs_widget = widgets.Tab(children=children)
            for i, title in enumerate(titles):
                ctx.tabs_widget.set_title(i, title)
            if selected_index is not None:
                ctx.tabs_widget.selected_index = selected_index
        else:
            hand_over_tabs(ctx.tabs_widget, children, titles, selected_index)
        ctx.tab_indices = {title: i for i, title in enumerate(titles)}
        _install_tab_observer()

    ctx.rebuild_dashboard_tabs = _rebuild_dashboard_tabs

    # Disable spellcheck in all textareas (browser-level red underlines).
    disable_spellcheck_global(ctx)
    # Text boxes marked 'delfin-debounced' report their value once typing
    # stops, so a handler that renders on every keystroke runs once instead.
    debounce_input(ctx)
    # A rebuilt widget does not move the page the reader is on.
    keep_scroll_position(ctx)
    logo_data_uri = _load_logo_data_uri()
    apply_branding(ctx, title='DELFIN Dashboard', favicon_data_uri=logo_data_uri)

    _rebuild_dashboard_tabs()
    tabs = ctx.tabs_widget

    # -- header bar --------------------------------------------------------
    # One row: the branch switcher, Test, the running previews, Pull, More
    # and the status hint; HOME, the session switch and the busy state on
    # the right. Rare actions (Clean up, roll back) wait behind More.
    # Messages go to one area (header_ui.Toasts) instead of an output
    # under the bar that grew until the page was reloaded.
    from delfin.dashboard import header_ui as _hu
    header_toasts = _hu.Toasts()
    pull_delfin_btn = _hu.icon_button(
        '↓ Pull', 'Update main from origin (fast-forward only)', 'delfin-btn')
    git_hint_btn = _hu.icon_button('…', '', 'delfin-btn', 'delfin-hint')
    git_hint_fold = _hu.icon_button('‹', 'Show only the icon', 'delfin-icon-btn')
    _hint_state = {'compact': bool(_hu.load_prefs().get('hint_compact')),
                   'hint': {'text': '', 'icon': '', 'title': '', 'state': ''}}
    # Branches. Everyone works from main: the shared checkout stays on it,
    # and a branch is tried in a preview -- its own worktree and server --
    # instead of switching this checkout under every person and agent on
    # it (branch_panel / branch_preview). Inside a preview the header shows
    # which code runs instead of the branch tools.
    from delfin.dashboard import branch_preview as _branch_preview
    _in_preview = _branch_preview.in_preview()
    _preview_banner = None
    branch_panel = None
    rollback_delfin_btn = _hu.icon_button(
        '↶ Roll back 1 commit',
        'Move the shared checkout back by one commit (git reset --hard HEAD~1); '
        'refused when it has uncommitted changes', 'delfin-btn', 'is-danger')
    if _in_preview:
        from delfin.dashboard.branch_panel import preview_bar as _preview_bar
        try:
            _preview_banner = _preview_bar(header_toasts)
        except Exception as _banner_err:
            _preview_banner = widgets.HTML(
                '<div class="delfin-preview-banner" style="background:#b71c1c;'
                'color:#fff;padding:6px;">'
                f'PREVIEW (banner failed: {html.escape(str(_banner_err))})</div>')
    else:
        if ctx.repo_dir:
            os.environ.setdefault(_branch_preview.SHARED_CHECKOUT_ENV,
                                  str(ctx.repo_dir))
        from delfin.dashboard.branch_panel import BranchPanel
        branch_panel = BranchPanel(
            ctx.repo_dir, pull=pull_delfin_btn, more_extra=[rollback_delfin_btn],
            hint=widgets.HBox([git_hint_btn, git_hint_fold],
                              layout=widgets.Layout(align_items='center')),
            toasts=header_toasts)

    home_usage_label = _build_home_usage_widget(home)

    def _draw_hint():
        h = _hint_state['hint']
        compact = _hint_state['compact']
        git_hint_btn.description = h['icon'] if compact else h['text']
        git_hint_btn.tooltip = h['title']
        for cls in ('is-behind', 'is-warn', 'is-ok'):
            git_hint_btn.remove_class(cls)
        if h['state']:
            git_hint_btn.add_class(f"is-{h['state']}")
        git_hint_fold.description = '›' if compact else '‹'
        git_hint_fold.tooltip = 'Show the text' if compact else 'Show only the icon'

    def _fold_hint(_b=None):
        _hint_state['compact'] = not _hint_state['compact']
        _hu.save_prefs(hint_compact=_hint_state['compact'])
        _draw_hint()

    git_hint_fold.on_click(_fold_hint)

    def refresh_git_status_label():
        rd = ctx.repo_dir
        if not rd or not Path(rd).exists():
            _hint_state['hint'] = {'text': 'no repository', 'icon': '⚠',
                                   'title': 'No DELFIN checkout found for this dashboard',
                                   'state': 'warn'}
            pull_delfin_btn.disabled = True
            rollback_delfin_btn.disabled = True
            if branch_panel is not None:
                branch_panel.row.layout.display = 'none'
            _draw_hint()
            return
        try:
            _hint_state['hint'] = _git_status_hint(rd)
        except Exception as exc:
            _hint_state['hint'] = {'text': 'git status unavailable', 'icon': '⚠',
                                   'title': str(exc), 'state': 'warn'}
        lock_present = bool(_git_index_lock_message(rd))
        pull_delfin_btn.disabled = lock_present
        rollback_delfin_btn.disabled = lock_present
        _draw_hint()

    refresh_git_status_label()
    if branch_panel is not None:
        branch_panel.on_change = refresh_git_status_label
        branch_panel.refresh()

    def handle_pull_delfin(button=None):
        rd = ctx.repo_dir
        if not rd or not Path(rd).exists():
            header_toasts.show('No DELFIN checkout found for this dashboard.',
                               level='error', details=f'path: {rd}')
            return
        header_toasts.show('Updating main …', level='busy', key='pull')
        ctx.set_busy(True)
        out = []
        level, text = 'error', "Couldn't update main."
        try:
            with _dashboard_git_write_lock(rd, 'git pull'):
                if _ensure_git_index_lock_clear(rd, 'Pull'):
                    out.append('Removed stale Git index lock before pull.')
                # Detached at a branch (one another worktree holds):
                # git pull has no branch to merge into, so move to the
                # branch's current commit, fast-forward only.
                from delfin.agent import held_branch as _hb
                moved = _hb.update_detached(rd)
                if moved is not None:
                    out.append(moved.message)
                    level, text = 'ok', moved.message.splitlines()[0][:160]
                    return
                # Fast-forward main only: the shared checkout is the
                # main checkout, and a merge commit made here would be
                # a commit nobody wrote on a branch everybody runs.
                current = _get_current_git_branch(rd)
                if current != 'main':
                    level = 'error'
                    text = (f'Pull only updates main; the shared checkout is on '
                            f'{current or "a detached commit"}.')
                    out.append('Commit or move its changes, then run "git switch '
                               'main" in a terminal. To try a branch, use Test.')
                    return
                before = _run_git_capture(rd, 'rev-parse', 'HEAD').stdout.strip()
                result = _run_git_write(rd, 'pull', '--ff-only', 'origin', 'main')
                out.append(result.stdout.strip() or '(no output)')
                if result.returncode == 0:
                    n = _run_git_capture(rd, 'rev-list', '--count',
                                         f'{before}..HEAD').stdout.strip()
                    level = 'ok'
                    text = ('main is up to date.' if n in ('', '0')
                            else f'Updated main: {n} new commit(s).')
                else:
                    text = "Couldn't update main (fast-forward only)."
        except RuntimeError as e:
            out.append(str(e))
        except Exception as e:
            out.append(f'Error: {e}')
        finally:
            refresh_git_status_label()
            ctx.set_busy(False)
            header_toasts.show(text, level=level, details='\n'.join(out), key='pull')

    pull_delfin_btn.on_click(handle_pull_delfin)
    git_hint_btn.on_click(handle_pull_delfin)

    def handle_rollback_delfin(button=None):
        rd = ctx.repo_dir
        if not rd or not Path(rd).exists():
            header_toasts.show('No DELFIN checkout found for this dashboard.',
                               level='error', details=f'path: {rd}')
            return
        ctx.set_busy(True)
        out = []
        level, text = 'error', "Couldn't roll back."
        try:
            with _dashboard_git_write_lock(rd, 'git reset --hard HEAD~1'):
                if _ensure_git_index_lock_clear(rd, 'Rollback'):
                    out.append('Removed stale Git index lock before rollback.')

                status_result = _run_git_capture(rd, 'status', '--porcelain')
                if status_result.returncode != 0:
                    out.append(status_result.stdout.strip() or 'git status failed')
                    return
                if status_result.stdout.strip():
                    text = 'Not rolled back: the shared checkout has uncommitted changes.'
                    out.append('Please commit, stash, or discard local changes first.')
                    return

                old_head = _run_git_capture(rd, 'rev-parse', '--short', 'HEAD')
                target_head = _run_git_capture(rd, 'rev-parse', '--short', 'HEAD~1')
                if target_head.returncode != 0:
                    text = 'Not rolled back: there is no earlier commit.'
                    out.append(target_head.stdout.strip())
                    return

                out.append(
                    'Rolling DELFIN back one commit: '
                    f"{old_head.stdout.strip() or 'HEAD'} -> {target_head.stdout.strip()}"
                )
                result = _run_git_write(rd, 'reset', '--hard', 'HEAD~1')
                out.append(result.stdout.strip() or '(no output)')
                if result.returncode == 0:
                    level, text = 'ok', 'Rolled back one commit.'
                    _force_reload_delfin()
        except RuntimeError as e:
            out.append(str(e))
        except Exception as e:
            out.append(f'Error: {e}')
        finally:
            refresh_git_status_label()
            ctx.set_busy(False)
            header_toasts.show(text, level=level, details='\n'.join(out))

    rollback_delfin_btn.on_click(handle_rollback_delfin)

    # Esc, / and Ctrl+K for the branch menu; kept for a resumed page.
    _keys = _hu.keys_js()
    ctx.page_scripts.append(_keys)
    _append_js(ctx, _keys)

    # -- display -----------------------------------------------------------
    title_html = (
        '<div style="display:flex; align-items:center; gap:8px;">'
        + (
            f'<img src="{logo_data_uri}" alt="DELFIN logo" '
            'style="height:34px; width:34px; object-fit:contain;" />'
            if logo_data_uri else ''
        )
        # Wordmark in the logo's own blue (#628FBF, sampled from the dolphin
        # disc) so header and mark read as one unit.
        + '<h1 style="color:#628FBF; margin:0; font-size:28px; '
          'font-weight:700; line-height:1; letter-spacing:1px;">DELFIN</h1>'
        + '</div>'
    )

    # The session control. It sits in the header beside the other
    # controls that change what the dashboard IS rather than what it
    # shows, so a kept session is never invisible.
    _session_strip = _session.build_status_strip()

    # Landing while another session is still running.
    #
    # Offered where people actually land rather than behind a route of
    # its own: they open the address they always open, and it tells them.
    # Excluded by kernel id, so a resume -- which renders INTO the kept
    # kernel -- never offers to go back to the page you are on.
    _returning = _session.build_returning_banner()

    # One beacon, so the server can tell a page that was closed from a
    # connection that dropped. It renders nothing.
    try:
        from delfin.dashboard import window_close as _window_close
        _closing_beacon = _window_close.beacon_widget()
    except Exception:
        _closing_beacon = None

    # The right-hand group: what the dashboard is doing and holding.
    _status_group = widgets.HBox(
        [busy_indicator, ctx.agent_status_html, home_usage_label, _session_strip])
    _status_group.add_class('delfin-group')
    _spacer = widgets.HBox([], layout=widgets.Layout(flex='1 1 auto'))
    _bar = widgets.HBox(
        [widgets.HTML(title_html)]
        + ([_preview_banner] if _preview_banner is not None else [])
        + ([branch_panel.row] if branch_panel is not None else [])
        + [_spacer, _status_group],
        # Wrap rather than shrink: the session strip carries an address,
        # and a row that squeezes its items renders it unreadable before
        # it overflows.
        layout=widgets.Layout(width='100%'))
    _bar.add_class('delfin-topbar')
    _header_root = widgets.VBox(([_returning] if _returning else [])
                                + ([_closing_beacon] if _closing_beacon else [])
                                + [
        busy_css,
        create_page_css(),
        _hu.css_widget(),
        header_toasts.box,
        _bar,
    ]
        + ([branch_panel.box] if branch_panel is not None else []),
        layout=widgets.Layout(width='100%'))
    _header_root.add_class('delfin-hdr')
    _body_root = widgets.VBox([js_output, tabs])

    # The two roots are kept before they are shown, not after.
    #
    # Coming back to a session does not rebuild anything: the widget
    # objects still exist in this kernel with every value and callback
    # they had, so re-displaying THESE is the whole restore. Registering
    # them here is what lets the resume path find them without knowing
    # how the dashboard was assembled -- and it is why the work does not
    # grow when somebody adds a twentieth tab.
    _session.register_root(_header_root, _body_root)

    # The widgets come back by themselves; the scripts around them do
    # not. What the page needs before any viewer in it is worth looking
    # at -- the bundled 3Dmol, every tab's startup script, the editor's
    # bootstraps once sent -- is collected on ctx and handed to a resume
    # to run ahead of the roots.
    _session.register_bootstrap(ctx.resume_bootstrap_js, before=_body_root)

    display(_header_root)
    display(_body_root)

    return ctx


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get_notebook_dir():
    """Best-effort determination of the directory containing the notebook."""
    try:
        import IPython
        return Path(IPython.extract_module_locals()[1]['__vsc_ipynb_file__']).parent
    except Exception:
        return Path.cwd()


def _find_root_dir(start_dir):
    """Search upwards for a directory containing ``software/delfin``."""
    cur = start_dir
    while True:
        if (cur / 'software' / 'delfin').is_dir():
            return cur
        if cur == cur.parent:
            return None
        cur = cur.parent


def _find_delfin_root():
    """Find the DELFIN repo root (directory containing ``delfin/__init__.py``).

    Resolve it from the INSTALLED package location, NOT the cwd: the dashboard's
    repo / git / branch / HOME indicators must always reflect the DELFIN
    checkout, no matter where ``delfin-voila`` was launched. Only the agent tab
    follows the launch dir (via DELFIN_LAUNCH_CWD /
    ``_agent_workspace_from_launch``); everything else stays pinned to the repo.
    Launching from an agent workspace previously made this return ``None``
    (cwd had no ``delfin/``) → the header showed ``git: repo n/a``.
    """
    try:
        import delfin as _delfin
        pkg_root = Path(_delfin.__file__).resolve().parent.parent
        if (pkg_root / 'delfin' / '__init__.py').exists():
            return pkg_root
    except Exception:
        pass
    # Fallback: legacy cwd-based search (e.g. delfin not importable yet).
    for base in [Path.cwd(), *Path.cwd().parents]:
        if (base / 'delfin' / '__init__.py').exists():
            return base
    return None


def _load_logo_data_uri() -> str:
    """Return the packaged DELFIN logo as a data URI when available."""
    try:
        logo_dir = importlib.resources.files('delfin').joinpath('logo')
    except Exception:
        return ''

    try:
        preferred = logo_dir.joinpath('DELFIN_logo.png')
        if preferred.is_file():
            candidates = [preferred]
        else:
            candidates = sorted(
                [
                    entry for entry in logo_dir.iterdir()
                    if entry.is_file() and entry.name.lower().endswith('.png')
                ],
                key=lambda entry: entry.name.lower(),
            )
    except Exception:
        return ''

    if not candidates:
        return ''

    try:
        data = candidates[0].read_bytes()
    except Exception:
        return ''
    if not data:
        return ''
    encoded = base64.b64encode(data).decode('ascii')
    return f'data:image/png;base64,{encoded}'


def _ensure_delfin_importable():
    """Make sure the local DELFIN package is on ``sys.path``."""
    root = _find_delfin_root()
    if root and str(root) not in sys.path:
        sys.path.insert(0, str(root))


def _force_reload_delfin():
    """Reload core DELFIN modules so code changes are picked up."""
    try:
        import delfin.config
        import delfin.common.control_validator
        importlib.reload(delfin.common.control_validator)
        importlib.reload(delfin.config)
    except Exception:
        pass


def _run_git_capture(repo_dir, *args):
    """Run read-only git in the given repo without taking optional index locks."""
    env = os.environ.copy()
    env.setdefault('GIT_OPTIONAL_LOCKS', '0')
    return subprocess.run(
        ['git', '-C', str(repo_dir), *args],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
        env=env,
    )


def _git_dir(repo_dir):
    """This checkout's own git directory.

    ``<repo>/.git`` in a main checkout; in a ``git worktree`` checkout
    ``.git`` is a file and the directory (index, index.lock) is
    ``<common>/.git/worktrees/<name>``, which only git can name."""
    plain = Path(repo_dir) / '.git'
    if plain.is_dir():
        return plain
    result = _run_git_capture(repo_dir, 'rev-parse', '--absolute-git-dir')
    out = (result.stdout or '').strip()
    if result.returncode == 0 and out and Path(out).is_dir():
        return Path(out)
    return plain


@contextlib.contextmanager
def _dashboard_git_write_lock(repo_dir, action_name, timeout=15.0):
    """Serialize repo-changing dashboard Git actions across kernels/processes."""
    lock_path = _git_dir(repo_dir) / 'delfin-dashboard.lock'
    with open(lock_path, 'a+', encoding='utf-8') as handle:
        deadline = time.monotonic() + timeout
        while True:
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                handle.seek(0)
                handle.truncate()
                handle.write(f'{os.getpid()} {action_name}\n')
                handle.flush()
                break
            except BlockingIOError:
                if time.monotonic() >= deadline:
                    raise RuntimeError(
                        'Another DELFIN dashboard session is already changing the Git repo. '
                        'Please wait a moment and retry.'
                    )
                time.sleep(0.2)
        try:
            yield
        finally:
            try:
                handle.seek(0)
                handle.truncate()
                handle.flush()
            except Exception:
                pass
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            except Exception:
                pass


def _git_index_lock_path(repo_dir):
    """Return the Git index lock path for the repository."""
    return _git_dir(repo_dir) / 'index.lock'


def _git_index_lock_message(repo_dir):
    """Describe a present Git index lock, or return an empty string."""
    lock_path = _git_index_lock_path(repo_dir)
    if not lock_path.exists():
        return ''
    try:
        stamp = datetime.fromtimestamp(lock_path.stat().st_mtime).strftime('%Y-%m-%d %H:%M:%S')
        return f'Git index lock present: {lock_path} (mtime {stamp})'
    except Exception:
        return f'Git index lock present: {lock_path}'


def _git_index_lock_in_use(repo_dir):
    """Return whether another process currently holds the Git index lock."""
    lock_path = _git_index_lock_path(repo_dir)
    if not lock_path.exists():
        return False
    try:
        result = subprocess.run(
            ['lsof', str(lock_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            check=False,
        )
    except Exception:
        return False
    return result.returncode == 0


def _ensure_git_index_lock_clear(repo_dir, action_name, wait_seconds=10.0):
    """Wait for an active Git index lock or remove it when it is stale."""
    lock_path = _git_index_lock_path(repo_dir)
    removed_stale_lock = False
    deadline = time.monotonic() + wait_seconds
    while lock_path.exists():
        if _git_index_lock_in_use(repo_dir):
            if time.monotonic() >= deadline:
                raise RuntimeError(
                    f'{action_name} blocked: another git process is still using {lock_path}. '
                    'Please finish that operation or close the other DELFIN session.'
                )
            time.sleep(0.2)
            continue
        try:
            lock_path.unlink()
            removed_stale_lock = True
        except FileNotFoundError:
            break
        except Exception as exc:
            raise RuntimeError(
                f'{action_name} blocked: stale Git lock could not be removed ({exc}).'
            ) from exc
    return removed_stale_lock


def _git_output_mentions_index_lock(output):
    """Return whether Git output points to an index lock conflict."""
    message = (output or '').lower()
    return 'index.lock' in message or 'could not write index' in message


def _run_git_write(repo_dir, *args, retry_on_index_lock=True):
    """Run a repo-changing Git command with stale-index-lock recovery."""
    _ensure_git_index_lock_clear(repo_dir, 'Git write')
    result = subprocess.run(
        ['git', '-C', str(repo_dir), *args],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if result.returncode == 0 or not retry_on_index_lock:
        return result
    if not _git_output_mentions_index_lock(result.stdout):
        return result
    _ensure_git_index_lock_clear(repo_dir, 'Git write retry', wait_seconds=1.0)
    return _run_git_write(repo_dir, *args, retry_on_index_lock=False)


def _format_git_status_html(repo_dir):
    """Build a compact HTML summary of the current repo position."""
    result = _run_git_capture(repo_dir, 'status', '--porcelain=2', '--branch')
    if result.returncode != 0:
        message = result.stdout.strip() or 'git status failed'
        return f'<span style="color:#b71c1c;">git: {html.escape(message)}</span>'

    branch = '?'
    head_short = '?'
    ahead = 0
    behind = 0

    for line in result.stdout.splitlines():
        if line.startswith('# branch.head '):
            branch = line.removeprefix('# branch.head ').strip()
        elif line.startswith('# branch.ab '):
            parts = line.removeprefix('# branch.ab ').strip().split()
            for part in parts:
                if part.startswith('+'):
                    try:
                        ahead = int(part[1:])
                    except ValueError:
                        ahead = 0
                elif part.startswith('-'):
                    try:
                        behind = int(part[1:])
                    except ValueError:
                        behind = 0

    head_result = _run_git_capture(repo_dir, 'rev-parse', '--short', 'HEAD')
    if head_result.returncode == 0 and head_result.stdout.strip():
        head_short = head_result.stdout.strip()

    if branch == '(detached)':
        # A bare "(detached)" / "HEAD" hides that this is usually main's
        # commit run here because another worktree holds the name.
        try:
            from delfin.agent.held_branch import checkout_label
            label = checkout_label(repo_dir)
        except Exception:
            label = ''
        parts = [label or f'{branch} @{head_short}']
    else:
        parts = [f'{branch} @{head_short}']
    if ahead:
        parts.append(f'ahead {ahead}')
    if behind:
        parts.append(f'behind {behind}')
    if not ahead and not behind and branch != '(detached)':
        parts.append('synced')
    if _git_index_lock_message(repo_dir):
        parts.append('index.lock present')

    color = '#b71c1c' if _git_index_lock_message(repo_dir) else '#455a64'
    text = ' | '.join(parts)
    return (
        f'<span style="font-family:monospace; color:{color}; '
        f'padding:2px 6px; border:1px solid #cfd8dc; border-radius:4px;">'
        f'{html.escape(text)}</span>'
    )

def _git_status_hint(repo_dir):
    """The header's status hint for the shared checkout, as data.

    Output: ``{'text', 'icon', 'title', 'state'}`` -- ``text`` a few words
    (``main ✓ up to date``, ``main ↓ 3 new commits``), ``icon`` the same
    in one or two characters for the compact form, ``title`` the full
    status line of :func:`_format_git_status_html` plus what a click does,
    ``state`` ``ok`` / ``behind`` / ``warn``. Ahead and behind are counted
    against the upstream as of the last fetch (a click pulls, which
    fetches).
    """
    result = _run_git_capture(repo_dir, 'status', '--porcelain=2', '--branch')
    if result.returncode != 0:
        msg = result.stdout.strip() or 'git status failed'
        return {'text': 'git status failed', 'icon': '⚠', 'title': msg,
                'state': 'warn'}
    branch, ahead, behind = '?', 0, 0
    for line in result.stdout.splitlines():
        if line.startswith('# branch.head '):
            branch = line.removeprefix('# branch.head ').strip()
        elif line.startswith('# branch.ab '):
            for part in line.removeprefix('# branch.ab ').split():
                try:
                    if part.startswith('+'):
                        ahead = int(part[1:])
                    elif part.startswith('-'):
                        behind = int(part[1:])
                except ValueError:
                    pass
    full = re.sub(r'<[^>]+>', '', _format_git_status_html(repo_dir))
    full = html.unescape(full)
    click = 'Click to update main (fast-forward only).'
    if _git_index_lock_message(repo_dir):
        return {'text': f'{branch} ⚠ git is busy', 'icon': '⚠',
                'title': f'{full}\nAnother git command holds the index lock.',
                'state': 'warn'}
    if branch == '(detached)':
        try:
            from delfin.agent.held_branch import checkout_label
            label = checkout_label(repo_dir).split(' @ ')[0] or 'detached'
        except Exception:
            label = 'detached'
        return {'text': f'{label} (detached)', 'icon': '⚠',
                'title': f'{full}\n{click}', 'state': 'warn'}
    if behind:
        return {'text': f'{branch} ↓ {behind} new commit' + ('s' if behind != 1 else ''),
                'icon': f'↓{behind}', 'title': f'{full}\n{click}', 'state': 'behind'}
    if ahead:
        return {'text': f'{branch} ↑ {ahead} local commit' + ('s' if ahead != 1 else ''),
                'icon': f'↑{ahead}', 'title': f'{full}\nCommits here that origin does '
                'not have.', 'state': 'warn'}
    return {'text': f'{branch} ✓ up to date', 'icon': '✓',
            'title': f'{full} (as of the last fetch)\n{click}', 'state': 'ok'}


def _get_current_git_branch(repo_dir):
    """Return the current local branch name, or ``None`` when detached/unknown."""
    result = _run_git_capture(repo_dir, 'branch', '--show-current')
    if result.returncode != 0:
        return None
    branch = result.stdout.strip()
    return branch or None


def _list_local_git_branches(repo_dir):
    """Return the set of local branch names available in the repo."""
    result = _run_git_capture(repo_dir, 'for-each-ref', '--format=%(refname:short)', 'refs/heads')
    if result.returncode != 0:
        return set()
    return {line.strip() for line in result.stdout.splitlines() if line.strip()}


def _compute_branch_options(repo_dir):
    """Branch dropdown contents as ``[(label, ref)]``: ``main`` first, then
    the user's pinned entries and local branches with commits not in
    ``origin/main`` (branch_list.branch_entries, git only; the panel adds
    the pull-request state in the background)."""
    from delfin.dashboard import branch_list as _bl
    try:
        return [(e.label, e.ref) for e in _bl.branch_entries(repo_dir, gh=None)]
    except Exception:
        return [('main', 'main')]


def _find_orca_candidates(notebook_dir):
    """Find local and system ORCA installations visible to the dashboard."""
    root_dir = _find_root_dir(notebook_dir)
    search_roots = [Path.home() / 'software', Path.home() / 'apps', Path.home() / 'local', Path('/opt')]
    if root_dir:
        search_roots.append(root_dir / 'software')
    try:
        return [Path(path) for path in discover_orca_installations(search_roots=search_roots)]
    except Exception:
        return []


def _find_only_goat_template(notebook_dir):
    """Locate the ``only_GOAT.txt`` template file."""
    root_dir = _find_root_dir(notebook_dir)
    candidates = []
    if root_dir:
        candidates.append(root_dir / 'software' / 'CONTROL_Templates' / 'only_GOAT.txt')
        candidates.append(root_dir / 'CONTROL_Templates' / 'only_GOAT.txt')
    candidates.append(notebook_dir / 'CONTROL_Templates' / 'only_GOAT.txt')
    for p in candidates:
        if p.exists():
            return p
    return candidates[0] if candidates else None


def _build_orca_version_widget(orca_base, orca_candidates):
    """Build an informational ORCA version widget for the dashboard header."""
    import ipywidgets as widgets

    def _label_for(path_value):
        version = Path(path_value).name.replace('orca_', 'ORCA ').replace('_', '.')
        if version == path_value or version == Path(path_value).name:
            version = 'ORCA'
        return f'{version} [{describe_installation_source(path_value)}]'

    if orca_candidates and len(orca_candidates) > 1:
        options = []
        for p in orca_candidates:
            p = Path(p)
            options.append((_label_for(str(p)), str(p)))
        return widgets.Dropdown(
            options=options,
            value=str(orca_base) if orca_base else str(orca_candidates[-1]),
            description='ORCA:',
            layout=widgets.Layout(width='250px'),
            style={'description_width': 'initial'},
            disabled=True,
        )
    if orca_base:
        version = _label_for(str(orca_base))
    else:
        version = 'ORCA (not found)'
    return widgets.Dropdown(
        options=[(version, orca_base or '')],
        value=orca_base or '',
        description='ORCA:',
        layout=widgets.Layout(width='250px'),
        style={'description_width': 'initial'},
        disabled=True,
    )


def _build_home_usage_widget(home_dir, warn_threshold_gb=400):
    """Build an asynchronously refreshed HOME usage widget."""
    import ipywidgets as widgets
    label = widgets.HTML(
        '<span class="delfin-home" title="Space used in your home directory">'
        '🏠 …</span>'
    )

    def _set_label(message, color='#455a64', *, short='', title=''):
        warn = color == '#b71c1c'
        tip = title or message
        label.value = (
            f'<span class="delfin-home{" is-warn" if warn else ""}" '
            f'title="{html.escape(tip)}">🏠 {html.escape(short or message)}</span>'
        )

    def _refresh():
        # Usage comes from quota accounting, never from a directory walk: on a
        # shared cluster file system ``du`` costs one metadata RPC per file and
        # degrades the machine for everyone. See delfin/quota.py.
        try:
            usage = home_usage(home_dir)
        except Exception:
            usage = None
        if usage is None:
            _set_label('HOME: n/a', color='#757575', short='n/a',
                       title='Home usage is not available on this machine')
            return
        used_gb = usage.used_bytes / (1024 ** 3)
        if usage.soft_limit_bytes:
            limit_gb = usage.soft_limit_bytes / (1024 ** 3)
            percent = 100.0 * usage.used_bytes / usage.soft_limit_bytes
            over = percent >= 85.0 or used_gb >= float(warn_threshold_gb)
            _set_label(
                f'HOME: {used_gb:.1f} / {limit_gb:.0f} GB ({percent:.0f}%)',
                color='#b71c1c' if over else '#455a64',
                short=f'{percent:.0f}%',
                title=(f'Home directory: {used_gb:.1f} of {limit_gb:.0f} GB used '
                       f'({percent:.0f}%), from the quota accounting'),
            )
        else:
            color = '#b71c1c' if used_gb >= float(warn_threshold_gb) else '#455a64'
            _set_label(f'HOME: {used_gb:.1f} GB', color=color,
                       short=f'{used_gb:.1f} GB',
                       title=f'Home directory: {used_gb:.1f} GB used')

    threading.Thread(target=_refresh, daemon=True).start()
    return label
