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
        DASHBOARD_TABS_CLASS,
        _append_js,
        self_running_html,
        apply_branding,
        create_busy_css,
        create_page_css,
        debounce_input,
        disable_spellcheck_global,
        hand_over_tabs,
        keep_dashboard_tab_in_view,
        keep_scroll_position,
        name_page_controls,
    )
    from .molecule_viewer import RIGHT_MOUSE_TRANSLATE_PATCH_JS, vendored_3dmol_js
    from .tab_defaults import apply_first_tabs, first_start, resolve_tab_preferences

    from . import tab_schedule

    from . import (
        tab_agent,
        tab_agent_activity,
        tab_archive_statistics,
        tab_home,
        tab_office,
        tab_calculations_browser,
        tab_literature,
        tab_settings,
        tab_job_status,
        tab_orca_builder,
        tab_recalc,
        tab_submit,
        tab_tools,
    )
    # tab_remote_archive is NOT imported here: its module load pulls the
    # delfin.smiles_converter -> delfin.manta.ml_tables -> stk chain (seconds
    # of process start). It imports lazily inside its spec 'build', which
    # only runs when the Remote Archive tab is first opened.

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
    # Every tab below the first visible one is deferred: its spec carries a
    # 'build' lambda and widget stays None until the reader opens the tab
    # (tab_schedule.open_spec). A cold start then pays only the first tab's
    # build cost -- the numbers in the lazy-tab task showed the other tabs
    # together cost several times the first tab's build.
    # The refs dicts are filled in by the lambdas on first build, in place,
    # the same way _clone_of_browser_now fills the Settings-held dicts: the
    # dicts are created here so Settings can capture them by value before
    # any tab is built.
    refs1, refs2, refs3, refs4, refs5, refs7, refs_tools = {}, {}, {}, {}, {}, {}, {}

    def _build_now(factory, refs, attribute=None):
        """Run one deferred tab factory, fill its refs dict in place."""
        made = factory(ctx)
        widget = made[0]
        if refs is not None:
            refs.update(made[1] or {})
        if attribute:
            setattr(ctx, attribute, refs)
        return widget
    def _deferred_orca_builder():
        """The ORCA Builder tab: its widgets are the PAL/maxcore source the
        Calculations tab reads, so the build also publishes them on ctx."""
        def build():
            widget = _build_now(
                tab_orca_builder.create_tab, refs4, 'orca_builder_refs')
            ctx.orca_pal_widget = refs4.get('orca_pal')
            ctx.orca_maxcore_widget = refs4.get('orca_maxcore')
            return widget
        return build

    def _deferred(id_, factory, refs, attribute=None):
        return lambda: _build_now(factory, refs, attribute)

    # The first visible tab (Submit Job, default_order 10) is built eagerly:
    # it is what the reader sees first and its refs are the startup
    # contract of the dashboard.
    tab1, refs1 = tab_submit.create_tab(ctx)
    ctx.submit_refs = refs1

    # Wire cross-tab widget references for agent dashboard control. The
    # ORCA builder, job status and recalc refs are filled on the tab's
    # first open through the same dicts the agent reads.
    ctx.submit_refs = refs1

    # TURBOMOLE tab (only for SLURM backends)
    if ctx.backend.supports_turbomole:
        from . import tab_turbomole_builder
        _tm_available = True
    else:
        _tm_available = False

    try:
        remote_archive_enabled = load_remote_archive_enabled()
    except Exception:
        remote_archive_enabled = False

    refs6, refs_off, refs_home = {}, {}, {}

    def _clone_of_browser_now(build, refs):
        """Build one of the browser clones and fill in the refs Settings holds."""
        widget, made = build(ctx)
        refs.update(made or {})
        return widget

    _agent_backend_available = bool(shutil.which("claude"))
    try:
        from delfin.agent.credentials import load_credential as _lc_anth
        _agent_backend_available = bool(
            shutil.which("claude") or _lc_anth("ANTHROPIC_API_KEY")
        )
    except Exception:
        _agent_backend_available = bool(
            shutil.which("claude") or os.environ.get("ANTHROPIC_API_KEY", "")
        )

    def _clone_of_browser(tab_id, build, refs):
        """The lazy build lambda for a browser clone.

        The clone is built through the spec's ``build`` on first open, like
        every tab that is not the first visible one -- a hidden one is
        therefore never built, and switching it on in Settings builds it
        then. (A clone hidden at start used to get no build at all and
        stayed off the bar until a reload.) ``refs`` is the dict the
        Settings tab is handed by value further down; it is filled in place
        when the build runs, so applying a new Archive or Office folder calls
        back through the same dict.
        """
        return lambda: _clone_of_browser_now(build, refs)

    tab6 = _clone_of_browser('archive', tab_archive_statistics.create_tab, refs6)
    tab_off = _clone_of_browser('office', tab_office.create_tab, refs_off)
    tab_hm = _clone_of_browser('home', tab_home.create_tab, refs_home)
    archive_build, office_build, home_build = tab6, tab_off, tab_hm

    # Archive, Office and Home build lazily through their specs' 'build'
    # lambdas on first open; the refs dicts are filled in place so Settings
    # keeps working when it applies a new folder.
    ctx.remote_archive_refs = refs7

    def _deferred_remote_archive():
        """Remote Archive: its module imports on first open (its load pulls
        the stk chain, seconds of process start), not when the spec is
        written."""
        def build():
            from . import tab_remote_archive
            return _build_now(
                tab_remote_archive.create_tab, refs7, 'remote_archive_refs')
        return build

    def _deferred_agent_tab():
        """The DELFIN Agent tab: the session list, or the single agent tab
        when the list itself fails to build."""
        def build():
            try:
                from . import agent_sessions
                return agent_sessions.create_tab(ctx)[0]
            except Exception:
                return tab_agent.create_tab(ctx)[0]
        return build

    def _deferred_tools_tab():
        """The Tools tab; a failing build yields None, never a startup break."""
        def build():
            try:
                widget, made = tab_tools.create_tab(ctx)
            except Exception:
                return None
            refs_tools.update(made or {})
            return widget
        return build

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
            'widget': None,
            'default_order': 20,
            'default_visible': True,
            'available': True,
            'fixed': False,
            'build': _deferred('recalc', tab_recalc.create_tab, refs2, 'recalc_refs'),
            'reason': '',
        },
        {
            'id': 'orca_builder',
            'title': 'ORCA Builder',
            'widget': None,
            'default_order': 30,
            'default_visible': True,
            'available': True,
            'fixed': False,
            'build': _deferred_orca_builder(),
            'reason': '',
        },
        {
            'id': 'turbomole_builder',
            'title': 'TURBOMOLE Builder',
            'widget': None,
            'default_order': 40,
            'default_visible': _tm_available,
            'available': _tm_available,
            'fixed': False,
            'build': _deferred('turbomole', tab_turbomole_builder.create_tab, None) if _tm_available else None,
            'reason': '' if _tm_available else 'Only available for supported SLURM backends.',
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
            'widget': None,
            'default_order': 55,
            'default_visible': _agent_backend_available,
            'available': True,
            'fixed': False,
            'build': _deferred_agent_tab(),
            'reason': '' if _agent_backend_available else 'Install the Anthropic CLI or set ANTHROPIC_API_KEY.',
        },
        {
            'id': 'agent_activity',
            'title': 'Agent Activity',
            'widget': None,
            'default_order': 56,
            'default_visible': _agent_backend_available,
            'available': True,
            'fixed': False,
            'build': _deferred('agent_activity', tab_agent_activity.create_tab, None),
            'reason': '' if _agent_backend_available else 'Install the Anthropic CLI or set ANTHROPIC_API_KEY.',
        },
        {
            'id': 'job_status',
            'title': 'Job Status',
            'widget': None,
            'default_order': 60,
            'default_visible': True,
            'available': True,
            'fixed': False,
            'build': _deferred('job_status', tab_job_status.create_tab, refs3, 'job_status_refs'),
            'reason': '',
        },
        {
            'id': 'calculations',
            'title': 'Calculations',
            'widget': None,
            'default_order': 70,
            'default_visible': True,
            'available': True,
            'fixed': False,
            'build': _deferred('calculations', tab_calculations_browser.create_tab, refs5, 'calc_browser_refs'),
            'reason': '',
        },
        {
            'id': 'archive',
            'title': 'Archive',
            'widget': None,
            'default_order': 80,
            'default_visible': True,
            'available': True,
            'fixed': False,
            'build': archive_build,
            'reason': '',
        },
        {
            'id': 'home',
            'title': 'Home',
            'widget': None,
            'default_order': 81,
            # Off in the view of before 2026-11-13, which an existing user
            # keeps; one of the five a new installation starts with
            # (tab_defaults.FIRST_TABS, first_start).
            'default_visible': False,
            'available': True,
            'fixed': False,
            'build': home_build,
            'reason': '',
        },
        {
            'id': 'office',
            'title': 'Office',
            'widget': None,
            'default_order': 82,
            'default_visible': True,
            'available': True,
            'fixed': False,
            'build': office_build,
            'reason': '',
        },
        {
            'id': 'literature',
            'title': 'Literature',
            'widget': None,
            'default_order': 85,
            'default_visible': True,
            'available': True,
            'fixed': False,
            'build': _deferred('literature', tab_literature.create_tab, None),
            'reason': '',
        },
        {
            'id': 'remote_archive',
            'title': 'Remote Archive',
            'widget': None,
            'default_order': 90,
            'default_visible': False,
            'available': remote_archive_enabled,
            'fixed': False,
            'build': _deferred_remote_archive() if remote_archive_enabled else None,
            'reason': '' if remote_archive_enabled else 'Enable Remote Archive in Settings first.',
        },
        {
            'id': 'tools',
            'title': 'Tools',
            'widget': None,
            'default_order': 9_000,
            'default_visible': True,
            'available': True,
            'fixed': False,
            'build': _deferred_tools_tab(),
            'reason': 'Tools panel unavailable if its build fails.',
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

    def _chemdarwin_available():
        """Import the ChemDarwin tab module once, to know it is installed.

        Its module load pulls the stk chain (seconds of process start), so
        the probe runs at first open through the spec's build, not at
        dashboard startup. Result is cached for later spec reads.
        """
        try:
            from . import tab_chemdarwin
        except Exception:
            return None
        return tab_chemdarwin

    _chemdarwin_state = {'module': None, 'probed': False}

    def _chemdarwin_build():
        """Build ChemDarwin on first open, probing its module then.

        The module load pulls the stk chain (seconds of process start), so
        the availability probe must not run at dashboard startup. A module
        that is not installed leaves the tab unavailable.
        """
        if not _chemdarwin_state['probed']:
            _chemdarwin_state['module'] = _chemdarwin_available()
            _chemdarwin_state['probed'] = True
        module = _chemdarwin_state['module']
        if module is None:
            return None
        return module.create_tab(ctx)[0]

    # Whether the tab can exist is decided without importing it: the module
    # and its one hard dependency beyond the dashboard's own are findable.
    # The import itself waits for the first open (_chemdarwin_build).
    import importlib.util as _ilu
    try:
        _cd_findable = all(_ilu.find_spec(name) is not None for name in
                           ('delfin.dashboard.tab_chemdarwin', 'rdkit'))
    except Exception:
        _cd_findable = False
    for spec in tab_specs:
        if spec['id'] == 'chemdarwin':
            spec['build'] = _chemdarwin_build
            spec['available'] = _cd_findable
            if _cd_findable:
                spec['default_visible'] = True
                spec['reason'] = ''
            break

    # Additively append dynamically-registered tabs (tab_registry). No-op when
    # nothing is registered; defensive so a registered tab can never break
    # dashboard startup. Built-in tabs above are untouched.
    try:
        from delfin.dashboard.tab_registry import registered_tab_specs
        tab_specs.extend(registered_tab_specs(ctx))
    except Exception:
        pass

    # A new installation starts with five tabs and Settings; anyone who used
    # DELFIN before keeps the tabs they had. Decided once, at the first
    # start, and saved as the person's arrangement before the Settings tab
    # reads it (tab_defaults.first_start); every other tab is switched on
    # in Settings > Dashboard Tabs.
    apply_first_tabs(tab_specs)
    try:
        from delfin.agent import session_store as _session_store
        _sessions_dir = getattr(_session_store, '_SESSIONS_DIR', None)
    except Exception:
        _sessions_dir = None
    try:
        first_start(tab_specs, calc_dir=ctx.calc_dir, archive_dir=ctx.archive_dir,
                    sessions_dir=_sessions_dir,
                    jobs_file=Path.home() / '.delfin_jobs.json')
    except Exception:
        pass

    ctx.tab_specs = tab_specs
    tab8, _refs8 = tab_settings.create_tab(ctx, calc_refs=refs5, archive_refs=refs6, office_refs=refs_off)
    for spec in ctx.tab_specs:
        if spec['id'] == 'settings':
            spec['widget'] = tab8
            break

    # Unnamed selects and slider handles take their widget's words as an
    # accessible name; every tab takes its full name as a tooltip.
    name_page_controls(ctx)
    # The tab strip, where it scrolls: faded edges, wheel, current tab in view.
    keep_dashboard_tab_in_view(ctx)

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
    # Without the bundle no viewer can draw (there is no network fallback),
    # so a broken installation says so in the browser console.
    from delfin.dashboard.login_session import keepalive_js as _keepalive_js
    _startup_parts = [
        vendored_3dmol_js() or (
            'console.error("DELFIN: the bundled 3Dmol '
            '(delfin/dashboard/static/3Dmol-min.js) is missing from this '
            'installation; molecule viewers cannot draw.");'),
        f'window.__delfinViewerPixelRatio = {_pixel_ratio};',
        # User input keeps the login alive (login_session.keepalive_js).
        _keepalive_js(),
        RIGHT_MOUSE_TRANSLATE_PATCH_JS,
        *ctx.init_js_parts,
    ]
    _calc_init = '\n'.join(_startup_parts)
    # The page's startup script runs from the page itself (an element in the
    # body, helpers.self_running_html), not from ctx.js_output: run_js
    # clears that output before it writes, and with every tab built at
    # start a tab's own run_js (Remote Archive's preview reset) cleared the
    # startup script before the browser ran it -- no 3Dmol (no molecule
    # viewer drew), no login keep-alive (measured in a real browser,
    # 11.11.). Kept in page_scripts for a resumed window as before.
    _startup_runner = None
    if _calc_init.strip():
        ctx.page_scripts.append(_calc_init)
        _startup_runner = widgets.HTML(
            self_running_html(*_startup_parts),
            layout=widgets.Layout(width='0px', height='0px', margin='0',
                                  overflow='hidden'))

    def _resolve_tab_preferences():
        try:
            prefs = (load_settings().get('ui', {}) or {}).get('tabs', {}) or {}
        except Exception:
            prefs = {}
        return resolve_tab_preferences(ctx.tab_specs or [], prefs)

    # A tab that drives other tabs' widgets (through ctx.submit_refs,
    # ctx.orca_builder_refs, ...) needs them built first: the agent fills
    # the ORCA Builder, submits and reads the job list whichever tab the
    # reader opened before it.
    _TAB_NEEDS = {
        'agent': ('submit', 'recalc', 'orca_builder', 'job_status', 'calculations'),
    }

    def _build_tab_now(spec):
        """Build a deferred tab: on its first open, or one Settings had hidden.

        Its startup script cannot go out through ctx.run_js: that clears the
        output before it writes, which would wipe the script every tab already
        on the page is running from.
        """
        build = spec.pop('build', None)
        if not callable(build):
            return
        for need in _TAB_NEEDS.get(spec.get('id'), ()):
            for other in ctx.tab_specs or []:
                if other.get('id') == need and other.get('widget') is None:
                    _build_tab_now(other)
        already = len(ctx.init_js_parts)
        try:
            spec['widget'] = build()
        except Exception:
            spec['widget'] = None
            return
        # The tab carries its own init scripts and runs them when it is
        # shown (helpers.self_running_html). Sent as a Javascript output
        # into ctx.js_output after the page was out, they never ran: every
        # tab built on its first open lost its splitter, explorer clicks,
        # upload and keys (Calculations, Remote Archive, ... since 09.11.).
        parts = [p for p in ctx.init_js_parts[already:] if p and str(p).strip()]
        if parts and spec.get('widget') is not None:
            runner = widgets.HTML(self_running_html(*parts), layout=widgets.Layout(
                width='0px', height='0px', margin='0', overflow='hidden'))
            built = spec['widget']
            if isinstance(built, widgets.Box):
                built.children = tuple(built.children) + (runner,)
            else:
                spec['widget'] = widgets.VBox([built, runner])

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

    refresh_job_list_holder = {'getter': None}
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

        def _on_tab_change(change):
            # A tab that is on the bar but not yet built builds here, on
            # its first open, and keeps its state from then on.
            index = change.get('new')
            if index is None or not ctx.tab_indices:
                return
            title = ctx.tab_indices_title_of(index)
            if title is None:
                return
            for spec in ctx.tab_specs or []:
                if spec.get('title') == title and spec.get('widget') is None:
                    _build_tab_now(spec)
                    if spec.get('widget') is not None:
                        ctx.rebuild_dashboard_tabs(spec['title'])
                    break
            if title == 'Job Status':
                getter = refresh_job_list_holder['getter']
                if callable(getter):
                    getter()

        ctx.tabs_widget.observe(_on_tab_change, names='selected_index')
        _tab_change_state['handler'] = _on_tab_change

    def _bar_widget(spec):
        """What the bar shows for *spec*: its built widget, or an empty
        placeholder while it waits for its first open; None for a tab that
        has neither (a build that failed)."""
        if spec.get('widget') is not None:
            return spec['widget']
        if not callable(spec.get('build')):
            return None
        if spec.get('placeholder') is None:
            spec['placeholder'] = widgets.VBox()
        return spec['placeholder']

    def _rebuild_dashboard_tabs(selected_title=None):
        specs = _sorted_visible_specs()
        # Only the tab the reader sees is built: at the first build the
        # first visible one, later the one being selected. Every other tab
        # stays on the bar as a placeholder until its first open (see
        # _on_tab_change) -- a cold start pays one tab's build cost instead
        # of the bar's, and no tab goes missing from the bar.
        if ctx.tabs_widget is None:
            to_build = tab_schedule.startup_specs_to_build(specs)
        else:
            to_build = [spec for spec in specs if spec['title'] == selected_title]
        for spec in to_build:
            if spec.get('widget') is None and callable(spec.get('build')):
                _build_tab_now(spec)
        shown = [(spec['title'], _bar_widget(spec)) for spec in specs]
        shown = [(title, widget) for title, widget in shown if widget is not None]
        children = [widget for _, widget in shown]
        titles = [title for title, _ in shown]
        selected_index = titles.index(selected_title) if selected_title in titles else None
        if ctx.tabs_widget is None:
            ctx.tabs_widget = widgets.Tab(children=children)
            ctx.tabs_widget.add_class(DASHBOARD_TABS_CLASS)
            for i, title in enumerate(titles):
                ctx.tabs_widget.set_title(i, title)
            if selected_index is not None:
                ctx.tabs_widget.selected_index = selected_index
        else:
            hand_over_tabs(ctx.tabs_widget, children, titles, selected_index)
        ctx.tab_indices = {title: i for i, title in enumerate(titles)}
        ctx.tab_indices_title_of = (
            lambda index: titles[index]
            if index is not None and 0 <= index < len(titles) else None)
        refresh_job_list_holder['getter'] = (
            ctx.job_status_refs or {}).get('refresh_job_list')
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
        'Move the shared checkout back by one commit (git reset --hard to the '
        'first parent) after a confirmation that names both commits; refused '
        'when it has uncommitted changes', 'delfin-btn', 'is-danger')
    from delfin.dashboard import rollback as _rollback
    rollback_confirm = _rollback.Confirm()
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
            more_below=[rollback_confirm.box],
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
        installed = _installed_copy_hint(rd)
        if installed is not None:
            _hint_state['hint'] = installed
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
        try:
            from delfin.dashboard.branch_list import pull_unavailable
            no_pull = pull_unavailable(rd)
        except Exception:
            no_pull = ''
        pull_delfin_btn.disabled = lock_present or bool(no_pull)
        pull_delfin_btn.tooltip = (no_pull or
                                   'Update main from origin (fast-forward only)')
        rollback_delfin_btn.disabled = lock_present
        _draw_hint()

    if branch_panel is not None:
        branch_panel.on_change = refresh_git_status_label
        branch_panel.on_more_closed(rollback_confirm.cancel)

    def _first_header_read():
        """The git hint and the branch list, read once beside the page.

        Both are git reads (a branch list of a few hundred branches takes
        seconds on a cluster file system); the page does not wait for
        them. Until they land the hint shows its placeholder and the
        branch menu reads the list again when it is opened. Pull and
        roll back check the checkout themselves on the click.
        """
        try:
            refresh_git_status_label()
        except Exception:
            pass
        if branch_panel is not None:
            try:
                branch_panel.refresh()
            except Exception:
                pass

    threading.Thread(target=_first_header_read, daemon=True,
                     name='delfin-header-first-read').start()
    if branch_panel is not None:
        branch_panel.start_polling(10.0)

    def handle_pull_delfin(button=None):
        rd = ctx.repo_dir
        if not rd or not Path(rd).exists():
            header_toasts.show('No DELFIN checkout found for this dashboard.',
                               level='error', details=f'path: {rd}')
            return
        try:
            from delfin.dashboard.branch_list import pull_unavailable
            no_pull = pull_unavailable(rd)
        except Exception:
            no_pull = ''
        if no_pull:
            header_toasts.show("Pull is not available.", level='error',
                               details=no_pull, key='pull')
            return
        header_toasts.show('Updating main …', level='busy', key='pull')
        ctx.set_busy(True)
        out = []
        level, text = 'error', "Couldn't update main."
        # The pull itself lives in delfin/agent/branch_ops.pull_main
        # (UI-free, shared with the Pi's /pull); it takes the same write
        # lock and answers {'ok', 'text'}. The dashboard keeps its own
        # busy state and toast.
        try:
            from delfin.agent import branch_ops
            result = branch_ops.pull_main(rd)
            level = 'ok' if result.get('ok') else 'error'
            text = result.get('text') or ("Couldn't update main." if level == 'error'
                                          else 'main is up to date.')
            out.append(text)
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
        """Ask first (rollback.Confirm); the reset runs on the confirm."""
        rd = ctx.repo_dir
        if not rd or not Path(rd).exists():
            header_toasts.show('No DELFIN checkout found for this dashboard.',
                               level='error', details=f'path: {rd}')
            return
        try:
            plan, why = _rollback.plan(rd, run=_run_git_capture)
        except Exception as e:
            plan, why = None, f'Error: {e}'
        if plan is None:
            rollback_confirm.cancel()
            header_toasts.show(why, level='error')
            return
        rollback_confirm.ask(_rollback.describe(plan),
                             f'Roll back to {plan.target_short}',
                             lambda: _do_rollback_delfin(plan))

    def _do_rollback_delfin(plan):
        rd = ctx.repo_dir
        ctx.set_busy(True)
        out = []
        level, text = 'error', "Couldn't roll back."
        try:
            with _dashboard_git_write_lock(rd, f'git reset --hard {plan.target_short}'):
                if _ensure_git_index_lock_clear(rd, 'Rollback'):
                    out.append('Removed stale Git index lock before rollback.')
                ok, text, details = _rollback.execute(
                    rd, plan, read=_run_git_capture,
                    write=_run_git_write)
                out.extend(details)
                if ok:
                    level = 'ok'
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
    class _PageRoot(widgets.VBox):
        """One of the page's two root boxes, with a short plain-text form.

        display() asks every object for its text/plain form beside the
        widget view; a widget's default is the repr of its whole tree,
        every tab's widgets and the startup scripts included -- built on
        every display and then cut to 110 characters. The page shows the
        widget view, never that text.
        """

        def __repr__(self):
            return f'{type(self).__name__}(DELFIN dashboard)'

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
    # In a branch preview the strip makes no sense: a preview is its own
    # worktree and server, judged and then thrown away -- keeping a
    # session there would save the preview, not the real dashboard.
    _session_strip = None
    if not _in_preview:
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

    # Logging out ends the login on the server (login_session); the
    # control is the last thing in the bar.
    from delfin.dashboard import login_session as _login_session
    _logout = widgets.HTML(_login_session.logout_control_html())

    # The right-hand group: what the dashboard is holding. What it is doing
    # (the busy spinner, the agent's working line) sits just left of it, in
    # the spacer, where its changing text cannot wrap the bar (top_bar).
    _status_group = widgets.HBox(
        [home_usage_label]
        + ([_session_strip] if _session_strip is not None else [])
        + [_logout])
    # Wrap rather than shrink: the session strip carries an address, and a
    # row that squeezes its items renders it unreadable before it overflows.
    _bar = _hu.top_bar(
        [widgets.HTML(title_html)]
        + ([_preview_banner] if _preview_banner is not None else [])
        + ([branch_panel.row] if branch_panel is not None else []),
        [busy_indicator, ctx.agent_status_html],
        _status_group)
    _header_root = _PageRoot(([_returning] if _returning else [])
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
    _body_root = _PageRoot(
        [js_output] + ([_startup_runner] if _startup_runner is not None else [])
        + [tabs])

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


# The git write/capture helpers and their lock live in
# delfin/dashboard/git_write.py (UI-free): the agent's branch operations
# (delfin/agent/branch_ops.py) take the same lock and clear the same
# stale index.lock, so a dashboard Pull and a service pull cannot
# interleave. The names stay: the dashboard's callers and the tests
# pinning them use the underscore spellings.
from delfin.dashboard.git_write import (
    run_git_capture as _run_git_capture,
    git_dir as _git_dir,
    dashboard_git_write_lock as _dashboard_git_write_lock,
    git_index_lock_path as _git_index_lock_path,
    git_index_lock_message as _git_index_lock_message,
    git_index_lock_in_use as _git_index_lock_in_use,
    ensure_git_index_lock_clear as _ensure_git_index_lock_clear,
    git_output_mentions_index_lock as _git_output_mentions_index_lock,
    run_git_write as _run_git_write,
)



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

def _installed_copy_hint(repo_dir):
    """The header's hint when DELFIN runs from an installed copy; else None.

    The header follows the installed package (_find_delfin_root), and a
    package is not always a checkout: a pip install lives in site-packages,
    a deploy (scripts/pi) is an unpacked archive. Asked for its git status,
    such a directory answered "git status failed" -- a warning about
    nothing, on every page load -- with Pull and Rollback aimed at a
    repository that is not there.
    """
    if (Path(repo_dir) / '.git').exists():
        return None
    try:
        from delfin import __version__
        name = f'DELFIN {__version__}'
    except Exception:
        name = 'DELFIN'
    return {'text': 'installed', 'icon': '◆',
            'title': (f'{name} runs from an installed copy ({repo_dir}), not '
                      'a git checkout: there is nothing to pull or roll back '
                      'here. Update it the way it was installed.'),
            'state': 'ok'}


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
