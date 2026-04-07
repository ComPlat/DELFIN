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
    from IPython.display import clear_output, display

    from .constants import DEFAULT_CONTROL, ONLY_GOAT_TEMPLATE
    from .context import DashboardContext
    from .helpers import apply_branding, create_busy_css, disable_spellcheck_global
    from .molecule_viewer import RIGHT_MOUSE_TRANSLATE_PATCH_JS

    from . import (
        tab_agent,
        tab_archive_statistics,
        tab_calculations_browser,
        tab_remote_archive,
        tab_settings,
        tab_job_status,
        tab_orca_builder,
        tab_recalc,
        tab_submit,
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

    if calc_dir is None:
        calc_dir = configured_paths.get('calculations_dir') or default_calc_dir
    calc_dir = Path(calc_dir).expanduser()
    calc_dir.mkdir(parents=True, exist_ok=True)
    archive_dir = configured_paths.get('archive_dir') or default_archive_dir
    archive_dir = Path(archive_dir).expanduser()
    archive_dir.mkdir(parents=True, exist_ok=True)

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
            tool_binaries=runtime_settings.get('tool_binaries', {}) or {},
        )

        slurm_profile = (runtime_settings.get('slurm', {}) or {}).get('profile', '')
        backend_obj = SlurmJobBackend(
            submit_templates_dir=submit_templates_dir,
            orca_base=orca_base,
            tool_binaries=runtime_settings.get('tool_binaries', {}) or {},
            slurm_profile=slurm_profile,
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
        primary_calc_dir=calc_dir,
        default_calc_dir=default_calc_dir,
        default_archive_dir=default_archive_dir,
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

    # Wire cross-tab ORCA widget references
    ctx.orca_pal_widget = refs4.get('orca_pal')
    ctx.orca_maxcore_widget = refs4.get('orca_maxcore')

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
    tab6, refs6 = tab_archive_statistics.create_tab(ctx)
    tab7, refs7 = (tab_remote_archive.create_tab(ctx) if remote_archive_enabled else (None, {}))
    tab_ag, refs_ag = tab_agent.create_tab(ctx)
    _agent_backend_available = bool(
        shutil.which("claude") or os.environ.get("ANTHROPIC_API_KEY", "")
    )

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
            'reason': '' if _agent_backend_available else 'Install Claude Code CLI or set ANTHROPIC_API_KEY.',
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

    ctx.tab_specs = tab_specs
    tab8, _refs8 = tab_settings.create_tab(ctx, calc_refs=refs5, archive_refs=refs6)
    for spec in ctx.tab_specs:
        if spec['id'] == 'settings':
            spec['widget'] = tab8
            break

    # Run both calc-browser init scripts in ONE ctx.run_js() call.
    # If called separately, the second call's clear_output() would wipe the
    # first tab's splitter-init JS before the browser ever executes it.
    _calc_init = (
        RIGHT_MOUSE_TRANSLATE_PATCH_JS
        + '\n'
        + refs4.get('init_js', '')
        + '\n'
        + refs5.get('init_js', '')
        + '\n'
        + refs6.get('init_js', '')
        + '\n'
        + refs7.get('init_js', '')
        + '\n'
        + refs_ag.get('init_js', '')
    )
    if _calc_init.strip():
        ctx.run_js(_calc_init)

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
        children = [spec['widget'] for spec in specs if spec.get('widget') is not None]
        titles = [spec['title'] for spec in specs if spec.get('widget') is not None]
        if ctx.tabs_widget is None:
            ctx.tabs_widget = widgets.Tab(children=children)
        else:
            ctx.tabs_widget.children = tuple(children)
        for i, title in enumerate(titles):
            ctx.tabs_widget.set_title(i, title)
        ctx.tab_indices = {title: i for i, title in enumerate(titles)}
        _install_tab_observer()
        if selected_title and selected_title in ctx.tab_indices:
            ctx.tabs_widget.selected_index = ctx.tab_indices[selected_title]

    ctx.rebuild_dashboard_tabs = _rebuild_dashboard_tabs

    # Disable spellcheck in all textareas (browser-level red underlines).
    disable_spellcheck_global(ctx)
    logo_data_uri = _load_logo_data_uri()
    apply_branding(ctx, title='DELFIN Dashboard', favicon_data_uri=logo_data_uri)

    _rebuild_dashboard_tabs()
    tabs = ctx.tabs_widget

    # -- header bar --------------------------------------------------------
    backend_label = 'Local' if backend == 'local' else 'BwUniCluster'

    pull_delfin_btn = widgets.Button(
        description='PULL DELFIN', button_style='info',
        layout=widgets.Layout(width='150px'),
    )
    branch_options = [
        ('main', 'main'),
        ('tools-and-workflows', 'tools-and-workflows'),
        ('GUPPY', 'GUPPY'),
    ]
    branch_switch_dropdown = widgets.Dropdown(
        options=branch_options,
        value='main',
        description='Branch:',
        layout=widgets.Layout(width='320px'),
        style={'description_width': 'initial'},
    )
    switch_branch_btn = widgets.Button(
        description='SWITCH', button_style='',
        layout=widgets.Layout(width='110px'),
    )
    rollback_delfin_btn = widgets.Button(
        description='HEAD -1', button_style='warning',
        layout=widgets.Layout(width='110px'),
    )
    git_status_label = widgets.HTML(value='')
    pull_delfin_output = widgets.Output()

    home_usage_label = _build_home_usage_widget(home)

    def refresh_git_status_label():
        rd = ctx.repo_dir
        if not rd or not Path(rd).exists():
            git_status_label.value = '<span style="color:#666;">git: repo n/a</span>'
            pull_delfin_btn.disabled = True
            branch_switch_dropdown.disabled = True
            switch_branch_btn.disabled = True
            rollback_delfin_btn.disabled = True
            return
        git_status_label.value = _format_git_status_html(rd)
        lock_present = bool(_git_index_lock_message(rd))
        pull_delfin_btn.disabled = lock_present
        branch_switch_dropdown.disabled = lock_present
        switch_branch_btn.disabled = lock_present
        rollback_delfin_btn.disabled = lock_present
        current_branch = _get_current_git_branch(rd)
        if current_branch in {value for _label, value in branch_options}:
            branch_switch_dropdown.value = current_branch

    refresh_git_status_label()

    def handle_pull_delfin(button):
        with pull_delfin_output:
            clear_output()
            rd = ctx.repo_dir
            if not rd or not Path(rd).exists():
                print(f'Repo path not found: {rd}')
                return
            ctx.set_busy(True)
            try:
                with _dashboard_git_write_lock(rd, 'git pull'):
                    if _ensure_git_index_lock_clear(rd, 'Pull'):
                        print('Removed stale Git index lock before pull.')
                    print(f'Running git pull in {rd} ...')
                    result = _run_git_write(rd, 'pull')
                    print(result.stdout.strip() or '(no output)')
            except RuntimeError as e:
                print(str(e))
            except Exception as e:
                print(f'Error: {e}')
            finally:
                refresh_git_status_label()
                ctx.set_busy(False)

    pull_delfin_btn.on_click(handle_pull_delfin)

    def handle_switch_branch(button):
        with pull_delfin_output:
            clear_output()
            rd = ctx.repo_dir
            if not rd or not Path(rd).exists():
                print(f'Repo path not found: {rd}')
                return

            target_branch = str(branch_switch_dropdown.value or '').strip()
            current_branch = _get_current_git_branch(rd)
            if not target_branch:
                print('No target branch selected.')
                return
            if target_branch == current_branch:
                print(f'Already on {target_branch}.')
                refresh_git_status_label()
                return

            ctx.set_busy(True)
            try:
                with _dashboard_git_write_lock(rd, f'branch-switch:{current_branch or "detached"}->{target_branch}'):
                    if _ensure_git_index_lock_clear(rd, 'Branch switch'):
                        print('Removed stale Git index lock before branch switch.')

                    stash_name = ''
                    status_result = _run_git_capture(rd, 'status', '--porcelain')
                    if status_result.returncode != 0:
                        print(status_result.stdout.strip() or 'git status failed')
                        return
                    if status_result.stdout.strip():
                        stash_name = (
                            'dashboard-switch-'
                            f'{datetime.utcnow().strftime("%Y%m%d-%H%M%S")}-'
                            f'{current_branch or "detached"}-to-{target_branch}'
                        )
                        print('Working tree is not clean.')

                        backup_dir = Path(rd).parent / 'delfin_branch_switch_backups' / stash_name
                        try:
                            copied = _backup_status_entries(rd, status_result.stdout, backup_dir)
                            if copied:
                                print(f'File backup: {copied} file(s) saved to {backup_dir}')
                            else:
                                print(f'Warning: no files could be backed up (directory created: {backup_dir})')
                        except Exception as backup_err:
                            print(f'Warning: file backup failed ({backup_err}), continuing with git stash only.')

                        print(f'Creating local safety stash: {stash_name}')
                        stash_result = _run_git_write(
                            rd, 'stash', 'push',
                            '--include-untracked',
                            '--message', stash_name,
                        )
                        print(stash_result.stdout.strip() or '(no output)')
                        if stash_result.returncode != 0:
                            print('Branch switch aborted: could not create safety stash.')
                            lock_message = _git_index_lock_message(rd)
                            if lock_message:
                                print(lock_message)
                            return

                    available_branches = _list_local_git_branches(rd)
                    if target_branch not in available_branches:
                        print(f'Branch not available locally: {target_branch}')
                        if stash_name:
                            print(f'Your changes are preserved locally in stash "{stash_name}".')
                        return

                    print(f'Switching DELFIN branch: {current_branch or "?"} -> {target_branch}')
                    result = _run_git_write(rd, 'switch', target_branch)
                    print(result.stdout.strip() or '(no output)')
                    if result.returncode == 0:
                        if stash_name:
                            print(
                                'Local changes were preserved in git stash AND as file backup. '
                                f'Backup: {backup_dir}'
                            )
                        try:
                            stash_list = _run_git_capture(rd, 'stash', 'list')
                            for sline in (stash_list.stdout or '').strip().splitlines():
                                if f'-{target_branch}-to-' in sline:
                                    stash_ref = sline.split(':')[0]
                                    print(f'Found stash from {target_branch}: {stash_ref} – restoring...')
                                    pop_result = _run_git_write(rd, 'stash', 'pop', stash_ref)
                                    print(pop_result.stdout.strip() or '(no output)')
                                    if pop_result.returncode != 0:
                                        print('Auto-restore had conflicts. Check your files.')
                                    break
                        except Exception:
                            pass
                        refresh_git_status_label()
                        print('Reloading dashboard to pick up code from the selected branch...')
                        ctx.run_js('window.location.reload();')
                    elif stash_name:
                        print(
                            f'Branch switch failed, but your changes are preserved in stash "{stash_name}" '
                            f'and as file backup in {backup_dir}.'
                        )
            except RuntimeError as e:
                print(str(e))
            except Exception as e:
                print(f'Error: {e}')
            finally:
                refresh_git_status_label()
                ctx.set_busy(False)

    switch_branch_btn.on_click(handle_switch_branch)

    def handle_rollback_delfin(button):
        with pull_delfin_output:
            clear_output()
            rd = ctx.repo_dir
            if not rd or not Path(rd).exists():
                print(f'Repo path not found: {rd}')
                return
            ctx.set_busy(True)
            try:
                with _dashboard_git_write_lock(rd, 'git reset --hard HEAD~1'):
                    if _ensure_git_index_lock_clear(rd, 'Rollback'):
                        print('Removed stale Git index lock before rollback.')

                    status_result = _run_git_capture(rd, 'status', '--porcelain')
                    if status_result.returncode != 0:
                        print(status_result.stdout.strip() or 'git status failed')
                        return
                    if status_result.stdout.strip():
                        print('Rollback aborted: working tree is not clean.')
                        print('Please commit, stash, or discard local changes first.')
                        return

                    old_head = _run_git_capture(rd, 'rev-parse', '--short', 'HEAD')
                    target_head = _run_git_capture(rd, 'rev-parse', '--short', 'HEAD~1')
                    if target_head.returncode != 0:
                        print(target_head.stdout.strip() or 'No previous commit available.')
                        return

                    print(
                        'Rolling DELFIN back one commit: '
                        f"{old_head.stdout.strip() or 'HEAD'} -> {target_head.stdout.strip()}"
                    )
                    result = _run_git_write(rd, 'reset', '--hard', 'HEAD~1')
                    print(result.stdout.strip() or '(no output)')
                    if result.returncode == 0:
                        _force_reload_delfin()
            except RuntimeError as e:
                print(str(e))
            except Exception as e:
                print(f'Error: {e}')
            finally:
                refresh_git_status_label()
                ctx.set_busy(False)

    rollback_delfin_btn.on_click(handle_rollback_delfin)

    # -- display -----------------------------------------------------------
    title_html = (
        '<div style="display:flex; align-items:center; gap:12px;">'
        + (
            f'<img src="{logo_data_uri}" alt="DELFIN logo" '
            'style="height:46px; width:46px; object-fit:contain;" />'
            if logo_data_uri else ''
        )
        + f'<h2 style="color:#1976d2; margin:0;">DELFIN ({backend_label})</h2>'
        + '</div>'
    )

    display(widgets.VBox([
        busy_css,
        widgets.HBox([
            widgets.HTML(title_html),
            widgets.HBox(
                [
                    busy_indicator,
                    home_usage_label,
                    git_status_label,
                    branch_switch_dropdown,
                    switch_branch_btn,
                    pull_delfin_btn,
                    rollback_delfin_btn,
                ],
                layout=widgets.Layout(
                    margin='0 0 0 12px', align_items='center', gap='8px',
                ),
            ),
        ], layout=widgets.Layout(
            align_items='center', justify_content='space-between', width='100%',
        )),
        pull_delfin_output,
    ], layout=widgets.Layout(width='100%')))
    display(widgets.VBox([js_output, tabs]))

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
    """Find the DELFIN repo root (directory containing ``delfin/__init__.py``)."""
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


@contextlib.contextmanager
def _dashboard_git_write_lock(repo_dir, action_name, timeout=15.0):
    """Serialize repo-changing dashboard Git actions across kernels/processes."""
    lock_path = Path(repo_dir) / '.git' / 'delfin-dashboard.lock'
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
    return Path(repo_dir) / '.git' / 'index.lock'


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


def _backup_status_entries(repo_dir, status_output, backup_dir):
    """Copy dirty files and directories from git status output into a backup directory."""
    backup_dir = Path(backup_dir)
    backup_dir.mkdir(parents=True, exist_ok=True)
    copied = 0

    for line in (status_output or '').strip().splitlines():
        raw_entry = line.lstrip(' MADRCU?!').lstrip().strip()
        entries = []
        for candidate in raw_entry.split(' -> '):
            candidate = candidate.strip().strip('"')
            if candidate and candidate not in entries:
                entries.append(candidate)
        for entry in entries:
            src = Path(repo_dir) / entry
            if not src.exists():
                continue
            dst = backup_dir / entry
            if src.is_dir():
                shutil.copytree(str(src), str(dst), dirs_exist_ok=True)
                copied += 1
            elif src.is_file():
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(str(src), str(dst))
                copied += 1

    return copied


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

    parts = [f'{branch} @{head_short}']
    if ahead:
        parts.append(f'ahead {ahead}')
    if behind:
        parts.append(f'behind {behind}')
    if not ahead and not behind:
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
        '<span style="font-family:monospace; color:#455a64; '
        'padding:2px 6px; border:1px solid #cfd8dc; border-radius:4px;">'
        'HOME: measuring...</span>'
    )

    def _set_label(message, color='#455a64'):
        label.value = (
            f'<span style="font-family:monospace; color:{color}; '
            f'padding:2px 6px; border:1px solid #cfd8dc; border-radius:4px;">'
            f'{html.escape(message)}</span>'
        )

    def _measure_bytes():
        commands = (
            ['du', '-sb', str(home_dir)],
            ['du', '-sk', str(home_dir)],
        )
        for cmd in commands:
            try:
                result = subprocess.run(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False,
                    timeout=1800,
                )
            except Exception:
                continue
            if result.returncode != 0 or not result.stdout.strip():
                continue
            try:
                value = int(result.stdout.split()[0])
            except Exception:
                continue
            if cmd[1] == '-sk':
                value *= 1024
            return value
        raise RuntimeError('du failed')

    def _refresh():
        try:
            used_bytes = _measure_bytes()
            used_gb = used_bytes / (1024 ** 3)
            color = '#b71c1c' if used_gb >= float(warn_threshold_gb) else '#455a64'
            _set_label(f'HOME: {used_gb:.1f} GB', color=color)
        except Exception:
            _set_label('HOME: n/a', color='#757575')

    threading.Thread(target=_refresh, daemon=True).start()
    return label
