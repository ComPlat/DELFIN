"""DELFIN Dashboard: modular, backend-agnostic dashboard for Jupyter / Voila.

Usage::

    from delfin.dashboard import create_dashboard
    ctx = create_dashboard(backend='auto')
"""

import importlib
import html
import shutil
import subprocess
import sys
import threading
from pathlib import Path

import ipywidgets as widgets
from IPython.display import clear_output, display

from .constants import DEFAULT_CONTROL, ONLY_GOAT_TEMPLATE
from .context import DashboardContext
from .helpers import create_busy_css, disable_spellcheck_global
from .molecule_viewer import RIGHT_MOUSE_TRANSLATE_PATCH_JS
from delfin.runtime_setup import (
    apply_runtime_environment,
    detect_local_runtime_limits,
    get_packaged_submit_templates_dir,
    resolve_backend_choice,
    resolve_orca_base,
    resolve_submit_templates_dir,
)
from delfin.user_settings import load_remote_archive_enabled, load_settings

from . import (
    tab_archive_statistics,
    tab_calculations_browser,
    tab_remote_archive,
    tab_settings,
    tab_job_status,
    tab_orca_builder,
    tab_recalc,
    tab_submit,
)


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
            _find_submit_templates_dir(notebook_dir),
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
        )

        backend_obj = SlurmJobBackend(
            submit_templates_dir=submit_templates_dir,
            orca_base=orca_base,
        )
        only_goat_template_path = _find_only_goat_template(notebook_dir)
    else:
        from .backend_local import LocalJobBackend

        run_script = _find_run_script(notebook_dir)
        detected_local_cores, detected_local_ram_mb = detect_local_runtime_limits()
        orca_base = resolve_orca_base(
            orca_base,
            runtime_settings,
            'local',
            auto_candidates=[],
            local_default='/opt/orca',
        )
        orca_candidates = []
        apply_runtime_environment(
            qm_tools_root=runtime_settings.get('qm_tools_root', ''),
            orca_base=orca_base,
        )
        backend_obj = LocalJobBackend(
            run_script=run_script,
            orca_base=orca_base,
            max_cores=int((runtime_settings.get('local', {}) or {}).get('max_cores', detected_local_cores)),
            max_ram_mb=int((runtime_settings.get('local', {}) or {}).get('max_ram_mb', detected_local_ram_mb)),
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
        run_script=getattr(backend_obj, 'run_script', None),
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
    tab8, _refs8 = tab_settings.create_tab(ctx, calc_refs=refs5, archive_refs=refs6)

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
    )
    if _calc_init.strip():
        ctx.run_js(_calc_init)

    # -- assemble tabs -----------------------------------------------------
    children = [tab1, tab2, tab4]
    titles = ['Submit Job', 'Recalc', 'ORCA Builder']

    if tab_tm is not None:
        children.append(tab_tm)
        titles.append('TURBOMOLE Builder')

    # ChemDarwin tab (optional, local only)
    tab_cd = None
    try:
        from . import tab_chemdarwin
        tab_cd, _ = tab_chemdarwin.create_tab(ctx)
    except Exception:
        tab_cd = None

    if tab_cd is not None:
        children.append(tab_cd)
        titles.append('ChemDarwin')

    children.append(tab3)
    titles.append('Job Status')
    children.append(tab5)
    titles.append('Calculations')
    children.append(tab6)
    titles.append('Archive')
    if tab7 is not None:
        children.append(tab7)
        titles.append('Remote Archive')
    children.append(tab8)
    titles.append('Settings')

    # Disable spellcheck in all textareas (browser-level red underlines).
    disable_spellcheck_global(ctx)

    tabs = widgets.Tab(children=children)
    for i, title in enumerate(titles):
        tabs.set_title(i, title)
    ctx.tabs_widget = tabs
    ctx.tab_indices = {title: i for i, title in enumerate(titles)}

    # Auto-refresh job list when switching to Job Status tab
    job_status_idx = titles.index('Job Status')
    refresh_job_list = refs3.get('refresh_job_list')

    def _on_tab_change(change):
        if change['new'] == job_status_idx and refresh_job_list:
            refresh_job_list()

    tabs.observe(_on_tab_change, names='selected_index')

    # -- header bar --------------------------------------------------------
    backend_label = 'Local' if backend == 'local' else 'BwUniCluster'

    pull_delfin_btn = widgets.Button(
        description='PULL DELFIN', button_style='info',
        layout=widgets.Layout(width='150px'),
    )
    rollback_delfin_btn = widgets.Button(
        description='HEAD -1', button_style='warning',
        layout=widgets.Layout(width='110px'),
    )
    git_status_label = widgets.HTML(value='')
    pull_delfin_output = widgets.Output()

    orca_version_label = _build_orca_version_widget(orca_base, orca_candidates)
    home_usage_label = _build_home_usage_widget(home)

    def refresh_git_status_label():
        rd = ctx.repo_dir
        if not rd or not Path(rd).exists():
            git_status_label.value = '<span style="color:#666;">git: repo n/a</span>'
            return
        git_status_label.value = _format_git_status_html(rd)

    refresh_git_status_label()

    def handle_pull_delfin(button):
        with pull_delfin_output:
            clear_output()
            rd = ctx.repo_dir
            if not rd or not Path(rd).exists():
                print(f'Repo path not found: {rd}')
                return
            print(f'Running git pull in {rd} ...')
            ctx.set_busy(True)
            try:
                result = subprocess.run(
                    ['git', '-C', str(rd), 'pull'],
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, check=False,
                )
                print(result.stdout.strip() or '(no output)')
            except Exception as e:
                print(f'Error: {e}')
            finally:
                refresh_git_status_label()
                ctx.set_busy(False)

    pull_delfin_btn.on_click(handle_pull_delfin)

    def handle_rollback_delfin(button):
        with pull_delfin_output:
            clear_output()
            rd = ctx.repo_dir
            if not rd or not Path(rd).exists():
                print(f'Repo path not found: {rd}')
                return
            ctx.set_busy(True)
            try:
                status_result = subprocess.run(
                    ['git', '-C', str(rd), 'status', '--porcelain'],
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, check=False,
                )
                if status_result.returncode != 0:
                    print(status_result.stdout.strip() or 'git status failed')
                    return
                if status_result.stdout.strip():
                    print('Rollback aborted: working tree is not clean.')
                    print('Please commit, stash, or discard local changes first.')
                    return

                old_head = subprocess.run(
                    ['git', '-C', str(rd), 'rev-parse', '--short', 'HEAD'],
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, check=False,
                )
                target_head = subprocess.run(
                    ['git', '-C', str(rd), 'rev-parse', '--short', 'HEAD~1'],
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, check=False,
                )
                if target_head.returncode != 0:
                    print(target_head.stdout.strip() or 'No previous commit available.')
                    return

                print(
                    'Rolling DELFIN back one commit: '
                    f"{old_head.stdout.strip() or 'HEAD'} -> {target_head.stdout.strip()}"
                )
                result = subprocess.run(
                    ['git', '-C', str(rd), 'reset', '--hard', 'HEAD~1'],
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, check=False,
                )
                print(result.stdout.strip() or '(no output)')
                if result.returncode == 0:
                    _force_reload_delfin()
            except Exception as e:
                print(f'Error: {e}')
            finally:
                refresh_git_status_label()
                ctx.set_busy(False)

    rollback_delfin_btn.on_click(handle_rollback_delfin)

    # -- display -----------------------------------------------------------
    display(widgets.VBox([
        busy_css,
        widgets.HBox([
            widgets.HTML(
                f'<h2 style="color:#1976d2; margin:0;">'
                f'DELFIN ({backend_label})</h2>'
            ),
            widgets.HBox(
                [
                    busy_indicator,
                    home_usage_label,
                    git_status_label,
                    pull_delfin_btn,
                    rollback_delfin_btn,
                    orca_version_label,
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
    """Run git in the given repo and return the completed process."""
    return subprocess.run(
        ['git', '-C', str(repo_dir), *args],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
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

    parts = [f'{branch} @{head_short}']
    if ahead:
        parts.append(f'ahead {ahead}')
    if behind:
        parts.append(f'behind {behind}')
    if not ahead and not behind:
        parts.append('synced')

    color = '#455a64'
    text = ' | '.join(parts)
    return (
        f'<span style="font-family:monospace; color:{color}; '
        f'padding:2px 6px; border:1px solid #cfd8dc; border-radius:4px;">'
        f'{html.escape(text)}</span>'
    )


def _find_submit_templates_dir(notebook_dir):
    """Locate the ``submit_sh`` directory for SLURM templates."""
    root_dir = _find_root_dir(notebook_dir)
    candidates = []
    if root_dir:
        candidates.append(
            root_dir / 'software' / 'delfin' / 'examples'
            / 'example_Job_Submission_Scripts' / 'BwUniCluster' / 'submit_sh'
        )
    candidates.append(notebook_dir / 'submit_sh')
    candidates.append(get_packaged_submit_templates_dir())
    for p in candidates:
        if p.exists():
            return p
    return candidates[0] if candidates else notebook_dir / 'submit_sh'


def _find_orca_candidates(notebook_dir):
    """Find all ``orca_*`` directories under the software folder."""
    root_dir = _find_root_dir(notebook_dir)
    if not root_dir:
        return []
    try:
        orca_root = root_dir / 'software'
        if orca_root.exists():
            return sorted(p for p in orca_root.glob('orca_*') if p.is_dir())
    except Exception:
        pass
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


def _find_run_script(notebook_dir):
    """Locate ``run_local.sh`` for the local backend."""
    candidates = [notebook_dir / 'run_local.sh', Path.cwd() / 'run_local.sh']
    # Look in repo examples
    repo_dir = _find_delfin_root()
    if repo_dir:
        candidates.append(
            repo_dir / 'examples' / 'example_Job_Submission_Scripts'
            / 'LocalServer' / 'run_local.sh'
        )
    # Look in software/delfin if installed under a root_dir
    root_dir = _find_root_dir(notebook_dir)
    if root_dir:
        candidates.append(
            root_dir / 'software' / 'delfin' / 'examples'
            / 'example_Job_Submission_Scripts' / 'LocalServer' / 'run_local.sh'
        )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return notebook_dir / 'run_local.sh'


def _build_orca_version_widget(orca_base, orca_candidates):
    """Build an ORCA version dropdown or label widget."""
    if orca_candidates and len(orca_candidates) > 1:
        options = []
        for p in orca_candidates:
            p = Path(p)
            label = p.name.replace('orca_', 'ORCA ').replace('_', '.')
            options.append((label, str(p)))
        return widgets.Dropdown(
            options=options,
            value=str(orca_base) if orca_base else str(orca_candidates[-1]),
            description='ORCA:',
            layout=widgets.Layout(width='170px'),
            style={'description_width': 'initial'},
        )
    # Single ORCA or none: static label
    if orca_base:
        version = Path(orca_base).name.replace('orca_', 'ORCA ').replace('_', '.')
        if version == orca_base or version == Path(orca_base).name:
            version = 'ORCA'
    else:
        version = 'ORCA (not found)'
    return widgets.Dropdown(
        options=[(version, orca_base or '')],
        value=orca_base or '',
        description='ORCA:',
        layout=widgets.Layout(width='170px'),
        style={'description_width': 'initial'},
        disabled=True,
    )


def _build_home_usage_widget(home_dir, warn_threshold_gb=400):
    """Build an asynchronously refreshed HOME usage widget."""
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
