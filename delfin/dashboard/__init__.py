"""DELFIN Dashboard: modular, backend-agnostic dashboard for Jupyter / Voila.

Usage::

    from delfin.dashboard import create_dashboard
    ctx = create_dashboard(backend='auto')
"""

import importlib
import shutil
import subprocess
import sys
from pathlib import Path

import ipywidgets as widgets
from IPython.display import clear_output, display

from .constants import DEFAULT_CONTROL, ONLY_GOAT_TEMPLATE
from .context import DashboardContext
from .helpers import create_busy_css, disable_spellcheck_global

from . import (
    tab_calculations_browser,
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

    if calc_dir is None:
        # SLURM: look for a root_dir/calc pattern; local: ~/calc
        root_dir = _find_root_dir(notebook_dir)
        if root_dir and (root_dir / 'calc').exists():
            calc_dir = root_dir / 'calc'
        else:
            calc_dir = home / 'calc'
    calc_dir = Path(calc_dir)
    calc_dir.mkdir(parents=True, exist_ok=True)

    repo_dir = _find_delfin_root()

    # -- auto-detect backend -----------------------------------------------
    if backend == 'auto':
        backend = 'slurm' if shutil.which('sbatch') else 'local'

    # -- build backend object ----------------------------------------------
    if backend == 'slurm':
        from .backend_slurm import SlurmJobBackend

        submit_templates_dir = _find_submit_templates_dir(notebook_dir)
        orca_candidates = _find_orca_candidates(notebook_dir)
        if orca_base is None and orca_candidates:
            # Prefer 6_1_1 if available, else latest
            for p in orca_candidates:
                if '6_1_1' in p.name:
                    orca_base = str(p)
                    break
            if orca_base is None:
                orca_base = str(orca_candidates[-1])

        backend_obj = SlurmJobBackend(
            submit_templates_dir=submit_templates_dir,
            orca_base=orca_base,
        )
        only_goat_template_path = _find_only_goat_template(notebook_dir)
    else:
        from .backend_local import LocalJobBackend

        run_script = _find_run_script(notebook_dir)
        if orca_base is None:
            orca_base = '/opt/orca'
        orca_candidates = []
        backend_obj = LocalJobBackend(
            run_script=run_script,
            orca_base=orca_base,
        )
        only_goat_template_path = None

    # -- shared widgets ----------------------------------------------------
    js_output = widgets.Output()
    busy_indicator = widgets.HTML(value='')
    busy_css = create_busy_css()

    # -- build context -----------------------------------------------------
    ctx = DashboardContext(
        calc_dir=calc_dir,
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

    tab5, refs5 = tab_calculations_browser.create_tab(ctx)

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

    # Disable spellcheck in all textareas (browser-level red underlines).
    disable_spellcheck_global(ctx)

    tabs = widgets.Tab(children=children)
    for i, title in enumerate(titles):
        tabs.set_title(i, title)

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
    pull_delfin_output = widgets.Output()

    orca_version_label = _build_orca_version_widget(orca_base, orca_candidates)

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
                ctx.set_busy(False)

    pull_delfin_btn.on_click(handle_pull_delfin)

    # -- display -----------------------------------------------------------
    display(widgets.VBox([
        busy_css,
        widgets.HBox([
            widgets.HTML(
                f'<h2 style="color:#1976d2; margin:0;">'
                f'DELFIN Dashboard ({backend_label})</h2>'
            ),
            widgets.HBox(
                [busy_indicator, pull_delfin_btn, orca_version_label],
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
