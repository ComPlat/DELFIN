"""Dashboard Settings tab backed by a user-local settings file."""

import html
import logging
import os
import shlex
import shutil
from pathlib import Path

import ipywidgets as widgets

from delfin.runtime_setup import (
    apply_runtime_environment,
    collect_bwunicluster_verification,
    collect_runtime_diagnostics,
    describe_installation_source,
    detect_local_runtime_limits,
    describe_orca_installation,
    discover_orca_installations,
    prepare_bwunicluster_user_setup,
    run_bwunicluster_installer,
    get_packaged_submit_templates_dir,
    get_user_qm_tools_dir,
    get_user_csp_tools_dir,
    get_user_mlp_tools_dir,
    rebuild_venv,
    run_pip_install_editable,
    run_qm_tools_installer,
    run_csp_tools_installer,
    run_mlp_tools_installer,
    run_analysis_tools_installer,
    get_user_analysis_tools_dir,
    stage_packaged_analysis_tools,
    resolve_backend_choice,
    resolve_orca_base,
    resolve_submit_templates_dir,
    stage_packaged_qm_tools,
    stage_packaged_csp_tools,
    stage_packaged_mlp_tools,
)
from delfin.qm_runtime import canonical_tool_name, discover_tool_installations, settings_selectable_tools
from delfin.user_settings import (
    get_settings_path,
    load_settings,
    normalize_choice_setting,
    normalize_local_directory_setting,
    normalize_positive_int_setting,
    normalize_ssh_transfer_settings,
    save_settings,
)


def create_tab(ctx, calc_refs=None, archive_refs=None):
    """Create the dashboard Settings tab."""
    settings_path = get_settings_path()
    status_html = widgets.HTML(value='')
    runtime_diagnostics_html = widgets.HTML(value='')
    qm_tools_log = widgets.Textarea(
        value='',
        disabled=True,
        layout=widgets.Layout(width='100%', height='160px'),
    )
    show_sensitive_btn = widgets.ToggleButton(
        value=False,
        description='Show Values',
        icon='eye',
        layout=widgets.Layout(width='118px', height='28px'),
    )
    reload_btn = widgets.Button(
        description='Reload',
        layout=widgets.Layout(width='90px', height='28px'),
    )
    validate_runtime_btn = widgets.Button(
        description='Validate Setup',
        button_style='info',
        layout=widgets.Layout(width='120px', height='28px'),
    )
    prepare_qm_tools_btn = widgets.Button(
        description='Prepare qm_tools',
        layout=widgets.Layout(width='135px', height='28px'),
    )
    install_qm_tools_btn = widgets.Button(
        description='Install qm_tools',
        button_style='warning',
        layout=widgets.Layout(width='125px', height='28px'),
    )
    update_qm_tools_btn = widgets.Button(
        description='Update qm_tools',
        button_style='info',
        layout=widgets.Layout(width='130px', height='28px'),
    )
    # Individual qm_tool install buttons
    _qm_tool_btn_layout = widgets.Layout(width='100px', height='26px')
    install_xtb_btn = widgets.Button(description='xtb', button_style='warning', layout=_qm_tool_btn_layout)
    install_crest_btn = widgets.Button(description='crest', button_style='warning', layout=_qm_tool_btn_layout)
    install_dftbplus_btn = widgets.Button(description='dftb+', button_style='warning', layout=_qm_tool_btn_layout)
    install_stda_btn = widgets.Button(description='xtb4stda', button_style='warning', layout=_qm_tool_btn_layout)
    install_std2_btn = widgets.Button(description='std2', button_style='warning', layout=_qm_tool_btn_layout)
    install_micromamba_btn = widgets.Button(
        description='Install micromamba',
        button_style='',
        layout=widgets.Layout(width='145px', height='26px'),
    )
    csp_tools_log = widgets.Textarea(
        value='',
        disabled=True,
        layout=widgets.Layout(width='100%', height='160px'),
    )
    csp_tools_root_input = widgets.Text(
        placeholder='Optional csp_tools root override',
        layout=widgets.Layout(width='100%', min_width='280px', height='28px'),
    )
    install_csp_tools_btn = widgets.Button(
        description='Install csp_tools',
        button_style='warning',
        layout=widgets.Layout(width='130px', height='28px'),
    )
    update_csp_tools_btn = widgets.Button(
        description='Update csp_tools',
        button_style='info',
        layout=widgets.Layout(width='130px', height='28px'),
    )
    mlp_tools_log = widgets.Textarea(
        value='',
        disabled=True,
        layout=widgets.Layout(width='100%', height='160px'),
    )
    mlp_tools_root_input = widgets.Text(
        placeholder='Optional mlp_tools root override',
        layout=widgets.Layout(width='100%', min_width='280px', height='28px'),
    )
    install_mlp_tools_btn = widgets.Button(
        description='Install mlp_tools',
        button_style='warning',
        layout=widgets.Layout(width='130px', height='28px'),
    )
    update_mlp_tools_btn = widgets.Button(
        description='Update mlp_tools',
        button_style='info',
        layout=widgets.Layout(width='130px', height='28px'),
    )
    analysis_tools_log = widgets.Textarea(
        value='',
        disabled=True,
        layout=widgets.Layout(width='100%', height='160px'),
    )
    install_analysis_tools_btn = widgets.Button(
        description='Install analysis_tools',
        button_style='warning',
        layout=widgets.Layout(width='160px', height='28px'),
    )
    update_analysis_tools_btn = widgets.Button(
        description='Update analysis_tools',
        button_style='info',
        layout=widgets.Layout(width='160px', height='28px'),
    )
    qm_status_box = widgets.VBox(layout=widgets.Layout(width='100%'))
    refresh_qm_status_btn = widgets.Button(
        description='Refresh QM status',
        button_style='info',
        layout=widgets.Layout(width='155px', height='28px'),
    )
    csp_status_box = widgets.VBox(layout=widgets.Layout(width='100%'))
    refresh_csp_status_btn = widgets.Button(
        description='Refresh CSP status',
        button_style='info',
        layout=widgets.Layout(width='155px', height='28px'),
    )
    analysis_status_html = widgets.HTML(value='')
    analysis_status_box = widgets.VBox(layout=widgets.Layout(width='100%'))
    refresh_analysis_status_btn = widgets.Button(
        description='Refresh status',
        button_style='info',
        layout=widgets.Layout(width='130px', height='28px'),
    )
    mlp_status_html = widgets.HTML(value='')
    mlp_status_box = widgets.VBox(layout=widgets.Layout(width='100%'))
    refresh_mlp_status_btn = widgets.Button(
        description='Refresh MLP status',
        button_style='info',
        layout=widgets.Layout(width='155px', height='28px'),
    )
    ai_status_html = widgets.HTML(value='')
    ai_status_box = widgets.VBox(layout=widgets.Layout(width='100%'))
    refresh_ai_status_btn = widgets.Button(
        description='Refresh AI status',
        button_style='info',
        layout=widgets.Layout(width='155px', height='28px'),
    )
    tool_install_log = widgets.Textarea(
        value='',
        disabled=True,
        layout=widgets.Layout(width='100%', height='120px'),
    )
    pip_install_log = widgets.Textarea(
        value='',
        disabled=True,
        layout=widgets.Layout(width='100%', height='160px'),
    )
    pip_install_btn = widgets.Button(
        description='pip install -e .',
        button_style='warning',
        layout=widgets.Layout(width='140px', height='28px'),
    )
    rebuild_venv_btn = widgets.Button(
        description='Rebuild .venv',
        button_style='danger',
        layout=widgets.Layout(width='140px', height='28px'),
    )
    detect_local_resources_btn = widgets.Button(
        description='Detect local resources',
        layout=widgets.Layout(width='165px', height='28px'),
    )
    setup_bwunicluster_btn = widgets.Button(
        description='Setup bwUniCluster',
        button_style='success',
        layout=widgets.Layout(width='150px', height='28px'),
    )
    verify_bwunicluster_btn = widgets.Button(
        description='Verify bwUniCluster',
        button_style='info',
        layout=widgets.Layout(width='155px', height='28px'),
    )
    full_install_bwunicluster_btn = widgets.Button(
        description='Full bwUni install',
        button_style='warning',
        layout=widgets.Layout(width='145px', height='28px'),
    )
    save_btn = widgets.Button(
        description='Save Settings',
        button_style='primary',
        disabled=True,
        layout=widgets.Layout(width='120px', height='28px'),
    )
    tabs_status_html = widgets.HTML(value='')
    tabs_rows_box = widgets.VBox(layout=widgets.Layout(width='100%', gap='6px'))
    remote_archive_toggle = widgets.Checkbox(
        value=False,
        description='Enabled',
        indent=False,
        layout=widgets.Layout(width='110px', height='28px'),
    )
    calc_path_input = widgets.Text(
        placeholder=str(ctx.default_calc_dir),
        layout=widgets.Layout(width='100%', min_width='280px', height='28px'),
    )
    archive_path_input = widgets.Text(
        placeholder=str(ctx.default_archive_dir),
        layout=widgets.Layout(width='100%', min_width='280px', height='28px'),
    )
    backend_dropdown = widgets.Dropdown(
        options=[
            ('Auto', 'auto'),
            ('Local', 'local'),
            ('SLURM', 'slurm'),
        ],
        value='auto',
        layout=widgets.Layout(width='160px', height='28px'),
    )
    global_orca_input = widgets.Text(
        placeholder='Prefer PATH / auto-detect',
        layout=widgets.Layout(width='100%', min_width='280px', height='28px'),
    )
    scan_orca_btn = widgets.Button(
        description='Scan ORCA',
        layout=widgets.Layout(width='105px', height='28px'),
    )
    detected_orca_dropdown = widgets.Dropdown(
        options=[('Auto-detect / PATH', '')],
        value='',
        disabled=True,
        layout=widgets.Layout(width='320px', min_width='250px', height='28px'),
    )
    _tool_display_names = {
        'gaussian': 'Gaussian',
        'turbomole': 'Turbomole',
        'xtb': 'xTB',
        'crest': 'CREST',
        'std2': 'STD2',
        'stda': 'STDA',
        'xtb4stda': 'xtb4stda',
        'dftb+': 'DFTB+',
    }
    _selectable_tool_names = tuple(
        tool_name for tool_name in settings_selectable_tools() if tool_name in _tool_display_names
    )
    tool_binary_inputs = {}
    tool_scan_buttons = {}
    tool_detected_dropdowns = {}
    for _tool_name in _selectable_tool_names:
        tool_binary_inputs[_tool_name] = widgets.Text(
            placeholder='Prefer PATH / auto-detect',
            layout=widgets.Layout(width='100%', min_width='280px', height='28px'),
        )
        tool_scan_buttons[_tool_name] = widgets.Button(
            description='Scan',
            layout=widgets.Layout(width='80px', height='28px'),
        )
        tool_detected_dropdowns[_tool_name] = widgets.Dropdown(
            options=[('Auto-detect / PATH', '')],
            value='',
            disabled=True,
            layout=widgets.Layout(width='320px', min_width='250px', height='28px'),
        )
    qm_tools_root_input = widgets.Text(
        placeholder='Optional qm_tools root override',
        layout=widgets.Layout(width='100%', min_width='280px', height='28px'),
    )
    local_orca_input = widgets.Text(
        placeholder='Optional local ORCA override',
        layout=widgets.Layout(width='100%', min_width='280px', height='28px'),
    )
    _init_cores, _init_ram_mb = detect_local_runtime_limits()
    local_max_cores_input = widgets.BoundedIntText(
        value=_init_cores,
        min=1,
        max=1_000_000,
        step=1,
        layout=widgets.Layout(width='120px', min_width='120px', height='28px'),
    )
    local_max_ram_input = widgets.BoundedIntText(
        value=_init_ram_mb,
        min=1,
        max=100_000_000,
        step=1000,
        layout=widgets.Layout(width='140px', min_width='140px', height='28px'),
    )
    slurm_orca_input = widgets.Text(
        placeholder='Optional SLURM ORCA override',
        layout=widgets.Layout(width='100%', min_width='280px', height='28px'),
    )
    slurm_templates_input = widgets.Text(
        placeholder='Optional submit template directory',
        layout=widgets.Layout(width='100%', min_width='280px', height='28px'),
    )
    slurm_profile_input = widgets.Text(
        placeholder='Optional site profile label',
        layout=widgets.Layout(width='220px', min_width='200px', height='28px'),
    )
    state = {
        'remote_archive_enabled': False,
        'calculations_dir': str(ctx.calc_dir),
        'archive_dir': str(ctx.archive_dir),
        'tab_prefs': {
            'order': [],
            'hidden': [],
        },
    }

    host_hidden = widgets.Password(
        placeholder='Host or IP',
        layout=widgets.Layout(width='230px', min_width='180px', height='28px'),
    )
    host_visible = widgets.Text(
        placeholder='Host or IP',
        layout=widgets.Layout(
            width='230px',
            min_width='180px',
            height='28px',
            display='none',
        ),
    )
    user_hidden = widgets.Password(
        placeholder='Remote user',
        layout=widgets.Layout(width='180px', min_width='150px', height='28px'),
    )
    user_visible = widgets.Text(
        placeholder='Remote user',
        layout=widgets.Layout(
            width='180px',
            min_width='150px',
            height='28px',
            display='none',
        ),
    )
    path_hidden = widgets.Password(
        placeholder='Remote target folder',
        layout=widgets.Layout(width='100%', min_width='280px', height='28px'),
    )
    path_visible = widgets.Text(
        placeholder='Remote target folder',
        layout=widgets.Layout(
            width='100%',
            min_width='280px',
            height='28px',
            display='none',
        ),
    )
    port_input = widgets.BoundedIntText(
        value=22,
        min=1,
        max=65535,
        step=1,
        layout=widgets.Layout(width='96px', min_width='96px', height='28px'),
    )

    def _set_status(message, color='#455a64'):
        status_html.value = f'<span style="color:{color};">{message}</span>'

    def _slurm_available():
        return shutil.which('sbatch') is not None

    def _runtime_auto_candidates():
        candidates = []
        if ctx.orca_base:
            candidates.append(str(ctx.orca_base))
        for candidate in ctx.orca_candidates or []:
            candidate_text = str(candidate)
            if candidate_text not in candidates:
                candidates.append(candidate_text)
        return candidates

    def _orca_search_roots():
        roots = [Path.home() / 'software', Path.home() / 'apps', Path.home() / 'local', Path('/opt')]
        if getattr(ctx, 'repo_dir', None):
            roots.append(Path(ctx.repo_dir).parent / 'software')
        return roots

    _path_source_labels = {
        'auto': 'Auto',
        'local': 'Local',
        'cluster': 'System',
        'system': 'System',
        'external': 'External',
    }
    _turbomole_program_names = ('ridft', 'dscf', 'define')

    def _settings_source_label_from_path(path_value: str) -> str:
        source = describe_installation_source(path_value)
        return _path_source_labels.get(source, str(source or '').strip().title() or 'External')

    def _module_name_from_source(source: str) -> str:
        normalized_source = str(source or '').strip()
        if normalized_source.startswith('module:'):
            return normalized_source.split(':', 1)[1]
        return ''

    def _candidate_ancestors(path_value: str):
        text = str(path_value or '').strip()
        if not text:
            return []
        path = Path(text).expanduser()
        try:
            path = path.resolve()
        except Exception:
            pass
        return list(path.parents)

    def _infer_turbomole_install_root(path_value: str):
        for ancestor in _candidate_ancestors(path_value):
            if ancestor.name == 'bin':
                return ancestor.parent
        for ancestor in _candidate_ancestors(path_value):
            bin_dir = ancestor / 'bin'
            if not bin_dir.is_dir():
                continue
            if any((bin_dir / name).is_file() for name in _turbomole_program_names):
                return ancestor
            try:
                subdirs = [item for item in bin_dir.iterdir() if item.is_dir()]
            except Exception:
                subdirs = []
            for subdir in subdirs:
                if any((subdir / name).is_file() for name in _turbomole_program_names):
                    return ancestor
        return None

    def _summarize_turbomole_install(path_value: str, source: str) -> str:
        module_name = _module_name_from_source(source)
        if module_name:
            return module_name
        root = _infer_turbomole_install_root(path_value)
        if root is None:
            return Path(path_value).name or str(path_value)
        root_name = root.name.strip()
        if root_name and root_name.lower() not in {'turbomole', 'tmole'}:
            return root_name
        parent_name = root.parent.name.strip()
        if parent_name and parent_name.lower() not in {'chem', 'common', 'opt'}:
            return f'{parent_name}/{root_name}'
        return str(root)

    def _prefer_turbomole_candidate(current, candidate):
        priority = {'ridft': 0, 'dscf': 1, 'define': 2}
        current_rank = priority.get(Path(current.path).name.lower(), 99)
        candidate_rank = priority.get(Path(candidate.path).name.lower(), 99)
        return candidate if candidate_rank < current_rank else current

    def _dropdown_resolved_tools(tool_name):
        resolved_tools = discover_tool_installations(tool_name)
        if tool_name != 'turbomole':
            return resolved_tools

        grouped = {}
        ordered_keys = []
        for resolved in resolved_tools:
            install_root = _infer_turbomole_install_root(resolved.path)
            group_key = str(install_root) if install_root is not None else str(Path(resolved.path).expanduser())
            if group_key not in grouped:
                grouped[group_key] = resolved
                ordered_keys.append(group_key)
                continue
            grouped[group_key] = _prefer_turbomole_candidate(grouped[group_key], resolved)

        return [grouped[key] for key in ordered_keys]

    def _format_orca_option_label(path_value, *, configured=False):
        source = _settings_source_label_from_path(path_value)
        base_label = describe_orca_installation(path_value)
        if configured:
            return f'{base_label} [{source}, configured]'
        return f'{base_label} [{source}]'

    def _refresh_orca_dropdown(selected_value=None):
        current_value = str(selected_value if selected_value is not None else global_orca_input.value or '').strip()
        options = [('Auto-detect / PATH', '')]
        for candidate in discover_orca_installations(
            seed_candidates=_runtime_auto_candidates(),
            search_roots=_orca_search_roots(),
        ):
            options.append((_format_orca_option_label(candidate), candidate))

        if current_value and all(current_value != value for _label, value in options):
            options.append((_format_orca_option_label(current_value, configured=True), current_value))

        detected_orca_dropdown.options = options
        allowed_values = {value for _label, value in options}
        detected_orca_dropdown.value = current_value if current_value in allowed_values else ''
        detected_orca_dropdown.disabled = len(options) <= 1

    def _tool_source_label(source: str, path_value: str) -> str:
        if _module_name_from_source(source):
            return 'Module'
        normalized_source = str(source or '').strip()
        if normalized_source == 'qm_tools':
            return 'qm_tools'
        if normalized_source.startswith('env:'):
            return 'Environment'
        return _settings_source_label_from_path(path_value)

    def _format_tool_option_label(tool_name, path_value, source, *, configured=False):
        display_name = _tool_display_names.get(tool_name, tool_name.upper())
        module_name = _module_name_from_source(source)
        if tool_name == 'turbomole':
            label = f'{display_name} ({_summarize_turbomole_install(path_value, source)})'
        else:
            path_obj = Path(path_value)
            label = f'{display_name} ({path_obj.name})'
            if module_name:
                label += f' @ {module_name}'
            else:
                parent_name = path_obj.parent.name
                if parent_name and parent_name not in {'bin', path_obj.name}:
                    label += f' @ {parent_name}'
        source_label = _tool_source_label(source, path_value)
        if configured:
            source_label = f'{source_label}, configured'
        return f'{label} [{source_label}]'

    def _refresh_tool_dropdown(tool_name, selected_value=None):
        current_value = str(
            selected_value if selected_value is not None else tool_binary_inputs[tool_name].value or ''
        ).strip()
        options = [('Auto-detect / PATH', '')]
        for resolved in _dropdown_resolved_tools(tool_name):
            options.append((
                _format_tool_option_label(tool_name, resolved.path, resolved.source),
                resolved.path,
            ))

        if current_value and all(current_value != value for _label, value in options):
            options.append((
                _format_tool_option_label(tool_name, current_value, 'explicit', configured=True),
                current_value,
            ))

        dropdown = tool_detected_dropdowns[tool_name]
        dropdown.options = options
        allowed_values = {value for _label, value in options}
        dropdown.value = current_value if current_value in allowed_values else ''
        dropdown.disabled = len(options) <= 1

    def _fallback_submit_templates_dir():
        current = getattr(ctx, 'submit_templates_dir', None)
        if current:
            return Path(current)
        return get_packaged_submit_templates_dir()

    def _resolve_runtime_effective_values(runtime_payload):
        effective_backend = resolve_backend_choice(
            None,
            runtime_payload.get('backend', 'auto'),
            slurm_available=_slurm_available(),
        )
        effective_orca_base = resolve_orca_base(
            None,
            runtime_payload,
            effective_backend,
            auto_candidates=_runtime_auto_candidates(),
            local_default='/opt/orca' if effective_backend == 'local' else '',
        )
        submit_templates_dir = resolve_submit_templates_dir(
            runtime_payload,
            _fallback_submit_templates_dir(),
        )
        return effective_backend, effective_orca_base, submit_templates_dir

    def _render_runtime_diagnostics(runtime_payload, *, reload_required=False):
        effective_backend, effective_orca_base, submit_templates_dir = _resolve_runtime_effective_values(
            runtime_payload
        )
        diagnostics = collect_runtime_diagnostics(
            runtime_payload,
            backend=effective_backend,
            effective_orca_base=effective_orca_base,
            submit_templates_dir=submit_templates_dir if effective_backend == 'slurm' else None,
        )

        summary_parts = [
            f'Effective backend: <code>{html.escape(effective_backend)}</code>',
            (
                f'ORCA: <code>{html.escape(effective_orca_base)}</code>'
                if effective_orca_base
                else 'ORCA: <code>PATH / auto-detect</code>'
            ),
        ]
        if runtime_payload.get('qm_tools_root'):
            summary_parts.append(
                f'qm_tools: <code>{html.escape(str(runtime_payload.get("qm_tools_root", "")))}</code>'
            )
        if runtime_payload.get('csp_tools_root'):
            summary_parts.append(
                f'csp_tools: <code>{html.escape(str(runtime_payload.get("csp_tools_root", "")))}</code>'
            )
        if runtime_payload.get('mlp_tools_root'):
            summary_parts.append(
                f'mlp_tools: <code>{html.escape(str(runtime_payload.get("mlp_tools_root", "")))}</code>'
            )
        if effective_backend == 'slurm':
            summary_parts.append(
                f'Submit templates: <code>{html.escape(str(submit_templates_dir))}</code>'
            )
        if reload_required:
            summary_parts.append(
                f'<span style="color:#ef6c00;">Reload required to switch from '
                f'<code>{html.escape(ctx.runtime_backend)}</code> to '
                f'<code>{html.escape(effective_backend)}</code>.</span>'
            )

        rows = []
        for item in diagnostics:
            status = item.get('status', 'missing')
            if status == 'ok':
                status_color, status_label = '#2e7d32', 'OK'
            elif status == 'system':
                status_color, status_label = '#1565c0', 'System'
            elif status == 'module':
                status_color, status_label = '#e65100', 'Module'
            else:
                status_color, status_label = '#d32f2f', 'Missing'
            rows.append(
                '<tr>'
                f'<td style="padding:4px 8px;"><code>{html.escape(str(item.get("name", "")))}</code></td>'
                f'<td style="padding:4px 8px; color:{status_color}; font-weight:600;">{status_label}</td>'
                f'<td style="padding:4px 8px;"><code>{html.escape(str(item.get("detail", "")))}</code></td>'
                '</tr>'
            )

        runtime_diagnostics_html.value = (
            '<div style="border:1px solid #d9dee3; border-radius:6px; padding:10px; background:#fafbfc;">'
            '<b>Runtime validation</b><br>'
            + ' | '.join(summary_parts)
            + '<table style="margin-top:8px; border-collapse:collapse; width:100%;">'
            '<thead>'
            '<tr>'
            '<th style="text-align:left; padding:4px 8px;">Check</th>'
            '<th style="text-align:left; padding:4px 8px;">Status</th>'
            '<th style="text-align:left; padding:4px 8px;">Detail</th>'
            '</tr>'
            '</thead>'
            '<tbody>'
            + ''.join(rows)
            + '</tbody></table>'
            + '<div style="margin-top:8px; color:#546e7a; font-size:12px;">'
            + 'System = found system-wide. Module = available via Environment Modules.'
            + '</div></div>'
        )
        return effective_backend, effective_orca_base, submit_templates_dir

    def _persist_runtime_payload(runtime_payload):
        settings_payload = load_settings()
        settings_payload['runtime'] = runtime_payload
        settings_payload = save_settings(settings_payload, settings_path)
        _set_runtime_widgets(settings_payload)
        backend_switch_required, effective_backend, effective_orca_base = _apply_runtime_settings(
            settings_payload.get('runtime', {}) or runtime_payload
        )
        _render_runtime_diagnostics(
            settings_payload.get('runtime', {}) or runtime_payload,
            reload_required=backend_switch_required,
        )
        return backend_switch_required, effective_backend, effective_orca_base

    def _set_transfer_widgets(payload):
        host_value = str((payload or {}).get('host') or '')
        user_value = str((payload or {}).get('user') or '')
        path_value = str((payload or {}).get('remote_path') or '')
        port_value = int((payload or {}).get('port') or 22)

        host_hidden.value = host_value
        user_hidden.value = user_value
        path_hidden.value = path_value
        port_input.value = port_value
        if show_sensitive_btn.value:
            host_visible.value = host_value
            user_visible.value = user_value
            path_visible.value = path_value
        else:
            host_visible.value = ''
            user_visible.value = ''
            path_visible.value = ''

    def _set_path_widgets(settings_payload):
        paths_payload = ((settings_payload or {}).get('paths') or {})
        calc_path_input.value = str(paths_payload.get('calculations_dir') or '')
        archive_path_input.value = str(paths_payload.get('archive_dir') or '')

    def _set_runtime_widgets(settings_payload):
        detected_local_cores, detected_local_ram_mb = detect_local_runtime_limits()
        runtime_payload = ((settings_payload or {}).get('runtime') or {})
        local_payload = runtime_payload.get('local', {}) or {}
        slurm_payload = runtime_payload.get('slurm', {}) or {}

        backend_dropdown.value = runtime_payload.get('backend', 'auto')
        global_orca_input.value = str(runtime_payload.get('orca_base') or '')
        qm_tools_root_input.value = str(runtime_payload.get('qm_tools_root') or '')
        csp_tools_root_input.value = str(runtime_payload.get('csp_tools_root') or '')
        mlp_tools_root_input.value = str(runtime_payload.get('mlp_tools_root') or '')
        local_orca_input.value = str(local_payload.get('orca_base') or '')
        local_max_cores_input.value = int(local_payload.get('max_cores', detected_local_cores))
        local_max_ram_input.value = int(local_payload.get('max_ram_mb', detected_local_ram_mb))
        slurm_orca_input.value = str(slurm_payload.get('orca_base') or '')
        slurm_templates_input.value = str(slurm_payload.get('submit_templates_dir') or '')
        slurm_profile_input.value = str(slurm_payload.get('profile') or '')
        tool_binaries_payload = runtime_payload.get('tool_binaries', {}) or {}
        for tool_name in _selectable_tool_names:
            tool_binary_inputs[tool_name].value = str(tool_binaries_payload.get(tool_name) or '')
            _refresh_tool_dropdown(tool_name, tool_binary_inputs[tool_name].value)
        _refresh_orca_dropdown(global_orca_input.value)

    def _effective_paths_from_widgets():
        calc_override = normalize_local_directory_setting(
            calc_path_input.value,
            'Calculations path',
        )
        archive_override = normalize_local_directory_setting(
            archive_path_input.value,
            'Archive path',
        )
        effective_calc_dir = Path(calc_override) if calc_override else Path(ctx.default_calc_dir)
        effective_archive_dir = (
            Path(archive_override) if archive_override else Path(ctx.default_archive_dir)
        )
        return calc_override, archive_override, effective_calc_dir, effective_archive_dir

    def _runtime_payload_from_widgets():
        tool_binaries = {}
        for tool_name in _selectable_tool_names:
            configured_path = normalize_local_directory_setting(
                tool_binary_inputs[tool_name].value,
                f'{_tool_display_names.get(tool_name, tool_name)} binary',
            )
            if configured_path:
                tool_binaries[canonical_tool_name(tool_name)] = configured_path
        return {
            'backend': normalize_choice_setting(
                backend_dropdown.value,
                'Runtime backend',
                {'auto', 'local', 'slurm'},
                'auto',
            ),
            'orca_base': normalize_local_directory_setting(
                global_orca_input.value,
                'Global ORCA path',
            ),
            'qm_tools_root': normalize_local_directory_setting(
                qm_tools_root_input.value,
                'qm_tools root',
            ),
            'csp_tools_root': normalize_local_directory_setting(
                csp_tools_root_input.value,
                'csp_tools root',
            ),
            'mlp_tools_root': normalize_local_directory_setting(
                mlp_tools_root_input.value,
                'mlp_tools root',
            ),
            'tool_binaries': tool_binaries,
            'local': {
                'orca_base': normalize_local_directory_setting(
                    local_orca_input.value,
                    'Local ORCA path',
                ),
                'max_cores': normalize_positive_int_setting(
                    local_max_cores_input.value,
                    'Local max cores',
                    384,
                ),
                'max_ram_mb': normalize_positive_int_setting(
                    local_max_ram_input.value,
                    'Local max RAM (MB)',
                    1_400_000,
                ),
            },
            'slurm': {
                'orca_base': normalize_local_directory_setting(
                    slurm_orca_input.value,
                    'SLURM ORCA path',
                ),
                'submit_templates_dir': normalize_local_directory_setting(
                    slurm_templates_input.value,
                    'SLURM submit templates path',
                ),
                'profile': str(slurm_profile_input.value or '').strip(),
            },
        }

    def _apply_workspace_paths(effective_calc_dir, effective_archive_dir):
        effective_calc_dir.mkdir(parents=True, exist_ok=True)
        effective_archive_dir.mkdir(parents=True, exist_ok=True)
        ctx.calc_dir = effective_calc_dir
        ctx.archive_dir = effective_archive_dir
        state['calculations_dir'] = str(effective_calc_dir)
        state['archive_dir'] = str(effective_archive_dir)

        if calc_refs and callable(calc_refs.get('calc_set_root')):
            calc_refs['calc_set_root'](effective_calc_dir)
        if archive_refs and callable(archive_refs.get('calc_set_root')):
            archive_refs['calc_set_root'](effective_archive_dir)
        if archive_refs and callable(archive_refs.get('calc_set_primary_root')):
            archive_refs['calc_set_primary_root'](effective_calc_dir)

    def _apply_runtime_settings(runtime_payload):
        effective_backend, effective_orca_base, submit_templates_dir = _resolve_runtime_effective_values(
            runtime_payload
        )
        backend_switch_required = effective_backend != ctx.runtime_backend

        ctx.runtime_settings = runtime_payload
        apply_runtime_environment(
            qm_tools_root=runtime_payload.get('qm_tools_root', ''),
            orca_base=effective_orca_base,
            csp_tools_root=runtime_payload.get('csp_tools_root', ''),
            tool_binaries=runtime_payload.get('tool_binaries', {}) or {},
        )

        if not backend_switch_required:
            ctx.orca_base = effective_orca_base
            if getattr(ctx, 'backend', None) is not None and hasattr(ctx.backend, 'orca_base'):
                ctx.backend.orca_base = effective_orca_base
            if getattr(ctx, 'backend', None) is not None and hasattr(ctx.backend, 'tool_binaries'):
                ctx.backend.tool_binaries = dict(runtime_payload.get('tool_binaries', {}) or {})
            if effective_backend == 'local':
                local_settings = runtime_payload.get('local', {}) or {}
                if hasattr(ctx.backend, 'max_cores'):
                    ctx.backend.max_cores = int(local_settings.get('max_cores', ctx.backend.max_cores))
                if hasattr(ctx.backend, 'max_ram_mb'):
                    ctx.backend.max_ram_mb = int(local_settings.get('max_ram_mb', ctx.backend.max_ram_mb))
            if effective_backend == 'slurm':
                ctx.submit_templates_dir = Path(submit_templates_dir)
                if hasattr(ctx.backend, 'submit_templates_dir'):
                    ctx.backend.submit_templates_dir = Path(submit_templates_dir)

        return backend_switch_required, effective_backend, effective_orca_base

    def _transfer_payload_from_widgets():
        if show_sensitive_btn.value:
            host_value = host_visible.value
            user_value = user_visible.value
            path_value = path_visible.value
            host_hidden.value = host_value
            user_hidden.value = user_value
            path_hidden.value = path_value
            return host_value, user_value, path_value, port_input.value
        return host_hidden.value, user_hidden.value, path_hidden.value, port_input.value

    def _apply_sensitive_visibility():
        visible = bool(show_sensitive_btn.value)
        if visible:
            host_visible.value = host_hidden.value
            user_visible.value = user_hidden.value
            path_visible.value = path_hidden.value
        else:
            host_hidden.value = host_visible.value
            user_hidden.value = user_visible.value
            path_hidden.value = path_visible.value
            host_visible.value = ''
            user_visible.value = ''
            path_visible.value = ''

        host_hidden.layout.display = 'none' if visible else 'inline-flex'
        user_hidden.layout.display = 'none' if visible else 'inline-flex'
        path_hidden.layout.display = 'none' if visible else 'inline-flex'
        host_visible.layout.display = 'inline-flex' if visible else 'none'
        user_visible.layout.display = 'inline-flex' if visible else 'none'
        path_visible.layout.display = 'inline-flex' if visible else 'none'
        show_sensitive_btn.icon = 'eye-slash' if visible else 'eye'
        show_sensitive_btn.description = 'Hide Values' if visible else 'Show Values'

    def _set_remote_archive_widget(settings_payload):
        features = ((settings_payload or {}).get('features') or {})
        enabled = bool(features.get('remote_archive_enabled', False))
        remote_archive_toggle.value = enabled
        state['remote_archive_enabled'] = enabled

    def _available_tab_specs():
        return list(getattr(ctx, 'tab_specs', []) or [])

    def _tab_pref_lookup():
        prefs = state.get('tab_prefs', {}) or {}
        order = [str(item) for item in (prefs.get('order', []) or []) if str(item).strip()]
        hidden = {str(item) for item in (prefs.get('hidden', []) or []) if str(item).strip()}
        return order, hidden

    def _normalized_tab_prefs():
        specs = _available_tab_specs()
        spec_map = {spec.get('id'): spec for spec in specs}
        order, hidden = _tab_pref_lookup()
        movable_defaults = [
            spec['id']
            for spec in sorted(specs, key=lambda item: item.get('default_order', 10_000))
            if not spec.get('fixed')
        ]
        normalized_order = []
        seen = set()
        for tab_id in order + movable_defaults:
            if tab_id in seen or tab_id not in spec_map or spec_map[tab_id].get('fixed'):
                continue
            seen.add(tab_id)
            normalized_order.append(tab_id)
        saved_order = set(order)
        normalized_hidden = [
            tab_id
            for tab_id in normalized_order
            if (
                (
                    tab_id in hidden
                    if tab_id in saved_order
                    else not bool(spec_map[tab_id].get('default_visible', True))
                )
                and not spec_map[tab_id].get('fixed')
            )
        ]
        return {
            'order': normalized_order,
            'hidden': normalized_hidden,
        }

    def _set_tab_preferences(settings_payload):
        ui_payload = ((settings_payload or {}).get('ui') or {})
        tabs_payload = (ui_payload.get('tabs') or {})
        state['tab_prefs'] = {
            'order': list(tabs_payload.get('order', []) or []),
            'hidden': list(tabs_payload.get('hidden', []) or []),
        }
        _rebuild_tab_rows()

    def _move_tab(tab_id, direction):
        prefs = _normalized_tab_prefs()
        order = list(prefs['order'])
        try:
            idx = order.index(tab_id)
        except ValueError:
            return
        new_idx = idx + int(direction)
        if new_idx < 0 or new_idx >= len(order):
            return
        order[idx], order[new_idx] = order[new_idx], order[idx]
        state['tab_prefs'] = {
            'order': order,
            'hidden': list(prefs['hidden']),
        }
        _rebuild_tab_rows()
        _mark_dirty()

    def _toggle_tab_hidden(tab_id, visible):
        prefs = _normalized_tab_prefs()
        hidden = set(prefs['hidden'])
        if visible:
            hidden.discard(tab_id)
        else:
            hidden.add(tab_id)
        state['tab_prefs'] = {
            'order': list(prefs['order']),
            'hidden': [item for item in prefs['order'] if item in hidden],
        }
        _rebuild_tab_rows()
        _mark_dirty()

    def _rebuild_tab_rows():
        specs = _available_tab_specs()
        if not specs:
            tabs_rows_box.children = (
                widgets.HTML('<span style="color:#616161;">No dashboard tabs registered.</span>'),
            )
            tabs_status_html.value = ''
            return

        prefs = _normalized_tab_prefs()
        state['tab_prefs'] = prefs
        spec_map = {spec['id']: spec for spec in specs}
        rows = []
        for idx, tab_id in enumerate(prefs['order']):
            spec = spec_map.get(tab_id)
            if not spec:
                continue
            is_hidden = tab_id in set(prefs['hidden'])
            available = bool(spec.get('available'))
            fixed = bool(spec.get('fixed'))
            visible_checkbox = widgets.Checkbox(
                value=(not is_hidden) if not fixed else True,
                description='Visible',
                indent=False,
                disabled=fixed or False,
                layout=widgets.Layout(width='85px'),
            )
            up_btn = widgets.Button(
                description='Up',
                layout=widgets.Layout(width='56px', height='26px'),
                disabled=fixed or idx == 0,
            )
            down_btn = widgets.Button(
                description='Down',
                layout=widgets.Layout(width='64px', height='26px'),
                disabled=fixed or idx == len(prefs['order']) - 1,
            )
            badge_parts = []
            if fixed:
                badge_parts.append('<span style="color:#1565c0;">fixed</span>')
            elif is_hidden:
                badge_parts.append('<span style="color:#ef6c00;">hidden</span>')
            else:
                badge_parts.append('<span style="color:#2e7d32;">visible</span>')
            if not available:
                badge_parts.append('<span style="color:#9e9e9e;">unavailable</span>')
            detail = spec.get('reason', '')
            label_html = (
                f'<b>{html.escape(str(spec.get("title", tab_id)))}</b> '
                + ' | '.join(badge_parts)
            )
            if detail:
                label_html += (
                    f'<br><span style="color:#616161; font-size:0.92em;">'
                    f'{html.escape(str(detail))}</span>'
                )
            label = widgets.HTML(label_html, layout=widgets.Layout(flex='1 1 auto'))
            visible_checkbox.observe(
                lambda change, _tab_id=tab_id: (
                    _toggle_tab_hidden(_tab_id, bool(change.get('new')))
                    if change.get('name') == 'value'
                    else None
                ),
                names='value',
            )
            up_btn.on_click(lambda _button, _tab_id=tab_id: _move_tab(_tab_id, -1))
            down_btn.on_click(lambda _button, _tab_id=tab_id: _move_tab(_tab_id, 1))
            rows.append(
                widgets.HBox(
                    [label, visible_checkbox, up_btn, down_btn],
                    layout=widgets.Layout(
                        width='100%',
                        gap='8px',
                        align_items='center',
                        border='1px solid #e0e0e0',
                        padding='8px 10px',
                    ),
                )
            )

        fixed_specs = [
            spec for spec in sorted(specs, key=lambda item: item.get('default_order', 10_000))
            if spec.get('fixed')
        ]
        for spec in fixed_specs:
            if spec['id'] in prefs['order']:
                continue
            label = widgets.HTML(
                f'<b>{html.escape(str(spec.get("title", spec["id"])))}</b> '
                '<span style="color:#1565c0;">fixed</span><br>'
                '<span style="color:#616161; font-size:0.92em;">Always visible.</span>',
                layout=widgets.Layout(flex='1 1 auto'),
            )
            rows.append(
                widgets.HBox(
                    [
                        label,
                        widgets.Checkbox(
                            value=True,
                            description='Visible',
                            indent=False,
                            disabled=True,
                            layout=widgets.Layout(width='85px'),
                        ),
                        widgets.Button(
                            description='Up',
                            disabled=True,
                            layout=widgets.Layout(width='56px', height='26px'),
                        ),
                        widgets.Button(
                            description='Down',
                            disabled=True,
                            layout=widgets.Layout(width='64px', height='26px'),
                        ),
                    ],
                    layout=widgets.Layout(
                        width='100%',
                        gap='8px',
                        align_items='center',
                        border='1px solid #e0e0e0',
                        padding='8px 10px',
                    ),
                )
            )

        hidden_count = len(prefs['hidden'])
        tabs_status_html.value = (
            f'<span style="color:#455a64;">'
            f'{len(prefs["order"])} configurable tabs, '
            f'{hidden_count} currently hidden. '
            'Settings stays visible and pinned at the end.'
            '</span>'
        )
        tabs_rows_box.children = tuple(rows)

    def _load_settings_to_widgets(set_status=True):
        try:
            settings_payload = load_settings()
        except Exception as exc:
            _set_transfer_widgets({})
            _set_path_widgets({})
            _set_runtime_widgets({})
            _set_remote_archive_widget({})
            _set_tab_preferences({})
            _set_status(
                (
                    f'Could not load <code>{html.escape(str(settings_path))}</code>: '
                    f'{html.escape(str(exc))}'
                ),
                color='#d32f2f',
            )
            return None

        payload = settings_payload.get('transfer', {}) or {}
        _set_transfer_widgets(payload)
        _set_path_widgets(settings_payload)
        _set_runtime_widgets(settings_payload)
        _set_remote_archive_widget(settings_payload)
        _set_tab_preferences(settings_payload)
        runtime_payload = settings_payload.get('runtime', {}) or {}
        effective_backend, effective_orca_base, _submit_templates_dir = _render_runtime_diagnostics(
            runtime_payload,
            reload_required=False,
        )
        if set_status:
            current_calc_dir = html.escape(str(ctx.calc_dir))
            current_archive_dir = html.escape(str(ctx.archive_dir))
            if payload:
                _set_status(
                    (
                        f'Loaded local settings from '
                        f'<code>{html.escape(str(settings_path))}</code>. '
                        'Sensitive values stay masked until revealed. '
                        f'Current paths: <code>{current_calc_dir}</code> and '
                        f'<code>{current_archive_dir}</code>. '
                        f'Runtime resolves to <code>{html.escape(effective_backend)}</code> '
                        f'with ORCA <code>{html.escape(effective_orca_base or "PATH / auto-detect")}</code>.'
                    ),
                    color='#2e7d32',
                )
            else:
                _set_status(
                    (
                        f'No transfer target saved yet. Settings will be created in '
                        f'<code>{html.escape(str(settings_path))}</code>. '
                        f'Current paths: <code>{current_calc_dir}</code> and '
                        f'<code>{current_archive_dir}</code>. '
                        f'Runtime resolves to <code>{html.escape(effective_backend)}</code>.'
                    ),
                    color='#ef6c00',
                )
        return settings_payload

    _long_running_buttons = [
        install_qm_tools_btn, update_qm_tools_btn,
        install_csp_tools_btn, update_csp_tools_btn,
        install_mlp_tools_btn, update_mlp_tools_btn,
        install_analysis_tools_btn, update_analysis_tools_btn,
        pip_install_btn, rebuild_venv_btn,
        setup_bwunicluster_btn, verify_bwunicluster_btn,
        full_install_bwunicluster_btn,
    ]

    def _with_buttons_disabled(func):
        """Wrapper that disables long-running buttons during execution."""
        def wrapper(button):
            for btn in _long_running_buttons:
                btn.disabled = True
            try:
                func(button)
            finally:
                for btn in _long_running_buttons:
                    btn.disabled = False
        return wrapper

    def _mark_dirty(change=None):
        """Enable the Save button when any setting widget changes."""
        save_btn.disabled = False

    def _on_toggle_sensitive(change):
        if change.get('name') != 'value':
            return
        _apply_sensitive_visibility()

    def _on_reload(button):
        _load_settings_to_widgets(set_status=True)
        save_btn.disabled = True

    def _on_validate_runtime(button):
        try:
            runtime_payload = _runtime_payload_from_widgets()
            effective_backend, _effective_orca_base, _submit_templates_dir = _resolve_runtime_effective_values(
                runtime_payload
            )
            _render_runtime_diagnostics(
                runtime_payload,
                reload_required=(effective_backend != ctx.runtime_backend),
            )
            _set_status(
                'Runtime validation finished. Review the checks below before saving.',
                color='#2e7d32',
            )
        except Exception as exc:
            _set_status(
                f'Runtime validation failed: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _on_scan_orca(button):
        try:
            _refresh_orca_dropdown(global_orca_input.value)
            options_count = max(len(detected_orca_dropdown.options) - 1, 0)
            if options_count:
                _set_status(
                    f'Found {options_count} ORCA installation(s). Select one below or keep auto-detect.',
                    color='#2e7d32',
                )
            else:
                _set_status(
                    'No ORCA installation was found in PATH, configured roots, or detected system locations. Keep auto-detect or enter a path manually.',
                    color='#ef6c00',
                )
        except Exception as exc:
            _set_status(
                f'ORCA scan failed: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _on_detect_local_resources(button):
        try:
            detected_local_cores, detected_local_ram_mb = detect_local_runtime_limits()
            local_max_cores_input.value = detected_local_cores
            local_max_ram_input.value = detected_local_ram_mb
            _set_status(
                (
                    'Detected local resources from the current system: '
                    f'<code>{detected_local_cores}</code> cores and '
                    f'<code>{detected_local_ram_mb}</code> MB RAM.'
                ),
                color='#2e7d32',
            )
        except Exception as exc:
            _set_status(
                f'Local resource detection failed: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _on_select_detected_orca(change):
        if change.get('name') != 'value':
            return
        selected = str(change.get('new') or '')
        global_orca_input.value = selected

    def _on_change_global_orca(change):
        if change.get('name') != 'value':
            return
        _refresh_orca_dropdown(change.get('new'))

    def _on_prepare_qm_tools(button):
        try:
            target = stage_packaged_qm_tools()
            runtime_payload = _runtime_payload_from_widgets()
            runtime_payload['qm_tools_root'] = str(target)
            backend_switch_required, effective_backend, effective_orca_base = _persist_runtime_payload(
                runtime_payload
            )
            qm_tools_log.value = (
                f'Prepared bundled qm_tools under {target}\n'
                'Only the packaged DELFIN bundle was copied. '
                'For optional downloads / conda-based extras use "Install qm_tools" or "Update qm_tools".'
            )
            backend_hint = (
                f' Reload DELFIN to switch execution backend from '
                f'{ctx.runtime_backend} to {effective_backend}.'
                if backend_switch_required
                else ''
            )
            _set_status(
                (
                    f'Bundled qm_tools prepared in <code>{html.escape(str(target))}</code>. '
                    f'Runtime now uses <code>{html.escape(str(target))}</code>. '
                    f'Effective backend: <code>{html.escape(effective_backend)}</code>; '
                    f'ORCA: <code>{html.escape(effective_orca_base or "PATH / auto-detect")}</code>.'
                    f'{backend_hint}'
                ),
                color='#2e7d32',
            )
        except Exception as exc:
            qm_tools_log.value = ''
            _set_status(
                f'Preparing qm_tools failed: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _on_install_qm_tools(button):
        try:
            target, result = run_qm_tools_installer(
                extra_env={'INSTALL_STD2_FROM_SOURCE': '1'},
            )
            runtime_payload = _runtime_payload_from_widgets()
            runtime_payload['qm_tools_root'] = str(target)
            backend_switch_required, effective_backend, effective_orca_base = _persist_runtime_payload(
                runtime_payload
            )
            qm_tools_log.value = result.stdout or '(no installer output)'
            if result.returncode == 0:
                backend_hint = (
                    f' Reload DELFIN to switch execution backend from '
                    f'{ctx.runtime_backend} to {effective_backend}.'
                    if backend_switch_required
                    else ''
                )
                _set_status(
                    (
                        f'qm_tools installer completed in <code>{html.escape(str(target))}</code>. '
                        f'Runtime now uses <code>{html.escape(str(target))}</code>. '
                        f'Effective ORCA: <code>{html.escape(effective_orca_base or "PATH / auto-detect")}</code>.'
                        f'{backend_hint}'
                    ),
                    color='#2e7d32',
                )
            else:
                _set_status(
                    (
                        f'qm_tools installer failed with exit code {result.returncode}. '
                        f'Inspect the log below. Runtime still points to '
                        f'<code>{html.escape(str(target))}</code>.'
                    ),
                    color='#d32f2f',
                )
        except Exception as exc:
            qm_tools_log.value = ''
            _set_status(
                f'Running qm_tools installer failed: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _on_update_qm_tools(button):
        try:
            target, result = run_qm_tools_installer(
                extra_env={
                    'FORCE_REDOWNLOAD': '1',
                    'FORCE_CONDA_UPDATE': '1',
                }
            )
            runtime_payload = _runtime_payload_from_widgets()
            runtime_payload['qm_tools_root'] = str(target)
            backend_switch_required, effective_backend, effective_orca_base = _persist_runtime_payload(
                runtime_payload
            )
            qm_tools_log.value = result.stdout or '(no updater output)'
            if result.returncode == 0:
                backend_hint = (
                    f' Reload DELFIN to switch execution backend from '
                    f'{ctx.runtime_backend} to {effective_backend}.'
                    if backend_switch_required
                    else ''
                )
                _set_status(
                    (
                        f'qm_tools updated in <code>{html.escape(str(target))}</code>. '
                        'Bundled downloads were refreshed and DELFIN-managed conda envs were updated when used. '
                        'External PATH/.venv tools are still reused, not modified. '
                        f'Effective ORCA: <code>{html.escape(effective_orca_base or "PATH / auto-detect")}</code>.'
                        f'{backend_hint}'
                    ),
                    color='#2e7d32',
                )
            else:
                _set_status(
                    (
                        f'qm_tools update failed with exit code {result.returncode}. '
                        f'Inspect the log below. Runtime still points to '
                        f'<code>{html.escape(str(target))}</code>.'
                    ),
                    color='#d32f2f',
                )
        except Exception as exc:
            qm_tools_log.value = ''
            _set_status(
                f'Updating qm_tools failed: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _make_single_qm_tool_handler(tool_name):
        def _handler(button):
            try:
                extra = {'INSTALL_STD2_FROM_SOURCE': '1'} if tool_name == 'std2' else None
                target, result = run_qm_tools_installer(tools=[tool_name], extra_env=extra)
                runtime_payload = _runtime_payload_from_widgets()
                runtime_payload['qm_tools_root'] = str(target)
                _persist_runtime_payload(runtime_payload)
                qm_tools_log.value = result.stdout or '(no installer output)'
                if result.returncode == 0:
                    _set_status(
                        f'<code>{html.escape(tool_name)}</code> installed in '
                        f'<code>{html.escape(str(target))}</code>.',
                        color='#2e7d32',
                    )
                else:
                    _set_status(
                        f'Installing <code>{html.escape(tool_name)}</code> failed '
                        f'(exit code {result.returncode}). Check log below.',
                        color='#d32f2f',
                    )
            except Exception as exc:
                qm_tools_log.value = ''
                _set_status(
                    f'Installing {html.escape(tool_name)} failed: {html.escape(str(exc))}',
                    color='#d32f2f',
                )
        return _handler

    def _on_install_micromamba(button):
        try:
            import os, shutil, subprocess
            if shutil.which('micromamba') or shutil.which('conda'):
                _set_status('micromamba/conda is already available.', color='#2e7d32')
                return
            local_bin = Path.home() / '.local' / 'bin'
            local_bin.mkdir(parents=True, exist_ok=True)
            result = subprocess.run(
                ['bash', '-c',
                 'URL=https://micro.mamba.pm/api/micromamba/linux-64/latest; '
                 'DST="$HOME/.local/bin/"; '
                 '{ curl -fsSL "$URL" || curl -fsSL --insecure "$URL"; } '
                 '| tar -xvj -C "$DST" --strip-components=1 bin/micromamba'],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
                timeout=120,
            )
            qm_tools_log.value = result.stdout or ''
            mamba_path = local_bin / 'micromamba'
            if result.returncode == 0 and mamba_path.is_file():
                os.environ['PATH'] = f'{local_bin}:{os.environ.get("PATH", "")}'
                _set_status(
                    f'micromamba installed to <code>{html.escape(str(mamba_path))}</code>. '
                    'You can now install xtb, crest, and dftb+.',
                    color='#2e7d32',
                )
            else:
                _set_status(
                    f'micromamba install failed (exit code {result.returncode}). Check log below.',
                    color='#d32f2f',
                )
        except Exception as exc:
            qm_tools_log.value = ''
            _set_status(
                f'Installing micromamba failed: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _on_install_csp_tools(button):
        try:
            target, result = run_csp_tools_installer()
            runtime_payload = _runtime_payload_from_widgets()
            runtime_payload['csp_tools_root'] = str(target)
            backend_switch_required, effective_backend, effective_orca_base = _persist_runtime_payload(
                runtime_payload
            )
            csp_tools_log.value = result.stdout or '(no installer output)'
            if result.returncode == 0:
                _set_status(
                    (
                        f'csp_tools (Genarris) installed in <code>{html.escape(str(target))}</code>. '
                        'Crystal structure prediction tools are now available.'
                    ),
                    color='#2e7d32',
                )
            else:
                _set_status(
                    (
                        f'csp_tools installer failed with exit code {result.returncode}. '
                        'Inspect the CSP tools log below. '
                        'Common issues: missing SWIG, MPI mismatch (OpenMPI vs Intel MPI), '
                        'or missing mpi4py.'
                    ),
                    color='#d32f2f',
                )
        except Exception as exc:
            csp_tools_log.value = ''
            _set_status(
                f'Running csp_tools installer failed: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _on_update_csp_tools(button):
        try:
            target, result = run_csp_tools_installer(
                extra_env={'FORCE_REINSTALL': '1'}
            )
            runtime_payload = _runtime_payload_from_widgets()
            runtime_payload['csp_tools_root'] = str(target)
            backend_switch_required, effective_backend, effective_orca_base = _persist_runtime_payload(
                runtime_payload
            )
            csp_tools_log.value = result.stdout or '(no updater output)'
            if result.returncode == 0:
                _set_status(
                    (
                        f'csp_tools (Genarris) updated in <code>{html.escape(str(target))}</code>. '
                        'Genarris was rebuilt from latest source.'
                    ),
                    color='#2e7d32',
                )
            else:
                _set_status(
                    (
                        f'csp_tools update failed with exit code {result.returncode}. '
                        'Inspect the CSP tools log below.'
                    ),
                    color='#d32f2f',
                )
        except Exception as exc:
            csp_tools_log.value = ''
            _set_status(
                f'Updating csp_tools failed: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _on_install_mlp_tools(button):
        try:
            target, result = run_mlp_tools_installer()
            runtime_payload = _runtime_payload_from_widgets()
            runtime_payload['mlp_tools_root'] = str(target)
            backend_switch_required, effective_backend, effective_orca_base = _persist_runtime_payload(
                runtime_payload
            )
            mlp_tools_log.value = result.stdout or '(no installer output)'
            if result.returncode == 0:
                _set_status(
                    (
                        f'mlp_tools installed in <code>{html.escape(str(target))}</code>. '
                        'ML potential backends (ANI-2x, AIMNet2) are now available.'
                    ),
                    color='#2e7d32',
                )
            else:
                _set_status(
                    (
                        f'mlp_tools installer failed with exit code {result.returncode}. '
                        'Inspect the MLP tools log below.'
                    ),
                    color='#d32f2f',
                )
        except Exception as exc:
            mlp_tools_log.value = ''
            _set_status(
                f'Running mlp_tools installer failed: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _on_update_mlp_tools(button):
        try:
            target, result = run_mlp_tools_installer(
                extra_env={'FORCE_REINSTALL': '1'}
            )
            runtime_payload = _runtime_payload_from_widgets()
            runtime_payload['mlp_tools_root'] = str(target)
            backend_switch_required, effective_backend, effective_orca_base = _persist_runtime_payload(
                runtime_payload
            )
            mlp_tools_log.value = result.stdout or '(no updater output)'
            if result.returncode == 0:
                _set_status(
                    (
                        f'mlp_tools updated in <code>{html.escape(str(target))}</code>. '
                        'All MLP backends were reinstalled.'
                    ),
                    color='#2e7d32',
                )
            else:
                _set_status(
                    (
                        f'mlp_tools update failed with exit code {result.returncode}. '
                        'Inspect the MLP tools log below.'
                    ),
                    color='#d32f2f',
                )
        except Exception as exc:
            mlp_tools_log.value = ''
            _set_status(
                f'Updating mlp_tools failed: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _on_install_analysis_tools(button):
        try:
            target, result = run_analysis_tools_installer()
            analysis_tools_log.value = result.stdout or '(no installer output)'
            if result.returncode == 0:
                _set_status(
                    (
                        f'analysis_tools installed in <code>{html.escape(str(target))}</code>. '
                        'morfeus and CENSO are now available.'
                    ),
                    color='#2e7d32',
                )
            else:
                _set_status(
                    (
                        f'analysis_tools installer failed with exit code {result.returncode}. '
                        'Inspect the log below.'
                    ),
                    color='#d32f2f',
                )
            _refresh_analysis_status()
        except Exception as exc:
            analysis_tools_log.value = ''
            _set_status(
                f'Running analysis_tools installer failed: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _on_update_analysis_tools(button):
        try:
            target, result = run_analysis_tools_installer(
                extra_env={'FORCE_REINSTALL': '1'}
            )
            analysis_tools_log.value = result.stdout or '(no updater output)'
            if result.returncode == 0:
                _set_status(
                    (
                        f'analysis_tools updated in <code>{html.escape(str(target))}</code>.'
                    ),
                    color='#2e7d32',
                )
            else:
                _set_status(
                    (
                        f'analysis_tools update failed with exit code {result.returncode}. '
                        'Inspect the log below.'
                    ),
                    color='#d32f2f',
                )
            _refresh_analysis_status()
        except Exception as exc:
            analysis_tools_log.value = ''
            _set_status(
                f'Updating analysis_tools failed: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _refresh_qm_status(button=None):
        """Refresh the QM tools status box with checkmarks/crosses and install buttons."""
        try:
            runtime_payload = _runtime_payload_from_widgets()
            effective_backend, effective_orca_base, submit_templates_dir = _resolve_runtime_effective_values(
                runtime_payload
            )
            diagnostics = collect_runtime_diagnostics(
                runtime_payload,
                backend=effective_backend,
                effective_orca_base=effective_orca_base,
                submit_templates_dir=submit_templates_dir if effective_backend == 'slurm' else None,
            )

            _QM_TOOL_NAMES = {'xtb', 'crest', 'std2', 'stda', 'xtb4stda', 'dftb+'}
            _QM_DESCRIPTIONS = {
                'xtb': 'Extended tight-binding semi-empirical method (GFN-xTB)',
                'crest': 'Conformer-Rotamer Ensemble Sampling Tool',
                'std2': 'Simplified TD-DFT v2 for UV/Vis spectra',
                'stda': 'Simplified Tamm-Dancoff approximation for UV/Vis',
                'xtb4stda': 'xTB optimized for sTDA input generation',
                'dftb+': 'Density-functional tight-binding method',
            }
            # Maps diagnostic name to install_qm_tools.sh argument
            _QM_INSTALL_ARGS = {
                'xtb': 'xtb',
                'crest': 'crest',
                'dftb+': 'dftb+',
                'stda': 'xtb4stda',
                'xtb4stda': 'xtb4stda',
                'std2': 'std2',
            }

            rows = []
            for item in diagnostics:
                name = item.get('name', '')
                if name not in _QM_TOOL_NAMES:
                    continue
                installed = item.get('status') == 'ok'
                detail = item.get('detail', '')
                desc = _QM_DESCRIPTIONS.get(name, '')
                if installed and detail:
                    desc = f'{desc} — <code>{html.escape(detail)}</code>' if desc else html.escape(detail)

                icon = '&#x2705;' if installed else '&#x274C;'
                color = '#455a64' if installed else '#9e9e9e'
                label = widgets.HTML(
                    f'<span>{icon} <b>{html.escape(name)}</b>'
                    f' &mdash; <span style="color:{color};font-size:0.9em;">'
                    f'{desc}</span></span>',
                    layout=widgets.Layout(min_width='300px'),
                )

                install_arg = _QM_INSTALL_ARGS.get(name)
                if not installed and install_arg:
                    btn = widgets.Button(
                        description='Install',
                        button_style='warning',
                        layout=widgets.Layout(width='80px', height='24px'),
                    )
                    btn.on_click(lambda b, tool=install_arg: _with_buttons_disabled(
                        _make_single_qm_tool_handler(tool)
                    )(b))
                    rows.append(widgets.HBox(
                        [label, btn],
                        layout=widgets.Layout(width='100%', margin='1px 0', align_items='center'),
                    ))
                else:
                    rows.append(widgets.HBox(
                        [label],
                        layout=widgets.Layout(width='100%', margin='1px 0'),
                    ))

            qm_status_box.children = rows if rows else [
                widgets.HTML('<span style="color:#616161;">No QM tool status available.</span>')
            ]
        except Exception as exc:
            logging.exception("Failed to refresh QM tools status")
            qm_status_box.children = [widgets.HTML(
                f'<span style="color:#d32f2f;">Could not load QM tools status: {html.escape(str(exc))}</span>'
            )]

    def _refresh_csp_status(button=None):
        """Refresh the CSP tools status box with checkmarks/crosses."""
        try:
            runtime_payload = _runtime_payload_from_widgets()
            effective_backend, effective_orca_base, submit_templates_dir = _resolve_runtime_effective_values(
                runtime_payload
            )
            diagnostics = collect_runtime_diagnostics(
                runtime_payload,
                backend=effective_backend,
                effective_orca_base=effective_orca_base,
                submit_templates_dir=submit_templates_dir if effective_backend == 'slurm' else None,
            )

            _CSP_NAMES = {'genarris', 'gnrs-cli', 'cgenarris'}
            _CSP_DESCRIPTIONS = {
                'genarris': 'Crystal structure prediction with random sampling',
                'gnrs-cli': 'Genarris command-line interface',
                'cgenarris': 'Genarris C extension for fast structure generation',
            }

            rows = []
            for item in diagnostics:
                name = item.get('name', '')
                if name not in _CSP_NAMES:
                    continue
                installed = item.get('status') == 'ok'
                detail = item.get('detail', '')
                desc = _CSP_DESCRIPTIONS.get(name, '')
                if installed and detail:
                    desc = f'{desc} — <code>{html.escape(detail)}</code>' if desc else html.escape(detail)

                icon = '&#x2705;' if installed else '&#x274C;'
                color = '#455a64' if installed else '#9e9e9e'
                rows.append(widgets.HBox(
                    [widgets.HTML(
                        f'<span>{icon} <b>{html.escape(name)}</b>'
                        f' &mdash; <span style="color:{color};font-size:0.9em;">'
                        f'{desc}</span></span>',
                        layout=widgets.Layout(min_width='300px'),
                    )],
                    layout=widgets.Layout(width='100%', margin='1px 0'),
                ))

            csp_status_box.children = rows if rows else [
                widgets.HTML('<span style="color:#616161;">No CSP tool status available.</span>')
            ]
        except Exception as exc:
            logging.exception("Failed to refresh CSP tools status")
            csp_status_box.children = [widgets.HTML(
                f'<span style="color:#d32f2f;">Could not load CSP tools status: {html.escape(str(exc))}</span>'
            )]

    def _refresh_analysis_status(button=None):
        try:
            from delfin.analysis_tools import collect_analysis_summary
            info = collect_analysis_summary()

            if not _outdated_cache:
                _check_outdated_packages()

            # (pip_install_cmd, pip_pkg_name, conda_cmd)
            _ANALYSIS_PKGS = {
                'Multiwfn': ('', '', ''),                    # manual binary
                'CENSO':    ('censo', 'censo', ''),
                'morfeus':  ('morfeus-ml', 'morfeus-ml', ''),
                'cclib':    ('cclib', 'cclib', ''),
                'nglview':  ('nglview', 'nglview', ''),
                'Packmol':  ('', '', 'packmol'),             # conda only
                'xtb-python': ('', '', 'xtb-python'),        # conda only
            }

            rows = []
            for t in info['tools']:
                cmd, pkg, conda = _ANALYSIS_PKGS.get(
                    t['name'], (t['name'].lower(), t['name'].lower(), ''),
                )
                rows.append(_make_tool_row(
                    t['name'], t['installed'], t.get('version', ''),
                    t['description'], cmd, _refresh_analysis_status,
                    pip_pkg_name=pkg, conda_cmd=conda,
                ))

            analysis_status_box.children = rows
        except Exception as exc:
            logging.exception("Failed to refresh analysis tools status")
            analysis_status_box.children = [widgets.HTML(
                f'<span style="color:#d32f2f;">Could not load analysis tools status: {html.escape(str(exc))}</span>'
            )]

    # ── conda/mamba environment detection ─────────────────────────────────

    def _detect_conda():
        """Detect conda/mamba executable and active environment name.

        Returns (exe, env_name) or (None, None) if not in a conda env.
        """
        import os
        import shutil

        env_name = os.environ.get('CONDA_DEFAULT_ENV', '')
        if not env_name:
            return None, None

        for cmd in ('micromamba', 'mamba', 'conda'):
            path = shutil.which(cmd)
            if path:
                return path, env_name
        return None, None

    # ── pip helpers for install / update ──────────────────────────────────

    _outdated_cache: dict = {}  # package_name_lower -> latest_version

    def _check_outdated_packages():
        """Query pip for outdated packages once; cache the result."""
        import json
        import subprocess
        import sys

        try:
            result = subprocess.run(
                [sys.executable, '-m', 'pip', 'list', '--outdated', '--format=json'],
                capture_output=True, text=True, timeout=60,
            )
            if result.returncode == 0 and result.stdout.strip():
                data = json.loads(result.stdout)
                _outdated_cache.clear()
                for pkg in data:
                    _outdated_cache[pkg['name'].lower()] = pkg['latest_version']
        except Exception:
            pass  # non-critical, just won't show update badges

    def _pip_install_tool(pip_cmd, label, refresh_fn):
        """Run pip install in background and refresh status afterwards."""
        import subprocess
        import sys

        tool_install_log.value = f'Installing {label}...\n'
        try:
            result = subprocess.run(
                [sys.executable, '-m', 'pip', 'install'] + shlex.split(pip_cmd),
                capture_output=True, text=True, timeout=600,
            )
            tool_install_log.value += result.stdout or ''
            if result.stderr:
                tool_install_log.value += result.stderr
            if result.returncode == 0:
                _set_status(f'{label} installed successfully.', color='#2e7d32')
            else:
                _set_status(f'{label} installation failed (exit {result.returncode}).', color='#d32f2f')
        except Exception as exc:
            tool_install_log.value += f'\nError: {exc}'
            _set_status(f'{label} installation error: {exc}', color='#d32f2f')
        # Invalidate cache so next refresh re-checks
        _outdated_cache.clear()
        refresh_fn()

    def _pip_update_tool(pip_pkg, label, refresh_fn):
        """Run pip install --upgrade for a specific package."""
        import subprocess
        import sys

        tool_install_log.value = f'Updating {label}...\n'
        try:
            result = subprocess.run(
                [sys.executable, '-m', 'pip', 'install', '--upgrade'] + shlex.split(pip_pkg),
                capture_output=True, text=True, timeout=600,
            )
            tool_install_log.value += result.stdout or ''
            if result.stderr:
                tool_install_log.value += result.stderr
            if result.returncode == 0:
                _set_status(f'{label} updated successfully.', color='#2e7d32')
            else:
                _set_status(f'{label} update failed (exit {result.returncode}).', color='#d32f2f')
        except Exception as exc:
            tool_install_log.value += f'\nError: {exc}'
            _set_status(f'{label} update error: {exc}', color='#d32f2f')
        _outdated_cache.clear()
        refresh_fn()

    def _conda_install_tool(conda_spec, label, refresh_fn, channel='conda-forge'):
        """Run conda/mamba install in the active environment."""
        import subprocess

        conda_exe, env_name = _detect_conda()
        if not conda_exe:
            tool_install_log.value = (
                f'Cannot install {label}: no conda/mamba found '
                'or not running inside a conda environment.\n'
            )
            _set_status(
                f'{label}: no conda/mamba environment detected.',
                color='#d32f2f',
            )
            return

        tool_install_log.value = (
            f'Installing {label} via {conda_exe} '
            f'into environment "{env_name}"...\n'
        )
        try:
            result = subprocess.run(
                [conda_exe, 'install', '-n', env_name, '-c', channel, '-y']
                + shlex.split(conda_spec),
                capture_output=True, text=True, timeout=600,
            )
            tool_install_log.value += result.stdout or ''
            if result.stderr:
                tool_install_log.value += result.stderr
            if result.returncode == 0:
                _set_status(f'{label} installed successfully.', color='#2e7d32')
            else:
                _set_status(
                    f'{label} installation failed (exit {result.returncode}).',
                    color='#d32f2f',
                )
        except Exception as exc:
            tool_install_log.value += f'\nError: {exc}'
            _set_status(f'{label} installation error: {exc}', color='#d32f2f')
        _outdated_cache.clear()
        refresh_fn()

    def _conda_update_tool(conda_spec, label, refresh_fn, channel='conda-forge'):
        """Run conda/mamba update for a specific package."""
        import subprocess

        conda_exe, env_name = _detect_conda()
        if not conda_exe:
            tool_install_log.value = (
                f'Cannot update {label}: no conda/mamba found '
                'or not running inside a conda environment.\n'
            )
            _set_status(
                f'{label}: no conda/mamba environment detected.',
                color='#d32f2f',
            )
            return

        tool_install_log.value = (
            f'Updating {label} via {conda_exe} '
            f'in environment "{env_name}"...\n'
        )
        try:
            result = subprocess.run(
                [conda_exe, 'update', '-n', env_name, '-c', channel, '-y']
                + shlex.split(conda_spec),
                capture_output=True, text=True, timeout=600,
            )
            tool_install_log.value += result.stdout or ''
            if result.stderr:
                tool_install_log.value += result.stderr
            if result.returncode == 0:
                _set_status(f'{label} updated successfully.', color='#2e7d32')
            else:
                _set_status(
                    f'{label} update failed (exit {result.returncode}).',
                    color='#d32f2f',
                )
        except Exception as exc:
            tool_install_log.value += f'\nError: {exc}'
            _set_status(f'{label} update error: {exc}', color='#d32f2f')
        _outdated_cache.clear()
        refresh_fn()

    def _make_tool_row(name, installed, version, detail_text, install_cmd, refresh_fn,
                       pip_pkg_name='', conda_cmd='', conda_channel='conda-forge'):
        """Create a widget row with Install or Update button.

        Parameters
        ----------
        pip_pkg_name : the pip package name used to check for updates
                       (e.g. 'cclib', 'torchani'). If empty, uses install_cmd.
        conda_cmd : conda/mamba package spec (e.g. 'xtb-python'). When set,
                     the Install/Update buttons use conda/mamba instead of pip.
        conda_channel : conda channel (default 'conda-forge').
        """
        use_conda = bool(conda_cmd)
        _raw = pip_pkg_name or install_cmd or ''
        pkg_key = shlex.split(_raw)[0].lower() if _raw.strip() else ''

        if installed:
            ver = f' v{version}' if version else ''
            # Check if an update is available (pip-based check only)
            latest = _outdated_cache.get(pkg_key, '') if not use_conda else ''

            if latest and version and latest != version:
                # Update available — show update button
                label = widgets.HTML(
                    f'<span>&#x2705; <b>{html.escape(name)}</b>{html.escape(ver)}'
                    f' <span style="color:#ef6c00;font-size:0.85em;">'
                    f'&#x2192; v{html.escape(latest)} available</span>'
                    f' &mdash; <span style="color:#455a64;font-size:0.9em;">'
                    f'{html.escape(detail_text)}</span></span>',
                    layout=widgets.Layout(min_width='300px'),
                )
                if use_conda:
                    btn = widgets.Button(
                        description='Update',
                        button_style='info',
                        layout=widgets.Layout(width='80px', height='24px'),
                        tooltip=f'conda update {conda_cmd}',
                    )
                    btn.on_click(lambda b, spec=conda_cmd, lbl=name, fn=refresh_fn, ch=conda_channel:
                                 _conda_update_tool(spec, lbl, fn, channel=ch))
                else:
                    btn = widgets.Button(
                        description='Update',
                        button_style='info',
                        layout=widgets.Layout(width='80px', height='24px'),
                        tooltip=f'pip install --upgrade {pkg_key}',
                    )
                    btn.on_click(lambda b, pkg=pip_pkg_name or install_cmd, lbl=name, fn=refresh_fn:
                                 _pip_update_tool(pkg, lbl, fn))
                return widgets.HBox(
                    [label, btn],
                    layout=widgets.Layout(width='100%', margin='1px 0', align_items='center'),
                )
            else:
                # Up to date
                label = widgets.HTML(
                    f'<span>&#x2705; <b>{html.escape(name)}</b>{html.escape(ver)}'
                    f' &mdash; <span style="color:#455a64;font-size:0.9em;">'
                    f'{html.escape(detail_text)}</span></span>',
                    layout=widgets.Layout(min_width='300px'),
                )
                return widgets.HBox([label], layout=widgets.Layout(width='100%', margin='1px 0'))
        else:
            # Not installed — show install button
            label = widgets.HTML(
                f'<span>&#x274C; <b>{html.escape(name)}</b>'
                f' &mdash; <span style="color:#9e9e9e;font-size:0.9em;">'
                f'{html.escape(detail_text)}</span></span>',
                layout=widgets.Layout(min_width='300px'),
            )
            has_cmd = conda_cmd if use_conda else install_cmd
            if has_cmd:
                if use_conda:
                    conda_exe_name = (_detect_conda()[0] or 'conda').rsplit('/', 1)[-1]
                    btn = widgets.Button(
                        description='Install',
                        button_style='warning',
                        layout=widgets.Layout(width='80px', height='24px'),
                        tooltip=f'{conda_exe_name} install {conda_cmd}',
                    )
                    btn.on_click(lambda b, spec=conda_cmd, lbl=name, fn=refresh_fn, ch=conda_channel:
                                 _conda_install_tool(spec, lbl, fn, channel=ch))
                else:
                    btn = widgets.Button(
                        description='Install',
                        button_style='warning',
                        layout=widgets.Layout(width='80px', height='24px'),
                        tooltip=f'pip install {install_cmd}',
                    )
                    btn.on_click(lambda b, cmd=install_cmd, lbl=name, fn=refresh_fn:
                                 _pip_install_tool(cmd, lbl, fn))
                return widgets.HBox(
                    [label, btn],
                    layout=widgets.Layout(width='100%', margin='1px 0', align_items='center'),
                )
            else:
                return widgets.HBox(
                    [label],
                    layout=widgets.Layout(width='100%', margin='1px 0'),
                )

    def _refresh_mlp_status(button=None):
        try:
            from delfin.mlp_tools import collect_mlp_summary
            info = collect_mlp_summary()

            if not _outdated_cache:
                _check_outdated_packages()

            # install_cmd -> pip_pkg_name mapping
            _MLP_PKGS = {
                'ANI-2x':     ('torchani',    'torchani'),
                'AIMNet2':    ('aimnet2calc',  'aimnet2calc'),
                'MACE-OFF':   ('mace-torch',  'mace-torch'),
                'CHGNet':     ('chgnet',       'chgnet'),
                'M3GNet':     ('matgl',        'matgl'),
                'SchNetPack': ('schnetpack',   'schnetpack'),
                'NequIP':     ('nequip',       'nequip'),
                'ALIGNN':     ('alignn',       'alignn'),
            }

            def _format_mlp_elements(name, elements):
                if not elements:
                    return 'trainable (any elements)'
                if name in {'CHGNet', 'M3GNet'}:
                    return f'{len(elements)} supported elements (Materials Project coverage)'
                if len(elements) > 20:
                    return f'{len(elements)} supported elements'
                return ', '.join(elements)

            rows = []
            for b in info['backends']:
                elems = _format_mlp_elements(b['name'], b['elements'])
                cmd, pkg = _MLP_PKGS.get(b['name'], (b['name'].lower(), b['name'].lower()))
                rows.append(_make_tool_row(
                    b['name'], b['installed'], b.get('version', ''),
                    elems, cmd, _refresh_mlp_status, pip_pkg_name=pkg,
                ))

            # PyTorch / CUDA info
            torch_line = ''
            if info['torch_version']:
                cuda_badge = (
                    '<span style="color:#2e7d32;font-weight:bold;">CUDA available</span>'
                    if info['cuda']
                    else '<span style="color:#ef6c00;">CPU only</span>'
                )
                torch_line = f'PyTorch {html.escape(info["torch_version"])} &mdash; {cuda_badge}'
            else:
                torch_line = '<span style="color:#d32f2f;">PyTorch not installed</span>'

            gpu_line = ''
            if info['gpu_partition']:
                gpu_line = (
                    f'<br>SLURM GPU partition: '
                    f'<code>{html.escape(info["gpu_partition"])}</code>'
                )

            footer = widgets.HTML(
                f'<div style="margin-top:4px;font-size:0.9em;">{torch_line}{gpu_line}</div>'
            )

            mlp_status_box.children = rows + [footer]
        except Exception as exc:
            logging.exception("Failed to refresh MLP status")
            mlp_status_box.children = [widgets.HTML(
                f'<span style="color:#d32f2f;">Could not load MLP status: {html.escape(str(exc))}</span>'
            )]

    def _refresh_ai_status(button=None):
        try:
            from delfin.ai_tools import collect_ai_summary
            info = collect_ai_summary()

            if not _outdated_cache:
                _check_outdated_packages()

            children = []
            # Group by category
            categories = {}
            for t in info['tools']:
                categories.setdefault(t['category'], []).append(t)

            for cat, tools in categories.items():
                children.append(widgets.HTML(
                    f'<div style="margin-top:6px;font-weight:bold;color:#37474f;">{html.escape(cat)}</div>'
                ))
                for t in tools:
                    hint = t.get('install_hint', '')
                    # Extract pip package spec from install_hint
                    if hint.startswith('pip install '):
                        cmd = hint[len('pip install '):]
                    elif hint:
                        cmd = hint
                    else:
                        cmd = ''
                    pkg = shlex.split(cmd)[0] if cmd else ''
                    children.append(_make_tool_row(
                        t['name'], t['installed'], t.get('version', ''),
                        t['description'] if t['installed'] else hint,
                        cmd, _refresh_ai_status, pip_pkg_name=pkg,
                    ))

            n_installed = sum(1 for t in info['tools'] if t['installed'])
            n_total = len(info['tools'])
            children.append(widgets.HTML(
                f'<div style="margin-top:6px;color:#546e7a;">'
                f'{n_installed}/{n_total} AI tools installed</div>'
            ))

            ai_status_box.children = children
        except Exception as exc:
            logging.exception("Failed to refresh AI tools status")
            ai_status_box.children = [widgets.HTML(
                f'<span style="color:#d32f2f;">Could not load AI tools status: {html.escape(str(exc))}</span>'
            )]

    def _on_pip_install_editable(button):
        import sys

        try:
            result = run_pip_install_editable()
            pip_install_log.value = result.stdout or '(no output)'
            if result.returncode == 0:
                _set_status(
                    (
                        f'<code>pip install -e .</code> completed successfully '
                        f'using <code>{html.escape(sys.executable)}</code>. '
                        'New CLI commands are now available.'
                    ),
                    color='#2e7d32',
                )
            else:
                _set_status(
                    (
                        f'<code>pip install -e .</code> failed with exit code {result.returncode}. '
                        'Inspect the pip install log below.'
                    ),
                    color='#d32f2f',
                )
        except Exception as exc:
            pip_install_log.value = ''
            _set_status(
                f'pip install -e . failed: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _on_rebuild_venv(button):
        try:
            voila_port_str = os.environ.get('DELFIN_VOILA_PORT', '')
            voila_port = int(voila_port_str) if voila_port_str.isdigit() else None
            log_file = rebuild_venv(
                repo_dir=getattr(ctx, 'repo_dir', None),
                voila_port=voila_port,
            )
            if voila_port:
                restart_msg = (
                    f'The dashboard will restart automatically on port {voila_port}. '
                    '<b>Reload this page</b> after ~2 minutes.'
                )
            else:
                restart_msg = (
                    'Could not detect the Voilà port — run '
                    '<code>delfin-voila</code> manually after the rebuild.'
                )
            pip_install_log.value = (
                f'Background rebuild started.\n'
                f'Log: {log_file}\n\n'
                'The dashboard will stop responding shortly because the\n'
                'old .venv is being deleted.  This is expected.\n\n'
                + (f'Voilà will restart on port {voila_port} when done.\n'
                   'Just reload this browser tab.\n'
                   if voila_port else
                   'Run delfin-voila again after ~2 minutes.\n')
            )
            _set_status(restart_msg, color='#ef6c00')
        except Exception as exc:
            pip_install_log.value = ''
            _set_status(
                f'venv rebuild failed to launch: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _on_setup_bwunicluster(button):
        try:
            calc_override, archive_override, effective_calc_dir, effective_archive_dir = _effective_paths_from_widgets()
            runtime_payload = _runtime_payload_from_widgets()
            prepared = prepare_bwunicluster_user_setup(
                repo_dir=getattr(ctx, 'repo_dir', None),
                calc_dir=effective_calc_dir,
                archive_dir=effective_archive_dir,
                orca_base=runtime_payload.get('slurm', {}).get('orca_base')
                or runtime_payload.get('orca_base', ''),
                qm_tools_root=runtime_payload.get('qm_tools_root', ''),
                install_qm_tools=True,
            )
            prepared_runtime = prepared.get('runtime', {}) or {}
            prepared_runtime['tool_binaries'] = dict(runtime_payload.get('tool_binaries', {}) or {})
            prepared['runtime'] = prepared_runtime

            settings_payload = load_settings()
            paths_payload = {}
            if calc_override or prepared["paths"].get("calculations_dir"):
                paths_payload["calculations_dir"] = str(prepared["paths"]["calculations_dir"])
            if archive_override or prepared["paths"].get("archive_dir"):
                paths_payload["archive_dir"] = str(prepared["paths"]["archive_dir"])
            settings_payload["paths"] = paths_payload
            settings_payload["runtime"] = prepared["runtime"]
            settings_payload = save_settings(settings_payload, settings_path)

            _apply_workspace_paths(Path(prepared["paths"]["calculations_dir"]), Path(prepared["paths"]["archive_dir"]))
            _set_runtime_widgets(settings_payload)
            backend_switch_required, effective_backend, effective_orca_base = _apply_runtime_settings(
                settings_payload.get('runtime', {}) or prepared["runtime"]
            )
            _render_runtime_diagnostics(
                settings_payload.get('runtime', {}) or prepared["runtime"],
                reload_required=backend_switch_required,
            )

            log_lines = [
                f"Prepared bwUniCluster runtime profile.",
                f"Submit templates: {prepared['submit_templates_dir']}",
                f"qm_tools root: {prepared['qm_tools_root']}",
                f"ORCA: {prepared['orca_base'] or 'not detected'}",
                f"Environment file: {prepared['env_file']}",
            ]
            if prepared.get("venv_tarball"):
                log_lines.append(f"Venv tarball: {prepared['venv_tarball']}")
            if prepared.get("runtime_cache_dir"):
                log_lines.append(f"Runtime cache: {prepared['runtime_cache_dir']}")
            if prepared.get("shell_files"):
                log_lines.append("Shell rc updated: " + ", ".join(prepared["shell_files"]))
            installer_output = str(prepared.get("qm_tools_installer_output") or "").strip()
            if installer_output:
                log_lines.append("")
                log_lines.append(installer_output)
            qm_tools_log.value = "\n".join(log_lines)

            backend_hint = (
                f' Reload DELFIN to switch execution backend from '
                f'{ctx.runtime_backend} to {effective_backend}.'
                if backend_switch_required
                else ''
            )
            _set_status(
                (
                    'bwUniCluster setup completed. '
                    f'Runtime now uses <code>slurm</code> with profile <code>bwunicluster3</code>, '
                    f'ORCA <code>{html.escape(effective_orca_base or "not detected")}</code>, '
                    f'qm_tools <code>{html.escape(str(prepared["qm_tools_root"]))}</code>, and '
                    f'submit templates <code>{html.escape(str(prepared["submit_templates_dir"]))}</code>. '
                    f'New shells will source <code>{html.escape(str(prepared["env_file"]))}</code>.'
                    f'{backend_hint}'
                ),
                color='#2e7d32',
            )
        except Exception as exc:
            _set_status(
                f'bwUniCluster setup failed: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _on_verify_bwunicluster(button):
        try:
            _calc_override, _archive_override, effective_calc_dir, effective_archive_dir = _effective_paths_from_widgets()
            runtime_payload = _runtime_payload_from_widgets()
            checks = collect_bwunicluster_verification(
                repo_dir=getattr(ctx, 'repo_dir', None),
                orca_base=runtime_payload.get('slurm', {}).get('orca_base')
                or runtime_payload.get('orca_base', ''),
                calc_dir=effective_calc_dir,
                archive_dir=effective_archive_dir,
            )
            rows = []
            missing = 0
            for item in checks:
                status = item.get("status", "missing")
                if status not in ("ok", "module"):
                    missing += 1
                if status == "ok":
                    status_color, status_label = "#2e7d32", "OK"
                elif status == "module":
                    status_color, status_label = "#e65100", "Module"
                else:
                    status_color, status_label = "#d32f2f", "Missing"
                rows.append(
                    "<tr>"
                    f"<td style='padding:4px 8px;'><code>{html.escape(str(item.get('name', '')))}</code></td>"
                    f"<td style='padding:4px 8px; color:{status_color}; font-weight:600;'>{status_label}</td>"
                    f"<td style='padding:4px 8px;'><code>{html.escape(str(item.get('detail', '')))}</code></td>"
                    "</tr>"
                )
            qm_tools_log.value = (
                "=== bwUniCluster verification ===\n"
                + "\n".join(
                    f"{item['name']}: {item['status']} - {item['detail']}"
                    for item in checks
                )
            )
            if missing:
                _set_status(
                    f'Verify bwUniCluster found {missing} missing item(s). Review the log below. This check does not modify your system.',
                    color='#ef6c00',
                )
            else:
                _set_status(
                    'Verify bwUniCluster passed. No missing items were detected. This check did not modify your system.',
                    color='#2e7d32',
                )
        except Exception as exc:
            _set_status(
                f'Verify bwUniCluster failed: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _on_full_install_bwunicluster(button):
        try:
            _set_status(
                'Running the full bwUniCluster installer script. This can take a while because OpenMPI may be built and the repo venv may be recreated.',
                color='#ef6c00',
            )
            _on_detect_local_resources(button)
            _refresh_orca_dropdown(global_orca_input.value)
            _calc_override, _archive_override, effective_calc_dir, effective_archive_dir = _effective_paths_from_widgets()
            runtime_payload = _runtime_payload_from_widgets()
            result = run_bwunicluster_installer(
                repo_dir=getattr(ctx, 'repo_dir', None),
                orca_base=runtime_payload.get('slurm', {}).get('orca_base')
                or runtime_payload.get('orca_base', ''),
                calc_dir=effective_calc_dir,
                archive_dir=effective_archive_dir,
            )
            qm_tools_log.value = result.stdout or '(no installer output)'
            if result.returncode != 0:
                _set_status(
                    (
                        f'Full bwUniCluster install failed with exit code {result.returncode}. '
                        'Inspect the log below. The installer script is '
                        f'<code>{html.escape(str(Path(ctx.repo_dir or ".") / "scripts" / "install_delfin_bwu.sh"))}</code>.'
                    ),
                    color='#d32f2f',
                )
                return

            settings_payload = _load_settings_to_widgets(set_status=False) or {}
            runtime_payload = settings_payload.get('runtime', {}) or {}
            paths_payload = settings_payload.get('paths', {}) or {}
            effective_calc_dir = Path(paths_payload.get('calculations_dir') or ctx.default_calc_dir)
            effective_archive_dir = Path(paths_payload.get('archive_dir') or ctx.default_archive_dir)
            _apply_workspace_paths(effective_calc_dir, effective_archive_dir)
            backend_switch_required, effective_backend, effective_orca_base = _apply_runtime_settings(
                runtime_payload
            )
            _render_runtime_diagnostics(
                runtime_payload,
                reload_required=backend_switch_required,
            )
            backend_hint = (
                f' Reload DELFIN to switch execution backend from {ctx.runtime_backend} to {effective_backend}.'
                if backend_switch_required
                else ''
            )
            _set_status(
                (
                    'Full bwUniCluster install completed via '
                    f'<code>{html.escape(str(Path(ctx.repo_dir or ".") / "scripts" / "install_delfin_bwu.sh"))}</code>. '
                    f'Runtime resolves to <code>{html.escape(effective_backend)}</code> with ORCA '
                    f'<code>{html.escape(effective_orca_base or "PATH / auto-detect")}</code>.'
                    f'{backend_hint}'
                ),
                color='#2e7d32',
            )
        except Exception as exc:
            _set_status(
                f'Full bwUniCluster install failed: {html.escape(str(exc))}',
                color='#d32f2f',
            )

    def _on_save(button):
        try:
            host, user, remote_path, port = _transfer_payload_from_widgets()
            host = str(host or '').strip()
            user = str(user or '').strip()
            remote_path = str(remote_path or '').strip()
            calc_override, archive_override, effective_calc_dir, effective_archive_dir = _effective_paths_from_widgets()
            runtime_payload = _runtime_payload_from_widgets()
            settings_payload = load_settings()
            if host or user or remote_path:
                host, user, remote_path, port = normalize_ssh_transfer_settings(
                    host,
                    user,
                    remote_path,
                    port,
                )
                settings_payload['transfer'] = {
                    'host': host,
                    'user': user,
                    'remote_path': remote_path,
                    'port': port,
                }
            else:
                settings_payload['transfer'] = {}
            paths_payload = {}
            if calc_override:
                paths_payload['calculations_dir'] = calc_override
            if archive_override:
                paths_payload['archive_dir'] = archive_override
            settings_payload['paths'] = paths_payload
            settings_payload['runtime'] = runtime_payload
            settings_payload.setdefault('features', {})
            settings_payload['features']['remote_archive_enabled'] = bool(remote_archive_toggle.value)
            settings_payload.setdefault('ui', {})
            settings_payload['ui']['tabs'] = _normalized_tab_prefs()
            settings_payload = save_settings(settings_payload, settings_path)
            _apply_workspace_paths(effective_calc_dir, effective_archive_dir)
            backend_switch_required, effective_backend, effective_orca_base = _apply_runtime_settings(
                settings_payload.get('runtime', {}) or runtime_payload
            )
        except Exception as exc:
            _set_status(
                (
                    f'Saving settings failed: {html.escape(str(exc))}. '
                    'No passwords are stored; use SSH keys or ssh-agent for transfers.'
                ),
                color='#d32f2f',
            )
            return

        toggle_changed = bool(remote_archive_toggle.value) != bool(state.get('remote_archive_enabled'))
        _set_transfer_widgets(settings_payload.get('transfer', {}) or {})
        _set_path_widgets(settings_payload)
        _set_runtime_widgets(settings_payload)
        _set_remote_archive_widget(settings_payload)
        _set_tab_preferences(settings_payload)
        _render_runtime_diagnostics(
            settings_payload.get('runtime', {}) or {},
            reload_required=backend_switch_required,
        )
        reload_hint = ' Reload DELFIN to apply Remote Archive tab/button visibility.' if toggle_changed else ''
        tabs_hint = ''
        if callable(getattr(ctx, 'rebuild_dashboard_tabs', None)):
            current_title = None
            try:
                selected_idx = ctx.tabs_widget.selected_index
                current_title = next(
                    title for title, index in (ctx.tab_indices or {}).items()
                    if index == selected_idx
                )
            except Exception:
                current_title = None
            try:
                ctx.rebuild_dashboard_tabs(selected_title=current_title or 'Settings')
                tabs_hint = ' Tab visibility/order applied to the current session.'
            except Exception:
                tabs_hint = ' Reload DELFIN to apply tab visibility/order changes.'
        backend_hint = (
            f' Reload DELFIN to switch execution backend from '
            f'<code>{html.escape(ctx.runtime_backend)}</code> to '
            f'<code>{html.escape(effective_backend)}</code>.'
            if backend_switch_required
            else ' Runtime settings applied to the current session.'
        )
        _set_status(
            (
                f'Saved local settings to '
                f'<code>{html.escape(str(settings_path))}</code>. '
                f'Calculations now point to <code>{html.escape(str(effective_calc_dir))}</code>; '
                f'Archive now points to <code>{html.escape(str(effective_archive_dir))}</code>. '
                f'Runtime resolves to <code>{html.escape(effective_backend)}</code> with ORCA '
                f'<code>{html.escape(effective_orca_base or "PATH / auto-detect")}</code>.'
                f'{backend_hint}{reload_hint}{tabs_hint}'
            ),
            color='#2e7d32',
        )
        save_btn.disabled = True

    # ── dirty-tracking: enable Save button when any setting widget changes ──
    _settings_widgets_to_watch = [
        calc_path_input, archive_path_input,
        backend_dropdown, global_orca_input, qm_tools_root_input,
        csp_tools_root_input, mlp_tools_root_input,
        *[tool_binary_inputs[name] for name in _selectable_tool_names],
        local_orca_input, local_max_cores_input, local_max_ram_input,
        slurm_orca_input, slurm_templates_input, slurm_profile_input,
        host_hidden, host_visible, user_hidden, user_visible,
        path_hidden, path_visible, port_input,
        remote_archive_toggle,
    ]
    for _w in _settings_widgets_to_watch:
        _w.observe(_mark_dirty, names='value')

    show_sensitive_btn.observe(_on_toggle_sensitive, names='value')
    reload_btn.on_click(_on_reload)
    validate_runtime_btn.on_click(_on_validate_runtime)
    scan_orca_btn.on_click(_on_scan_orca)
    detect_local_resources_btn.on_click(_on_detect_local_resources)
    prepare_qm_tools_btn.on_click(_on_prepare_qm_tools)
    install_qm_tools_btn.on_click(_with_buttons_disabled(_on_install_qm_tools))
    update_qm_tools_btn.on_click(_with_buttons_disabled(_on_update_qm_tools))
    install_xtb_btn.on_click(_with_buttons_disabled(_make_single_qm_tool_handler('xtb')))
    install_crest_btn.on_click(_with_buttons_disabled(_make_single_qm_tool_handler('crest')))
    install_dftbplus_btn.on_click(_with_buttons_disabled(_make_single_qm_tool_handler('dftb+')))
    install_stda_btn.on_click(_with_buttons_disabled(_make_single_qm_tool_handler('xtb4stda')))
    install_std2_btn.on_click(_with_buttons_disabled(_make_single_qm_tool_handler('std2')))
    install_micromamba_btn.on_click(_with_buttons_disabled(_on_install_micromamba))
    install_csp_tools_btn.on_click(_with_buttons_disabled(_on_install_csp_tools))
    update_csp_tools_btn.on_click(_with_buttons_disabled(_on_update_csp_tools))
    install_mlp_tools_btn.on_click(_with_buttons_disabled(_on_install_mlp_tools))
    update_mlp_tools_btn.on_click(_with_buttons_disabled(_on_update_mlp_tools))
    install_analysis_tools_btn.on_click(_with_buttons_disabled(_on_install_analysis_tools))
    update_analysis_tools_btn.on_click(_with_buttons_disabled(_on_update_analysis_tools))
    refresh_qm_status_btn.on_click(_refresh_qm_status)
    refresh_csp_status_btn.on_click(_refresh_csp_status)
    refresh_analysis_status_btn.on_click(_refresh_analysis_status)
    pip_install_btn.on_click(_with_buttons_disabled(_on_pip_install_editable))
    rebuild_venv_btn.on_click(_with_buttons_disabled(_on_rebuild_venv))
    refresh_mlp_status_btn.on_click(_refresh_mlp_status)
    refresh_ai_status_btn.on_click(_refresh_ai_status)
    setup_bwunicluster_btn.on_click(_with_buttons_disabled(_on_setup_bwunicluster))
    verify_bwunicluster_btn.on_click(_with_buttons_disabled(_on_verify_bwunicluster))
    full_install_bwunicluster_btn.on_click(_with_buttons_disabled(_on_full_install_bwunicluster))
    save_btn.on_click(_on_save)
    detected_orca_dropdown.observe(_on_select_detected_orca, names='value')
    global_orca_input.observe(_on_change_global_orca, names='value')
    for _tool_name in _selectable_tool_names:
        tool_detected_dropdowns[_tool_name].observe(
            lambda change, tool_name=_tool_name: (
                None if change.get('name') != 'value' else tool_binary_inputs[tool_name].set_trait('value', str(change.get('new') or ''))
            ),
            names='value',
        )
        tool_binary_inputs[_tool_name].observe(
            lambda change, tool_name=_tool_name: (
                None if change.get('name') != 'value' else _refresh_tool_dropdown(tool_name, change.get('new'))
            ),
            names='value',
        )
        tool_scan_buttons[_tool_name].on_click(
            lambda _button, tool_name=_tool_name: _refresh_tool_dropdown(tool_name)
        )

    info_box = widgets.HTML(
        value=(
            f'<div style="color:#455a64;">'
            f'<b>Local settings:</b> <code>{html.escape(str(settings_path))}</code><br>'
            'This file lives in your HOME, outside the git repo. '
            '<code>git pull</code> will not overwrite it. Missing future settings are '
            'added without replacing your existing values.<br>'
            f'Leave path fields empty to use the defaults '
            f'<code>{html.escape(str(ctx.default_calc_dir))}</code> and '
            f'<code>{html.escape(str(ctx.default_archive_dir))}</code>.<br>'
            'Leave runtime path fields empty to prefer environment variables, '
            '<code>PATH</code>, or DELFIN auto-detection. '
            'Passwords are not stored here. Resumable SSH transfers use SSH keys or '
            '<code>ssh-agent</code>.'
            '</div>'
        )
    )
    # ── helper: standard row layout ─────────────────────────────────────
    _row_layout = widgets.Layout(width='100%', gap='8px', align_items='center', flex_flow='row wrap')
    _section_layout = widgets.Layout(width='100%', gap='10px')

    # ══════════════════════════════════════════════════════════════════════
    # Top bar: always-visible buttons + status
    # ══════════════════════════════════════════════════════════════════════
    top_bar = widgets.VBox(
        [
            widgets.HBox(
                [save_btn, reload_btn, show_sensitive_btn, validate_runtime_btn],
                layout=widgets.Layout(gap='8px', align_items='center', flex_flow='row wrap'),
            ),
            status_html,
        ],
        layout=_section_layout,
    )

    # ══════════════════════════════════════════════════════════════════════
    # Section 1: Workspace Paths
    # ══════════════════════════════════════════════════════════════════════
    workspace_section = widgets.VBox(
        [
            widgets.HBox(
                [widgets.HTML('<b>Calculations</b>'), calc_path_input],
                layout=_row_layout,
            ),
            widgets.HBox(
                [widgets.HTML('<b>Archive</b>'), archive_path_input],
                layout=_row_layout,
            ),
        ],
        layout=_section_layout,
    )

    tabs_section = widgets.VBox(
        [
            widgets.HTML(
                '<div style="color:#455a64;">'
                'Choose which dashboard tabs are visible and in which order they appear. '
                '<code>Settings</code> is always visible and stays at the end.'
                '</div>'
            ),
            tabs_status_html,
            tabs_rows_box,
        ],
        layout=_section_layout,
    )

    # ══════════════════════════════════════════════════════════════════════
    # Section 2: Runtime Backend
    # ══════════════════════════════════════════════════════════════════════
    tool_selector_rows = []
    if _selectable_tool_names:
        tool_selector_rows.extend([
            widgets.HTML(
                '<hr style="border:none; border-top:1px solid #e0e0e0; margin:4px 0;">'
                '<b style="color:#455a64;">Preferred Program Binaries</b>'
            )
        ])
        for tool_name in _selectable_tool_names:
            display_name = _tool_display_names.get(tool_name, tool_name.upper())
            tool_selector_rows.append(
                widgets.HBox(
                    [widgets.HTML(f'<b>{display_name} binary</b>'), tool_binary_inputs[tool_name], tool_scan_buttons[tool_name]],
                    layout=_row_layout,
                )
            )
            tool_selector_rows.append(
                widgets.HBox(
                    [
                        widgets.HTML(f'<b>Detected {display_name}</b>'),
                        tool_detected_dropdowns[tool_name],
                        widgets.HTML(
                            '<span style="color:#616161;">Choose a detected installation to fill the preferred binary automatically.</span>'
                        ),
                    ],
                    layout=_row_layout,
                )
            )

    runtime_section = widgets.VBox(
        [
            widgets.HBox(
                [
                    widgets.HTML('<b>Backend</b>'),
                    backend_dropdown,
                    widgets.HTML(
                        '<span style="color:#616161;">Auto chooses SLURM when <code>sbatch</code> is available.</span>'
                    ),
                ],
                layout=_row_layout,
            ),
            widgets.HBox(
                [widgets.HTML('<b>ORCA path</b>'), global_orca_input, scan_orca_btn],
                layout=_row_layout,
            ),
            widgets.HBox(
                [
                    widgets.HTML('<b>Detected ORCA</b>'),
                    detected_orca_dropdown,
                    widgets.HTML(
                        '<span style="color:#616161;">Choose a detected installation to fill the ORCA path automatically.</span>'
                    ),
                ],
                layout=_row_layout,
            ),
            *tool_selector_rows,
            widgets.HTML(
                '<hr style="border:none; border-top:1px solid #e0e0e0; margin:4px 0;">'
                '<b style="color:#455a64;">Local overrides</b>'
            ),
            widgets.HBox(
                [widgets.HTML('<b>Local ORCA</b>'), local_orca_input],
                layout=_row_layout,
            ),
            widgets.HBox(
                [
                    widgets.HTML('<b>Local max cores</b>'),
                    local_max_cores_input,
                    widgets.HTML('<b>Local max RAM (MB)</b>'),
                    local_max_ram_input,
                    detect_local_resources_btn,
                ],
                layout=_row_layout,
            ),
            widgets.HTML(
                '<hr style="border:none; border-top:1px solid #e0e0e0; margin:4px 0;">'
                '<b style="color:#455a64;">SLURM / Cluster</b>'
            ),
            widgets.HBox(
                [
                    setup_bwunicluster_btn,
                    verify_bwunicluster_btn,
                    full_install_bwunicluster_btn,
                    widgets.HTML(
                        '<span style="color:#616161;">'
                        '<b>Setup</b> prepares an existing DELFIN install for bwUniCluster. '
                        '<b>Verify</b> is read-only. '
                        '<b>Full install</b> runs the packaged installer.'
                        '</span>'
                    ),
                ],
                layout=_row_layout,
            ),
            widgets.HBox(
                [widgets.HTML('<b>SLURM ORCA</b>'), slurm_orca_input],
                layout=_row_layout,
            ),
            widgets.HBox(
                [widgets.HTML('<b>Submit templates</b>'), slurm_templates_input],
                layout=_row_layout,
            ),
            widgets.HBox(
                [widgets.HTML('<b>Site profile</b>'), slurm_profile_input],
                layout=_row_layout,
            ),
            runtime_diagnostics_html,
        ],
        layout=_section_layout,
    )

    # ══════════════════════════════════════════════════════════════════════
    # Section 3: Tool Installation (sub-accordion)
    # ══════════════════════════════════════════════════════════════════════

    # -- QM Tools --
    qm_tools_section = widgets.VBox(
        [
            widgets.HBox(
                [widgets.HTML('<b>qm_tools root</b>'), qm_tools_root_input],
                layout=_row_layout,
            ),
            widgets.HBox(
                [
                    prepare_qm_tools_btn,
                    install_qm_tools_btn,
                    update_qm_tools_btn,
                    refresh_qm_status_btn,
                    widgets.HTML(
                        f'<span style="color:#616161;">'
                        f'Installs into <code>{html.escape(str(get_user_qm_tools_dir()))}</code>. '
                        'ORCA stays external; PATH/.venv tools are reused.'
                        f'</span>'
                    ),
                ],
                layout=_row_layout,
            ),
            qm_status_box,
            qm_tools_log,
        ],
        layout=_section_layout,
    )

    # -- CSP Tools --
    csp_tools_section = widgets.VBox(
        [
            widgets.HBox(
                [widgets.HTML('<b>csp_tools root</b>'), csp_tools_root_input],
                layout=_row_layout,
            ),
            widgets.HBox(
                [
                    install_csp_tools_btn,
                    update_csp_tools_btn,
                    refresh_csp_status_btn,
                    widgets.HTML(
                        f'<span style="color:#616161;">'
                        f'Installs Genarris into <code>{html.escape(str(get_user_csp_tools_dir()))}</code>. '
                        'Requires SWIG, MPI (OpenMPI), and mpi4py.'
                        f'</span>'
                    ),
                ],
                layout=_row_layout,
            ),
            csp_status_box,
            csp_tools_log,
        ],
        layout=_section_layout,
    )

    # -- MLP Tools --
    mlp_tools_section = widgets.VBox(
        [
            widgets.HBox(
                [widgets.HTML('<b>mlp_tools root</b>'), mlp_tools_root_input],
                layout=_row_layout,
            ),
            widgets.HBox(
                [
                    install_mlp_tools_btn,
                    update_mlp_tools_btn,
                    refresh_mlp_status_btn,
                    widgets.HTML(
                        f'<span style="color:#616161;">'
                        f'Installs ANI-2x + AIMNet2 into <code>{html.escape(str(get_user_mlp_tools_dir()))}</code>. '
                        'Requires PyTorch.'
                        f'</span>'
                    ),
                ],
                layout=_row_layout,
            ),
            mlp_status_box,
            mlp_tools_log,
        ],
        layout=_section_layout,
    )

    # -- Analysis Tools --
    analysis_tools_section = widgets.VBox(
        [
            widgets.HBox(
                [
                    install_analysis_tools_btn,
                    update_analysis_tools_btn,
                    refresh_analysis_status_btn,
                    widgets.HTML(
                        f'<span style="color:#616161;">'
                        f'Installs morfeus + CENSO into <code>{html.escape(str(get_user_analysis_tools_dir()))}</code>. '
                        'Multiwfn requires manual binary download.'
                        f'</span>'
                    ),
                ],
                layout=_row_layout,
            ),
            analysis_status_box,
            analysis_tools_log,
        ],
        layout=_section_layout,
    )

    # -- AI/ML Tools --
    ai_tools_section = widgets.VBox(
        [
            widgets.HBox(
                [
                    refresh_ai_status_btn,
                    widgets.HTML(
                        '<span style="color:#616161;">'
                        'Optional AI tools — install individually via pip as needed. '
                        'See install hints in the status panel.'
                        '</span>'
                    ),
                ],
                layout=_row_layout,
            ),
            ai_status_box,
            tool_install_log,
        ],
        layout=_section_layout,
    )

    # Sub-accordion for tool categories
    tools_accordion = widgets.Accordion(
        children=[
            qm_tools_section,
            csp_tools_section,
            mlp_tools_section,
            analysis_tools_section,
            ai_tools_section,
        ],
    )
    tools_accordion.set_title(0, 'QM Tools (xtb, crest, dftb+, stda, std2)')
    tools_accordion.set_title(1, 'CSP Tools (Crystal Structure Prediction)')
    tools_accordion.set_title(2, 'MLP Tools (Machine Learning Potentials)')
    tools_accordion.set_title(3, 'Analysis Tools (Multiwfn, CENSO, morfeus, cclib, nglview, Packmol)')
    tools_accordion.set_title(4, 'AI/ML Tools (Foundation Models, Generative, Retrosynthesis, ADMET)')
    tools_accordion.selected_index = None  # all collapsed by default

    tools_section = widgets.VBox(
        [tools_accordion],
        layout=_section_layout,
    )

    # ══════════════════════════════════════════════════════════════════════
    # Section 4: Developer
    # ══════════════════════════════════════════════════════════════════════
    developer_section = widgets.VBox(
        [
            widgets.HBox(
                [
                    rebuild_venv_btn,
                    widgets.HTML(
                        '<span style="color:#616161;">'
                        'Deletes <code>.venv</code>, recreates it, runs '
                        '<code>pip install -e .</code>, and packages '
                        '<code>delfin_venv.tar</code> for job staging. '
                        'Use after upgrading Python or when the venv is broken.'
                        '</span>'
                    ),
                ],
                layout=_row_layout,
            ),
            widgets.HBox(
                [
                    pip_install_btn,
                    widgets.HTML(
                        '<span style="color:#616161;">'
                        'Runs <code>pip install -e .</code> in the DELFIN repo root. '
                        'Registers new CLI entry-points without restarting.'
                        '</span>'
                    ),
                ],
                layout=_row_layout,
            ),
            pip_install_log,
            widgets.HTML(
                '<hr style="border:none; border-top:1px solid #e0e0e0; margin:4px 0;">'
                '<b style="color:#455a64;">Prerequisites</b>'
            ),
            widgets.HBox(
                [
                    install_micromamba_btn,
                    widgets.HTML(
                        '<span style="color:#616161;">'
                        'Required for conda-based QM tools (xtb, crest, dftb+).'
                        '</span>'
                    ),
                ],
                layout=_row_layout,
            ),
        ],
        layout=_section_layout,
    )

    # ══════════════════════════════════════════════════════════════════════
    # Section 5: SSH Transfer & Remote Archive
    # ══════════════════════════════════════════════════════════════════════
    ssh_help_box = widgets.Accordion(
        children=[
            widgets.HTML(
                value=(
                    '<div style="color:#455a64; padding:6px;">'
                    '<b>Option A &ndash; SSH key (simple servers):</b><br>'
                    'If the remote server accepts key-based login without 2FA:<br>'
                    '<pre style="margin:8px 0 0 0; white-space:pre-wrap;">'
                    'ssh-keygen -t ed25519\n'
                    'ssh-copy-id USER@HOST\n'
                    'ssh USER@HOST'
                    '</pre>'
                    'When running <code>ssh-copy-id</code>, enter the target account password once. '
                    'If the final <code>ssh</code> test succeeds, DELFIN will work afterward.<br><br>'
                    '<b>Option B &ndash; SSH ControlMaster (servers with 2FA/TOTP):</b><br>'
                    'Some HPC systems require a one-time password (TOTP) at every login. '
                    'SSH ControlMaster keeps one authenticated session open so that all '
                    'subsequent transfers run without re-entering the code.<br><br>'
                    '<b>Step 1 &ndash; Create an SSH key (skip if you already have one):</b><br>'
                    '<pre style="margin:4px 0; white-space:pre-wrap;">'
                    'ssh-keygen -t ed25519'
                    '</pre>'
                    'You will be asked three questions &ndash; just press <b>Enter</b> each time:<br>'
                    '&nbsp;&nbsp;1. <code>Enter file in which to save the key (...)</code> &rarr; Enter<br>'
                    '&nbsp;&nbsp;2. <code>Enter passphrase (empty for no passphrase)</code> &rarr; Enter<br>'
                    '&nbsp;&nbsp;3. <code>Enter same passphrase again</code> &rarr; Enter<br><br>'
                    '<b>Step 2 &ndash; Create or edit ~/.ssh/config:</b><br>'
                    'Open the file in a terminal editor and add the block below. '
                    'Replace <code>REMOTE_ALIAS</code>, <code>REMOTE_HOST</code>, '
                    'and <code>YOUR_USERNAME</code> with your values.<br>'
                    '<pre style="margin:4px 0; white-space:pre-wrap;">'
                    'nano ~/.ssh/config'
                    '</pre>'
                    'Paste this block at the end of the file:<br>'
                    '<pre style="margin:4px 0; white-space:pre-wrap;">'
                    'Host REMOTE_ALIAS\n'
                    '    HostName REMOTE_HOST\n'
                    '    User YOUR_USERNAME\n'
                    '    ControlMaster auto\n'
                    '    ControlPath ~/.ssh/cm-delfin-%r@%h:%p\n'
                    '    ControlPersist yes'
                    '</pre>'
                    'Save with <code>Ctrl+O</code>, <code>Enter</code>, <code>Ctrl+X</code>. '
                    'Then set the correct permissions:<br>'
                    '<pre style="margin:4px 0; white-space:pre-wrap;">'
                    'chmod 600 ~/.ssh/config'
                    '</pre><br>'
                    '<b>Step 3 &ndash; Open the master connection:</b><br>'
                    'Connect once manually. Enter your TOTP code and password when prompted.<br>'
                    '<pre style="margin:4px 0; white-space:pre-wrap;">'
                    'ssh REMOTE_ALIAS'
                    '</pre>'
                    'You should now be logged into the remote server. '
                    'You can type <code>exit</code> &ndash; the master session stays open '
                    'in the background.<br><br>'
                    '<b>Step 4 &ndash; Verify:</b><br>'
                    'This should connect <b>instantly</b> without asking for a password:<br>'
                    '<pre style="margin:4px 0; white-space:pre-wrap;">'
                    'ssh REMOTE_ALIAS hostname'
                    '</pre>'
                    'If it prints the remote hostname without a prompt, everything works.<br><br>'
                    '<b>Step 5 &ndash; DELFIN settings:</b><br>'
                    'In the transfer settings above, enter <code>REMOTE_ALIAS</code> as the '
                    '<b>Host</b> (the same name you used in ~/.ssh/config).<br><br>'
                    'The master session stays open until the remote server disconnects or '
                    'you close it manually with <code>ssh -O exit REMOTE_ALIAS</code>.<br>'
                    'If your local key has a passphrase, load it with '
                    '<code>ssh-add ~/.ssh/id_ed25519</code>.'
                    '</div>'
                )
            ),
        ],
    )
    ssh_help_box.set_title(0, 'SSH Setup Guide (click to expand)')
    ssh_help_box.selected_index = None  # collapsed by default

    transfer_section = widgets.VBox(
        [
            widgets.HBox(
                [
                    widgets.HTML('<b>Remote Archive</b>'),
                    remote_archive_toggle,
                    widgets.HTML(
                        '<span style="color:#616161;">'
                        'Controls the <code>Remote Archive</code> tab and remote buttons in <code>Archive</code>.'
                        '</span>'
                    ),
                ],
                layout=_row_layout,
            ),
            widgets.HBox(
                [
                    widgets.HTML('<b>Host</b>'),
                    widgets.Box([host_hidden, host_visible], layout=widgets.Layout(width='230px')),
                    widgets.HTML('<b>User</b>'),
                    widgets.Box([user_hidden, user_visible], layout=widgets.Layout(width='180px')),
                    widgets.HTML('<b>Port</b>'),
                    port_input,
                ],
                layout=_row_layout,
            ),
            widgets.HBox(
                [
                    widgets.HTML('<b>Remote Path</b>'),
                    widgets.Box(
                        [path_hidden, path_visible],
                        layout=widgets.Layout(flex='1 1 420px', min_width='280px'),
                    ),
                ],
                layout=_row_layout,
            ),
            ssh_help_box,
        ],
        layout=_section_layout,
    )

    # ══════════════════════════════════════════════════════════════════════
    # Main accordion
    # ══════════════════════════════════════════════════════════════════════
    main_accordion = widgets.Accordion(
        children=[
            workspace_section,
            tabs_section,
            runtime_section,
            tools_section,
            developer_section,
            transfer_section,
        ],
    )
    main_accordion.set_title(0, 'Workspace Paths')
    main_accordion.set_title(1, 'Dashboard Tabs')
    main_accordion.set_title(2, 'Runtime Backend')
    main_accordion.set_title(3, 'Tool Installation')
    main_accordion.set_title(4, 'Developer')
    main_accordion.set_title(5, 'SSH Transfer & Remote Archive')
    main_accordion.selected_index = 0  # Workspace open by default

    tab = widgets.VBox(
        [
            info_box,
            top_bar,
            main_accordion,
        ],
        layout=widgets.Layout(width='100%', gap='12px', padding='10px'),
    )

    _apply_sensitive_visibility()
    _load_settings_to_widgets(set_status=True)
    _refresh_qm_status()
    _refresh_csp_status()
    _refresh_mlp_status()
    _refresh_analysis_status()
    _refresh_ai_status()
    save_btn.disabled = True  # clean state after initial load

    return tab, {'reload_settings': _load_settings_to_widgets}
