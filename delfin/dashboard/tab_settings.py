"""Dashboard Settings tab backed by a user-local settings file."""

import html
import shutil
from pathlib import Path

import ipywidgets as widgets

from delfin.runtime_setup import (
    apply_runtime_environment,
    collect_bwunicluster_verification,
    collect_runtime_diagnostics,
    detect_local_runtime_limits,
    describe_orca_installation,
    discover_orca_installations,
    prepare_bwunicluster_user_setup,
    run_bwunicluster_installer,
    get_packaged_submit_templates_dir,
    get_user_qm_tools_dir,
    run_qm_tools_installer,
    resolve_backend_choice,
    resolve_orca_base,
    resolve_submit_templates_dir,
    stage_packaged_qm_tools,
)
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
        layout=widgets.Layout(width='120px', height='28px'),
    )
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
    qm_tools_root_input = widgets.Text(
        placeholder='Optional qm_tools root override',
        layout=widgets.Layout(width='100%', min_width='280px', height='28px'),
    )
    local_orca_input = widgets.Text(
        placeholder='Optional local ORCA override',
        layout=widgets.Layout(width='100%', min_width='280px', height='28px'),
    )
    local_max_cores_input = widgets.BoundedIntText(
        value=detect_local_runtime_limits()[0],
        min=1,
        max=1_000_000,
        step=1,
        layout=widgets.Layout(width='120px', min_width='120px', height='28px'),
    )
    local_max_ram_input = widgets.BoundedIntText(
        value=detect_local_runtime_limits()[1],
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

    def _refresh_orca_dropdown(selected_value=None):
        current_value = str(selected_value if selected_value is not None else global_orca_input.value or '').strip()
        options = [('Auto-detect / PATH', '')]
        for candidate in discover_orca_installations(
            seed_candidates=_runtime_auto_candidates(),
            search_roots=_orca_search_roots(),
        ):
            options.append((describe_orca_installation(candidate), candidate))

        if current_value and all(current_value != value for _label, value in options):
            options.append((f'Configured: {current_value}', current_value))

        detected_orca_dropdown.options = options
        allowed_values = {value for _label, value in options}
        detected_orca_dropdown.value = current_value if current_value in allowed_values else ''
        detected_orca_dropdown.disabled = len(options) <= 1

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
            status_color = '#2e7d32' if status == 'ok' else '#d32f2f'
            status_label = 'OK' if status == 'ok' else 'Missing'
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
            + '</tbody></table></div>'
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
        local_orca_input.value = str(local_payload.get('orca_base') or '')
        local_max_cores_input.value = int(local_payload.get('max_cores', detected_local_cores))
        local_max_ram_input.value = int(local_payload.get('max_ram_mb', detected_local_ram_mb))
        slurm_orca_input.value = str(slurm_payload.get('orca_base') or '')
        slurm_templates_input.value = str(slurm_payload.get('submit_templates_dir') or '')
        slurm_profile_input.value = str(slurm_payload.get('profile') or '')
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
        )

        if not backend_switch_required:
            ctx.orca_base = effective_orca_base
            if getattr(ctx, 'backend', None) is not None and hasattr(ctx.backend, 'orca_base'):
                ctx.backend.orca_base = effective_orca_base
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

    def _load_settings_to_widgets(set_status=True):
        try:
            settings_payload = load_settings()
        except Exception as exc:
            _set_transfer_widgets({})
            _set_path_widgets({})
            _set_runtime_widgets({})
            _set_remote_archive_widget({})
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

    def _on_toggle_sensitive(change):
        if change.get('name') != 'value':
            return
        _apply_sensitive_visibility()

    def _on_reload(button):
        _load_settings_to_widgets(set_status=True)

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
                    'No ORCA installation was found in PATH or the standard search roots. Keep auto-detect or enter a path manually.',
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
            target, result = run_qm_tools_installer()
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
                if status != "ok":
                    missing += 1
                status_color = "#2e7d32" if status == "ok" else "#d32f2f"
                status_label = "OK" if status == "ok" else "Missing"
                rows.append(
                    "<tr>"
                    f"<td style='padding:4px 8px;'><code>{html.escape(str(item.get('name', '')))}</code></td>"
                    f"<td style='padding:4px 8px; color:{status_color}; font-weight:600;'>{status_label}</td>"
                    f"<td style='padding:4px 8px;'><code>{html.escape(str(item.get('detail', '')))}</code></td>"
                    "</tr>"
                )
            qm_tools_log.value = (
                "<bwUniCluster verification>\n"
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
        _render_runtime_diagnostics(
            settings_payload.get('runtime', {}) or {},
            reload_required=backend_switch_required,
        )
        reload_hint = ' Reload DELFIN to apply Remote Archive tab/button visibility.' if toggle_changed else ''
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
                f'{backend_hint}{reload_hint}'
            ),
            color='#2e7d32',
        )

    show_sensitive_btn.observe(_on_toggle_sensitive, names='value')
    reload_btn.on_click(_on_reload)
    validate_runtime_btn.on_click(_on_validate_runtime)
    scan_orca_btn.on_click(_on_scan_orca)
    detect_local_resources_btn.on_click(_on_detect_local_resources)
    prepare_qm_tools_btn.on_click(_on_prepare_qm_tools)
    install_qm_tools_btn.on_click(_on_install_qm_tools)
    update_qm_tools_btn.on_click(_on_update_qm_tools)
    setup_bwunicluster_btn.on_click(_on_setup_bwunicluster)
    verify_bwunicluster_btn.on_click(_on_verify_bwunicluster)
    full_install_bwunicluster_btn.on_click(_on_full_install_bwunicluster)
    save_btn.on_click(_on_save)
    detected_orca_dropdown.observe(_on_select_detected_orca, names='value')
    global_orca_input.observe(_on_change_global_orca, names='value')

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
    workspace_header = widgets.HTML('<h3 style="margin:0;">Workspace Paths</h3>')
    workspace_rows = widgets.VBox(
        [
            widgets.HBox(
                [
                    widgets.HTML('<b>Calculations</b>'),
                    calc_path_input,
                ],
                layout=widgets.Layout(
                    width='100%',
                    gap='8px',
                    align_items='center',
                    flex_flow='row wrap',
                ),
            ),
            widgets.HBox(
                [
                    widgets.HTML('<b>Archive</b>'),
                    archive_path_input,
                ],
                layout=widgets.Layout(
                    width='100%',
                    gap='8px',
                    align_items='center',
                    flex_flow='row wrap',
                ),
            ),
        ],
        layout=widgets.Layout(width='100%', gap='10px'),
    )
    runtime_header = widgets.HTML('<h3 style="margin:0;">Runtime / Execution</h3>')
    runtime_rows = widgets.VBox(
        [
            widgets.HBox(
                [
                    widgets.HTML('<b>Backend</b>'),
                    backend_dropdown,
                    widgets.HTML(
                        '<span style="color:#616161;">Auto chooses SLURM when <code>sbatch</code> is available.</span>'
                    ),
                ],
                layout=widgets.Layout(
                    width='100%',
                    gap='8px',
                    align_items='center',
                    flex_flow='row wrap',
                ),
            ),
            widgets.HBox(
                [
                    widgets.HTML('<b>ORCA path</b>'),
                    global_orca_input,
                    scan_orca_btn,
                ],
                layout=widgets.Layout(
                    width='100%',
                    gap='8px',
                    align_items='center',
                    flex_flow='row wrap',
                ),
            ),
            widgets.HBox(
                [
                    widgets.HTML('<b>Detected ORCA</b>'),
                    detected_orca_dropdown,
                    widgets.HTML(
                        '<span style="color:#616161;">Choose a detected installation to fill the ORCA path automatically.</span>'
                    ),
                ],
                layout=widgets.Layout(
                    width='100%',
                    gap='8px',
                    align_items='center',
                    flex_flow='row wrap',
                ),
            ),
            widgets.HBox(
                [
                    widgets.HTML('<b>qm_tools root</b>'),
                    qm_tools_root_input,
                ],
                layout=widgets.Layout(
                    width='100%',
                    gap='8px',
                    align_items='center',
                    flex_flow='row wrap',
                ),
            ),
            widgets.HBox(
                [
                    prepare_qm_tools_btn,
                    install_qm_tools_btn,
                    update_qm_tools_btn,
                    widgets.HTML(
                        (
                            f'<span style="color:#616161;">'
                            f'Installs into <code>{html.escape(str(get_user_qm_tools_dir()))}</code>, '
                            'not into the Python package directory. '
                            'ORCA and site modules stay external; PATH/.venv tools are reused.'
                            f'</span>'
                        )
                    ),
                ],
                layout=widgets.Layout(
                    width='100%',
                    gap='8px',
                    align_items='center',
                    flex_flow='row wrap',
                ),
            ),
            qm_tools_log,
            widgets.HTML('<b>Local overrides</b>'),
            widgets.HBox(
                [
                    widgets.HTML('<b>Local ORCA</b>'),
                    local_orca_input,
                ],
                layout=widgets.Layout(
                    width='100%',
                    gap='8px',
                    align_items='center',
                    flex_flow='row wrap',
                ),
            ),
            widgets.HBox(
                [
                    widgets.HTML('<b>Local max cores</b>'),
                    local_max_cores_input,
                    widgets.HTML('<b>Local max RAM (MB)</b>'),
                    local_max_ram_input,
                    detect_local_resources_btn,
                ],
                layout=widgets.Layout(
                    width='100%',
                    gap='8px',
                    align_items='center',
                    flex_flow='row wrap',
                ),
            ),
            widgets.HTML('<b>SLURM / Cluster</b>'),
            widgets.HBox(
                [
                    setup_bwunicluster_btn,
                    verify_bwunicluster_btn,
                    full_install_bwunicluster_btn,
                    widgets.HTML(
                        '<span style="color:#616161;">'
                        '<b>Setup</b> prepares an existing DELFIN install for bwUniCluster. '
                        '<b>Verify</b> is read-only and only reports what is present or missing. '
                        '<b>Full install</b> runs the packaged or repo installer for OpenMPI, repo/venv setup, runtime defaults, and environment wiring. '
                        'ORCA itself still stays external, but an existing ORCA install or tarball is reused when found.'
                        '</span>'
                    ),
                ],
                layout=widgets.Layout(
                    width='100%',
                    gap='8px',
                    align_items='center',
                    flex_flow='row wrap',
                ),
            ),
            widgets.HBox(
                [
                    widgets.HTML('<b>SLURM ORCA</b>'),
                    slurm_orca_input,
                ],
                layout=widgets.Layout(
                    width='100%',
                    gap='8px',
                    align_items='center',
                    flex_flow='row wrap',
                ),
            ),
            widgets.HBox(
                [
                    widgets.HTML('<b>Submit templates</b>'),
                    slurm_templates_input,
                ],
                layout=widgets.Layout(
                    width='100%',
                    gap='8px',
                    align_items='center',
                    flex_flow='row wrap',
                ),
            ),
            widgets.HBox(
                [
                    widgets.HTML('<b>Site profile</b>'),
                    slurm_profile_input,
                ],
                layout=widgets.Layout(
                    width='100%',
                    gap='8px',
                    align_items='center',
                    flex_flow='row wrap',
                ),
            ),
            runtime_diagnostics_html,
        ],
        layout=widgets.Layout(width='100%', gap='10px'),
    )
    transfer_header = widgets.HTML('<h3 style="margin:0;">Transfer Target</h3>')
    remote_archive_row = widgets.HBox(
        [
            widgets.HTML('<b>Remote Archive</b>'),
            remote_archive_toggle,
            widgets.HTML(
                '<span style="color:#616161;">'
                'Controls the <code>Remote Archive</code> tab and remote buttons in <code>Archive</code>.'
                '</span>'
            ),
        ],
        layout=widgets.Layout(
            width='100%',
            gap='8px',
            align_items='center',
            flex_flow='row wrap',
        ),
    )
    transfer_rows = widgets.VBox(
        [
            remote_archive_row,
            widgets.HBox(
                [
                    widgets.HTML('<b>Host</b>'),
                    widgets.Box([host_hidden, host_visible], layout=widgets.Layout(width='230px')),
                    widgets.HTML('<b>User</b>'),
                    widgets.Box([user_hidden, user_visible], layout=widgets.Layout(width='180px')),
                    widgets.HTML('<b>Port</b>'),
                    port_input,
                ],
                layout=widgets.Layout(
                    width='100%',
                    gap='8px',
                    align_items='center',
                    flex_flow='row wrap',
                ),
            ),
            widgets.HBox(
                [
                    widgets.HTML('<b>Remote Path</b>'),
                    widgets.Box(
                        [path_hidden, path_visible],
                        layout=widgets.Layout(flex='1 1 420px', min_width='280px'),
                    ),
                ],
                layout=widgets.Layout(
                    width='100%',
                    gap='8px',
                    align_items='center',
                    flex_flow='row wrap',
                ),
            ),
            widgets.HBox(
                [show_sensitive_btn, reload_btn, validate_runtime_btn, save_btn],
                layout=widgets.Layout(gap='8px', align_items='center', flex_flow='row wrap'),
            ),
            status_html,
        ],
        layout=widgets.Layout(width='100%', gap='10px'),
    )
    future_box = widgets.HTML(
        value=(
            '<div style="color:#616161;">'
            '<b>Future settings:</b> Additional DELFIN options can be added here later '
            'without moving user-specific data into the repository.'
            '</div>'
        )
    )
    ssh_help_box = widgets.HTML(
        value=(
            '<div style="color:#455a64; border:1px solid #d9dee3; border-radius:6px; '
            'padding:10px; background:#fafbfc;">'
            '<b>Set up SSH once in the terminal:</b><br>'
            'To allow <code>SSH Transfer</code> to run as a background job without a password prompt, '
            'set up the SSH key once in the terminal for the target account:<br>'
            '<pre style="margin:8px 0 0 0; white-space:pre-wrap;">'
            'ssh-keygen -t ed25519\n'
            'ssh-copy-id USER@HOST\n'
            'ssh USER@HOST'
            '</pre>'
            'When running <code>ssh-copy-id</code>, enter the target account password once. '
            'If the final <code>ssh</code> test succeeds, DELFIN will work afterward. '
            'If your local key has a passphrase, you can load it when needed with '
            '<code>ssh-add ~/.ssh/id_ed25519</code>.'
            '</div>'
        )
    )

    tab = widgets.VBox(
        [
            info_box,
            workspace_header,
            workspace_rows,
            transfer_header,
            transfer_rows,
            ssh_help_box,
            runtime_header,
            runtime_rows,
            future_box,
        ],
        layout=widgets.Layout(width='100%', gap='12px', padding='10px'),
    )

    _apply_sensitive_visibility()
    _load_settings_to_widgets(set_status=True)

    return tab, {'reload_settings': _load_settings_to_widgets}
