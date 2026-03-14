"""Dashboard Settings tab backed by a user-local settings file."""

import html

import ipywidgets as widgets

from delfin.user_settings import (
    get_settings_path,
    load_transfer_settings,
    save_transfer_settings,
)


def create_tab(ctx):
    """Create the dashboard Settings tab."""
    settings_path = get_settings_path()
    status_html = widgets.HTML(value='')
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
    save_btn = widgets.Button(
        description='Save Settings',
        button_style='primary',
        layout=widgets.Layout(width='120px', height='28px'),
    )

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
        status_html.value = (
            f'<span style="color:{color};">{message}</span>'
        )

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

    def _load_settings_to_widgets(set_status=True):
        try:
            payload = load_transfer_settings()
        except Exception as exc:
            _set_transfer_widgets({})
            _set_status(
                (
                    f'Could not load <code>{html.escape(str(settings_path))}</code>: '
                    f'{html.escape(str(exc))}'
                ),
                color='#d32f2f',
            )
            return None

        _set_transfer_widgets(payload or {})
        if set_status:
            if payload:
                _set_status(
                    (
                        f'Loaded local settings from '
                        f'<code>{html.escape(str(settings_path))}</code>. '
                        'Sensitive values stay masked until revealed.'
                    ),
                    color='#2e7d32',
                )
            else:
                _set_status(
                    (
                        f'No transfer target saved yet. Settings will be created in '
                        f'<code>{html.escape(str(settings_path))}</code>.'
                    ),
                    color='#ef6c00',
                )
        return payload

    def _on_toggle_sensitive(change):
        if change.get('name') != 'value':
            return
        _apply_sensitive_visibility()

    def _on_reload(button):
        _load_settings_to_widgets(set_status=True)

    def _on_save(button):
        try:
            host, user, remote_path, port = _transfer_payload_from_widgets()
            payload = save_transfer_settings(host, user, remote_path, port)
        except Exception as exc:
            _set_status(
                (
                    f'Saving settings failed: {html.escape(str(exc))}. '
                    'No passwords are stored; use SSH keys or ssh-agent for transfers.'
                ),
                color='#d32f2f',
            )
            return

        _set_transfer_widgets(payload)
        _set_status(
            (
                f'Saved transfer target to '
                f'<code>{html.escape(str(payload["settings_path"]))}</code>.'
            ),
            color='#2e7d32',
        )

    show_sensitive_btn.observe(_on_toggle_sensitive, names='value')
    reload_btn.on_click(_on_reload)
    save_btn.on_click(_on_save)

    info_box = widgets.HTML(
        value=(
            f'<div style="color:#455a64;">'
            f'<b>Local settings:</b> <code>{html.escape(str(settings_path))}</code><br>'
            'This file lives in your HOME, outside the git repo. '
            '<code>git pull</code> will not overwrite it. Missing future settings are '
            'added without replacing your existing values.<br>'
            'Passwords are not stored here. Resumable SSH transfers use SSH keys or '
            '<code>ssh-agent</code>.'
            '</div>'
        )
    )
    transfer_header = widgets.HTML('<h3 style="margin:0;">Transfer Target</h3>')
    transfer_rows = widgets.VBox(
        [
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
                [show_sensitive_btn, reload_btn, save_btn],
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
        [info_box, transfer_header, transfer_rows, ssh_help_box, future_box],
        layout=widgets.Layout(width='100%', gap='12px', padding='10px'),
    )

    _apply_sensitive_visibility()
    _load_settings_to_widgets(set_status=True)

    return tab, {'reload_settings': _load_settings_to_widgets}
