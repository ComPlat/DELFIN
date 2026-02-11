"""Shared helper / utility functions for the DELFIN Dashboard."""

import ipywidgets as widgets

from .constants import JOB_TIME_LIMITS


def resolve_time_limit(toggle_widget, custom_widget, default='48:00:00'):
    """Resolve the effective time-limit string from the toggle/custom widgets."""
    toggle_val = toggle_widget.value
    if toggle_val == 'Custom':
        try:
            hours = int(custom_widget.value)
            if hours < 1:
                raise ValueError
            return f'{hours}:00:00'
        except (ValueError, TypeError):
            return default
    return JOB_TIME_LIMITS.get(toggle_val, default)


def create_busy_css():
    """Return an HTML widget containing the busy-spinner CSS."""
    return widgets.HTML(
        "<style>"
        ".delfin-busy { display:inline-block; width:14px; height:14px; "
        "border:2px solid #1976d2; border-top-color: transparent; "
        "border-radius:50%; animation: delfin-spin 0.9s linear infinite; }"
        "@keyframes delfin-spin { from { transform: rotate(0deg); } "
        "to { transform: rotate(360deg); } }"
        ".widget-textarea textarea { resize: none !important; }"
        "textarea { resize: none !important; }"
        "</style>"
    )


def create_time_limit_widgets(style=None):
    """Create and return ``(toggle, custom)`` time-limit widgets.

    The *custom* widget is hidden by default and shown when the user
    selects 'Custom' in the toggle.
    """
    if style is None:
        style = {'description_width': 'initial'}

    toggle = widgets.ToggleButtons(
        options=['12h', '24h', '36h', '48h', '60h', 'Custom'],
        value='48h',
        description='Time Limit:',
        style=style,
        button_style='',
        layout=widgets.Layout(width='550px'),
    )

    custom = widgets.BoundedIntText(
        value=72, min=1, max=720, step=1,
        description='Hours:',
        layout=widgets.Layout(width='180px', display='none'),
        style=style,
    )

    def _toggle_visibility(change):
        custom.layout.display = '' if change['new'] == 'Custom' else 'none'

    toggle.observe(_toggle_visibility, names='value')

    return toggle, custom


def parse_time_to_seconds(time_str):
    """Convert an ``HH:MM:SS`` string to total seconds."""
    parts = time_str.strip().split(':')
    try:
        if len(parts) == 3:
            return int(parts[0]) * 3600 + int(parts[1]) * 60 + int(parts[2])
        elif len(parts) == 2:
            return int(parts[0]) * 60 + int(parts[1])
        else:
            return int(parts[0]) * 3600
    except (ValueError, IndexError):
        return 48 * 3600
