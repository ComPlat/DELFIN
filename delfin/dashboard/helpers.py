"""Shared helper / utility functions for the DELFIN Dashboard."""

import base64
import csv
import io
import re

import ipywidgets as widgets

from .constants import JOB_TIME_LIMITS


_NEB_TRAJECTORY_SECTION_RE = re.compile(r'^\s*(Images|Interp\.)\s*:\s*(.*)$')
_NEB_TRAJECTORY_ROW_RE = re.compile(
    r'^\s*([-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?)'
    r'\s+([-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?)'
    r'\s+([-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?)\s*$'
)
_BOHR_TO_ANGSTROM = 0.529177210903
_HARTREE_TO_KCAL_MOL = 627.5094740631


def parse_neb_final_interp(text):
    """Parse a ``*.final.interp`` file into Images / Interp. sections."""
    sections = {}
    current_key = None
    for raw_line in str(text or '').splitlines():
        match = _NEB_TRAJECTORY_SECTION_RE.match(raw_line)
        if match:
            current_key = match.group(1)
            sections[current_key] = []
            continue
        if current_key is None:
            continue
        row_match = _NEB_TRAJECTORY_ROW_RE.match(raw_line)
        if not row_match:
            if raw_line.strip():
                current_key = None
            continue
        try:
            progress = float(row_match.group(1).replace('D', 'E').replace('d', 'e'))
            distance_bohr = float(row_match.group(2).replace('D', 'E').replace('d', 'e'))
            energy_eh = float(row_match.group(3).replace('D', 'E').replace('d', 'e'))
        except ValueError:
            continue
        sections.setdefault(current_key, []).append(
            {
                'progress': progress,
                'distance_bohr': distance_bohr,
                'distance_angstrom': distance_bohr * _BOHR_TO_ANGSTROM,
                'energy_eh': energy_eh,
                'energy_kcal_mol': energy_eh * _HARTREE_TO_KCAL_MOL,
            }
        )
    return sections


def _build_neb_trajectory_figure(text, title='Trajectory Plot', energy_unit='kcal/mol'):
    """Build a matplotlib figure for ``*.final.interp`` content."""
    sections = parse_neb_final_interp(text)
    interp_rows = sections.get('Interp.', [])
    image_rows = sections.get('Images', [])
    rows = interp_rows or image_rows
    if not rows:
        raise ValueError('No trajectory data found in .final.interp file.')
    energy_unit_norm = str(energy_unit or 'kcal/mol').strip().lower()
    if energy_unit_norm in {'eh', 'hartree'}:
        energy_key = 'energy_eh'
        energy_label = 'Energy (Hartree)'
    else:
        energy_key = 'energy_kcal_mol'
        energy_label = 'Energy (kcal/mol)'

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.4, 4.8), dpi=180)

    if interp_rows:
        ax.plot(
            [row['distance_angstrom'] for row in interp_rows],
            [row[energy_key] for row in interp_rows],
            color='#1565c0',
            linewidth=2.2,
            label='Interp.',
        )
    if image_rows:
        ax.scatter(
            [row['distance_angstrom'] for row in image_rows],
            [row[energy_key] for row in image_rows],
            color='#ef6c00',
            edgecolors='white',
            linewidths=0.7,
            s=38,
            zorder=3,
            label='Images',
        )
        ax.plot(
            [row['distance_angstrom'] for row in image_rows],
            [row[energy_key] for row in image_rows],
            color='#ef6c00',
            linewidth=1.0,
            alpha=0.45,
            zorder=2,
        )

    ax.set_xlabel('Distance (Angstrom)')
    ax.set_ylabel(energy_label)
    ax.set_title(str(title or 'Trajectory Plot'))
    ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.7)
    if interp_rows and image_rows:
        ax.legend(frameon=False)

    peak_row = max(rows, key=lambda row: row[energy_key])
    end_row = rows[-1]
    y_values = [row[energy_key] for row in rows]
    y_min = min(y_values)
    y_max = max(y_values)
    y_span = y_max - y_min
    if y_span <= 0:
        y_span = max(abs(y_max), 1.0) * 0.12
    top_padding = max(y_span * 0.16, 0.02)
    bottom_padding = max(y_span * 0.06, 0.01)
    ax.set_ylim(y_min - bottom_padding, y_max + top_padding)
    peak_label = f"{peak_row[energy_key]:.2f}"
    end_label = f"{end_row[energy_key]:.2f}"

    ax.annotate(
        peak_label,
        xy=(peak_row['distance_angstrom'], peak_row[energy_key]),
        xytext=(0, 10),
        textcoords='offset points',
        fontsize=9,
        color='#b71c1c',
        bbox={'boxstyle': 'round,pad=0.18', 'facecolor': 'white', 'edgecolor': '#ef9a9a', 'alpha': 0.95},
        ha='center',
        va='bottom',
    )
    ax.annotate(
        end_label,
        xy=(end_row['distance_angstrom'], end_row[energy_key]),
        xytext=(0, 10),
        textcoords='offset points',
        fontsize=9,
        color='#1b5e20',
        bbox={'boxstyle': 'round,pad=0.18', 'facecolor': 'white', 'edgecolor': '#a5d6a7', 'alpha': 0.95},
        ha='center',
        va='bottom',
    )
    fig.tight_layout()

    return fig


def save_neb_trajectory_plot_png(text, output_path, title='Trajectory Plot', energy_unit='kcal/mol'):
    """Save ``*.final.interp`` content as a PNG plot."""
    fig = _build_neb_trajectory_figure(text, title=title, energy_unit=energy_unit)
    try:
        fig.savefig(output_path, format='png', bbox_inches='tight')
    finally:
        import matplotlib.pyplot as plt
        plt.close(fig)


def _format_csv_number(value, decimal=','):
    text = f'{float(value):.10f}'.rstrip('0').rstrip('.')
    if decimal == ',':
        return text.replace('.', ',')
    return text


def save_neb_trajectory_csv(text, output_path, decimal='.'):
    """Save parsed ``*.final.interp`` trajectory data as CSV."""
    sections = parse_neb_final_interp(text)
    rows_written = 0
    with open(output_path, 'w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.writer(handle, delimiter=';')
        writer.writerow([
            'section',
            'progress',
            'distance_bohr',
            'distance_angstrom',
            'energy_eh',
            'energy_kcal_mol',
        ])
        for section_name in ('Images', 'Interp.'):
            for row in sections.get(section_name, []):
                writer.writerow([
                    section_name,
                    _format_csv_number(row['progress'], decimal=decimal),
                    _format_csv_number(row['distance_bohr'], decimal=decimal),
                    _format_csv_number(row['distance_angstrom'], decimal=decimal),
                    _format_csv_number(row['energy_eh'], decimal=decimal),
                    _format_csv_number(row['energy_kcal_mol'], decimal=decimal),
                ])
                rows_written += 1
    if rows_written <= 0:
        raise ValueError('No trajectory data found in .final.interp file.')


def render_neb_trajectory_plot_html(text, title='Trajectory Plot'):
    """Render ``*.final.interp`` content as an inline PNG plot and return HTML."""
    fig = _build_neb_trajectory_figure(text, title=title)
    buffer = io.BytesIO()
    fig.savefig(buffer, format='png', bbox_inches='tight')
    import matplotlib.pyplot as plt
    plt.close(fig)
    b64 = base64.b64encode(buffer.getvalue()).decode('ascii')
    return (
        "<div style='width:100%; border:1px solid #d9dee3; border-radius:6px; background:#fff; "
        "padding:10px; box-sizing:border-box;'>"
        f"<img src='data:image/png;base64,{b64}' "
        "style='display:block; width:100%; height:auto; max-width:100%;' />"
        "</div>"
    )


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
        ".delfin-time-limit-toggle .widget-toggle-button:nth-of-type(n+4), "
        ".delfin-time-limit-toggle.widget-toggle-buttons "
        ".widget-toggle-button:nth-of-type(n+4) { "
        "margin-top: 6px !important; }"
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
    toggle.add_class('delfin-time-limit-toggle')

    custom = widgets.BoundedIntText(
        value=72, min=1, max=720, step=1,
        description='Hours:',
        layout=widgets.Layout(width='180px', display='none', margin='6px 0 0 0'),
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


def _append_js(ctx, script):
    """Append JavaScript to the shared output without clearing prior JS."""
    if not script:
        return
    from IPython.display import Javascript, display
    with ctx.js_output:
        display(Javascript(script))


def disable_spellcheck(ctx, class_name='delfin-nospell'):
    """Disable spellcheck for textareas with the given CSS class.

    Note: ipywidgets attaches custom classes to the widget container, not
    the underlying form controls, so we target ".class textarea" and
    ".class input".
    """
    script = f"""
    (function() {{
        window.__delfinSpellcheckClasses = window.__delfinSpellcheckClasses || [];
        if (window.__delfinSpellcheckClasses.indexOf('{class_name}') === -1) {{
            window.__delfinSpellcheckClasses.push('{class_name}');
        }}
        function disableInputLike(el) {{
            if (!el) return;
            const tag = (el.tagName || '').toUpperCase();
            if (tag !== 'TEXTAREA' && tag !== 'INPUT') return;
            el.setAttribute('spellcheck', 'false');
            el.spellcheck = false;
            el.setAttribute('autocorrect', 'off');
            el.setAttribute('autocapitalize', 'off');
            el.setAttribute('autocomplete', 'off');
            el.setAttribute('data-gramm', 'false');
            el.setAttribute('data-gramm_editor', 'false');
            el.setAttribute('data-enable-grammarly', 'false');
            el.setAttribute('data-lt-active', 'false');
            el.setAttribute('data-grammarly', 'false');
            el.setAttribute('data-ms-editor', 'false');
            if (!el.hasAttribute('lang')) el.setAttribute('lang', '');
        }}
        function scan(root) {{
            if (!root || !root.querySelectorAll) return;
            window.__delfinSpellcheckClasses.forEach(function(cls) {{
                root.querySelectorAll('.' + cls + ' textarea, .' + cls + ' input')
                    .forEach(disableInputLike);
            }});
        }}
        scan(document);
        if (!window.__delfinSpellcheckFocusBound) {{
            window.__delfinSpellcheckFocusBound = true;
            document.addEventListener('focusin', function(ev) {{
                const t = ev && ev.target;
                disableInputLike(t);
            }}, true);
        }}
        if (window.__delfinSpellcheckObserver && window.__delfinSpellcheckObserver.disconnect) {{
            try {{ window.__delfinSpellcheckObserver.disconnect(); }} catch (_e) {{}}
        }}
        const obs = new MutationObserver(function(muts) {{
            muts.forEach(function(m) {{
                if (m.type === 'attributes') {{
                    disableInputLike(m.target);
                    return;
                }}
                m.addedNodes && m.addedNodes.forEach(function(node) {{
                    if (!node || node.nodeType !== 1) return;
                    if (node.matches && (node.matches('textarea') || node.matches('input'))) {{
                        disableInputLike(node);
                    }}
                    scan(node);
                }});
            }});
        }});
        obs.observe(document.documentElement, {{
            childList: true,
            subtree: true,
            attributes: true,
            attributeFilter: ['spellcheck', 'autocorrect', 'autocapitalize', 'autocomplete'],
        }});
        window.__delfinSpellcheckObserver = obs;
    }})();
    """
    _append_js(ctx, script)


def disable_spellcheck_global(ctx):
    """Disable spellcheck for all textareas in the dashboard."""
    script = """
    (function() {
        if (window.__delfinSpellcheckObserverAll) return;
        function disableTextarea(el) {
            if (!el || el.tagName !== 'TEXTAREA') return;
            el.setAttribute('spellcheck', 'false');
            el.spellcheck = false;
        }
        function scan(root) {
            if (!root || !root.querySelectorAll) return;
            root.querySelectorAll('textarea').forEach(disableTextarea);
        }
        scan(document);
        const obs = new MutationObserver(function(muts) {
            muts.forEach(function(m) {
                m.addedNodes && m.addedNodes.forEach(function(node) {
                    if (!node || node.nodeType !== 1) return;
                    if (node.matches && node.matches('textarea')) {
                        disableTextarea(node);
                    }
                    scan(node);
                });
            });
        });
        obs.observe(document.documentElement, { childList: true, subtree: true });
        window.__delfinSpellcheckObserverAll = obs;
    })();
    """
    _append_js(ctx, script)


def apply_branding(ctx, title='DELFIN Dashboard', favicon_data_uri=''):
    """Set the browser title and favicon for the dashboard page."""
    safe_title = (
        str(title)
        .replace("\\", "\\\\")
        .replace("'", "\\'")
        .replace("\n", " ")
    )
    safe_icon = (
        str(favicon_data_uri)
        .replace("\\", "\\\\")
        .replace("'", "\\'")
        .replace("\n", "")
    )
    script = f"""
    (function() {{
        document.title = '{safe_title}';
        if (!window.__delfinBrandingObserver) {{
            window.__delfinBrandingObserver = new MutationObserver(function() {{
                if (document.title !== '{safe_title}') {{
                    document.title = '{safe_title}';
                }}
            }});
            window.__delfinBrandingObserver.observe(document.head || document.documentElement, {{
                childList: true,
                subtree: true,
            }});
        }}
        if (!window.__delfinSetFavicon) {{
            window.__delfinSetFavicon = function() {{
                if (!'{safe_icon}') return;
                var head = document.head || document.getElementsByTagName('head')[0];
                if (!head) return;
                var link = document.querySelector("link[rel='icon'][data-delfin-favicon='1']");
                if (!link) {{
                    link = document.createElement('link');
                    link.setAttribute('rel', 'icon');
                    link.setAttribute('data-delfin-favicon', '1');
                    head.appendChild(link);
                }}
                if (link.getAttribute('href') !== '{safe_icon}') {{
                    link.setAttribute('href', '{safe_icon}');
                    link.setAttribute('type', 'image/png');
                }}
            }};
        }}
        window.__delfinSetFavicon();
    }})();
    """
    _append_js(ctx, script)
