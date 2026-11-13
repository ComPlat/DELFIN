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


def js_string_literal(value):
    """A JavaScript string literal that stays a string, whatever it holds.

    File content (an .xyz comment line, a molblock, a .cube header) ends
    up inside ``<script>`` blocks the dashboard renders. ``json.dumps``
    alone is not enough there: it leaves ``<`` and ``>`` alone, so a
    ``</script>`` inside the content closes the script element early and
    whatever follows is parsed as markup. A template literal is worse
    still: a backtick or ``${`` in the content is executed as code.

    Every ``<``, ``>`` and ``&`` is written as a unicode escape (backslash-u and four hex digits),
    which JavaScript reads back as the same character while the HTML
    parser never sees it. ``ensure_ascii`` (the default) already covers
    U+2028/U+2029 and everything else outside ASCII.
    """
    import json
    return (json.dumps(str(value if value is not None else ""))
            .replace("<", "\\u003c")
            .replace(">", "\\u003e")
            .replace("&", "\\u0026"))


_HREF_RE = re.compile(r'(\bhref\s*=\s*)(["\'])(.*?)\2', re.IGNORECASE | re.DOTALL)
_SAFE_SCHEMES = ('http://', 'https://', 'mailto:', '#')
_IMAGE_TYPE_RE = re.compile(r'^image/[A-Za-z0-9.+-]+$')


def neutralise_link_targets(markup):
    """Replace every ``href`` that is not http(s), mailto or a fragment by ``#``.

    For HTML another program converted from a user's document (mammoth's
    .docx conversion): text there is escaped, but a hyperlink target is
    copied as written, and a ``javascript:`` target would run on the
    dashboard's origin when clicked. The target is compared after HTML
    entity decoding and with whitespace and control characters removed,
    the forms browsers ignore inside a scheme.
    """
    import html as _html_mod

    def _fix(match):
        target = _html_mod.unescape(match.group(3))
        target = ''.join(ch for ch in target if ch > ' ').lower()
        if target.startswith(_SAFE_SCHEMES):
            return match.group(0)
        return f'{match.group(1)}{match.group(2)}#{match.group(2)}'

    return _HREF_RE.sub(_fix, str(markup or ''))


def image_data_uri(content_type, b64):
    """A ``data:`` image URI, or '' when *content_type* is not an image type."""
    content_type = str(content_type or '').strip()
    if not _IMAGE_TYPE_RE.match(content_type):
        return ''
    return f'data:{content_type};base64,{b64}'


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
            if hasattr(custom_widget, 'hours_widget'):
                hours = int(custom_widget.hours_widget.value)
                minutes = int(custom_widget.minutes_widget.value)
                seconds = int(custom_widget.seconds_widget.value)
            else:
                # Backward-compat: bare hours-only widget
                hours = int(custom_widget.value)
                minutes = 0
                seconds = 0
            if hours < 0 or not (0 <= minutes <= 59) or not (0 <= seconds <= 59):
                raise ValueError
            if hours == 0 and minutes == 0 and seconds == 0:
                raise ValueError
            return f'{hours}:{minutes:02d}:{seconds:02d}'
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


#: The class on the dashboard's top-level Tab widget (create_page_css).
DASHBOARD_TABS_CLASS = 'delfin-dashboard-tabs'


def create_page_css():
    """Return an HTML widget with the page-level layout CSS.

    Voila pads the notebook and every rendered cell vertically. Those few
    dozen pixels sit on top of the tab bodies, which size themselves as
    ``calc(100vh - 145px)``, so the page as a whole ends up taller than the
    viewport and grows its own thin scrollbar. Dropping the vertical padding
    lifts the header to the top edge and lets the tabs fit inside the window.
    """
    return widgets.HTML(
        "<style>"
        # Notebook / cell container: no vertical padding, keep the horizontal
        # padding so the content is not glued to the window edge.
        ".jp-Notebook, #rendered_cells { "
        "padding-top: 0 !important; padding-bottom: 0 !important; }"
        # Every rendered cell (header cell and tab cell).
        "#rendered_cells .jp-Cell { "
        "padding-top: 0 !important; padding-bottom: 0 !important; "
        "margin-top: 0 !important; margin-bottom: 0 !important; }"
        # JupyterLab puts a 4px gap above each output area.
        "#rendered_cells .jp-Cell-outputWrapper { margin-top: 0 !important; }"
        # An inactive tab's label is the theme's primary text colour, as in
        # JupyterLab's own tab bars. The widgets' secondary colour on the
        # bar's grey falls below the 4.5:1 contrast WCAG AA asks of text.
        ".jupyter-widgets.jupyter-widget-tab > .lm-TabBar "
        ".lm-TabBar-tab:not(.lm-mod-current) { "
        "color: var(--jp-ui-font-color1); }"
        + DASHBOARD_TABS_CSS
        + DASHBOARD_FORM_CSS
        + "</style>"
    )


#: One look for the widgets' own controls on every tab: 6 px corners like
#: the header's (8 px there, at 30 px height), the accent border and a soft
#: ring on focus, a visible keyboard focus on buttons. The widgets draw
#: square fields and selects, 2 px buttons, and mark focus by a border
#: colour only.
#:
#: Precedence, by construction: the selectors are the widgets' own, so the
#: specificity is theirs (focus is added through ``:where``, which counts
#: nothing), and nothing is ``!important``. The widgets' styles sit in the
#: document head, before this; every tab's own CSS sits after it (this is
#: in the header cell), so a rule written for one control -- the agent's
#: composer and chips, the session list, the structure editor -- wins at
#: equal specificity and keeps its look.
#:
#: Rejected: per-widget layouts (dozens of inline widths and radii to keep
#: in step); one shared class added to every control (each tab would have
#: to remember it, and the ones that forget look different again).
_FIELDS = (".widget-text input[type='text'], .widget-text input[type='number'], "
           ".widget-text input[type='password'], .widget-textarea textarea, "
           ".widget-dropdown > select")
_FOCUSED = ", ".join(s.strip() + ":where(:focus)" for s in _FIELDS.split(","))
DASHBOARD_FORM_CSS = (
    f"{_FIELDS} {{ border-radius: 6px; }}"
    f"{_FOCUSED} {{ border-color: var(--delfin-accent, #3d74b6); "
    "box-shadow: 0 0 0 2px rgba(61, 116, 182, 0.18); }"
    ".jupyter-button { border-radius: 6px; }"
    ".jupyter-button:where(:focus-visible) { outline: 2px solid var(--delfin-accent, #3d74b6); "
    "outline-offset: 1px; }"
)


_BAR = f".jupyter-widgets.jupyter-widget-tab.{DASHBOARD_TABS_CLASS} > .lm-TabBar"

#: The dashboard's own tab bar (DASHBOARD_TABS_CLASS): one flat row.
#:
#: * One row, every tab as wide as its name. The widgets' default shrank
#:   the tabs and clipped the labels ("Submit Jo"); wrapping onto a second
#:   row (2026-11-12) left "Settings" alone below the rest at 1366 px and
#:   moved every tab body down by a row. Where the names do not fit, the
#:   strip scrolls sideways (wheel, trackpad, keyboard focus), its
#:   scrollbar hidden and the side that has more faded
#:   (keep_dashboard_tab_in_view sets the two classes).
#: * The padding follows the window: clamp(6px, 0.6vw, 14px) fits the
#:   eighteen labels (about 1020 px of text at 13 px) in one row at
#:   1366 px and gives them room at 1920 px.
#: * Flat tabs as in current IDEs and web apps: no grey boxes; inactive
#:   labels in the theme's secondary text colour (4.6:1 on the white bar,
#:   WCAG AA -- the 2026-11-12 contrast finding was the secondary colour on
#:   the boxes' grey), the current one in the primary colour with an accent
#:   underline that lies on the body's top border, a light background on
#:   hover. The current tab is not bold: a bolder label would widen it and
#:   shift every tab after it.
#: * 28 px high, as the old single row (26 px plus the overlap): the tab
#:   bodies are calc(100vh - 145px) tall.
#:
#: Rejected: an overflow menu ("More") for the tabs that do not fit. The
#: tab bar is Lumino's, re-rendered from its own state; hiding tabs in it
#: and mirroring them in a menu would fight that renderer on every update.
DASHBOARD_TABS_CSS = (
    ":root { --delfin-accent: #3d74b6; }"
    f"{_BAR} {{ min-height: 0; }}"
    f"{_BAR} > .lm-TabBar-content {{ flex-wrap: nowrap; align-items: flex-end; "
    "gap: 2px; overflow-x: auto; overflow-y: hidden; scrollbar-width: none; "
    "overscroll-behavior-x: contain; position: relative; z-index: 1; "
    "margin-bottom: -1px; }"
    f"{_BAR} > .lm-TabBar-content::-webkit-scrollbar {{ display: none; }}"
    f"{_BAR} .lm-TabBar-tab {{ flex: 0 0 auto; min-width: 0; max-width: none; "
    "margin: 0; padding: 0 clamp(6px, 0.6vw, 14px); height: 28px; min-height: 28px; "
    "box-sizing: border-box; border: none; border-bottom: 2px solid transparent; "
    "border-radius: 6px 6px 0 0; background: transparent; font-size: 13px; "
    "cursor: pointer; transition: color .12s, background-color .12s; }"
    f"{_BAR} .lm-TabBar-tabLabel {{ overflow: visible; line-height: 26px; }}"
    f"{_BAR} .lm-TabBar-tab:not(.lm-mod-current) {{ "
    "color: var(--jp-ui-font-color2); }"
    f"{_BAR} .lm-TabBar-tab:hover:not(.lm-mod-current) {{ "
    "background: var(--jp-layout-color2); color: var(--jp-ui-font-color1); }"
    f"{_BAR} .lm-TabBar-tab.lm-mod-current {{ color: var(--jp-ui-font-color0); "
    "background: transparent; border-bottom-color: var(--delfin-accent, #3d74b6); "
    "min-height: 28px; transform: none; }"
    f"{_BAR} .lm-TabBar-tab.lm-mod-current::before {{ display: none; }}"
    f"{_BAR} .lm-TabBar-tab:focus-visible {{ outline: 2px solid var(--delfin-accent, #3d74b6); "
    "outline-offset: -2px; }"
    f"{_BAR}.delfin-tabs-more-right > .lm-TabBar-content {{ "
    "-webkit-mask-image: linear-gradient(to right, #000 calc(100% - 40px), transparent); "
    "mask-image: linear-gradient(to right, #000 calc(100% - 40px), transparent); }"
    f"{_BAR}.delfin-tabs-more-left > .lm-TabBar-content {{ "
    "-webkit-mask-image: linear-gradient(to left, #000 calc(100% - 40px), transparent); "
    "mask-image: linear-gradient(to left, #000 calc(100% - 40px), transparent); }"
    f"{_BAR}.delfin-tabs-more-left.delfin-tabs-more-right > .lm-TabBar-content {{ "
    "-webkit-mask-image: linear-gradient(to right, transparent, #000 40px, "
    "#000 calc(100% - 40px), transparent); "
    "mask-image: linear-gradient(to right, transparent, #000 40px, "
    "#000 calc(100% - 40px), transparent); }"
)


def keep_dashboard_tab_in_view(ctx):
    """The dashboard's tab strip, where it scrolls: edges, wheel, current tab.

    Input: the page. Output, for each tab bar of class DASHBOARD_TABS_CLASS
    (also one rendered later, e.g. after a session comes back):

    * ``delfin-tabs-more-left`` / ``delfin-tabs-more-right`` on the bar
      while tabs lie beyond that edge (DASHBOARD_TABS_CSS fades it);
    * a vertical wheel over the strip scrolls it sideways, only while it
      overflows -- otherwise the wheel stays the page's;
    * the current tab is scrolled into view when it changes (a tab opened
      from code, e.g. the agent switching to Calculations).

    Events and observers only, no timer. Rejected: ``scrollIntoView`` --
    it also scrolls every scrollable ancestor, the page included.
    """
    script = """
    (function() {
        if (window.__delfinTabStrip) return;
        window.__delfinTabStrip = true;
        var BAR = '.jupyter-widget-tab.%(cls)s > .lm-TabBar';
        function strip(bar) { return bar.querySelector(':scope > .lm-TabBar-content'); }
        function edges(bar) {
            var c = strip(bar);
            if (!c) return;
            var max = c.scrollWidth - c.clientWidth;
            bar.classList.toggle('delfin-tabs-more-left', max > 1 && c.scrollLeft > 1);
            bar.classList.toggle('delfin-tabs-more-right', max > 1 && c.scrollLeft < max - 1);
        }
        function reveal(bar) {
            var c = strip(bar);
            if (!c) return;
            var cur = c.querySelector('.lm-TabBar-tab.lm-mod-current');
            if (cur && c.scrollWidth > c.clientWidth) {
                var left = cur.offsetLeft - c.offsetLeft, right = left + cur.offsetWidth;
                if (left < c.scrollLeft) c.scrollLeft = Math.max(0, left - 40);
                else if (right > c.scrollLeft + c.clientWidth)
                    c.scrollLeft = right - c.clientWidth + 40;
            }
            edges(bar);
        }
        var watched = new WeakSet();
        var changes = new MutationObserver(function(records) {
            var done = new Set();
            records.forEach(function(r) {
                var bar = r.target.closest && r.target.closest(BAR);
                if (bar && !done.has(bar)) { done.add(bar); reveal(bar); }
            });
        });
        function watch(bar) {
            var c = strip(bar);
            if (!c || watched.has(c)) return;
            watched.add(c);
            c.addEventListener('scroll', function() { edges(bar); }, {passive: true});
            c.addEventListener('wheel', function(ev) {
                if (c.scrollWidth <= c.clientWidth) return;
                if (Math.abs(ev.deltaY) <= Math.abs(ev.deltaX)) return;
                c.scrollLeft += ev.deltaY;
                ev.preventDefault();
            }, {passive: false});
            changes.observe(c, {attributes: true, attributeFilter: ['class'],
                                childList: true, subtree: true});
            reveal(bar);
        }
        function scan(root) {
            if (!root || root.nodeType !== 1) return;
            if (root.matches(BAR)) watch(root);
            root.querySelectorAll(BAR).forEach(watch);
        }
        scan(document.body);
        new MutationObserver(function(records) {
            for (var i = 0; i < records.length; i++) {
                var added = records[i].addedNodes;
                for (var j = 0; j < added.length; j++) scan(added[j]);
            }
        }).observe(document.body, {childList: true, subtree: true});
        window.addEventListener('resize', function() {
            document.querySelectorAll(BAR).forEach(reveal);
        }, {passive: true});
    })();
    """ % {"cls": DASHBOARD_TABS_CLASS}
    ctx.add_init_js(script)


def name_page_controls(ctx):
    """Give the page's unnamed controls an accessible name, and tabs a tooltip.

    A widget's description is a ``<label>`` beside the control, and only the
    dropdown links it to its ``<select>`` -- and then only when the
    description is not empty. A slider's handle (``role="slider"``) is never
    linked to anything. Screen readers announced them as "combo box" and
    "slider" with no name, and an accessibility check (axe, WCAG 2 AA)
    flagged both. Each control without an ``aria-label`` takes the words its
    widget shows: the description, else the widget's tooltip.

    A tab bar clips a long name, and the widgets give a tab no tooltip. Each
    tab gets its full name as ``title``, as JupyterLab's tabs have; a tab
    renamed in place gets it again when the pointer or the focus reaches it.

    It works on what is added to the page as it is added (one observer, no
    timer), so a tab built on its first open is covered like the first one.
    """
    script = """
    (function() {
        if (window.__delfinControlNames) return;
        window.__delfinControlNames = true;
        var CONTROLS = 'select:not([aria-label]):not([aria-labelledby]), '
            + '.noUi-handle:not([aria-label]):not([aria-labelledby])';
        var WIDGET = '.widget-slider, .widget-dropdown, .widget-select, .jupyter-widgets';
        var TAB = '.lm-TabBar-tab';
        function words(el) {
            var widget = el.closest(WIDGET), label = null;
            if (widget) {
                for (var i = 0; i < widget.children.length; i++) {
                    if (widget.children[i].classList.contains('widget-label')) {
                        label = widget.children[i];
                        break;
                    }
                }
            }
            var shown = label ? (label.textContent || '').trim().replace(/:$/, '') : '';
            return shown || (label && label.title) || el.getAttribute('title')
                || (widget && widget.getAttribute('title')) || '';
        }
        function nameControl(el) {
            var name = words(el);
            if (name) el.setAttribute('aria-label', name);
        }
        function titleTab(tab) {
            var label = tab.querySelector('.lm-TabBar-tabLabel');
            var name = label ? (label.textContent || '').trim() : '';
            if (name && tab.getAttribute('title') !== name) tab.setAttribute('title', name);
        }
        function scan(root) {
            if (!root || root.nodeType !== 1) return;
            if (root.matches(CONTROLS)) nameControl(root);
            if (root.matches(TAB)) titleTab(root);
            root.querySelectorAll(CONTROLS).forEach(nameControl);
            root.querySelectorAll(TAB).forEach(titleTab);
        }
        scan(document.body);
        new MutationObserver(function(records) {
            for (var i = 0; i < records.length; i++) {
                var added = records[i].addedNodes;
                for (var j = 0; j < added.length; j++) scan(added[j]);
            }
        }).observe(document.body, {childList: true, subtree: true});
        function refresh(ev) {
            var tab = ev.target && ev.target.closest && ev.target.closest(TAB);
            if (tab) titleTab(tab);
        }
        document.addEventListener('pointerover', refresh, {capture: true, passive: true});
        document.addEventListener('focusin', refresh, true);
    })();
    """
    ctx.add_init_js(script)


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

    hours_widget = widgets.BoundedIntText(
        value=72, min=0, max=720, step=1,
        description='Hours:',
        layout=widgets.Layout(width='150px'),
        style=style,
    )
    minutes_widget = widgets.BoundedIntText(
        value=0, min=0, max=59, step=1,
        description='Min:',
        layout=widgets.Layout(width='120px'),
        style=style,
    )
    seconds_widget = widgets.BoundedIntText(
        value=0, min=0, max=59, step=1,
        description='Sec:',
        layout=widgets.Layout(width='120px'),
        style=style,
    )

    custom = widgets.HBox(
        [hours_widget, minutes_widget, seconds_widget],
        layout=widgets.Layout(display='none', margin='6px 0 0 0'),
    )
    # Expose the sub-widgets so resolve_time_limit and resets can reach them
    # while keeping the (toggle, custom) two-tuple contract intact.
    custom.hours_widget = hours_widget
    custom.minutes_widget = minutes_widget
    custom.seconds_widget = seconds_widget

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


def self_running_html(*scripts):
    """HTML that runs *scripts* when the browser shows it.

    Each script sits as raw text in an inert ``<script type="text/x-delfin">``
    (raw text: no '<' in the code is read as a tag) and an <img> whose load
    fails runs each one on its own through its onerror -- the way the agent
    tab's chat tags run their code. One script that throws stops no other.

    For a tab built after the page went out: its init scripts were sent as a
    Javascript output into the shared JS output widget, and since tabs build
    on first open (09.11.) those outputs were never executed -- measured: the
    Calculations tab's resolvers stayed undefined, its splitter and explorer
    were never bound.
    """
    out = []
    for script in scripts:
        if not script or not str(script).strip():
            continue
        body = str(script).replace("</script", "<\\/script")
        out.append('<script type="text/x-delfin" class="delfin-tab-js">'
                   + body + "</script>")
    if not out:
        return ""
    out.append('<img src="" alt="" style="display:none" onerror="'
               "this.parentNode.querySelectorAll('script.delfin-tab-js')"
               ".forEach(function(t){try{(0,eval)(t.textContent)}"
               "catch(e){console.error('DELFIN tab script',e)}});"
               "this.remove();"
               '">')
    return "".join(out)


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
            // Set an attribute only when it differs: the observer below
            // watches these attributes, and setAttribute writes a mutation
            // record even for the value already there -- unguarded, one
            // change elsewhere fed the observer its own writes forever and
            // froze the page (the Calculations tab, built on first open).
            function put(name, value) {{
                if (el.getAttribute(name) !== value) el.setAttribute(name, value);
            }}
            put('spellcheck', 'false');
            put('autocorrect', 'off');
            put('autocapitalize', 'off');
            put('autocomplete', 'off');
            put('data-gramm', 'false');
            put('data-gramm_editor', 'false');
            put('data-enable-grammarly', 'false');
            put('data-lt-active', 'false');
            put('data-grammarly', 'false');
            put('data-ms-editor', 'false');
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


def debounce_input(ctx, class_name='delfin-debounced', delay_ms=350):
    """Let a text box report its value once typing stops, instead of per keystroke.

    A widget with ``continuous_update=True`` sends every keystroke to the kernel,
    and a handler that renders a 3D viewer or re-parses a file then runs dozens
    of times while a line is being typed -- the box lags behind the keyboard and
    the work is thrown away as soon as the next character arrives.  Setting
    ``continuous_update=False`` alone would go too far the other way: nothing at
    all happens until the box loses focus.

    So the widget is declared ``continuous_update=False`` and this script fires
    its change event *delay_ms* after the last keystroke.  Typing stays smooth,
    the handler runs once, and pasting still updates without clicking elsewhere.
    Add ``class_name`` to the widget to opt in.
    """
    script = f"""
    (function() {{
        window.__delfinDebounceClasses = window.__delfinDebounceClasses || [];
        if (window.__delfinDebounceClasses.indexOf('{class_name}') === -1) {{
            window.__delfinDebounceClasses.push('{class_name}');
        }}
        window.__delfinDebounceDelay = {int(delay_ms)};
        function arm(el) {{
            if (!el || el.__delfinDebounced) return;
            const tag = (el.tagName || '').toUpperCase();
            if (tag !== 'TEXTAREA' && tag !== 'INPUT') return;
            el.__delfinDebounced = true;
            el.addEventListener('input', function() {{
                if (el.__delfinDebounceTimer) clearTimeout(el.__delfinDebounceTimer);
                el.__delfinDebounceTimer = setTimeout(function() {{
                    el.__delfinDebounceTimer = 0;
                    el.dispatchEvent(new Event('change', {{ bubbles: true }}));
                }}, window.__delfinDebounceDelay);
            }});
            /* Leaving the box must not wait for the timer. */
            el.addEventListener('blur', function() {{
                if (el.__delfinDebounceTimer) {{
                    clearTimeout(el.__delfinDebounceTimer);
                    el.__delfinDebounceTimer = 0;
                }}
            }});
        }}
        function scan(root) {{
            if (!root || !root.querySelectorAll) return;
            window.__delfinDebounceClasses.forEach(function(cls) {{
                root.querySelectorAll('.' + cls + ' textarea, .' + cls + ' input')
                    .forEach(arm);
            }});
        }}
        scan(document);
        if (window.__delfinDebounceObserver) return;
        const obs = new MutationObserver(function(muts) {{
            muts.forEach(function(m) {{
                m.addedNodes && m.addedNodes.forEach(function(node) {{
                    if (!node || node.nodeType !== 1) return;
                    scan(node);
                }});
            }});
        }});
        obs.observe(document.documentElement, {{ childList: true, subtree: true }});
        window.__delfinDebounceObserver = obs;
    }})();
    """
    _append_js(ctx, script)


def keep_scroll_position(ctx):
    """Keep the page where it is when widgets are rebuilt under the reader.

    When a container's ``children`` change, the browser removes the old
    widget views at once and renders the new ones a moment later. In that
    moment the page is shorter, the scroll position is clamped, and when
    the content is back the reader is somewhere else -- in Settings a single
    tick moved the page up by 760 px, and Save sent it to the top.

    Measured in the real page: the shrink never reaches a painted frame, so
    a ResizeObserver does not see it; the removal of widget views does, and
    so does a widget view being hidden, which is how a tab switch the reader
    did not ask for shrinks the page. So a scroll upwards that follows either
    within half a second, with no wheel, touch, scroll key or held mouse
    button behind it, is undone as soon as the content is tall enough again.

    Everything the reader does wins: any of those inputs cancels a pending
    correction, and a click on a tab or a section header is left to the
    browser, which is where the page is meant to change. The structure
    editor's scroll hold learnt that a clamp is worse than the jump, so this
    corrects once and never holds.
    """
    script = """
    (function() {
        if (window.__delfinScrollKeeper) return;
        window.__delfinScrollKeeper = true;
        var INTENT_MS = 400, NAV_MS = 1500, REBUILD_MS = 500, WAIT_MS = 3000;
        var NAV = '.lm-TabBar-tab, .lm-Collapse-header, .jupyter-widget-Collapse-header, summary';
        var SCROLL_KEYS = ['PageUp', 'PageDown', 'Home', 'End', 'ArrowUp', 'ArrowDown', ' ', 'Tab'];
        var s = {pane: null, lastTop: 0, held: false, intentAt: -1e9, navAt: -1e9,
                 rebuildAt: -1e9, pending: null};
        function now() { return performance.now(); }
        function busy() {
            return s.held || now() - s.intentAt < INTENT_MS || now() - s.navAt < NAV_MS;
        }
        function intent() { s.intentAt = now(); s.pending = null; }
        function holdsView(node) {
            return node.nodeType === 1 && (node.classList.contains('lm-Widget')
                || !!(node.querySelector && node.querySelector('.lm-Widget')));
        }
        function wasHidden(m) {
            var el = m.target;
            if (!el.classList || !el.classList.contains('lm-Widget')) return false;
            var before = m.oldValue || '';
            if (m.attributeName === 'class') {
                return el.classList.contains('lm-mod-hidden') && before.indexOf('lm-mod-hidden') === -1;
            }
            return el.style.display === 'none' && !/display:\\s*none/.test(before);
        }
        function restore(target) {
            var pane = s.pane, until = now() + WAIT_MS;
            s.pending = target;
            (function step() {
                if (s.pending !== target) return;
                if (busy() || now() > until) {
                    s.pending = null;
                    s.lastTop = pane.scrollTop;
                    return;
                }
                if (pane.scrollHeight - pane.clientHeight >= target - 1) {
                    s.pending = null;
                    pane.scrollTop = target;
                    s.lastTop = target;
                    return;
                }
                requestAnimationFrame(step);
            })();
        }
        function onScroll() {
            var top = s.pane.scrollTop;
            if (top < s.lastTop - 2 && !busy() && now() - s.rebuildAt < REBUILD_MS) {
                if (s.pending === null) restore(s.lastTop);
                return;
            }
            if (s.pending === null) s.lastTop = top;
        }
        function attach(tries) {
            var cells = document.getElementById('rendered_cells');
            if (!cells && tries > 0) { setTimeout(function() { attach(tries - 1); }, 500); return; }
            var pane = cells || document.scrollingElement;
            if (!pane) return;
            s.pane = pane;
            s.lastTop = pane.scrollTop;
            (cells ? pane : window).addEventListener('scroll', onScroll, {passive: true});
            new MutationObserver(function(records) {
                for (var i = 0; i < records.length; i++) {
                    var m = records[i];
                    if (m.type === 'attributes') {
                        if (wasHidden(m)) { s.rebuildAt = now(); return; }
                        continue;
                    }
                    for (var j = 0; j < m.removedNodes.length; j++) {
                        if (holdsView(m.removedNodes[j])) { s.rebuildAt = now(); return; }
                    }
                }
            }).observe(cells ? pane : document.body, {
                childList: true, subtree: true, attributes: true,
                attributeFilter: ['class', 'style'], attributeOldValue: true});
        }
        var opts = {capture: true, passive: true};
        window.addEventListener('wheel', intent, opts);
        window.addEventListener('touchstart', intent, opts);
        window.addEventListener('touchmove', intent, opts);
        window.addEventListener('keydown', function(ev) {
            /* On a control these keys type, choose or press; they scroll
               only the page itself. Tab moves the focus and may scroll. */
            var el = ev.target, tag = (el && el.tagName || '').toUpperCase();
            var control = tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT'
                || tag === 'BUTTON' || (el && el.isContentEditable);
            if (SCROLL_KEYS.indexOf(ev.key) !== -1 && (!control || ev.key === 'Tab')) intent();
        }, opts);
        window.addEventListener('pointerdown', function(ev) {
            s.held = true;
            s.pending = null;
            if (ev.target && ev.target.closest && ev.target.closest(NAV)) s.navAt = now();
        }, opts);
        ['pointerup', 'pointercancel', 'dragend'].forEach(function(type) {
            window.addEventListener(type, function() { s.held = false; }, opts);
        });
        window.addEventListener('blur', function() { s.held = false; });
        attach(20);
    })();
    """
    _append_js(ctx, script)


def hand_over_tabs(tabs_widget, children, titles, selected_index=None):
    """Give *tabs_widget* new children without the selection drifting away.

    The browser rebuilds every tab from the first position that changed,
    and until the rebuilt ones have rendered the selected index can point
    at a view that is about to be removed; the first rebuilt tab to render
    then takes the selection and reports it back. Settings is pinned last,
    so saving a hidden or moved tab sent the reader to another tab.

    So the selection is first parked on the first tab, which the browser
    keeps whenever the two lists share it, and only moved to
    *selected_index* once the new list is in place. When the first tab
    changed too, the first answer from the browser within a few seconds
    that is not *selected_index* is corrected once.
    """
    old = tuple(tabs_widget.children)
    new = tuple(children)
    shared = 0
    while shared < min(len(old), len(new)) and old[shared] is new[shared]:
        shared += 1
    if new != old:
        current = tabs_widget.selected_index
        if shared and current is not None and current >= shared:
            tabs_widget.selected_index = 0
        tabs_widget.children = new
    for i, title in enumerate(titles):
        tabs_widget.set_title(i, title)
    if selected_index is None:
        return
    tabs_widget.selected_index = selected_index
    if new == old or shared:
        return
    import time

    deadline = time.monotonic() + 3.0

    def _correct_once(change):
        tabs_widget.unobserve(_correct_once, names='selected_index')
        if time.monotonic() < deadline and change.get('new') != selected_index:
            tabs_widget.selected_index = selected_index

    tabs_widget.observe(_correct_once, names='selected_index')


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
