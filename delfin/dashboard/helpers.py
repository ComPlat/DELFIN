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
