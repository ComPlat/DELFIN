"""The structure editor's toolbar, measured in a browser at four widths.

Reported by the user: "Knöpfe schauen teils raus aus Bildschirm" -- with a scan
armed, the last buttons of the internal-coordinate row (Path from here, Find
the path, Path to saddle) were past the right edge of the window with no way to
reach them.

Nothing about that is visible from the source, and jsdom does no layout, so it
is measured: the real Submit tab is built, its widget tree is written out as
the DOM ipywidgets gives the browser -- the element and the classes each view
builds, read out of the shipped ``@jupyter-widgets`` bundle rather than
guessed, and the Layout traits written on as inline style the way
``LayoutView.handleChange`` writes them -- and chromium is asked where every
control ended up.

What was measured before the fix, at a 1280 px window with a scan armed: ten
controls past the right edge of the toolbar and nine of them not the topmost
element at their own centre, embedded; eight and seven in the fullscreen
overlay.  No container between the toolbar and the document scrolls sideways
(they all hold their overflow with ``overflow-x: hidden``), so those controls
could not be reached by scrolling either.  The cause was one nested row: the
internal-coordinate group said ``flex-flow: row nowrap`` and ``flex: 0 0 auto``
inside a toolbar that wraps, and a flexbox breaks between its items and never
inside one.

The assertion is the acceptance test the user's complaint implies, and nothing
narrower: at every width, in both places, with the controls the method and the
hand leave on screen and again with every control a user can be shown, no
control's right edge is past the toolbar's and nothing is covered.

The three buttons that reached a transition state have since become one press
with two boxes beside it, and the checkbox that only revealed the scan's end
field has become the third answer in the box that asks where the walk goes.
Measured here again with a scan armed, before and after: 306/272 px at 1920
embedded, 406/372 at 1536, 472/436 at 1280, 604/568 at 1024, and 204/170,
236/204, 236/236, 304/272 in the fullscreen overlay -- shorter at seven of the
eight, the same at the eighth, and nothing past the edge in any of them.
"""

from __future__ import annotations

import html as _html
import pathlib
import sys
import tempfile

import pytest


WIDTHS = (1920, 1536, 1280, 1024)

#: The channels the page and the kernel talk through.  They live in the
#: toolbar because a widget outside the tree is outside the DOM; they are never
#: shown, so "every control a user can be shown" is measured without them.
CHANNELS = ('submit-pick-sync', 'submit-cmd-sync', 'submit-manip-sync',
            'submit-gfn-frame', 'submit-gfn-wall', 'submit-fs-row-break')

#: Every ``css_properties`` key of ipywidgets' LayoutModel: what a Layout can
#: actually put on an element.  Anything else a caller passes -- ``gap``,
#: ``flex_wrap``, ``box_sizing`` -- is dropped by traitlets and never reaches
#: the browser, so a harness that honoured it would be measuring a page the
#: user never sees.
_LAYOUT_TRAITS = (
    'align_content', 'align_items', 'align_self', 'border_top', 'border_right',
    'border_bottom', 'border_left', 'bottom', 'display', 'flex', 'flex_flow',
    'height', 'justify_content', 'justify_items', 'left', 'margin',
    'max_height', 'max_width', 'min_height', 'min_width', 'overflow', 'order',
    'padding', 'right', 'top', 'visibility', 'width', 'object_fit',
    'object_position', 'grid_auto_columns', 'grid_auto_flow', 'grid_auto_rows',
    'grid_gap', 'grid_template_rows', 'grid_template_columns',
    'grid_template_areas', 'grid_row', 'grid_column', 'grid_area',
)

_BUTTON_STYLE = {'primary': 'mod-primary', 'success': 'mod-success',
                 'info': 'mod-info', 'warning': 'mod-warning',
                 'danger': 'mod-danger'}


def _widget_stylesheet():
    """The real ipywidgets CSS, out of the extension that ships it.

    The extension inlines its stylesheets as template literals for webpack's
    style loader, so they are on disk but not as ``.css``.  Reading them is
    what makes this a measurement of the running dashboard rather than of a
    stand-in: the toolbar's only spacing is the 2 px margin that sheet gives
    every ``.jupyter-widgets`` element -- the ``gap`` the toolbar asks for is
    not a Layout trait and never arrives -- and a hand-written approximation
    would not have it.
    """
    import re

    static = (pathlib.Path(sys.prefix) / 'share' / 'jupyter' / 'labextensions'
              / '@jupyter-widgets' / 'jupyterlab-manager' / 'static')
    if not static.is_dir():
        pytest.skip('the @jupyter-widgets lab extension is not installed')
    blocks = []
    for path in sorted(static.glob('*.js')):
        source = path.read_text(errors='replace')
        for match in re.finditer(
                r"___CSS_LOADER_EXPORT___\.push\(\[module\.id, `", source):
            start = match.end()
            at = start
            while True:
                end = source.index('`', at)
                if source[end - 1] != '\\':
                    break
                at = end + 1
            block = source[start:end].replace('\\`', '`').replace('\\\\', '\\')
            blocks.append(re.sub(r"\$\{[^}]*\}", 'none', block))
    if not blocks:
        pytest.skip('the widget stylesheet is not in the shipped bundle')
    return '\n'.join(blocks)


def _style(widget):
    return '; '.join(
        f"{trait.replace('_', '-')}: {getattr(widget.layout, trait)}"
        for trait in _LAYOUT_TRAITS
        if getattr(widget.layout, trait, None) not in (None, ''))


def _attrs(widget, *classes):
    names = ' '.join(list(classes) + list(getattr(widget, '_dom_classes', ())))
    return (f'class="{_html.escape(names)}" '
            f'style="{_html.escape(_style(widget))}"')


def _button(widget, kind):
    classes = ['jupyter-widgets', 'jupyter-button', kind]
    style = _BUTTON_STYLE.get(getattr(widget, 'button_style', ''))
    if style:
        classes.append(style)
    if kind == 'widget-toggle-button' and getattr(widget, 'value', False):
        classes.append('mod-active')
    icon = (getattr(widget, 'icon', '') or '').strip()
    label = getattr(widget, 'description', '') or ''
    body = ''
    if icon:
        body += f'<i class="fa fa-{_html.escape(icon)}'
        body += '' if label else ' center'
        body += '"></i>'
    body += _html.escape(label) or ('' if icon else '&nbsp;')
    return f'<button {_attrs(widget, *classes)}>{body}</button>'


def _render_labelled(widget, classes, inner):
    """A control ipywidgets lays out as label-plus-body in an inline row."""
    label = getattr(widget, 'description', '') or ''
    tag = (f'<label class="widget-label">{_html.escape(label)}</label>'
           if label else
           '<label class="widget-label" style="display: none"></label>')
    names = ['jupyter-widgets', 'widget-inline-hbox'] + classes.split()
    return f'<div {_attrs(widget, *names)}>{tag}{inner}</div>'


def _render(widget):
    """The element (and children) this widget's view builds."""
    kind = type(widget).__name__
    if kind in ('Box', 'HBox', 'VBox', 'GridBox'):
        classes = ['jupyter-widgets', 'widget-container', 'widget-box']
        if kind == 'HBox':
            classes.append('widget-hbox')
        elif kind == 'VBox':
            classes.append('widget-vbox')
        elif kind == 'GridBox':
            classes = ['jupyter-widgets', 'widget-container', 'widget-gridbox']
        inner = ''.join(_render(child) for child in widget.children)
        return f'<div {_attrs(widget, *classes)}>{inner}</div>'
    if kind == 'Button':
        return _button(widget, 'widget-button')
    if kind == 'ToggleButton':
        return _button(widget, 'widget-toggle-button')
    if kind in ('Dropdown', 'Select', 'SelectMultiple'):
        options = ''.join(
            '<option>%s</option>'
            % _html.escape(str(o[0] if isinstance(o, tuple) else o))
            for o in (widget.options or ()))
        return _render_labelled(widget, 'widget-dropdown',
                                f'<select>{options}</select>')
    if kind in ('IntText', 'FloatText', 'BoundedIntText', 'BoundedFloatText'):
        return _render_labelled(widget, 'widget-text',
                                f'<input type="number" value="{widget.value}">')
    if kind in ('Text', 'Password', 'Combobox'):
        return _render_labelled(widget, 'widget-text',
                                '<input type="text" class="widget-input">')
    if kind == 'Textarea':
        return _render_labelled(
            widget, 'widget-textarea',
            '<textarea class="widget-input" rows="5"></textarea>')
    if kind == 'Checkbox':
        label = _html.escape(getattr(widget, 'description', '') or '')
        attrs = _attrs(widget, 'jupyter-widgets', 'widget-inline-hbox',
                       'widget-checkbox', 'widget-label-basic')
        return (f'<div {attrs}><label class="widget-label-basic">'
                f'<input type="checkbox">{label}</label></div>')
    if kind in ('IntSlider', 'FloatSlider', 'IntRangeSlider',
                'FloatRangeSlider', 'SelectionSlider', 'FloatLogSlider'):
        return _render_labelled(
            widget, 'widget-slider widget-hslider',
            '<div class="slider-container">'
            '<div class="slider noUi-target noUi-horizontal"></div></div>'
            f'<div class="widget-readout">{widget.value}</div>')
    if kind == 'HTML':
        return _render_labelled(
            widget, 'widget-html',
            f'<div class="widget-html-content">{widget.value}</div>')
    if kind == 'Label':
        return (f'<div {_attrs(widget, "jupyter-widgets", "widget-label")}>'
                f'{_html.escape(widget.value)}</div>')
    # An Output holds a picture that never arrives without a kernel, and every
    # widget not named above still occupies a box of its own.
    classes = ['jupyter-widgets']
    if kind == 'Output':
        classes += ['widget-output', 'jupyter-widgets-output-area']
    inner = ''.join(_render(c) for c in getattr(widget, 'children', ()))
    return f'<div {_attrs(widget, *classes)}>{inner}</div>'


def _build_tab():
    from delfin.dashboard.context import DashboardContext
    from delfin.dashboard import tab_submit

    folder = pathlib.Path(tempfile.mkdtemp())
    for name in ('calc', 'archive', 'office'):
        (folder / name).mkdir()
    ctx = DashboardContext(calc_dir=folder / 'calc',
                           archive_dir=folder / 'archive',
                           office_dir=folder / 'office')
    ctx.run_js = lambda _script: None
    return tab_submit.create_tab(ctx)


def _arm_a_scan(exports):
    """The state the complaint was made in: xtb, a pulling hand, a scan set up.

    Driven through the widgets the tab exports, so the tab's own observers --
    ``_refresh_method_controls`` and ``_refresh_hand_controls`` -- decide what
    is on screen, which is the thing the layout has to survive.
    """
    gfn = [value for _label, value in exports['submit_ff_dd'].options
           if str(value).lower().startswith('gfn2')]
    if gfn:
        exports['submit_ff_dd'].value = gfn[0]
    exports['submit_hand_dd'].value = 'pull'
    exports['submit_thermal_btn'].value = True
    for name in ('submit_scan_way', 'submit_scan_steps', 'submit_scan_to',
                 'submit_scan_dd', 'submit_scan_del',
                 'submit_scan_whole', 'submit_scan_how', 'submit_scan_energy',
                 'submit_scan_run_btn', 'submit_saddle_from',
                 'submit_saddle_how', 'submit_constraint_dd',
                 'submit_constraint_del'):
        exports[name].layout.display = ''


_MEASURE_JS = r"""
(sel) => {
  const bar = document.querySelector(sel);
  if (!bar) return {error: 'no toolbar'};
  const cs = getComputedStyle(bar);
  const box = bar.getBoundingClientRect();
  const right = box.right - parseFloat(cs.borderRightWidth)
                          - parseFloat(cs.paddingRight);
  const controls = [];
  bar.querySelectorAll('.jupyter-widgets').forEach((el) => {
    if (el.classList.contains('widget-box')) return;
    const s = getComputedStyle(el);
    if (s.display === 'none' || s.visibility === 'hidden') return;
    if (el.parentElement.closest('.jupyter-widgets:not(.widget-box)')) return;
    const r = el.getBoundingClientRect();
    controls.push({el: el,
                   name: (el.textContent || '').trim().slice(0, 20)
                         || el.className.slice(0, 40),
                   right: Math.round(r.right), width: Math.round(r.width),
                   past: r.right > right + 0.5});
  });
  // Reachable is asked of the browser rather than inferred. The page is
  // scrolled vertically to bring each control's centre to the middle of the
  // window -- which is all a user can do, because every container between the
  // toolbar and the document holds its overflow -- and the browser is then
  // asked what is at that point.
  const was = scrollY;
  controls.forEach((c) => {
    const first = c.el.getBoundingClientRect();
    scrollTo(0, scrollY + first.top + first.height / 2 - innerHeight / 2);
    const now = c.el.getBoundingClientRect();
    const x = now.left + now.width / 2, y = now.top + now.height / 2;
    const hit = (x >= 0 && x < innerWidth && y >= 0 && y < innerHeight)
        ? document.elementFromPoint(x, y) : null;
    c.covered = !(hit && (hit === c.el || c.el.contains(hit)
                          || hit.contains(c.el)));
    delete c.el;
  });
  scrollTo(0, was);
  return {height: Math.round(box.height),
          barRight: Math.round(right),
          scrolls: bar.scrollWidth > bar.clientWidth + 1,
          controls: controls};
}
"""

_SHOW_EVERYTHING_JS = r"""
([sel, channels]) => {
  document.querySelector(sel).querySelectorAll('.jupyter-widgets')
    .forEach((el) => {
      if (channels.some((c) => el.classList.contains(c))) return;
      if (el.style.display === 'none') el.style.removeProperty('display');
    });
}
"""

#: The manipulate toolbar, which is the first strip of either name in the
#: document -- it precedes the isomer stepper and the copy row in the tab, and
#: it is the first member the overlay takes.  Both names, so that a toolbar
#: that has lost the one the layout rules hang off is still found and still
#: measured, and fails on where its controls are rather than on a missing
#: class.
_TOOLBAR = '.delfin-structure-toolbar, .delfin-structure-fs-toolbar'


def test_no_control_in_the_toolbar_can_end_up_off_the_row():
    """Every control is inside the toolbar and clickable, at four widths.

    Both places the toolbar lives: embedded in the tab, where it shares the
    window with the form on the left, and in the body-level fullscreen overlay,
    which is as wide as the screen and outside the tab sheet that used to be
    the only thing keeping a control inside its column.

    And both extremes of what is on screen: what the chosen method and hand
    leave showing with a scan armed, which is what the user reported, and every
    control a user can be shown at once, which no single combination reaches
    but which no combination can exceed either.
    """
    pytest.importorskip('ipywidgets')
    playwright = pytest.importorskip(
        'playwright.sync_api', reason='needs a browser to lay the toolbar out')
    stylesheet = _widget_stylesheet()

    from delfin.dashboard.molecule_viewer import (
        structure_viewer_fullscreen_bootstrap_js,
        structure_viewer_fullscreen_kind_js,
    )

    tab_widget, exports = _build_tab()
    # The toolbar shows itself once there is a structure to act on; there is no
    # kernel here to load one.
    exports['submit_manip_toolbar'].layout.display = 'flex'
    _arm_a_scan(exports)
    document = (
        '<!doctype html><html><head><meta charset="utf-8"><style>'
        + stylesheet + '\nhtml, body { margin: 0; padding: 0; }\n'
        + '</style></head><body>' + _render(tab_widget)
        + '<script>' + structure_viewer_fullscreen_bootstrap_js() + '\n'
        + structure_viewer_fullscreen_kind_js(
            'submit', 'submit-scope-', ['_submitMolViewerByScope'])
        + '</script></body></html>')

    trouble = []
    with playwright.sync_playwright() as engine:
        try:
            browser = engine.chromium.launch()
        except Exception as exc:                     # no browser on this box
            pytest.skip(f'chromium unavailable: {exc}')
        try:
            for width in WIDTHS:
                page = browser.new_page(viewport={'width': width,
                                                  'height': 900})
                page.set_content(document)
                page.wait_for_timeout(120)
                for everything in (False, True):
                    if everything:
                        page.evaluate(_SHOW_EVERYTHING_JS,
                                      [_TOOLBAR, list(CHANNELS)])
                        page.wait_for_timeout(60)
                    for where in ('embedded', 'fullscreen'):
                        if where == 'fullscreen':
                            page.evaluate(
                                "() => document.querySelector("
                                "'.delfin-structure-fullscreen-btn').click()")
                            page.wait_for_timeout(120)
                        got = page.evaluate(_MEASURE_JS, _TOOLBAR)
                        where_said = (f'{width}px {where}'
                                      + (' everything' if everything else ''))
                        assert 'error' not in got, f'{where_said}: {got}'
                        assert got['controls'], f'{where_said}: no controls'
                        for control in got['controls']:
                            if control['past'] or control['covered']:
                                trouble.append(
                                    f"{where_said}: {control['name']!r} ends "
                                    f"at {control['right']}, the toolbar at "
                                    f"{got['barRight']}"
                                    + (' (covered)' if control['covered']
                                       else ''))
                        if got['scrolls']:
                            trouble.append(
                                f'{where_said}: the toolbar scrolls sideways')
                        if where == 'fullscreen':
                            page.keyboard.press('Escape')
                            page.wait_for_timeout(120)
                page.close()
        finally:
            browser.close()

    assert not trouble, 'controls the user cannot reach:\n' + '\n'.join(trouble)


def test_every_row_inside_the_toolbar_can_wrap_and_can_give_way():
    """The property the measurement above rests on, checked without a browser.

    A flexbox breaks between its items and never inside one, so a nested row
    that says ``nowrap`` puts its whole content on a single line however narrow
    the toolbar is -- which is exactly how the internal-coordinate row came to
    hang off the right of the screen.  A row that additionally says
    ``flex: 0 0 auto`` cannot even be shrunk back onto the line.

    This is the widget state, not the source text: it is asked of the Layout
    objects the tab hands the browser, so it stays true however the toolbar is
    assembled.  It runs where the browser test cannot, which is every machine
    without chromium.
    """
    pytest.importorskip('ipywidgets')
    import ipywidgets as widgets

    _tab_widget, exports = _build_tab()
    toolbar = exports['submit_manip_toolbar']

    def rows(box):
        for child in box.children:
            if isinstance(child, widgets.Box):
                yield child
                yield from rows(child)

    for row in [toolbar] + list(rows(toolbar)):
        if not row.children:
            continue                        # the overlay's line break
        flow = (row.layout.flex_flow or 'row').split()
        assert 'wrap' in flow, (
            f'a row of {len(row.children)} controls cannot wrap: '
            f'flex_flow={row.layout.flex_flow!r}')
        shrink = ((row.layout.flex or '0 1 auto').split() + ['1'])[1]
        assert shrink != '0', (
            f'a row of {len(row.children)} controls cannot be shrunk onto '
            f'its line: flex={row.layout.flex!r}')
