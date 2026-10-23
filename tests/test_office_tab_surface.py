"""The Office tab is the file browser without the chemistry surface.

Office mirrors the calculations browser on its own root, so it inherits
every control the browser has -- including the ones that report on
calculations: a keyword menu of ORCA section headings, the DELFIN report,
the 3D structure toggle, the move into the calculation archive, and the
running SSH transfers. In a folder of spreadsheets and letters none of
them has anything to act on.

They are retired in CSS rather than by setting display on the widgets,
because the per-file-type code enables and re-shows those buttons as
files are selected and would undo a flag set once at build time. These
tests pin the rule and, just as importantly, that Calculations still has
all of it.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard import tab_archive_statistics, tab_office
from delfin.dashboard.context import DashboardContext

# Controls that report on a calculation, not on a document.
_CHEMISTRY_CONTROLS = (
    'calc-search-suggest',      # ORCA section keywords
    'calc-report-btn',          # DELFIN report
    'calc-view-toggle',         # 3D structure
    'calc-move-archive-btn',    # into the calculation archive
    'calc-ssh-transfer-btn',    # send to the compute cluster
    'calc-transfer-jobs-btn',   # running SSH transfers
)


@pytest.fixture(autouse=True)
def every_optional_control_is_present(monkeypatch):
    """Build the full control surface, not the one this machine happens to have.

    calc-ssh-transfer-btn is only placed when the remote archive feature is
    switched on in the user settings, which is off by default. A developer
    who runs with it on sees the whole surface; CI, with default settings,
    sees one control fewer -- and the test that pins the hiding rule then
    fails for a reason that has nothing to do with Office. Turning the
    feature on here keeps both assertions strong: the rule must name a
    control that exists, and Calculations must keep all of them.
    """
    monkeypatch.setattr(browser, 'load_remote_archive_enabled', lambda: True)


@pytest.fixture
def ctx(tmp_path):
    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    return DashboardContext(
        calc_dir=tmp_path / 'calc',
        archive_dir=tmp_path / 'archive',
        office_dir=tmp_path / 'office',
    )


def _classes(widget):
    """Every CSS class in the widget tree, and the CSS it carries."""
    found: set[str] = set()
    css: list[str] = []

    def walk(node):
        for name in getattr(node, '_dom_classes', ()) or ():
            found.add(name)
        value = getattr(node, 'value', None)
        if isinstance(value, str) and '<style>' in value:
            css.append(value)
        for child in getattr(node, 'children', ()) or ():
            walk(child)

    walk(widget)
    return found, '\n'.join(css)


def test_office_marks_itself_and_calculations_does_not(ctx):
    office, _ = tab_office.create_tab(ctx)
    calc, _ = browser.create_tab(ctx)

    assert 'calc-office' in _classes(office)[0]
    assert 'calc-office' not in _classes(calc)[0]


def test_the_archive_mirror_is_not_mistaken_for_office(ctx):
    """Both mirrors clone the context; only the root tells them apart."""
    archive, _ = tab_archive_statistics.create_tab(ctx)
    assert 'calc-office' not in _classes(archive)[0]


def test_the_chemistry_controls_are_retired_in_office(ctx):
    office, _ = tab_office.create_tab(ctx)
    classes, css = _classes(office)

    for control in _CHEMISTRY_CONTROLS:
        assert control in classes, (
            f'{control} is not on any widget, so the rule that hides it in '
            'Office points at nothing')
        assert f'.calc-office .{control}' in css, (
            f'{control} is still shown in Office')


def test_calculations_keeps_every_one_of_them(ctx):
    """The rule is scoped to .calc-office. Calculations must be untouched."""
    calc, _ = browser.create_tab(ctx)
    classes, _css = _classes(calc)
    for control in _CHEMISTRY_CONTROLS:
        assert control in classes


def test_the_rule_hides_rather_than_merely_dimming(ctx):
    office, _ = tab_office.create_tab(ctx)
    _classes_, css = _classes(office)
    rule_start = css.index('.calc-office .calc-search-suggest')
    rule = css[rule_start:rule_start + 400]
    assert 'display:none !important' in rule, (
        'a disabled-looking button still reads as a feature that is broken '
        'rather than one that does not apply here')


def test_office_says_what_it_is(ctx):
    office, _ = tab_office.create_tab(ctx)
    calc, _ = browser.create_tab(ctx)

    def heading(widget):
        for child in widget.children:
            value = getattr(child, 'value', '')
            if isinstance(value, str) and value.startswith('<h3>'):
                return value
        raise AssertionError('no heading')

    assert 'Office' in heading(office)
    assert 'Calculations' not in heading(office)
    assert 'Calculations Browser' in heading(calc)


def test_the_root_is_the_only_thing_that_decides():
    """One definition, asked once. Two places asking the same question in
    different words is how the cycle inspector ended up in a mode that has
    no cycle to inspect."""
    source = pathlib.Path(browser.__file__).read_text(encoding='utf-8')
    assert source.count('_is_office_tab = ') == 1
    # The root comparison happens in that definition and nowhere else.
    assert source.count('ctx.office_dir') == 1, (
        'the Office root is consulted somewhere besides the definition of '
        '_is_office_tab')
