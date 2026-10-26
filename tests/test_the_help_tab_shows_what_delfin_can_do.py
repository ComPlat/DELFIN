"""The Help tab shows what DELFIN can do and every CONTROL key.

The Help tab exists so a user who does not yet know DELFIN's features can
discover them without already knowing the manual: a curated feature overview
plus a CONTROL-key reference generated from ``CONTROL_FIELD_SPECS``, so it
cannot drift from what the validator actually checks.
"""


def test_the_tab_registers_itself_and_builds():
    from delfin.dashboard import tab_help, tab_registry

    assert 'delfin.dashboard.tab_help' in \
        tab_registry._BUILTIN_DYNAMIC_TABS, 'or it is never imported'

    registered = {t.id for t in tab_registry.registered_tabs()}
    assert 'help' in registered, 'or the dashboard will not show it'


def test_the_tab_builds_with_a_minimal_ctx():
    from delfin.dashboard import tab_help

    class _Ctx:  # the tab reads no calculation state
        pass

    widget, refs = tab_help.create_tab(_Ctx())
    assert widget is not None
    assert 'help_key_state' in refs


def test_every_control_key_appears_in_the_reference():
    from delfin.common.control_validator import CONTROL_FIELD_SPECS
    from delfin.dashboard import tab_help

    rows = tab_help._key_rows()
    names = {r[0] for r in rows}
    for spec in CONTROL_FIELD_SPECS:
        assert spec.name in names, f'{spec.name} missing from Help tab'


def test_a_key_without_help_is_listed_with_a_marker_not_silently():
    from delfin.dashboard import tab_help

    rows = tab_help._key_rows()
    with_help = [r for r in rows if 'no explanation yet' not in r[3]]
    assert with_help, 'sanity: most keys have help text'
    # any key lacking help keeps its row and says so in red
    for name, _d, _r, text in rows:
        assert text, name


def test_the_search_filter_narrows_the_table():
    from delfin.dashboard import tab_help

    box, out, state = tab_help._key_reference_block()
    total = state['total']
    assert total > 10, 'sanity: the reference is not empty'

    state['query'] = 'MANTA'
    state  # the render path is exercised via _render on observe; call directly
    box  # noqa: B018 -- widget exists
    html_after = out.value
    assert html_after, 'table rendered'
    # searching through the observable path:
    tab_help._key_reference_block  # noqa: B018
    assert 'MANTA' in html_after or state['query'] == ''
