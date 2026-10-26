"""The Help tab shows what DELFIN can do and every CONTROL key.

The Help tab exists so a user who does not yet know DELFIN's features can
discover them without already knowing the manual: generalised functions as
an accordion of sub-functions, plus a CONTROL-key reference generated from
``CONTROL_FIELD_SPECS`` so it cannot drift from what the validator checks.
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


def test_feature_groups_have_subfunctions_and_keys():
    from delfin.dashboard import tab_help

    assert tab_help._FEATURE_GROUPS, 'or the overview is empty'
    for title, summary, subs in tab_help._FEATURE_GROUPS:
        assert title and summary
        assert subs, f'{title} has no sub-functions'
        for sub, what, keys in subs:
            assert sub and what, f'{sub} in {title} is incomplete'
            assert keys, f'{sub} in {title} names no keys'


def test_the_key_reference_covers_every_spec_key():
    from delfin.common.control_validator import CONTROL_FIELD_SPECS
    from delfin.dashboard import tab_help

    widget, refs = tab_help.create_tab(type('_Ctx', (), {})())
    html = refs['help_key_state'] and ''
    box, out, state = tab_help._key_reference_block()
    rendered = out.value
    for spec in CONTROL_FIELD_SPECS:
        assert spec.name in rendered, (
            f'{spec.name} is validated by CONTROL_FIELD_SPECS but missing '
            'from the Help tab key reference'
        )


def test_the_key_reference_groups_know_the_families():
    from delfin.dashboard import tab_help

    assert tab_help._group_of('MANTA_OPT') == 'MANTA'
    assert tab_help._group_of('GUPPY_RUNS') == 'GUPPY (legacy)'
    assert tab_help._group_of('OCCUPIER_compare') == 'OCCUPIER'
    assert tab_help._group_of('co2_coordination') == 'CO2'
    assert tab_help._group_of('functional') == 'General'


def test_every_spec_key_has_a_help_text_now():
    from delfin.common.control_validator import CONTROL_FIELD_SPECS

    missing = [s.name for s in CONTROL_FIELD_SPECS
               if not (s.help or '').strip()]
    assert not missing, f'no explanation for: {missing}'
