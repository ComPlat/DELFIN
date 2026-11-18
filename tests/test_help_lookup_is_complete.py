"""Every answerable thing has an entry: the help_lookup completeness test.

The lookup answers from CONTROL_FIELD_SPECS and _TABS_HELP; this test
pins that those two tables are complete, so the lookup (and therefore
every key/tab answer, the Help tab and the header search) never reports
a gap for something DELFIN actually has.
"""

from delfin.dashboard.help_lookup import lookup


def test_every_control_key_has_a_help_text():
    from delfin.common.control_validator import CONTROL_FIELD_SPECS

    thin = [
        s.name for s in CONTROL_FIELD_SPECS
        if len((s.help or "").strip()) < 30
    ]
    assert not thin, f"{len(thin)} keys have no real explanation: {thin[:8]}"


def test_every_control_key_is_findable_and_returns_its_name():
    from delfin.common.control_validator import CONTROL_FIELD_SPECS

    missing = []
    for spec in CONTROL_FIELD_SPECS:
        r = lookup(spec.name)
        if not (r["found"] and any(
                e["kind"] == "key" and e["name"] == spec.name
                for e in r["entries"])):
            missing.append(spec.name)
    assert not missing, f"lookup cannot find: {missing[:8]}"


def test_every_control_key_states_what_it_accepts():
    from delfin.common.control_validator import CONTROL_FIELD_SPECS

    no_allowed = [
        s.name for s in CONTROL_FIELD_SPECS
        if not (getattr(s, "allowed", "") or "").strip()
    ]
    assert not no_allowed, (
        f"{len(no_allowed)} keys never state their allowed values; "
        f"add allowed= to each: {no_allowed[:10]}"
    )


def test_every_dashboard_tab_has_a_tabs_help_entry():
    """The fixed tabs create_dashboard always builds must all be described."""
    from delfin.dashboard.tab_help import _TABS_HELP

    described = {t for t, *_rest in _TABS_HELP}
    # The names create_dashboard gives its fixed tabs (id -> title):
    required = {
        "Submit Job", "Recalc", "ORCA Builder", "TURBOMOLE Builder",
        "Job Status", "Calculations", "Archive", "Home", "Office",
        "Literature", "Remote Archive", "Tools", "Settings",
        "DELFIN Agent", "Agent Activity",
    }
    missing = required - described
    assert not missing, f"tabs without a guide entry: {sorted(missing)}"
    for _t, d, _w, _f, _i in _TABS_HELP:
        assert len(d.strip()) >= 25, "a tab entry is too thin"


def test_lookup_reports_no_entry_for_something_that_is_not_there():
    r = lookup("definitely_not_a_key_or_tab_xyz")
    assert not r["found"] and r["entries"] == []


def test_lookup_scores_exact_names_first():
    r = lookup("maxcore")
    assert r["entries"][0]["name"] == "maxcore"
    r = lookup("Submit Job")
    assert r["entries"][0]["kind"] == "tab"
    assert r["entries"][0]["name"] == "Submit Job"
