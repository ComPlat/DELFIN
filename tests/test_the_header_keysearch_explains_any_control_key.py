"""The header key search explains any CONTROL key without leaving the tab.

Data source is CONTROL_FIELD_SPECS (same as the Help tab and the
validator), the suggestions filter on name/help/allowed, and picking one
renders an explanation card with default and allowed values.  The widget
must degrade to absent, never break dashboard startup.
"""

import ipywidgets as W


def test_the_search_builds_and_lists_suggestions():
    from delfin.dashboard import header_keysearch

    wrapper = header_keysearch.build_key_search()
    assert wrapper is not None
    # _specs reads the real table; the reducer is the assertion here.
    entries = header_keysearch._specs()
    assert len(entries) > 100, 'the specs table shrank unexpectedly'
    names = [e['name'] for e in entries]
    assert 'maxcore' in names and 'reduction_steps' in names, (
        'the keys the user asked about must be searchable'
    )


def test_a_picked_suggestion_shows_an_explanation_card():
    from delfin.dashboard import header_keysearch

    wrapper = header_keysearch.build_key_search(max_suggestions=3)
    search, dropdown = wrapper.children[0].children
    card = wrapper.children[1]

    search.value = 'reduction'
    dropdown.children[0].click()  # first suggestion button
    html = card.value
    assert 'reduction_steps' in html
    assert 'default' in html.lower()
    assert 'Allowed' in html
    assert card.layout.display != 'none'


def test_no_matches_close_cleanly_and_empty_query_hides_everything():
    from delfin.dashboard import header_keysearch

    wrapper = header_keysearch.build_key_search()
    search, dropdown = wrapper.children[0].children
    card = wrapper.children[1]

    search.value = 'zzz_no_such_key'
    assert dropdown.children and 'no key matches' in dropdown.children[0].value
    search.value = ''
    assert dropdown.layout.display == 'none'
    assert card.layout.display == 'none'
