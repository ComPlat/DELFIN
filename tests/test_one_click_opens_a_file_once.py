"""One press, two reports, and a guard that never matched.

Clicking a file in the calculations browser is reported twice: by the browser's
click bridge and by the selection change the same press causes. Only one of
them is supposed to open the file; the other sees that it has just been opened
and stands down.

It never did. The two do not spell the label alike. The bridge keeps the
non-breaking space the list is drawn with and the selection change hands over
a plain one:

    bridge            '\U0001f52c\\xa0second.xyz'
    selection change  '\U0001f52c second.xyz'

so the comparison was False every time and both reports opened the file. Every
open builds a whole 3Dmol viewer. Measured in a browser on a single-frame xyz,
per click:

    viewers built     2-3      ->  1
    main thread       232-309 ms  ->  90-164 ms
    what survives     1 stage, 1 canvas -- unchanged

The key is the file name now, which is what the path is built from a few lines
above, so the two reports agree by construction. Clicking the same file again
still re-opens it: that goes through the bridge, which has no guard.
"""

import pathlib
import re

import pytest

from delfin.dashboard import tab_calculations_browser as browser


@pytest.fixture(scope='module')
def source():
    return pathlib.Path(browser.__file__).read_text(encoding='utf-8')


def test_the_two_spellings_really_differ(source):
    """The reason the guard failed, kept as a fact rather than a story."""
    from_bridge = '\U0001f52c\xa0second.xyz'
    from_selection = '\U0001f52c second.xyz'

    assert from_bridge != from_selection

    # And the name behind them is the same, which is why it is the key:
    # _calc_label_to_name flattens the non-breaking space before it splits.
    namer = source.split('def _calc_label_to_name')[1].split('\n    def ')[0]
    assert "replace('\\xa0', ' ')" in namer


def test_the_guard_is_keyed_on_the_name(source):
    opener = source.split('def _calc_open_file_label')[1].split('\n    def ')[0]
    assert "state['last_opened_label'] = (state['current_path'], name)" in opener
    assert "(state['current_path'], label)" not in opener, 'that is the label again'


def test_both_sides_compare_the_same_thing(source):
    selection = source.split('def calc_on_selection_change')[1].split('\n    def ')[0]
    assert '_calc_label_to_name(labels[0])' in selection
    assert "state.get('last_opened_label')" in selection


def test_the_bridge_still_opens_without_a_guard(source):
    """Clicking the same file again is a request to open it again. The bridge
    is what serves that, so it must not be guarded."""
    bridge = source.split('def calc_on_open_request')[1].split('\n    def ')[0]
    assert '_calc_open_file_label(label)' in bridge
    assert 'last_opened_label' not in bridge


def test_a_file_name_is_still_scoped_to_its_folder(source):
    """The same name in another directory has to open."""
    opener = source.split('def _calc_open_file_label')[1].split('\n    def ')[0]
    assert "state['current_path']" in opener
    stored = re.search(r"state\['last_opened_label'\] = \((.*?)\)", opener)
    assert stored and 'current_path' in stored.group(1)
