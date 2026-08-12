"""Opening DRAW threw the page down, seconds after it was opened.

The frame appears below the button, so nothing moves at first -- and then two
to three seconds later Ketcher finishes loading inside it, takes the focus, and
the browser scrolls the frame into view. Measured in a browser:

    before    DRAW button 538 px down the viewport, pane at 0
    +1.2 s    unchanged, the 560 px frame already in place
    +3 s      pane at 654, button at -116, off the top of the screen

So the first scroll within six seconds of opening is undone, once, and the
correction then takes itself away. After the fix, across +1 to +7 s the button
stays at 538 and the pane at 0, and a scroll of the user's own to 400 stands.

A clamp was tried first -- hold the position, let go on the user's next
gesture -- and it is worse than the jump it fixes: it also swallowed the
user's own scrolling, and whether their wheel released it depended on what the
pointer happened to be over. One correction cannot do that. At most one
deliberate scroll inside those six seconds is undone, and the next one stands.
"""

from delfin.dashboard import tab_submit


def _hold_js():
    source = open(tab_submit.__file__, encoding='utf-8').read()
    return source.split('_KETCHER_SCROLL_HOLD_JS = """')[1].split('"""')[0]


def test_it_corrects_once_and_then_gets_out_of_the_way():
    js = _hold_js()

    # The scroll handler releases before it corrects, so it cannot fire twice.
    body = js.split('function onScroll()')[1].split('}')[0]
    assert 'release();' in body
    assert body.index('release();') < body.index('pane.scrollTop = keep')


def test_it_is_not_a_clamp():
    """No gesture listeners on window: a clamp needs them to know when to let
    go, and getting that wrong locks the user out of scrolling."""
    js = _hold_js()

    assert 'window.addEventListener' not in js
    assert "'wheel'" not in js and 'pointerdown' not in js


def test_it_lets_go_by_itself_even_if_nothing_scrolls():
    js = _hold_js()

    assert 'setTimeout(release, 6000)' in js
    assert 'clearTimeout(timer)' in js


def test_a_second_opening_replaces_the_first():
    js = _hold_js()

    assert 'if (window.__delfinDrawScrollRelease) window.__delfinDrawScrollRelease();' in js
    assert 'window.__delfinDrawScrollRelease = release;' in js
    assert 'window.__delfinDrawScrollRelease = null;' in js


def test_it_finds_the_pane_that_actually_scrolls():
    """Not document.body: the tab lives in a Lumino box that scrolls itself."""
    js = _hold_js()

    assert '/(auto|scroll)/.test(how.overflowY)' in js
    assert 'pane.scrollHeight > pane.clientHeight + 4' in js
    assert 'document.scrollingElement' in js


def test_it_runs_when_the_editor_is_opened():
    source = open(tab_submit.__file__, encoding='utf-8').read()
    opener = source.split('def on_submit_draw_open')[1].split('\n    def ')[0]

    assert '_run_manip_js(_KETCHER_SCROLL_HOLD_JS)' in opener
    # Only on the way in. Folding it shut moves nothing.
    assert opener.index('if not submit_draw_open_btn.value:') < opener.index(
        '_run_manip_js(_KETCHER_SCROLL_HOLD_JS)')
