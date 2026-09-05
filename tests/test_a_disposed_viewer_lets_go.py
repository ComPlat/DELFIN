"""3Dmol's viewer constructor takes three listeners and never gives them back.

One ``window resize``, one ``body mouseup``, one ``body touchend``, each
``.bind()``-ed so no handle survives to remove them by, and the vendored bundle
contains no matching ``removeEventListener``. They accumulate for as long as
the tab is open.

Measured against the real bundle in a browser, one viewer:

    live listeners     window:resize, body:mouseup, body:touchend  -> +3
    after dispose      0 remain
    addEventListener   no own property left shadowing the inherited one

And what it was costing before: over 12 structure opens, handling a single
window resize event went from 0.14 ms to 2.44 ms and kept climbing. It only
bites while the window or the splitter is being dragged -- which is exactly
when a picture is supposed to follow the hand.

The shadow is the trap worth naming. ``addEventListener`` is inherited from
``EventTarget.prototype``; putting one back by assignment leaves a permanent
own property in front of it. The descriptor is saved and restored the same way
``__delfinWithPixelRatio`` handles ``devicePixelRatio``.
"""

from delfin.dashboard import molecule_viewer


PATCH = molecule_viewer.RIGHT_MOUSE_TRANSLATE_PATCH_JS


def test_what_the_constructor_adds_is_noted():
    noting = PATCH.split('window.__delfinNotingListeners')[1].split(
        'window.__delfinCreateViewer')[0]

    assert '[window, document.body]' in noting
    assert 'noted.push([t, type, fn, opts])' in noting
    # The original is still called, or the viewer loses its own events.
    assert 'original.call(t, type, fn, opts)' in noting


def test_the_patch_is_put_back_by_descriptor_and_not_by_assignment():
    """An assignment would leave an own property shadowing the inherited
    addEventListener for the life of the page."""
    noting = PATCH.split('window.__delfinNotingListeners')[1].split(
        'window.__delfinCreateViewer')[0]

    assert 'Object.getOwnPropertyDescriptor(t, "addEventListener")' in noting
    assert 'Object.defineProperty(t, "addEventListener", saved[i])' in noting
    assert 'else delete t.addEventListener' in noting
    # In a finally, so a constructor that throws does not leave it patched.
    assert '} finally {' in noting


def test_the_viewer_carries_its_own_list():
    create = PATCH.split('window.__delfinCreateViewer = function')[1].split(
        'var bootstrapAttempt')[0]

    assert '__delfinNotingListeners' in create
    assert 'viewer.__delfinNotedListeners = (caught.noted || [])' in create
    # The pixel-ratio wrapper still wraps the construction, not the other way
    # round: it has to be in place while 3Dmol reads devicePixelRatio.
    assert create.index('__delfinWithPixelRatio') < create.index(
        '__delfinNotingListeners')


def test_disposing_takes_them_off():
    dispose = PATCH.split('window.__delfinDisposeViewer = function')[1].split(
        'window.__delfinWithPixelRatio')[0]

    assert 'viewer.__delfinNotedListeners' in dispose
    assert 'removeEventListener(noted[n][1],noted[n][2],noted[n][3])' in dispose
    # A viewer built past the funnel simply has none, and must not throw.
    assert 'Array.isArray(noted)' in dispose
    assert 'viewer.__delfinNotedListeners=null' in dispose


def test_it_happens_before_the_context_is_thrown_away():
    """Removing a listener from a viewer whose WebGL context has already been
    lost still works, but the order says what was meant."""
    dispose = PATCH.split('window.__delfinDisposeViewer = function')[1].split(
        'window.__delfinWithPixelRatio')[0]

    assert dispose.index('removeEventListener') < dispose.index('loseContext')
