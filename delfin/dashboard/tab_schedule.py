"""Startup/lazy scheduling for dashboard tabs.

Holds the pure part of the lazy-build contract so it is testable without a
running dashboard:

* :func:`startup_specs_to_build` -- which specs are built before the page
  goes out: only the first visible one, so the first tab is interactive
  before any other tab's build cost is paid.
* :func:`open_spec` -- build one spec's widget the first time it is opened
  and return the same widget on every later open, so a built tab keeps its
  state.
"""


def first_visible_spec(specs):
    """The first spec in *specs* that is available and not fixed."""
    for spec in specs:
        if spec.get('available') and not spec.get('fixed'):
            return spec
    return None


def startup_specs_to_build(specs):
    """The specs whose widget must exist before the page goes out.

    Input: the ordered visible specs (as ``_sorted_visible_specs`` returns
    them). Output: the one spec the reader sees at startup, or [] when the
    bar is empty. Every other visible spec is built lazily on first open by
    :func:`open_spec`.
    """
    first = first_visible_spec(specs)
    return [] if first is None else [first]


def open_spec(spec):
    """Return ``spec['widget']``, building it once through ``spec['build']``.

    A spec that is already built is returned unchanged -- the same widget
    object, so nothing on the page is re-created and the tab keeps its
    state. A build that raises leaves the spec unbuilt and yields None, so
    a broken tab never breaks the bar; the next open retries the build.
    """
    if spec.get('widget') is not None:
        return spec['widget']
    build = spec.get('build')
    if not callable(build):
        return None
    try:
        spec['widget'] = build()
    except Exception:
        spec['widget'] = None
        return None
    return spec['widget']
