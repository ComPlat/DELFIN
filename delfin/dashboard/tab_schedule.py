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


#: Seconds after the page went out before the first idle build, and
#: between two. Each build holds the kernel for its own duration; a click
#: in that moment waits at most for one tab, never for the whole bar.
WARM_FIRST_S = 1.0
WARM_GAP_S = 0.4


def _call_later(delay, fn):
    """Run *fn* on the kernel's event loop after *delay* s; False when
    there is no running loop (a test, a plain script)."""
    try:
        import asyncio
        loop = asyncio.get_event_loop()
        if not loop.is_running():
            return False
        loop.call_later(delay, fn)
        return True
    except Exception:
        return False


def warm_up(ctx, build_now, order, *, call_later=None, env=None):
    """Build the tabs not yet built, one per idle moment, in *order*.

    Input: the dashboard context, its ``_build_tab_now`` and the specs in
    the order they should become ready. Each step builds the first spec of
    *order* that is still unbuilt (and buildable), then schedules the next
    step; a spec a click already built is skipped. Off with
    ``DELFIN_DASHBOARD_WARM=0``. Returns whether the first step was
    scheduled.
    """
    import os
    env = os.environ if env is None else env
    if str(env.get('DELFIN_DASHBOARD_WARM', '1')).strip() in ('0', 'false', 'no', 'off'):
        return False
    call_later = call_later or _call_later
    queue = []
    seen = set()
    for spec in order or []:
        if id(spec) in seen:
            continue
        seen.add(id(spec))
        queue.append(spec)

    def step():
        while queue:
            spec = queue.pop(0)
            if spec.get('widget') is None and callable(spec.get('build')):
                try:
                    build_now(spec)
                except Exception:
                    pass
                break
        if queue:
            call_later(WARM_GAP_S, step)

    return bool(call_later(WARM_FIRST_S, step))
