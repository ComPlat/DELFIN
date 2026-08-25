"""The energy profile of one route through a reaction network.

A step diagram, not a curve.  A scan walks a coordinate and its picture is a
line through the points it computed; a mechanism has no coordinate -- between
an intermediate and the saddle after it there is nothing anybody calculated,
and drawing a smooth line there would be drawing a path that was never walked.
So every species is a flat bar at its own energy and the bars are joined by
dashes, which is the convention every paper in the field uses and means
exactly "these two are connected and what is between them is not to scale".

Three things go on the picture that a table of numbers cannot say at a glance:

* **Where the highest point of the whole route is**, which is what decides the
  rate -- and it is frequently not the tallest single barrier, because a
  barrier is measured from the state before it while the route is measured
  from where it started.
* **The thermal ceiling** as a horizontal line.  Above it the route does not go
  at that temperature within that window; below it, it does.  It is the same
  instrument the editor holds a drag against, drawn across a whole mechanism.
* **What is missing.**  A point of the route with no record at this level is
  not drawn and is named under the picture instead.  A profile with a gap
  silently interpolated is a picture of a reaction nobody computed.

Drawn like :mod:`delfin.dashboard.scan_profile` -- a matplotlib figure on the
Agg backend through the object API rather than ``pyplot``, saved to PNG and
carried in as a data URI.  The object API touches no global state, which
matters because this is drawn while other tabs may be drawing too.
"""

from __future__ import annotations

import base64
import html
import io
from typing import Any, Dict, List, Optional, Sequence

#: Wide rather than tall.  A profile is read left to right and a mechanism of
#: six species needs the width; the height only has to separate the levels.
_FIGSIZE = (9.0, 4.6)
_DPI = 110

#: How wide a species' bar is, as a fraction of the space between two of them.
#: Wide enough to read a label on, narrow enough that the dashes between two
#: bars are visibly longer than the bars themselves -- which is the whole
#: visual claim: the flat part is computed, the sloped part is not.
_BAR = 0.52

_STATE_COLOUR = '#1f4e79'
_SADDLE_COLOUR = '#b7472a'
_JOIN_COLOUR = '#90a4ae'
_CEILING_COLOUR = '#2e7d32'


def _drawable(points: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    return [one for one in points if one.get('kcal') is not None]


def figure(points: Sequence[Dict[str, Any]], *, level: str,
           title: str = '', ceiling: Optional[float] = None,
           ceiling_label: str = '', missing: Sequence[str] = ()):
    """The profile as a matplotlib figure.

    *points* are what :func:`reaction_graph.profile` returns: dicts with
    ``ref``, ``label``, ``kind`` and ``kcal``, in reading order, already
    measured against the first state of the route.

    Returned rather than saved, so a test can read the axes instead of a PNG.
    """
    from matplotlib.figure import Figure

    fig = Figure(figsize=_FIGSIZE, dpi=_DPI)
    fig.set_facecolor('white')
    ax = fig.add_subplot(111)

    shown = _drawable(points)
    values = [float(one['kcal']) for one in shown]
    if ceiling is not None:
        values.append(float(ceiling))
    low = min(values) if values else 0.0
    high = max(values) if values else 1.0
    span = max(high - low, 1.0)
    ax.set_ylim(low - 0.18 * span, high + 0.22 * span)
    ax.set_xlim(-0.6, max(len(points) - 1, 1) + 0.6)

    # The bars, and the dashes between the ones that are there.  A gap in the
    # route breaks the chain rather than being bridged: two species joined
    # across a saddle nobody computed would read as a step that has been
    # established.
    previous = None
    for n, one in enumerate(points):
        if one.get('kcal') is None:
            previous = None
            ax.text(n, low - 0.10 * span, '?', ha='center', va='center',
                    fontsize=13, color='#9e9e9e')
            continue
        y = float(one['kcal'])
        colour = _SADDLE_COLOUR if one['kind'] == 'saddle' else _STATE_COLOUR
        ax.plot([n - _BAR / 2, n + _BAR / 2], [y, y],
                color=colour, linewidth=3.0, solid_capstyle='butt', zorder=3)
        if previous is not None:
            ax.plot([previous[0] + _BAR / 2, n - _BAR / 2],
                    [previous[1], y], color=_JOIN_COLOUR, linewidth=1.0,
                    linestyle='--', zorder=2)
        previous = (n, y)

    for n, one in enumerate(points):
        if one.get('kcal') is None:
            continue
        y = float(one['kcal'])
        above = one['kind'] == 'saddle'
        ax.annotate(f"{one['label']}\n{y:+.1f}", (n, y),
                    textcoords='offset points',
                    xytext=(0, 9 if above else -22), ha='center',
                    fontsize=8.5,
                    color=_SADDLE_COLOUR if above else _STATE_COLOUR)

    if ceiling is not None:
        ax.axhline(float(ceiling), color=_CEILING_COLOUR, linewidth=1.1,
                   linestyle=':', zorder=1)
        ax.text(ax.get_xlim()[1], float(ceiling), ' ' + (ceiling_label or ''),
                va='bottom', ha='right', fontsize=8, color=_CEILING_COLOUR)

    ax.set_ylabel('kcal/mol from the first state')
    ax.set_xticks([])
    ax.set_title(title or f'Reaction profile at {level}', fontsize=11)
    for edge in ('top', 'right', 'bottom'):
        ax.spines[edge].set_visible(False)
    ax.grid(axis='y', color='#eceff1', linewidth=0.8)
    ax.set_axisbelow(True)
    fig.tight_layout()
    return fig


def profile_png(points: Sequence[Dict[str, Any]], **drawn) -> bytes:
    """The same figure as PNG bytes."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    fig = figure(points, **drawn)
    FigureCanvasAgg(fig)
    out = io.BytesIO()
    fig.savefig(out, format='png', facecolor='white')
    return out.getvalue()


def summit(points: Sequence[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """The highest point of the whole route, which is what decides the rate.

    Frequently not the tallest single barrier: a barrier is measured from the
    state before it, and the route is measured from where it started.  A
    mechanism whose second step has the bigger barrier and whose first
    intermediate sits far uphill is governed by the second saddle *above the
    educt*, and quoting the bigger barrier there is the commonest way of
    getting a mechanism's rate wrong on paper.
    """
    shown = _drawable(points)
    if not shown:
        return None
    return max(shown, key=lambda one: float(one['kcal']))


def said_about(points: Sequence[Dict[str, Any]], *, missing: Sequence[str],
               ceiling: Optional[float], window: str) -> str:
    """One line under the picture, saying what the picture cannot."""
    said: List[str] = []
    top = summit(points)
    if top is not None:
        said.append(f"Highest point of the route: {top['label']} at "
                    f"{float(top['kcal']):+.1f} kcal/mol.")
        if ceiling is not None:
            said.append('The whole route is open within '
                        f'{window}.' if float(top['kcal']) <= float(ceiling)
                        else f'That is past what {window} pays for, so the '
                             'route is closed there.')
    if missing:
        said.append('Not drawn, because nothing is priced at this level: '
                    + ', '.join(str(one) for one in missing) + '.')
    return ' '.join(said)


def profile_html(points: Sequence[Dict[str, Any]], *, note: str = '',
                 **drawn) -> str:
    """The picture as a block of HTML, ready for an ``ipywidgets.HTML``."""
    encoded = base64.b64encode(profile_png(points, **drawn)).decode('ascii')
    said = (f"<div style='font-size:11px; color:#546e7a; margin-top:2px;'>"
            f"{html.escape(str(note))}</div>" if note else '')
    title = html.escape(str(drawn.get('title') or 'reaction profile'))
    return (
        "<div style='width:100%; border:1px solid #d9dee3; border-radius:6px; "
        "background:#fff; padding:6px; box-sizing:border-box;'>"
        f"<img alt='{title}' src='data:image/png;base64,{encoded}' "
        "style='display:block; width:100%; height:auto; max-width:100%;' />"
        f"{said}</div>"
    )


__all__ = ['figure', 'profile_png', 'profile_html', 'summit', 'said_about']
