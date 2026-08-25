"""The picture of the path a relaxed scan just walked.

A scan reports itself in a sentence -- the highest point, where it ended, the
free energies when they were asked for -- and a sentence is the wrong shape for
the one thing a person reads off a profile at a glance: whether that is one
barrier or three, whether the walk ended in a minimum or was still climbing,
and how flat the top is.  Nothing here computes anything: every number on the
picture was already in hand when the verdict was written.

Drawn the way the rest of the dashboard draws -- a matplotlib figure on the Agg
backend, saved to a PNG and carried into an ``ipywidgets.HTML`` as a data URI,
which is what :func:`delfin.dashboard.helpers.render_neb_trajectory_plot_html`
does for a nudged elastic band's trajectory.  A profile is the same object as
that one, so it looks like it.

Two deliberate differences from that file, both measured:

* The figure is built through the object API (``Figure`` and
  ``FigureCanvasAgg``) rather than through ``pyplot``.  This is drawn from the
  scan's own background thread while the dashboard's other tabs may be drawing
  from theirs, and ``pyplot`` keeps one global figure registry that threads
  share; the object API touches no global state and needs no ``plt.close`` to
  avoid leaking figures.
* The points a scan computed are drawn *as points*, with only a thin line
  joining them in the order they were walked.  A scan has eight to forty points
  and a smooth curve through them would claim a resolution nothing measured --
  the difference matters exactly at the top, which is where a reader looks.

The marked points are the three the scan puts into the undo history (where it
started, the highest point it crossed, where it was left), so the picture names
the same three places Undo steps back through.
"""

from __future__ import annotations

import base64
import html
import io

#: The house colours, taken from the NEB trajectory plot so that two profiles
#: in one dashboard look like two of the same thing: the walk in blue, the
#: second series in orange, and the two landmarks in the red and green that
#: file already labels a peak and an end point with.
_WALK = '#1565c0'
_FREE = '#ef6c00'
_TOP = '#b71c1c'
_END = '#1b5e20'
_CEILING = '#78909c'
#: Text stays ink whatever it is about; the marker beside it carries the
#: identity.  Red-on-red and green-on-green read as two more series, and the
#: red/green pair is the one a colour-blind reader cannot separate -- which is
#: why the landmarks differ in shape as well, and are labelled.
_INK = '#263238'

#: How large the picture is made.  Wide and short, because it lives in a column
#: under the structure and a profile is read along its x axis; at 130 dpi the
#: PNG is 40 to 66 kB over the walks measured here, which is one write of one
#: widget at the end of a walk that took minutes.  The height carries a row
#: of legend under the axes: put inside them, the legend lands wherever
#: matplotlib finds a hole, and on a profile the hole is by the barrier.
_FIGSIZE = (7.2, 3.0)
_DPI = 130


def _bounds(values, headroom=0.16, floor=0.06):
    """The y range the data needs, with room over the top for its labels."""
    low, high = min(values), max(values)
    span = high - low
    if span <= 0:
        span = max(abs(high), 1.0) * 0.12
    return low - max(span * floor, 0.01), high + max(span * headroom, 0.02)


def figure(points, *, x_label, y_label, title, started=None, top=None,
           ended=None, ended_label='where it ended', free=(), ceiling=None,
           ceiling_label='', kept=None):
    """The profile as a matplotlib figure, for whoever wants to look at it.

    *points* are ``(coordinate, kcal/mol)`` in the order they were walked.
    *started*, *top* and *ended* are three of them, marked because they are
    the three the scan puts into the undo history: the picture then names the
    same places Undo steps back through, and the two the verdict names are
    among them.  *free* is ``[(coordinate, kcal/mol, label), ...]`` -- the
    three places a free energy was taken, which is all there ever are, drawn on
    the same axis as the electronic path because both are kcal/mol above where
    the walk started.  *ceiling* is a height on that axis, not a temperature:
    the caller has already added the thermal budget to the minimum the rise is
    measured from, because that is the comparison the verdict makes.  *kept* is
    the point handed back into the box when the walk ended somewhere the
    temperature cannot pay for.

    Returned rather than saved so that a test can read the axes instead of a
    PNG.
    """
    from matplotlib.figure import Figure

    fig = Figure(figsize=_FIGSIZE, dpi=_DPI)
    fig.set_facecolor('white')
    ax = fig.add_subplot(111)
    xs = [float(x) for x, _y in points]
    ys = [float(y) for _x, y in points]

    ax.plot(xs, ys, color=_WALK, linewidth=1.0, alpha=0.55, zorder=2)
    ax.scatter(xs, ys, s=26, color=_WALK, edgecolors='white', linewidths=0.6,
               zorder=3, label='what was computed')

    heights = list(ys)
    if free:
        heights += [float(y) for _x, y, _name in free]
        ax.plot([float(x) for x, _y, _n in free],
                [float(y) for _x, y, _n in free],
                color=_FREE, linewidth=1.0, linestyle='--', alpha=0.7,
                zorder=2)
        ax.scatter([float(x) for x, _y, _n in free],
                   [float(y) for _x, y, _n in free],
                   s=44, marker='D', color=_FREE, edgecolors='white',
                   linewidths=0.7, zorder=4, label='free energy')

    low, high = _bounds(heights)
    # The ceiling only where it is on the picture the data needs.  Forced into
    # the range, a 22 kcal/mol budget over a torsion whose whole profile is 4
    # kcal/mol flattens the profile into a line -- and the profile is what the
    # picture is for.  Where it does not fit, the sentence beside the picture
    # is what says the path is open.
    if ceiling is not None and low <= float(ceiling) <= high:
        ax.axhline(float(ceiling), color=_CEILING, linewidth=1.0,
                   linestyle=(0, (5, 3)), zorder=1,
                   label=ceiling_label or 'the temperature\'s ceiling')

    def _mark(point, colour, marker, name):
        """One landmark, with its number where the number stays readable.

        Above the marker unless a free energy sits there -- G is taken at
        exactly the three places that are marked, so its diamond shares the
        vertical, and at the top of a barrier it is usually the higher of the
        two.  Measured on the first Diels-Alder picture drawn here: the +19.2
        sat inside the orange diamond above it, and the -63.8 of a walk that
        ended in the corner hung off the bottom of the axes into the tick
        labels, which is why the side and the alignment are both chosen rather
        than fixed.
        """
        if not point:
            return
        x, y = float(point[0]), float(point[1])
        near = [float(g) for fx, g, _n in free
                if abs(float(fx) - x) < 1e-9] if free else []
        # Below only when something is in the way above it *and* there is
        # somewhere to go: a walk that ends at the bottom of its own axis has
        # no room underneath, and the number went out through the floor into
        # the tick labels.
        above = not (near and max(near) > y
                     and (y - low) > 0.12 * ((high - low) or 1.0))
        edge = (x - min(xs)) / ((max(xs) - min(xs)) or 1.0)
        if xs[0] > xs[-1]:
            edge = 1.0 - edge
        ax.scatter([x], [y], s=70, marker=marker, color=colour,
                   edgecolors='white', linewidths=0.8, zorder=5, label=name)
        ax.annotate(
            f'{y:+.1f}', xy=(x, y),
            xytext=(6 if edge < 0.08 else -6 if edge > 0.92 else 0,
                    9 if above else -14),
            textcoords='offset points', fontsize=8, color=_INK,
            ha='left' if edge < 0.08 else 'right' if edge > 0.92 else 'center',
            va='bottom' if above else 'top',
            bbox={'boxstyle': 'round,pad=0.18', 'facecolor': 'white',
                  'edgecolor': colour, 'alpha': 0.95, 'linewidth': 0.7})

    # And where it started, in the walk's own colour and without a number on
    # it: it is the zero the y axis is quoted against, so "+0.0" beside it
    # would be noise.  Marked at all because it is the third station Undo
    # steps back to, and the line under the picture says so.
    if started:
        ax.scatter([float(started[0])], [float(started[1])], s=54, marker='o',
                   facecolors='white', edgecolors=_WALK, linewidths=1.4,
                   zorder=5, label='where it started')

    # One marker where there is one point.  A walk that was interrupted, and
    # a walk that ran to its last step still climbing, ends *at* its highest
    # point -- measured on a Diels-Alder approach stopped after 8 of 20
    # points, where the triangle and the square were drawn on the same pixel
    # with their two labels on top of each other.  What that point is, is that
    # the walk was left there; the sentence beside the picture is what says
    # the top is an interruption rather than a barrier.
    same = (top and ended and abs(float(top[0]) - float(ended[0])) < 1e-9
            and abs(float(top[1]) - float(ended[1])) < 1e-9)
    if not same:
        _mark(top, _TOP, '^', 'the highest point')
    _mark(ended, _END, 's', ended_label)
    if kept:
        ax.scatter([float(kept[0])], [float(kept[1])], s=70, marker='o',
                   facecolors='none', edgecolors=_END, linewidths=1.4,
                   zorder=5, label='what the box was given')

    pad = 0.03 * ((max(xs) - min(xs)) or 1.0)
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    # Read in the order it was walked, whichever way the coordinate ran.
    #
    # A scan that closes a bond drives its coordinate downwards, and drawn on
    # an ascending axis the walk reads backwards: measured on the Diels-Alder
    # approach from 3.26 to 1.55 A, the product stood on the left and the
    # structure the walk started from on the right, against every reaction
    # profile ever drawn and against the sentence beside it.  So the axis
    # follows the walk: where it started is on the left, where it was left is
    # on the right, and the numbers descend if that is what the walk did.
    if xs[0] > xs[-1]:
        ax.invert_xaxis()
    ax.set_ylim(low, high)
    ax.set_xlabel(x_label, fontsize=9)
    ax.set_ylabel(y_label, fontsize=9)
    ax.set_title(title, fontsize=9.5, color=_INK)
    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.22, linestyle='--', linewidth=0.7)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    # A legend as soon as there is more than the walk itself on the picture:
    # what a mark means must never be carried by its colour alone.
    handles, labels = ax.get_legend_handles_labels()
    if len(labels) > 1:
        ax.legend(handles, labels, frameon=False, fontsize=7.5,
                  loc='upper center', bbox_to_anchor=(0.5, -0.24),
                  ncol=min(len(labels), 3), labelcolor=_INK,
                  handletextpad=0.4, columnspacing=1.4, borderaxespad=0.0)
    fig.tight_layout(pad=0.6)
    return fig


def profile_png(points, **drawn):
    """The same figure as PNG bytes."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    fig = figure(points, **drawn)
    FigureCanvasAgg(fig)
    out = io.BytesIO()
    fig.savefig(out, format='png', facecolor='white')
    return out.getvalue()


def profile_html(points, note='', **drawn):
    """The picture as a block of HTML, ready for an ``ipywidgets.HTML``.

    Sized in percent rather than in pixels: the column it sits in is the one
    the structure is drawn in, and that column is a different width in the
    tab, in the ORCA Builder and in the fullscreen overlay.

    *note* is one line under the picture, for what the picture cannot say
    itself -- which coordinates a multi-leg walk drove together, and that Undo
    steps back through the marked points.
    """
    encoded = base64.b64encode(profile_png(points, **drawn)).decode('ascii')
    said = (f"<div style='font-size:11px; color:#546e7a; margin-top:2px;'>"
            f"{html.escape(str(note))}</div>" if note else '')
    return (
        "<div style='width:100%; border:1px solid #d9dee3; border-radius:6px; "
        "background:#fff; padding:6px; box-sizing:border-box;'>"
        f"<img alt='{html.escape(str(drawn.get('title') or 'scan profile'))}' "
        f"src='data:image/png;base64,{encoded}' "
        "style='display:block; width:100%; height:auto; max-width:100%;' />"
        f"{said}</div>"
    )
