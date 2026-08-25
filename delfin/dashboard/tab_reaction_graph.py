"""The Reactions tab: a network you can read, before it is a network you can look at.

The document underneath is :mod:`delfin.dashboard.reaction_graph`.  This is the
first thing put on top of it, and it is deliberately an indented list rather
than a drawing.  A drawing is what the graph is *for* and it is also the part
that hides mistakes: a layout algorithm that overlaps two nodes, an edge drawn
to the wrong end, a barrier printed beside the wrong arrow -- all of them look
like a picture either way.  A list is checkable by reading it, so everything
underneath can be got right while it is still cheap to be wrong.

What is on the screen
---------------------

Three questions across the top, because they are the three that change what
every number below means:

* **Which graph.**  A folder, chosen from the calculation directory.
* **At which level.**  Filled from the document, because a network can only be
  read at a level something in it actually has -- and because the alternative,
  a list of every method DELFIN can run, invites drawing a network at a level
  nothing in it has been computed at.
* **At which temperature, over how long.**  The same two halves of "possible"
  the editor holds a drag against, applied here to a whole network.

Then the network on the left and whatever is selected on the right.  Selection
is a :class:`ipywidgets.Select` rather than something clicked in the drawing,
and that is not a placeholder: it is keyboard-navigable, it survives a rebuild
of the list, and it works in a headless test -- which is the only way this
gets tested at all, because there is no browser on the machine it is built on.

What is *not* on the screen is the part worth stating.  A button appears when
it does something and not before: opening a structure in the editor, sending
one to a builder and drawing a profile are the next three steps, so their
buttons are not here yet.  A row of controls that raise "not implemented" is a
worse lie than an empty toolbar.
"""

from __future__ import annotations

import html
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from delfin.dashboard import reaction_graph as rg
from delfin.dashboard import reaction_profile
from delfin.dashboard.thermal import thermal_ceiling
from delfin.dashboard import run_results

try:
    import ipywidgets as widgets
    HAS_WIDGETS = True
except ImportError:                                  # pragma: no cover
    widgets = None
    HAS_WIDGETS = False


#: How long a barrier is given, as the windows a chemist actually means.  The
#: editor settled on one hour and made it a constant rather than a control --
#: between a second and a year the ceiling moves ten kcal/mol, while the gap
#: between chemistry and nonsense is twenty against a hundred.  Here it is a
#: control again, and the reason is that the question has changed: a drag is a
#: gesture and the answer has to arrive before the hand stops, while a network
#: is read deliberately, and "does this route open if I leave it overnight" is
#: exactly the question somebody brings to a mechanism.
WINDOWS: Tuple[Tuple[str, float], ...] = (
    ('a second', 1.0),
    ('a minute', 60.0),
    ('an hour', 3600.0),
    ('a day', 86400.0),
    ('a week', 604800.0),
    ('a year', 31557600.0),
)

#: Where a drawn profile is kept, inside the graph beside the calculations
#: it was made from.  A figure a person can only look at is one they will
#: redraw by hand for the paper.
FIGURES = 'figures'

#: Enough of a geometry to recognise it, and not enough to fill the panel.
_PREVIEW_ROWS = 12

#: What the editor hands over, as a plain dict.
#:
#: The editor is a workbench and knows nothing about networks; it knows what is
#: on screen, what made it, and what the last press left lying about.  So the
#: hand-over is one-way and untyped: the editor describes what it has, and this
#: side decides what -- if anything -- can be done with it.  Anything else puts
#: the graph's rules inside a 15,000-line module that has no business holding
#: them, and puts a second copy of them there to drift.
#:
#:     {'xyz':          the structure on screen,
#:      'level':        what a caption would call the method, e.g. 'GFN2-xTB',
#:      'method':       the short machine word, 'gfn2',
#:      'charge':       int,        'multiplicity': int,
#:      'energy':       Hartree or None,
#:      'free_energy':  Hartree or None,
#:      'imaginary':    int or None,   'frequency': cm-1 or None,
#:      'ends':         (xyz, xyz) or None -- the two a walk left,
#:      'gesture':      what produced it, for the history}
OFFER_KEYS = ('xyz', 'level', 'method', 'charge', 'multiplicity', 'energy',
              'free_energy', 'imaginary', 'frequency', 'ends', 'gesture')


def graphs_in(folder: Any) -> List[Path]:
    """Every reaction graph directly inside *folder*, by name.

    One level down and no deeper.  A graph carries a ``runs/`` full of ordinary
    job folders, and a recursive search would offer those as though a
    calculation were a network -- so the rule is "a folder holding a document",
    looked for exactly where graphs are made.
    """
    root = Path(folder or '.')
    out: List[Path] = []
    try:
        entries = sorted(root.iterdir(), key=lambda p: p.name.lower())
    except OSError:
        return out
    for entry in entries:
        try:
            if entry.is_dir() and (entry / rg.DOCUMENT).is_file():
                out.append(entry)
        except OSError:
            continue
    return out


# ---------------------------------------------------------------------------
# Saying what a number is, in words that do not overstate it
# ---------------------------------------------------------------------------

#: How a geometry came about, said in one word.
#:
#: There are two ways a structure changes in DELFIN and they are not the
#: same kind of fact.  One is the editor: a hand on an atom, a drawing, an
#: optimisation run from the workbench -- interactive, quick, and as good
#: as the method in the box.  The other is a calculation: an input
#: written, a job submitted, hours or days on a cluster, an output parsed.
#:
#: The document has always recorded which one (Record.source['kind']).
#: What it did not do was *show* it, and a geometry somebody dragged into
#: shape read exactly like one a two-day optimisation produced.  They are
#: not interchangeable: the first is where a calculation starts, the
#: second is what a paper quotes.
_CAME_FROM = {
    'editor': 'by hand',
    'scan': 'from a scan',
    'run': 'computed',
    'job': 'computed',
}


def said_origin(record: Optional[rg.Record]) -> str:
    """Which of the two kinds of change made this geometry."""
    if record is None:
        return ''
    return _CAME_FROM.get(str((record.source or {}).get('kind') or ''), '')


def said_energy(record: Optional[rg.Record]) -> str:
    """One record as a line, with the quantity it actually is.

    ``G`` and ``E`` are printed under their own names rather than both as
    "energy", because which one a diagram was drawn from is the difference
    between a free-energy profile and an electronic one, and nobody can tell
    them apart from the picture afterwards.

    And what nobody has asked is said as that.  ``imaginary`` is ``None``
    until a press takes the frequencies, which is a third state and not
    the same as zero -- a structure nobody has checked is not a minimum,
    it is a structure.
    """
    if record is None:
        return '--'
    parts: List[str] = []
    if record.free_energy is not None:
        parts.append(f'G {record.free_energy:.6f}')
    if record.energy is not None:
        parts.append(f'E {record.energy:.6f}')
    if not parts:
        parts.append('no energy')
    if record.imaginary == 0:
        parts.append('minimum')
    elif record.imaginary == 1:
        parts.append('saddle'
                     + (f' {record.frequency:.0f}i' if record.frequency else ''))
    elif record.imaginary is not None and record.imaginary > 1:
        parts.append(f'{record.imaginary} imaginary')
    else:
        parts.append('not checked')
    came = said_origin(record)
    if came:
        parts.append(came)
    if record.conformer:
        parts.append(f'conformer {record.conformer}')
    return ', '.join(parts)


def said_gate(verdict: Optional[bool], rise: Optional[float],
              window: str, *, why: str = '') -> str:
    """Whether the temperature pays for that barrier, in three answers.

    The third is the one that has to be said rather than implied.  A network
    drawn with the uncomputed edges greyed out like the closed ones reports a
    conclusion about calculations nobody has run, and the reader has no way to
    tell which grey is which.

    And the third answer says *which* thing is missing.  "Not at this
    level" sent a reader off to compute the saddle when what was unpriced
    was the state the step leaves -- an answer that is true and points the
    wrong way.
    """
    if verdict is None or rise is None:
        return why or 'nothing priced here'
    if verdict:
        return f'open within {window}'
    return f'closed within {window}'


def why_no_barrier(graph: rg.Graph, edge: rg.Edge, level: str) -> str:
    """Which end of a barrier is missing, named so it can be gone and got."""
    gaps = []
    if rg.priced(edge, level) is None:
        gaps.append('the saddle')
    foot = graph.node(edge.source)
    if rg.priced(foot, level) is None:
        gaps.append((foot.label or edge.source) if foot else edge.source)
    if not gaps:
        return 'nothing priced here'
    return 'no energy yet for ' + ' and '.join(gaps)


def network_lines(graph: rg.Graph, level: str, *,
                  temperature: Optional[float] = None,
                  seconds: Optional[float] = None,
                  window: str = 'an hour') -> List[Tuple[str, str]]:
    """The network as ``(id, line)`` in reading order, states then their steps.

    Every state, each followed by the transitions that leave it.  A state that
    several routes converge on appears once and is pointed at more than once,
    which is the thing a tree could not do and the reason the document is a
    graph.

    Returned as pairs rather than as text so that the selector and the display
    are built from one list: two lists in the same order is a thing that stays
    in the same order right up until it does not.
    """
    gates = rg.open_at(graph, level, temperature=temperature, seconds=seconds)
    out: List[Tuple[str, str]] = []
    for node in graph.nodes:
        record = rg.best(node, level)
        name = node.label or node.id
        charge = '' if node.charge == 0 else f' {node.charge:+d}'
        spin = '' if node.multiplicity == 1 else f' m={node.multiplicity}'
        out.append((node.id,
                    f'{node.id}  {name}{charge}{spin}   {said_energy(record)}'))
        for edge in graph.edges_from(node.id):
            rise = rg.barrier(graph, edge.id, level)
            target = graph.node(edge.target)
            to = (target.label or edge.target) if target else edge.target
            height = '   --  ' if rise is None else f'{rise:+7.1f}'
            mark = '' if edge.confirmed else '   (unconfirmed)'
            out.append((edge.id,
                        f'   {edge.id} -> {to}   {height} kcal/mol   '
                        f'{said_gate(gates.get(edge.id), rise, window, why=why_no_barrier(graph, edge, level))}{mark}'))
    return out


def _escape(text: Any) -> str:
    return html.escape(str(text or ''))


def _many(count: int, one: str, many: str = '') -> str:
    """``1 state`` and ``2 states``, because "1 transitions" reads as a bug.

    It is a small thing and it is the kind of small thing that tells a reader
    the numbers beside it were not looked at either.
    """
    return f'{count} {one if count == 1 else (many or one + "s")}'


def detail_html(graph: rg.Graph, ref: str, level: str) -> str:
    """Everything known about one state or transition, as the panel shows it.

    Every record, not only the one at the chosen level, because the question
    the panel answers is "what has been done to this" -- and a level that is
    missing here is the thing a user is deciding whether to start a calculation
    about.
    """
    holder = graph.holder(ref)
    if holder is None:
        return '<i>nothing selected</i>'
    rows: List[str] = []
    is_edge = isinstance(holder, rg.Edge)
    if is_edge:
        source = graph.node(holder.source)
        target = graph.node(holder.target)
        head = (f'{holder.id} &middot; {_escape(holder.label or "transition")}'
                f' &mdash; {_escape((source.label if source else holder.source))}'
                f' &rarr; {_escape((target.label if target else holder.target))}')
        if holder.confirmed:
            how = _escape(holder.confirmed.get('how') or 'confirmed')
            rows.append(f'<div>Confirmed to join these two ({how}).</div>')
        else:
            rows.append(
                '<div>Not confirmed to join these two. A saddle optimiser '
                'returns a saddle whatever it was given; only following the '
                'imaginary mode down both ways says which reaction it is the '
                'transition state of.</div>')
    else:
        head = f'{holder.id} &middot; {_escape(holder.label or "state")}'
        rows.append(f'<div>Charge {holder.charge}, multiplicity '
                    f'{holder.multiplicity}.</div>')
    if holder.note:
        rows.append(f'<div><i>{_escape(holder.note)}</i></div>')

    rows.append('<div style="margin-top:8px"><b>Records</b></div>')
    if not holder.records:
        rows.append('<div>None yet.</div>')
    for record in holder.records:
        here = ' &larr;' if record.level == level else ''
        rows.append(
            f'<div style="font-family:monospace">'
            f'{_escape(record.level)}{here} &nbsp; {_escape(said_energy(record))}'
            f'</div>')
        origin = record.source.get('kind') if record.source else ''
        where = record.source.get('run') if record.source else ''
        if origin:
            rows.append(
                f'<div style="margin-left:14px; opacity:.75">from {_escape(origin)}'
                + (f', {_escape(where)}' if where else '')
                + (f', {_escape(record.at)}' if record.at else '') + '</div>')
    if rg.best(holder, level) is None:
        rows.append(f'<div style="margin-top:6px">Nothing at '
                    f'<b>{_escape(level)}</b> yet.</div>')
    elif rg.priced(holder, level) is None:
        rows.append(
            f'<div style="margin-top:6px">There is a geometry at '
            f'<b>{_escape(level)}</b> and no energy for it, so nothing can '
            f'be measured against it yet. It is what a calculation would '
            f'start from.</div>')

    waiting = [p for p in graph.pending if p.on == ref]
    if waiting:
        rows.append('<div style="margin-top:8px"><b>Running</b></div>')
        for entry in waiting:
            rows.append(
                f'<div style="font-family:monospace">{_escape(entry.level)} '
                f'&nbsp; {_escape(entry.run)} &nbsp; since '
                f'{_escape(entry.submitted)}</div>')

    geometry = rg.geometry(graph, ref, level)
    if geometry:
        lines = geometry.splitlines()[:_PREVIEW_ROWS]
        if len(geometry.splitlines()) > _PREVIEW_ROWS:
            lines.append('...')
        rows.append('<div style="margin-top:8px"><b>Geometry at this level'
                    '</b></div>')
        rows.append('<pre style="margin:2px 0; font-size:12px">'
                    + _escape('\n'.join(lines)) + '</pre>')
    return (f'<div style="font-family:sans-serif; font-size:13px">'
            f'<div style="font-size:14px"><b>{head}</b></div>{"".join(rows)}</div>')


def summary_html(graph: rg.Graph, level: str, *, window: str = 'an hour',
                 temperature: Optional[float] = None) -> str:
    """One line about the whole network, and what it is not able to say.

    The gaps are named here rather than left to be noticed.  An energy diagram
    is a selection out of a document that may be complete at one level and not
    at another, and the sentence that makes it honest is "these four are not at
    this level yet" -- said where the level is chosen, not in a footnote.
    """
    missing = rg.missing_at(graph, level)
    gates = rg.open_at(graph, level, temperature=temperature)
    opened = sum(1 for one in gates.values() if one is True)
    closed = sum(1 for one in gates.values() if one is False)
    unknown = sum(1 for one in gates.values() if one is None)
    T = float(graph.temperature if temperature is None else temperature)
    said = [f'{_many(len(graph.nodes), "state")}, '
            f'{_many(len(graph.edges), "transition")}.']
    if graph.edges:
        said.append(f'At {T:g} K within {window}: {opened} open, '
                    f'{closed} closed'
                    + (f', {unknown} not priced yet.' if unknown else '.'))
    if missing:
        shown = ', '.join(missing[:8]) + ('...' if len(missing) > 8 else '')
        said.append(f'No energy at {level} for: {shown}')
    if graph.pending:
        said.append(f'{_many(len(graph.pending), "calculation")} running.')
    return ('<div style="font-family:sans-serif; font-size:13px">'
            + ' '.join(_escape(one) for one in said) + '</div>')


# ---------------------------------------------------------------------------
# The tab
# ---------------------------------------------------------------------------

class ReactionGraphPanel:
    """One Reactions tab, over the graphs in a calculation directory."""

    def __init__(self, calc_dir: Any = None, ctx: Any = None):
        self.calc_dir = Path(calc_dir) if calc_dir else Path.home() / 'calc'
        self.ctx = ctx
        self.graph: Optional[rg.Graph] = None
        self._building = False
        if not HAS_WIDGETS:                          # pragma: no cover
            self.widget = None
            return

        self.graph_dd = widgets.Dropdown(options=[], description='Graph:',
                                         layout=widgets.Layout(width='320px'))
        self.new_name = widgets.Text(placeholder='name for a new graph',
                                     layout=widgets.Layout(width='220px'))
        self.new_btn = widgets.Button(description='New', button_style='info',
                                      layout=widgets.Layout(width='80px'))
        self.reload_btn = widgets.Button(description='Reload',
                                         layout=widgets.Layout(width='90px'))
        #: Look at every calculation this graph is waiting for.  Pressed
        #: by hand as well as run on opening, because a person who has
        #: just watched a job leave the queue wants to ask now rather
        #: than close the tab and open it again.
        self.harvest_btn = widgets.Button(
            description='Check calculations',
            layout=widgets.Layout(width='auto', display='none'))
        self.harvest_btn.on_click(lambda _b: self._on_harvest())

        self.level_dd = widgets.Dropdown(options=[], description='Level:',
                                         layout=widgets.Layout(width='300px'))
        self.temperature = widgets.BoundedFloatText(
            value=298.15, min=1.0, max=5000.0, step=5.0, description='T (K):',
            layout=widgets.Layout(width='180px'))
        self.window_dd = widgets.Dropdown(
            options=[(name, name) for name, _ in WINDOWS], value='an hour',
            description='within:', layout=widgets.Layout(width='200px'))

        self.network = widgets.Select(options=[], rows=18,
                                      layout=widgets.Layout(width='100%'))
        self.detail = widgets.HTML(value='<i>nothing selected</i>')
        self.summary = widgets.HTML(value='')
        self.status = widgets.HTML(value='')
        #: The profile of one route, and the boxes that choose it.  Two
        #: ends and then which way between them, because a network with a
        #: branch in it has several and the one a person means is the one
        #: they are arguing about.
        self.from_dd = widgets.Dropdown(options=[], description='From:',
                                        layout=widgets.Layout(width='260px'))
        self.to_dd = widgets.Dropdown(options=[], description='to:',
                                      layout=widgets.Layout(width='240px'))
        self.route_dd = widgets.Dropdown(
            options=[], description='via:',
            layout=widgets.Layout(width='420px', display='none'))
        self.draw_btn = widgets.Button(
            description='Draw the profile', button_style='info',
            layout=widgets.Layout(width='auto', display='none'))
        self.picture = widgets.HTML(value='')
        self.from_dd.observe(self._on_ends, names='value')
        self.to_dd.observe(self._on_ends, names='value')
        self.draw_btn.on_click(lambda _b: self._on_draw())

        self.label_box = widgets.Text(placeholder='name',
                                      layout=widgets.Layout(width='220px'))
        self.rename_btn = widgets.Button(description='Rename',
                                         layout=widgets.Layout(width='100px'))
        self.note_box = widgets.Text(placeholder='note',
                                     layout=widgets.Layout(width='320px'))
        self.note_btn = widgets.Button(description='Save note',
                                       layout=widgets.Layout(width='110px'))
        #: Back to the workbench.  Absent where there is no workbench to go to
        #: -- this tab is usable on its own, and a button that would reach a
        #: tab that is not there is the row claiming something untrue.
        self.to_editor_btn = widgets.Button(
            description='Open in the editor', button_style='info',
            tooltip='Put this geometry in the Submit tab and go there',
            layout=widgets.Layout(width='auto', display='none'))
        self.to_editor_btn.on_click(self._on_to_editor)
        #: And on to a calculation.  The Submit tab is where a job is set
        #: up -- every field it has, every workflow it knows -- and this
        #: does not reproduce any of that.  It carries the geometry and a
        #: name across and notes in the graph that it went, so what comes
        #: back has somewhere to come back to.
        self.to_submit_btn = widgets.Button(
            description='Set up a calculation',
            tooltip='Take this geometry to the Submit tab and note it here',
            layout=widgets.Layout(width='auto', display='none'))
        self.to_submit_btn.on_click(self._on_to_submit)

        self.graph_dd.observe(self._on_graph, names='value')
        self.level_dd.observe(self._on_view, names='value')
        self.temperature.observe(self._on_view, names='value')
        self.window_dd.observe(self._on_view, names='value')
        self.network.observe(self._on_select, names='value')
        self.new_btn.on_click(self._on_new)
        self.reload_btn.on_click(lambda _b: self.refresh_graphs())
        self.rename_btn.on_click(self._on_rename)
        self.note_btn.on_click(self._on_note)

        self.widget = widgets.VBox([
            widgets.HBox([self.graph_dd, self.new_name, self.new_btn,
                          self.reload_btn, self.harvest_btn],
                         layout=widgets.Layout(flex_flow='row wrap')),
            widgets.HBox([self.level_dd, self.temperature, self.window_dd],
                         layout=widgets.Layout(flex_flow='row wrap')),
            self.summary,
            widgets.HBox([self.from_dd, self.to_dd, self.route_dd,
                          self.draw_btn],
                         layout=widgets.Layout(flex_flow='row wrap')),
            self.picture,
            widgets.HBox([
                widgets.VBox([self.network],
                             layout=widgets.Layout(width='52%')),
                widgets.VBox([self.detail,
                              widgets.HBox([self.label_box, self.rename_btn],
                                           layout=widgets.Layout(
                                               flex_flow='row wrap')),
                              widgets.HBox([self.note_box, self.note_btn],
                                           layout=widgets.Layout(
                                               flex_flow='row wrap')),
                              widgets.HBox([self.to_editor_btn,
                                            self.to_submit_btn],
                                           layout=widgets.Layout(
                                               flex_flow='row wrap'))],
                             layout=widgets.Layout(width='46%'))],
                layout=widgets.Layout(flex_flow='row wrap')),
            self.status,
        ])
        self.refresh_graphs()

    # -- what the boxes are showing --------------------------------------

    @property
    def window_seconds(self) -> float:
        want = str(self.window_dd.value)
        return next((s for name, s in WINDOWS if name == want), 3600.0)

    def _say(self, text: str) -> None:
        self.status.value = (
            f'<div style="font-family:sans-serif; font-size:12px">'
            f'{_escape(text)}</div>' if text else '')

    def refresh_graphs(self) -> None:
        """Look again at the calculation directory, keeping what was chosen."""
        found = graphs_in(self.calc_dir)
        want = self.graph_dd.value
        self._building = True
        try:
            self.graph_dd.options = [(p.name, str(p)) for p in found]
            if want and any(str(p) == want for p in found):
                self.graph_dd.value = want
            elif found:
                self.graph_dd.value = str(found[0])
        finally:
            self._building = False
        self._open(self.graph_dd.value)
        if not found:
            self._say(f'No reaction graphs in {self.calc_dir}. '
                      f'Name one and press New.')

    def _open(self, folder: Any) -> None:
        if not folder:
            self.graph = None
            self.network.options = []
            self.detail.value = '<i>nothing selected</i>'
            self.summary.value = ''
            return
        try:
            self.graph = rg.load(folder)
        except Exception as exc:                     # noqa: BLE001
            self.graph = None
            self._say(f'{Path(str(folder)).name} did not open: {exc}')
            return
        self.temperature.value = float(self.graph.temperature)
        # Opening the graph is when a person finds out what came back
        # while they were away, which is the whole point of a document
        # that outlives the session.  Quiet when nothing was waiting.
        #
        # Before the level box is filled, not after: what landed is very
        # often at a level the graph did not have until this moment -- that
        # is what the calculation was started for -- and filled first, the
        # box offered every level except the new one and the user had to
        # reload to see the answer they had been waiting two days for.
        landed = self.harvest(say=False)
        self._refresh_levels()
        self.refresh_view()
        if landed:
            self._say(' '.join(landed))

    def _refresh_levels(self) -> None:
        """The level box, filled from the document and nowhere else."""
        if self.graph is None:
            return
        levels = self.graph.levels()
        want = self.level_dd.value
        self._building = True
        try:
            self.level_dd.options = levels
            if want in levels:
                self.level_dd.value = want
            elif levels:
                self.level_dd.value = levels[0]
        finally:
            self._building = False

    def refresh_view(self) -> None:
        """Rebuild the list and the panel around whatever is selected."""
        if self.graph is None:
            return
        level = str(self.level_dd.value or '')
        window = str(self.window_dd.value)
        chosen = self.network.value
        lines = network_lines(self.graph, level,
                              temperature=float(self.temperature.value),
                              seconds=self.window_seconds, window=window)
        self._building = True
        try:
            self.network.options = [(line, ref) for ref, line in lines]
            if chosen and any(ref == chosen for ref, _ in lines):
                self.network.value = chosen
            elif lines:
                self.network.value = lines[0][0]
        finally:
            self._building = False
        self.harvest_btn.layout.display = ('' if self.graph.pending
                                           else 'none')
        self._refresh_ends()
        self.summary.value = summary_html(
            self.graph, level, window=window,
            temperature=float(self.temperature.value))
        self._refresh_detail()

    def _refresh_detail(self) -> None:
        if self.graph is None or not self.network.value:
            self.detail.value = '<i>nothing selected</i>'
            return
        ref = str(self.network.value)
        self.detail.value = detail_html(self.graph, ref,
                                        str(self.level_dd.value or ''))
        holder = self.graph.holder(ref)
        if holder is not None:
            self.label_box.value = getattr(holder, 'label', '') or ''
            self.note_box.value = getattr(holder, 'note', '') or ''
        # There has to be a geometry at the chosen level to send, and a Submit
        # tab to send it to.  Either missing and the press would be a promise
        # this tab cannot keep.
        can = bool(rg.geometry(self.graph, ref, str(self.level_dd.value or ''))
                   and self._submit_refs())
        self.to_editor_btn.layout.display = '' if can else 'none'
        self.to_submit_btn.layout.display = '' if can else 'none'

    # -- what the editor may hand over ------------------------------------

    def offer_label(self, offer: Dict[str, Any]) -> Optional[str]:
        """What pressing "Put in graph" would do, or None for nothing.

        The editor asks this to decide whether to show the button at all, and
        what to write on it.  Three answers and a silence:

        * A walk left **two ends**, so the press lays down two states and the
          step between them.  That is the shape a scan produces and it is the
          most valuable hand-over there is -- a barrier arrives with both the
          structures it is a barrier between.
        * The structure **matches exactly one state** already in the graph, so
          it is another geometry of that species at this level.
        * Otherwise it is a **new state**.

        The silence is a saddle with no two ends.  A transition is an edge, and
        an edge needs the two states it joins; a saddle put in as a state would
        be a minimum in the document that is not one.  There is nothing correct
        to do yet, so there is no button -- and standing right beside it is
        ``Down to both ends``, which appears exactly when there is a saddle and
        is what produces the ends.  The absence is the instruction.
        """
        if self.graph is None:
            return None
        xyz = str(offer.get('xyz') or '')
        if not rg.fingerprint(xyz):
            return None
        name = self.graph.name or Path(self.graph.folder).name
        ends = offer.get('ends')
        if ends and len(ends) == 2 and all(rg.fingerprint(one) for one in ends):
            return f'Put both ends and the step in {name}'
        if offer.get('imaginary') == 1:
            return None
        same = rg.looks_like(self.graph, xyz,
                             charge=offer.get('charge'),
                             multiplicity=offer.get('multiplicity'))
        level = str(offer.get('level') or 'this method')
        if len(same) == 1:
            return f'Add as a {level} record on {same[0].label or same[0].id}'
        return f'Put in {name} as a new state'

    def take(self, offer: Dict[str, Any]) -> str:
        """Do it, and hand back the sentence the editor puts on its status row.

        Everything that goes in carries where it came from, because a record
        that cannot say what produced it is a number in a document with no
        account behind it -- and the account is the whole reason the graph is
        on disk rather than in a notebook.
        """
        label = self.offer_label(offer)
        if label is None or self.graph is None:
            return ''
        graph = self.graph
        gesture = str(offer.get('gesture') or 'the editor')
        source = {'kind': 'editor', 'gesture': gesture}

        def _record(**over):
            fields = dict(
                level=str(offer.get('level') or 'unnamed method'),
                method=str(offer.get('method') or ''),
                charge=offer.get('charge'), multiplicity=offer.get('multiplicity'),
                energy=offer.get('energy'), free_energy=offer.get('free_energy'),
                imaginary=offer.get('imaginary'), frequency=offer.get('frequency'),
                source=source)
            fields.update(over)
            return rg.Record(**fields)

        ends = offer.get('ends')
        try:
            if ends and len(ends) == 2:
                said = self._take_a_step(offer, ends, _record)
            elif label.startswith('Add as a'):
                same = rg.looks_like(graph, str(offer['xyz']),
                                     charge=offer.get('charge'),
                                     multiplicity=offer.get('multiplicity'))
                rg.add_record(graph, same[0].id, str(offer['xyz']), _record())
                said = (f'{same[0].label or same[0].id} now has a '
                        f'{offer.get("level")} geometry as well.')
            else:
                node = rg.add_state(
                    graph, str(offer['xyz']), _record(),
                    label=str(offer.get('name') or
                              rg.formula(offer['xyz'])),
                    charge=int(offer.get('charge') or 0),
                    multiplicity=int(offer.get('multiplicity') or 1))
                said = f'{node.id} is in {graph.name}.'
            rg.save(graph)
        except Exception as exc:                     # noqa: BLE001
            return f'It did not go into the graph: {exc}'
        self._reopen()
        return said

    def _take_a_step(self, offer, ends, _record) -> str:
        """Two ends and the structure between them, as two states and an edge.

        The two ends go in as states only where they are not already there --
        a walk run twice from the same educt must not put that educt in twice,
        or the network grows a duplicate whose energies quietly disagree with
        its twin.
        """
        graph = self.graph
        made: List[str] = []

        def _end(xyz, what):
            same = rg.looks_like(graph, xyz, charge=offer.get('charge'),
                                 multiplicity=offer.get('multiplicity'))
            if len(same) == 1:
                return same[0].id
            node = rg.add_state(
                graph, xyz,
                _record(imaginary=None, frequency=None, energy=None,
                        free_energy=None,
                        source={'kind': 'editor',
                                'gesture': f'{offer.get("gesture")}, {what}'}),
                label=rg.formula(xyz), charge=int(offer.get('charge') or 0),
                multiplicity=int(offer.get('multiplicity') or 1))
            made.append(node.id)
            return node.id

        first = _end(ends[0], 'where it started')
        second = _end(ends[1], 'where it came to')
        if first == second:
            raise ValueError('both ends are the same structure, so there is '
                             'nothing between them')
        edge = rg.add_transition(graph, str(offer['xyz']), _record(),
                                 source=first, target=second,
                                 confirmed=offer.get('confirmed'))
        grew = (f' ({", ".join(made)} '
                f'{"is" if len(made) == 1 else "are"} new)') if made else ''
        return f'{edge.id} joins {first} to {second}{grew}.'

    def _reopen(self) -> None:
        """Read the document back and redraw, so the tab shows what landed."""
        if self.graph is None:
            return
        try:
            self.graph = rg.load(self.graph.folder)
        except Exception:                            # noqa: BLE001
            return
        self._refresh_levels()
        self.refresh_view()

    # -- what the presses do ---------------------------------------------

    def _on_graph(self, change) -> None:
        if not self._building:
            self._open(change.get('new'))

    def _on_view(self, _change) -> None:
        if not self._building:
            self.refresh_view()

    def _on_select(self, _change) -> None:
        if not self._building:
            self._refresh_detail()

    # -- the profile of one route -----------------------------------------

    def _refresh_ends(self) -> None:
        """The two end boxes, from the states the network actually has."""
        if self.graph is None:
            return
        names = [(f'{n.label or n.id} ({n.id})', n.id)
                 for n in self.graph.nodes]
        self._building = True
        try:
            for box in (self.from_dd, self.to_dd):
                want = box.value
                box.options = names
                if want in [one for _label, one in names]:
                    box.value = want
        finally:
            self._building = False
        self._on_ends(None)

    def _on_ends(self, _change) -> None:
        """Which ways there are between the two ends that are chosen.

        A branch is the ordinary case here, so the third box appears when
        there is more than one way and stays away when there is one -- the
        rule the rest of this dashboard follows.  Where there is no way at
        all the press goes too, because a profile of nothing is not a
        picture a person should be offered.
        """
        if self._building or self.graph is None:
            return
        start, end = str(self.from_dd.value or ''), str(self.to_dd.value or '')
        found = (rg.routes_between(self.graph, start, end)
                 if start and end and start != end else [])
        self._routes = found
        self._building = True
        try:
            self.route_dd.options = [
                (' -> '.join((self.graph.node(one).label or one)
                             for one in way), n)
                for n, way in enumerate(found)]
        finally:
            self._building = False
        self.route_dd.layout.display = '' if len(found) > 1 else 'none'
        self.draw_btn.layout.display = '' if found else 'none'

    def _on_draw(self) -> None:
        """Draw the chosen route, and keep the picture with the evidence.

        Written into the graph's own folder as well as shown.  A figure a
        person can only look at is a figure they will redraw by hand for
        the paper; one that is beside the calculations it was made from,
        under a name saying which route and which level, is one they can
        take -- and it travels with the folder, like everything else here.
        """
        if self.graph is None or not getattr(self, '_routes', None):
            return
        which = int(self.route_dd.value or 0)
        way = self._routes[min(which, len(self._routes) - 1)]
        level = str(self.level_dd.value or '')
        found = rg.profile(self.graph, way, level)
        if not [one for one in found['points'] if one['kcal'] is not None]:
            self.picture.value = ''
            self._say(f'Nothing on that route is priced at {level} yet.')
            return
        window = str(self.window_dd.value)
        T = float(self.temperature.value)
        ceiling = thermal_ceiling(T, self.window_seconds)
        title = (f'{self.graph.name}: '
                 + ' -> '.join((self.graph.node(one).label or one)
                               for one in way)
                 + f'  at {level}')
        note = reaction_profile.said_about(
            found['points'], missing=found['missing'], ceiling=ceiling,
            window=window)
        drawn = dict(level=level, title=title, ceiling=ceiling,
                     ceiling_label=f'{T:g} K, {window}',
                     missing=found['missing'])
        self.picture.value = reaction_profile.profile_html(
            found['points'], note=note, **drawn)
        kept = self._keep_the_picture(way, level, found, drawn)
        self._say(note + (f' Kept as {kept}.' if kept else ''))

    def _keep_the_picture(self, way, level, found, drawn) -> str:
        """The PNG into the graph folder, and a line in the history."""
        stem = rg.safe_name('-'.join(way) + '-' + level)
        relative = f'{FIGURES}/{stem}.png'
        try:
            (self.graph.folder / FIGURES).mkdir(parents=True, exist_ok=True)
            self.graph.path(relative).write_bytes(
                reaction_profile.profile_png(found['points'], **drawn))
        except Exception as exc:                 # noqa: BLE001
            self._say(f'The picture could not be kept: {exc}')
            return ''
        rg.remember(self.graph, 'profile drawn', route=list(way),
                    level=level, figure=relative,
                    missing=list(found['missing']))
        return relative

    # -- what came back while nobody was looking ---------------------------

    def harvest(self, *, say: bool = True) -> List[str]:
        """Look at every calculation this graph is waiting for.

        Asked of the folders and of nothing else.  A scheduler forgets a
        job an hour after it ends, a process is gone, and a handle in
        memory never survived the night -- but the folder is still there,
        and the document has been pointing at it the whole time.

        A calculation that is still running is left alone.  One that
        finished becomes a record; one that failed becomes a line in the
        history saying so, which is a result and not a gap: a graph that
        quietly showed nothing would invite the same two days to be spent
        again next week by whoever read that nothing.
        """
        if self.graph is None:
            return []
        told: List[str] = []
        for entry in list(self.graph.pending):
            said = run_results.what_a_run_says(self.graph.path(entry.run))
            if said['state'] in ('nothing', 'running'):
                continue
            if said['state'] == 'failed':
                rg.settle_pending(self.graph, entry.id,
                                  failed=said['why'] or 'it stopped')
                told.append(f'{entry.on}: {entry.level} stopped. '
                            f'{said["why"]}')
                continue
            xyz = said.get('xyz')
            if not xyz:
                # It finished and wrote no geometry.  A single point does
                # that, and so does an optimisation that stopped at the
                # last moment -- so the structure that went in is what the
                # numbers are about, and that is what the record carries.
                xyz = (rg.geometry(self.graph, entry.on, entry.level)
                       or self._sent_geometry(entry))
            if not xyz:
                rg.settle_pending(
                    self.graph, entry.id,
                    failed='it finished and left no geometry to attach')
                told.append(f'{entry.on}: {entry.level} finished with no '
                            f'structure to attach.')
                continue
            record = rg.Record(
                level=entry.level, method='orca',
                energy=said.get('energy'),
                free_energy=said.get('free_energy'),
                imaginary=said.get('imaginary'),
                frequency=said.get('frequency'),
                source={'kind': 'run', 'run': entry.run,
                        'output': said.get('output')})
            rg.settle_pending(self.graph, entry.id, xyz=xyz, record=record)
            told.append(f'{entry.on}: {entry.level} landed.')
        if told:
            try:
                rg.save(self.graph)
            except Exception as exc:             # noqa: BLE001
                told.append(f'and it could not be saved: {exc}')
        if told and say:
            self._say(' '.join(told))
        return told

    def _sent_geometry(self, entry) -> str:
        """The structure that was sent for this calculation, if it is there.

        Written into the run folder when the job was set up, so a run that
        produced numbers and no geometry still has something true to hang
        them on -- which is what a single point is.
        """
        try:
            return (self.graph.path(entry.run) / 'from_graph.xyz').read_text(
                encoding='utf-8')
        except OSError:
            return ''

    def _on_harvest(self) -> None:
        told = self.harvest()
        self._reopen()
        if not told:
            waiting = len(self.graph.pending) if self.graph else 0
            self._say(f'Nothing has come back yet; {waiting} still '
                      f'running.' if waiting else 'Nothing is running.')

    # -- back to the workbench --------------------------------------------

    def _submit_refs(self) -> Dict[str, Any]:
        """The Submit tab's widgets, or an empty dict when there is no tab.

        Asked each time.  This panel is built while the dashboard is being
        assembled and is also used on its own in tests, so "the Submit tab"
        is a thing that may or may not exist and the answer is allowed to be
        no.
        """
        return (getattr(self.ctx, 'submit_refs', None) or {})

    def _on_to_editor(self, _button=None) -> None:
        """This geometry into the Submit tab's box, and go and look at it.

        Written rather than handed over: the editor reads the box, and every
        other route into it -- a conversion, a drawing, a file -- writes there
        too, so a structure arriving from the graph is a structure like any
        other and the editor needs to know nothing about where it came from.

        The graph is not told.  What comes back later comes back through the
        hand-over the other way, with whatever the editor did to it, and a
        document that recorded "sent to the editor" would be recording an
        intention rather than a result.
        """
        if self.graph is None or not self.network.value:
            return
        ref = str(self.network.value)
        level = str(self.level_dd.value or '')
        text = rg.geometry(self.graph, ref, level)
        refs = self._submit_refs()
        box = refs.get('coords_widget')
        if not text or box is None:
            self._say('There is no geometry at this level to open.')
            return
        holder = self.graph.holder(ref)
        name = (getattr(holder, 'label', '') or ref) if holder else ref
        rows = text.splitlines()
        body = '\n'.join(rows[2:]) if len(rows) > 2 else ''
        box.value = (f'{len(rg.fingerprint(text).split("|"))}\n'
                     f'{name} at {level}, from the reaction graph\n{body}\n')
        self._say(f'{ref} is in the Submit tab at {level}.')
        self._go_to_tab('Submit Job')

    def _on_to_submit(self, _button=None) -> None:
        """This geometry into the Submit tab, and the graph notes that it went.

        The job is set up over there.  Every field the Submit tab has, and
        every workflow it knows, is what a calculation needs; a second
        copy of that here would be a smaller, staler one.  So this carries
        across what the graph knows -- the geometry, the charge, the spin
        and a name saying what the job is for -- and stops.

        What it does do is write the run folder into the graph before
        anything is submitted, so the job has somewhere to come back to.
        The folder is inside the graph: a calculation belongs to the
        reaction it was run for, and a graph whose evidence is scattered
        through a shared calculation directory cannot be handed to anyone.
        """
        if self.graph is None or not self.network.value:
            return
        ref = str(self.network.value)
        level = str(self.level_dd.value or '')
        text = rg.geometry(self.graph, ref, level)
        refs = self._submit_refs()
        box = refs.get('coords_widget')
        if not text or box is None:
            self._say('There is no geometry at this level to send.')
            return
        holder = self.graph.holder(ref)
        name = (getattr(holder, 'label', '') or ref) if holder else ref
        relative, real = rg.run_folder(self.graph, ref, level)
        rows = text.splitlines()
        body = chr(10).join(rows[2:]) if len(rows) > 2 else ''
        atoms = len(rg.fingerprint(text).split('|'))
        box.value = (f'{atoms}' + chr(10) +
                     f'{name}, for a calculation from the reaction graph'
                     + chr(10) + body + chr(10))
        # The name the job will carry, where the tab has a box for it.
        job_box = refs.get('job_name_widget')
        if job_box is not None:
            try:
                job_box.value = Path(relative).name
            except Exception:                    # noqa: BLE001
                pass
        try:
            real.mkdir(parents=True, exist_ok=True)
            (real / 'from_graph.xyz').write_text(text, encoding='utf-8')
            rg.remember(self.graph, 'sent for a calculation', ref=ref,
                        run=relative, level=level)
        except OSError as exc:
            self._say(f'The run folder could not be made: {exc}')
            return
        self._say(f'{ref} is in the Submit tab. Its run folder is '
                  f'{relative}, inside this graph.')
        self._go_to_tab('Submit Job')

    def _go_to_tab(self, title: str) -> None:
        """Show that tab, where there is one to show."""
        tabs = getattr(self.ctx, 'tabs_widget', None)
        where = (getattr(self.ctx, 'tab_indices', None) or {}).get(title)
        if tabs is not None and where is not None:
            try:
                tabs.selected_index = where
            except Exception:                        # noqa: BLE001
                pass

    def _on_new(self, _button=None) -> None:
        name = str(self.new_name.value or '').strip()
        if not name:
            self._say('Give the graph a name first.')
            return
        folder = Path(self.calc_dir) / rg.safe_name(name)
        try:
            rg.create(folder, name=name,
                      temperature=float(self.temperature.value),
                      seconds=self.window_seconds)
        except FileExistsError:
            self._say(f'{folder.name} already holds a graph.')
            return
        except Exception as exc:                     # noqa: BLE001
            self._say(f'Could not make it: {exc}')
            return
        self.new_name.value = ''
        self.refresh_graphs()
        self.graph_dd.value = str(folder)
        self._say(f'{folder.name} started, in {folder.parent}.')

    def _on_rename(self, _button=None) -> None:
        self._change('label', str(self.label_box.value or ''), 'renamed')

    def _on_note(self, _button=None) -> None:
        self._change('note', str(self.note_box.value or ''), 'note written')

    def _change(self, field: str, value: str, what: str) -> None:
        """Write one field of the selected thing, and say so in the history."""
        if self.graph is None or not self.network.value:
            return
        ref = str(self.network.value)
        holder = self.graph.holder(ref)
        if holder is None:
            return
        setattr(holder, field, value)
        try:
            rg.save(self.graph)
        except Exception as exc:                     # noqa: BLE001
            self._say(f'Not saved: {exc}')
            return
        rg.remember(self.graph, what, ref=ref, **{field: value})
        self.refresh_view()
        self._say(f'{ref}: {what}.')


def create_tab(ctx: Any):
    """Build the Reactions tab.  Returns ``(widget, refs)``.

    The refs are published onto *ctx* here rather than by the dashboard.  The
    hardcoded tabs are assigned by ``create_dashboard`` after it builds them;
    a registered tab goes through the registry, which keeps the widget and
    drops the refs -- so a tab that wants to be reachable from another one has
    to hand itself over, and this is where.
    """
    calc_dir = getattr(ctx, 'calc_dir', None) or (Path.home() / 'calc')
    panel = ReactionGraphPanel(calc_dir, ctx=ctx)
    refs = {'reaction_graph': panel}
    try:
        ctx.reaction_graph_refs = refs
    except Exception:                                # noqa: BLE001
        pass
    return panel.widget, refs


# Registered additively, so nothing built in changes and a tab that fails to
# build cannot take the dashboard down with it -- the registry catches that and
# marks it unavailable.  Named in tab_registry._BUILTIN_DYNAMIC_TABS, or the
# module is never imported and the registration never runs.
try:                                                 # pragma: no cover
    from delfin.dashboard.tab_registry import register_tab
    register_tab('reactions', 'Reactions', create_tab, order=8100)
except Exception:                                    # pragma: no cover
    pass


__all__ = ['WINDOWS', 'graphs_in', 'said_energy', 'said_gate', 'said_origin',
           'why_no_barrier',
           'network_lines', 'detail_html', 'summary_html',
           'ReactionGraphPanel', 'create_tab']
