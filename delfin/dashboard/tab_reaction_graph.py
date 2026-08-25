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

#: Enough of a geometry to recognise it, and not enough to fill the panel.
_PREVIEW_ROWS = 12


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

def said_energy(record: Optional[rg.Record]) -> str:
    """One record as a line, with the quantity it actually is.

    ``G`` and ``E`` are printed under their own names rather than both as
    "energy", because which one a diagram was drawn from is the difference
    between a free-energy profile and an electronic one, and nobody can tell
    them apart from the picture afterwards.
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
    if record.conformer:
        parts.append(f'conformer {record.conformer}')
    return ', '.join(parts)


def said_gate(verdict: Optional[bool], rise: Optional[float],
              window: str) -> str:
    """Whether the temperature pays for that barrier, in three answers.

    The third is the one that has to be said rather than implied.  A network
    drawn with the uncomputed edges greyed out like the closed ones reports a
    conclusion about calculations nobody has run, and the reader has no way to
    tell which grey is which.
    """
    if verdict is None:
        return 'not at this level'
    if rise is None:
        return 'not at this level'
    if verdict:
        return f'open within {window}'
    return f'closed within {window}'


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
                        f'{said_gate(gates.get(edge.id), rise, window)}{mark}'))
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
                    + (f', {unknown} not at this level.' if unknown
                       else '.'))
    if missing:
        shown = ', '.join(missing[:8]) + ('...' if len(missing) > 8 else '')
        said.append(f'Nothing at {level} for: {shown}')
    if graph.pending:
        said.append(f'{_many(len(graph.pending), "calculation")} running.')
    return ('<div style="font-family:sans-serif; font-size:13px">'
            + ' '.join(_escape(one) for one in said) + '</div>')


# ---------------------------------------------------------------------------
# The tab
# ---------------------------------------------------------------------------

class ReactionGraphPanel:
    """One Reactions tab, over the graphs in a calculation directory."""

    def __init__(self, calc_dir: Any = None):
        self.calc_dir = Path(calc_dir) if calc_dir else Path.home() / 'calc'
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

        self.label_box = widgets.Text(placeholder='name',
                                      layout=widgets.Layout(width='220px'))
        self.rename_btn = widgets.Button(description='Rename',
                                         layout=widgets.Layout(width='100px'))
        self.note_box = widgets.Text(placeholder='note',
                                     layout=widgets.Layout(width='320px'))
        self.note_btn = widgets.Button(description='Save note',
                                       layout=widgets.Layout(width='110px'))

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
                          self.reload_btn],
                         layout=widgets.Layout(flex_flow='row wrap')),
            widgets.HBox([self.level_dd, self.temperature, self.window_dd],
                         layout=widgets.Layout(flex_flow='row wrap')),
            self.summary,
            widgets.HBox([
                widgets.VBox([self.network],
                             layout=widgets.Layout(width='52%')),
                widgets.VBox([self.detail,
                              widgets.HBox([self.label_box, self.rename_btn],
                                           layout=widgets.Layout(
                                               flex_flow='row wrap')),
                              widgets.HBox([self.note_box, self.note_btn],
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
        self._refresh_levels()
        self.temperature.value = float(self.graph.temperature)
        self.refresh_view()

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
    """Build the Reactions tab.  Returns ``(widget, refs)``."""
    calc_dir = getattr(ctx, 'calc_dir', None) or (Path.home() / 'calc')
    panel = ReactionGraphPanel(calc_dir)
    return panel.widget, {'reaction_graph': panel}


# Registered additively, so nothing built in changes and a tab that fails to
# build cannot take the dashboard down with it -- the registry catches that and
# marks it unavailable.  Named in tab_registry._BUILTIN_DYNAMIC_TABS, or the
# module is never imported and the registration never runs.
try:                                                 # pragma: no cover
    from delfin.dashboard.tab_registry import register_tab
    register_tab('reactions', 'Reactions', create_tab, order=8100)
except Exception:                                    # pragma: no cover
    pass


__all__ = ['WINDOWS', 'graphs_in', 'said_energy', 'said_gate',
           'network_lines', 'detail_html', 'summary_html',
           'ReactionGraphPanel', 'create_tab']
