"""A reaction network as a document on disk: what it is, and how it is kept.

The structure editor is a workbench.  One structure, now, and when the page is
closed it is gone.  That is right for a workbench and wrong for a reaction: a
mechanism is worked on for weeks, its calculations run for days, and at the end
of it somebody has to be able to hand the whole thing to a colleague and have
them arrive at the same picture.  So the network is a *document*, and this
module is the document -- read, written, and asked questions about, with no
widget anywhere in it.

What a graph is
---------------

A directed graph, and deliberately not a tree.  Reactions branch, and they also
converge: two educts into one saddle, two routes to the same product.  Stored
as a tree the same product exists twice and its two copies disagree about their
own energy, which is the one thing a network is for.  Drawn as a tree is fine.

**A node is a state**, a minimum: educt, intermediate, product.  **An edge is a
transition**, a saddle -- and it carries a geometry of its own, because a
barrier is a property of the connection and not a third thing standing in the
row between two structures.  This is what Chemoton, AutoMeKin and YARP each
settled on, and the editor already tells the two apart: ``Optimise`` reaches
one, ``To the saddle`` the other.

**A node is a species, not a geometry.**  A GFN2-optimised and a DFT-optimised
structure are two geometries of the same species, so a node carries several
:class:`Record` -- each one a geometry with the method that made it, what it
cost, and what was measured on it.  Three things follow, and all three are the
reason this is a document rather than a dictionary:

* "The structure at this node" is not a question until a method is named.
* A barrier is a record on an edge, never a number.  Measured on one
  sixteen-atom Diels-Alder the same barrier is +6.31 kcal/mol under GFN2 and
  +21.67 under g-xTB -- the same reaction, told twice, and averaging them or
  quoting whichever arrived last would be worse than either.
* An energy diagram is therefore a *selection* -- "draw this network at
  r2SCAN-3c and say where it is missing" -- rather than something derived.
  Silently mixing levels is the commonest quiet error in exactly this kind of
  figure, and :func:`missing_at` exists so that the answer can be "here is the
  network, and these four points are not at that level yet".

The folder
----------

A graph is a folder, not a file, and everything inside it is referred to by a
path relative to that folder.  Copy the folder and the whole reproducible
record travels: the structures, the calculations that produced them, and the
account of what was done in what order.

    <name>/
        graph.json          the document
        structures/         every geometry, once, as .xyz, never rewritten
        runs/               the calculations, one folder each
        history.jsonl       append-only: every change, and what it produced

Geometries are 10 to 100 kB of text apiece.  Inside the JSON the document would
be unreadable and every save a rewrite of the whole thing; written once to
``structures/`` and never edited, the folder stays diffable, and a save that is
interrupted cannot take a structure with it.

``runs/`` is inside the graph and not out in the calculation directory because
a calculation belongs to the reaction it was run for.  A job folder there is an
ordinary job folder -- the same layout every other DELFIN workflow writes, so
the calculations browser opens it, and the agent reads it, without knowing that
a graph put it there.

``history.jsonl`` is the reproducibility half.  One JSON object per line,
appended and never rewritten: what was done, to what, at what time, and what
came out.  A record in ``graph.json`` says the answer; the history says how it
was arrived at, which is the half a reviewer asks for.

Identity
--------

The one question this module refuses to answer on its own is whether two
structures are the same species.  :func:`looks_like` returns candidates and
merges nothing.  Two conformers of one species share an element column and a
bond graph and are usually the same node -- and so does a stereoisomer.  Merged
wrongly and silently, every energy downstream of that node is about something
else, and nothing anywhere says so.  The graph looks; the chemist decides.
"""

from __future__ import annotations

import json
import os
import re
import tempfile
from dataclasses import asdict, dataclass, field, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from delfin.dashboard.thermal import _THERMAL_SECONDS, thermal_ceiling

#: The document format.  Written into every graph, checked on every load, and
#: the reason a graph written today can still be opened when the shape has
#: moved on: a loader that meets a number it does not know says so instead of
#: reading half a file and inventing the rest.
FORMAT_VERSION = 1

#: What the folder is called inside, and nothing else may be assumed about it.
#: A graph is identified by its folder, so these are relative names throughout
#: -- an absolute path written into ``graph.json`` would break the moment the
#: folder was copied anywhere, which is the one thing it is designed for.
DOCUMENT = 'graph.json'
STRUCTURES = 'structures'
RUNS = 'runs'
HISTORY = 'history.jsonl'

#: Hartree to kcal/mol, so a barrier can be quoted in the unit a chemist reads.
HARTREE_TO_KCAL = 627.5094740631

_ID = re.compile(r'^([a-z]+)(\d+)$')
_UNSAFE = re.compile(r'[^A-Za-z0-9._-]+')


# ---------------------------------------------------------------------------
# What is in the document
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Record:
    """One geometry of one species, with the method that made it.

    *level* is the text a caption would carry -- ``GFN2-xTB``,
    ``r2SCAN-3c/CPCM(THF)`` -- and it is what the network is drawn *at*.  It is
    written out rather than assembled from parts because what makes two numbers
    comparable is not a method name: it is the method, the basis, the solvent
    model and the dispersion together, and the only thing that reliably knows
    all four is whoever set the calculation up.

    *method* is the short machine word beside it (``gfn2``, ``orca``), for the
    one thing code needs to do with it: decide whether it can run something
    else the same way.

    *energy* and *free_energy* are Hartree, absolute, exactly as the program
    reported them.  Nothing here converts or shifts: a relative number is only
    meaningful against a stated zero, and which zero depends on the question
    being asked, so it is worked out where the question is.

    *imaginary* and *frequency* are how a structure answers "what are you".
    Zero imaginary modes is a minimum, one is a saddle, and the frequency is
    the imaginary one in cm-1.  ``None`` means nobody has asked, which is a
    third state and not the same as zero.
    """

    level: str
    method: str = ''
    structure: str = ''          # graph-relative, e.g. 'structures/n01-gfn2.xyz'
    energy: Optional[float] = None
    free_energy: Optional[float] = None
    temperature: Optional[float] = None
    imaginary: Optional[int] = None
    frequency: Optional[float] = None
    charge: Optional[int] = None
    multiplicity: Optional[int] = None
    #: Where this came from, kept verbatim so that the history can be read
    #: without the code that wrote it: ``{'kind': 'editor', 'gesture': ...}``,
    #: ``{'kind': 'run', 'run': 'runs/n01_r2scan'}``, ``{'kind': 'scan', ...}``.
    source: Dict[str, Any] = field(default_factory=dict)
    at: str = ''
    note: str = ''
    #: A conformer of the same species at the same level, rather than a
    #: different level.  Kept as a mark on a record instead of a third layer
    #: between species and geometry: a conformer set is a handful of records
    #: with this set, and everything that reads the node picks the lowest.
    conformer: str = ''

    def is_minimum(self) -> bool:
        return self.imaginary == 0

    def is_saddle(self) -> bool:
        return self.imaginary == 1


@dataclass
class Node:
    """A species, with every geometry anyone has made of it."""

    id: str
    label: str = ''
    charge: int = 0
    multiplicity: int = 1
    records: List[Record] = field(default_factory=list)
    note: str = ''


@dataclass
class Edge:
    """A transition between two states, with the saddle geometries for it.

    *confirmed* is what says this edge is about these two states rather than
    about some other pair the saddle also connects.  A saddle optimiser returns
    a saddle whatever it was given; only pushing it down its imaginary mode
    both ways and landing on the two ends says which reaction it is the
    transition state *of* -- see :func:`climb.follow_the_mode_down`, which is
    the editor's answer to the same question.  ``None`` means nobody has
    checked, and a diagram that quotes an unconfirmed barrier ought to say so.
    """

    id: str
    source: str                  # node id
    target: str                  # node id
    records: List[Record] = field(default_factory=list)
    confirmed: Optional[Dict[str, Any]] = None
    label: str = ''
    note: str = ''


@dataclass
class Pending:
    """A calculation that has been started and has not answered yet.

    On disk, because that is the whole point: a DFT optimisation of a
    forty-atom complex runs for two days, the dashboard is restarted twice in
    that time, and the graph still has to know that the number it is waiting
    for is coming and where from.  Keyed by the run folder rather than by a job
    id, because a job id belongs to a scheduler and a folder belongs to the
    graph -- a queue that forgets a job leaves the folder, and the folder is
    what holds the answer.
    """

    id: str
    on: str                      # node id or edge id
    run: str                     # graph-relative, e.g. 'runs/n01_r2scan'
    level: str = ''              # what the record will be called when it lands
    job_name: str = ''
    job_id: str = ''
    backend: str = ''
    submitted: str = ''
    note: str = ''


@dataclass
class Graph:
    """One reaction network, and the folder it lives in."""

    folder: Path
    name: str = ''
    version: int = FORMAT_VERSION
    #: The temperature the network is read at, in kelvin.  A property of the
    #: document rather than of a control, because "which routes are open" is a
    #: question about this network and the answer belongs with it.
    temperature: float = 298.15
    #: How long a barrier is given, in seconds -- the other half of "possible".
    seconds: float = _THERMAL_SECONDS
    nodes: List[Node] = field(default_factory=list)
    edges: List[Edge] = field(default_factory=list)
    pending: List[Pending] = field(default_factory=list)
    note: str = ''
    #: The highest number handed out for each prefix, kept in the document.
    #: Not derived from what is in the graph now: an id is a name, and a node
    #: deleted whose number came round again would make two different species
    #: share a line in the history -- the one file that has to stay readable
    #: years later.  See :func:`_next_id`.
    counters: Dict[str, int] = field(default_factory=dict)

    # -- finding things ---------------------------------------------------

    def node(self, node_id: str) -> Optional[Node]:
        return next((n for n in self.nodes if n.id == node_id), None)

    def edge(self, edge_id: str) -> Optional[Edge]:
        return next((e for e in self.edges if e.id == edge_id), None)

    def holder(self, ref: str):
        """The node or edge called *ref*, whichever it is.

        One lookup rather than two, because almost everything that acts on a
        record -- adding one, reading one, choosing a level -- does not care
        which of the two it is standing on.
        """
        return self.node(ref) or self.edge(ref)

    def edges_from(self, node_id: str) -> List[Edge]:
        return [e for e in self.edges if e.source == node_id]

    def edges_into(self, node_id: str) -> List[Edge]:
        return [e for e in self.edges if e.target == node_id]

    def path(self, relative: str) -> Path:
        """A graph-relative name as a real path."""
        return self.folder / str(relative)

    def levels(self) -> List[str]:
        """Every level anything in this graph has been computed at, sorted.

        What the level box on a tab is filled from, and it comes from the
        document rather than from a list of methods somebody maintained: a
        network can only be drawn at a level something in it actually has.
        """
        seen = {r.level for holder in (*self.nodes, *self.edges)
                for r in holder.records if r.level}
        return sorted(seen)


# ---------------------------------------------------------------------------
# Making, reading, writing
# ---------------------------------------------------------------------------

def safe_name(text: Any) -> str:
    """A folder name that cannot escape, collide with a shell, or be empty."""
    cleaned = _UNSAFE.sub('_', str(text or '').strip()).strip('._-')
    return cleaned[:96] or 'graph'


def now() -> str:
    """The time, in the one format that sorts and parses everywhere."""
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def create(folder: Any, *, name: str = '', temperature: float = 298.15,
           seconds: float = _THERMAL_SECONDS) -> Graph:
    """Start a new graph in *folder*, which must not already hold one."""
    root = Path(folder)
    if (root / DOCUMENT).exists():
        raise FileExistsError(f'{root} already holds a graph')
    root.mkdir(parents=True, exist_ok=True)
    (root / STRUCTURES).mkdir(exist_ok=True)
    (root / RUNS).mkdir(exist_ok=True)
    graph = Graph(folder=root, name=str(name or root.name),
                  temperature=float(temperature), seconds=float(seconds))
    save(graph)
    remember(graph, 'created', name=graph.name)
    return graph


def load(folder: Any) -> Graph:
    """Read the graph in *folder*.

    A document from a future version is refused rather than half-read: the
    fields this version does not know are exactly the ones whose absence would
    make the answers wrong, and a network that quietly loses its edges is worse
    than one that will not open.
    """
    root = Path(folder)
    raw = json.loads((root / DOCUMENT).read_text(encoding='utf-8'))
    version = int(raw.get('version') or 0)
    if version > FORMAT_VERSION:
        raise ValueError(
            f'{root / DOCUMENT} is version {version}; this DELFIN reads up to '
            f'{FORMAT_VERSION}. Update DELFIN rather than opening it here.')
    return Graph(
        folder=root,
        name=str(raw.get('name') or root.name),
        version=version or FORMAT_VERSION,
        temperature=float(raw.get('temperature') or 298.15),
        seconds=float(raw.get('seconds') or _THERMAL_SECONDS),
        note=str(raw.get('note') or ''),
        nodes=[_node_from(one) for one in (raw.get('nodes') or ())],
        edges=[_edge_from(one) for one in (raw.get('edges') or ())],
        pending=[Pending(**_only(one, Pending))
                 for one in (raw.get('pending') or ())],
        counters={str(k): int(v)
                  for k, v in (raw.get('counters') or {}).items()},
    )


def save(graph: Graph) -> None:
    """Write the document, whole or not at all.

    Into a temporary file in the same directory and then renamed over the old
    one, which on every filesystem this runs on is atomic.  A graph is worked
    on for weeks; a save interrupted by a full disk or a closed laptop must
    leave last week's document standing rather than half of this one.
    """
    root = Path(graph.folder)
    root.mkdir(parents=True, exist_ok=True)
    body = {
        'version': FORMAT_VERSION,
        'name': graph.name,
        'temperature': graph.temperature,
        'seconds': graph.seconds,
        'note': graph.note,
        'nodes': [_node_to(one) for one in graph.nodes],
        'edges': [_edge_to(one) for one in graph.edges],
        'pending': [asdict(one) for one in graph.pending],
        'counters': dict(graph.counters),
    }
    text = json.dumps(body, indent=2, ensure_ascii=False, sort_keys=False)
    handle, temporary = tempfile.mkstemp(dir=str(root), prefix='.graph-',
                                         suffix='.json')
    try:
        with os.fdopen(handle, 'w', encoding='utf-8') as out:
            out.write(text + '\n')
            out.flush()
            os.fsync(out.fileno())
        os.replace(temporary, root / DOCUMENT)
    except BaseException:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def remember(graph: Graph, what: str, **fields: Any) -> None:
    """Append one line to the history, and never rewrite one.

    The document says what the network is; this says how it got that way.  It
    is the half a reviewer asks for and the half nobody writes down, so it is
    written by the code that makes the change rather than by the person, and
    it is append-only so that it cannot be tidied up afterwards.

    Failure here does not fail the change.  A history that cannot be written
    is worth reporting, but a graph that refuses to add a structure because a
    log line would not go down is a graph that has confused its two jobs.
    """
    line = {'at': now(), 'what': str(what)}
    line.update({k: _plain(v) for k, v in fields.items()})
    try:
        with open(Path(graph.folder) / HISTORY, 'a', encoding='utf-8') as out:
            out.write(json.dumps(line, ensure_ascii=False) + '\n')
    except OSError:
        pass


def history(graph: Graph) -> List[Dict[str, Any]]:
    """Everything that has been done to this graph, oldest first.

    A line that cannot be parsed is skipped rather than raising: the history is
    read to be shown, and one bad line from a half-written append should not
    take the account of a month's work with it.
    """
    out: List[Dict[str, Any]] = []
    try:
        text = (Path(graph.folder) / HISTORY).read_text(encoding='utf-8')
    except OSError:
        return out
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            got = json.loads(line)
        except ValueError:
            continue
        if isinstance(got, dict):
            out.append(got)
    return out


# ---------------------------------------------------------------------------
# Putting things in
# ---------------------------------------------------------------------------

def add_state(graph: Graph, xyz: str, record: Record, *, label: str = '',
              charge: int = 0, multiplicity: int = 1,
              note: str = '') -> Node:
    """A new species, with the geometry that is being put in as its first."""
    node = Node(id=_next_id(graph, 'n'), label=str(label or ''),
                charge=int(charge), multiplicity=int(multiplicity),
                note=str(note or ''))
    graph.nodes.append(node)
    add_record(graph, node.id, xyz, record)
    remember(graph, 'state added', ref=node.id, label=node.label,
             level=record.level)
    return node


def add_transition(graph: Graph, xyz: str, record: Record, *,
                   source: str, target: str, label: str = '',
                   confirmed: Optional[Dict[str, Any]] = None,
                   note: str = '') -> Edge:
    """A transition from one state to another, with its saddle geometry.

    Both ends have to exist.  An edge to a state nobody has put in is a
    reaction to something unnamed, and a network that permits it stops being
    able to answer any question about a route.
    """
    for end in (source, target):
        if graph.node(end) is None:
            raise KeyError(f'no state {end!r} in this graph')
    edge = Edge(id=_next_id(graph, 'e'), source=source, target=target,
                label=str(label or ''), confirmed=confirmed,
                note=str(note or ''))
    graph.edges.append(edge)
    add_record(graph, edge.id, xyz, record)
    # ``from_state`` and ``to_state`` rather than ``source`` and ``target``.
    # A record's line already uses ``source`` for where the number came from,
    # and one key meaning provenance on one line and a node id on the next is
    # how a log stops being readable -- which is the only thing this file is
    # for.  The dataclass keeps its own names; the history is a different
    # document with a different reader.
    remember(graph, 'transition added', ref=edge.id, from_state=source,
             to_state=target, level=record.level,
             confirmed=bool(confirmed))
    return edge


def add_record(graph: Graph, ref: str, xyz: str, record: Record) -> Record:
    """One more geometry of the thing called *ref*, at its own level.

    The structure is written to ``structures/`` under a name of its own and the
    record points at it.  Nothing is ever overwritten: a second record at the
    same level is a second file, because the first one is what an earlier
    number was measured on and a document that quietly replaced it could no
    longer account for its own history.
    """
    holder = graph.holder(ref)
    if holder is None:
        raise KeyError(f'no state or transition {ref!r} in this graph')
    relative = _write_structure(graph, ref, record.level, xyz)
    stored = replace(record, structure=relative, at=record.at or now())
    holder.records.append(stored)
    remember(graph, 'record added', ref=ref, level=stored.level,
             structure=relative, energy=stored.energy,
             free_energy=stored.free_energy, imaginary=stored.imaginary,
             source=stored.source)
    return stored


def add_pending(graph: Graph, ref: str, *, run: str, level: str,
                job_name: str = '', job_id: str = '', backend: str = '',
                note: str = '') -> Pending:
    """Say that a calculation has been started for *ref*, and where."""
    if graph.holder(ref) is None:
        raise KeyError(f'no state or transition {ref!r} in this graph')
    entry = Pending(id=_next_id(graph, 'p'), on=ref, run=str(run),
                    level=str(level), job_name=str(job_name),
                    job_id=str(job_id), backend=str(backend),
                    submitted=now(), note=str(note or ''))
    graph.pending.append(entry)
    remember(graph, 'calculation started', ref=ref, run=entry.run,
             level=entry.level, job_name=entry.job_name, job_id=entry.job_id)
    return entry


def settle_pending(graph: Graph, pending_id: str, *, xyz: str = '',
                   record: Optional[Record] = None,
                   failed: str = '') -> Optional[Record]:
    """A calculation has answered.  Take the record, or take the refusal.

    Both endings are recorded.  A job that failed is not the absence of an
    answer -- it is a fact about this structure at this level, and a graph that
    silently dropped it would invite the same calculation to be started again
    next week by somebody reading a gap.
    """
    entry = next((p for p in graph.pending if p.id == pending_id), None)
    if entry is None:
        raise KeyError(f'no pending calculation {pending_id!r}')
    graph.pending = [p for p in graph.pending if p.id != pending_id]
    if failed or record is None:
        remember(graph, 'calculation failed', ref=entry.on, run=entry.run,
                 level=entry.level, why=str(failed or 'nothing came back'))
        return None
    stored = add_record(graph, entry.on, xyz, record)
    remember(graph, 'calculation landed', ref=entry.on, run=entry.run,
             level=stored.level)
    return stored


def run_folder(graph: Graph, ref: str, level: str) -> Tuple[str, Path]:
    """Where a calculation for *ref* at *level* goes, as name and real path.

    Inside the graph, so that the folder carries its own calculations, and
    named after what it is for so that the folder is readable without the
    document.  A name already taken gets a number rather than being reused: a
    calculation is evidence, and two of them sharing a directory means neither
    can be pointed at afterwards.
    """
    stem = safe_name(f'{ref}_{level}') or safe_name(ref)
    root = Path(graph.folder) / RUNS
    relative, candidate, n = f'{RUNS}/{stem}', root / stem, 1
    while candidate.exists():
        n += 1
        relative, candidate = f'{RUNS}/{stem}_{n}', root / f'{stem}_{n}'
    return relative, candidate


# ---------------------------------------------------------------------------
# Asking it things
# ---------------------------------------------------------------------------

def best(holder, level: str) -> Optional[Record]:
    """The record of *holder* at *level* to quote, or None.

    The lowest one, which is how a conformer set answers for its species: a
    handful of records at the same level marked as conformers, and the number
    that belongs in a diagram is the bottom of them.  Free energy decides where
    there is one, because that is what the diagram is about; electronic energy
    where there is not, so a network computed without frequencies still draws.
    """
    if holder is None:
        return None
    same = [r for r in holder.records if r.level == level]
    if not same:
        return None
    priced = [r for r in same if r.free_energy is not None]
    if priced:
        return min(priced, key=lambda r: r.free_energy)
    priced = [r for r in same if r.energy is not None]
    if priced:
        return min(priced, key=lambda r: r.energy)
    return same[0]


def _value(record: Optional[Record]) -> Optional[float]:
    """What a record is worth, in Hartree: free energy if it has one."""
    if record is None:
        return None
    if record.free_energy is not None:
        return float(record.free_energy)
    if record.energy is not None:
        return float(record.energy)
    return None


def barrier(graph: Graph, edge_id: str, level: str) -> Optional[float]:
    """The forward barrier over *edge_id* at *level*, in kcal/mol.

    Against the state the edge leaves, at the same level, and ``None`` unless
    all three of the saddle, the state and the level are there.  Refusing is
    the point: a barrier assembled out of two levels is a number with no
    meaning that looks exactly like one with meaning.
    """
    edge = graph.edge(edge_id)
    if edge is None:
        return None
    top = _value(best(edge, level))
    foot = _value(best(graph.node(edge.source), level))
    if top is None or foot is None:
        return None
    return (top - foot) * HARTREE_TO_KCAL


def reaction_energy(graph: Graph, edge_id: str,
                    level: str) -> Optional[float]:
    """What the step itself costs or releases at *level*, in kcal/mol."""
    edge = graph.edge(edge_id)
    if edge is None:
        return None
    foot = _value(best(graph.node(edge.source), level))
    head = _value(best(graph.node(edge.target), level))
    if foot is None or head is None:
        return None
    return (head - foot) * HARTREE_TO_KCAL


def open_at(graph: Graph, level: str, *, temperature: Optional[float] = None,
            seconds: Optional[float] = None) -> Dict[str, Optional[bool]]:
    """Which transitions the temperature pays for, by edge id.

    ``True`` open, ``False`` closed, ``None`` unanswerable at this level -- and
    the third is not a detail.  A network drawn with the unknown ones shown as
    closed is a network that reports a conclusion about calculations nobody has
    run, which is the failure this whole module is arranged against.

    The ceiling is :func:`thermal.thermal_ceiling`, the same instrument the
    editor holds a drag against; here it is held against a week's work.
    """
    T = float(graph.temperature if temperature is None else temperature)
    window = float(graph.seconds if seconds is None else seconds)
    ceiling = thermal_ceiling(T, window)
    out: Dict[str, Optional[bool]] = {}
    for edge in graph.edges:
        rise = barrier(graph, edge.id, level)
        out[edge.id] = None if rise is None else bool(rise <= ceiling)
    return out


def missing_at(graph: Graph, level: str) -> List[str]:
    """Everything in the network with no record at *level*, by id.

    What makes an energy diagram honest.  The question a network is drawn to
    answer is comparative, so the useful statement is not "here is the diagram"
    but "here is the diagram, and these four points are not at this level yet".
    """
    return [holder.id for holder in (*graph.nodes, *graph.edges)
            if best(holder, level) is None]


def looks_like(graph: Graph, xyz: str, *, charge: Optional[int] = None,
               multiplicity: Optional[int] = None) -> List[Node]:
    """States this structure might be, and it decides nothing.

    Same element column, same charge and multiplicity where they were given.
    That is a shortlist and never an answer: two conformers of a species match
    on all three and are usually one node, and so do two stereoisomers, which
    are never one.  Merged wrongly and silently, every energy downstream is
    about a different molecule and nothing says so -- so this returns the
    candidates and a person picks, or picks none.
    """
    want = fingerprint(xyz)
    if not want:
        return []
    out: List[Node] = []
    for node in graph.nodes:
        if charge is not None and node.charge != int(charge):
            continue
        if multiplicity is not None and node.multiplicity != int(multiplicity):
            continue
        for record in node.records:
            try:
                text = graph.path(record.structure).read_text(encoding='utf-8')
            except OSError:
                continue
            if fingerprint(text) == want:
                out.append(node)
                break
    return out


def fingerprint(xyz: Any) -> str:
    """The element column of an XYZ block -- what makes it the same molecule.

    The same rule the editor identifies a structure by, and for the same
    reason: a row counts only where a symbol is followed by three numbers, so a
    comment line of free text cannot enter the fingerprint and make one
    molecule look like two.
    """
    out: List[str] = []
    for line in str(xyz or '').splitlines():
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            [float(one) for one in parts[1:4]]
        except ValueError:
            continue
        out.append(parts[0].strip().capitalize())
    return '|'.join(out)


def geometry(graph: Graph, ref: str, level: str) -> Optional[str]:
    """The XYZ text of *ref* at *level*, ready for the editor or an input."""
    record = best(graph.holder(ref), level)
    if record is None or not record.structure:
        return None
    try:
        return graph.path(record.structure).read_text(encoding='utf-8')
    except OSError:
        return None


def route(graph: Graph, node_ids: Iterable[str]) -> List[Edge]:
    """The edges joining those states in order, for a diagram along one path.

    Raises where two named states are not joined, rather than skipping: a
    profile drawn over a gap it did not mention is a picture of a reaction that
    does not happen.
    """
    ids = [str(one) for one in node_ids]
    out: List[Edge] = []
    for first, second in zip(ids, ids[1:]):
        edge = next((e for e in graph.edges
                     if e.source == first and e.target == second), None)
        if edge is None:
            raise KeyError(f'no transition from {first!r} to {second!r}')
        out.append(edge)
    return out


# ---------------------------------------------------------------------------
# The small things underneath
# ---------------------------------------------------------------------------

def _write_structure(graph: Graph, ref: str, level: str, xyz: str) -> str:
    """One geometry into ``structures/``, under a name nothing else has."""
    text = str(xyz or '')
    if not fingerprint(text):
        raise ValueError('that is not a structure: no atom rows in it')
    root = Path(graph.folder) / STRUCTURES
    root.mkdir(parents=True, exist_ok=True)
    stem = safe_name(f'{ref}-{level}') or safe_name(ref)
    name, n = f'{stem}.xyz', 1
    while (root / name).exists():
        n += 1
        name = f'{stem}-{n}.xyz'
    (root / name).write_text(text if text.endswith('\n') else text + '\n',
                             encoding='utf-8')
    return f'{STRUCTURES}/{name}'


def _next_id(graph: Graph, prefix: str) -> str:
    """The next id with that prefix, and it is never handed out twice.

    From a counter kept in the document rather than from what is in the graph
    now, because an id is a name.  Derived from the highest one present, a node
    deleted would free its number, the next node would take it, and the two
    would share every line the history has about either -- and the history is
    the one file that has to still be readable years later.

    What is present is consulted as well, and only ever to raise the counter: a
    file edited by hand, or merged from two copies of the same graph, can carry
    ids the counter has not reached, and a name that collides is worse than a
    gap in the numbering.
    """
    highest = int(graph.counters.get(prefix) or 0)
    for one in (*graph.nodes, *graph.edges, *graph.pending):
        got = _ID.match(str(one.id))
        if got and got.group(1) == prefix:
            highest = max(highest, int(got.group(2)))
    graph.counters[prefix] = highest + 1
    return f'{prefix}{highest + 1:02d}'


def _only(raw: Dict[str, Any], kind) -> Dict[str, Any]:
    """The keys of *raw* that *kind* actually has.

    A document written by a newer DELFIN that added a field is refused at the
    version check; this is for the other direction -- a hand-edited file with
    something extra in it opens, without the extra silently becoming an
    attribute nothing reads.
    """
    known = {f for f in kind.__dataclass_fields__}
    return {k: v for k, v in (raw or {}).items() if k in known}


def _record_from(raw: Dict[str, Any]) -> Record:
    return Record(**_only(raw, Record))


def _node_from(raw: Dict[str, Any]) -> Node:
    fields = _only(raw, Node)
    fields['records'] = [_record_from(one) for one in (raw.get('records') or ())]
    return Node(**fields)


def _edge_from(raw: Dict[str, Any]) -> Edge:
    fields = _only(raw, Edge)
    fields['records'] = [_record_from(one) for one in (raw.get('records') or ())]
    return Edge(**fields)


def _node_to(node: Node) -> Dict[str, Any]:
    out = asdict(node)
    out['records'] = [asdict(one) for one in node.records]
    return out


def _edge_to(edge: Edge) -> Dict[str, Any]:
    out = asdict(edge)
    out['records'] = [asdict(one) for one in edge.records]
    return out


def _plain(value: Any) -> Any:
    """Whatever it is, as something json can write."""
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(one) for one in value]
    return str(value)
