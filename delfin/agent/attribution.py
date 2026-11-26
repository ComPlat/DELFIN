"""What a model call was FOR — the attribution the metric rows were missing.

``turn_metrics`` and ``agent_metrics`` record how a call went (tokens, ttft,
wall time) and who answered it (provider, model, role). ``budget_ledger``
charges a delegate's spend to the calling session. None of them records **what
the call served**: a session that spent 200 turns over three molecules and two
workflow runs is one undifferentiated block in the store, and "what did this
molecule cost?" has no answer in the data — only in somebody's memory.

This module holds that context: an immutable :class:`Attribution` of
``molecule_id``, ``workflow_id``, ``task_id``, ``role``, ``subagent_name`` and
``parent_session``, set where the work is known (a workflow run, a calculation
job, a delegated subagent, a swarm member) and read by the metric writers on
every call made under it.

Two mechanisms, because the callers live on both sides of a thread boundary:

* A :class:`contextvars.ContextVar` — the inherited path. A nested generator,
  a tool call, or a model request issued deep inside the client all see the
  caller's context without any of them taking a parameter. This is how a call
  a TOOL triggers gets attributed: its writer never heard of attribution and
  nothing passes it a value.
* A value bound on the client object — the path that survives a thread.
  DELFIN spawns a background sub-agent with ``threading.Thread`` and calls the
  runner *inside* that thread (``api_client.py``: ``_bg_run`` ->
  ``subagent_runner``), and a new thread starts with an EMPTY context — that is
  what ``ContextVar`` means, not an implementation detail. So an ambient-only
  design would silently record ``None`` for exactly the delegated calls this
  feature exists to attribute. :func:`bind` therefore writes the client too,
  :func:`inherit` reads it back in the child, and :func:`for_client` prefers it.

Three rules the data depends on:

**Never guess an id.** A field nobody set stays ``None``, and ``None`` is
written as ``None`` rather than omitted: an absent key means "this recorder had
no such field", which is a different statement from "nothing was known" — the
same distinction ``turn_metrics`` draws for ``first_event_ms`` and ``retries``.
An invented id is worse than a null, because it bills one molecule for another's
compute and nothing downstream can tell the two apart.

**Identifiers only.** The metric stores are kept for months and are read by
other people. A molecule *name*, a SMILES string, a filesystem path or an
account must never reach them — so a value that does not look like an opaque
identifier is dropped, not truncated. Truncating leaks anyway: a home path
like ``/home/user/...`` cut to 20 characters still names a user's directory.

**Inherit, do not re-supply.** :func:`child_of` copies the parent's ids and
changes only who is running, so a delegation cannot lose the ids by forgetting
to repeat them, and the parent's record is never mutated on the way.
"""

from __future__ import annotations

import contextvars
import dataclasses
import os
import re
from typing import Any

__all__ = [
    "Attribution",
    "EMPTY",
    "COLUMNS",
    "attribution_scope",
    "current",
    "bind",
    "bound_to",
    "inherit",
    "for_client",
    "child_of",
    "columns",
    "clean_id",
    "open_scope",
    "close_scope",
    "from_env",
    "adopt_process",
    "reset_process",
    "turn_scope",
    "ENV_RUN",
    "ENV_LABEL",
]

#: The two variables ``swarm/wave.py`` puts in every member process
#: (``_launch``), and the only identity such a process is given: its own
#: run and its own member id. ``ENV_LABEL`` is ``"<run>/<member>"``. They are
#: read rather than re-derived because a member is a SEPARATE process with a
#: fresh context: without this, every model call a member makes would record
#: six nulls, which is the same defect this module exists to close, one
#: process further out.
ENV_RUN = "DELFIN_SWARM_RUN"
ENV_LABEL = "DELFIN_AGENT_LABEL"

#: The six columns this adds to a metric row, in the order they are written.
#: ``role`` is published as ``agent_role`` because ``role`` is already a
#: turn_metrics key with its own meaning (the engine's routed role).
COLUMNS: tuple[str, ...] = (
    "molecule_id", "workflow_id", "task_id",
    "agent_role", "subagent_name", "parent_session",
)

#: An identifier is a short opaque token: an account name, a calc folder, a
#: run id. Everything DELFIN actually uses (``mol_2026-12``, ``esd-chain``,
#: ``t3a``, ``calc.4``, ``run:12``) is in this set; a path (``/``, ``\``), a
#: SMILES string (``(``, ``=``, ``)``), an e-mail (``@``) or a name with a
#: space (``benzene dimer complex``) is not.
_ID_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:+#-]{0,79}")

#: A leading directory separator is the commonest leak and the one a bare
#: ``^[A-Za-z0-9]`` already rejects; named so the intent survives reading.
_MAX_ID_LEN = 80

#: "this client carried no binding", distinct from one carrying ``None``.
_MISSING = object()


def clean_id(value: Any) -> str | None:
    """``value`` as a storable identifier, or ``None`` when it is not one.

    Deliberately strict and deliberately lossy in one direction only: a value
    that is not an identifier becomes ``None``, never a shortened version of
    itself. A caller that passes a molecule name gets a null column and can
    find out it passed a name; it does not get its path into a store that is
    kept for months.
    """
    if not isinstance(value, str):
        return None
    text = value.strip()
    if not text or len(text) > _MAX_ID_LEN:
        return None
    if _ID_RE.fullmatch(text) is None:
        return None
    # A path is rejected by the charset; this catches a Windows drive or a
    # UNC root that slipped through as a bare token.
    if ":" in text and not text[:1].isalpha():
        return None
    return text


@dataclasses.dataclass(frozen=True)
class Attribution:
    """What one model call served. Immutable, all fields optional.

    Frozen because the record is inherited: a delegate that "just sets its own
    role" must not be able to rewrite what its parent's other calls mean.
    """

    molecule_id: str | None = None
    workflow_id: str | None = None
    task_id: str | None = None
    role: str | None = None
    subagent_name: str | None = None
    parent_session: str | None = None

    def is_empty(self) -> bool:
        """Nothing known at all — the state a row with six nulls records."""
        return not any((self.molecule_id, self.workflow_id, self.task_id,
                        self.role, self.subagent_name, self.parent_session))


EMPTY = Attribution()

#: The attribute a bound client carries. Named for what it is: an internal
#: attribute, like the client's other per-session state (``metrics_session``).
CLIENT_ATTR = "_delfin_attribution"

_CURRENT: "contextvars.ContextVar[Attribution]" = contextvars.ContextVar(
    "delfin_attribution", default=EMPTY)

#: Who this PROCESS works for, when DELFIN started it rather than a person.
#: A plain module global, deliberately NOT a ``ContextVar``: a member's
#: identity is a property of the process, and a new thread starts with an
#: empty context. ``repl.py`` runs every turn on a fresh ``threading.Thread``,
#: so an identity adopted in ``main``'s context alone was measured to vanish
#: before the turn wrote its row (``workflow_id=None`` for a member process).
_PROCESS: Attribution = EMPTY


def _effective() -> Attribution:
    """The attribution in force here: this context's, on top of the process's.

    A context value wins over the process value field by field, so a scope
    inside a member's turn can narrow the run it belongs to without having to
    restate it.
    """
    here = _CURRENT.get() or EMPTY
    if _PROCESS.is_empty():
        return here
    keep = dataclasses.asdict(_PROCESS)
    for name, value in dataclasses.asdict(here).items():
        if value is not None:
            keep[name] = value
    return Attribution(**keep)


def current() -> Attribution:
    """The attribution in force here, or ``None`` when nothing is known.

    Returns ``None`` rather than :data:`EMPTY` so a caller cannot read an empty
    record as "attributed to nothing" and skip the column: the two are the same
    to a writer, and different to a reader asking "was this ever attributed?".
    """
    got = _effective()
    return None if got is None or got.is_empty() else got


def adopt_process(**fields: Any) -> Attribution:
    """Make ``fields`` the identity of this process, for every thread in it.

    Called once from the process entry point (``cli.main``) with what the
    environment states: a swarm member's run and member. Threads started later
    inherit it through :func:`_effective`, a context :func:`set_attribution`
    would not reach.
    """
    global _PROCESS
    _PROCESS = _merged(_PROCESS, **fields)
    return _PROCESS


def reset_process() -> None:
    """Forget the process identity. For a test, and for a caller that ends one
    run's work in a process it keeps."""
    global _PROCESS
    _PROCESS = EMPTY


def _merged(base: Attribution, **fields: Any) -> Attribution:
    """``base`` with the given fields replaced, every value cleaned."""
    keep = dataclasses.asdict(base)
    for name, value in fields.items():
        if value is None:
            continue                    # absent says nothing; do not erase
        cleaned = clean_id(value)
        if cleaned is None:
            continue                    # rejected: keep the parent's, never guess
        keep[name] = cleaned
    return Attribution(**keep)


class attribution_scope:
    """Attribute every model call made inside this block.

    Used by the code that knows the work — a workflow run, a calculation job, a
    swarm member::

        with attribution.attribution_scope(molecule_id=calc.molecule_id,
                                           workflow_id=run_id):
            engine.stream_response(...)

    Nesting replaces: an inner scope sees the outer's ids through
    :func:`current` and can merge them (``child_of``), and the outer one is
    restored on exit — a molecule's ids never survive into unrelated work.

    A value that is not an identifier is dropped by :func:`clean_id`, so this
    never raises and never stores a path: a caller that passes a name gets a
    null column.
    """

    __slots__ = ("_new", "_token", "_active")

    def __init__(self, *, molecule_id: Any = None, workflow_id: Any = None,
                 task_id: Any = None, role: Any = None,
                 subagent_name: Any = None, parent_session: Any = None,
                 inherit: bool = True) -> None:
        base = _effective() if inherit else EMPTY
        if base is None:
            base = EMPTY
        self._new = _merged(
            base, molecule_id=molecule_id, workflow_id=workflow_id,
            task_id=task_id, role=role, subagent_name=subagent_name,
            parent_session=parent_session)
        self._token: Any = None
        self._active = False

    def __enter__(self) -> Attribution:
        self._token = _CURRENT.set(self._new)
        self._active = True
        return self._new

    def __exit__(self, exc_type, exc, tb) -> bool:
        if self._active:
            _CURRENT.reset(self._token)
            self._active = False
        return False


def set_attribution(**fields: Any) -> Attribution:
    """Set the attribution for the rest of THIS context (no restore).

    The escape hatch for a caller that has no ``with`` block to put it in — a
    long-lived worker that starts one run per task. Prefer
    :class:`attribution_scope`: a scope restores, and a leaked id silently
    bills the next unrelated call.
    """
    base = _CURRENT.get() or EMPTY
    new = _merged(base, **fields)
    _CURRENT.set(new)
    return new


def open_scope(client: Any = None, **fields: Any):
    """A restorable scope for a caller with a ``try``/``finally``, not a ``with``.

    The engine's turn is that shape. Restoring matters more than it looks: a
    workflow run that leaves its ``molecule_id`` in a long-lived dashboard
    thread's context would bill every later call on that thread to a molecule
    it has nothing to do with, and the row would look authoritative.

    With ``client`` the scope also binds it, so a call made on that client from
    another thread still finds the ids.
    """
    scope = attribution_scope(**fields)
    got = scope.__enter__()
    had = _MISSING
    if client is not None:
        had = getattr(client, CLIENT_ATTR, _MISSING)
        try:
            setattr(client, CLIENT_ATTR, got)
        except Exception:
            had = _MISSING
    return scope, got, client, had


def close_scope(handle) -> None:
    """Undo :func:`open_scope`, the client binding included. Never raises.

    The client binding has to go back, not just stay: the engine holds ONE
    client for a whole session, so a binding left by turn 1 would report turn
    2 — which may know nothing — as belonging to turn 1's molecule. Restoring a
    client that had no binding removes the attribute entirely, which reads as
    "nothing known", the honest state.
    """
    try:
        scope, _got, client, had = handle
    except Exception:
        return
    try:
        scope.__exit__(None, None, None)
    except Exception:
        pass
    if client is None:
        return
    try:
        if had is _MISSING:
            delattr(client, CLIENT_ATTR)
        else:
            setattr(client, CLIENT_ATTR, had)
    except Exception:
        pass


def bind(client: Any, **fields: Any) -> Attribution:
    """Attribute this client's calls, on the thread AND on the client.

    Called where the turn starts, beside the existing
    ``client.metrics_session`` assignment: the engine knows its role and
    session, and the ids it merges in come from whoever wrapped the turn.

    Writing the client is what makes a delegation inheritable. ``run_subagent``
    takes the parent's client by shallow copy, so the attribute travels with
    it into a new thread, where the context variable does not exist yet.
    """
    base = _CURRENT.get() or EMPTY
    new = _merged(base, **fields)
    _CURRENT.set(new)
    try:
        setattr(client, CLIENT_ATTR, new)
    except Exception:
        pass                            # a client that refuses attributes
    return new


def bound_to(client: Any) -> Attribution | None:
    """The attribution bound to ``client``, or ``None``."""
    got = getattr(client, CLIENT_ATTR, None)
    if isinstance(got, Attribution) and not got.is_empty():
        return got
    return None


def inherit(client: Any) -> Attribution | None:
    """What the calls under ``client`` should carry: its binding, else here.

    A child thread has no context, so a delegate's ids come off the client it
    was copied from. Falling back to the ambient context keeps a foreground
    delegation that never called :func:`bind` working.
    """
    return bound_to(client) or current()


def for_client(client: Any) -> Attribution | None:
    """The attribution to record for a call made on ``client``.

    The client's own binding wins: it says what the client was handed to work
    on, which is what a row should say, and it is the value that survives a
    thread. ``None`` means nothing is known — the writer records nulls.
    """
    return bound_to(client) or current()


def child_of(parent: Attribution | None, *, role: Any = "subagent",
             subagent_name: Any = "", parent_session: Any = "") -> Attribution:
    """The delegate's record: the parent's ids, a new owner.

    Copying rather than starting fresh is the point — a child row that lost its
    parent's ``molecule_id`` would put the delegate's spend outside the answer
    to "what did this molecule cost?", which is the hole this whole module
    fills. Only ``role``, ``subagent_name`` and ``parent_session`` change, and
    the parent's record is not touched.
    """
    base = parent if isinstance(parent, Attribution) else EMPTY
    return _merged(base, role=role, subagent_name=subagent_name,
                   parent_session=parent_session)


def columns(attribution: Attribution | None = None,
            client: Any = None) -> dict[str, str | None]:
    """The six metric columns for a call, always all six, null when unknown.

    Read by the writers, which is why nothing here can raise: a metrics write
    that throws must not break the turn it describes, and these functions are
    called from inside those writes.
    """
    got = attribution
    if got is None and client is not None:
        got = for_client(client)
    if got is None:
        got = current()
    if not isinstance(got, Attribution):
        got = EMPTY
    return {
        "molecule_id": got.molecule_id,
        "workflow_id": got.workflow_id,
        "task_id": got.task_id,
        "agent_role": got.role,
        "subagent_name": got.subagent_name,
        "parent_session": got.parent_session,
    }


def turn_scope(engine: Any):
    """The attribution for ONE turn of ``engine``, as a context manager.

    Called at the turn ENTRY POINTS (``cli._run_once``, ``repl.run_turn``)
    rather than inside ``engine.stream_response``: the engine and the client
    are on the protected list, so a change there needs the patch route, and
    the entry point is where a turn is actually started anyway.

    What it contributes is only what the engine itself KNOWS: the role it is
    running as and the session that pays. The molecule / workflow / task come
    from whoever wrapped the turn (:class:`attribution_scope`), are merged in
    by :func:`_merged`, and stay null when nobody set them — a guessed id
    bills one molecule for another's compute.

    The binding goes on the CLIENT as well as in the context, because a
    shallow copy of that client carries it into a delegate's thread, where a
    ``ContextVar`` does not exist. It is taken back off on exit, so the next
    turn of the same long-lived client does not inherit this one's molecule.

    Never raises, and works whether or not ``engine`` has a usable client: a
    turn that loses its attribution is a smaller loss than a turn that fails
    to start.
    """
    return _TurnScope(engine)


class _TurnScope:
    """See :func:`turn_scope`. Kept a class so ``with`` guarantees the reset."""

    __slots__ = ("_engine", "_handle")

    def __init__(self, engine: Any) -> None:
        self._engine = engine
        self._handle: Any = None

    def __enter__(self):
        try:
            self._handle = open_scope(
                getattr(self._engine, "client", None),
                role=str(getattr(self._engine, "current_role", "") or ""),
                parent_session=str(_session_of(self._engine) or ""))
        except Exception:
            self._handle = None
        return self._handle

    def __exit__(self, exc_type, exc, tb) -> bool:
        if self._handle is not None:
            close_scope(self._handle)
            self._handle = None
        return False


def _session_of(engine: Any) -> str:
    """The engine's session key, by whichever accessor it exposes."""
    for name in ("trace_session", "session_id"):
        try:
            got = getattr(engine, name, None)
            value = got() if callable(got) else got
            if value:
                return str(value)
        except Exception:
            continue
    return ""


def from_env(environ: "dict | None" = None) -> Attribution:
    """Attribution for a process DELFIN itself started, or :data:`EMPTY`.

    A swarm member is its own process: ``swarm/wave.py:_launch`` hands it
    ``DELFIN_SWARM_RUN`` and ``DELFIN_AGENT_LABEL`` (``"<run>/<member>"``) and
    nothing else. Those two are its whole identity, and a fresh process starts
    with an empty context, so reading them here is what keeps a member's model
    calls attributable to the run it works for. Called once at session start.

    An ordinary session sets neither variable and gets :data:`EMPTY` — all six
    columns stay null, which is the honest row, not a guess.

    ``run`` is published as ``workflow_id`` and the member as
    ``subagent_name``: in a wave the run IS the unit of work whose cost
    somebody asks about, and the member is the thing that ran.
    """
    env = os.environ if environ is None else environ
    run = clean_id(env.get(ENV_RUN) or "") or ""
    label = str(env.get(ENV_LABEL) or "").strip()
    member = ""
    if label:
        # "<run>/<member>"; a label without a separator is the member alone.
        head, _, tail = label.partition("/")
        member = clean_id(tail) if tail else clean_id(head) or ""
        if not run and tail:
            run = clean_id(head) or ""
    if not (run or member):
        return EMPTY
    return Attribution(workflow_id=run or None, subagent_name=member or None)
