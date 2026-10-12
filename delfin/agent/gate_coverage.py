"""gate_coverage -- static scan for ungated filesystem sinks in tool executors.

The V2-type review holes of DELFIN waves 13/14 shared one shape: a ``_execute_*``
tool executor that opens or writes a file without routing the path through the
read/write gate (``_check_read_access`` / ``_gate_write_path``), or that lets a
``Path`` reach a filesystem call the gate never saw.  Each was found by human
review, late.

This module performs a pure AST scan (it never imports :mod:`api_client`, so it
can run against any source text, including test fixtures) and answers, per
``_execute_*`` method:

* which filesystem sinks it touches, directly or through a named helper,
* whether a gate call dominates each sink (appears earlier in the same method).

A sink is considered *gated* only when a gate call appears at an earlier line in
the same function that holds the sink.  A gate called *after* the open does not
protect it, and neither does a gate in a sibling method that never runs before
the sink.

Executors with an ungated sink are findings unless the executor's name is listed
in :data:`EXEMPTIONS`.  Every exemption entry names its reason and the source
line it authorises; an entry that no longer matches any executor is *stale* and
is itself a failure (the list must not rot).

The public entry points are :func:`analyze` (parse + scan one source string) and
the methods of :class:`CoverageReport`.  Nothing here imports the DELFIN
runtime, touches files, or changes state.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass, field, replace
from typing import Iterable, Optional

# --------------------------------------------------------------------------- 
# Vocabulary of sinks and gates
# ---------------------------------------------------------------------------

#: Method-call names that touch the filesystem: an attribute call ``x.<name>(...)``
#: on any receiver (a Path, an open file handle, a str).  These are the sink
#: carriers the task names -- open, read_text, write_text -- plus the byte and
#: line siblings a careless rewrite introduces without changing the shape.
FILE_SINK_METHODS: frozenset[str] = frozenset({
    "open", "read_text", "read_bytes", "write_text", "write_bytes",
    "read", "readlines", "write", "writelines",
})

#: ``os.<name>(...)`` calls that create, mutate or destroy a path.
OS_SINK_METHODS: frozenset[str] = frozenset({
    "open", "makedirs", "mkdir", "rename", "replace", "remove", "unlink",
})

#: ``subprocess.<name>(...)`` calls, flagged as sinks only when a Path-like
#: argument is present (a bare ``subprocess.run(["ls"])`` has no file sink).
SUBPROCESS_SINK_METHODS: frozenset[str] = frozenset({
    "run", "Popen", "call", "check_call", "check_output", "check_call",
})

#: Gate-call names, in precedence order.  A call to one of these dominates a
#: sink when it appears earlier in the same function.  ``_resolve_in_workspace``
#: is write_file's write gate (it resolves the path against the writable roots
#: and returns an error when the write is out of scope); ``_gate_write_path`` is
#: the richer variant edit_file/multi_edit use; ``_run_permission_gate`` is the
#: policy+callback gate apply_patch and the pinned-executors call (it applies
#: the sandbox + deny-list + calc-confirm + self-mod guard, so a diff's pre/post
#: image reads under it are gated).  All count as gates.
GATE_METHODS: frozenset[str] = frozenset({
    "_check_read_access", "_gate_write_path", "_resolve_in_workspace",
    "_run_permission_gate",
})

#: Name prefixes treated as gate calls as well (so future gates are recognised
#: without an edit here).  An executor that only *calls* a gate by a renamed
#: variant is still protected.
GATE_PREFIXES: tuple[str, ...] = ("_gate_", "_check_read_", "_check_write_")

# --------------------------------------------------------------------------- 
# Data model
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Sink:
    """One filesystem access in an executor or its helper.

    ``owner`` names the function/method whose body contains the call -- the
    executor itself for an inline open, or the helper it delegates to.  It is
    this owner's gate list (plus any gate an ancestor caller runs before
    reaching the sink) that decides whether the sink is protected.
    """

    line: int
    col: int
    kind: str           # "open" | "read_text" | "os.open" | "subprocess" | ...
    receiver: Optional[str]  # the object whose method is called, or None
    arg: Optional[str]  # the first path-ish argument, best-effort; "" when none
    owner: Optional[str] = None  # the enclosing function name
    #: Control-flow identity of the statement carrying this sink, set by
    #: ``_index_statements``: ``preorder`` is the statement's depth-first order
    #: in the method, ``forkpath`` the tuple of ``(id(fork), id(branch_list))``
    #: frames for every conditional fork (if/for/while/try/match) that statically
    #: contains it.  A gate and a sink dominate one another only along
    #: fork-compatible paths (see :func:`_covers`) -- a gate inside ``if``/``else``
    #: does not cover a sink the fork can reach without it, while a gate in an
    #: enclosing body *does* cover a sink nested in a ``with`` (whose body always
    #: runs, so it is not a fork).
    preorder: Optional[int] = None
    forkpath: Optional[tuple] = None
    #: The concrete path expressions this sink touches, as canonical strings
    #: (:meth:`_expr_key`).  A gate covers the sink only when its own path
    #: argument appears in this set (a gate on ``a`` does not cover ``open(b)``).
    path_exprs: tuple[str, ...] = ()

    @property
    def key(self) -> tuple:
        return (self.line, self.col, self.kind)


@dataclass(frozen=True)
class GateCall:
    """One gate call found inside a method.

    ``path_arg`` is the gate's concrete path argument as canonical text, or
    ``None`` when the gate takes the path at an unknown position (e.g. the
    policy gate, whose per-path tier lives elsewhere) -- in that case the gate
    is treated as covering any path (conservative, keeps same-block behaviour
    for gates that hand the whole arguments dict to the gate chain).
    """

    line: int
    name: str
    #: Control-flow identity of the gate's statement, exactly as on :class:`Sink`
    #: -- a gate only dominates a sink on fork-compatible paths, before it.  A
    #: gate inside ``if force:`` is not on the path of a sink that follows the
    #: ``if``, so it does not cover it.
    preorder: Optional[int] = None
    forkpath: Optional[tuple] = None
    #: Canonical text of the path argument the gate checks, or None = permissive.
    path_arg: Optional[str] = None


@dataclass(frozen=True)
class SinkResult:
    """A sink resolved against the gates of its own enclosing method."""

    sink: Sink
    gated: bool
    gate_line: Optional[int]
    gate_name: Optional[str]


@dataclass
class ExecutorFinding:
    """Everything the scan learned about one ``_execute_*`` method."""

    name: str
    line: int
    sinks: list[Sink] = field(default_factory=list)
    gates: list[GateCall] = field(default_factory=list)
    results: list[SinkResult] = field(default_factory=list)

    @property
    def ungated_sinks(self) -> list[SinkResult]:
        return [r for r in self.results if not r.gated]


@dataclass
class CoverageReport:
    """The outcome of :func:`analyze` over one source text."""

    module_path: str
    executors: list[ExecutorFinding] = field(default_factory=list)

    def names(self) -> list[str]:
        return [e.name for e in self.executors]

    def ungated(self) -> list[ExecutorFinding]:
        """Executors carrying at least one ungated sink."""
        return [e for e in self.executors if e.ungated_sinks]

    def stale_exemptions(
        self, exemptions: Optional[dict[str, str]] = None
    ) -> list[str]:
        """Exemption entries that name no executor in this source (they rot)."""
        live = set(self.names())
        table = EXEMPTIONS if exemptions is None else exemptions
        return [name for name in table if name not in live]


#: Executors whose ungated sink is intentional, with the reason and the source
#: line it authorises.  Every key must name an executor that still exists --
#: :meth:`CoverageReport.stale_exemptions` flags entries that no longer match
#: any code, so the list cannot silently rot as executors are renamed or removed.
#:
#: Kept empty on purpose for the initial build: any ungated sink the scanner
#: finds on today's main belongs in ``.gate/h1d_findings.md`` for the operator,
#: not in an exemption written by the builder.  Add an entry here only with a
#: concrete reason and line, reviewed like any finding.
EXEMPTIONS: dict[str, str] = {}


# --------------------------------------------------------------------------- 
# Sink detection
# ---------------------------------------------------------------------------


def _path_like(arg: ast.expr) -> bool:
    """Best-effort: is this expression a path? A Path(...) call, an f-string
    holding a path variable, or any identifier/attribute (assume string paths
    flow through string variables)."""
    if isinstance(arg, ast.Call) and isinstance(arg.func, ast.Name) \
            and arg.func.id == "Path":
        return True
    if isinstance(arg, ast.JoinedStr):
        return True
    if isinstance(arg, (ast.Name, ast.Attribute, ast.Constant)):
        return True
    return False


def _expr_key(node: ast.expr) -> str:
    """A canonical string for comparing two path expressions.

    ``arguments['cwd']`` and ``arguments['cwd']`` compare equal; the constant
    ``'a'`` and the name ``a`` do not collapse together.  Falls back to
    ``ast.unparse`` so any expression is comparable.
    """
    try:
        return ast.unparse(node).strip()
    except Exception:  # pragma: no cover - unparse exists on Py>=3.9
        return getattr(node, "id", "") or repr(node)[:40]


def _gate_path_slot(name: str) -> Optional[int]:
    """The argument position holding the path for a known gate, else None.

    ``None`` means the gate's path position is unknown (the policy gate hands
    the whole arguments dict through the tier chain) -- callers treat that as
    permissive: it covers any same-block sink.
    """
    if name in ("_check_read_access", "_gate_write_path"):
        # _check_read_access(perms, path, ...) ; _gate_write_path(path, ...)
        # Both are documented with the path in a fixed, known slot.
        return 1 if name == "_check_read_access" else 0
    if name in ("_resolve_in_workspace",):
        # _resolve_in_workspace(path_or_dst, ...) -- the path is the first arg.
        return 0
    # A _gate_*/_check_read_*/_check_write_* variant we do not special-case:
    # path position unknown -> permissive.
    return None


def _sink_path_exprs(sink_kind: str, call: ast.Call) -> tuple[str, ...]:
    """The concrete path expressions a sink touches, as canonical text.

    ``open(path)`` / ``os.open(path)`` use the first positional arg;
    ``subprocess.*`` collects any Path-like positional arg *and* the ``cwd=`` /
    ``executable=`` keyword values (a path arriving as a keyword is still a
    filesystem sink).  The file-method sinks (``f.read()``, ``path.write_text()``)
    act on a receiver that is usually a *handle derived from a gated open* or a
    freshly built Path; their path identity is not reliably comparable to the
    gate's argument, so their coverage is judged by the branch/order rules alone
    (returning ``()`` = undeterminable = permissive).  The gate-relevant sinks
    -- the ones that *open or create* a path -- keep a strict path match.
    """
    if sink_kind in ("open", "os.open"):
        if call.args:
            return (_expr_key(call.args[0]),)
    if sink_kind.startswith("subprocess."):
        exprs = [_expr_key(a) for a in call.args if _path_like(a)]
        for kw in call.keywords:
            if kw.arg in ("cwd", "executable") and kw.value is not None:
                exprs.append(_expr_key(kw.value))
        return tuple(exprs)
    # file-method sink: permissive (see docstring).
    return ()


def _first_arg_text(call: ast.Call) -> Optional[str]:
    if not call.args:
        return None
    first = call.args[0]
    if isinstance(first, ast.Constant) and isinstance(first.value, str):
        return first.value
    if isinstance(first, ast.Name):
        return first.id
    if isinstance(first, ast.Call) and isinstance(first.func, ast.Name):
        return first.func.id + "(...)"
    return ast.unparse(first)[:40] if ast.unparse else None


def _receiver_text(node: ast.expr) -> Optional[str]:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return ast.unparse(node) if ast.unparse else None
    if isinstance(node, ast.Call):
        return ast.unparse(node.func) if ast.unparse else None
    return None


def _is_temp_source_call(node) -> bool:
    """Is this expression a tempfile handle-creation source?

    ``tempfile.mkstemp / NamedTemporaryFile / TemporaryFile`` and
    ``os.fdopen`` (of a mkstemp fd) produce transient staging handles whose
    reads/writes are not persistent user-data accesses; a sink reached through
    them is excluded, not an ungated write.
    """
    if isinstance(node, ast.Call):
        f = node.func
        if isinstance(f, ast.Attribute):
            mod = f.value
            if isinstance(mod, ast.Name) and mod.id == "tempfile":
                return f.attr in ("mkstemp", "NamedTemporaryFile",
                                  "TemporaryFile", "TemporaryDirectory")
            if isinstance(mod, ast.Attribute) and isinstance(mod.value, ast.Name) \
                    and mod.value.id == "_tf":
                return f.attr in ("mkstemp", "NamedTemporaryFile", "TemporaryFile")
            if isinstance(mod, ast.Name) and mod.id == "os" \
                    and f.attr == "fdopen":
                return True
    return False


def _detect_sink(call: ast.Call) -> Optional[Sink]:
    """Return a Sink if ``call`` touches the filesystem, else None."""
    func = call.func

    # os.<name>(...) -- os.open, os.makedirs, ...
    if isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name) \
            and func.value.id == "os" and func.attr in OS_SINK_METHODS:
        return Sink(
            line=call.lineno, col=call.col_offset,
            kind=f"os.{func.attr}", receiver="os",
            arg=_first_arg_text(call),
            path_exprs=_sink_path_exprs("os.open", call),
        )

    # subprocess.<name>(...) -- a sink only when a path arg is present.
    if isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name) \
            and func.value.id == "subprocess" and func.attr in SUBPROCESS_SINK_METHODS:
        paths = _sink_path_exprs(f"subprocess.{func.attr}", call)
        if paths:
            return Sink(
                line=call.lineno, col=call.col_offset,
                kind=f"subprocess.{func.attr}", receiver="subprocess",
                arg=_first_arg_text(call),
                path_exprs=paths,
            )
        return None

    # obj.<name>(...) where name is a file method (p.read_text, f.write, ...)
    if isinstance(func, ast.Attribute) and func.attr in FILE_SINK_METHODS:
        return Sink(
            line=call.lineno, col=call.col_offset,
            kind=func.attr, receiver=_receiver_text(func.value),
            arg=_first_arg_text(call),
            path_exprs=_sink_path_exprs("file-method", call),
        )

    # bare open(...)
    if isinstance(func, ast.Name) and func.id == "open":
        return Sink(
            line=call.lineno, col=call.col_offset,
            kind="open", receiver=None, arg=_first_arg_text(call),
            path_exprs=_sink_path_exprs("open", call),
        )

    return None


def _gate_path_arg(call: ast.Call, name: str) -> Optional[str]:
    """The canonical text of the gate's path argument, or None if permissive."""
    slot = _gate_path_slot(name)
    if slot is None:
        return None
    args = call.args
    if slot < len(args):
        return _expr_key(args[slot])
    # A keyword could carry the path; check the well-known ones.
    for kw in call.keywords:
        if kw.arg in ("path", "path_arg", "rel", "full", "dst"):
            return _expr_key(kw.value)
    return None


def _is_gate_call(call: ast.Call) -> Optional[GateCall]:
    func = call.func
    # self._check_read_access(...) / self._gate_write_path(...)
    if isinstance(func, ast.Attribute) and func.attr in GATE_METHODS:
        return GateCall(line=call.lineno, name=func.attr,
                        path_arg=_gate_path_arg(call, func.attr))
    # Any call whose function name carries a gate prefix.
    name = None
    if isinstance(func, ast.Name):
        name = func.id
    elif isinstance(func, ast.Attribute):
        name = func.attr
    if name and any(name.startswith(p) for p in GATE_PREFIXES):
        return GateCall(line=call.lineno, name=name,
                        path_arg=_gate_path_arg(call, name))
    return None


# --------------------------------------------------------------------------- 
# Per-method scan
# ---------------------------------------------------------------------------


# Conditional forks: a statement inside one of these may be skipped or choose a
# different branch, so a gate in one branch never covers a sink the fork can
# reach without it.  ``with`` is deliberately NOT a fork: its body always runs.
_FORKS: tuple = (ast.If, ast.For, ast.While, ast.Try, ast.Match)


def _index_statements(tree: ast.AST) -> None:
    """Attach to every statement its control-flow facts for dominance.

    For each statement node this attaches:

    * ``.delfin_parent`` -- its AST parent (set on every node),
    * ``.delfin_preorder`` -- a depth-first order number within the module,
    * ``.delfin_forkpath`` -- the tuple of ``(id(fork), id(branch_list))``
      frames for every conditional fork that statically contains it, outer
      first.

    A gate covers a sink only on fork-compatible paths (same frames for all
    forks both descend through, or the gate sits in an enclosing body and
    precedes it).  ``with`` bodies contribute no frame, so a gate in the
    function body still covers a sink inside ``with open(...) as f:``.
    """
    for parent_node in ast.walk(tree):
        for child in ast.iter_child_nodes(parent_node):
            child.delfin_parent = parent_node

    counter = 0
    forks: list[tuple[int, int]] = []

    def step(stmt: ast.stmt) -> None:
        nonlocal counter
        counter += 1
        stmt.delfin_preorder = counter
        stmt.delfin_forkpath = tuple(forks)

    def walk_list(body_list: object) -> None:
        if not isinstance(body_list, list):
            return
        for stmt in body_list:
            if not isinstance(stmt, ast.stmt):
                continue
            step(stmt)
            descend(stmt)

    def descend(stmt: ast.stmt) -> None:
        # Plain containers (with/class/module): their body always runs -> no frame.
        if isinstance(stmt, ast.With):
            walk_list(stmt.body)
            return
        # Conditional forks: push a frame per branch, then recurse.
        if isinstance(stmt, ast.If):
            forks.append((id(stmt), id(stmt.body)))
            walk_list(stmt.body)
            forks.pop()
            forks.append((id(stmt), id(stmt.orelse)))
            walk_list(stmt.orelse)
            forks.pop()
            return
        if isinstance(stmt, (ast.For, ast.While)):
            forks.append((id(stmt), id(stmt.body)))
            walk_list(stmt.body)
            forks.pop()
            forks.append((id(stmt), id(stmt.orelse)))
            walk_list(stmt.orelse)
            forks.pop()
            return
        if isinstance(stmt, ast.Try):
            for bl in (stmt.body, stmt.orelse, stmt.finalbody):
                forks.append((id(stmt), id(bl)))
                walk_list(bl)
                forks.pop()
            for handler in getattr(stmt, "handlers", None) or ():
                forks.append((id(stmt), id(handler.body)))
                walk_list(handler.body)
                forks.pop()
            return
        if isinstance(stmt, ast.Match):
            for case in getattr(stmt, "cases", None) or ():
                forks.append((id(stmt), id(case.body)))
                walk_list(case.body)
                forks.pop()
            return
        # Any other compound statement (function/class) body: always run bodies
        # are walked without a frame; nested statements inside lists are handled
        # by walking every statement-list attribute recursively.
        for node in ast.iter_child_nodes(stmt):
            if isinstance(node, ast.stmt):
                step(node)
                descend(node)
            elif isinstance(node, ast.excepthandler):
                for body_st in (getattr(node, "body", None) or ()):
                    if isinstance(body_st, ast.stmt):
                        step(body_st)
                        descend(body_st)

    walk_list(tree.body)


def _stmt_ctx(node: ast.AST) -> tuple[Optional[int], Optional[tuple]]:
    """The (preorder, forkpath) of the statement enclosing ``node``."""
    cur = node
    while cur is not None and not isinstance(cur, ast.stmt):
        cur = getattr(cur, "delfin_parent", None)
    if cur is None:
        return (None, None)
    return (getattr(cur, "delfin_preorder", None),
            getattr(cur, "delfin_forkpath", None))


def _fork_compatible(fg: Optional[tuple], fs: Optional[tuple]) -> bool:
    """Is the gate's fork-path compatible with the sink's, i.e. is the gate on a
    branch of every fork the sink can flow through?

    Where the gate and the sink diverge at a *shared* fork frame (both sit
    inside the same ``if``/``try``/loop but in different bodies -- the gate in
    ``if body``, the sink in ``else``), the sink can run on a path the gate is
    not on, so the gate never covers it: incompatible.  By contrast a gate that
    is *deeper* inside a conditional with the sink after that whole conditional
    is the DELFIN *guard* idiom -- ``if perms is not None: check_read(...)``
    followed by the read -- and is judged by the path argument instead, so it is
    compatible here.  ``None`` on either side (no CF facts, e.g. a synthesised
    edge sink) is treated as compatible so the pre-order test below decides.
    """
    if fg is None or fs is None:
        return True
    i = 0
    fg, fs = tuple(fg), tuple(fs)
    while i < len(fg) and i < len(fs) and fg[i] == fs[i]:
        i += 1
    if i < len(fg) and i < len(fs):
        # Both still have frames at the same depth but from different bodies:
        # a gate in one branch and a sink in a sibling branch of the same fork.
        # The sink can run on a path the gate's branch does not take -- reject.
        return False
    # One side exhausted: the gate is deeper inside a conditional than the
    # sink (sink after the fork) -- the guard idiom, path-matching decides.
    return True


def _resolve_receiver(receiver: Optional[str], bindings) -> Optional[str]:
    """Follow a file-method sink's receiver to the path it was opened/bound from.

    ``with open(x) as f: f.read()`` resolves ``f`` to the path ``x`` via the
    method's binding table, so a gate on ``x`` covers the read.  Returns None
    when the receiver is not a plain bound name (permissive).
    """
    if receiver is None:
        return None
    if not isinstance(receiver, str):
        return None
    if receiver not in bindings:
        return None
    return bindings[receiver]


def _covers(gate: GateCall, sink: Sink, bindings) -> bool:
    """Does ``gate`` dominate ``sink``?

    Every requirement must hold:

    * the gate's fork-path is compatible with the sink's (see
      :func:`_fork_compatible`) -- a gate inside ``if force:`` is not on the
      path of a sink that follows the ``if``, a gate in one branch does not
      cover an ``open`` in the sibling ``else``, yet a gate in the body *does*
      cover a sink nested in a ``with`` (whose body always runs), and
    * the gate statement precedes the sink statement (pre-order), and
    * the gate's path argument covers a path the sink touches -- directly for
      ``open``/``os.*``/``subprocess`` sinks, through handle binding for a file
      method (``f.read()`` where ``with open(x) as f``), else permissive.  A
      gate on ``a`` never covers ``open(b)``.
    """
    if not _fork_compatible(gate.forkpath, sink.forkpath):
        return False
    if gate.preorder is None or sink.preorder is None:
        return False
    if gate.preorder >= sink.preorder:
        return False
    # A gate that sits *deeper* inside a conditional than the sink (the DELFIN
    # guard idiom, ``if perms is not None: check_read(...)`` then the read) only
    # protects the sink when it demonstrably cites the same path -- the
    # permission check is a real gate, not a license to open arbitrary files.
    guard_idiom = bool(gate.forkpath) and len(gate.forkpath) > len(
        sink.forkpath or ())
    if gate.path_arg is None:
        return not guard_idiom  # permissive elsewhere; not for the guard idiom
    if sink.path_exprs:
        # ``open(path)`` / ``os.*`` / ``subprocess`` sinks carry concrete paths;
        # a gate on a different path does not cover them.
        return any(gate.path_arg == e for e in sink.path_exprs)
    # File-method sink (receiver-based, no concrete path): resolve the receiver
    # to the path it was opened from so ``with open(x) as f: f.read()`` is
    # covered by a gate on ``x``; an unresolvable handle is permissive *except*
    # under the guard idiom, where a gate must cite the path it allows.
    resolved = _resolve_receiver(sink.receiver, bindings)
    if resolved is None:
        return not guard_idiom  # cannot prove a mismatch (elsewhere)
    return gate.path_arg == resolved


def _any_covering_gate(
    gates: Iterable[GateCall], sink: Sink, bindings
) -> Optional[GateCall]:
    """The first gate (by statement order) covering ``sink``, or None."""
    cands = [g for g in gates if _covers(g, sink, bindings)]
    if not cands:
        return None
    return min(cands, key=lambda g: (g.preorder, g.line))


class _MethodScan:
    """Sinks, gates and outbound helper edges inside one function/method body."""

    __slots__ = ("name", "line", "sinks", "gates", "helper_calls",
                   "helper_edges", "bindings", "param_names", "temp_sources")

    def __init__(self, name: str, line: int):
        self.name = name
        self.line = line
        self.sinks: list[Sink] = []
        self.gates: list[GateCall] = []
        self.helper_calls: list[str] = []  # same-module defs this method calls
        # (def_name, call_line, call_stmt_preorder, call_stmt_forkpath,
        #  call_pos_arg_texts) -- the last is the canonical texts of the call's
        #  positional path-ish arguments, so a callee param can be bound back to
        #  the caller's (possibly gated) expression when a sink is reached here.
        self.helper_edges: list[
            tuple[str, int, Optional[int], Optional[tuple], tuple[str, ...]]] = []
        # variable name -> canonical path expression it binds (open handles)
        self.bindings: dict[str, str] = {}
        # variable/handle names bound to a tempfile handle (mkstemp /
        # NamedTemporaryFile wrapped through os.fdopen / with): reads/writes on
        # them are transient staging, not ungated user-data writes.
        self.temp_sources: set[str] = set()
        # positional parameter names, for interprocedural argument->param binding
        self.param_names: list[str] = []


class _BodyVisitor(ast.NodeVisitor):
    """Collects sinks, gates and same-module helper-call edges in one body."""

    def __init__(self) -> None:
        self.scan = _MethodScan("_", 0)
        self._local_defs: set[str] = set()
        self._owner_stack: list[str] = []  # innermost enclosing def name

    def visit_FunctionDef(self, node):
        self._owner_stack.append(node.name)
        self.generic_visit(node)
        self._owner_stack.pop()

    def visit_arguments(self, node):
        # Record the positional parameter names of the innermost enclosing def so
        # a gate's dominance walk can bind a helper's parameters back to the
        # caller's argument expressions at a call edge.  kwonlyargs/vararg are
        # not path-carrying positions, so only posonlyargs + args are kept.
        names = [a.arg for a in getattr(node, "posonlyargs", ())]
        names += [a.arg for a in getattr(node, "args", ())]
        self.scan.param_names = names

    visit_AsyncFunctionDef = visit_FunctionDef

    def visit_Call(self, node):
        owner = self._owner_stack[-1] if self._owner_stack else None
        pre, fpath = _stmt_ctx(node)
        sink = _detect_sink(node)
        if sink is not None:
            # A sink on a tempfile staging handle/name (mkstemp bound through
            # os.fdopen/with / NamedTemporaryFile) is transient, not an ungated
            # user-data write -- exclude it (its receiver handle OR any path
            # argument that is such a temp source).
            if sink.receiver is not None and sink.receiver in self.scan.temp_sources:
                sink = None
            elif sink.path_exprs and any(
                    e in self.scan.temp_sources for e in sink.path_exprs):
                sink = None
        if sink is not None:
            self.scan.sinks.append(replace(
                sink, owner=owner, preorder=pre, forkpath=fpath))
        gate = _is_gate_call(node)
        if gate is not None:
            self.scan.gates.append(replace(gate, preorder=pre, forkpath=fpath))
        # A call to a same-module helper: record the edge (name, call line and
        # the call statement's own block) so the executor's sinks are
        # attributed and its gates dominate through the helper chain.
        func = node.func
        target = None
        if isinstance(func, ast.Name):
            target = func.id
        elif isinstance(func, ast.Attribute):
            target = func.attr
        if target in self._local_defs:
            self.scan.helper_calls.append(target)
            arg_texts = tuple(_expr_key(a) for a in node.args[:3])
            self.scan.helper_edges.append(
                (target, node.lineno, pre, fpath, arg_texts))
        self.generic_visit(node)

    def visit_With(self, node):
        # Bind the exported handle name to the opened path expression so a later
        # file-method sink on that handle resolves to the underlying path.  For
        # ``with open(p) as f`` the path is the open's first argument, not the
        # whole ``open(p)`` call.  A tempfile-context handle (``NamedTemporaryFile``
        # / ``os.fdopen`` of a mkstemp fd) is a transient staging handle.
        for item in node.items:
            if item.optional_vars is not None and isinstance(
                    item.optional_vars, ast.Name):
                ctx = item.context_expr
                if _is_temp_source_call(ctx):
                    self.scan.temp_sources.add(item.optional_vars.id)
                elif (isinstance(ctx, ast.Call) and isinstance(ctx.func, ast.Name)
                        and ctx.func.id == "open" and ctx.args):
                    self.scan.bindings[item.optional_vars.id] = _expr_key(
                        ctx.args[0])
                else:
                    self.scan.bindings[item.optional_vars.id] = _expr_key(ctx)
        self.generic_visit(node)

    def visit_Assign(self, node):
        # `var = open(path)` (or a Path/result the gate checks) makes a handle.
        # `fd, name = tempfile.mkstemp(...)` makes a tempfile handle AND name;
        # both are transient staging, marked so their later sinks are excluded.
        if len(node.targets) == 1 and isinstance(node.value, ast.Call):
            val = node.value
            if _is_temp_source_call(val):
                t = node.targets[0]
                if isinstance(t, ast.Name):
                    self.scan.temp_sources.add(t.id)
                elif isinstance(t, ast.Tuple):
                    for elt in t.elts:
                        if isinstance(elt, ast.Name):
                            self.scan.temp_sources.add(elt.id)
            elif isinstance(node.targets[0], ast.Name):
                recv = _receiver_text(val.func) if isinstance(
                    val.func, ast.Attribute) else None
                if recv == "open":
                    self.scan.bindings[node.targets[0].id] = _expr_key(
                        val.args[0]) if val.args else ""
        self.generic_visit(node)


# --------------------------------------------------------------------------- 
# Public entry point
# ---------------------------------------------------------------------------


def analyze(source: str, module_path: str = "") -> CoverageReport:
    """Parse ``source`` and scan every ``_execute_*`` method for file sinks.

    Returns a :class:`CoverageReport`.  Raises ``SyntaxError`` on unparsable
    input.  Does not import the DELFIN runtime and touches nothing on disk.
    """
    tree = ast.parse(source)
    report = CoverageReport(module_path=module_path)
    _index_statements(tree)

    # Every function/method name in the module, used to recognise a call to a
    # same-module helper (for example ``self._read_all(...)``, a sibling
    # method, or a module-level ``read_text_file(...)``) so its sinks and gates
    # are attributed back to the executor that names it.
    all_defs = {
        d.name for d in ast.walk(tree)
        if isinstance(d, (ast.FunctionDef, ast.AsyncFunctionDef))
    }

    # First pass: map every method by name -> its scan, cloned per definition.
    methods: dict[str, _MethodScan] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            vs = _BodyVisitor()
            vs._local_defs = all_defs
            vs.scan = _MethodScan(node.name, node.lineno)
            vs.visit(node)
            methods[node.name] = vs.scan

    # Second pass: build the reverse call graph (child -> parent+call_stmt) and
    # transitive reachability (child -> set of methods it can reach), so a sink's
    # protection can be followed up the path that reaches it from its executor.
    inbound: dict[str, list[tuple[str, int, Optional[int], Optional[int]]]] = {}
    for name, m in methods.items():
        for child, call_line, blk, idx, _args in m.helper_edges:
            inbound.setdefault(child, []).append((name, call_line, blk, idx))
    reaches = _reachability(methods)

    # Third pass: for each _execute_* method, attribute the sinks of every
    # helper it can reach, then decide each sink's gating by walking the call
    # path from this executor to the method that owns the sink.
    for name, scan in methods.items():
        if not name.startswith("_execute_"):
            continue
        finding = ExecutorFinding(name=name, line=scan.line)
        seen: set = set()
        stack = [scan]
        while stack:
            cur = stack.pop()
            if id(cur) in seen:
                continue
            seen.add(id(cur))
            finding.sinks.extend(cur.sinks)
            finding.gates.extend(cur.gates)
            for child, _line, _blk, _idx, _args in cur.helper_edges:
                hscan = methods.get(child)
                if hscan is not None:
                    stack.append(hscan)
        finding.sinks.sort(key=lambda s: (s.line, s.col))
        finding.results = [
            _dominated_result(s, name, methods, inbound, reaches)
            for s in finding.sinks
        ]
        report.executors.append(finding)

    report.executors.sort(key=lambda e: e.line)
    return report


def _reachability(methods: dict[str, _MethodScan]) -> dict[str, set[str]]:
    """For each method name, the set of method names reachable via helper edges."""
    out: dict[str, set[str]] = {}
    for name, m in methods.items():
        seen: set[str] = set()
        stack = list(m.helper_calls)
        while stack:
            cur = stack.pop()
            if cur in seen or cur not in methods:
                continue
            seen.add(cur)
            stack.extend(methods[cur].helper_calls)
        out[name] = seen
    return out


def _dominated_result(
    sink: Sink,
    executor: str,
    methods: dict[str, _MethodScan],
    inbound: dict[str, list[tuple[str, int, Optional[int], Optional[int]]]],
    reaches: dict[str, set[str]],
) -> SinkResult:
    """Whether a gate dominates ``sink`` along the caller path from executor.

    Walk outward from the executor following the helper edges that lead to the
    sink's owner.  A gate in the owner covering the sink protects it; otherwise
    a gate in any method on the path *before* the edge it uses to descend
    protects it.  The path is the chain of inbound calls back to ``executor``,
    so a gate in an unrelated caller of a shared helper never grants protection
    here passes a gate that does not dominate in control-flow terms.

    Coverage is judged by :func:`_covers` -- same statement block, gate before
    the sink, and a matching path argument -- never by line number alone.
    """
    owner = sink.owner or executor
    if owner == executor:
        return _sink_result_with(sink, _any_covering_gate(
            methods[owner].gates, sink, methods[owner].bindings))
    if owner in reaches.get(executor, ()):
        # Walk the path from the executor down to the owner.
        gate = _dominant_gate_on_path(
            executor, owner, sink, methods, reaches, set())
        return _sink_result_with(sink, gate)
    # Owner not reachable from this executor; the sink was attributed by a
    # mis-parse.  Fall back to the owner's own gates only.
    own = methods.get(owner)
    gate = _any_covering_gate(
        own.gates if own else [], sink, own.bindings if own else {})
    return _sink_result_with(sink, gate)


def _sink_result_with(sink: Sink, gate: Optional[GateCall]) -> SinkResult:
    if gate is None:
        return SinkResult(sink, False, None, None)
    return SinkResult(sink, True, gate.line, gate.name)


def _subst_expr(expr: str, sub: dict[str, str]) -> str:
    """Rewrite a canonical path expression through one param->arg binding.

    A bare identifier matching a parameter name is replaced by its bound
    argument text; anything else is returned unchanged.  The binding maps a
    callee parameter name to the caller's (already caller-scoped) expression.
    """
    if expr.isidentifier() and expr in sub:
        return sub[expr]
    return expr


def _dominant_gate_on_path(
    here: str, target: str, sink: Sink,
    methods: dict[str, _MethodScan], reaches: dict[str, set[str]],
    path: set[str],
) -> Optional[GateCall]:
    """Walk the call-path from ``here`` to ``target``; return the protecting gate.

    ``here`` descends into the child whose subtree contains ``target``; if that
    child *is* ``target``, a gate in it covering the sink protects it; otherwise
    a gate in ``here`` covering the *call statement* that descends protects it
    and we stop.  Coverage always goes through :func:`_covers`, so a gate in a
    sibling branch, after its sink, or on another path never protects.

    Interprocedural argument binding: when a helper's parameter is the path a
    sink inside it uses, that parameter is bound back to the caller's argument
    text at the descending edge, so a gate in the caller on that argument is
    credited for the deep sink (``def h(p): open(p)`` called as ``h(gated)``).
    """
    if here in path:
        return None
    path = path | {here}
    scan = methods.get(here)
    if scan is None:
        return None
    if here == target:
        return _any_covering_gate(scan.gates, sink, scan.bindings)
    # Gate in `here` that covers the descending call statement, before it.
    for child, _line, pre, fpath, arg_texts in scan.helper_edges:
        if child == target or target in reaches.get(child, ()):
            child_scan = methods.get(child)
            child_params = child_scan.param_names if child_scan else []
            # Bind the child's path-carrying params to this call's arguments, so
            # the edge_sink's path is expressed in ``here`` scope for the gate.
            sub = {child_params[k]: arg_texts[k]
                   for k in range(min(len(child_params), len(arg_texts)))}
            # The gate in ``here`` must compare against the sink's path expressed
            # in ``here`` scope.  Concrete open/os/subprocess paths substitute
            # directly; a file-method sink (empty paths, e.g. ``fh.read()``)
            # resolves its receiver through the *owner's* bindings (``fh`` ->
            # the path it was opened from, in the owner's scope) first, then that
            # owner-scope path substitutes through the edge to caller scope.
            if sink.path_exprs:
                src_paths = sink.path_exprs
            elif sink.receiver is not None:
                owner_scan = methods.get(target)
                owner_bind = owner_scan.bindings if owner_scan else {}
                res = _resolve_receiver(sink.receiver, owner_bind)
                src_paths = (res,) if res else ()
            else:
                src_paths = ()
            bound_paths = tuple(
                _subst_expr(e, sub) for e in src_paths) if sub else src_paths
            edge_sink = replace(
                sink, preorder=pre, forkpath=fpath, path_exprs=bound_paths)
            before = _any_covering_gate(scan.gates, edge_sink, scan.bindings)
            if before:
                return before
            return _dominant_gate_on_path(
                child, target, sink, methods, reaches, path)
    return None


# --------------------------------------------------------------------------- 
# Convenience helper for callers / tests
# ---------------------------------------------------------------------------


def format_report(report: CoverageReport) -> str:
    """A compact human-readable rendering, line per ungated sink."""
    if not report.ungated():
        return f"{report.module_path or '<source>'}: no ungated file sinks."
    lines = [f"{report.module_path or '<source>'}: {len(report.ungated())} "
             f"executor(s) with ungated file sink(s):"]
    for ex in report.ungated():
        for r in ex.ungated_sinks:
            s = r.sink
            where = ex.name
            lines.append(
                f"  {where}@{s.line}:{s.col}  {s.kind}"
                f"{f'({s.receiver})' if s.receiver else ''}"
                f"{f'  arg={s.arg!r}' if s.arg else ''}  UNGATED"
            )
    return "\n".join(lines)
