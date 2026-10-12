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
from dataclasses import dataclass, field
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

    @property
    def key(self) -> tuple:
        return (self.line, self.col, self.kind)


@dataclass(frozen=True)
class GateCall:
    """One gate call found inside a method."""

    line: int
    name: str


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
        )

    # subprocess.<name>(...) -- a sink only when a path arg is present.
    if isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name) \
            and func.value.id == "subprocess" and func.attr in SUBPROCESS_SINK_METHODS:
        if any(_path_like(a) for a in call.args):
            return Sink(
                line=call.lineno, col=call.col_offset,
                kind=f"subprocess.{func.attr}", receiver="subprocess",
                arg=_first_arg_text(call),
            )
        return None

    # obj.<name>(...) where name is a file method (p.read_text, f.write, ...)
    if isinstance(func, ast.Attribute) and func.attr in FILE_SINK_METHODS:
        return Sink(
            line=call.lineno, col=call.col_offset,
            kind=func.attr, receiver=_receiver_text(func.value),
            arg=_first_arg_text(call),
        )

    # bare open(...)
    if isinstance(func, ast.Name) and func.id == "open":
        return Sink(
            line=call.lineno, col=call.col_offset,
            kind="open", receiver=None, arg=_first_arg_text(call),
        )

    return None


def _is_gate_call(call: ast.Call) -> Optional[GateCall]:
    func = call.func
    # self._check_read_access(...) / self._gate_write_path(...)
    if isinstance(func, ast.Attribute) and func.attr in GATE_METHODS:
        return GateCall(line=call.lineno, name=func.attr)
    # Any call whose function name carries a gate prefix.
    name = None
    if isinstance(func, ast.Name):
        name = func.id
    elif isinstance(func, ast.Attribute):
        name = func.attr
    if name and any(name.startswith(p) for p in GATE_PREFIXES):
        return GateCall(line=call.lineno, name=name)
    return None


# --------------------------------------------------------------------------- 
# Per-method scan
# ---------------------------------------------------------------------------


class _MethodScan:
    """Sinks, gates and outbound helper edges inside one function/method body."""

    __slots__ = ("name", "line", "sinks", "gates", "helper_calls",
                 "helper_edges")

    def __init__(self, name: str, line: int):
        self.name = name
        self.line = line
        self.sinks: list[Sink] = []
        self.gates: list[GateCall] = []
        self.helper_calls: list[str] = []      # same-module defs this method calls
        self.helper_edges: list[tuple[str, int]] = []  # (def_name, call_line)


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

    visit_AsyncFunctionDef = visit_FunctionDef

    def visit_Call(self, node):
        owner = self._owner_stack[-1] if self._owner_stack else None
        sink = _detect_sink(node)
        if sink is not None:
            self.scan.sinks.append(
                Sink(line=sink.line, col=sink.col, kind=sink.kind,
                     receiver=sink.receiver, arg=sink.arg, owner=owner))
        gate = _is_gate_call(node)
        if gate is not None:
            self.scan.gates.append(gate)
        # A call to a same-module helper: record the edge (name, call line) so
        # the executor's sinks are attributed and its gates dominate through
        # the helper chain.
        func = node.func
        target = None
        if isinstance(func, ast.Name):
            target = func.id
        elif isinstance(func, ast.Attribute):
            target = func.attr
        if target in self._local_defs:
            self.scan.helper_calls.append(target)
            self.scan.helper_edges.append((target, node.lineno))
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

    # Second pass: build the reverse call graph (child -> parent+call_line) and
    # transitive reachability (child -> set of methods it can reach), so a sink's
    # protection can be followed up the path that reaches it from its executor.
    inbound: dict[str, list[tuple[str, int]]] = {}
    for name, m in methods.items():
        for child, call_line in m.helper_edges:
            inbound.setdefault(child, []).append((name, call_line))
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
            for child, _line in cur.helper_edges:
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
    inbound: dict[str, list[tuple[str, int]]],
    reaches: dict[str, set[str]],
) -> SinkResult:
    """Whether a gate dominates ``sink`` along the caller path from executor.

    Walk outward from the executor following the helper edges that lead to the
    sink's owner.  A gate in the owner before the sink protects it; otherwise a
    gate in any method on the path *before* the edge it uses to descend protects
    it.  The path is the chain of inbound calls back to ``executor``, so a gate
    in an unrelated caller of a shared helper never grants protection here.
    """
    owner = sink.owner or executor
    if owner == executor or owner in reaches.get(executor, ()):
        # Walk the path from the executor down to the owner.
        gate = _dominant_gate_on_path(
            executor, owner, sink.line, methods, reaches, set())
        return SinkResult(sink, gate is not None,
                          gate.line if gate else None,
                          gate.name if gate else None)
    # Owner not reachable from this executor; the sink was attributed by a
    # mis-parse.  Fall back to the owner's own gates only.
    own = methods.get(owner)
    earlier = [g for g in (own.gates if own else []) if g.line < sink.line]
    gate = min(earlier, key=lambda g: g.line) if earlier else None
    return SinkResult(sink, gate is not None,
                      gate.line if gate else None, gate.name if gate else None)


def _dominant_gate_on_path(
    here: str, target: str, sink_line: int,
    methods: dict[str, _MethodScan], reaches: dict[str, set[str]],
    path: set[str],
) -> Optional[GateCall]:
    """Walk the call-path from ``here`` to ``target``; return the protecting gate.

    ``here`` descends into the child whose subtree contains ``target``; if that
    child *is* ``target``, a gate in it before ``sink_line`` protects the sink;
    otherwise a gate in ``here`` before the chosen edge protects it and we stop.
    """
    if here in path:
        return None
    path = path | {here}
    scan = methods.get(here)
    if scan is None:
        return None
    if here == target:
        earlier = [g for g in scan.gates if g.line < sink_line]
        return min(earlier, key=lambda g: g.line) if earlier else None
    for child, call_line in scan.helper_edges:
        if child == target or target in reaches.get(child, ()):
            before = [g for g in scan.gates if g.line < call_line]
            if before:
                return min(before, key=lambda g: g.line)
            return _dominant_gate_on_path(
                child, target, sink_line, methods, reaches, path)
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
