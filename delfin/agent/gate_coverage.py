"""gate_coverage -- static scan for ungated filesystem sinks in tool executors.

Input: the source text of a module (``api_client.py`` in production, fixture
text in tests).  Output: a :class:`CoverageReport` that lists, for every
``_execute_*`` method, the filesystem sinks it reaches (in its own body or in
a same-module helper it calls) and whether a gate dominates each of them.

Pure AST: the module is parsed, never imported, and nothing on disk is read
or written.

A sink is *gated* when a gate call

* runs before it (statement pre-order),
* on every path that reaches it: a gate inside ``if``/``else``/``try``/loop
  does not cover a sink the fork can reach without it.  The one conditional
  that is accepted is ``if perms is not None:`` (``permissions`` alike): with
  no permissions object there is no sandbox to consult, and the dispatcher
  refuses such calls for every tool that needs one,
* and names the same path: each path a sink touches must be the path
  argument of a covering gate (``str(x)``/``Path(x)``/``os.fspath(x)`` are
  read as ``x``; the path a gate returns, ``resolved, err =
  self._resolve_in_workspace(p, ...)``, counts as checked).  A gate whose
  path position is unknown (the policy gate, the bash gates, which take the
  whole command) covers every path.

What the scan cannot see: file access inside another module (an image
loader, an office reader, the change journal) unless the call is listed in
:data:`XMOD_SINKS`, and control-flow correlation (``if not check_only:`` around
both the gate and the read).  An executor with no sink is therefore *not
shown to be file-free*; it is only free of the sinks named here.

Rejected: crediting a gate by line number (a gate below the open, or in a
sibling branch, would count); crediting any ``_gate_*``-prefixed call (the
MCP tool gates are not path gates and covered every later open).
"""

from __future__ import annotations

import ast
from dataclasses import dataclass, field, replace
from typing import Iterable, Optional

# ---------------------------------------------------------------------------
# Vocabulary of sinks and gates
# ---------------------------------------------------------------------------

#: Methods on a path object; the receiver is the path (``p.unlink()``).  The
#: second set also takes a second path as the first argument.
PATH_METHODS: frozenset[str] = frozenset({
    "open", "read_text", "read_bytes", "write_text", "write_bytes",
    "unlink", "mkdir", "rmdir", "touch", "chmod", "iterdir", "glob", "rglob",
    "rename", "symlink_to", "hardlink_to",
})
PATH_METHODS_TWO: frozenset[str] = frozenset(
    {"rename", "symlink_to", "hardlink_to"})
#: Methods on an open handle; the path is the one the handle was opened from.
#: (``replace`` is left out of both sets: ``str.replace`` is everywhere.)
HANDLE_METHODS: frozenset[str] = frozenset(
    {"read", "readlines", "write", "writelines"})

#: ``os.<name>`` / ``shutil.<name>`` -> how many leading positional
#: arguments are paths.
MODULE_SINKS: dict[str, dict[str, int]] = {
    "os": {"open": 1, "makedirs": 1, "mkdir": 1, "remove": 1, "unlink": 1,
           "rmdir": 1, "removedirs": 1, "listdir": 1, "scandir": 1, "walk": 1,
           "chmod": 1, "truncate": 1, "rename": 2, "replace": 2,
           "symlink": 2, "link": 2},
    "shutil": {"copy": 2, "copy2": 2, "copyfile": 2, "copytree": 2,
               "copymode": 2, "copystat": 2, "move": 2, "rmtree": 1},
}

#: File access in another module of the package, by import alias: the first
#: argument is the path.  Only the readers/writers the executors call by
#: these aliases today; any other cross-module access is invisible.
XMOD_SINKS: dict[str, frozenset[str]] = {
    "_text_files": frozenset({"read_text_file", "write_text_file"}),
    "_nb": frozenset({"read_cells", "apply_edit"}),
}

#: ``subprocess.<name>(...)``: a sink when any argument is not a constant.
SUBPROCESS_SINK_METHODS: frozenset[str] = frozenset({
    "run", "Popen", "call", "check_call", "check_output",
})

#: Gate name -> position of its path argument, or None when the gate takes
#: the path elsewhere (the whole arguments dict, a command string) and so
#: covers any path.
GATES: dict[str, Optional[int]] = {
    "_check_read_access": 1,
    "_gate_write_path": 0,
    "_resolve_in_workspace": 0,
    "_worktree_path_refusal": 1,
    "_run_permission_gate": None,
    "_gate_bash_read_paths": None,
    "_gate_bash_write_targets": None,
}
#: Gates that return the checked path as the first element of their result.
_GATES_RETURNING_PATH: frozenset[str] = frozenset({"_resolve_in_workspace"})

#: The parameter names of the permissions object (see the module docstring).
_PERMS_NAMES: frozenset[str] = frozenset({"perms", "permissions"})

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Sink:
    """One filesystem access in an executor or in a helper it reaches.

    ``owner`` is the function whose body holds the call.  ``preorder`` and
    ``forkpath`` place the enclosing statement in the control flow (see
    :func:`_index_statements`).  ``path_exprs`` are the canonical texts of the
    paths touched; empty for a handle method whose handle is resolved later
    through the owner's bindings (``receiver``).
    """

    line: int
    col: int
    kind: str
    receiver: Optional[str]
    arg: Optional[str]
    owner: Optional[str] = None
    preorder: Optional[int] = None
    forkpath: Optional[tuple] = None
    path_exprs: tuple[str, ...] = ()


@dataclass(frozen=True)
class GateCall:
    """One gate call.  ``path_args`` is None for a gate that covers any path,
    else the paths it checks (its argument, and the name it returns the
    resolved path into)."""

    line: int
    name: str
    preorder: Optional[int] = None
    forkpath: Optional[tuple] = None
    path_args: Optional[tuple[str, ...]] = None


@dataclass(frozen=True)
class SinkResult:
    sink: Sink
    gated: bool
    gate_line: Optional[int]
    gate_name: Optional[str]


@dataclass
class ExecutorFinding:
    name: str
    line: int
    end_line: int = 0
    sinks: list[Sink] = field(default_factory=list)
    gates: list[GateCall] = field(default_factory=list)
    results: list[SinkResult] = field(default_factory=list)

    @property
    def ungated_sinks(self) -> list[SinkResult]:
        return [r for r in self.results if not r.gated]


@dataclass
class CoverageReport:
    module_path: str
    executors: list[ExecutorFinding] = field(default_factory=list)

    def names(self) -> list[str]:
        return [e.name for e in self.executors]

    def ungated(self) -> list[ExecutorFinding]:
        """Executors carrying at least one ungated sink."""
        return [e for e in self.executors if e.ungated_sinks]

    def stale_exemptions(self, exemptions: Iterable[str]) -> list[str]:
        """Names in ``exemptions`` that are no executor of this source."""
        live = set(self.names())
        return [name for name in exemptions if name not in live]


# ---------------------------------------------------------------------------
# Sink and gate detection
# ---------------------------------------------------------------------------


def _strip_path_wrappers(node: ast.expr) -> ast.expr:
    """``str(x)``, ``Path(x)``, ``pathlib.Path(x)``, ``os.fspath(x)`` -> ``x``."""
    while (isinstance(node, ast.Call) and len(node.args) == 1
           and not node.keywords):
        f = node.func
        name = f.id if isinstance(f, ast.Name) else (
            f.attr if isinstance(f, ast.Attribute) else None)
        if name not in ("str", "Path", "fspath"):
            break
        node = node.args[0]
    return node


def _expr_key(node: ast.expr) -> str:
    """Canonical text of a path expression, wrappers stripped."""
    return ast.unparse(_strip_path_wrappers(node)).strip()


def _short(call: ast.Call) -> Optional[str]:
    return ast.unparse(call.args[0])[:40] if call.args else None


def _is_temp_source_call(node, temp_names: set[str]) -> bool:
    """A ``tempfile`` creation call (``tempfile``/``_tf`` alias), or
    ``os.fdopen`` of a descriptor that one returned."""
    if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)):
        return False
    mod, attr = node.func.value.id, node.func.attr
    if mod in ("tempfile", "_tf"):
        return attr in ("mkstemp", "NamedTemporaryFile", "TemporaryFile")
    if mod == "os" and attr == "fdopen":
        return (bool(node.args) and isinstance(node.args[0], ast.Name)
                and node.args[0].id in temp_names)
    return False


def _from_arguments(node: ast.expr) -> bool:
    """``arguments[...]`` / ``arguments.get(...)``: the model's tool input."""
    if isinstance(node, ast.Subscript):
        node = node.value
    elif (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
            and node.func.attr == "get"):
        node = node.func.value
    return isinstance(node, ast.Name) and node.id == "arguments"


def _argv_path(elt: ast.expr) -> bool:
    """An argv element that is a path: wrapped in ``str``/``Path``/
    ``os.fspath``, or taken from the tool arguments."""
    return _strip_path_wrappers(elt) is not elt or _from_arguments(elt)


def _detect_sink(call: ast.Call) -> Optional[Sink]:
    """Return a Sink if ``call`` touches the filesystem, else None."""
    func = call.func

    def make(kind, receiver, paths):
        return Sink(line=call.lineno, col=call.col_offset, kind=kind,
                    receiver=receiver, arg=_short(call),
                    path_exprs=tuple(paths))

    if isinstance(func, ast.Name):
        if func.id == "open" and call.args:
            return make("open", None, [_expr_key(call.args[0])])
        return None
    if not isinstance(func, ast.Attribute):
        return None
    mod = func.value.id if isinstance(func.value, ast.Name) else None

    if mod in MODULE_SINKS and func.attr in MODULE_SINKS[mod]:
        n = MODULE_SINKS[mod][func.attr]
        return make(f"{mod}.{func.attr}", mod,
                    [_expr_key(a) for a in call.args[:n]])
    if mod in XMOD_SINKS and func.attr in XMOD_SINKS[mod] and call.args:
        return make(f"{mod}.{func.attr}", mod, [_expr_key(call.args[0])])
    if mod == "subprocess":
        if func.attr not in SUBPROCESS_SINK_METHODS:
            return None
        # A whole command held in a variable is a path carrier; inside an
        # argv list only the elements that are visibly paths or model input
        # are, so a git ref or a flag variable is not.
        paths: list[str] = []
        for a in call.args:
            if isinstance(a, (ast.List, ast.Tuple)):
                paths += [_expr_key(e) for e in a.elts if _argv_path(e)]
            elif not isinstance(a, ast.Constant):
                paths.append(_expr_key(a))
        paths += [_expr_key(kw.value) for kw in call.keywords
                  if kw.arg in ("cwd", "executable")
                  and not isinstance(kw.value, ast.Constant)]
        return make(f"subprocess.{func.attr}", "subprocess", paths) \
            if paths else None
    if func.attr in PATH_METHODS:
        paths = [_expr_key(func.value)]
        if func.attr in PATH_METHODS_TWO and call.args:
            paths.append(_expr_key(call.args[0]))
        return make(func.attr, ast.unparse(func.value), paths)
    if func.attr in HANDLE_METHODS:
        return make(func.attr, ast.unparse(func.value), [])
    return None


def _gate_name(call: ast.Call) -> Optional[str]:
    f = call.func
    name = f.attr if isinstance(f, ast.Attribute) else (
        f.id if isinstance(f, ast.Name) else None)
    return name if name in GATES else None


def _gate_paths(call: ast.Call, name: str) -> Optional[tuple[str, ...]]:
    slot = GATES[name]
    if slot is None:
        return None
    if slot < len(call.args):
        paths = [_expr_key(call.args[slot])]
    else:
        kw = next((k for k in call.keywords
                   if k.arg in ("path", "rel", "full", "dst")), None)
        if kw is None:
            return None
        paths = [_expr_key(kw.value)]
    parent = getattr(call, "delfin_parent", None)
    if (name in _GATES_RETURNING_PATH and isinstance(parent, ast.Assign)
            and len(parent.targets) == 1):
        t = parent.targets[0]
        first = t.elts[0] if isinstance(t, ast.Tuple) and t.elts else t
        if isinstance(first, ast.Name):
            paths.append(first.id)
    return tuple(paths)


# ---------------------------------------------------------------------------
# Control flow
# ---------------------------------------------------------------------------


def _is_perms_guard(test: ast.expr) -> bool:
    """``perms is not None`` / ``permissions is not None``."""
    return (isinstance(test, ast.Compare) and len(test.ops) == 1
            and isinstance(test.ops[0], ast.IsNot)
            and isinstance(test.left, ast.Name)
            and test.left.id in _PERMS_NAMES
            and isinstance(test.comparators[0], ast.Constant)
            and test.comparators[0].value is None)


def _index_statements(tree: ast.AST) -> None:
    """Attach ``delfin_parent`` to every node, and to every statement its
    ``delfin_preorder`` (depth-first order) and ``delfin_forkpath``: one frame
    ``(id(fork), id(branch), is_perms_guard_body)`` per conditional fork that
    contains it, outer first.  ``with`` and ``def`` bodies add no frame."""
    for parent in ast.walk(tree):
        for child in ast.iter_child_nodes(parent):
            child.delfin_parent = parent

    counter = 0
    forks: list[tuple] = []

    def walk(body) -> None:
        nonlocal counter
        for stmt in body or ():
            if not isinstance(stmt, ast.stmt):
                continue
            counter += 1
            stmt.delfin_preorder = counter
            stmt.delfin_forkpath = tuple(forks)
            descend(stmt)

    def branch(fork, body, guard=False) -> None:
        forks.append((id(fork), id(body), guard))
        walk(body)
        forks.pop()

    def descend(stmt: ast.stmt) -> None:
        if isinstance(stmt, ast.If):
            branch(stmt, stmt.body, _is_perms_guard(stmt.test))
            branch(stmt, stmt.orelse)
        elif isinstance(stmt, (ast.For, ast.While)):
            branch(stmt, stmt.body)
            branch(stmt, stmt.orelse)
        elif isinstance(stmt, ast.Try):
            for bl in (stmt.body, stmt.orelse, stmt.finalbody):
                branch(stmt, bl)
            for handler in stmt.handlers:
                branch(stmt, handler.body)
        elif isinstance(stmt, ast.Match):
            for case in stmt.cases:
                branch(stmt, case.body)
        elif hasattr(stmt, "body"):
            walk(stmt.body)  # with / def / class: always-run bodies

    walk(tree.body)


def _stmt_ctx(node: ast.AST) -> tuple[Optional[int], Optional[tuple]]:
    cur = node
    while cur is not None and not isinstance(cur, ast.stmt):
        cur = getattr(cur, "delfin_parent", None)
    if cur is None:
        return (None, None)
    return (getattr(cur, "delfin_preorder", None),
            getattr(cur, "delfin_forkpath", None))


def _on_every_path(gate: GateCall, sink: Sink) -> bool:
    """Does the gate run on every path that reaches the sink, before it?

    The common prefix of the two fork paths is shared.  Two different
    branches of one fork at the first difference reject.  Otherwise the
    gate's remaining frames are a conditional that closed before the sink,
    and each must be a ``perms is not None`` body; the sink's remaining
    frames (a later ``if`` around the sink) do not matter."""
    if gate.preorder is None or sink.preorder is None:
        return False
    if gate.preorder >= sink.preorder:
        return False
    fg, fs = gate.forkpath or (), sink.forkpath or ()
    i = 0
    while i < len(fg) and i < len(fs) and fg[i] == fs[i]:
        i += 1
    if i < len(fg) and i < len(fs) and fg[i][0] == fs[i][0]:
        return False
    return all(frame[2] for frame in fg[i:])


def _sink_paths(sink: Sink, bindings: dict[str, str]) -> Optional[tuple]:
    """The paths a sink touches, or None when they cannot be named (a
    handle method on a handle not bound to an ``open``)."""
    if sink.path_exprs:
        return sink.path_exprs
    if sink.receiver in bindings:
        return (bindings[sink.receiver],)
    return None


def _covering_gate(
    gates: Iterable[GateCall], sink: Sink, bindings: dict[str, str],
) -> Optional[GateCall]:
    """The first gate on every path to ``sink`` such that, together with the
    other such gates, every path the sink touches is checked; else None."""
    ok = sorted((g for g in gates if _on_every_path(g, sink)),
                key=lambda g: (g.preorder, g.line))
    if not ok:
        return None
    if any(g.path_args is None for g in ok):
        return ok[0]
    paths = _sink_paths(sink, bindings)
    if paths is None:
        return ok[0]  # an unnamed handle: no mismatch can be shown
    checked = {p for g in ok for p in g.path_args}
    if all(p in checked for p in paths):
        return next(g for g in ok if set(g.path_args) & set(paths))
    return None


# ---------------------------------------------------------------------------
# Per-function scan
# ---------------------------------------------------------------------------


@dataclass
class _MethodScan:
    name: str
    line: int
    end_line: int = 0
    sinks: list[Sink] = field(default_factory=list)
    gates: list[GateCall] = field(default_factory=list)
    #: (callee, call preorder, call forkpath, positional argument texts)
    helper_edges: list[tuple] = field(default_factory=list)
    #: handle name -> path it was opened from
    bindings: dict[str, str] = field(default_factory=dict)
    temp_sources: set[str] = field(default_factory=set)
    param_names: list[str] = field(default_factory=list)


class _BodyVisitor(ast.NodeVisitor):
    """Sinks, gates and same-module call edges of one function.  Sinks of a
    nested function are kept with that function as owner; its gates and
    calls are not the outer function's (it may never run)."""

    def __init__(self, scan: _MethodScan, local_defs: set[str]) -> None:
        self.scan = scan
        self._local_defs = local_defs
        self._owners: list[str] = []

    def visit_FunctionDef(self, node):
        self._owners.append(node.name)
        self.generic_visit(node)
        self._owners.pop()

    visit_AsyncFunctionDef = visit_FunctionDef

    def visit_Call(self, node):
        owner = self._owners[-1] if self._owners else self.scan.name
        own = owner == self.scan.name
        pre, fpath = _stmt_ctx(node)
        sink = _detect_sink(node)
        temp = self.scan.temp_sources
        if sink is not None and not (
                sink.receiver in temp or any(e in temp for e in sink.path_exprs)):
            self.scan.sinks.append(
                replace(sink, owner=owner, preorder=pre, forkpath=fpath))
        name = _gate_name(node) if own else None
        if name is not None:
            self.scan.gates.append(GateCall(
                line=node.lineno, name=name, preorder=pre, forkpath=fpath,
                path_args=_gate_paths(node, name)))
        func = node.func
        target = func.id if isinstance(func, ast.Name) else (
            func.attr if isinstance(func, ast.Attribute) else None)
        if own and target in self._local_defs:
            self.scan.helper_edges.append(
                (target, pre, fpath,
                 tuple(_expr_key(a) for a in node.args[:3])))
        self.generic_visit(node)

    def visit_With(self, node):
        # ``with open(p) as f`` / ``with p.open() as f`` bind f to p; a
        # tempfile context marks f as staging.
        for item in node.items:
            var, ctx = item.optional_vars, item.context_expr
            if not isinstance(var, ast.Name):
                continue
            if _is_temp_source_call(ctx, self.scan.temp_sources):
                self.scan.temp_sources.add(var.id)
            elif (isinstance(ctx, ast.Call) and isinstance(ctx.func, ast.Name)
                    and ctx.func.id == "open" and ctx.args):
                self.scan.bindings[var.id] = _expr_key(ctx.args[0])
            elif (isinstance(ctx, ast.Call)
                    and isinstance(ctx.func, ast.Attribute)
                    and ctx.func.attr == "open"):
                self.scan.bindings[var.id] = _expr_key(ctx.func.value)
        self.generic_visit(node)

    def visit_Assign(self, node):
        # ``fd, name = tempfile.mkstemp()`` marks both names as staging, and
        # so does ``p = Path(name)``; ``fh = open(p)`` binds fh to p.
        val = node.value
        if len(node.targets) == 1 and isinstance(val, ast.Call):
            t = node.targets[0]
            temp = self.scan.temp_sources
            if _is_temp_source_call(val, temp) or (
                    isinstance(t, ast.Name) and _expr_key(val) in temp):
                elts = t.elts if isinstance(t, ast.Tuple) else [t]
                temp.update(e.id for e in elts if isinstance(e, ast.Name))
            elif (isinstance(t, ast.Name) and isinstance(val.func, ast.Name)
                    and val.func.id == "open" and val.args):
                self.scan.bindings[t.id] = _expr_key(val.args[0])
        self.generic_visit(node)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def analyze(source: str, module_path: str = "") -> CoverageReport:
    """Parse ``source`` and scan every ``_execute_*`` function.

    Raises ``SyntaxError`` on unparsable input.  A call from one executor to
    another (the dispatcher) is not followed: each executor is reported on
    its own, and following it charged every tool's sinks to the dispatcher.
    """
    tree = ast.parse(source)
    _index_statements(tree)
    defs = [d for d in ast.walk(tree) if isinstance(d, ast.FunctionDef)]
    all_defs = {d.name for d in defs}

    methods: dict[str, _MethodScan] = {}
    for node in defs:
        scan = _MethodScan(
            node.name, node.lineno, node.end_lineno, param_names=[
                a.arg for a in node.args.posonlyargs + node.args.args])
        _BodyVisitor(scan, all_defs).generic_visit(node)
        methods[node.name] = scan

    def callees(name: str) -> list[str]:
        return [c for c, *_ in methods[name].helper_edges if _followed(c)]

    def _followed(callee: str) -> bool:
        # Not into another executor (scanned on its own) and not into a
        # gate (its body is the check; its own reads, such as the change
        # preview of the path it has just resolved, are not the tool's).
        return (callee in methods and not callee.startswith("_execute_")
                and callee not in GATES)

    reaches: dict[str, set[str]] = {}
    for name in methods:
        seen: set[str] = set()
        stack = callees(name)
        while stack:
            cur = stack.pop()
            if cur not in seen:
                seen.add(cur)
                stack.extend(callees(cur))
        reaches[name] = seen

    report = CoverageReport(module_path=module_path)
    for name, scan in methods.items():
        if not name.startswith("_execute_"):
            continue
        finding = ExecutorFinding(name=name, line=scan.line,
                                  end_line=scan.end_line)
        for m in [scan] + [methods[c] for c in sorted(reaches[name])]:
            finding.sinks.extend(m.sinks)
            finding.gates.extend(m.gates)
        finding.sinks.sort(key=lambda s: (s.line, s.col))
        finding.results = [
            _result(s, _gate_for(s, name, methods, reaches))
            for s in finding.sinks]
        report.executors.append(finding)
    report.executors.sort(key=lambda e: e.line)
    return report


def _result(sink: Sink, gate: Optional[GateCall]) -> SinkResult:
    if gate is None:
        return SinkResult(sink, False, None, None)
    return SinkResult(sink, True, gate.line, gate.name)


def _gate_for(sink, executor, methods, reaches) -> Optional[GateCall]:
    owner = sink.owner or executor
    if owner == executor or owner not in reaches[executor]:
        # The executor itself, or a nested function it passes on (a
        # callback): only that function's own gates can cover the sink.
        own = methods.get(owner) or methods[executor]
        return _covering_gate(own.gates, sink, own.bindings)
    return _gate_on_call_path(executor, owner, sink, methods, reaches,
                              frozenset())


def _gate_on_call_path(here, target, sink, methods, reaches, visited):
    """Walk the calls from ``here`` down to ``target`` (the sink's owner).

    At each call edge on the way, a gate in ``here`` that covers the call
    statement, with the callee's parameters bound back to this call's
    arguments (``def h(p): open(p)`` called as ``h(gated)``), protects the
    sink; otherwise descend.  In ``target`` its own gates decide."""
    if here in visited:
        return None
    visited = visited | {here}
    scan = methods[here]
    if here == target:
        return _covering_gate(scan.gates, sink, scan.bindings)
    for child, pre, fpath, arg_texts in scan.helper_edges:
        if child not in reaches[here]:
            continue
        if child != target and target not in reaches[child]:
            continue
        params = methods[child].param_names
        if params[:1] in (["self"], ["cls"]):
            params = params[1:]
        sub = dict(zip(params, arg_texts))
        paths = _sink_paths(sink, methods[target].bindings) or ()
        edge_sink = replace(sink, preorder=pre, forkpath=fpath,
                            path_exprs=tuple(sub.get(p, p) for p in paths))
        gate = _covering_gate(scan.gates, edge_sink, {})
        if gate is not None:
            return gate
        return _gate_on_call_path(child, target, sink, methods, reaches,
                                  visited)
    return None


def format_report(report: CoverageReport) -> str:
    """One line per ungated sink."""
    head = report.module_path or "<source>"
    if not report.ungated():
        return f"{head}: no ungated file sinks."
    lines = [f"{head}: {len(report.ungated())} executor(s) with ungated "
             f"file sink(s):"]
    for ex in report.ungated():
        for r in ex.ungated_sinks:
            s = r.sink
            lines.append(
                f"  {ex.name}@{s.line}:{s.col}  {s.kind}"
                f"{f'({s.receiver})' if s.receiver else ''}"
                f"{f'  arg={s.arg!r}' if s.arg else ''}  UNGATED")
    return "\n".join(lines)
