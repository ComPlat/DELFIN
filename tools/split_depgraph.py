"""Map the internal structure of a monolith module.

Usage: python tools/split_depgraph.py <module.py> <depgraph.json>

Top-level statements, the names each top-level definition references, the
function-level dependency graph, its strongly connected components and the
section headers the author left in the file.  Read-only; prints a report.
"""
import ast
import sys
import json
from collections import defaultdict

path = sys.argv[1]
src = open(path).read()
tree = ast.parse(src)
lines = src.splitlines()

# ---- top-level statements ----------------------------------------------
tops = []  # (kind, names, lineno, end_lineno, node)
defined = {}  # name -> index into tops


def _names_of_target(t):
    if isinstance(t, ast.Name):
        return [t.id]
    if isinstance(t, (ast.Tuple, ast.List)):
        out = []
        for e in t.elts:
            out += _names_of_target(e)
        return out
    return []


def _stmt_defs(node):
    """Names a top-level statement binds (recursing into try/if bodies)."""
    names = []
    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
        names.append(node.name)
    elif isinstance(node, ast.Assign):
        for t in node.targets:
            names += _names_of_target(t)
    elif isinstance(node, (ast.AnnAssign, ast.AugAssign)):
        names += _names_of_target(node.target)
    elif isinstance(node, (ast.Import, ast.ImportFrom)):
        for a in node.names:
            names.append((a.asname or a.name).split(".")[0])
    elif isinstance(node, (ast.Try, ast.If, ast.With, ast.For, ast.While)):
        for field in ("body", "orelse", "finalbody", "handlers"):
            for sub in getattr(node, field, []) or []:
                if isinstance(sub, ast.ExceptHandler):
                    for s2 in sub.body:
                        names += _stmt_defs(s2)
                else:
                    names += _stmt_defs(sub)
    return names


for i, node in enumerate(tree.body):
    kind = type(node).__name__
    names = _stmt_defs(node)
    tops.append((kind, names, node.lineno, node.end_lineno, node))
    for n in names:
        defined.setdefault(n, []).append(i)

module_names = set(defined)


class FreeNames(ast.NodeVisitor):
    """Names loaded in a definition that are not bound locally (approx)."""

    def __init__(self):
        self.loads = set()
        self.stores = set()
        self.globals_decl = set()

    def visit_Name(self, node):
        if isinstance(node.ctx, ast.Load):
            self.loads.add(node.id)
        else:
            self.stores.add(node.id)

    def visit_Global(self, node):
        self.globals_decl.update(node.names)

    def visit_Attribute(self, node):
        self.generic_visit(node)


refs = {}
glob_writers = {}
for i, (kind, names, l0, l1, node) in enumerate(tops):
    v = FreeNames()
    v.visit(node)
    # any loaded name that is a module-level name counts as a dependency
    # (locals shadowing a module name are rare and only add a false edge)
    deps = (v.loads & module_names) - set(names)
    refs[i] = deps
    if v.globals_decl:
        glob_writers[i] = v.globals_decl

# ---- graph on definitions -----------------------------------------------
idx_of = {}
for n, ids in defined.items():
    idx_of[n] = ids[-1]
edges = defaultdict(set)
for i, deps in refs.items():
    for d in deps:
        j = idx_of[d]
        if j != i:
            edges[i].add(j)

# Tarjan SCC
sys.setrecursionlimit(100000)
index = {}
low = {}
onstack = set()
stack = []
sccs = []
counter = [0]


def strongconnect(v):
    index[v] = low[v] = counter[0]
    counter[0] += 1
    stack.append(v)
    onstack.add(v)
    for w in edges[v]:
        if w not in index:
            strongconnect(w)
            low[v] = min(low[v], low[w])
        elif w in onstack:
            low[v] = min(low[v], index[w])
    if low[v] == index[v]:
        comp = []
        while True:
            w = stack.pop()
            onstack.discard(w)
            comp.append(w)
            if w == v:
                break
        sccs.append(comp)


for v in range(len(tops)):
    if v not in index:
        strongconnect(v)

big = [c for c in sccs if len(c) > 1]
print("top-level statements:", len(tops))
print("definitions (names):", len(module_names))
print("edges:", sum(len(e) for e in edges.values()))
print("SCCs with >1 member:", len(big))
for c in sorted(big, key=len, reverse=True):
    print("  SCC size", len(c), ":", [tops[i][1][0] if tops[i][1] else tops[i][0] for i in sorted(c)][:40])
print("global writers:", {tops[i][1][0]: sorted(g) for i, g in glob_writers.items()})

# ---- section headers ----------------------------------------------------
print("\n== section headers (comment banners at column 0) ==")
for ln, line in enumerate(lines, 1):
    if line.startswith("# ===") or line.startswith("# ---") or line.startswith("# ###") or line.startswith("#####"):
        nxt = lines[ln] if ln < len(lines) else ""
        print(f"{ln:6d}: {line[:100]}")
        if nxt.startswith("# ") and not nxt.startswith("# ==="):
            print(f"        {nxt[:100]}")

# ---- non-def top-level statements ---------------------------------------
print("\n== top-level statements that are not def/class/import ==")
for i, (kind, names, l0, l1, node) in enumerate(tops):
    if kind in ("FunctionDef", "ClassDef", "Import", "ImportFrom"):
        continue
    print(f"{l0:6d}-{l1:6d} {kind:10s} {names[:6]}")

# dump machine-readable
out = {
    "tops": [
        {"i": i, "kind": k, "names": n, "l0": l0, "l1": l1,
         "deps": sorted(tops[j][1][0] if tops[j][1] else "?" for j in edges[i])}
        for i, (k, n, l0, l1, _) in enumerate(tops)
    ]
}
json.dump(out, open(sys.argv[2], "w"), indent=1)
