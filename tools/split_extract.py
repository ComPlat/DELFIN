"""Mechanical extractor for the MANTA monolith split.

Moves whole top-level units (a def/class/assignment plus the unnamed statements
glued to it: attribute docstrings, module-level calls, loops that fill a table)
verbatim from the shim module into a new module, writes the import block the
new module needs, and replaces the moved units in the shim by a re-export
import of every name the new module defines.  No statement text is edited.

    python tools/split_extract.py --shim delfin/smiles_converter.py \
        --module delfin/manta/NAME.py --names a,b,c --doc "..." [--dry-run]

    python tools/split_extract.py --shim ... --plan plan.json   # closure check only

    python tools/split_extract.py --shim ... --plan-step plan.json NAME \
        --registry split_registry.json [--logger-name delfin.smiles_converter]

A registry (name -> module) persists which names already left the shim.
"""
from __future__ import annotations

import argparse
import ast
import json
import symtable
import sys
from pathlib import Path

HEADER_IMPORTS = {}   # name -> import line, filled from the shim's own header
LOGGER_LINE = None    # set from --logger-name


def header_imports(tree):
    """Every name the shim binds through a top-level import, with its line."""
    table = {}
    for node in tree.body:
        if isinstance(node, ast.ImportFrom):
            if node.module == "__future__":
                continue
            for a in node.names:
                bound = a.asname or a.name
                table[bound] = "from %s%s import %s" % (
                    "." * node.level, node.module or "",
                    a.name + (" as " + a.asname if a.asname else ""))
        elif isinstance(node, ast.Import):
            for a in node.names:
                bound = (a.asname or a.name).split(".")[0]
                table[bound] = "import " + a.name + (" as " + a.asname if a.asname else "")
    return table


BUILTINS = set(dir(__builtins__)) | {"__name__", "__file__", "__doc__", "logger"}


# ----------------------------------------------------------------------------
# parsing the shim into units
# ----------------------------------------------------------------------------
def _names_of_target(t):
    if isinstance(t, ast.Name):
        return [t.id]
    if isinstance(t, (ast.Tuple, ast.List)):
        out = []
        for e in t.elts:
            out += _names_of_target(e)
        return out
    return []


def stmt_defs(node):
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
                        names += stmt_defs(s2)
                else:
                    names += stmt_defs(sub)
    return names


class Unit:
    def __init__(self, idx, names, stmts):
        self.idx = idx
        self.names = names          # names defined by the leading statement
        self.stmts = stmts          # ast nodes (leading + glued)
        self.l0 = stmts[0].lineno   # first line of the leading statement
        self.l1 = stmts[-1].end_lineno
        self.start = None           # text start line (after the previous unit)
        self.text = None
        self.refs = set()


def parse_units(src: str):
    tree = ast.parse(src)
    lines = src.splitlines(keepends=True)
    units = []
    pending_unnamed = []
    for node in tree.body:
        names = stmt_defs(node)
        is_named = bool(names) or isinstance(node, (ast.FunctionDef, ast.ClassDef, ast.If))
        if isinstance(node, ast.If):
            # an If with a test on __name__ is its own (unnamed) unit kept in the shim
            is_named = True
        if is_named:
            units.append(Unit(len(units), names, [node]))
        elif units:
            units[-1].stmts.append(node)
            units[-1].l1 = node.end_lineno
        else:
            pending_unnamed.append(node)
    # text slices: from the end of the previous unit to the end of this unit
    prev_end = 0
    if pending_unnamed:
        prev_end = pending_unnamed[-1].end_lineno
    for u in units:
        u.start = prev_end + 1
        u.text = "".join(lines[prev_end:u.l1])
        prev_end = u.l1
    tail = "".join(lines[prev_end:])
    head = "".join(lines[: (pending_unnamed[-1].end_lineno if pending_unnamed else 0)])
    return tree, units, head, tail


def _annotation_refs(node):
    """Names used in annotations anywhere inside the definition.

    With ``from __future__ import annotations`` they are strings at run time and
    symtable does not report them, but pyflakes/ruff resolve them, so the
    module needs the typing imports.
    """
    out = set()
    parts = []
    for n in ast.walk(node):
        if isinstance(n, ast.arg) and n.annotation is not None:
            parts.append(n.annotation)
        elif isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.returns is not None:
            parts.append(n.returns)
        elif isinstance(n, ast.AnnAssign):
            parts.append(n.annotation)
    for p in parts:
        for n in ast.walk(p):
            if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load):
                out.add(n.id)
    return out


def _decl_refs(node):
    """Names evaluated in the enclosing scope: decorators, defaults, bases."""
    out = set()
    parts = []
    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
        parts += node.decorator_list
        a = node.args
        parts += [d for d in a.defaults if d is not None]
        parts += [d for d in a.kw_defaults if d is not None]
    elif isinstance(node, ast.ClassDef):
        parts += node.decorator_list + node.bases + [k.value for k in node.keywords]
    for p in parts:
        for n in ast.walk(p):
            if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load):
                out.add(n.id)
    return out


def _table_globals(tab, acc):
    for sym in tab.get_symbols():
        if sym.is_global() and (sym.is_referenced() or sym.is_assigned()):
            acc.add(sym.get_name())
        if sym.is_declared_global():
            acc.add(sym.get_name())
    for child in tab.get_children():
        _table_globals(child, acc)


def compute_refs(src, units):
    """Free module-level names each unit references."""
    top = symtable.symtable(src, "<shim>", "exec")
    child_by_line = {}
    for ch in top.get_children():
        child_by_line.setdefault(ch.get_lineno(), []).append(ch)
    for u in units:
        refs = set()
        for st in u.stmts:
            if isinstance(st, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                for ch in child_by_line.get(st.lineno, []):
                    if ch.get_name() == st.name:
                        _table_globals(ch, refs)
                refs |= _decl_refs(st)
                refs |= _annotation_refs(st)
            else:
                # module-level statement: every loaded name (comprehension
                # variables fall out later: they are not module-level names)
                for n in ast.walk(st):
                    if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load):
                        refs.add(n.id)
                    elif isinstance(n, (ast.FunctionDef, ast.Lambda)):
                        pass
                # lambdas / nested defs inside module-level statements
                for n in ast.walk(st):
                    if isinstance(n, (ast.Lambda, ast.FunctionDef)):
                        for ch in child_by_line.get(n.lineno, []):
                            _table_globals(ch, refs)
        u.refs = refs - set(u.names)
    return units


# ----------------------------------------------------------------------------
# planning / closure check
# ----------------------------------------------------------------------------
def closure(units, group_names, registry, defined_in_shim):
    """Which referenced names the group needs from where."""
    group = [u for u in units if u.names and u.names[0] in group_names]
    group_defined = set()
    for u in group:
        group_defined |= set(u.names)
    missing_from_shim = []
    need_registry = {}
    need_header = set()
    for u in group:
        for r in sorted(u.refs):
            if r in group_defined or r in BUILTINS:
                continue
            if r in registry:
                need_registry.setdefault(registry[r], set()).add(r)
            elif r in HEADER_IMPORTS:
                need_header.add(r)
            elif r in defined_in_shim:
                missing_from_shim.append((u.names[0], r))
            else:
                # not a module-level name anywhere (comprehension var, nested
                # local the symtable walk reported as global, ...): ignore
                pass
    return group, group_defined, need_registry, need_header, missing_from_shim


def render_module(doc, group, need_registry, need_header, origin):
    out = []
    out.append(f'"""{doc}\n\nMoved verbatim from {origin} (MANTA split, 2026-10); every\nstatement is the original text, only the imports below are new.\n"""\n')
    out.append("from __future__ import annotations\n\n")
    lines = {HEADER_IMPORTS[n] for n in need_header if n != "get_logger"}
    plain = sorted(s for s in lines if s.startswith("import "))
    froms = {}
    for s_ in lines:
        if s_.startswith("from "):
            mod, name = s_[5:].split(" import ", 1)
            froms.setdefault(mod, []).append(name)
    std_plain = [s for s in plain if not s.startswith("import delfin")]
    loc_plain = [s for s in plain if s.startswith("import delfin")]
    std_from = {m: v for m, v in froms.items() if not m.startswith("delfin")}
    loc_from = {m: v for m, v in froms.items() if m.startswith("delfin")}
    for s_ in std_plain:
        out.append(s_ + "\n")
    for m in sorted(std_from):
        out.append("from %s import %s\n" % (m, ", ".join(sorted(std_from[m]))))
    out.append("\n")
    if LOGGER_LINE:
        out.append("from delfin.common.logging import get_logger\n")
    for s_ in loc_plain:
        out.append(s_ + "\n")
    for m in sorted(loc_from):
        out.append("from %s import %s\n" % (m, ", ".join(sorted(loc_from[m]))))
    for mod in sorted(need_registry):
        names = sorted(need_registry[mod])
        out.append(_import_block(mod, names))
    if LOGGER_LINE:
        out.append("\n" + LOGGER_LINE + "\n")
    for u in group:
        body = u.text.lstrip("\n")
        out.append("\n\n" + body.rstrip("\n") + "\n")
    return "".join(out)


def _import_block(mod, names, noqa=False):
    line = f"from {mod} import ("
    body = ",\n".join("    " + n for n in names)
    tail = ")" + ("  # noqa: F401" if noqa else "")
    return f"{line}\n{body},\n{tail}\n"


def rewrite_shim(src, units, head, tail, group, modname, group_defined, insert_after_unit):
    """Shim without the group's units, with a re-export import inserted."""
    keep = [u for u in units if u not in group]
    out = [head]
    inserted = False
    for u in keep:
        out.append(u.text)
        if not inserted and u.idx >= insert_after_unit:
            out.append("\n" + _import_block(modname, sorted(group_defined), noqa=True))
            inserted = True
    out.append(tail)
    return "".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--shim", required=True)
    ap.add_argument("--module", help="new module file path (delfin/manta/x.py)")
    ap.add_argument("--names", help="comma-separated names of the leading statements")
    ap.add_argument("--names-file", help="file with one name per line")
    ap.add_argument("--doc", default="MANTA constructor module.")
    ap.add_argument("--registry", default="split_registry.json")
    ap.add_argument("--plan", help="JSON [{module, names|names_file, doc}] closure check")
    ap.add_argument("--plan-step", nargs=2, metavar=("PLAN", "MODULE"),
                    help="take module, names and doc of one step from the plan")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--logger-name", default=None,
                    help="emit logger = get_logger(NAME) in every module")
    ap.add_argument("--origin", default="delfin/smiles_converter.py")
    ap.add_argument("--insert-after", default="logger",
                    help="shim unit name after which the re-export import goes")
    args = ap.parse_args(argv)

    shim_path = Path(args.shim)
    src = shim_path.read_text()
    tree, units, head, tail = parse_units(src)
    compute_refs(src, units)
    global LOGGER_LINE
    HEADER_IMPORTS.update(header_imports(tree))
    if args.logger_name:
        HEADER_IMPORTS["get_logger"] = "from delfin.common.logging import get_logger"
        LOGGER_LINE = 'logger = get_logger("%s")' % args.logger_name
    reg_path = Path(args.registry)
    registry = json.loads(reg_path.read_text()) if reg_path.exists() else {}
    defined_in_shim = set()
    for u in units:
        defined_in_shim |= set(u.names)

    def _load_names(spec_names, spec_file):
        names = []
        if spec_names:
            names += [n.strip() for n in spec_names.split(",") if n.strip()]
        if spec_file:
            names += [ln.strip() for ln in Path(spec_file).read_text().splitlines()
                      if ln.strip() and not ln.startswith("#")]
        return names

    if args.plan:
        plan = json.loads(Path(args.plan).read_text())
        sim_registry = dict(registry)
        sim_shim = set(defined_in_shim)
        ok = True
        for step in plan:
            names = _load_names(step.get("names"), step.get("names_file"))
            modname = Path(step["module"]).with_suffix("").as_posix().replace("/", ".")
            unknown = [n for n in names if n not in sim_shim]
            group, gdef, need_reg, need_hdr, missing = closure(
                units, set(names), sim_registry, sim_shim)
            nlines = sum(u.l1 - u.start + 1 for u in group)
            print(f"== {modname}: {len(group)} units, {nlines} lines, needs "
                  f"{ {m.rsplit('.',1)[-1]: len(v) for m, v in need_reg.items()} }")
            if unknown:
                print("   UNKNOWN names:", unknown)
                ok = False
            if missing:
                ok = False
                print("   NOT YET EXTRACTED (still in shim):")
                by = {}
                for u, r in missing:
                    by.setdefault(r, []).append(u)
                for r, us in sorted(by.items()):
                    print(f"      {r}  <- {', '.join(sorted(set(us)))}")
            for n in gdef:
                sim_registry[n] = modname
                sim_shim.discard(n)
        rest = sum(u.l1 - u.start + 1 for u in units if u.names and u.names[0] in sim_shim)
        print(f"== shim keeps {rest} lines in named units; plan {'OK' if ok else 'HAS GAPS'}")
        return 0 if ok else 1

    if args.plan_step:
        plan = json.loads(Path(args.plan_step[0]).read_text())
        step = [s for s in plan if Path(s["module"]).stem == args.plan_step[1]]
        if len(step) != 1:
            print("plan step not found:", args.plan_step[1])
            return 2
        args.module = step[0]["module"]
        args.names_file = step[0].get("names_file")
        args.names = step[0].get("names")
        args.doc = step[0]["doc"]
    names = _load_names(args.names, args.names_file)
    modname = Path(args.module).with_suffix("").as_posix().replace("/", ".")
    unknown = [n for n in names if n not in defined_in_shim]
    if unknown:
        print("unknown names:", unknown)
        return 2
    group, gdef, need_reg, need_hdr, missing = closure(units, set(names), registry,
                                                       defined_in_shim)
    if missing:
        print("group references names still in the shim:")
        for u, r in missing:
            print(f"   {u} -> {r}")
        return 3
    text = render_module(args.doc, group, need_reg, need_hdr, args.origin)
    # the shim: where to insert the import
    after = [u.idx for u in units if u.names and u.names[0] == args.insert_after]
    # after the last re-export block already in the shim, so the blocks read in
    # extraction order
    prior = [u.idx for u in units
             if isinstance(u.stmts[0], ast.ImportFrom)
             and (u.stmts[0].module or "").startswith("delfin.manta")]
    insert_after_unit = max(after + prior) if (after or prior) else 0
    new_shim = rewrite_shim(src, units, head, tail, group, modname, gdef, insert_after_unit)
    nlines = sum(u.l1 - u.start + 1 for u in group)
    print(f"{modname}: {len(group)} units, {nlines} lines; imports from "
          f"{sorted(m.rsplit('.', 1)[-1] for m in need_reg)}; header {sorted(need_hdr)}")
    if args.dry_run:
        return 0
    Path(args.module).write_text(text)
    shim_path.write_text(new_shim)
    for n in gdef:
        registry[n] = modname
    reg_path.write_text(json.dumps(registry, indent=0, sort_keys=True))
    ast.parse(text)
    ast.parse(new_shim)
    print(f"wrote {args.module} and rewrote {shim_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
