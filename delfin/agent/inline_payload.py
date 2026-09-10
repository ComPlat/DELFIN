"""Reading a ``python -c`` payload instead of giving up on it.

An interpreter is kept off the bash auto-allow list because a scanner that
decides on the command string cannot see what the command will do:
``xargs`` reads its targets from a file, ``make`` runs whatever the
Makefile says, ``base64 -d | bash`` carries its program encoded. For all
of those the objection is exact -- the text is genuinely not there.

``python -c`` is the one member of that family where it IS there. The
program is a literal in the command line, and the two checks that the
interpreter rule stands in for can both be run on it:

* the CONTENT scan, which every executed script file already gets
  (:meth:`_scan_bash_script_payloads`) and which an inline payload
  escaped only because :meth:`_referenced_script_paths` has no file to
  open;
* the WRITE gate, which reads paths out of a command
  (:func:`_bash_write_targets` knows cp, mv, tee, sed -i, dd) and knew
  nothing about ``open(path, 'w')``.

So this module answers one question about a payload: **what would it
write, and is there anything in it I cannot account for?** The two
outcomes are deliberately asymmetric.

``opaque`` is the default. Every construct that is not on the small list
below -- an unparsed payload, a computed path, ``exec``, ``subprocess``,
``os.system``, an attribute this module does not model -- yields a reason
string, and the caller keeps doing exactly what it does today: refuse,
and name the allowed way. Nothing that is refused now starts running
because a construct was overlooked; a construct that is overlooked stays
refused.

``writes`` is the earned outcome. A payload that only reads, prints and
computes has no write targets and nothing opaque, so the write gate has
nothing to gate and the interpreter objection has nothing left to stand
for. A payload that writes a LITERAL path yields that path, and it goes
through the same gate as ``tee out.txt`` -- which is a containment gain,
because today a write through python is invisible to that gate wherever
the command runs at all.

The local-import rule is the third half of it. ``import mymod`` executes
``mymod.py`` from the working directory, and that file may be one this
session wrote a minute ago. That is not an argument for refusing the
payload -- ``python3 mymod.py`` runs the same file and is auto-allowed --
it is an argument for scanning the same file. So imports that resolve
inside the working directory come back as :attr:`Effects.local_modules`
and the caller feeds them to the content scan it already has.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass, field
from pathlib import Path

__all__ = ["Effects", "analyze_payload", "extract_c_payloads"]


# Modules whose surface is read-and-compute. Importing one of these
# executes stdlib code only, so nothing here can be a file this session
# wrote. The list is short on purpose: it covers what the recorded
# refusals actually asked for (json, csv, a workbook reader, arithmetic)
# and grows only when a real refusal shows a real gap.
_SAFE_MODULES = frozenset({
    "ast", "base64", "binascii", "bisect", "calendar", "cmath",
    "collections", "colorsys", "copy", "csv", "datetime", "decimal",
    "difflib", "enum", "fnmatch", "fractions", "functools", "gzip",
    "hashlib", "heapq", "html", "itertools", "json", "keyword", "locale",
    "math", "numbers", "operator", "pprint", "random", "re", "statistics",
    "string", "textwrap", "time", "tomllib", "types", "typing",
    "unicodedata", "uuid", "warnings", "zoneinfo",
    # `sys` is here for one idiom that appears in the refusal log five
    # times: `sys.path.insert(0, '.')` in front of a local import. The
    # module itself reads and computes; what it can reach that matters --
    # `sys.stdout.write` -- is judged by the attribute rules below like
    # any other write, and the local import is judged as a local import.
    "sys",
    # Scientific and document readers. Each is a third-party package the
    # scientist's shell one-liner is actually made of; each can also
    # write, which is why the WRITE side below models their writers by
    # name rather than trusting the import.
    "numpy", "pandas", "scipy", "openpyxl", "docx", "pypdf", "PyPDF2",
    "matplotlib", "yaml", "h5py", "ase", "rdkit", "cclib",
})

# Modules that reach the filesystem, the process table or the network in
# ways this module does not model. Importing one is not a crime; it makes
# the payload opaque, which means "ask", not "refuse forever".
_OPAQUE_MODULES = frozenset({
    "os", "shutil", "subprocess", "pathlib", "glob", "tempfile",
    "socket", "urllib", "requests", "http", "ftplib", "smtplib",
    "importlib", "ctypes", "pickle", "shelve", "sqlite3", "multiprocessing",
    "threading", "signal", "pty", "resource", "mmap", "webbrowser",
})

# Callables that execute text or a name resolved at run time. The payload
# stops being readable at the point one of these appears.
_OPAQUE_BUILTINS = frozenset({
    "exec", "eval", "compile", "__import__", "input", "breakpoint",
    "globals", "locals", "vars", "setattr", "delattr",
})

# Attribute calls that write. Matched by NAME, without knowing the
# receiver's type -- deliberately: `x.write_text(p)` is a write whether x
# is a Path, a stub or something this module never heard of. The cost of
# the loose match is that a same-named harmless method reads as a write
# and the payload asks; that is the safe direction.
_WRITE_ATTRS = frozenset({
    "write", "writelines", "writerow", "writerows", "write_text",
    "write_bytes", "truncate", "unlink", "rmdir", "mkdir", "makedirs",
    "removedirs", "rename", "replace", "touch", "chmod", "chown",
    "symlink_to", "hardlink_to", "rmtree", "copy", "copy2", "copyfile",
    "copytree", "move", "save", "savefig", "to_csv", "to_excel",
    "to_json", "to_parquet", "to_pickle", "to_hdf", "dump", "system",
    "popen", "spawn", "spawnl", "spawnv", "fork", "execv", "execve",
})

# Write-attribute calls whose FIRST positional argument is the path
# written. `wb.save('out.xlsx')`, `df.to_csv('out.csv')`,
# `Path.rename(dst)` -- knowing this turns an opaque call into a named
# target the write gate can actually judge.
_WRITE_ATTR_PATH_ARG0 = frozenset({
    "save", "savefig", "to_csv", "to_excel", "to_json", "to_parquet",
    "to_pickle", "to_hdf", "rename", "replace", "symlink_to",
    "hardlink_to",
})

# ... and the ones whose destination is the SECOND argument.
_WRITE_ATTR_PATH_ARG1 = frozenset({
    "copy", "copy2", "copyfile", "copytree", "move",
})

# `open(p, mode)`: the modes that create or change a file. A missing mode
# is 'r'.
_WRITE_MODE_CHARS = frozenset("wax+")

# Names that mean a write on one type and something ordinary on another,
# told apart by how many positional arguments they were given. Without
# this, `s.replace(',', '.')` -- how every German-written amount is
# parsed, and the single commonest line in the recorded refusals -- read
# as `Path.replace(target)` and reported a write to the file ".".
#
# The rule is (name -> the arities that mean a write). An arity outside
# the set is not a write and not opaque either: it is the other method.
_WRITE_ARITY: dict[str, frozenset[int]] = {
    "replace": frozenset({1}),        # Path.replace(t) vs str.replace(a, b)
    "copy": frozenset({2, 3}),        # shutil.copy(s, d) vs df.copy()
    "copy2": frozenset({2, 3}),
    "copyfile": frozenset({2}),
    "copytree": frozenset({2}),
    "move": frozenset({2}),
    "update": frozenset(),            # never a write here; dict.update
}

# Dotted calls that write to a stream the shell already owns. Printing to
# stdout is not a filesystem effect and must not read as one.
_STREAM_WRITES = frozenset({
    "sys.stdout.write", "sys.stderr.write",
    "sys.stdout.writelines", "sys.stderr.writelines",
    "sys.stdout.flush", "sys.stderr.flush",
})


@dataclass
class Effects:
    """What a payload was found to do.

    ``opaque`` non-empty means: this module could not account for the
    payload, and the caller must fall back to asking. ``writes`` and
    ``local_modules`` are only meaningful when ``opaque`` is empty.
    """

    writes: list[str] = field(default_factory=list)
    # Literal paths the payload OPENS for reading. Collected for the same
    # reason as ``writes``: the gate that stops `cat /etc/passwd` reads
    # the arguments of known content-dumping commands, and a payload's
    # `open(p).read()` is that command with the path one level in.
    reads: list[str] = field(default_factory=list)
    local_modules: list[str] = field(default_factory=list)
    opaque: list[str] = field(default_factory=list)

    @property
    def readable(self) -> bool:
        return not self.opaque


def _literal_str(node: ast.AST) -> str | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    return None


def _dotted(node: ast.AST) -> str | None:
    """`os.path.join` -> "os.path.join"; anything computed -> None."""
    parts: list[str] = []
    cur = node
    while isinstance(cur, ast.Attribute):
        parts.append(cur.attr)
        cur = cur.value
    if not isinstance(cur, ast.Name):
        return None
    parts.append(cur.id)
    return ".".join(reversed(parts))


class _Walker(ast.NodeVisitor):
    def __init__(self, cwd: Path | str | None) -> None:
        # Coerced, not merely annotated. The signature said Path and
        # every real caller passes a string -- the workspace, a base
        # directory -- so `self.cwd / f"{top}.py"` raised TypeError on
        # any payload containing an import, analyze_payload caught it,
        # and the whole payload came back opaque.
        #
        # It failed SAFE, which is why it was invisible: an opaque
        # payload is refused, exactly as before this module existed. But
        # `import math; print(math.exp(-1.5))` is the shape of nearly
        # every scientific one-liner, so the capability was dead for its
        # main case while the tests -- which pass tmp_path, a Path --
        # stayed green.
        try:
            self.cwd = Path(cwd) if cwd else None
        except TypeError:
            self.cwd = None
        self.eff = Effects()

    # -- imports ---------------------------------------------------------

    def _note_module(self, dotted: str, where: ast.AST) -> None:
        top = (dotted or "").split(".", 1)[0]
        if not top:
            self.eff.opaque.append("a relative import")
            return
        local = self._resolve_local(top)
        if local is not None:
            # A file in the working directory. Runs on import exactly as
            # `python3 <that file>` would, and that form is auto-allowed
            # WITH a content scan, so the caller gets the same file to
            # scan rather than a refusal.
            self.eff.local_modules.append(str(local))
            return
        if top in _SAFE_MODULES:
            return
        if top in _OPAQUE_MODULES:
            self.eff.opaque.append(f"imports {top}")
            return
        self.eff.opaque.append(f"imports {top}, which is not a known-read module")

    def _resolve_local(self, top: str) -> Path | None:
        if self.cwd is None:
            return None
        try:
            for cand in (self.cwd / f"{top}.py", self.cwd / top / "__init__.py"):
                if cand.is_file():
                    return cand
        except OSError:
            return None
        return None

    def visit_Import(self, node: ast.Import) -> None:
        for alias in node.names:
            self._note_module(alias.name, node)
        self.generic_visit(node)

    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        if node.level:
            self.eff.opaque.append("a relative import")
        else:
            self._note_module(node.module or "", node)
        self.generic_visit(node)

    # -- calls -----------------------------------------------------------

    def visit_Call(self, node: ast.Call) -> None:
        func = node.func
        if isinstance(func, ast.Name):
            self._plain_call(func.id, node)
        elif isinstance(func, ast.Attribute):
            self._attr_call(func, node)
        self.generic_visit(node)

    def _plain_call(self, name: str, node: ast.Call) -> None:
        if name in _OPAQUE_BUILTINS:
            self.eff.opaque.append(f"calls {name}()")
            return
        if name == "open":
            self._open_call(node)

    def _open_call(self, node: ast.Call) -> None:
        mode = "r"
        if len(node.args) >= 2:
            lit = _literal_str(node.args[1])
            if lit is None:
                self.eff.opaque.append("open() with a computed mode")
                return
            mode = lit
        for kw in node.keywords:
            if kw.arg == "mode":
                lit = _literal_str(kw.value)
                if lit is None:
                    self.eff.opaque.append("open() with a computed mode")
                    return
                mode = lit
            elif kw.arg is None:
                self.eff.opaque.append("open() with **kwargs")
                return
        path = _literal_str(node.args[0]) if node.args else None
        if not (set(mode) & _WRITE_MODE_CHARS):
            if path is not None:
                self.eff.reads.append(path)
            return
        if path is None:
            self.eff.opaque.append("open() writes a computed path")
            return
        self.eff.writes.append(path)

    def _attr_call(self, func: ast.Attribute, node: ast.Call) -> None:
        attr = func.attr
        dotted = _dotted(func)
        if dotted and dotted.split(".", 1)[0] in _OPAQUE_MODULES:
            self.eff.opaque.append(f"calls {dotted}()")
            return
        if dotted in _STREAM_WRITES:
            return
        if attr not in _WRITE_ATTRS:
            return
        if attr in _WRITE_ARITY and len(node.args) not in _WRITE_ARITY[attr]:
            return                        # the same name on another type
        if attr in _WRITE_ATTR_PATH_ARG0:
            path = _literal_str(node.args[0]) if node.args else None
            if path is not None:
                self.eff.writes.append(path)
                return
        if attr in _WRITE_ATTR_PATH_ARG1:
            path = _literal_str(node.args[1]) if len(node.args) > 1 else None
            if path is not None:
                self.eff.writes.append(path)
                return
        # A write whose destination this module cannot name. `f.write(x)`
        # on a handle from an open() already accounted for above is the
        # common case, and it is still opaque here on purpose: pairing a
        # handle to its open() is dataflow, and getting that wrong is a
        # write nobody gated.
        self.eff.opaque.append(f"calls .{attr}(), whose target is not a literal")

    # -- assignment to a filesystem-ish attribute ------------------------

    def visit_With(self, node: ast.With) -> None:
        self.generic_visit(node)


def analyze_payload(source: str, cwd: Path | None = None) -> Effects:
    """Read *source* as a Python program and report what it would do.

    *cwd* is the directory the payload would run in; it decides whether an
    import names a local file. Passing None means "cannot tell", and then
    a non-stdlib import is opaque rather than local.

    Never raises.
    """
    try:
        tree = ast.parse(source or "")
    except SyntaxError as exc:
        return Effects(opaque=[f"does not parse as Python ({exc.msg})"])
    except (ValueError, RecursionError, MemoryError) as exc:
        return Effects(opaque=[f"could not be read ({type(exc).__name__})"])
    walker = _Walker(cwd)
    try:
        walker.visit(tree)
    except Exception as exc:                                  # pragma: no cover
        return Effects(opaque=[f"could not be analysed ({type(exc).__name__})"])
    eff = walker.eff
    # Handles from a write-mode open() are the one dataflow this module
    # does do, and it does it by not doing it: the open() already recorded
    # the target, and the .write() that follows recorded an opaque entry.
    # Drop that entry only when EVERY write-mode open() in the payload
    # named a literal path, so a payload that writes one known and one
    # computed file still asks.
    if eff.writes and all(
            not o.startswith("open()") for o in eff.opaque):
        eff.opaque = [o for o in eff.opaque
                      if not o.startswith("calls .write")]
    _dedupe(eff.writes)
    _dedupe(eff.reads)
    _dedupe(eff.local_modules)
    _dedupe(eff.opaque)
    return eff


def _dedupe(items: list[str]) -> None:
    seen: set[str] = set()
    keep: list[str] = []
    for it in items:
        if it not in seen:
            seen.add(it)
            keep.append(it)
    items[:] = keep


def extract_c_payloads(cmd: str) -> list[str] | None:
    """The ``-c`` payloads in *cmd*, or None when it cannot be split.

    None and [] mean different things: None is "this command line could
    not be tokenised, so do not claim to have read it", [] is "tokenised
    fine, no inline payload in it".
    """
    import re
    import shlex

    try:
        toks = shlex.split(cmd or "", posix=True)
    except ValueError:
        return None
    out: list[str] = []
    interp = re.compile(r"(?:[\w./~+-]*/)?python[0-9.]*$")
    i = 0
    while i < len(toks):
        if interp.match(toks[i]):
            j = i + 1
            while j < len(toks) and toks[j].startswith("-") and toks[j] != "-c":
                j += 1
            if j < len(toks) and toks[j] == "-c":
                if j + 1 >= len(toks):
                    return None
                out.append(toks[j + 1])
                i = j + 2
                continue
        i += 1
    return out
