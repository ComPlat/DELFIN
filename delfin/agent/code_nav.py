"""Code navigation: find_definition / find_references.

Backed by `jedi`_ for Python (in-process, no LSP server needed) and
a grep-with-heuristics fallback for everything else. Both return a
stable JSON shape::

    {
        "symbol": "compute_mean",
        "language": "python",
        "matches": [
            {"path": "mylib/stats.py", "line": 14, "col": 4,
             "kind": "function", "preview": "def compute_mean(values):"}
        ],
        "backend": "jedi"|"grep",
        "truncated": false
    }

For Python the agent gets near-LSP-quality navigation:
inheritance-aware, follows imports, distinguishes definitions from
references. For non-Python languages we fall back to a regex-grep
that's right "most of the time" and at least gives the agent
candidates to read.

.. _jedi: https://github.com/davidhalter/jedi
"""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

_MAX_MATCHES = 50
_MAX_FILES_GREP = 5000


@dataclass
class NavMatch:
    path: str
    line: int
    col: int
    kind: str
    preview: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "path": self.path, "line": self.line, "col": self.col,
            "kind": self.kind, "preview": self.preview,
        }


def _try_jedi():
    try:
        import jedi   # noqa: F401
        return jedi
    except ImportError:
        return None


def _read_line(path: Path, lineno: int) -> str:
    try:
        with path.open(encoding="utf-8", errors="replace") as f:
            for i, line in enumerate(f, start=1):
                if i == lineno:
                    return line.rstrip("\n")[:200]
                if i > lineno:
                    break
    except OSError:
        return ""
    return ""


def _files_in(workspace: Path, *, ext: tuple[str, ...] | None = None) -> list[Path]:
    out: list[Path] = []
    for path in workspace.rglob("*"):
        if not path.is_file():
            continue
        if ext is not None and path.suffix not in ext:
            continue
        # Skip vendored / cache dirs; doesn't need to be exhaustive.
        parts = set(path.parts)
        if parts & {"__pycache__", ".git", "node_modules", ".venv", "venv"}:
            continue
        out.append(path)
        if len(out) >= _MAX_FILES_GREP:
            break
    return out


def _jedi_call(
    workspace: Path, symbol: str, kind: str,
    *, file_hint: str = "", line: int | None = None,
) -> list[NavMatch]:
    """``kind`` is 'definition' or 'reference'."""
    jedi = _try_jedi()
    if jedi is None:
        return []
    matches: list[NavMatch] = []
    seen: set[tuple[str, int]] = set()
    targets: list[Path]
    if file_hint:
        target = (workspace / file_hint).resolve()
        if target.is_file():
            targets = [target]
        else:
            return []
    else:
        targets = _files_in(workspace, ext=(".py",))
    project = None
    try:
        project = jedi.Project(path=str(workspace))
    except Exception:
        project = None
    for f in targets:
        try:
            src = f.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        # Find the symbol's offset within the file (any line, any col).
        for m in re.finditer(rf"\b{re.escape(symbol)}\b", src):
            offset = m.start()
            ln = src.count("\n", 0, offset) + 1
            line_start = src.rfind("\n", 0, offset) + 1
            col = offset - line_start
            try:
                script = jedi.Script(
                    code=src, path=str(f),
                    project=project,
                ) if project else jedi.Script(code=src, path=str(f))
            except Exception:
                continue
            try:
                if kind == "definition":
                    refs = script.goto(line=ln, column=col,
                                        follow_imports=True)
                else:
                    refs = script.get_references(line=ln, column=col)
            except Exception:
                continue
            for r in refs:
                rpath = (Path(str(r.module_path))
                         if r.module_path else Path("<inline>"))
                # Only keep workspace-internal hits
                try:
                    rpath.relative_to(workspace)
                except ValueError:
                    continue
                key = (str(rpath), int(r.line))
                if key in seen:
                    continue
                seen.add(key)
                matches.append(NavMatch(
                    path=str(rpath.relative_to(workspace)),
                    line=int(r.line),
                    col=int(r.column),
                    kind=str(r.type or kind),
                    preview=_read_line(rpath, r.line),
                ))
                if len(matches) >= _MAX_MATCHES:
                    return matches
            if matches:
                # One symbol position per file usually suffices.
                break
    return matches


_DEF_RE_BY_LANG: dict[str, list[re.Pattern[str]]] = {
    ".js": [re.compile(rf"\b(function|const|let|var|class)\s+{{name}}\b")],
    ".ts": [re.compile(rf"\b(function|const|let|var|class|interface|type)\s+{{name}}\b")],
    ".go": [re.compile(rf"\bfunc\s+(\([^)]*\)\s+)?{{name}}\b")],
    ".rs": [re.compile(rf"\b(fn|struct|enum|trait|impl)\s+{{name}}\b")],
    ".java": [re.compile(rf"\b(class|interface)\s+{{name}}\b"),
              re.compile(rf"\b\w[\w<>,\s]*\s+{{name}}\s*\(")],
    ".cpp": [re.compile(rf"\b(class|struct)\s+{{name}}\b"),
             re.compile(rf"\b\w[\w*&\s]*\s+{{name}}\s*\(")],
    ".c": [re.compile(rf"\b\w[\w*\s]*\s+{{name}}\s*\(")],
}


def _grep_definition(
    workspace: Path, symbol: str,
) -> list[NavMatch]:
    """Best-effort definition hunt across non-Python languages."""
    matches: list[NavMatch] = []
    sym_re = re.compile(rf"\b{re.escape(symbol)}\b")
    for f in _files_in(workspace):
        ext = f.suffix
        patterns = _DEF_RE_BY_LANG.get(ext)
        if not patterns:
            continue
        try:
            text = f.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for ln_no, line in enumerate(text.splitlines(), start=1):
            if sym_re.search(line) and any(
                p.pattern.replace("{name}", re.escape(symbol)) and re.search(
                    p.pattern.replace("{name}", re.escape(symbol)), line,
                )
                for p in patterns
            ):
                matches.append(NavMatch(
                    path=str(f.relative_to(workspace)),
                    line=ln_no, col=line.find(symbol),
                    kind="definition?",
                    preview=line.strip()[:200],
                ))
                if len(matches) >= _MAX_MATCHES:
                    return matches
    return matches


def _grep_references(
    workspace: Path, symbol: str,
) -> list[NavMatch]:
    """Plain word-boundary grep across the workspace. Best-effort."""
    matches: list[NavMatch] = []
    sym_re = re.compile(rf"\b{re.escape(symbol)}\b")
    for f in _files_in(workspace):
        if f.suffix in {".png", ".jpg", ".gif", ".pdf", ".webp", ".bin"}:
            continue
        try:
            text = f.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for ln_no, line in enumerate(text.splitlines(), start=1):
            if sym_re.search(line):
                matches.append(NavMatch(
                    path=str(f.relative_to(workspace)),
                    line=ln_no, col=line.find(symbol),
                    kind="reference",
                    preview=line.strip()[:200],
                ))
                if len(matches) >= _MAX_MATCHES:
                    return matches
    return matches


def find_definition(
    workspace: Path | str, symbol: str, *,
    file_hint: str = "",
    language: str = "auto",
) -> dict[str, Any]:
    workspace = Path(workspace).expanduser().resolve()
    if not workspace.is_dir():
        return {"error": f"workspace not a directory: {workspace}"}
    if not symbol or not re.match(r"^[A-Za-z_][\w]*$", symbol):
        return {"error": f"symbol must be a plain identifier, got {symbol!r}"}
    backend = "grep"
    matches: list[NavMatch] = []
    use_python = (
        language in ("auto", "python")
        and (not file_hint or file_hint.endswith(".py"))
    )
    if use_python and _try_jedi() is not None:
        matches = _jedi_call(
            workspace, symbol, kind="definition", file_hint=file_hint,
        )
        backend = "jedi"
    if not matches:
        matches = _grep_definition(workspace, symbol)
        backend = "grep"
    return {
        "symbol": symbol,
        "language": "python" if backend == "jedi" else "any",
        "backend": backend,
        "matches": [m.to_dict() for m in matches[:_MAX_MATCHES]],
        "truncated": len(matches) >= _MAX_MATCHES,
    }


def find_references(
    workspace: Path | str, symbol: str, *,
    file_hint: str = "",
    language: str = "auto",
) -> dict[str, Any]:
    workspace = Path(workspace).expanduser().resolve()
    if not workspace.is_dir():
        return {"error": f"workspace not a directory: {workspace}"}
    if not symbol or not re.match(r"^[A-Za-z_][\w]*$", symbol):
        return {"error": f"symbol must be a plain identifier, got {symbol!r}"}
    backend = "grep"
    matches: list[NavMatch] = []
    use_python = (
        language in ("auto", "python")
        and (not file_hint or file_hint.endswith(".py"))
    )
    if use_python and _try_jedi() is not None:
        matches = _jedi_call(
            workspace, symbol, kind="reference", file_hint=file_hint,
        )
        backend = "jedi"
    if not matches:
        matches = _grep_references(workspace, symbol)
        backend = "grep"
    return {
        "symbol": symbol,
        "language": "python" if backend == "jedi" else "any",
        "backend": backend,
        "matches": [m.to_dict() for m in matches[:_MAX_MATCHES]],
        "truncated": len(matches) >= _MAX_MATCHES,
    }


__all__ = ["find_definition", "find_references"]
