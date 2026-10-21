"""Compact, task-scoped repository map for prompt injection.

Input: the session's workspace directory and the task text.
Output: ``Repo map`` lines -- for a foreign workspace a one-line layout
(top-level directories, top-level files, source files counted by
extension), then up to three task-ranked files with their top-level
symbols and the tests beside them.

Two index shapes. A tree with a top-level ``delfin/`` directory (the
DELFIN checkout) keeps the original index: ``delfin/`` and ``tests/``,
symbols from ``ast``. Every other workspace is walked whole, pruned at
vendored/build/hidden directories, capped at ``_MAX_FILES`` files, and
symbols come from one regex family per language (Python still uses
``ast``). No new dependency: tree-sitter would read more languages more
exactly and is not installed.
"""

from __future__ import annotations

import ast
import os
import re
from collections import Counter
from functools import lru_cache
import json
from pathlib import Path
from typing import Any


_SOURCE_EXTS = {".py", ".md", ".json", ".yaml", ".yml"}
_IGNORE_DIRS = {
    ".git",
    ".pytest_cache",
    ".mypy_cache",
    ".ruff_cache",
    "__pycache__",
    ".venv",
    "venv",
    "node_modules",
}
# Pruned on the foreign walk on top of _IGNORE_DIRS: build output and
# vendored trees, whose files are not the project's source. Any directory
# whose name starts with "." is pruned as well (.git, .delfin, .idea, ...).
_FOREIGN_IGNORE_DIRS = _IGNORE_DIRS | {
    "target", "dist", "build", "out", "vendor", "coverage", "bin", "obj",
    "env", "site-packages", "bower_components", "Pods", "DerivedData",
}
# Source files of the languages a code repository is commonly written in.
_CODE_EXTS = {
    ".py", ".ts", ".tsx", ".js", ".jsx", ".mjs", ".cjs", ".go", ".rs",
    ".java", ".kt", ".kts", ".scala", ".c", ".h", ".cc", ".cpp", ".cxx",
    ".hpp", ".cs", ".rb", ".php", ".swift", ".m", ".sh", ".lua", ".jl",
    ".ex", ".exs", ".erl", ".hs", ".ml", ".zig", ".dart", ".vue",
    ".svelte",
}
_FOREIGN_EXTS = _CODE_EXTS | {".md", ".json", ".yaml", ".yml", ".toml"}
# A walk over a monorepo must not turn a prompt build into a stat storm.
_MAX_FILES = 3000
# Symbols are read from the head of a file only.
_MAX_SYMBOL_READ = 64_000

_PATH_HINT_RE = re.compile(
    r"[\w./-]+\.(?:py|md|json|ya?ml|toml|tsx?|jsx?|mjs|cjs|go|rs|java|kts?"
    r"|scala|c|h|cc|cpp|hpp|cs|rb|php|swift|sh|vue|svelte)\b")

_JS_SYMBOL_RES = (
    re.compile(r"^\s*(?:export\s+)?(?:default\s+)?(?:async\s+)?function\s*\*?"
               r"\s*([A-Za-z_$][\w$]*)", re.M),
    re.compile(r"^\s*(?:export\s+)?(?:default\s+)?(?:abstract\s+)?class\s+"
               r"([A-Za-z_$][\w$]*)", re.M),
    re.compile(r"^\s*export\s+(?:const|let)\s+([A-Za-z_$][\w$]*)\s*=", re.M),
    re.compile(r"^\s*export\s+(?:interface|type|enum)\s+([A-Za-z_$][\w$]*)",
               re.M),
)
_SYMBOL_RES: dict[str, tuple[re.Pattern, ...]] = {
    ".go": (
        re.compile(r"^func\s+(?:\([^)]*\)\s*)?([A-Za-z_]\w*)", re.M),
        re.compile(r"^type\s+([A-Za-z_]\w*)", re.M),
    ),
    ".rs": (
        re.compile(r"^\s*(?:pub(?:\([^)]*\))?\s+)?(?:const\s+)?(?:async\s+)?"
                   r"(?:unsafe\s+)?fn\s+([A-Za-z_]\w*)", re.M),
        re.compile(r"^\s*(?:pub(?:\([^)]*\))?\s+)?(?:struct|enum|trait)\s+"
                   r"([A-Za-z_]\w*)", re.M),
    ),
    ".rb": (re.compile(r"^\s*(?:def|class|module)\s+([A-Za-z_][\w.?!]*)",
                       re.M),),
    ".php": (re.compile(r"^\s*(?:abstract\s+|final\s+)?(?:class|function)\s+"
                        r"([A-Za-z_]\w*)", re.M),),
    ".sh": (re.compile(r"^\s*(?:function\s+)?([A-Za-z_]\w*)\s*\(\)", re.M),),
}
for _ext in (".ts", ".tsx", ".js", ".jsx", ".mjs", ".cjs", ".vue",
             ".svelte"):
    _SYMBOL_RES[_ext] = _JS_SYMBOL_RES
_JVM_LIKE = (
    re.compile(r"^\s*(?:(?:public|private|protected|internal|static|final|"
               r"abstract|open|data|sealed|partial)\s+)*(?:class|interface|"
               r"enum|record|struct|object|trait)\s+([A-Za-z_]\w*)", re.M),
    re.compile(r"^\s*(?:(?:public|private|protected|internal|static|"
               r"override|suspend|inline)\s+)*(?:fun|def|func)\s+"
               r"([A-Za-z_]\w*)", re.M),
)
for _ext in (".java", ".kt", ".kts", ".scala", ".cs", ".swift"):
    _SYMBOL_RES[_ext] = _JVM_LIKE
_C_LIKE = (re.compile(r"^\s*(?:typedef\s+)?(?:struct|class|enum)\s+"
                      r"([A-Za-z_]\w*)\s*\{", re.M),)
for _ext in (".c", ".h", ".cc", ".cpp", ".cxx", ".hpp"):
    _SYMBOL_RES[_ext] = _C_LIKE

_TEST_NAME_RE = re.compile(
    r"(?:^test_|_test\.|\.test\.|\.spec\.|_spec\.|Tests?\.(?:java|kt|cs|"
    r"swift|scala)$)", re.I)
_TEST_DIR_NAMES = {"tests", "test", "__tests__", "spec", "specs"}
_CODE_HINTS = (
    ".py",
    "pytest",
    "test_",
    "traceback",
    "stack trace",
    "function",
    "class ",
    "module",
    "code",
    "repo",
    "diff",
    "implement",
    "fix",
    "refactor",
    "edit",
    "change",
)
_STOPWORDS = {
    "the",
    "and",
    "for",
    "with",
    "this",
    "that",
    "from",
    "into",
    "your",
    "about",
    "task",
    "module",
    "file",
    "code",
    "repo",
    "tests",
    "test",
    "change",
    "changes",
    "fix",
}


def _tokenize(text: str) -> set[str]:
    return {
        token
        for token in re.split(r"[^a-z0-9]+", (text or "").lower())
        if len(token) >= 3 and token not in _STOPWORDS
    }


def _explicit_path_hints(task_text: str) -> set[str]:
    hints = set()
    for match in _PATH_HINT_RE.findall(task_text or ""):
        normalized = match.strip().lstrip("./")
        if normalized:
            hints.add(normalized.lower())
            hints.add(Path(normalized).name.lower())
    return hints


def _looks_like_code_task(task_text: str) -> bool:
    lower = (task_text or "").lower()
    if not lower.strip():
        return False
    return any(hint in lower for hint in _CODE_HINTS)


def _normalize_catalog_module(repo_dir: Path, module_name: str) -> str:
    raw = (module_name or "").strip().replace("\\", "/")
    if not raw:
        return ""
    if raw.startswith("test_"):
        candidates = [f"tests/{raw}", raw]
    else:
        candidates = [
            raw,
            f"delfin/{raw}" if not raw.startswith(("delfin/", "tests/")) else raw,
        ]
    for candidate in candidates:
        if (repo_dir / candidate).exists():
            return candidate
    return candidates[-1]


def _extract_symbols(path: Path) -> list[str]:
    if path.suffix != ".py":
        return _extract_symbols_by_regex(path)
    if not path.exists():
        return []
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except (OSError, SyntaxError, UnicodeDecodeError, ValueError):
        return []
    symbols: list[str] = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            symbols.append(node.name)
    return symbols


def _extract_symbols_by_regex(path: Path) -> list[str]:
    """Top-level definitions of a non-Python source file, in file order.

    One pattern family per language (see ``_SYMBOL_RES``); an unknown
    extension yields nothing. Reads at most ``_MAX_SYMBOL_READ`` bytes.
    """
    patterns = _SYMBOL_RES.get(path.suffix.lower())
    if not patterns:
        return []
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            text = fh.read(_MAX_SYMBOL_READ)
    except OSError:
        return []
    found: list[tuple[int, str]] = []
    for pattern in patterns:
        for m in pattern.finditer(text):
            found.append((m.start(), m.group(1)))
    out: list[str] = []
    for _, name in sorted(found):
        if name not in out:
            out.append(name)
    return out


def _is_delfin_tree(repo_dir: Path) -> bool:
    """The DELFIN checkout keeps its original, narrower index."""
    try:
        return (repo_dir / "delfin").is_dir()
    except OSError:
        return False


def _is_test_path(rel_path: str) -> bool:
    parts = rel_path.split("/")
    if any(p.lower() in _TEST_DIR_NAMES for p in parts[:-1]):
        return True
    return bool(_TEST_NAME_RE.search(parts[-1]))


def _walk_foreign(repo_dir: Path) -> list[str]:
    """Workspace-relative source files of a foreign repository.

    Pruned at hidden and vendored/build directories, hidden files
    skipped, extensions limited to ``_FOREIGN_EXTS``, at most
    ``_MAX_FILES`` paths. Symlinked directories are not followed.
    """
    out: list[str] = []
    for dirpath, dirnames, filenames in os.walk(repo_dir):
        dirnames[:] = sorted(
            d for d in dirnames
            if not d.startswith(".") and d not in _FOREIGN_IGNORE_DIRS)
        rel_dir = Path(dirpath).relative_to(repo_dir).as_posix()
        for name in sorted(filenames):
            if name.startswith("."):
                continue
            if Path(name).suffix.lower() not in _FOREIGN_EXTS:
                continue
            out.append(name if rel_dir == "." else f"{rel_dir}/{name}")
            if len(out) >= _MAX_FILES:
                return out
    return out


@lru_cache(maxsize=8)
def _build_repo_index(repo_dir_str: str) -> dict[str, dict[str, Any]]:
    repo_dir = Path(repo_dir_str)
    index: dict[str, dict[str, Any]] = {}
    if not _is_delfin_tree(repo_dir):
        for rel in _walk_foreign(repo_dir):
            index[rel] = {"symbols": _extract_symbols(repo_dir / rel)}
        return index
    for path in repo_dir.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix not in _SOURCE_EXTS:
            continue
        rel = path.relative_to(repo_dir).as_posix()
        parts = rel.split("/")
        if any(part in _IGNORE_DIRS for part in parts[:-1]):
            continue
        if not rel.startswith(("delfin/", "tests/")):
            continue
        index[rel] = {
            "symbols": _extract_symbols(path),
        }
    return index


def _catalog_test_mapping(
    repo_dir: Path,
    catalog: dict[str, Any],
) -> dict[str, list[str]]:
    raw = catalog.get("codebase_map", {}).get("test_mapping", {})
    if not isinstance(raw, dict):
        return {}
    mapped: dict[str, list[str]] = {}
    for module_name, tests in raw.items():
        normalized = _normalize_catalog_module(repo_dir, module_name)
        if not normalized:
            continue
        rendered: list[str] = []
        for test_name in tests or []:
            test_path = _normalize_catalog_module(repo_dir, test_name)
            rendered.append(test_path)
        mapped[normalized] = rendered
    return mapped


def _foreign_candidate_tests(index: dict[str, Any], rel_path: str) -> list[str]:
    """Test files in the index whose name carries the source file's stem:
    ``math.test.ts``, ``math_test.go``, ``test_math.py``, ``MathTest.java``."""
    stem = Path(rel_path).name.split(".")[0].lower()
    if not stem:
        return []
    matches = [rel for rel in sorted(index)
               if rel != rel_path and _is_test_path(rel)
               and stem in Path(rel).name.lower()]
    return matches[:3]


def _layout_line(repo_dir: Path, index: dict[str, Any]) -> str:
    """``Layout: <dirs>/ | files: <top-level files> | sources: .ext n, ...``"""
    dirs: list[str] = []
    files: list[str] = []
    try:
        for entry in sorted(repo_dir.iterdir(), key=lambda p: p.name):
            if entry.name.startswith("."):
                continue
            if entry.is_dir():
                if entry.name not in _FOREIGN_IGNORE_DIRS:
                    dirs.append(entry.name + "/")
            elif entry.is_file():
                files.append(entry.name)
    except OSError:
        pass
    counts = Counter(Path(rel).suffix.lower() for rel in index
                     if Path(rel).suffix.lower() in _CODE_EXTS)
    parts: list[str] = []
    if dirs:
        more = f" (+{len(dirs) - 12})" if len(dirs) > 12 else ""
        parts.append("Layout: " + ", ".join(dirs[:12]) + more)
    if files:
        more = f" (+{len(files) - 10})" if len(files) > 10 else ""
        parts.append("files: " + ", ".join(files[:10]) + more)
    if counts:
        parts.append("sources: " + ", ".join(
            f"{ext} {n}" for ext, n in sorted(
                counts.items(), key=lambda kv: (-kv[1], kv[0]))[:6]))
    return " | ".join(parts)


def _candidate_tests(repo_dir: Path, rel_path: str) -> list[str]:
    stem = Path(rel_path).stem
    tests_dir = repo_dir / "tests"
    if not tests_dir.exists():
        return []
    matches: list[str] = []
    for path in sorted(tests_dir.glob("test_*.py")):
        name = path.name.lower()
        if stem.lower() in name:
            matches.append(path.relative_to(repo_dir).as_posix())
    return matches[:3]


def _score_entry(
    rel_path: str,
    symbols: list[str],
    task_text: str,
    task_tokens: set[str],
    explicit_hints: set[str],
) -> int:
    lower_rel = rel_path.lower()
    name = Path(rel_path).name.lower()
    stem = Path(rel_path).stem.lower()
    score = 0
    if lower_rel in explicit_hints or name in explicit_hints:
        score += 50
    if any(hint and hint in lower_rel for hint in explicit_hints):
        score += 25
    overlap = task_tokens & _tokenize(f"{lower_rel} {stem}")
    score += len(overlap) * 5
    if stem and re.search(rf"(?<![a-z0-9_]){re.escape(stem)}(?![a-z0-9_])", task_text.lower()):
        score += 12
    symbol_hits = 0
    for symbol in symbols[:6]:
        sym = symbol.lower()
        if sym and re.search(rf"(?<![a-z0-9_]){re.escape(sym)}(?![a-z0-9_])", task_text.lower()):
            symbol_hits += 1
    score += symbol_hits * 3
    if rel_path.startswith("tests/") or _is_test_path(rel_path):
        score -= 10
    return score


def format_repo_map_context(
    repo_dir: Path,
    task_text: str,
    *,
    provider: str = "claude",
    profile_path: Path | None = None,
    playbooks_path: Path | None = None,
    max_modules: int = 3,
    max_symbols_per_module: int = 4,
) -> str:
    """Return a compact, task-scoped repo map."""
    if not _looks_like_code_task(task_text):
        return ""

    try:
        from delfin.agent.provider_profile import _load_playbook_catalog
    except Exception:
        _load_playbook_catalog = None  # type: ignore[assignment]

    repo_dir = Path(repo_dir)
    foreign = not _is_delfin_tree(repo_dir)
    # A copy: the catalog entries added below must not leak into the
    # cached index of the next call.
    index = dict(_build_repo_index(str(repo_dir)))
    task_tokens = _tokenize(task_text)
    explicit_hints = _explicit_path_hints(task_text)

    catalog: dict[str, Any] = {}
    if _load_playbook_catalog is not None:
        try:
            profile_data = None
            if profile_path is not None:
                try:
                    profile_data = json.loads(profile_path.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError):
                    profile_data = None
            catalog = _load_playbook_catalog(provider, playbooks_path, profile_data)
        except Exception:
            catalog = {}
    test_mapping = _catalog_test_mapping(repo_dir, catalog)

    for module_name in catalog.get("codebase_map", {}).get("modules", {}):
        rel = _normalize_catalog_module(repo_dir, module_name)
        # The catalog names DELFIN's modules; in a foreign workspace only
        # a file that exists there may enter the map.
        if foreign and rel and not (repo_dir / rel).is_file():
            continue
        if rel and rel not in index:
            index[rel] = {"symbols": []}

    layout = _layout_line(repo_dir, index) if foreign else ""

    ranked: list[tuple[int, str]] = []
    for rel_path, meta in index.items():
        score = _score_entry(
            rel_path,
            list(meta.get("symbols", [])),
            task_text,
            task_tokens,
            explicit_hints,
        )
        if score > 0:
            ranked.append((score, rel_path))

    if not ranked:
        return f"Repo map\n{layout}" if layout else ""

    ranked.sort(key=lambda item: (-item[0], item[1]))
    selected = [rel_path for _, rel_path in ranked
                if not rel_path.startswith("tests/")
                and not (foreign and _is_test_path(rel_path))]
    if not selected:
        selected = [rel_path for _, rel_path in ranked]
    selected = selected[:max_modules]

    lines = ["Repo map"]
    if layout:
        lines.append(layout)
    for rel_path in selected:
        meta = index.get(rel_path, {})
        symbols = list(meta.get("symbols", []))[:max_symbols_per_module]
        tests = test_mapping.get(rel_path) or (
            _foreign_candidate_tests(index, rel_path) if foreign
            else _candidate_tests(repo_dir, rel_path))
        detail_parts = []
        if symbols:
            detail_parts.append("symbols: " + ", ".join(symbols))
        if tests:
            detail_parts.append("tests: " + ", ".join(tests[:2]))
        if detail_parts:
            lines.append(f"- {rel_path} -> " + " | ".join(detail_parts))
        else:
            lines.append(f"- {rel_path}")
    return "\n".join(lines)
