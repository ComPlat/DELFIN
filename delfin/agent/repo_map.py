"""Compact, task-scoped repository map for prompt injection."""

from __future__ import annotations

import ast
import re
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
_PATH_HINT_RE = re.compile(r"[\w./-]+\.(?:py|md|json|ya?ml)")
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
    if path.suffix != ".py" or not path.exists():
        return []
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except (OSError, SyntaxError, UnicodeDecodeError):
        return []
    symbols: list[str] = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            symbols.append(node.name)
    return symbols


@lru_cache(maxsize=8)
def _build_repo_index(repo_dir_str: str) -> dict[str, dict[str, Any]]:
    repo_dir = Path(repo_dir_str)
    index: dict[str, dict[str, Any]] = {}
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
    if rel_path.startswith("tests/"):
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
    index = _build_repo_index(str(repo_dir))
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
        if rel and rel not in index:
            index[rel] = {"symbols": []}

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
        return ""

    ranked.sort(key=lambda item: (-item[0], item[1]))
    selected = [rel_path for _, rel_path in ranked if not rel_path.startswith("tests/")]
    if not selected:
        selected = [rel_path for _, rel_path in ranked]
    selected = selected[:max_modules]

    lines = ["Repo map"]
    for rel_path in selected:
        meta = index.get(rel_path, {})
        symbols = list(meta.get("symbols", []))[:max_symbols_per_module]
        tests = test_mapping.get(rel_path) or _candidate_tests(repo_dir, rel_path)
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
