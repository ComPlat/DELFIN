from __future__ import annotations

import os
import re
import shlex
import subprocess
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence


@dataclass(frozen=True)
class SystemToolCandidate:
    path: str
    source: str
    module_name: str = ""


_MODULE_LINE_RE = re.compile(r"^[A-Za-z0-9_.+-]+(?:/[A-Za-z0-9_.+-]*)+$")
_LUA_CALL_RE = re.compile(
    r'(?P<func>prepend_path|prepend-path|setenv)\(\s*"(?P<key>[^"]+)"\s*,\s*"(?P<value>[^"]+)"'
)
_TCL_CALL_RE = re.compile(
    r'(?P<func>prepend-path|prepend_path|setenv)\s+(?P<key>[A-Za-z0-9_]+)\s+(?P<value>.+)$'
)


def _run_login_shell(command: str, *, timeout: int = 15) -> subprocess.CompletedProcess[str] | None:
    try:
        return subprocess.run(
            ["bash", "-lc", command],
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except Exception:
        return None


def _normalize_module_line(line: str) -> str:
    candidate = str(line or "").strip()
    if not candidate or candidate.startswith("Lmod has detected"):
        return ""
    if candidate.startswith("No module"):
        return ""
    if not _MODULE_LINE_RE.match(candidate):
        return ""
    return candidate.rstrip("/")


@lru_cache(maxsize=1)
def available_modules() -> tuple[str, ...]:
    result = _run_login_shell("module -t avail 2>&1", timeout=20)
    if result is None:
        return ()

    modules: list[str] = []
    seen: set[str] = set()
    for line in (result.stdout + result.stderr).splitlines():
        normalized = _normalize_module_line(line)
        if not normalized:
            continue
        lowered = normalized.lower()
        if lowered in seen:
            continue
        seen.add(lowered)
        modules.append(normalized)
    return tuple(modules)


def matching_modules(patterns: Sequence[str]) -> list[str]:
    normalized_patterns = [
        str(pattern or "").strip().lower().rstrip("/")
        for pattern in patterns
        if str(pattern or "").strip()
    ]
    if not normalized_patterns:
        return []

    matches: list[str] = []
    for module_name in available_modules():
        lowered = module_name.lower()
        if any(lowered == pattern or lowered.startswith(f"{pattern}/") for pattern in normalized_patterns):
            matches.append(module_name)
    return sorted(matches, reverse=True)


def _clean_module_value(raw_value: str) -> str:
    value = str(raw_value or "").strip()
    if not value:
        return ""
    if "#" in value:
        value = value.split("#", 1)[0].strip()
    if value.startswith("{") and value.endswith("}"):
        value = value[1:-1].strip()
    if (value.startswith('"') and value.endswith('"')) or (value.startswith("'") and value.endswith("'")):
        value = value[1:-1]
    return value.strip()


@lru_cache(maxsize=128)
def module_show_paths(module_name: str) -> tuple[tuple[str, ...], tuple[tuple[str, str], ...]]:
    quoted = shlex.quote(module_name)
    result = _run_login_shell(f"module show {quoted} 2>&1", timeout=20)
    if result is None:
        return (), ()

    path_entries: list[str] = []
    env_entries: dict[str, str] = {}
    for raw_line in (result.stdout + result.stderr).splitlines():
        line = raw_line.strip()
        if not line:
            continue

        lua_match = _LUA_CALL_RE.search(line)
        if lua_match:
            key = lua_match.group("key").strip()
            value = _clean_module_value(lua_match.group("value"))
            if key == "PATH" and value:
                path_entries.append(value)
            elif value:
                env_entries[key] = value
            continue

        tcl_match = _TCL_CALL_RE.search(line)
        if tcl_match:
            key = tcl_match.group("key").strip()
            value = _clean_module_value(tcl_match.group("value"))
            if key == "PATH" and value:
                path_entries.append(value)
            elif value:
                env_entries[key] = value

    deduped_paths: list[str] = []
    seen_paths: set[str] = set()
    for value in path_entries:
        normalized = str(Path(value).expanduser())
        if normalized in seen_paths:
            continue
        seen_paths.add(normalized)
        deduped_paths.append(normalized)

    return tuple(deduped_paths), tuple(sorted(env_entries.items()))


def _iter_binary_matches(
    base_dir: Path,
    binary_names: Sequence[str],
    *,
    max_depth: int = 6,
) -> Iterable[str]:
    if not base_dir.is_dir():
        return

    wanted = {name for name in binary_names if name}

    def _walk(current: Path, depth: int) -> Iterable[str]:
        try:
            entries = sorted(current.iterdir(), key=lambda item: item.name.lower(), reverse=True)
        except (OSError, PermissionError):
            return

        child_dirs: list[Path] = []
        for entry in entries:
            if entry.is_file() and entry.name in wanted and os.access(str(entry), os.X_OK):
                yield str(entry.resolve())
            elif entry.is_dir() and depth < max_depth:
                child_dirs.append(entry)

        prioritized = [path for path in child_dirs if path.name == "bin"]
        prioritized.extend(path for path in child_dirs if path.name != "bin")
        for child in prioritized:
            yield from _walk(child, depth + 1)

    yield from _walk(base_dir, 0)


def _module_candidate_dirs(
    module_name: str,
    module_env_hints: Sequence[str],
) -> list[Path]:
    path_entries, env_entries = module_show_paths(module_name)
    env_map = {key: value for key, value in env_entries}

    candidates: list[Path] = []
    seen: set[str] = set()

    def _add_candidate(raw_value: str) -> None:
        text = str(raw_value or "").strip()
        if not text:
            return
        path = Path(text).expanduser()
        normalized = str(path)
        if normalized in seen:
            return
        seen.add(normalized)
        candidates.append(path)

    for value in path_entries:
        _add_candidate(value)

    for env_key in module_env_hints:
        value = env_map.get(env_key)
        if value:
            _add_candidate(value)

    for env_key, value in env_map.items():
        upper_key = env_key.upper()
        if upper_key.endswith(("ROOT", "HOME", "DIR", "PATH", "BINDIR", "BIN_DIR")):
            _add_candidate(value)

    return candidates


def discover_system_tool_candidates(
    binary_names: Sequence[str],
    *,
    base_dirs: Sequence[str | Path] = (),
    module_patterns: Sequence[str] = (),
    module_env_hints: Sequence[str] = (),
    max_depth: int = 6,
) -> list[SystemToolCandidate]:
    matches: list[SystemToolCandidate] = []
    seen_paths: set[str] = set()

    def _add_match(path_value: str, *, source: str, module_name: str = "") -> None:
        normalized = str(Path(path_value).expanduser())
        if normalized in seen_paths:
            return
        seen_paths.add(normalized)
        matches.append(SystemToolCandidate(path=normalized, source=source, module_name=module_name))

    for raw_dir in base_dirs:
        base_dir = Path(str(raw_dir)).expanduser()
        for found in _iter_binary_matches(base_dir, binary_names, max_depth=max_depth):
            _add_match(found, source="system")

    for module_name in matching_modules(module_patterns):
        for candidate_dir in _module_candidate_dirs(module_name, module_env_hints):
            if candidate_dir.is_file() and candidate_dir.name in binary_names and os.access(str(candidate_dir), os.X_OK):
                _add_match(str(candidate_dir.resolve()), source="module", module_name=module_name)
                continue
            if candidate_dir.is_dir():
                for binary_name in binary_names:
                    direct = candidate_dir / binary_name
                    if direct.is_file() and os.access(str(direct), os.X_OK):
                        _add_match(str(direct.resolve()), source="module", module_name=module_name)
                for found in _iter_binary_matches(candidate_dir, binary_names, max_depth=max_depth):
                    _add_match(found, source="module", module_name=module_name)

    return matches
