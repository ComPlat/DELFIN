"""One-call project introspection.

Returns everything the agent typically reads file-by-file at the
start of a session: existing venv(s), Python version of each,
installed packages, test framework, source layout, dependency
manifests, git state. The agent calls this once and decides
which workflow to use:

  - fresh project           → create venv + install + write code
  - existing venv           → activate + add deps + write code
  - extending existing code → just write/edit files, run tests

Nothing prescriptive — the report is a snapshot, the agent picks
the path.

Output schema::

    {
        "workspace": "/path/...",
        "venvs": [
            {"path": ".venv", "python": "3.11.4",
             "interpreter": ".venv/bin/python",
             "site_packages": "...",
             "installed_pkg_count": 142},
            ...
        ],
        "python_versions_on_path": ["3.11.4", "3.10.12"],
        "manifests": {
            "pyproject.toml": {"present": true, "has_pep621": true,
                                "deps": [...]},
            "requirements.txt": {"present": true, "lines": 12},
            "setup.py": {"present": false},
            "Pipfile": {"present": false},
            "poetry.lock": {"present": false},
            "uv.lock": {"present": true},
        },
        "test_framework": "pytest" | "unittest" | "" ,
        "source_layout": "src" | "flat" | "package" | "unknown",
        "top_level_dirs": ["src", "tests", "docs"],
        "git": {"is_repo": true, "branch": "main",
                "dirty": true, "ahead": 0, "behind": 0,
                "last_commit": "abc1234 fix: foo"},
        "summary_one_liner": "Python 3.11 project with .venv, "
                              "pyproject.toml, pytest, src/ layout, "
                              "git: 2 modified files",
    }

Failures are tolerated: each section catches its own errors and
returns what it can. ``summary_one_liner`` is always populated.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any


_VENV_DIR_PATTERNS = (".venv", "venv", "env", ".env")
_VENV_DIR_GLOB = (".venv*", "venv*", "env*")
_MANIFEST_FILES = (
    "pyproject.toml", "setup.py", "setup.cfg",
    "requirements.txt", "requirements-dev.txt", "requirements_dev.txt",
    "Pipfile", "poetry.lock", "uv.lock", "tox.ini",
    "Makefile", "package.json", "Cargo.toml", "go.mod",
)
_TEST_DIR_NAMES = ("tests", "test", "specs", "spec")


def _looks_like_venv(path: Path) -> bool:
    if not path.is_dir():
        return False
    return (path / "bin" / "python").is_file() or \
           (path / "Scripts" / "python.exe").is_file() or \
           (path / "pyvenv.cfg").is_file()


def _venv_python(path: Path) -> Path:
    cand = path / "bin" / "python"
    if cand.is_file():
        return cand
    cand = path / "Scripts" / "python.exe"
    if cand.is_file():
        return cand
    return cand


def _python_version(interpreter: Path) -> str:
    try:
        out = subprocess.run(
            [str(interpreter), "-c",
             "import sys; print('%d.%d.%d' % sys.version_info[:3])"],
            capture_output=True, text=True, timeout=4,
        )
        return (out.stdout or "").strip()
    except (FileNotFoundError, subprocess.SubprocessError):
        return ""


def _installed_pkg_count(interpreter: Path) -> int:
    try:
        out = subprocess.run(
            [str(interpreter), "-m", "pip", "list", "--format=freeze"],
            capture_output=True, text=True, timeout=8,
        )
        return sum(1 for ln in (out.stdout or "").splitlines() if ln.strip())
    except (FileNotFoundError, subprocess.SubprocessError):
        return 0


def _site_packages(interpreter: Path) -> str:
    try:
        out = subprocess.run(
            [str(interpreter), "-c",
             "import sysconfig; print(sysconfig.get_paths()['purelib'])"],
            capture_output=True, text=True, timeout=4,
        )
        return (out.stdout or "").strip()
    except (FileNotFoundError, subprocess.SubprocessError):
        return ""


def _discover_venvs(workspace: Path) -> list[dict[str, Any]]:
    seen: set[Path] = set()
    out: list[dict[str, Any]] = []
    candidates: list[Path] = []
    for name in _VENV_DIR_PATTERNS:
        candidates.append(workspace / name)
    # Also pick up *.venv-* style (dotted suffix variants like .venv-myproj)
    for pattern in _VENV_DIR_GLOB:
        try:
            candidates.extend(workspace.glob(pattern))
        except OSError:
            continue
    for cand in candidates:
        if cand in seen:
            continue
        seen.add(cand)
        if not _looks_like_venv(cand):
            continue
        py = _venv_python(cand)
        out.append({
            "path": str(cand.relative_to(workspace)) if cand.is_relative_to(workspace) else str(cand),
            "interpreter": str(py),
            "python": _python_version(py),
            "site_packages": _site_packages(py),
            "installed_pkg_count": _installed_pkg_count(py),
        })
    return out


def _read_pyproject_deps(path: Path) -> tuple[bool, list[str]]:
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        return False, []
    deps: list[str] = []
    has_pep621 = "[project]" in text
    # Cheap line-based dep extraction (don't pull tomllib for a hint)
    in_deps = False
    for line in text.splitlines():
        line_stripped = line.strip()
        if line_stripped.startswith("dependencies"):
            in_deps = True
            # inline list?
            inline = re.search(r"\[(.*?)\]", line_stripped)
            if inline:
                deps.extend(_split_dep_items(inline.group(1)))
                in_deps = False
            continue
        if in_deps:
            if line_stripped.startswith("]"):
                in_deps = False
                continue
            cleaned = line_stripped.rstrip(",").strip(" \"'")
            if cleaned and not cleaned.startswith("#"):
                deps.append(cleaned)
    return has_pep621, [d for d in deps if d]


def _split_dep_items(raw: str) -> list[str]:
    return [
        item.strip().strip("\"'")
        for item in raw.split(",")
        if item.strip().strip("\"'")
    ]


def _scan_manifests(workspace: Path) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for name in _MANIFEST_FILES:
        path = workspace / name
        info: dict[str, Any] = {"present": path.is_file()}
        if not info["present"]:
            out[name] = info
            continue
        try:
            stat = path.stat()
            info["size_bytes"] = stat.st_size
        except OSError:
            pass
        if name == "pyproject.toml":
            has_pep621, deps = _read_pyproject_deps(path)
            info["has_pep621"] = has_pep621
            info["deps"] = deps[:50]
            info["deps_truncated"] = len(deps) > 50
        elif name in {"requirements.txt", "requirements-dev.txt",
                      "requirements_dev.txt"}:
            try:
                lines = [
                    ln.strip() for ln in path.read_text(
                        encoding="utf-8"
                    ).splitlines() if ln.strip() and not ln.strip().startswith("#")
                ]
                info["lines"] = len(lines)
                info["packages"] = lines[:50]
            except OSError:
                pass
        out[name] = info
    return out


def _detect_test_framework(workspace: Path, manifests: dict) -> str:
    # pytest in any manifest?
    py = manifests.get("pyproject.toml") or {}
    if any("pytest" in d.lower() for d in (py.get("deps") or [])):
        return "pytest"
    req = manifests.get("requirements.txt") or {}
    if any("pytest" in p.lower() for p in (req.get("packages") or [])):
        return "pytest"
    # Hint files
    if (workspace / "pytest.ini").is_file():
        return "pytest"
    if (workspace / "tox.ini").is_file():
        return "tox"
    # Heuristic: a tests/ directory with files that import pytest
    for tdir in _TEST_DIR_NAMES:
        d = workspace / tdir
        if d.is_dir():
            for f in d.glob("test_*.py"):
                try:
                    if "pytest" in f.read_text(encoding="utf-8", errors="replace"):
                        return "pytest"
                    if "unittest" in f.read_text(encoding="utf-8", errors="replace"):
                        return "unittest"
                except OSError:
                    continue
            return "pytest"   # default for repos with a tests dir
    return ""


def _detect_source_layout(workspace: Path) -> tuple[str, list[str]]:
    top: list[str] = []
    try:
        for child in workspace.iterdir():
            if child.is_dir() and not child.name.startswith(".") \
                    and child.name not in {"node_modules", "__pycache__"}:
                top.append(child.name)
    except OSError:
        return "unknown", []
    top.sort()
    if "src" in top:
        return "src", top
    # Heuristic: a top-level dir with __init__.py = "package layout"
    for d in top:
        if d in _TEST_DIR_NAMES:
            continue
        if (workspace / d / "__init__.py").is_file():
            return "package", top
    if any((workspace / d).is_file() for d in top if d.endswith(".py")):
        return "flat", top
    if any(p.suffix == ".py" for p in workspace.glob("*.py")):
        return "flat", top
    return "unknown", top


def _git_state(workspace: Path) -> dict[str, Any]:
    info: dict[str, Any] = {"is_repo": False}
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--git-dir"],
            cwd=str(workspace), capture_output=True, text=True, timeout=4,
        )
        if out.returncode != 0:
            return info
    except (FileNotFoundError, subprocess.SubprocessError):
        return info
    info["is_repo"] = True
    try:
        info["branch"] = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            cwd=str(workspace), text=True, timeout=4,
        ).strip()
    except subprocess.SubprocessError:
        pass
    try:
        st = subprocess.check_output(
            ["git", "status", "--porcelain"],
            cwd=str(workspace), text=True, timeout=4,
        )
        modified = [ln for ln in st.splitlines() if ln.strip()]
        info["dirty"] = bool(modified)
        info["modified_count"] = len(modified)
    except subprocess.SubprocessError:
        pass
    try:
        info["last_commit"] = subprocess.check_output(
            ["git", "log", "-1", "--oneline"],
            cwd=str(workspace), text=True, timeout=4,
        ).strip()[:120]
    except subprocess.SubprocessError:
        pass
    try:
        ab = subprocess.check_output(
            ["git", "rev-list", "--left-right", "--count", "HEAD...@{u}"],
            cwd=str(workspace), text=True, timeout=4,
            stderr=subprocess.DEVNULL,
        ).strip().split()
        if len(ab) == 2:
            info["ahead"] = int(ab[0])
            info["behind"] = int(ab[1])
    except subprocess.SubprocessError:
        pass
    return info


def _python_versions_on_path() -> list[str]:
    out: list[str] = []
    for cmd in ("python3", "python"):
        try:
            res = subprocess.run(
                [cmd, "-c",
                 "import sys; print('%d.%d.%d' % sys.version_info[:3])"],
                capture_output=True, text=True, timeout=3,
            )
            if res.returncode == 0 and res.stdout.strip():
                out.append(res.stdout.strip())
        except (FileNotFoundError, subprocess.SubprocessError):
            continue
    return list(dict.fromkeys(out))   # dedupe preserving order


def _summary(report: dict[str, Any]) -> str:
    parts: list[str] = []
    venvs = report.get("venvs", [])
    if venvs:
        v0 = venvs[0]
        parts.append(f"Python {v0.get('python', '?')} venv "
                     f"at {v0.get('path')} ({v0.get('installed_pkg_count', '?')} pkgs)")
    else:
        parts.append("no venv detected")
    manifests = report.get("manifests", {})
    if (manifests.get("pyproject.toml") or {}).get("present"):
        parts.append("pyproject.toml")
    if (manifests.get("requirements.txt") or {}).get("present"):
        parts.append("requirements.txt")
    tf = report.get("test_framework")
    if tf:
        parts.append(tf)
    layout = report.get("source_layout")
    if layout and layout != "unknown":
        parts.append(f"{layout} layout")
    git = report.get("git", {})
    if git.get("is_repo"):
        s = f"git:{git.get('branch', '?')}"
        if git.get("dirty"):
            s += f"({git.get('modified_count', '?')} modified)"
        parts.append(s)
    return ", ".join(parts) or "empty workspace"


def introspect(workspace: Path | str) -> dict[str, Any]:
    """Return a one-shot snapshot of the project's state."""
    workspace = Path(workspace).expanduser().resolve()
    if not workspace.is_dir():
        return {"error": f"workspace not a directory: {workspace}"}

    venvs = _discover_venvs(workspace)
    manifests = _scan_manifests(workspace)
    test_framework = _detect_test_framework(workspace, manifests)
    layout, top_dirs = _detect_source_layout(workspace)
    git = _git_state(workspace)

    report: dict[str, Any] = {
        "workspace": str(workspace),
        "venvs": venvs,
        "python_versions_on_path": _python_versions_on_path(),
        "manifests": manifests,
        "test_framework": test_framework,
        "source_layout": layout,
        "top_level_dirs": top_dirs,
        "git": git,
    }
    report["summary_one_liner"] = _summary(report)
    return report


__all__ = ["introspect"]
