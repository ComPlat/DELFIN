"""Project bootstrap — detect tech stack + write AGENTS.md scaffold.

When the user runs ``/init`` in a fresh project the agent shouldn't
have to guess at the build/test/lint commands every time. This module
scans the repo for the usual indicators (pyproject.toml, package.json,
Cargo.toml, Makefile, …) and writes an ``AGENTS.md`` describing the
detected commands + a small workflow recipe.

The result is read back automatically by the prompt-loader's project-
memory layer (``project_memory.AGENTS.md`` is in the auto-load list)
so the agent has the context from the very next turn.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class ProjectProfile:
    """Snapshot of detected project metadata."""
    name: str = ""
    description: str = ""
    language: str = ""
    package_manager: str = ""
    test_command: str = ""
    build_command: str = ""
    lint_command: str = ""
    run_command: str = ""
    extras: dict[str, str] = field(default_factory=dict)


# --- Detectors --------------------------------------------------------------


def _read_json(path: Path) -> dict | None:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _detect_python(root: Path) -> ProjectProfile | None:
    pyproject = root / "pyproject.toml"
    setup_py = root / "setup.py"
    if not (pyproject.is_file() or setup_py.is_file()):
        return None
    name = root.name
    desc = ""
    if pyproject.is_file():
        try:
            text = pyproject.read_text(encoding="utf-8")
            m = re.search(r'^name\s*=\s*["\']([^"\']+)["\']', text, re.MULTILINE)
            if m:
                name = m.group(1)
            d = re.search(r'^description\s*=\s*["\']([^"\']+)["\']', text, re.MULTILINE)
            if d:
                desc = d.group(1)
        except OSError:
            pass
    pm = "uv" if (root / "uv.lock").is_file() else (
        "poetry" if (root / "poetry.lock").is_file() else "pip"
    )
    return ProjectProfile(
        name=name, description=desc,
        language="Python", package_manager=pm,
        test_command="pytest -q",
        build_command=f"{pm} build" if pm in ("uv", "poetry") else "",
        lint_command="ruff check .",
        run_command=f"python -m {name.replace('-', '_')}",
    )


def _detect_node(root: Path) -> ProjectProfile | None:
    pj = root / "package.json"
    if not pj.is_file():
        return None
    data = _read_json(pj) or {}
    name = str(data.get("name") or root.name)
    desc = str(data.get("description") or "")
    scripts = data.get("scripts") or {}
    pm = "pnpm" if (root / "pnpm-lock.yaml").is_file() else (
        "yarn" if (root / "yarn.lock").is_file() else "npm"
    )
    runner = f"{pm} run" if pm == "npm" else pm
    return ProjectProfile(
        name=name, description=desc,
        language="JavaScript/TypeScript",
        package_manager=pm,
        test_command=f"{runner} test" if "test" in scripts else "",
        build_command=f"{runner} build" if "build" in scripts else "",
        lint_command=(
            f"{runner} lint" if "lint" in scripts else "eslint ."
        ),
        run_command=f"{runner} start" if "start" in scripts else "",
    )


def _detect_rust(root: Path) -> ProjectProfile | None:
    cargo = root / "Cargo.toml"
    if not cargo.is_file():
        return None
    name = root.name
    desc = ""
    try:
        text = cargo.read_text(encoding="utf-8")
        m = re.search(r'^name\s*=\s*["\']([^"\']+)["\']', text, re.MULTILINE)
        if m:
            name = m.group(1)
        d = re.search(r'^description\s*=\s*["\']([^"\']+)["\']', text, re.MULTILINE)
        if d:
            desc = d.group(1)
    except OSError:
        pass
    return ProjectProfile(
        name=name, description=desc,
        language="Rust", package_manager="cargo",
        test_command="cargo test",
        build_command="cargo build --release",
        lint_command="cargo clippy --all-targets",
        run_command="cargo run",
    )


def _detect_go(root: Path) -> ProjectProfile | None:
    if not (root / "go.mod").is_file():
        return None
    return ProjectProfile(
        name=root.name, language="Go", package_manager="go",
        test_command="go test ./...",
        build_command="go build ./...",
        lint_command="go vet ./...",
        run_command="go run .",
    )


def _detect_make(root: Path) -> dict[str, str]:
    """Surface make targets even when another detector fires too."""
    mk = root / "Makefile"
    if not mk.is_file():
        return {}
    try:
        text = mk.read_text(encoding="utf-8")
    except OSError:
        return {}
    targets: dict[str, str] = {}
    for line in text.splitlines():
        m = re.match(r"^([a-zA-Z0-9_\-]+):", line)
        if m and not line.startswith("\t"):
            targets[m.group(1)] = f"make {m.group(1)}"
    return targets


_DETECTORS = (_detect_python, _detect_node, _detect_rust, _detect_go)


def detect_project(repo_root: Path | str) -> ProjectProfile:
    """Return a populated ``ProjectProfile`` for the given repo root.

    The first language detector that matches wins; Makefile targets are
    surfaced as extras regardless. Always returns a profile (never None)
    so the caller can still write a useful AGENTS.md even for empty
    repos — just with placeholders the user can fill in.
    """
    root = Path(repo_root).expanduser().resolve()
    for det in _DETECTORS:
        profile = det(root)
        if profile is not None:
            profile.extras.update(_detect_make(root))
            return profile
    return ProjectProfile(
        name=root.name, language="(undetected)",
        extras=_detect_make(root),
    )


# --- AGENTS.md renderer -----------------------------------------------------


def render_agents_md(profile: ProjectProfile, repo_root: Path | str) -> str:
    """Compose the markdown body for AGENTS.md from a profile."""
    root = Path(repo_root).expanduser().resolve()
    lines = [
        f"# AGENTS.md — {profile.name or root.name}",
        "",
    ]
    if profile.description:
        lines.append(profile.description.strip())
        lines.append("")
    lines.extend([
        "This file orients an agent working in this repo. It is auto-",
        "loaded into the agent's system prompt by ``project_memory`` —",
        "edit freely; the agent picks up changes on the next turn.",
        "",
        "## Stack",
        "",
        f"- **Language:** {profile.language or '(set me)'}",
    ])
    if profile.package_manager:
        lines.append(f"- **Package manager:** {profile.package_manager}")
    lines.extend([
        "",
        "## Commands",
        "",
    ])
    if profile.test_command:
        lines.append(f"- **Test:** `{profile.test_command}`")
    if profile.build_command:
        lines.append(f"- **Build:** `{profile.build_command}`")
    if profile.lint_command:
        lines.append(f"- **Lint:** `{profile.lint_command}`")
    if profile.run_command:
        lines.append(f"- **Run:** `{profile.run_command}`")
    if not (profile.test_command or profile.build_command
            or profile.lint_command or profile.run_command):
        lines.append(
            "- (fill in `test:` / `build:` / `lint:` / `run:` so the "
            "agent uses the right invocation)"
        )
    if profile.extras:
        lines.extend([
            "",
            "## Makefile targets",
            "",
        ])
        for name, cmd in sorted(profile.extras.items()):
            lines.append(f"- `{cmd}`")
    lines.extend([
        "",
        "## Workflow",
        "",
        "Standard incremental loop the agent uses:",
        "",
        "1. Read the failing test / target file before editing.",
        "2. Make the smallest change that satisfies the goal.",
        "3. Run the test command above; expect it to pass.",
        "4. Commit only after the user explicitly asks.",
        "",
        "## Don't touch",
        "",
        "- (list paths/files that should be left alone unless the user",
        "  explicitly asks — credentials, generated artifacts, lock",
        "  files for unrelated language stacks, …)",
        "",
    ])
    return "\n".join(lines)


# --- Main entrypoint --------------------------------------------------------


def init_project(
    repo_root: Path | str,
    *,
    overwrite: bool = False,
) -> dict:
    """Bootstrap a project: detect stack + write AGENTS.md + scaffold
    ``.delfin/settings.json``.

    Returns a summary dict with: ``profile`` (ProjectProfile), ``files``
    (list of paths that were created), ``skipped`` (existing files we
    didn't overwrite).
    """
    root = Path(repo_root).expanduser().resolve()
    profile = detect_project(root)
    body = render_agents_md(profile, root)

    files: list[str] = []
    skipped: list[str] = []

    agents = root / "AGENTS.md"
    if agents.exists() and not overwrite:
        skipped.append(str(agents))
    else:
        agents.write_text(body, encoding="utf-8")
        files.append(str(agents))

    settings_dir = root / ".delfin"
    settings_dir.mkdir(parents=True, exist_ok=True)
    settings_file = settings_dir / "settings.json"
    if settings_file.exists() and not overwrite:
        skipped.append(str(settings_file))
    else:
        scaffold = {
            "hooks": {},
            "_meta": {
                "created_by": "delfin /init",
                "profile_language": profile.language,
            },
        }
        settings_file.write_text(
            json.dumps(scaffold, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        files.append(str(settings_file))

    return {
        "profile": profile,
        "files": files,
        "skipped": skipped,
        "agents_path": str(agents),
        "settings_path": str(settings_file),
    }


__all__ = [
    "ProjectProfile",
    "detect_project",
    "render_agents_md",
    "init_project",
]
