"""swarm.yaml spec loading and validation (package SW1).

The spec is the contract the user edits. Validation errors are clear,
not tracebacks: a bad spec names the package and the glob that failed.

Rules enforced here (the operator's own swarm practice, made declarative):
- packages need disjoint write scopes: two packages that may write the
  same file race, and the later reviewer cannot tell which one moved it;
- no owned glob may cover a protected path (the safety frame is not a
  deliverable): see PROTECTED_GLOBS;
- the spec may not exceed the session limit;
- roles are a closed vocabulary: builder, reviewer, qs.
"""

from __future__ import annotations

from fnmatch import fnmatch
from pathlib import Path
from typing import Any

import yaml

SPEC_VERSION = 1

#: Roles a package member can hold. Closed on purpose: an unknown role in
#: a spec is a typo today and a permission widening tomorrow.
VALID_ROLES = frozenset({"builder", "reviewer", "qs"})

#: Path prefixes the safety frame lives in. A spec that claims one as a
#: write scope asks the agent to edit its own guard rails -- refused here,
#: at the earliest possible point, not by the permission gate mid-session.
PROTECTED_GLOBS: tuple[str, ...] = (
    "delfin/agent/sandbox.py",
    "delfin/agent/permission_rules.py",
    "delfin/agent/approval_guard.py",
    "delfin/agent/turn_outcome.py",
    "delfin/agent/api_client.py",
    "delfin/agent/engine.py",
    "delfin/agent/kit_settings.py",
    "delfin/agent/kit_confirm.py",
    "delfin/agent/terminal_confirm.py",
    "delfin/agent/file_confirm.py",
    "delfin/agent/mcp_isolation.py",
)

#: Upper bound on members per swarm. Each member is a paid session with a
#: worktree; more than this stops being a swarm and becomes a queue.
MAX_SESSIONS_LIMIT = 24


class SpecError(ValueError):
    """A swarm.yaml that cannot be started, with the reason in prose."""


def load_spec(path: Path | str) -> dict[str, Any]:
    """Read and parse a swarm.yaml file into a plain dict."""
    raw = Path(path).read_text(encoding="utf-8")
    try:
        data = yaml.safe_load(raw)
    except yaml.YAMLError as exc:
        raise SpecError(f"swarm.yaml is not valid YAML: {exc}") from exc
    if not isinstance(data, dict):
        raise SpecError("swarm.yaml must be a mapping at the top level")
    data.setdefault("spec_version", SPEC_VERSION)
    return data


def _glob_conflicts(a: str, b: str) -> bool:
    """Do two owned globs overlap?

    A pattern matches any path the other pattern's own text could name:
    the fine pattern's text is a valid concrete path for the coarse one
    (``delfin/agent/worktrees/**/sub/**`` matches ``delfin/agent/worktrees/**``),
    identical globs always conflict, and two sibling globs that differ in
    a literal segment do not.
    """
    if a == b:
        return True
    return fnmatch(b, a) or fnmatch(a, b)


def _check_protected(glob: str, pkg: str) -> None:
    for protected in PROTECTED_GLOBS:
        if _glob_conflicts(glob, protected):
            raise SpecError(
                f"package '{pkg}': owned glob '{glob}' covers the protected "
                f"path '{protected}' -- the safety frame is not a "
                f"deliverable")


def validate_spec(spec: dict[str, Any]) -> None:
    """Validate a parsed swarm spec in place; raise SpecError on any rule."""
    if not spec.get("goal"):
        raise SpecError("spec is missing 'goal'")
    packages = spec.get("packages")
    if not packages or not isinstance(packages, list):
        raise SpecError("spec is missing 'packages' (a non-empty list)")

    max_sessions = spec.get("max_sessions", MAX_SESSIONS_LIMIT)
    if not isinstance(max_sessions, int) or max_sessions < 1:
        raise SpecError("max_sessions must be a positive integer")
    if max_sessions > MAX_SESSIONS_LIMIT:
        raise SpecError(
            f"max_sessions is {max_sessions}; the limit is "
            f"{MAX_SESSIONS_LIMIT} -- each member is a paid session with "
            f"its own worktree")
    n_members = 0
    scopes: list[tuple[str, str]] = []

    for pkg in packages:
        if not isinstance(pkg, dict):
            raise SpecError("every package must be a mapping")
        name = pkg.get("name")
        if not name or not isinstance(name, str):
            raise SpecError("every package needs a 'name'")
        owned = pkg.get("owned")
        if not owned or not isinstance(owned, list):
            raise SpecError(
                f"package '{name}': missing 'owned' (a non-empty list of "
                f"write globs)")
        for glob in owned:
            if not isinstance(glob, str) or "*" not in glob and "/" not in glob:
                raise SpecError(
                    f"package '{name}': owned glob '{glob}' is not a path "
                    f"glob")
            _check_protected(glob, name)
            scopes.append((name, glob))

        roles = pkg.get("roles")
        if roles is None:
            # A package without roles is one builder: the common case.
            roles = {"builder": True}
        if not isinstance(roles, dict):
            raise SpecError(
                f"package '{name}': 'roles' must be a mapping "
                f"(role -> true)")
        for role in roles:
            if role not in VALID_ROLES:
                raise SpecError(
                    f"package '{name}': unknown role '{role}' "
                    f"(valid: {sorted(VALID_ROLES)})")
        n_members += len(roles)

    seen: dict[str, str] = {}
    for name, glob in scopes:
        for other_name, other in seen.items():
            if other_name != name and _glob_conflicts(glob, other):
                raise SpecError(
                    f"write scopes overlap: package '{other_name}' ({other}) "
                    f"and package '{name}' ({glob})")
        seen[name] = glob

    if n_members > max_sessions:
        raise SpecError(
            f"max_sessions is {max_sessions} but the spec declares "
            f"{n_members} member roles; raise max_sessions or cut packages")
