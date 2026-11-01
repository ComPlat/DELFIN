"""swarm launch: build the per-member launch plan (package SW1).

Read-only planning lives here; the actual start is a side effect and
lives in the same module behind ``execute=True`` on ``start_members`` --
one function per verb, JSON-able in and out, so the CLI, the SW2
dashboard and the SW3 lead call the same code.

The plan is the whole truth a member starts from: its own worktree and
branch (the existing worktree machinery), the write scope from the spec
(enforced by the existing permission gate -- nothing new is bypassed),
a permission mode that is never bypass, the role prompt built from the
built-in templates, and the task text. A dry run prints exactly this.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from typing import Any

from . import prompts as _prompts
from .spec import validate_spec

#: Directory the member worktrees are created under, relative to the
#: repo root. ``.swarm`` is git-ignored material, like node_modules.
WORKTREE_DIR = ".swarm"

ROLE_SUFFIX = {"builder": "builder", "reviewer": "reviewer", "qs": "qs"}


class LaunchError(RuntimeError):
    """A swarm cannot be started, with the reason in prose."""


def _git(repo: Path, *args: str) -> str:
    result = subprocess.run(["git", *args], cwd=repo, capture_output=True,
                            text=True)
    if result.returncode != 0:
        raise LaunchError(f"git {' '.join(args)} failed: "
                          f"{result.stderr.strip()}")
    return result.stdout


def plan_members(spec: dict[str, Any]) -> list[dict[str, Any]]:
    """One entry per member session, with everything start needs.

    JSON-able; no state changes. Branch and worktree names are derived
    deterministically from package and role so two planning runs agree.
    """
    validate_spec(spec)
    members: list[dict[str, Any]] = []
    for pkg in spec["packages"]:
        roles = pkg.get("roles") or {"builder": True}
        for role in roles:
            name = f"{pkg['name']}-{role}"
            branch = f"swarm/{pkg['name']}-{role}"
            worktree = f"{WORKTREE_DIR}/{pkg['name']}-{role}"
            members.append({
                "name": name,
                "package": pkg["name"],
                "role": role,
                "branch": branch,
                "worktree": worktree,
                "owned": list(pkg["owned"]),
                "task": pkg.get("task", ""),
                "goal": spec["goal"],
                "permission_mode": "default",
                "role_prompt": f"delfin/agent/swarm/prompts/{role}.md",
                "prompt_text": _prompts.render(role, {
                    "goal": spec["goal"],
                    "package": pkg["name"],
                    "task": pkg.get("task", ""),
                    "owned": ", ".join(pkg["owned"]),
                    "branch": branch,
                }),
                "memory_gate": True,
            })
    return members


def plan_launch(spec: dict[str, Any]) -> dict[str, Any]:
    """The full launch plan: members + start order + what the user sees."""
    members = plan_members(spec)
    return {
        "goal": spec["goal"],
        "members": members,
        "start_spacing_s": spec.get("start_spacing_s", 15),
        "memory_gate": spec.get("memory_gate", True),
        "dry_run": True,
    }


def start_members(spec: dict[str, Any], *, execute: bool = False) -> dict:
    """Start the swarm: one worktree + branch per member.

    The default refuses: starting paid sessions must be a deliberate act
    (``execute=True``), and the dry-run plan is what a user confirms.
    Real session processes are started by the caller (the CLI's launch
    path); this function creates the worktrees, branches and per-member
    settings, then returns the launch checklist.
    """
    plan = plan_launch(spec)
    if not execute:
        raise LaunchError(
            "refusing to start without execute=True -- pass --dry-run to "
            "see the plan, or execute to start it")
    repo = Path.cwd()
    started = []
    for member in plan["members"]:
        wt_path = repo / member["worktree"]
        if not wt_path.exists():
            _git(repo, "worktree", "add", "-b", member["branch"],
                 str(wt_path), "HEAD")
        write_member_settings(member, repo)
        started.append(member)
        from . import status as _status
        _status.emit_event({"kind": "phase_change",
                            "package": member["package"],
                            "member": member["name"],
                            "phase": "building"})
    return {"dry_run": False, "goal": plan["goal"], "started": started}


def render_dry_run(plan: dict[str, Any]) -> str:
    """The human-readable dry-run text: exactly what would start."""
    lines = [f"swarm: {plan['goal']}",
             f"{len(plan['members'])} member(s), spaced "
             f"{plan['start_spacing_s']}s, memory gate "
             f"{'on' if plan['memory_gate'] else 'off'}:",
             ""]
    for m in plan["members"]:
        lines.append(
            f"  {m['name']}: role={m['role']} branch={m['branch']} "
            f"worktree={m['worktree']}")
        lines.append(f"    writes: {m['owned']}")
        lines.append(f"    permission_mode: {m['permission_mode']} "
                     f"(never bypass)")
        lines.append(f"    task: {m['task']}")
        lines.append("")
    lines.append("Nothing started (dry run). Pass --execute to start.")
    return "\n".join(lines)


EXAMPLE_SPEC = """\
# swarm.yaml -- what agents run in this swarm. Edit, then:
#   delfin-agent swarm plan "<goal>"   # propose a spec draft from your tree
#   delfin-agent swarm start --dry-run # see exactly what would start
spec_version: 1

goal: "Ship the feature in three packages"

limits:
  max_sessions: 8        # every member is a paid session with a worktree
  memory_gate: true      # each start passes the memory gate first
  start_spacing_s: 15    # starts are spaced, not stampeded

packages:
  - name: core
    owned: ["delfin/agent/swarm/**"]   # write scope, enforced by the gate
    task: |
      Build the swarm core package. Red test first; commit per thought.
    roles:
      builder: true
      reviewer: true
  - name: ui
    owned: ["delfin/agent/dashboard/**"]
    task: |
      Build the dashboard view on the status API.
    roles:
      builder: true
"""


def write_example_spec(path: Path | str, *, force: bool = False) -> Path:
    """Write the commented example swarm.yaml; refuse to overwrite."""
    target = Path(path)
    if target.exists() and not force:
        raise LaunchError(f"{target} already exists -- not overwriting")
    target.write_text(EXAMPLE_SPEC, encoding="utf-8")
    return target


def write_member_settings(member: dict[str, Any], repo_dir: Path) -> Path:
    """Persist the member's write scope into ITS repo settings.

    Writes only to the member's own worktree settings file
    (``<worktree>/.delfin/settings.json``), never the user's. This is
    how the write scope from the spec reaches the existing permission
    gate: the gate reads ``kit`` write_allow_globs, so we record them
    where it looks.
    """
    from .. import kit_settings as _ks

    wt = repo_dir / member["worktree"]
    settings_path = wt / ".delfin" / "settings.json"
    settings_path.parent.mkdir(parents=True, exist_ok=True)
    settings = {
        "kit": {
            "write_allow_globs": list(member["owned"]),
            "default_mode": member["permission_mode"],
        },
    }
    settings_path.write_text(json.dumps(settings, indent=2),
                             encoding="utf-8")
    return settings_path
