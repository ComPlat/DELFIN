# Phase 2 of package SWP (nacht-s14).
#
# `write_member_settings` (delfin/agent/swarm/launch.py) must MERGE the
# member's role policy from delfin/agent/swarm/policy.py into the member's
# worktree settings file, and only ever TIGHTEN:
#
#   - never drop an existing deny_pattern (union, not replace);
#   - never widen an existing write_allow_globs (narrow-only);
#   - never write a bypass default_mode;
#   - effective policy == role policy, proven through kit_settings' own
#     merge (load on the worktree), not by re-reading the dict we wrote.
#
# Each case is RED on the current implementation (which overwrites the
# settings file with only member["owned"] + member["permission_mode"]).

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent.kit_settings import load as kit_load
from delfin.agent.swarm.launch import write_member_settings
from delfin.agent.swarm.policy import deny_patterns_for

HOME = Path("/home/alice")


def _member(*, role: str = "builder", owned: list[str] | None = None,
            permission_mode: str = "default") -> dict:
    return {
        "name": f"pkg-{role}",
        "role": role,
        "worktree": "member-worktree",
        "owned": owned or ["delfin/agent/swarm/"],
        "permission_mode": permission_mode,
    }


def _read_settings(worktree: Path) -> dict:
    with (worktree / ".delfin" / "settings.json").open("r", encoding="utf-8") as f:
        return json.load(f)


def test_merges_role_deny_into_existing_never_drops(tmp_path: Path) -> None:
    """An existing deny_pattern must survive, and the role's patterns must be
    added on top of it (unite, never replace)."""
    repo = tmp_path
    wt = repo / "member-worktree"
    (wt / ".delfin").mkdir(parents=True)
    with (wt / ".delfin" / "settings.json").open("w", encoding="utf-8") as f:
        json.dump({"kit": {"deny_patterns": ["custom-existing-pattern"]}}, f)

    write_member_settings(_member(), repo, home=HOME)

    deny = _read_settings(wt)["kit"]["deny_patterns"]
    expected = set(deny_patterns_for("builder", home=HOME))
    assert "custom-existing-pattern" in deny
    assert expected <= set(deny)


def test_never_widens_existing_write_allow_globs(tmp_path: Path) -> None:
    """A write_allow_globs that is already present must not grow. If the file
    already restricts writes to src/a.py and the member owns src/a.py AND
    docs/b.md, the result must not gain docs/b.md."""
    repo = tmp_path
    wt = repo / "member-worktree"
    (wt / ".delfin").mkdir(parents=True)
    with (wt / ".delfin" / "settings.json").open("w", encoding="utf-8") as f:
        json.dump({"kit": {"write_allow_globs": ["src/a.py"]}}, f)

    member = _member(owned=["src/a.py", "docs/b.md"])
    write_member_settings(member, repo, home=HOME)

    globs = _read_settings(wt)["kit"]["write_allow_globs"]
    assert "docs/b.md" not in globs
    assert set(globs) <= {"src/a.py"}


def test_writes_member_owned_as_write_allow_globs_when_none_exist(tmp_path: Path) -> None:
    """On a fresh member the owned scope becomes its write_allow_globs."""
    repo = tmp_path
    member = _member(owned=["delfin/agent/swarm/", "docs/"])
    write_member_settings(member, repo, home=HOME)

    wt = repo / "member-worktree"
    globs = _read_settings(wt)["kit"]["write_allow_globs"]
    assert set(globs) == {"delfin/agent/swarm/", "docs/"}


def test_never_writes_bypass_permission_mode(tmp_path: Path) -> None:
    """A requested bypass mode must be clamped, never persisted."""
    repo = tmp_path
    member = _member(permission_mode="bypassPermissions")
    write_member_settings(member, repo, home=HOME)

    wt = repo / "member-worktree"
    mode = _read_settings(wt)["kit"]["default_mode"]
    assert mode != "bypassPermissions"
    assert mode == "default"


def test_effective_policy_equals_role_policy_via_kit_settings_load(tmp_path: Path) -> None:
    """The real proof: read the member's settings back through kit_settings'
    own merge (load on the worktree) and check the effective policy equals
    the role policy. Not a re-read of the dict we just wrote."""
    repo = tmp_path
    write_member_settings(_member(role="reviewer"), repo, home=HOME)

    wt = repo / "member-worktree"
    user_path = tmp_path / "user" / "settings.json"  # empty/absent user tier
    merged = kit_load(repo_dir=wt, user_path=user_path)

    assert set(merged.deny_patterns) >= set(
        deny_patterns_for("reviewer", home=HOME)
    )
    assert set(merged.write_allow_globs) == {"delfin/agent/swarm/"}


def test_unknown_role_is_refused_before_writing(tmp_path: Path) -> None:
    """An unknown role must not silently produce an empty/partial policy."""
    repo = tmp_path
    wt = repo / "member-worktree"
    with pytest.raises(ValueError):
        write_member_settings(_member(role="tyrant"), repo, home=HOME)
    assert not (wt / ".delfin" / "settings.json").exists()
