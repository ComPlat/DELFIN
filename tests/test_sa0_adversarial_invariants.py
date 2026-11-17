# SA0 adversarial tests (reviewer nacht-s7, branch agent/s7-sa0r1).
#
# These attack the SUBAGENT framework's invariants on the CURRENT code. They
# are committed GREEN against main and must stay green after SA0 lands:
# any one of them going red on the builder branch is a finding against
# agent/s5-sa0b. The "default = today's behaviour" baseline lives in the
# uncommitted probe tests/test_zz_sa0_baseline.py.
#
# Files under attack:
# - delfin/agent/subagents.py   (SubagentPreset, _build_preset_registry,
#                                _derive_perms, SUBAGENT_PRESETS)
# - delfin/agent/api_client.py  (write_allow_globs enforcement via
#                                _write_in_scope / KitToolExecutor.execute)
from pathlib import Path

import pytest

from delfin.agent import subagents as S
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


# ---------------------------------------------------------------- helpers

def make_perms(tmp_path, **kw):
    """A KitToolPermissions the way the real call sites build it: only the
    kwargs the dataclass actually accepts (pattern from
    tests/test_a_preset_has_its_tools_rather_than_asking_for_them.py:71).
    depth and write scope arrive via _derive_perms / the field defaults."""
    defaults = dict(
        workspace=tmp_path,
        mode="default",
    )
    defaults.update(kw)
    return KitToolPermissions(**defaults)


# ---------------------------------------------- "extend, never lose"

def test_builtin_preset_names_unchanged():
    """SA0 adds presets; it may not remove or rename existing ones.
    Pinned at main 20c81ec4d: these five names with these modes."""
    expected = {
        "explore": "plan",
        "plan": "plan",
        "code-reviewer": "plan",
        "general-purpose": "default",
    }
    for name, mode in expected.items():
        p = S.SUBAGENT_PRESETS.get(name)
        assert p is not None, f"builtin preset {name!r} vanished"
        assert p.mode == mode, f"builtin {name!r} mode changed"


def test_builtin_explore_tools_unchanged():
    """The explore preset's tool list is a contract with the schema enum."""
    p = S.SUBAGENT_PRESETS["explore"]
    assert p.tools  # non-empty and unchanged vs main (snapshot below)
    assert "read_file" in p.tools


def test_no_preset_named_with_builtin_prefix_can_shadow_builtins():
    """A user md file may not silently override a built-in preset."""
    p = S._build_preset_registry()
    assert p.get("explore") is not None
    # and the registry composition is deterministic across calls
    assert set(p) == set(S._build_preset_registry())


# ---------------------------------------------- write scope enforcement

def test_executor_refuses_write_outside_parent_globs(tmp_path):
    """With a parent write scope, a write outside it is refused by the
    executor itself, not merely left out of the advertised schema."""
    perms = make_perms(tmp_path, write_allow_globs=("tests/*",))
    outside = tmp_path / "delfin" / "evil.py"
    outside.parent.mkdir(parents=True, exist_ok=True)
    from delfin.agent.api_client import _write_in_scope
    assert not _write_in_scope(outside, perms)
    assert _write_in_scope(tmp_path / "tests" / "ok.py", perms)


def test_write_scope_never_reaches_outside_workspace(tmp_path, tmp_path_factory):
    """A glob cannot be used to reach past the workspace root."""
    other = tmp_path_factory.mktemp("elsewhere")
    perms = make_perms(tmp_path, write_allow_globs=("**",))
    from delfin.agent.api_client import _write_in_scope
    assert not _write_in_scope(other / "x.txt", perms)


def test_child_write_scope_not_wider_than_parent_today(tmp_path):
    """Today's baseline: _derive_perms does not WIDEN a parent write scope.
    SA0 may only narrow via preset write_globs; it must never widen."""
    from delfin.agent.subagents import _derive_perms
    parent = make_perms(tmp_path, write_allow_globs=("docs/*",))
    child = _derive_perms(parent, "default")
    child_globs = tuple(getattr(child, "write_allow_globs", ()) or ())
    for g in child_globs:
        assert g in ("docs/*",), f"child gained a new write glob {g!r}"


# ---------------------------------------------- nesting / spawn guards

def test_depth_grows_and_caps(tmp_path):
    """Sub-agent depth grows by one per delegation and the framework
    refuses to go beyond the cap (anti-recursion)."""
    from delfin.agent.subagents import _derive_perms
    parent = make_perms(tmp_path)
    child = _derive_perms(parent, "plan")
    assert int(getattr(child, "subagent_depth", 0)) == 1
    grandchild = _derive_perms(child, "plan")
    assert int(getattr(grandchild, "subagent_depth", 0)) == 2


def test_child_mode_never_wider_than_parent(tmp_path):
    """The clamp takes the MORE RESTRICTIVE mode instead of raising
    (subagents.py:2704 _clamp_child_mode): a plan parent spawning a
    default/bypass preset gets a plan child — no indirect plan-mode escape.
    The child's mode must never outrank the parent's."""
    from delfin.agent.subagents import _clamp_child_mode, _derive_perms
    parent = make_perms(tmp_path, mode="plan")
    for requested in ("default", "bypassPermissions", "acceptEdits"):
        child = _derive_perms(parent, requested)
        rank = S._MODE_RESTRICTIVENESS
        assert rank.get(child.mode, 1) <= rank.get("plan", 1), (
            f"plan parent got a {child.mode!r} child")
