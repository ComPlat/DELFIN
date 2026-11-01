"""SA0 invariants: the security core of the sub-agent framework.

Each invariant goes through the REAL path, not a stand-in:
- a write outside the preset's write_globs is refused by the executor
  (api_client._write_in_scope on a KitToolPermissions built the way
  run_subagent builds it via _derive_perms), not just unadvertised;
- the main session's tool surface and every existing preset are
  unchanged by SA0 — only additions (measured against the preset
  registry's own invariants and the tool-surface budget tests);
- an expert cannot spawn; a manager cannot spawn outside
  may_delegate_to (the depth-2 guard's parse-side is here; the call
  side is the protected _execute_subagent hook).
"""
from __future__ import annotations

import dataclasses
from pathlib import Path

from delfin.agent import subagents as SA
from delfin.agent import api_client as A


# ---------------------------------------------------------------------------
# Executor refuses writes outside write_globs — real call site
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class _Perms:
    mode: str = "default"
    workspace: str = ""
    subagent_depth: int = 0
    write_allow_globs: tuple = ()
    session_allowed_tools: frozenset = frozenset()
    read_tracker: dict = None
    scope_locked: bool = False


def _child_perms(tmp_path, write_globs):
    parent = _Perms(workspace=str(tmp_path))
    return SA._derive_perms(
        parent, "default", preset_tools=(),
        preset_write_globs=write_globs)


def _kit(tmp_path, child):
    return A.KitToolPermissions(
        workspace=str(tmp_path), mode="default",
        write_allow_globs=child.write_allow_globs)


def test_write_outside_the_preset_scope_is_refused_by_the_executor(
        tmp_path):
    child = _child_perms(tmp_path, ("results/*.csv",))
    kp = _kit(tmp_path, child)
    assert A._write_in_scope(tmp_path / "results" / "table.csv", kp)
    assert not A._write_in_scope(tmp_path / "delfin" / "workflows" / "x.py",
                                 kp)
    assert not A._write_in_scope(tmp_path / "tests" / "x.py", kp)


def test_empty_preset_globs_change_nothing(tmp_path):
    """A preset without write_globs must not introduce a scope: today's
    behaviour for every existing preset."""
    child = _child_perms(tmp_path, ())
    kp = _kit(tmp_path, child)
    assert A._write_in_scope(tmp_path / "anything" / "goes.py", kp)


def test_run_subagent_passes_the_preset_globs(tmp_path, monkeypatch):
    """Wiring: run_subagent must hand the preset's write_globs to
    _derive_perms, not leave them parsed-but-unused."""
    captured = {}

    real = SA._derive_perms

    def spy(parent_perms, mode, workspace=None, preset_tools=(),
            preset_write_globs=(), preset_name=""):
        captured["globs"] = tuple(preset_write_globs)
        return real(parent_perms, mode, workspace=workspace,
                    preset_tools=preset_tools,
                    preset_write_globs=preset_write_globs,
                    preset_name=preset_name)

    monkeypatch.setattr(SA, "_derive_perms", spy)
    p = SA.SubagentPreset(
        name="wired-writer", description="d", system_prompt="You write.",
        mode="default", tools=("read_file",),
        write_globs=("results/*",))
    saved = dict(SA.SUBAGENT_PRESETS)
    SA.SUBAGENT_PRESETS.clear()
    SA.SUBAGENT_PRESETS["wired-writer"] = p
    try:
        # parent_perms=None keeps the run short: the wiring is what the
        # test pins, not the run itself.
        SA.run_subagent(
            subagent_type="wired-writer", description="d",
            prompt="Please write one results row, at least this long.",
            parent_client=type("C", (), {"_permissions": None,
                                         "model": "m"})(),
            parent_perms=None)
    finally:
        SA.SUBAGENT_PRESETS.clear()
        SA.SUBAGENT_PRESETS.update(saved)
    assert captured["globs"] == ("results/*",)


# ---------------------------------------------------------------------------
# Experts cannot spawn; whitelist parse-side
# ---------------------------------------------------------------------------

def test_every_builtin_expert_has_no_whitelist():
    for name in SA._BUILTIN_PRESETS:
        assert SA.SUBAGENT_PRESETS[name].may_delegate_to == (), name


def test_a_manager_whitelist_only_lists_known_experts(tmp_path,
                                                      monkeypatch):
    """A department manager's may_delegate_to names presets that exist;
    an unknown name in the frontmatter keeps the preset out of the
    registry (validate_preset refuses it — see
    tests/test_sa0_validate_preset.py) so the model never sees a type
    pointing at nothing."""
    problems = SA.validate_preset(
        "---\nname: broken-manager\ndescription: d\n"
        "may_delegate_to: ghost-expert\n---\nbody\n")
    assert any("ghost-expert" in p for p in problems)


# ---------------------------------------------------------------------------
# Extend, never lose: presets only gain, never lose or change
# ---------------------------------------------------------------------------

def test_builtins_keep_their_mode_and_prompt_length():
    """The built-ins' identity fields are untouched by SA0's additions —
    a rename or mode change would break every reference."""
    assert SA._BUILTIN_PRESETS["explore"].mode == "plan"
    assert SA._BUILTIN_PRESETS["plan"].mode == "plan"
    assert SA._BUILTIN_PRESETS["code-reviewer"].mode == "plan"
    assert SA._BUILTIN_PRESETS["general-purpose"].mode == "default"
    for name in SA._BUILTIN_PRESETS:
        assert SA._BUILTIN_PRESETS[name].system_prompt.strip()
        assert SA._BUILTIN_PRESETS[name].description.strip()


def test_new_fields_default_to_no_restriction_on_builtins():
    for name in SA._BUILTIN_PRESETS:
        p = SA.SUBAGENT_PRESETS[name]
        assert p.write_globs == ()
        assert p.output_schema is None
        assert p.department == ""
        assert p.enabled_setting == ""
        # An empty enabled_setting means always on: every builtin is
        # present in the live registry.
        assert name in SA.SUBAGENT_PRESETS
