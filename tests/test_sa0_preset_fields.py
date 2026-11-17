"""SA0 phase 1: new preset fields, default = today's behaviour.

Input: SubagentPreset plus the markdown-preset loader, _derive_perms and
run_subagent's output_schema handling. Output: parsed fields and the
refusals the wave-27 plan demands.

Semantics under test:
- ``write_globs`` narrows the child's write scope: the child perms carry
  them, and globs that cannot be PROVEN to lie inside the parent's scope
  refuse the delegation (never widen).
- ``output_schema`` from the preset is the default when the caller passes
  none; an explicit caller schema still wins.
- ``display_name`` is shown but ``name`` stays the stable key.
- ``enabled_setting`` gates a department: off => not in the registry =>
  not offered, refused when called (registry-level; the call-side refusal
  is api_client's existing unknown-subagent_type path).
- ``department`` / ``may_delegate_to`` are parsed; experts delegate to
  nobody and depth > 1 only along the whitelist (depth-2 enforcement
  itself lives in the protected _execute_subagent hook).
- "extend, never lose": the main tool list and all existing presets are
  unchanged by SA0 — only additions.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from delfin.agent import subagents as SA


# ---------------------------------------------------------------------------
# Parsed fields
# ---------------------------------------------------------------------------

def _preset(**over):
    kw = dict(name="xyz-editor", description="edits xyz files",
              system_prompt="You edit xyz.", mode="default",
              tools=("read_file", "edit_file"))
    kw.update(over)
    return SA.SubagentPreset(**kw)


def test_preset_fields_default_to_todays_behaviour():
    p = _preset()
    assert p.write_globs == ()
    assert p.output_schema is None
    assert p.department == ""
    assert p.may_delegate_to == ()
    assert p.display_name == ""
    assert p.enabled_setting == ""


def test_write_globs_output_schema_and_names_parse_from_frontmatter(
        tmp_path, monkeypatch):
    md = tmp_path / "xyz-editor_subagent.md"
    md.write_text(
        "---\n"
        "name: xyz-editor\n"
        "description: edits xyz files\n"
        "write_globs: delfin/workflows/xyz/*\n"
        "department: structure\n"
        "may_delegate_to: []\n"
        "display_name: Mole\n"
        "enabled_setting: agent.departments.structure\n"
        'output_schema: {"type": "object", "properties": '
        '{"verdict": {"type": "string"}}, "required": ["verdict"]}\n'
        "---\n"
        "You edit xyz files.\n",
        encoding="utf-8")
    monkeypatch.setattr(SA, "_md_preset_search_dirs", lambda: [tmp_path])
    monkeypatch.setattr(SA, "_BUILTIN_PRESETS", {})
    monkeypatch.setattr(SA, "_setting_on", lambda key: True)
    reg = SA._build_preset_registry()
    p = reg["xyz-editor"]
    assert p.write_globs == ("delfin/workflows/xyz/*",)
    assert p.department == "structure"
    assert p.may_delegate_to == ()
    assert p.display_name == "Mole"
    assert p.enabled_setting == "agent.departments.structure"
    assert p.output_schema == {
        "type": "object",
        "properties": {"verdict": {"type": "string"}},
        "required": ["verdict"],
    }
    # The stable key is the name, never the display name.
    assert "xyz-editor" in reg and "Mole" not in reg


# ---------------------------------------------------------------------------
# write_globs in _derive_perms
# ---------------------------------------------------------------------------

@dataclass
class _Perms:
    mode: str = "default"
    workspace: str = ""
    subagent_depth: int = 0
    write_allow_globs: tuple = field(default=())
    session_allowed_tools: frozenset = field(default=frozenset())
    read_tracker: dict = field(default_factory=dict)
    scope_locked: bool = False


def test_child_inherits_preset_write_globs(tmp_path):
    parent = _Perms(workspace=str(tmp_path))
    child = SA._derive_perms(
        parent, "default", preset_tools=(),
        preset_write_globs=("tests/*.py",))
    assert child.write_allow_globs == ("tests/*.py",)
    # The parent's own scope is untouched.
    assert parent.write_allow_globs == ()


def test_a_glob_outside_the_parent_scope_refuses(tmp_path):
    parent = _Perms(
        workspace=str(tmp_path), write_allow_globs=("tests/*.py",))
    with pytest.raises(SA.SubagentPermissionsError) as exc:
        SA._derive_perms(
            parent, "default", preset_tools=(),
            preset_write_globs=("delfin/**",))
    assert "write" in str(exc.value).lower()


def test_a_glob_inside_the_parent_scope_passes(tmp_path):
    parent = _Perms(
        workspace=str(tmp_path), write_allow_globs=("tests/*.py",))
    child = SA._derive_perms(
        parent, "default", preset_tools=(),
        preset_write_globs=("tests/test_sa0*.py",))
    assert child.write_allow_globs == ("tests/test_sa0*.py",)


def test_scope_proof_matches_a_real_path(tmp_path):
    """The refusal must be provable, not a string heuristic: the same
    executor regex rule must accept what a real path would match."""
    parent = _Perms(
        workspace=str(tmp_path), write_allow_globs=("tests/*.py",))
    # ``tests/*.py`` (executor: * = [^/]*) does NOT contain a deeper
    # directory glob: no glob can be proven contained there, so this
    # delegation refuses.
    with pytest.raises(SA.SubagentPermissionsError):
        SA._derive_perms(
            parent, "default", preset_tools=(),
            preset_write_globs=("tests/a/*.py",))
    # A single-segment child glob IS contained.
    child = SA._derive_perms(
        parent, "default", preset_tools=(),
        preset_write_globs=("tests/test_sa0*.py",))
    assert child.write_allow_globs == ("tests/test_sa0*.py",)
    # And the child's scope genuinely contains the file it needs...
    from delfin.agent import api_client as A
    kp = A.KitToolPermissions(
        workspace=str(tmp_path), mode="default",
        write_allow_globs=child.write_allow_globs)
    assert A._write_in_scope(tmp_path / "tests" / "test_sa0_x.py", kp)
    # ...and refuses one the parent could have written.
    assert not A._write_in_scope(tmp_path / "delfin" / "x.py", kp)
    # A parent glob with ** contains deep child globs.
    parent2 = _Perms(
        workspace=str(tmp_path), write_allow_globs=("tests/**",))
    child2 = SA._derive_perms(
        parent2, "default", preset_tools=(),
        preset_write_globs=("tests/a/*.py",))
    assert child2.write_allow_globs == ("tests/a/*.py",)


def test_tests_star_cannot_smuggle_a_deeper_directory(tmp_path):
    """The reviewer's attack: plain fnmatch would accept ``tests*`` under
    ``tests/*`` because fnmatch's ``*`` crosses slashes; the executor's
    does not, so containment must be refused."""
    parent = _Perms(
        workspace=str(tmp_path), write_allow_globs=("tests/*.py",))
    with pytest.raises(SA.SubagentPermissionsError):
        SA._derive_perms(
            parent, "default", preset_tools=(),
            preset_write_globs=("tests_evil/*.py",))


# ---------------------------------------------------------------------------
# output_schema default from the preset
# ---------------------------------------------------------------------------

_DEFAULT_SCHEMA = {
    "type": "object",
    "properties": {
        "verdict": {"type": "string"},
        "evidence": {"type": "array", "items": {"type": "string"}},
        "files_written": {"type": "array", "items": {"type": "string"}},
        "open_questions": {"type": "array", "items": {"type": "string"}},
    },
    "required": ["verdict", "evidence"],
}


def _client_round(text: str):
    """A fake client whose one streamed round emits *text*."""
    from types import SimpleNamespace

    def _fake_event(**kwargs):
        return SimpleNamespace(
            type=kwargs.get("type", "text_delta"),
            text=kwargs.get("text", ""),
            tool_name=kwargs.get("tool_name", ""),
            tool_input=kwargs.get("tool_input", ""),
            tool_output=kwargs.get("tool_output", ""),
            input_tokens=kwargs.get("input_tokens", 0),
            output_tokens=kwargs.get("output_tokens", 0),
        )

    client = MagicMock()
    client.stream_message = MagicMock(side_effect=[iter([
        _fake_event(type="text_delta", text=text),
        _fake_event(type="message_delta", input_tokens=10,
                    output_tokens=5),
    ])])
    client._permissions = None
    client.set_permissions = MagicMock()
    return client


def test_preset_schema_used_when_the_caller_passes_none():
    good = json.dumps({"verdict": "ok", "evidence": ["a.py"]})
    client = _client_round(good)
    p = _preset(output_schema=_DEFAULT_SCHEMA)
    with _registry_hook({"xyz-editor": p}):
        res = SA.run_subagent(
            subagent_type="xyz-editor", description="d",
            prompt="Report the verdict as JSON, at least this long.",
            parent_client=client, parent_perms=None)
    assert res.error == ""
    assert res.structured_output is not None
    assert "verdict" in res.structured_output


def test_caller_schema_wins_over_the_preset_default():
    caller_schema = {"type": "object",
                     "properties": {"n": {"type": "integer"}},
                     "required": ["n"]}
    good = json.dumps({"n": 7})
    client = _client_round(good)
    p = _preset(output_schema=_DEFAULT_SCHEMA)
    with _registry_hook({"xyz-editor": p}):
        res = SA.run_subagent(
            subagent_type="xyz-editor", description="d",
            prompt="Return n as JSON, at least this long please.",
            parent_client=client, parent_perms=None,
            output_schema=caller_schema)
    assert res.error == ""
    assert res.structured_output is not None
    assert "n" in res.structured_output


class _registry_hook:
    """Point run_subagent's lookup at a test registry."""

    def __init__(self, reg):
        self._reg = reg

    def __enter__(self):
        self._saved = dict(SA.SUBAGENT_PRESETS)
        SA.SUBAGENT_PRESETS.clear()
        SA.SUBAGENT_PRESETS.update(self._reg)
        return self

    def __exit__(self, *exc):
        SA.SUBAGENT_PRESETS.clear()
        SA.SUBAGENT_PRESETS.update(self._saved)
        return False


# ---------------------------------------------------------------------------
# Departments: enabled_setting, department, may_delegate_to
# ---------------------------------------------------------------------------

def test_a_department_whose_setting_is_off_is_not_offered(tmp_path,
                                                          monkeypatch):
    md = tmp_path / "xyz-editor_subagent.md"
    md.write_text(
        "---\n"
        "name: xyz-editor\n"
        "description: edits xyz files\n"
        "enabled_setting: agent.departments.structure\n"
        "---\n"
        "You edit xyz files.\n",
        encoding="utf-8")
    monkeypatch.setattr(SA, "_md_preset_search_dirs", lambda: [tmp_path])
    monkeypatch.setattr(SA, "_BUILTIN_PRESETS", {})
    # Setting off => absent.
    monkeypatch.setattr(SA, "_setting_on", lambda key: False)
    reg = SA._build_preset_registry()
    assert "xyz-editor" not in reg
    # Setting on => present.
    monkeypatch.setattr(SA, "_setting_on", lambda key: True)
    reg = SA._build_preset_registry()
    assert "xyz-editor" in reg


def test_preset_without_enabled_setting_is_always_on(tmp_path, monkeypatch):
    md = tmp_path / "plain_subagent.md"
    md.write_text(
        "---\nname: plain\ndescription: d\n---\nbody\n", encoding="utf-8")
    monkeypatch.setattr(SA, "_md_preset_search_dirs", lambda: [tmp_path])
    monkeypatch.setattr(SA, "_BUILTIN_PRESETS", {})
    monkeypatch.setattr(SA, "_setting_on", lambda key: False)
    reg = SA._build_preset_registry()
    assert "plain" in reg


def test_may_delegate_to_is_parsed(tmp_path, monkeypatch):
    md = tmp_path / "structure-manager_subagent.md"
    md.write_text(
        "---\n"
        "name: structure-manager\n"
        "description: department manager\n"
        "department: structure\n"
        "may_delegate_to: xyz-editor, smiles-builder\n"
        "---\n"
        "You plan structure work.\n",
        encoding="utf-8")
    monkeypatch.setattr(SA, "_md_preset_search_dirs", lambda: [tmp_path])
    monkeypatch.setattr(SA, "_BUILTIN_PRESETS", {})
    monkeypatch.setattr(SA, "_setting_on", lambda key: True)
    reg = SA._build_preset_registry()
    assert reg["structure-manager"].may_delegate_to == (
        "xyz-editor", "smiles-builder")
    assert reg["structure-manager"].department == "structure"


def test_experts_delegate_to_nobody():
    """The whitelist decides depth > 1: an expert without the field may
    not spawn at all. This is the subagents-side half of the rule; the
    call-side guard lives in the protected _execute_subagent hook."""
    reg = SA.SUBAGENT_PRESETS
    for name in ("explore", "plan", "code-reviewer"):
        assert getattr(reg[name], "may_delegate_to", ()) == ()


# ---------------------------------------------------------------------------
# Extend, never lose
# ---------------------------------------------------------------------------

def test_existing_presets_gain_no_restriction_and_keep_their_identity():
    for name in ("explore", "plan", "code-reviewer", "general-purpose"):
        p = SA.SUBAGENT_PRESETS[name]
        assert p.write_globs == ()
        assert p.output_schema is None
        assert p.department == ""
        assert p.enabled_setting == ""


def test_builtins_kebab_case_names_are_the_keys():
    for name in SA._BUILTIN_PRESETS:
        assert name == name.strip().lower()
        assert " " not in name and "_" not in name


def test_the_template_is_not_a_delegateable_type():
    """_template_subagent.md is a copy-and-fill skeleton, not an agent:
    an underscore-prefixed name is the convention for "not live" and it
    must never appear in the registry (it would be offered to the model
    and refuse validate_preset's kebab-case rule)."""
    assert "_template" not in SA.SUBAGENT_PRESETS
    assert "test_the_template" not in SA.SUBAGENT_PRESETS
