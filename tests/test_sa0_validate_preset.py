"""SA0 phase 2: validate_preset — the agent-builder expert's one check.

Input: a preset definition (a file path or its text), the live preset
registry, and the executor's real registered tool names. Output: a list
of problems (empty = valid), one per defect, with the field named.

Checked (operator brief, 2026-11-01):
- required fields present (name, description, prompt body);
- unique kebab-case id, not a built-in or pack preset name (an override
  of either is refused — user files may override, the agent-builder may
  not write one);
- tools only from the known catalogue (the REAL registered surface,
  built exactly the way api_client builds its T4 grounding list);
- write_globs inside the repository and never a protected path
  (api_client, engine, sandbox, permission_rules, approval_guard,
  turn_outcome, benchmark_runner, kit_* — the safety frame);
- may_delegate_to only names known presets;
- output_schema is a valid JSON schema (a conforming sample validates);
- every existing preset passes (test at the end).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import subagents as SA


def _write(tmp_path, text, name="my-expert_subagent.md"):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


def _md(**over):
    lines = ["---"]
    fields = {
        "name": "my-expert",
        "description": "does one expert thing",
        "tools": "read_file, grep_file",
    }
    fields.update(over)
    for k, v in fields.items():
        if v is None:
            continue
        lines.append(f"{k}: {v}")
    lines.append("---")
    lines.append("You are an expert. Do the one thing.")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Valid preset
# ---------------------------------------------------------------------------

def test_a_valid_preset_has_no_problems(tmp_path):
    path = _write(tmp_path, _md())
    assert SA.validate_preset(str(path)) == []


def test_text_form_is_accepted_too(tmp_path):
    problems = SA.validate_preset(_md())
    assert problems == []


# ---------------------------------------------------------------------------
# Required fields / identity
# ---------------------------------------------------------------------------

def test_missing_name_and_description_are_named(tmp_path):
    path = _write(tmp_path, _md(name=None, description=None))
    problems = SA.validate_preset(str(path))
    assert any("name" in p for p in problems)
    assert any("description" in p for p in problems)


def test_missing_body_is_a_problem(tmp_path):
    path = _write(tmp_path, "---\nname: my-expert\ndescription: d\n---\n   \n")
    problems = SA.validate_preset(str(path))
    assert any("prompt" in p or "body" in p for p in problems)


@pytest.mark.parametrize("bad", ["My_Expert", "has space", "-lead", "x_y"])
def test_non_kebab_case_names_are_refused(tmp_path, bad):
    path = _write(tmp_path, _md(name=bad))
    problems = SA.validate_preset(str(path))
    assert any("kebab" in p for p in problems)


# ---------------------------------------------------------------------------
# No override of a built-in or pack preset
# ---------------------------------------------------------------------------

def test_override_of_a_builtin_is_refused(tmp_path):
    path = _write(tmp_path, _md(name="explore"))
    problems = SA.validate_preset(str(path))
    assert any("explore" in p and ("override" in p or "exists" in p)
               for p in problems)


def test_override_of_a_pack_preset_is_refused(tmp_path):
    pack_names = [n for n in SA.SUBAGENT_PRESETS
                  if n not in SA._BUILTIN_PRESETS]
    assert pack_names, "the pack must ship presets for this test"
    target = pack_names[0]
    path = _write(tmp_path, _md(name=target))
    problems = SA.validate_preset(str(path))
    assert any(target in p and ("override" in p or "exists" in p)
               for p in problems)


def test_an_unused_id_is_accepted(tmp_path):
    # "my-expert" is not a real preset anywhere: no override problem.
    assert not any("override" in p or "exists" in p
                   for p in SA.validate_preset(_md()))


# ---------------------------------------------------------------------------
# Tools from the known catalogue
# ---------------------------------------------------------------------------

def test_unknown_tool_is_named(tmp_path):
    path = _write(tmp_path, _md(tools="read_file, not_a_real_tool"))
    problems = SA.validate_preset(str(path))
    assert any("not_a_real_tool" in p and "tool" in p for p in problems)


def test_known_tools_are_accepted(tmp_path):
    path = _write(tmp_path, _md(tools="read_file, grep_file, search_docs"))
    assert SA.validate_preset(str(path)) == []


# ---------------------------------------------------------------------------
# write_globs: inside the repo, never a protected path
# ---------------------------------------------------------------------------

def test_write_globs_outside_the_repo_are_refused(tmp_path):
    path = _write(tmp_path, _md(write_globs="/etc/*, ../secrets/*"))
    problems = SA.validate_preset(str(path))
    assert any("/etc/*" in p for p in problems)
    assert any("../secrets/*" in p for p in problems)


def test_protected_paths_are_refused(tmp_path):
    path = _write(
        tmp_path, _md(write_globs="delfin/agent/api_client.py, "
                                    "delfin/agent/sandbox.py, "
                                    "delfin/agent/pack/skills/*"))
    problems = SA.validate_preset(str(path))
    assert len([p for p in problems if "protected" in p]) >= 3


def test_an_ordinary_repo_glob_is_accepted(tmp_path):
    path = _write(tmp_path,
                  _md(write_globs="delfin/workflows/xyz/*, results/*"))
    problems = [p for p in SA.validate_preset(str(path)) if "glob" in p]
    assert problems == []


# ---------------------------------------------------------------------------
# may_delegate_to
# ---------------------------------------------------------------------------

def test_may_delegate_to_unknown_expert_is_refused(tmp_path):
    path = _write(tmp_path,
                  _md(may_delegate_to="not-a-preset, explore"))
    problems = SA.validate_preset(str(path))
    assert any("not-a-preset" in p and "delegate" in p for p in problems)


def test_may_delegate_to_known_builtin_is_accepted(tmp_path):
    path = _write(tmp_path, _md(may_delegate_to="explore"))
    problems = [p for p in SA.validate_preset(str(path))
                if "delegate" in p]
    assert problems == []


# ---------------------------------------------------------------------------
# output_schema
# ---------------------------------------------------------------------------

def test_invalid_output_schema_json_is_refused(tmp_path):
    path = _write(tmp_path, _md(output_schema='{"type": "object", "o'))
    problems = SA.validate_preset(str(path))
    assert any("schema" in p for p in problems)


def test_valid_output_schema_is_accepted(tmp_path):
    schema = json.dumps({
        "type": "object",
        "properties": {"verdict": {"type": "string"}},
        "required": ["verdict"],
    })
    path = _write(tmp_path, _md(output_schema=schema))
    problems = [p for p in SA.validate_preset(str(path))
                if "schema" in p]
    assert problems == []


def test_non_object_schema_is_refused(tmp_path):
    path = _write(tmp_path, _md(output_schema='["verdict"]'))
    problems = SA.validate_preset(str(path))
    assert any("schema" in p for p in problems)


# ---------------------------------------------------------------------------
# The check itself never explodes
# ---------------------------------------------------------------------------

def test_a_missing_file_is_one_problem_not_a_traceback(tmp_path):
    problems = SA.validate_preset(str(tmp_path / "nope.md"))
    assert len(problems) == 1 and "read" in problems[0].lower()


# ---------------------------------------------------------------------------
# Every existing preset passes
# ---------------------------------------------------------------------------

def test_every_existing_preset_passes_validation():
    """The validator's rules hold for the presets DELFIN ships: a rule
    that rejects a shipped preset would break every user's definitions."""
    for name, preset in SA.SUBAGENT_PRESETS.items():
        problems = SA.validate_preset(_md(
            name=preset.name,
            description=preset.description,
            tools=", ".join(preset.tools) if preset.tools else None,
            write_globs=", ".join(preset.write_globs)
            if preset.write_globs else None,
            department=preset.department or None,
            may_delegate_to=", ".join(preset.may_delegate_to)
            if preset.may_delegate_to else None,
            display_name=preset.display_name or None,
            enabled_setting=preset.enabled_setting or None,
            output_schema=json.dumps(preset.output_schema)
            if preset.output_schema else None,
        ))
        # The name of a shipped preset IS an existing id, so the
        # uniqueness rule must not flag it against itself.
        problems = [p for p in problems
                    if not (name in p and ("override" in p or "exists" in p))]
        assert problems == [], f"{name}: {problems}"
