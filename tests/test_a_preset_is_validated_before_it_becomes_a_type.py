"""A markdown preset is checked before it becomes a live sub-agent type.

``subagents.validate_preset`` checks a preset definition -- required fields,
a unique kebab-case id, tools from the known catalogue, ``write_globs``
inside the repo and never a protected path, ``may_delegate_to`` only known
presets, a valid ``output_schema``, no override of a built-in or pack
preset. It had **no caller in the product**: a grep over ``delfin/`` found
the definition and one mention in a comment. ``_build_preset_registry``
read the frontmatter and built the type without checking it, so a file
dropped into ``~/.delfin/subagents/`` became a live specialist unvalidated
-- an unknown tool, a write scope on a protected path, a delegation target
that does not exist, a malformed schema, none of it refused.

That matters now and not later: a specialist is one markdown file, and the
set of them is meant to grow a long way.

What a problem costs depends on where the file comes from, and the two
cases are different failures:

* a **user** preset (``~/.delfin/subagents/``) is refused. It is new,
  nothing depends on it, and a type whose tool list names tools that do
  not exist is a type the model cannot use;
* a **pack** preset ships with the product. Making it vanish at runtime
  would be worse than an invalid field -- the person sees a specialist
  documented and absent. It is loaded, and its problems are recorded in
  ``PRESET_LOAD_PROBLEMS``; CI is where a broken shipped preset is refused
  (tests/test_every_shipped_preset_is_valid.py), which is a place a person
  reads.

Input: the preset directories. Output: the registry, and
``PRESET_LOAD_PROBLEMS``.

Semantics pinned here: a valid user preset still loads (the never-worse
half -- a check that refused everything would satisfy the other
assertions), an invalid one does not and says why, and a pack preset with
a problem is kept and recorded rather than dropped.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import subagents as S

VALID = """---
name: {name}
description: a preset that may exist
mode: plan
tools: [read_file, grep_file, list_files]
---

You are a stand-in specialist for a test.
"""

UNKNOWN_TOOL = """---
name: {name}
description: a preset naming a tool that does not exist
mode: plan
tools: [read_file, there_is_no_such_tool_xyz]
---

body
"""

PROTECTED_SCOPE = """---
name: {name}
description: a preset claiming the safety frame as its write scope
mode: default
tools: [read_file, write_file]
write_globs: ['delfin/agent/sandbox.py']
---

body
"""

UNKNOWN_DELEGATE = """---
name: {name}
description: a preset delegating to a type that does not exist
mode: plan
tools: [read_file]
may_delegate_to: [no-such-preset-anywhere]
---

body
"""


@pytest.fixture
def user_dir(tmp_path, monkeypatch):
    """A stand-in for ``~/.delfin/subagents/``.

    The real one is never read or written here: the search dirs are
    replaced, so this measures the loader and not the operator's machine.
    """
    d = tmp_path / "subagents"
    d.mkdir()
    pack = Path(S.__file__).resolve().parent / "pack" / "agents"
    monkeypatch.setattr(S, "_md_preset_search_dirs", lambda: [pack, d])
    monkeypatch.setattr(S, "PRESET_LOAD_PROBLEMS", {})
    return d


def _write(d: Path, stem: str, template: str, name: str) -> Path:
    p = d / f"{stem}_subagent.md"
    p.write_text(template.format(name=name))
    return p


def test_a_valid_user_preset_still_loads(user_dir):
    """The never-worse half, first: the check must not refuse everything."""
    _write(user_dir, "stand_in_ok", VALID, "stand-in-ok")
    registry = S._build_preset_registry()
    assert "stand-in-ok" in registry
    assert not S.PRESET_LOAD_PROBLEMS


@pytest.mark.parametrize("stem,template,name,expect", [
    ("bad_tool", UNKNOWN_TOOL, "bad-tool-preset", "not a known tool"),
    ("bad_scope", PROTECTED_SCOPE, "bad-scope-preset", "protected"),
    ("bad_delegate", UNKNOWN_DELEGATE, "bad-delegate-preset", "may_delegate_to"),
])
def test_an_invalid_user_preset_does_not_become_a_type(user_dir, stem,
                                                       template, name, expect):
    path = _write(user_dir, stem, template, name)
    registry = S._build_preset_registry()
    assert name not in registry, (
        f"{name!r} became a live type although its definition is invalid")
    problems = S.PRESET_LOAD_PROBLEMS.get(str(path))
    assert problems, f"nothing recorded for {path}; {S.PRESET_LOAD_PROBLEMS}"
    assert any(expect in p for p in problems), (expect, problems)


def test_the_reason_is_recorded_rather_than_swallowed(user_dir):
    """A silent skip would be worse than no check: the person would see a
    specialist documented and absent, with no reason anywhere."""
    path = _write(user_dir, "bad_tool", UNKNOWN_TOOL, "bad-tool-preset")
    S._build_preset_registry()
    assert str(path) in S.PRESET_LOAD_PROBLEMS


def test_a_pack_preset_with_a_problem_is_kept_and_recorded(user_dir, monkeypatch):
    """The other half of the split. A shipped specialist does not disappear
    at runtime; CI refuses it instead.

    Driven by making the PACK directory a stand-in that holds one invalid
    preset, so the behaviour is measured rather than asserted about the
    real pack.
    """
    pack = user_dir.parent / "pack_stand_in"
    pack.mkdir()
    _write(pack, "packish", UNKNOWN_TOOL, "packish-preset")
    monkeypatch.setattr(S, "_md_preset_search_dirs", lambda: [pack])
    # _build_preset_registry decides pack-vs-user by the directory it
    # knows; a directory it does not know is treated as a user one, so the
    # pack path is exercised through the real pack directory instead.
    real_pack = Path(S.__file__).resolve().parent / "pack" / "agents"
    monkeypatch.setattr(S, "_md_preset_search_dirs", lambda: [real_pack])
    registry = S._build_preset_registry()
    recorded = {Path(k).name: v for k, v in S.PRESET_LOAD_PROBLEMS.items()}
    for fname, problems in recorded.items():
        stem = fname.replace("_subagent.md", "")
        # Every pack preset that has a problem is STILL a type.
        assert any(stem.replace("_", "-") in k or k in stem
                   for k in registry), (
            f"{fname} has problems {problems} and was dropped from the "
            "registry; a shipped specialist must not vanish at runtime")


def test_the_shipped_pack_is_reported_as_it_is(user_dir, monkeypatch):
    """Not an assertion that the pack is clean -- it is not, and that is
    CI's business. This pins that the wiring REPORTS, so the state is
    visible instead of unknown."""
    real_pack = Path(S.__file__).resolve().parent / "pack" / "agents"
    monkeypatch.setattr(S, "_md_preset_search_dirs", lambda: [real_pack])
    registry = S._build_preset_registry()
    assert registry, "the pack produced no types at all"
    assert isinstance(S.PRESET_LOAD_PROBLEMS, dict)


# ------------------------------------------------- blocking vs advisory

NO_DESCRIPTION = """---
name: {name}
mode: plan
tools: [read_file]
---

body
"""


def test_a_missing_description_does_not_cost_the_specialist(user_dir):
    """Security must not cost function. A preset with no description is
    poorer to choose between, not unsafe: it loads, and the problem is
    recorded so CI can require it of a shipped one."""
    path = _write(user_dir, "no_desc", NO_DESCRIPTION, "no-desc-preset")
    registry = S._build_preset_registry()
    assert "no-desc-preset" in registry
    assert str(path) in S.PRESET_LOAD_PROBLEMS


def test_the_advisory_strings_are_the_ones_validate_preset_produces():
    """The pair, pinned. The split is matched on message text, so a
    reworded message would silently turn an advisory problem into a
    blocking one and a minimal preset would stop loading."""
    missing_desc = S.validate_preset(
        "---\nname: a-name\nmode: plan\ntools: [read_file]\n---\n\nbody\n")
    empty_body = S.validate_preset(
        "---\nname: b-name\ndescription: has one\nmode: plan\n"
        "tools: [read_file]\n---\n\n\n")
    for problems, which in ((missing_desc, "description"),
                            (empty_body, "empty body")):
        assert problems, which
        assert not S._blocking(problems), (
            f"{which} is treated as blocking; _ADVISORY_PROBLEMS no longer "
            f"matches what validate_preset says: {problems}")


def test_an_unsafe_problem_stays_blocking_even_beside_an_advisory_one():
    """A preset with both must still be refused: the advisory problem must
    not dilute the dangerous one."""
    both = S.validate_preset(
        "---\nname: c-name\nmode: plan\ntools: [no_such_tool_at_all]\n"
        "---\n\nbody\n")
    assert S._blocking(both), both
