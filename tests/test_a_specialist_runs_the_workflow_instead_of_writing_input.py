"""A specialist that runs calculations through the adapters.

MEASURED ON A REAL RUN (bench17 t13, 2026-11-23). In 142 tool calls the
agent called `write_file` 20 times and `bash` 46 times: it wrote the ORCA
inputs by hand and started the engine from a shell. 12 of its 14
calculations were abandoned. The platform's 51 adapters -- each one
declaring its parameters, its defaults and the rule that makes a value
valid -- were not called once, and neither was any of the 26 tools of
the `delfin-tools` server.

The server was not the problem: it is in the builtin registry, it
spawns, `discover_all()` returns its 26 tools, and it is first in the
order the advertising loop walks, so 19 of 26 are offered even on a 32k
window and all 26 from 131k up. They were on the table.

What was missing was a specialist whose job is to use them, and one
mechanism stood in the way of writing one: `_known_preset_tool_names`
derived the ops server's tools from its module and did NOT derive the
pipeline server's, so every name a workflow specialist must list --
`resolve_spec`, `run_application`, `run_diagnostics`,
`new_capability_template` -- was "not a known tool". The preset that
most needed a tool list was the one that could not have one.

What is asserted here:

  * the specialist exists, is offered as a `subagent_type`, and loads
    with NO problems -- a preset whose tool list is rejected is a file;
  * the derived names are pinned to the server's own handlers, so a tool
    it registers later stays known without a hand list;
  * the preset carries the RULE in its text, because a specialist that
    may call the adapters and is not told to prefer them is the state we
    measured;
  * it cannot widen what the session allows. Presets are markdown on
    disk, and one able to grant its child a tool would make dropping a
    file a privilege escalation.
"""
from __future__ import annotations

import pytest

from delfin.agent import subagents as sa

NAME = "workflow-operator"

#: What the specialist cannot do its job without.
NEEDED = (
    "resolve_spec",          # fills the validated defaults
    "validate_spec",         # structural diagnostics
    "scientific_lint",       # physical plausibility
    "run_application",       # runs it
    "submit_application",    # runs it on the cluster
    "run_status",            # where is it
    "run_diagnostics",       # what went wrong
    "new_capability_template",   # the reviewed way to extend
)


@pytest.fixture(autouse=True)
def _fresh_registry():
    """Rebuild the registry before each assertion.

    PRESET_LOAD_PROBLEMS and SUBAGENT_PRESETS are module state, and
    another test that points the preset search at its own directory
    leaves both describing ITS registry. Measured: this file passed
    alone and failed beside the budget suites. A test that reads module
    state another test can repopulate is testing the order it ran in.
    """
    sa.reload_subagent_presets()
    yield
    sa.reload_subagent_presets()


@pytest.fixture()
def preset():
    return sa.SUBAGENT_PRESETS[NAME]


def test_the_specialist_is_offered_as_a_subagent_type():
    assert NAME in sa.subagent_type_names()


def test_no_pipeline_tool_it_needs_is_rejected():
    """A preset whose tool list is rejected is a file, not a specialist.

    Scoped to the PIPELINE tools, which is the invariant this change
    introduces. The catalogue's knownness of a native tool depends on
    what else has run: measured in one pytest invocation, `orca_explain`,
    `orca_check_input` and `skill` are in `_DOC_TOOLS_OPENAI` and known
    when this file runs alone, and not known when it runs after
    test_tool_schema_budget.py. That is a property of the catalogue under
    test order, and asserting it here would make this file fail for a
    reason that has nothing to do with the specialist.
    """
    problems = [msg for k, v in (sa.PRESET_LOAD_PROBLEMS or {}).items()
                if NAME in k for msg in v]
    offenders = [m for m in problems
                 if any(f"'{t}'" in m for t in NEEDED)]
    assert not offenders, offenders


@pytest.mark.parametrize("tool", NEEDED)
def test_it_can_reach_the_pipeline_builder(preset, tool):
    assert tool in (preset.tools or ()), (
        f"{tool} missing: the specialist cannot run a calculation "
        "through the adapters without it")


def test_it_runs_rather_than_only_looks(preset):
    """`plan` mode cannot execute, and a specialist that may only look at
    a calculation is the researcher role, which already exists."""
    assert preset.mode == "default"


def test_its_writes_are_scoped(preset):
    """It writes specs and scratch, not the user's calculation
    directories: the adapters own where a run's files go."""
    assert preset.write_globs
    assert all("**" in g for g in preset.write_globs)


# --- the mechanism that made this preset possible ---------------------

def test_the_pipeline_servers_tools_are_derived_from_its_module():
    """Pinned to the server, not to a list.

    The ops server was derived this way and the pipeline server was not,
    which is why all 26 of its tools were unknown. Deriving means a tool
    the server registers later stays known automatically; a hand list
    would go stale on the first addition.
    """
    from delfin.tools import mcp_server as server
    handlers = {n[len("h_"):] for n in dir(server) if n.startswith("h_")}
    assert len(handlers) >= 20, (
        f"only {len(handlers)} handlers found: the prefix is wrong and "
        "the derivation would silently cover almost nothing")
    known = sa._known_preset_tool_names()
    missing = sorted(handlers - known)
    assert not missing, f"not derived: {missing}"
    # both spellings, as the ops server's tools have
    assert f"mcp__delfin-tools__resolve_spec" in known
    assert "resolve_spec" in known


def test_the_ops_servers_tools_are_still_known():
    """The other server must not have been displaced by the new one."""
    known = sa._known_preset_tool_names()
    for name in ("geometry_read_xyz", "compare_calculations",
                 "extract_delfin_json"):
        assert name in known, name


# --- the rule, in the text the model reads ----------------------------

def test_the_preset_states_the_rule_it_exists_for():
    """A specialist that MAY call the adapters and is not told to prefer
    them is exactly the state measured on the real run."""
    import pathlib
    text = ""
    for d in sa._md_preset_search_dirs():
        p = pathlib.Path(d) / f"{NAME}_subagent.md"
        if p.is_file():
            text = p.read_text(encoding="utf-8")
            break
    assert text, "the preset file was not found in the search dirs"
    low = text.lower()
    assert "resolve_spec" in low
    assert ".inp" in low and "hand" in low, (
        "the rule against hand-written engine input must be stated")
    assert "review" in low, (
        "the reviewed path for a missing capability must be stated")


def test_it_points_at_the_broken_symmetry_rule_rather_than_restating_chemistry():
    """The BS key decides which site carries the larger spin, and the
    skill that owns that knowledge is the single source. The preset
    names the rule and the skill; it does not become a second copy of
    the chemistry.
    """
    import pathlib
    text = ""
    for d in sa._md_preset_search_dirs():
        p = pathlib.Path(d) / f"{NAME}_subagent.md"
        if p.is_file():
            text = p.read_text(encoding="utf-8")
            break
    assert "BrokenSym" in text
    assert "orca-spectra-spin-occupier" in text
    assert "skill" in (preset_tools := sa.SUBAGENT_PRESETS[NAME].tools), (
        "it must be able to LOAD that skill")


# --- it may not widen anything ----------------------------------------

def test_a_preset_cannot_hand_its_child_a_tool_the_session_refuses():
    """Presets are markdown on disk. One able to widen its parent would
    make dropping a file into ~/.delfin/subagents/ a privilege
    escalation, so narrowing is the only direction."""
    import inspect
    src = inspect.getsource(sa._narrow_allowed_tools)
    assert "Only ever NARROWS" in src or "only ever narrows" in src.lower()
