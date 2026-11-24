"""A program could not build what it was meant to orchestrate.

Code Mode exists so the model writes ONE program instead of one turn
per tool call. Its allow-list held seven tools -- search_calcs,
get_calc_info, search_docs, read_section, bash_status, bash_output,
list_docs -- and NONE of the pipeline server's 26, of which 21 are reads
by the executor's own classification (_MCP_READONLY_TOOL_BASES).

So a program could not DISCOVER a capability, RESOLVE a spec or LINT
one, which is most of what a multi-step calculation workflow is. The
case Code Mode pays off on -- the field's rule of thumb is multi-step
work with large intermediate data -- was the case it could not reach.
This is phase 1 as it was meant, not a widening of it: every tool added
here is read-only.

TWO SPELLINGS, ONE LIST. The executor dispatches
``mcp__delfin-tools__resolve_spec``; a Python function cannot be called
that, because the server name carries a hyphen. So the allow-list holds
the dispatch name and ``_py_name`` derives the short one for the stub.
A second list of short names would be the drift this repository keeps
producing.

Measured before adding them: no suffix here collides with a native
tool, with a delfin-ops tool or with the seven that were already
allowed. A test below asserts that, so a later addition that WOULD
collide fails there instead of shadowing something silently.
"""
from __future__ import annotations

import pytest

from delfin.agent import code_mode as cm

PIPE = "mcp__delfin-tools__"


# --- what a program can now reach -------------------------------------

@pytest.mark.parametrize("tool", [
    "get_guide", "catalog", "describe_capability", "describe_key",
    "compatible_successors", "resolve_spec", "validate_spec",
    "scientific_lint", "run_status", "run_diagnostics",
])
def test_a_read_only_pipeline_tool_is_allowed(tool):
    assert cm._check_allowed(PIPE + tool) is None


def test_the_seven_that_were_allowed_still_are():
    for tool in ("search_calcs", "get_calc_info", "search_docs",
                 "read_section", "bash_status", "bash_output",
                 "list_docs"):
        assert cm._check_allowed(tool) is None


def test_the_discover_resolve_lint_loop_is_complete():
    """The three steps a program needs before it can run anything."""
    for tool in ("get_guide", "resolve_spec", "validate_spec",
                 "scientific_lint"):
        assert PIPE + tool in cm.READ_ONLY_TOOLS


# --- and what it still cannot ----------------------------------------

@pytest.mark.parametrize("tool", [
    "register_module",                         # writes source and runs it
    "save_application",                        # creates an application
    "cancel_run",                              # ends somebody else's run
])
def test_a_side_effecting_pipeline_tool_is_still_refused(tool):
    """The three that are not calculations.

    run_application and submit_application left this list in part B,
    for one reason: the gate judges them. These three are not judged by
    judge_application_run and are not calculations -- the first writes
    source and runs it, the second creates an application, the third
    ends somebody else's run -- so they keep the direct route.
    """
    refusal = cm._check_allowed(PIPE + tool)
    assert refusal is not None
    assert tool in refusal, "the refusal must name the tool"


@pytest.mark.parametrize("tool", ["write_file", "bash", "apply_patch",
                                  "git_commit"])
def test_a_native_write_tool_is_still_refused(tool):
    assert cm._check_allowed(tool) is not None


def test_every_allowed_pipeline_tool_is_read_only_by_the_executor():
    """Not by this module's opinion: by the classification the gate
    itself uses. A tool that is side-effecting there and allowed here
    would be a way round the gate."""
    from delfin.agent.api_client import _MCP_READONLY_TOOL_BASES as RO
    for tool in cm.READ_ONLY_PIPELINE_TOOLS:
        base = tool.rsplit("__", 1)[-1]
        assert base in RO, f"{base} is not read-only to the executor"


# --- the two spellings ------------------------------------------------

def test_the_stub_emits_a_python_safe_name_and_dispatches_the_full_one():
    stub = cm.generate_stub()
    assert "def resolve_spec(**arguments):" in stub
    assert "_call('mcp__delfin-tools__resolve_spec'" in stub or \
           '_call("mcp__delfin-tools__resolve_spec"' in stub


def test_py_name_shortens_only_a_namespaced_name():
    assert cm._py_name(PIPE + "resolve_spec") == "resolve_spec"
    assert cm._py_name("search_calcs") == "search_calcs"


def test_no_two_allowed_tools_share_a_python_name():
    """The short name is what a program writes. Two tools sharing one
    would shadow each other, and which won would depend on sort order."""
    seen: dict[str, str] = {}
    for tool in cm.READ_ONLY_TOOLS:
        short = cm._py_name(tool)
        assert short not in seen, (
            f"{tool} and {seen.get(short)} both become {short}()")
        seen[short] = tool


def test_every_emitted_name_is_an_identifier():
    for tool in cm.READ_ONLY_TOOLS:
        short = cm._py_name(tool)
        assert short.isidentifier(), f"{tool} -> {short!r}"
