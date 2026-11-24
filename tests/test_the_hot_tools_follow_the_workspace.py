"""The schemas kept hot were chosen from sessions that were not chemistry.

The deferred surface sends 28 tool schemas on every request and names
the rest. The 28 were derived from 81,949 real tool calls -- and a
second review on 2026-11-23 measured WHOSE calls those were:

    sessions with tool calls              584
    of them, any chemistry tool at all     28   (4%)

So the sample was 96% coding work by swarm agents, and the hot set it
produced is the hot set of a coding agent: search_calcs, get_calc_info,
calc_summary, read_section, list_docs and list_sections were deferred,
and a chemist's first use of each cost an extra round. Nothing was
unreachable -- every one stayed executable and named -- but the
principle was wrong: a usage sample that does not look like the person
the product is for.

THE FIX REUSES A CLASSIFICATION THAT ALREADY EXISTED.
``_DELFIN_ONLY_TOOL_NAMES`` is the repository's own answer to "which
tools are chemistry", and the request path already strips them in a
workspace that is not a DELFIN one. In a DELFIN workspace they are now
hot. In any other workspace nothing changes, because they are not
offered there at all.

AND ONE ANSWER INSTEAD OF TWO. The hot set was consulted in two places
-- ``advertisable_tools(defer=True)``, which the budget and the report
measure, and an inline block in the request path, which is what ships.
Two copies of one answer drift, and this change would have made them
disagree on the first day. Both now ask ``_hot_tool_names``.
"""
from __future__ import annotations

import ast
import inspect
import json

from delfin.agent import api_client as ac

CHEMISTRY = ("search_calcs", "get_calc_info", "calc_summary",
             "read_section", "list_docs", "list_sections")


def _sent(ctx):
    return {t["function"]["name"]
            for t in ac.advertisable_tools(ac._DOC_TOOLS_OPENAI, ctx,
                                           defer=True)}


def test_in_a_chemistry_workspace_the_chemistry_schemas_are_sent():
    sent = _sent(ac.ToolSurfaceContext(role="solo_agent",
                                       chemistry_workspace=True))
    missing = [n for n in CHEMISTRY if n not in sent]
    assert not missing, f"deferred in a DELFIN workspace: {missing}"


def test_elsewhere_the_surface_is_what_it_was():
    """Outside a DELFIN workspace these tools are stripped before
    deferral, so the coding surface must not grow by a single schema."""
    sent = _sent(ac.ToolSurfaceContext(role="solo_agent",
                                       chemistry_workspace=False))
    assert not (set(CHEMISTRY) & sent), sorted(set(CHEMISTRY) & sent)


def test_the_hot_set_is_one_answer():
    """The request path and advertisable_tools ask the same function,
    and neither reads the frozen name set past it."""
    src = inspect.getsource(ac)
    request = src.split("# DEFERRAL, last of the native filters")[1][:3000]
    assert "_hot_tool_names(" in request
    assert "_hot_tool_names(" in inspect.getsource(ac.advertisable_tools)
    assert "_hot_tool_names(" in inspect.getsource(ac.deferred_tool_names)
    for fn in (ac.advertisable_tools, ac.deferred_tool_names):
        code = ast.unparse(ast.parse(inspect.getsource(fn)))
        assert "_HOT_TOOL_NAMES" not in code.split('"""')[-1], fn.__name__


def test_what_ships_in_a_chemistry_workspace_stays_within_budget():
    """The larger of the two surfaces is the one the ceiling must hold."""
    import tests.test_tool_schema_budget as B

    ctx = ac.ToolSurfaceContext(role="solo_agent", chemistry_workspace=True)
    shipped = ac.advertisable_tools(ac._DOC_TOOLS_OPENAI, ctx, defer=True)
    tokens = sum(len(json.dumps(t, ensure_ascii=False)) // 4
                 for t in shipped)
    tokens += len(ac.deferred_tool_listing("solo_agent")) // 4
    assert tokens <= B._SHIPPED_BUDGET, tokens
