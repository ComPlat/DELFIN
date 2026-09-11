"""Tool-schema budget + context-scoped advertising.

Every request re-sends the whole tool surface, so the schema block is a
per-request cost. These tests pin three things:

  * the schema stays inside its token budget (measured with the house
    ``chars // 4`` estimate, via ``tool_schema_token_report``);
  * advertising is a strict SUBSET of what the execution layer permits —
    advertising less must never grant more;
  * the compaction that bought the budget did not drop the phrases a model
    needs in order to call a tool correctly and safely.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import (
    _DASHBOARD_AGENT_ALLOWED_TOOLS,
    _DELFIN_ONLY_TOOL_NAMES,
    _DOC_TOOLS_OPENAI,
    _ROLE_EXEC_ALLOWLIST,
    _WEAK_MODEL_CORE_TOOLS,
    ToolSurfaceContext,
    advertisable_tools,
    estimate_schema_tokens,
    role_tool_surface_report,
    tool_schema_token_report,
    tool_unavailable_reason,
)


# Measured on the pre-compaction catalogue (60 tools, same estimator).
_BASELINE_TOKENS = 11_422
# Budget with headroom for a few future tools. The compaction target was a
# 35% cut of the baseline, which lands around 7.35k; the document tools
# were added on top of that and cost ~1.2k tokens for seven.
# Raised from 8,600 (measured surface 8,589) for the three PDF-assembly
# tools: merge_pdfs, split_pdf and create_pdf cost 298 tokens together,
# and about two thirds of that is the parameter contract (names, types,
# required), which cannot be written any shorter. Their prose was kept to
# when-to-call and the argument shapes; the caveats — page counts,
# refusals, what did not verify — are returned by the runtime, where they
# cost nothing until they apply. The remaining headroom is deliberate.
# Raised once, by exactly what draft_email measures at (125 tokens), when
# that capability was added. Not a relaxation: the catalogue had no slack
# left, and paying for a new tool by shortening descriptions that already
# work would trade clarity for capacity. The tool is also named in
# _POST_COMPACTION_TOOLS below, so the diet ratchet still measures the
# compacted surface on its own and cannot be quietly undone by additions.
#
# Raised a second time, 9_125 -> 9_133, for column paging: `start_col` on
# read_document. The reader used to say "showing 40 of 87 columns" and
# name nothing the caller could do — the slice always began at column 1,
# so columns 41 to 87 were unreachable through the tool and absent from
# the column profile without that being said. A limit announced without a
# remedy is worse than no limit; it tells the model something is missing
# and leaves it to answer from the part it has.
#
# Paid for as far as it can be: the new parameter costs 20 tokens, and
# read_document's own description plus two of its parameter texts were
# tightened to return 12. The ceiling moves by the measured remainder,
# 8, and not by a round number.
#
# Raised a third time, 9_133 -> 9_287, for sum_column: 154 tokens
# measured, after its description and parameter texts were cut to the
# shortest form that still says what the tool refuses and how to get past
# the refusal. Not a relaxation. Until now nothing in the catalogue
# TOTALLED anything: a sum was arithmetic the model did in its head over a
# rendered grid, so no tool result ever held the figure the answer stated,
# and the coverage ledger that checks an answer's figures had nothing to
# check a total against. Its own text pays for what it can — the six
# parameters are the file, the column, the sheet, the group, the header
# row and the convention, and none of them can be dropped without making
# the tool unusable on the layouts the users actually have.
#
# Raised a fourth time, 9_287 -> 9_368, for the row filter on sum_column:
# period, date_column and date_convention. 87 tokens measured, 6 of them
# returned by tightening the tool's own description and the period text,
# so the ceiling moves by 81. Two thirds of that is structure -- three
# names and three types -- which no wording can shrink. The prose that is
# left is the part a model cannot work without: the three ISO shapes it
# may pass, that a period without a date column is refused, and the word
# that settles a date ambiguity. Neither of the two existing escape
# hatches could carry it: convention settles the DECIMAL reading, and
# folding dates into it would change how money parses whenever a caller
# answered a date question. Before this, a total over one month was the
# model filtering rows in its head and adding them there -- the arithmetic
# the tool exists to take away from it, done on exactly the task where a
# quietly dropped row is hardest to see.
#
# Raised a fifth time, 9_368 -> 9_375, for notebook_read's `output`.
# Seven tokens, and it is the whole of what the field costs to advertise:
# the field itself is free (it rides in the RESULT, not the schema), the
# `max_output_chars` knob was measured at 24 and dropped from the
# catalogue rather than paid for -- the default is right for every
# recorded use, a capped output says so in its own text, and a knob is
# paid for by every request whether or not anyone turns it -- and the
# description was cut from 64 tokens to 41, below where it started.
#
# What the seven buy: notebook_read summarised its outputs, so a failed
# cell came back as `error(ValueError)` with no message, no line and no
# stack, and a cell that computed a number came back as the news that it
# had computed something. Its own description tells the agent to use it
# INSTEAD of read_file, which would at least have dumped the traceback as
# raw JSON. For a scientific agent the outputs of a notebook are the
# result, and "why did this cell fail" was unanswerable.
#
# Raised a sixth time, 9_375 -> 9_393, for enter_worktree's `base_ref`.
# 40 measured; 7 returned by cutting the tool's own description, and the
# parameter text cut from 26 tokens to 11 because the rule that says WHEN
# to reach for it belongs in the integrity addendum, not in a schema every
# request pays for. 13 of the 18 left are structure -- one name, one type
# -- which no wording shrinks.
#
# What it buys is the control. "My change broke this test" was a
# hypothesis the agent had no way to test: enter_worktree only ever
# branched from HEAD, so there was no way to run the same check against
# the state before the change. Getting that wrong costs the same in both
# directions -- assume the failure is yours and you revert work that was
# right, assume it is not and you ship the regression -- and it is the
# one step of the scientific loop the addendum described without giving
# the agent a way to perform it.
#
# Raised a seventh time, 9_393 -> 9_401, for list_files. Eight tokens,
# seven of them the description and one the new `path` property; the
# first draft cost twelve and was cut.
#
# Two things it corrects, both discovered by the required-argument sweep
# and by a test the sweep then broke. `pattern` was marked REQUIRED while
# the executor has always defaulted it to "*", so the schema described a
# contract the code did not have -- and a required argument that silently
# means "everything" is how a listing of the whole workspace became the
# answer to a narrower question. And `path` was passed by callers (this
# repository's own tool-surface test passes it) and silently ignored, so
# `list_files(path="src")` answered with every file in the workspace.
# Now it is optional, documented, and it works.
_TOKEN_BUDGET = 9_401
# Capability added after the compaction was measured. The diet ratchet
# below applies to the surface the diet was measured on — new tools have
# to justify their own cost (the per-tool cap and the budget above), but
# they must not be able to make a REGRESSION in the compacted surface look
# like growth that was paid for.
_POST_COMPACTION_TOOLS = frozenset({
    "read_document", "edit_sheet", "fill_pdf_form",
    "fill_docx_template", "create_docx", "compare_tables", "sum_column",
    "fill_series", "merge_pdfs", "split_pdf", "create_pdf", "draft_email",
})


def _catalogue_names() -> set[str]:
    return {t["function"]["name"] for t in _DOC_TOOLS_OPENAI}


def _schema(name: str) -> dict:
    return next(t["function"] for t in _DOC_TOOLS_OPENAI
                if t["function"]["name"] == name)


# ---------------------------------------------------------------------------
# 1. Measurement helper
# ---------------------------------------------------------------------------


def test_token_report_totals_match_per_tool_sum():
    report = tool_schema_token_report()
    assert report["count"] == len(_DOC_TOOLS_OPENAI)
    assert set(report["tools"]) == _catalogue_names()
    assert report["total_tokens"] == sum(
        e["total"] for e in report["tools"].values())
    # Each entry splits the cost without losing more than rounding.
    for name, entry in report["tools"].items():
        parts = (entry["description"] + entry["parameter_descriptions"]
                 + entry["structure"])
        assert abs(parts - entry["total"]) <= 3, name


def test_token_report_accepts_an_explicit_catalogue():
    subset = [t for t in _DOC_TOOLS_OPENAI
              if t["function"]["name"] in {"bash", "read_file"}]
    report = tool_schema_token_report(subset)
    assert report["count"] == 2
    assert report["total_tokens"] == sum(estimate_schema_tokens(t)
                                         for t in subset)


# ---------------------------------------------------------------------------
# 2. Schema budget
# ---------------------------------------------------------------------------


def test_tool_schema_stays_within_token_budget():
    total = tool_schema_token_report()["total_tokens"]
    assert total <= _TOKEN_BUDGET, (
        f"tool schemas grew to {total} tokens per request "
        f"(budget {_TOKEN_BUDGET})")


def test_tool_schema_is_at_least_35_percent_below_baseline():
    report = tool_schema_token_report()
    total = report["total_tokens"] - sum(
        entry["total"] for name, entry in report["tools"].items()
        if name in _POST_COMPACTION_TOOLS)
    assert total <= _BASELINE_TOKENS * 0.65, (
        f"{total} tokens (excluding tools added after the compaction) is "
        f"only {100 * (1 - total / _BASELINE_TOKENS):.1f}% below the "
        f"{_BASELINE_TOKENS}-token baseline")


def test_post_compaction_tools_are_named_and_present():
    """The exclusion list above must not outlive the tools it names."""
    assert _POST_COMPACTION_TOOLS <= _catalogue_names()


def test_no_single_tool_schema_is_oversized():
    """A single tool must not eat a disproportionate slice of the surface."""
    for name, entry in tool_schema_token_report()["tools"].items():
        assert entry["total"] <= 450, f"{name} schema is {entry['total']} tokens"


# ---------------------------------------------------------------------------
# 3. (a) No advertised tool is refused by that role's allow-list
# ---------------------------------------------------------------------------


def _all_roles() -> list[str]:
    return ["", *sorted(_ROLE_EXEC_ALLOWLIST)]


@pytest.mark.parametrize("role", _all_roles())
def test_advertised_tools_are_never_refused_by_the_role(role):
    """Advertising must be a subset of execution — for every role."""
    for tool in advertisable_tools(_DOC_TOOLS_OPENAI,
                                   ToolSurfaceContext(role=role)):
        name = tool["function"]["name"]
        assert not A._tool_denied_for_role(role, name), (
            f"role {role!r} is advertised {name} but may not execute it")


@pytest.mark.parametrize("role", _all_roles())
def test_advertised_surface_is_subset_of_execution_allowlist(role):
    allow = _ROLE_EXEC_ALLOWLIST.get(role)
    if allow is None:
        return
    advertised = {t["function"]["name"] for t in
                  advertisable_tools(_DOC_TOOLS_OPENAI,
                                     ToolSurfaceContext(role=role))}
    assert advertised <= set(allow)


def test_dashboard_agent_surface_is_exactly_its_allowlist():
    """The restricted role's surface is derived from the allow-list, not from
    a second hand-maintained name list that could drift away from it."""
    advertised = {t["function"]["name"] for t in
                  advertisable_tools(
                      _DOC_TOOLS_OPENAI,
                      ToolSurfaceContext(role="dashboard_agent"))}
    assert advertised == _DASHBOARD_AGENT_ALLOWED_TOOLS & _catalogue_names()
    # And it really is a read-only surface.
    for mutating in ("bash", "write_file", "edit_file", "multi_edit",
                     "apply_patch", "undo_changes", "notebook_edit"):
        assert mutating not in advertised


def test_role_surface_report_quantifies_the_saving():
    report = role_tool_surface_report()
    assert "" in report and "dashboard_agent" in report
    full = report[""]["total_tokens"]
    restricted = report["dashboard_agent"]["total_tokens"]
    assert restricted < full
    assert report["dashboard_agent"]["count"] < report[""]["count"]


# ---------------------------------------------------------------------------
# 4. (b) Every dispatchable tool is still advertised somewhere
# ---------------------------------------------------------------------------


def test_every_catalogue_tool_is_advertised_in_the_default_context():
    """The catalogue is what the executor treats as 'known' (it is the source
    of the unknown-tool near-miss hint), so nothing in it may become
    unreachable through the advertising filters."""
    advertised = {t["function"]["name"] for t in
                  advertisable_tools(_DOC_TOOLS_OPENAI, ToolSurfaceContext())}
    assert advertised == _catalogue_names()


def test_every_context_scoped_drop_is_advertised_in_another_context():
    contexts = [
        ToolSurfaceContext(),
        ToolSurfaceContext(subagent_depth=99),
        ToolSurfaceContext(has_doc_index=False, has_calc_index=False),
        *[ToolSurfaceContext(role=r) for r in _ROLE_EXEC_ALLOWLIST],
    ]
    seen: set[str] = set()
    for ctx in contexts:
        seen |= {t["function"]["name"]
                 for t in advertisable_tools(_DOC_TOOLS_OPENAI, ctx)}
    assert seen == _catalogue_names()


def test_name_sets_only_reference_tools_that_exist():
    """A typo in a filter name set would silently hide (or leak) a tool."""
    names = _catalogue_names()
    assert _DELFIN_ONLY_TOOL_NAMES <= names
    assert _WEAK_MODEL_CORE_TOOLS <= names
    # The role allow-list is compared by BASE name, so it may name a
    # read-only MCP tool (explain_delfin_feature) beside the built-ins.
    assert _DASHBOARD_AGENT_ALLOWED_TOOLS <= names | set(A._MCP_READONLY_TOOL_BASES)
    assert A._DOC_INDEX_TOOL_NAMES <= names
    assert A._CALC_INDEX_TOOL_NAMES <= names
    assert A._SUBAGENT_SPAWN_TOOL_NAMES <= names


# ---------------------------------------------------------------------------
# 5. Context scoping mirrors real execution refusals
# ---------------------------------------------------------------------------


def test_subagent_spawn_tools_are_dropped_at_the_nesting_cap():
    """``_execute_subagent`` / ``_execute_orchestrate`` hard-refuse at the cap,
    so advertising them to a nested sub-agent is pure waste."""
    ctx = ToolSurfaceContext(subagent_depth=A._max_subagent_depth())
    advertised = {t["function"]["name"]
                  for t in advertisable_tools(_DOC_TOOLS_OPENAI, ctx)}
    assert "subagent" not in advertised
    assert "orchestrate" not in advertised
    # Collecting a background run the PARENT started stays legal.
    assert "subagent_result" in advertised
    assert _catalogue_names() - advertised == {"subagent", "orchestrate"}


def test_top_level_agent_keeps_the_spawn_tools():
    advertised = {t["function"]["name"] for t in
                  advertisable_tools(_DOC_TOOLS_OPENAI,
                                     ToolSurfaceContext(subagent_depth=0))}
    assert {"subagent", "orchestrate"} <= advertised


def test_index_backed_tools_are_dropped_without_their_index():
    """Those executors return 'Doc index not available' / 'Calc index could
    not be built' for every call, so they are refusals, not tools."""
    ctx = ToolSurfaceContext(has_doc_index=False, has_calc_index=False)
    advertised = {t["function"]["name"]
                  for t in advertisable_tools(_DOC_TOOLS_OPENAI, ctx)}
    assert not (A._DOC_INDEX_TOOL_NAMES & advertised)
    assert not (A._CALC_INDEX_TOOL_NAMES & advertised)
    assert _catalogue_names() - advertised == (
        A._DOC_INDEX_TOOL_NAMES | A._CALC_INDEX_TOOL_NAMES)


def test_unavailable_reason_is_none_in_the_default_context():
    for name in _catalogue_names():
        assert tool_unavailable_reason(name, ToolSurfaceContext()) is None


def test_unavailable_reason_handles_namespaced_mcp_names():
    ctx = ToolSurfaceContext(role="dashboard_agent")
    # Base name on the allow-list -> allowed even when namespaced.
    assert tool_unavailable_reason("mcp__kit-coding__task_create", ctx) is None
    # Base name off the allow-list -> refused, matching _gate_mcp_tool.
    assert tool_unavailable_reason("mcp__kit-coding__bash", ctx) is not None
    assert tool_unavailable_reason("mcp__delfin-docs__read_file", ctx) is not None


# ---------------------------------------------------------------------------
# 6. (c) Safety-critical phrases survived the compaction
# ---------------------------------------------------------------------------

# (tool, dotted path inside the schema, required substring)
_SAFETY_PHRASES: list[tuple[str, str, str]] = [
    # Path handling: absolute-path rule + hard secret denial.
    ("read_file", "description", "ABSOLUTE"),
    ("read_file", "description", "Secret-deny"),
    ("write_file", "description", "ABSOLUTE"),
    # Read-before-write / read-before-edit ordering.
    ("write_file", "description", "read_file first"),
    ("edit_file", "description", "read with read_file first"),
    ("edit_file", "description", "EXACTLY once"),
    ("multi_edit", "description", "read with read_file first"),
    ("multi_edit", "description", "NOTHING is written"),
    ("notebook_edit", "description", "notebook_read FIRST"),
    # Shell safety: deny-list, audit trail, and the cd-defeats-the-gate rule.
    ("bash", "description", "rm -rf"),
    ("bash", "description", "rejected"),
    ("bash", "parameters.properties.description.description", "audit"),
    ("bash", "parameters.properties.cwd.description", "cd /path &&"),
    ("bash", "parameters.properties.cwd.description", "auto-allow gate"),
    ("bash_background", "description", "SAME safety gate"),
    ("bash_status", "description", "tight loop"),
    # Permission persistence must never happen unannounced.
    ("remember_permission", "description", "ALWAYS confirm with the user"),
    ("remember_permission_bundle", "description", "ALWAYS state in chat"),
    ("remember_permission_bundle", "description", "SINGLE confirm dialog"),
    # Plan-mode contract.
    ("exit_plan_mode", "description", "'plan' mode"),
    ("exit_plan_mode", "description", "blocked"),
    ("ask_user_question", "description", "exit_plan_mode"),
    # Delegation: the permission escalation and the no-shared-context rule.
    ("subagent", "description", "FULL"),
    ("subagent", "description", "read-only"),
    ("subagent", "parameters.properties.prompt.description", "NO"),
    # Destructive / irreversible operations.
    ("apply_patch", "description", "NO file is mutated"),
    ("undo_changes", "description", "hash"),
    ("undo_changes", "description", "conflict"),
    # Network egress limits.
    ("web_fetch", "description", "RFC1918"),
    ("remote_trigger", "description", "NOT chosen by the agent"),
    # Honesty / grounding rules the harness depends on.
    ("report_verdict", "parameters.properties.criteria.description",
     "never guess PASS"),
    ("report_verdict", "description", "ONCE"),
    ("remember", "description", "secrets"),
    ("check_environment", "description", "never values"),
    ("history_search", "description", "BEFORE"),
    ("list_changes_made", "description", "never from memory"),
    # Identifier contract that silently corrupts task state when ignored.
    ("task_create", "description", "`id`"),
    ("task_create", "description", "`seq`"),
    ("task_adopt", "description", "BEFORE"),
]


def _resolve(schema: dict, dotted: str) -> str:
    node = schema
    for part in dotted.split("."):
        node = node[part]
    assert isinstance(node, str)
    return node


@pytest.mark.parametrize("tool,path,phrase", _SAFETY_PHRASES,
                         ids=[f"{t}:{p}" for t, _, p in _SAFETY_PHRASES])
def test_safety_critical_phrase_survives_compaction(tool, path, phrase):
    assert phrase in _resolve(_schema(tool), path), (
        f"{tool}.{path} lost the safety-critical phrase {phrase!r}")


# ---------------------------------------------------------------------------
# 7. The contract itself (names / params / types) is untouched
# ---------------------------------------------------------------------------


def test_every_required_param_is_declared_in_properties():
    for tool in _DOC_TOOLS_OPENAI:
        fn = tool["function"]
        params = fn.get("parameters", {})
        props = params.get("properties", {})
        for req in params.get("required", []):
            assert req in props, f"{fn['name']}: required {req!r} not declared"


def test_every_tool_has_a_non_empty_description():
    for tool in _DOC_TOOLS_OPENAI:
        fn = tool["function"]
        assert tool["type"] == "function"
        assert fn["name"]
        assert fn.get("description", "").strip(), fn["name"]


@pytest.mark.parametrize("name,required", [
    ("read_file", ["path"]),
    ("write_file", ["path", "content"]),
    ("edit_file", ["path", "old_string", "new_string"]),
    ("multi_edit", ["path", "edits"]),
    ("bash", ["command", "description"]),
    ("subagent", ["subagent_type", "description", "prompt"]),
    ("report_verdict", ["status"]),
    ("ask_user_question", ["question", "options"]),
    ("remember_permission", ["kind", "value", "rationale"]),
    ("apply_patch", ["diff"]),
    ("undo_changes", ["scope"]),
])
def test_required_parameter_lists_are_unchanged(name, required):
    assert _schema(name)["parameters"]["required"] == required


def test_subagent_type_enum_is_still_resolved_dynamically():
    """User-extensible presets (pack/agents/*_subagent.md) must keep showing
    up in the enum — a hard-coded list would silently drop them."""
    from delfin.agent.subagents import subagent_type_names
    enum = _schema("subagent")["parameters"]["properties"]["subagent_type"]["enum"]
    assert set(enum) == set(subagent_type_names())
    assert {"explore", "plan", "code-reviewer", "general-purpose"} <= set(enum)


def test_ask_user_question_option_preview_is_still_declared():
    item = _schema("ask_user_question")["parameters"]["properties"]["options"]["items"]
    assert item["properties"]["preview"]["type"] == "string"
    assert "markdown" in item["properties"]["preview"]["description"].lower()


# ---------------------------------------------------------------------------
# The MCP surface rides on every request too
# ---------------------------------------------------------------------------

def test_mcp_schemas_have_a_budget():
    """The built-in catalogue is capped and sits twelve tokens under its
    ceiling. MCP schemas were appended after that check, uncapped: two
    servers with thirty tools each would silently double the largest single
    part of a request, with nothing to notice it."""
    from delfin.agent.api_client import _mcp_schema_budget_chars

    assert _mcp_schema_budget_chars() > 0


def test_the_mcp_budget_cannot_be_switched_off(monkeypatch):
    """A budget a typo disables is not a budget."""
    from delfin.agent import api_client as A

    for junk in (0, -1, "", None, "lots"):
        monkeypatch.setattr(
            "delfin.user_settings.load_settings",
            lambda *a, **kw: {"agent": {"mcp_schema_budget_chars": junk}})
        assert A._mcp_schema_budget_chars() >= 2000


def test_dropping_a_tool_is_recorded_not_silent():
    """A surface that shrinks in silence looks like a broken server to
    whoever debugs it next.

    Anchored on the assignment rather than on a call with an empty
    argument list: the budget takes the model's context window now, and a
    test that pins the spelling of the call fails on a change that alters
    nothing about what it guards.
    """
    import inspect

    from delfin.agent import api_client as A

    source = inspect.getsource(A.OpenAIClient.stream_message)
    block = source[source.index("_mcp_budget = _mcp_schema_budget_chars("):]
    block = block[:2000]
    assert "_mcp_dropped" in block
    assert "_record_security_event" in block
    assert "were not advertised" in block


def test_the_builtin_catalogue_is_still_measured_on_its_own():
    """The MCP budget is additive; it must not become an excuse to let the
    catalogue itself grow."""
    from delfin.agent.api_client import tool_schema_token_report

    assert tool_schema_token_report()["total_tokens"] <= _TOKEN_BUDGET
