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
_TOKEN_BUDGET = 9_133
# Capability added after the compaction was measured. The diet ratchet
# below applies to the surface the diet was measured on — new tools have
# to justify their own cost (the per-tool cap and the budget above), but
# they must not be able to make a REGRESSION in the compacted surface look
# like growth that was paid for.
_POST_COMPACTION_TOOLS = frozenset({
    "read_document", "edit_sheet", "fill_pdf_form",
    "fill_docx_template", "create_docx", "compare_tables", "fill_series",
    "merge_pdfs", "split_pdf", "create_pdf", "draft_email",
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
    assert _DASHBOARD_AGENT_ALLOWED_TOOLS <= names
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
    whoever debugs it next."""
    import pathlib
    from delfin.agent import api_client as A

    source = pathlib.Path(A.__file__).read_text(encoding="utf-8")
    block = source[source.index("_mcp_budget = _mcp_schema_budget_chars()"):]
    block = block[:2000]
    assert "_mcp_dropped" in block
    assert "_record_security_event" in block
    assert "were not advertised" in block


def test_the_builtin_catalogue_is_still_measured_on_its_own():
    """The MCP budget is additive; it must not become an excuse to let the
    catalogue itself grow."""
    from delfin.agent.api_client import tool_schema_token_report

    assert tool_schema_token_report()["total_tokens"] <= _TOKEN_BUDGET
