"""11,729 tokens of tool schema ride on every request; ten tools do 95%.

MEASURED on 81,949 real tool calls (578 trace files in
~/.delfin/tool_traces, read only):

    offered surface        96 tools   11,729 tokens   100%
    never called           18 tools    1,599 tokens    13%
    fewer than 10 calls    61 tools    6,799 tokens    57%
    the ten most called    10 tools    1,601 tokens    14%  <- 95% of calls

bash alone is 46,979 calls, 57% of everything; then read_file 14,350,
edit_file 4,379, grep_file 3,581, write_file 2,372.

So 1,601 tokens carry 95% of the work and 6,799 tokens carry tools with
fewer than ten calls each -- on every request, for months. Today's
budget ratchets are 8 to 29 tokens; this is three orders of magnitude
above them, and it is the budget every other capability has been
competing for.

The 18 never-called ones are not a random set: cluster_* (7), lab_* (5),
the PDF/office group (5) and view_image. Capabilities for sessions that
did not happen yet -- the Pi's lab control, the cluster runs. They are
not dead and must not be deleted; they must arrive when they are
needed.

THE PATTERN IS ALREADY IN THE HOUSE. ops_server has
``tool_list_tools(category, query)`` and ``tool_describe_tool(name)``
("Use AFTER list_tools to read the docstring for one specific tool"),
and solo_agent.md documents the loop. It covers the ops API and not
DELFIN's own 96 native schemas. This extends it rather than inventing
anything.

RED FIRST, and deliberately: three of these fail today. The fourth must
be green now and stay green, because deferral that breaks execution is
worse than the cost it saves.
"""
from __future__ import annotations

import json
import tempfile

import pytest

from delfin.agent import api_client as ac

#: A tool measured at ZERO calls in 81,949: one of the PDF group.
A_DEFERRED_TOOL = "create_pdf"

#: Measured AFTER the change, and it corrects the figure this file was
#: written with. The prediction recorded before building was 4,000, from
#: the hot set (3,992) alone. It was wrong by 10% for one reason: the
#: MECHANISM has its own cost, which the prediction did not carry --
#: describe_native_tool's schema is 136 tokens and the name listing is
#: 280. Measured end to end:
#:
#:     without deferral   97 tools   11,865 tokens
#:     hot surface        28 tools    4,128
#:     + the listing                    280
#:     = on every request              4,408   (63% less, 7,457 freed)
#:
#: The alternative that would have kept the prediction was a threshold
#: of 100 calls: 22 tools, 3,532 + 416 = 3,948, under 4,000. It was NOT
#: taken, because it costs a point of coverage (98% against 99%) and
#: one extra round trip costs far more than the 460 tokens it saves.
#: So the target is raised by the measured amount and the reason is
#: here, which is the only way a budget may move.
#
#: No longer a constant of its own. It was a second copy of
#: tests/test_tool_schema_budget.py::_SHIPPED_BUDGET, and on 2026-11-23
#: the two disagreed the first time the hot set changed -- the chemistry
#: tools became hot in a DELFIN workspace (+438 tokens, measured) and
#: this copy failed while the real ceiling had already moved. One
#: ceiling, read from where it is ratcheted.
from tests.test_tool_schema_budget import _SHIPPED_BUDGET as SURFACE_TARGET_TOKENS  # noqa: E402

#: What the whole mechanism costs per request: the hot schemas plus the
#: always-sent names. The same ceiling, so a later hot-set addition
#: cannot be paid for by quietly dropping the listing, or the reverse.
TOTAL_TARGET_TOKENS = SURFACE_TARGET_TOKENS


def _surface():
    # What SHIPS, which this file is about; the report's default is what
    # a role can call.
    return ac.role_tool_surface_report(["solo_agent"],
                                       defer=True)["solo_agent"]


# --- the cost ---------------------------------------------------------

def test_the_always_on_surface_is_a_fraction_of_what_it_was():
    rep = _surface()
    assert rep["total_tokens"] <= SURFACE_TARGET_TOKENS, (
        f"{rep['total_tokens']} tokens of schema on every request, for "
        f"{rep['count']} tools; ten of them are 95% of all calls")


def test_the_mechanism_pays_for_itself_including_its_own_cost():
    """A saving measured without the mechanism's own cost is not a
    saving. The prediction this file opened with left out 416 tokens --
    describe_native_tool's schema and the name listing -- and was wrong
    by 10% because of it."""
    import json

    rep = _surface()
    listing = ac.deferred_tool_listing("solo_agent")
    total = rep["total_tokens"] + len(listing) // 4
    assert total <= TOTAL_TARGET_TOKENS, (
        f"{rep['total_tokens']} tokens of schema + {len(listing)//4} of "
        f"listing = {total}")

    full = ac.advertisable_tools(
        ac._DOC_TOOLS_OPENAI,
        ac.ToolSurfaceContext(role="solo_agent"), defer=False)
    undeferred = sum(
        len(json.dumps(t, ensure_ascii=False)) // 4 for t in full)
    assert undeferred > 2 * total, (
        f"deferral must more than halve the surface: {undeferred} -> "
        f"{total}; if this ever fails the hot set has grown into the "
        f"cost it was meant to remove")


def test_the_hot_set_is_derived_and_not_curated():
    """The rule is in the module beside the data: at least 50 measured
    calls, plus every always-allowed tool and `subagent`. A name added
    for any other reason has no counter beside it."""
    import inspect
    src = inspect.getsource(ac)
    block = src.split("_HOT_TOOL_NAMES")[1][:4000]
    assert "81,949" in src or "81,949" in block
    for n in ("bash", "read_file", "subagent", "ask_user_question"):
        assert n in ac._HOT_TOOL_NAMES, n
    # and the two kept despite being rare are explained
    assert "UNDER-used" in src


# --- and the capability, which must not be lost ----------------------

def test_a_deferred_tool_is_still_executable():
    """GREEN TODAY AND MUST STAY GREEN. Deferral moves a schema out of
    the request; it must not move the tool out of the product."""
    ex = ac._DocToolExecutor()
    perms = ac.KitToolPermissions(workspace=tempfile.mkdtemp(),
                                  mode="default")
    out = ex.execute(A_DEFERRED_TOOL, {}, perms)
    assert "Unknown tool" not in str(out), out


def test_a_deferred_tool_is_named_somewhere_the_model_reads():
    """A capability with no cue is a capability nobody uses: measured
    here before, on the subagent table -- six roles callable, named
    nowhere, chosen zero times in 79,056 calls."""
    listing = ac.deferred_tool_listing("solo_agent")
    assert A_DEFERRED_TOOL in listing
    assert len(listing) <= 2_000, (
        f"{len(listing)} characters of listing is a second always-on cost")


def test_the_schema_arrives_on_demand():
    ex = ac._DocToolExecutor()
    perms = ac.KitToolPermissions(workspace=tempfile.mkdtemp(),
                                  mode="default")
    out = json.loads(ex.execute("describe_native_tool",
                                {"name": A_DEFERRED_TOOL}, perms))
    assert out.get("name") == A_DEFERRED_TOOL
    params = json.dumps(out.get("parameters") or out.get("schema") or {})
    assert "properties" in params, out


# --- and it has to reach the REQUEST, not just this module ------------

def test_deferral_is_the_last_native_filter_in_the_request_path():
    """A listing built before the later filters would name tools they
    then strip -- DELFIN-only tools in a generic repo, hidden groups,
    role denials. So the first call passes defer=False and deferral
    happens below it, and the listing is computed from what survived.

    This is asserted on the request path's own source because the
    alternative -- a listing that is right in this module and wrong in
    the request -- is the defect this file exists to prevent, and it
    has appeared three times in this repository in one day.
    """
    import inspect
    src = inspect.getsource(ac)
    head = src.split("advertisable_tools(\n            advertised_tools")[1][:200]
    assert "defer=False" in head, (
        "the request path must take the FULL surface first, or the "
        "listing names tools a later filter removes")
    assert "Deferred tools:" in src, (
        "nothing in the request path puts the names where the model "
        "reads them")


def test_the_listing_rides_on_the_tool_that_acts_on_it():
    """The house idiom, from the skill listing: the cue sits on the
    tool's own description, so no prompt-head budget is spent and the
    names are where the model reads about the mechanism."""
    import inspect
    src = inspect.getsource(ac)
    block = src.split("Deferred tools:")[0][-1500:]
    assert 'describe_native_tool' in block
    assert '"description"' in block


def test_a_failure_to_defer_keeps_the_whole_surface():
    """Fail OPEN on cost, never on capability: a broken deferral must
    send too many schemas, never too few with no listing. The skill
    listing was invisible for as long as it existed because a bare
    except hid an AttributeError."""
    import inspect
    src = inspect.getsource(ac)
    warn = src.index("tool schemas not deferred")
    # the guard is the nearest except above the warning, and it is not bare
    guard = src.rindex("except Exception as exc:", 0, warn)
    assert warn - guard < 800, "the warning belongs to another handler"
    # and the warning names the consequence, so a reader of the log knows
    # the turn paid tokens rather than lost tools
    assert "sending the full surface" in src


# --- a role with its own list ----------------------------------------

def _listed_roles():
    return sorted(r for r in ac._ROLE_EXEC_ALLOWLIST) + ["workflow-operator"]


import pytest as _pytest  # noqa: E402


@_pytest.mark.parametrize("role", _listed_roles())
def test_a_role_with_an_allow_list_loses_nothing_to_deferral(role):
    """A role with an execution allow-list ships its WHOLE list.

    Measured in a second review on 2026-11-23: describe_native_tool is
    not on any role's allow-list, so for office_agent, dashboard_agent
    and workflow-operator the role filter removed it -- and with it the
    listing and the only way to fetch a deferred schema. Their deferred
    tools were invisible and unreachable; the office agent lost
    read_document, its own core tool. Every test of the mechanism had
    used solo_agent.

    The hot set was derived from a generalist's calls and never measured
    for these roles; their allow-list already IS their curated surface.
    So deferral applies to a role without a list and not to one with.
    """
    allow = ac._role_allowlist(role)
    assert allow is not None, f"{role} has no allow-list; test is void"
    ctx = ac.ToolSurfaceContext(role=role)
    shipped = {t["function"]["name"] for t in
               ac.advertisable_tools(ac._DOC_TOOLS_OPENAI, ctx, defer=True)}
    available = {t["function"]["name"] for t in
                 ac.advertisable_tools(ac._DOC_TOOLS_OPENAI, ctx)}
    assert shipped == available, (
        f"{role} loses {sorted(available - shipped)} to deferral")
    assert ac.deferred_tool_names(ac._DOC_TOOLS_OPENAI, ctx) == []
