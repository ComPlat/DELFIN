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

#: Measured: the ten most-called tools cost 1,601 tokens and are 95% of
#: all calls; every name costs about 340. The target leaves room for the
#: hot set, every name, and a group or two a turn pulls in -- and it is
#: a THIRD of today's surface rather than a ratchet on it.
SURFACE_TARGET_TOKENS = 4_000


def _surface():
    return ac.role_tool_surface_report(["solo_agent"])["solo_agent"]


# --- the cost ---------------------------------------------------------

def test_the_always_on_surface_is_a_fraction_of_what_it_was():
    rep = _surface()
    assert rep["total_tokens"] <= SURFACE_TARGET_TOKENS, (
        f"{rep['total_tokens']} tokens of schema on every request, for "
        f"{rep['count']} tools; ten of them are 95% of all calls")


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
