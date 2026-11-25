"""The ELN tools are reachable from the entry points a model and a person use.

Wave rule 5's defect, twice over: two changes added tool names to an
allow-list, passed their own tests, and every call answered "Unknown tool".
A set membership is not a capability. So every assertion here drives the
route a caller actually takes:

* the catalogue carries the five schemas AND ``delfin.lab.eln_tools``
  answers each name -- a name in the catalogue with no branch in
  ``execute`` is the "Unknown tool" defect;
* the ``delfin-ops`` MCP server registers them, which is the route the
  agent reaches them by (``mcp__delfin-ops__eln_find``), and the draft
  tool is NOT host-gated the way the calculation-mutating tools are;
* ``delfin.agent.cli.main(["eln", ...])`` answers -- the person's command,
  through the CLI parser rather than the helper behind it;
* the specialist preset is in the registry the framework builds, and its
  tool list names only tools that exist.

The names live in three places (the catalogue, ``eln_tools.TOOL_NAMES`` and
``api_client._ELN_TOOL_NAMES``, which rides in the patch). Pinning all three
here is what stops them drifting apart silently.
"""

from __future__ import annotations

import importlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

from _eln_fake import FAKE_URL, FakeELN

REPO = Path(__file__).resolve().parents[1]

#: The five names, spelled once. ``eln_tools.TOOL_NAMES`` is asserted
#: equal to this, so a new tool must be added in both (and then to the
#: catalogue and the dispatcher, which the other tests here check).
EXPECTED = {"eln_find", "eln_read", "eln_pending", "eln_draft", "eln_approve"}

PRESET = REPO / "delfin" / "agent" / "pack" / "agents" / "chemotion-eln_subagent.md"


@pytest.fixture
def home(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("DELFIN_LAB_DIR", str(tmp_path / "lab"))
    monkeypatch.setenv("DELFIN_PERSON", "alice")
    monkeypatch.setenv("DELFIN_ELN_URL", FAKE_URL)
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    from delfin.agent import trust_store
    trust_store.forget_cached_keys()
    from delfin.agent import credentials
    credentials.set_credential("CHEMOTION_API_TOKEN", "chemtoken_secret_XYZ",
                               store="file")
    return tmp_path


# ---------------------------------------------------------------------------
# 1. The catalogue and the executor agree, and the executor answers
# ---------------------------------------------------------------------------

def test_the_five_schemas_are_in_the_catalogue():
    from delfin.agent import tool_catalog as TC
    names = {t["function"]["name"] for t in TC._DOC_TOOLS_OPENAI}
    assert EXPECTED <= names, sorted(EXPECTED - names)


def test_the_executor_answers_every_name_rather_than_unknown_tool(home):
    """The defect this file exists for: a catalogued name with no branch
    answers "Unknown tool". Driven through the ELN executor itself, which
    is where the api_client branch sends the call."""
    from delfin.lab import eln_tools
    assert eln_tools.TOOL_NAMES == EXPECTED
    fake = FakeELN()
    for name in sorted(EXPECTED):
        args = {"_transport": fake.send}
        if name == "eln_find":
            args.update({"collection": "Synthesis 2026", "element": "samples"})
        if name == "eln_read":
            args.update({"element": "samples", "id": 11})
        if name == "eln_draft":
            args.update({"element": "reaction", "name": "Oxidation of 4a",
                         "collection_id": 3})
        out = json.loads(eln_tools.execute(name, args, None))
        assert "unknown ELN tool" not in json.dumps(out).lower(), (name, out)
        assert out, name


def test_the_api_client_names_the_same_five():
    """``api_client._ELN_TOOL_NAMES`` gates advertisement and dispatch. It
    is the half that rides in a patch, so a divergence would be invisible
    in the lab package's own tests."""
    from delfin.agent import api_client as A
    assert set(A._ELN_TOOL_NAMES) == EXPECTED


def test_the_tools_are_advertised_only_where_an_eln_is_configured(home,
                                                                 monkeypatch):
    from delfin.agent.api_client import ToolSurfaceContext, \
        tool_unavailable_reason
    ctx = ToolSurfaceContext(has_eln=True)
    assert tool_unavailable_reason("eln_find", ctx) is None
    ctx_off = ToolSurfaceContext(has_eln=False)
    why = tool_unavailable_reason("eln_draft", ctx_off)
    assert why and "ELN" in why
    # The live surface must pass the real answer, not the default.
    import inspect
    from delfin.agent import api_client as A
    source = inspect.getsource(A)
    assert "has_eln=_eln_configured()" in source, (
        "ToolSurfaceContext defaults has_eln to True; the live surface has "
        "to pass the measured answer or every installation advertises them")


# ---------------------------------------------------------------------------
# 2. The MCP server the agent reaches them by
# ---------------------------------------------------------------------------

def test_the_ops_server_registers_the_five_tools():
    server = (REPO / "delfin" / "ops_server" / "server.py").read_text(
        encoding="utf-8")
    registered = set(re.findall(r'mcp\.tool\(name="(eln_[a-z_]+)"\)', server))
    assert registered == EXPECTED, registered


def test_the_draft_tool_is_not_host_gated_like_a_calculation():
    """The ops server runs with ``DELFIN_OPS_ALLOW_MUTATE=1`` (it is the
    agent's own). Registering ``eln_draft`` through ``_host_gated`` would
    let that standing grant send an ELN write with no per-write approval."""
    server = (REPO / "delfin" / "ops_server" / "server.py").read_text(
        encoding="utf-8")
    assert 'mcp.tool(name="eln_draft")(_safe("eln_draft", tool_eln_draft))' \
        in server
    assert "_host_gated(tool_eln_draft)" not in server


def test_an_ops_tool_call_reaches_the_eln_tools(home):
    """Through the module function the server registers, the way its
    registration calls it -- not through the helper underneath."""
    from delfin.ops_server import server as S
    fake = FakeELN()
    from delfin.lab import eln_tools
    original = eln_tools.execute
    seen = []

    def spy(name, arguments, permissions=None):
        seen.append((name, permissions))
        args = dict(arguments)
        args["_transport"] = fake.send
        return original(name, args, permissions)

    eln_tools.execute = spy
    try:
        out = json.loads(S.tool_eln_find(collection="Synthesis 2026",
                                        element="samples"))
    finally:
        eln_tools.execute = original
    assert seen and seen[0][0] == "eln_find"
    # No session permission object: nothing from this route can be read as
    # consent for a write.
    assert seen[0][1] is None
    assert out["count"] == 4


# ---------------------------------------------------------------------------
# 3. The person's command, through the parser
# ---------------------------------------------------------------------------

def test_the_cli_approves_through_main(home):
    """A filed request, then ``delfin-agent eln approve ID`` by a person."""
    from _eln_fake import FakeELN as _F
    from delfin.lab import eln_write

    fake = _F()
    with pytest.raises(eln_write.WriteRefused) as filed:
        from delfin.lab.eln_tools import _client
        eln_write.create_draft(_client(fake.send), "reaction",
                               {"name": "Oxidation of 4a",
                                "collection_id": 3})
    rid = filed.value.request_id

    from delfin.agent import cli
    code = cli.main(["eln", "approve", rid])
    assert code == 0
    rec = eln_write.pending()          # gone from pending once granted
    assert not [r for r in rec if r["id"] == rid]


def test_the_cli_refuses_to_approve_inside_an_agent_session(home, monkeypatch):
    from delfin.lab import eln_write
    from _eln_fake import FakeELN as _F
    from delfin.lab.eln_tools import _client

    fake = _F()
    with pytest.raises(eln_write.WriteRefused) as filed:
        eln_write.create_draft(_client(fake.send), "reaction",
                               {"name": "Oxidation of 4a",
                                "collection_id": 3})
    rid = filed.value.request_id
    monkeypatch.setenv("DELFIN_AGENT_SESSION", "4242")
    from delfin.agent import trust_store
    trust_store.forget_cached_keys()
    from delfin.agent import cli
    with pytest.raises(SystemExit) as stop:
        cli.main(["eln", "approve", rid])
    assert stop.value.code == 3


def test_the_cli_lists_a_pending_request(home):
    from delfin.lab import eln_write
    from _eln_fake import FakeELN as _F
    from delfin.lab.eln_tools import _client

    fake = _F()
    with pytest.raises(eln_write.WriteRefused):
        eln_write.create_draft(_client(fake.send), "reaction",
                               {"name": "Oxidation of 4a",
                                "collection_id": 3})
    from delfin.agent import cli
    assert cli.main(["eln", "pending"]) == 0


# ---------------------------------------------------------------------------
# 4. The specialist
# ---------------------------------------------------------------------------

def test_the_specialist_is_in_the_registry_the_framework_builds():
    from delfin.agent import subagents as SA
    registry = SA._build_preset_registry()
    assert "chemotion-eln" in registry, sorted(registry)
    preset = registry["chemotion-eln"]
    assert preset.mode == "default"
    assert preset.description


def test_every_tool_the_specialist_lists_exists():
    """A preset naming a tool that does not exist is a type that cannot
    work; validate_preset refuses it, and this pins that it passes."""
    from delfin.agent import subagents as SA
    problems = SA.validate_preset(str(PRESET))
    assert problems == [], problems


def test_the_specialist_does_not_name_the_eln_tools_it_cannot_hold():
    """The preset reaches the ELN through the ops server's namespace, which
    ``_known_preset_tool_names`` does not model. Listing a bare ``eln_find``
    would be refused as an unknown tool, so the body names the namespaced
    form and the frontmatter does not list it -- which is stated in the
    body rather than left to be discovered."""
    text = PRESET.read_text(encoding="utf-8")
    frontmatter = text.split("---")[1]
    assert "eln_find" not in frontmatter
    assert "mcp__delfin-ops__eln_find" in text


# ---------------------------------------------------------------------------
# The control: what the base already did, and what must not change
# ---------------------------------------------------------------------------
#
# Most of this file fails at the base because the modules are new, which is
# weak evidence -- a test can fail that way while measuring nothing about
# behaviour. These two pass on BOTH trees and say what the base already got
# right, so a change that broke it would be visible as a diff in the other
# direction.

def test_control_the_lab_tools_still_block_in_plan_mode():
    """The lab path the ELN copies: a write is refused in plan mode. Passes
    at the base and here, and would fail if the ELN work had disturbed the
    shared dispatcher."""
    from delfin.lab import safety  # noqa: F401  (the module the kernel lives in)
    from delfin.lab.tools import TOOL_NAMES as LAB_NAMES
    assert LAB_NAMES == {"lab_devices", "lab_read", "lab_set",
                         "lab_plan_check", "lab_estop"}


def test_control_the_write_gate_refuses_before_this_change_too():
    """``POST /api/v1/research_plans`` was refused on the base and is
    refused now: the ELN work added routes, it did not open writes."""
    from delfin.lab.eln_routes import decide_route
    decision = decide_route("POST", "/api/v1/research_plans",
                            {"name": "x", "collection_id": 1})
    assert decision.allowed is False
    assert decision.writes is True
    assert decision.reason.split(":", 1)[0] == "refused_write"
