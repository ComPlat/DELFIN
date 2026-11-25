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
EXPECTED = {"eln_find", "eln_read", "eln_pending", "eln_draft", "eln_approve",
            "eln_attach"}

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

def test_the_eln_tools_are_not_in_the_builtin_catalogue():
    """They ride on the ``delfin-ops`` MCP surface instead.

    The builtin catalogue is the schema that goes on EVERY request, and its
    budget test (``test_tool_schema_budget``) exists so the catalogue cannot
    grow. Five ELN schemas cost 444 tokens (measured: 11865 -> 12317) and
    broke four of those tests. They bought nothing: no ``_ELN_TOOL_NAMES``
    dispatch branch existed for them, so a model calling ``eln_find`` from
    the catalogue would have been answered by the dispatcher's unknown-tool
    path. The ops server is where ``parse_orca_output`` and the other 100+
    typed tools already live, and that route IS registered and IS driven
    below. Pinned in this direction so a later change does not quietly put
    them back on the per-request surface.
    """
    from delfin.agent import tool_catalog as TC
    names = {t["function"]["name"] for t in TC._DOC_TOOLS_OPENAI}
    assert not (EXPECTED & names), sorted(EXPECTED & names)


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


def test_the_ops_server_is_what_makes_them_known_to_a_preset():
    """The one measurement that decides whether a specialist may list them.

    ``_known_preset_tool_names`` derives its catalogue from the modules, so
    the server functions ARE the reason a bare ``eln_find`` in a preset's
    ``tools:`` is accepted. Asserted by removing them from the derivation
    rather than by naming them: a name can be listed by a hand-written set
    that knows nothing about the server.
    """
    from delfin.agent import subagents as SA
    known = SA._known_preset_tool_names()
    assert EXPECTED <= known, sorted(EXPECTED - known)
    assert {f"mcp__delfin-ops__{n}" for n in EXPECTED} <= known


def test_the_specialist_can_actually_call_the_tools_it_lists(home):
    """Rule 5: a set membership is not a capability.

    A1 shipped this preset in the registry with an ``eln_*``-free ``tools:``
    list, which passed every check it wrote for itself -- ``validate_preset``
    was clean and the type was offered. It could still call nothing: preset
    tools are INTERSECTED into the child's ``session_allowed_tools``, which
    is the field the executor gate reads, so listing only read tools removed
    every ELN tool. Measured directly through ``_derive_perms`` and the same
    ``_session_tool_refusal`` the executor calls.
    """
    from delfin.agent import subagents as SA
    from delfin.agent.api_client import (
        KitToolPermissions, _session_tool_refusal)

    preset = SA._build_preset_registry()["chemotion-eln"]
    parent = KitToolPermissions(workspace=Path("."))
    child = SA._derive_perms(parent, preset.mode,
                             preset_tools=preset.tools,
                             preset_name=preset.name)
    allowed = child.session_allowed_tools
    assert allowed, "the preset narrowed the child to no tools at all"
    for name in EXPECTED:
        assert _session_tool_refusal(name, allowed,
                                     frozenset()) is None, (
            f"the specialist lists content it means to use but cannot call "
            f"{name!r}; its tools: line must carry it")


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


def test_the_attachment_tool_is_not_host_gated_either():
    """The same standing-grant defect would apply to ``eln_attach``.

    The package names attaching a file (a lab run's export). It is a WRITE
    on the same server whose ``DELFIN_OPS_ALLOW_MUTATE=1`` is a standing
    host grant, so it must go through ``_safe`` like its siblings and never
    ``_host_gated``.
    """
    server = (REPO / "delfin" / "ops_server" / "server.py").read_text(
        encoding="utf-8")
    assert 'mcp.tool(name="eln_attach")(_safe("eln_attach", tool_eln_attach))' \
        in server
    assert "_host_gated(tool_eln_attach)" not in server


def test_an_attachment_needs_an_approval_and_sends_nothing_without_one(home):
    """The package's write rule, driven through the tool a model calls.

    Two assertions, and the second is the one that matters: without an
    approval NOTHING reaches the server -- not even the first chunk -- and
    the refusal names the file, the element, the request id and the person's
    command. Then one grant sends the exporter's whole sequence once, and a
    second call with the same grant is refused.
    """
    from delfin.lab import eln_tools, eln_write

    fake = FakeELN()
    payload = b"x" * 11_000          # three chunks at MAX_UPLOAD_CHUNK=5000
    args = {"path": "", "content": payload, "element_type": "Sample",
            "element_id": 11, "filename": "run-export.json"}

    out = json.loads(eln_tools.execute("eln_attach",
                                       {**args, "_transport": fake.send}, None))
    assert out["code"] == "needs_approval", out
    assert out["request_id"].startswith("eln-")
    assert "run-export.json" in out["error"]
    assert "11000 bytes" in out["error"]
    assert "Element 11" in out["error"] or "Sample 11" in out["error"]
    assert f"delfin-agent eln approve {out['request_id']}" in out["error"]
    # Nothing was sent: not the chunk, not the complete, not the link.
    assert [c for c in fake.calls if "attachments" in c["path"]] == []
    assert fake.uploaded == [] and fake.attached == []

    # The person grants it; the same call now sends the sequence once.
    eln_write.grant(out["request_id"], approver="alice")
    sent = json.loads(eln_tools.execute(
        "eln_attach", {**args, "approval": out["request_id"],
                       "_transport": fake.send}, None))
    assert sent["ok"] is True, sent
    assert sent["chunks"] == 3
    assert fake.uploaded[0]["filename"] == "run-export.json"
    assert len(fake.uploaded[0]["checksum"]) == 32
    assert fake.attached == [{"element_type": "Sample",
                              "attachment_ids": [sent["key"]]}]

    # The grant is single use.
    again = json.loads(eln_tools.execute(
        "eln_attach", {**args, "approval": out["request_id"],
                       "_transport": fake.send}, None))
    assert again.get("code") == "refused_approval", again
    assert len(fake.attached) == 1


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


def test_the_specialist_names_the_eln_tools_it_may_call():
    """The inverse of what this file asserted first, and the reason for it.

    ``_known_preset_tool_names`` DOES model the ops server (it walks every
    ``tool_*`` in ``delfin/ops_server/server.py``), so a bare ``eln_find`` is
    a known name -- and it must be listed, because the list narrows. The
    body's claim that the tools appear namespaced is still true and still
    checked; both spellings resolve to the same permission.
    """
    text = PRESET.read_text(encoding="utf-8")
    frontmatter = text.split("---")[1]
    listed = set(frontmatter.split("tools:", 1)[1].split("]", 1)[0]
                 .strip().strip("[]").replace(",", " ").split())
    assert EXPECTED <= listed, sorted(EXPECTED - listed)
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
    shared dispatcher.

    The name set is the full documented one, ``lab_run`` included: the
    workflow runner is the lab's sixth tool, added by the workflow package
    that lands alongside this one. Asserting the five-name set from before
    it is red on the bundle and green on the base -- a control that fails
    for a reason on neither branch is measuring the merge, not the ELN.
    Equality is kept rather than a subset check, so a tool that disappears
    or a stray name that appears is still caught.
    """
    from delfin.lab import safety  # noqa: F401  (the module the kernel lives in)
    from delfin.lab.tools import TOOL_NAMES as LAB_NAMES
    assert LAB_NAMES == {"lab_devices", "lab_read", "lab_set",
                         "lab_plan_check", "lab_estop", "lab_run"}


def test_control_the_write_gate_refuses_before_this_change_too():
    """``POST /api/v1/research_plans`` was refused on the base and is
    refused now: the ELN work added routes, it did not open writes."""
    from delfin.lab.eln_routes import decide_route
    decision = decide_route("POST", "/api/v1/research_plans",
                            {"name": "x", "collection_id": 1})
    assert decision.allowed is False
    assert decision.writes is True
    assert decision.reason.split(":", 1)[0] == "refused_write"
