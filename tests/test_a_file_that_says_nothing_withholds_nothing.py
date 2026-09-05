"""Two numbers about one settings file, and what each of them decides.

``Kind.count`` is how many definitions inside the file would *execute* if it
were honoured.  It is the number the notice quotes -- "3 hook commands were
NOT loaded" -- so a definition that is switched off has no business in it.

``Kind.declares`` is how many the file *names*, switched off or not.  It is a
different question and it decides a different thing: whether the file says
anything at all, and so whether there is a decision for the user to make.

They agree everywhere except at one point, and that point is a feature.  A
workspace's ``mcp_servers.json`` may switch a built-in server off -- disabling
is a tightening, which a project is entitled to do, while redefining a built-in
is refused.  Such a file runs nothing and is emphatically not silent.  Asked
with ``count``, it looked empty; the gate dropped it as offering nothing, the
loader never saw it, and the server it had turned off came back on.

The other direction is why the question is asked at all: a ``settings.json``
that merely exists, with an empty hooks block or none, met the user with a
paragraph beginning "0 hook commands ... were NOT loaded".  A warning about a
risk that is not there is what teaches people to scroll past the one that
matters.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import mcp_client as M
from delfin.agent import workspace_trust as WT


def _write(ws, rel, body):
    path = ws / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(body), encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Nothing to say
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("body", [
    {},
    {"hooks": {}},
    {"_meta": {"written": "2026-08-25"}},
    {"hooks": {"PreToolUse": []}},
])
def test_a_settings_file_with_no_hooks_in_it_asks_the_user_nothing(
        tmp_path, body):
    """It exists, and that is all it does."""
    ws = tmp_path / "proj"
    _write(ws, ".delfin/settings.json", body)
    decision = WT.gate(WT.KIND_HOOKS, ws)
    assert decision.state == WT.STATE_NOTHING_OFFERED
    assert not decision.notice
    assert not decision.short_notice


# ---------------------------------------------------------------------------
# Something to say, and nothing to run
# ---------------------------------------------------------------------------

def test_a_file_whose_only_entry_switches_a_server_off_is_not_silent(tmp_path):
    """Nothing in it would execute, and it is still a decision."""
    ws = tmp_path / "proj"
    _write(ws, ".delfin/mcp_servers.json",
           {"servers": {"delfin-tools": {"enabled": False}}})
    decision = WT.gate(WT.KIND_MCP_SERVERS, ws)
    assert decision.state != WT.STATE_NOTHING_OFFERED
    assert decision.withheld, 'the file is a statement and it was withheld'


def test_the_two_numbers_are_the_two_questions(tmp_path):
    """One file, both entries: what runs is one, what it names is two."""
    ws = tmp_path / "proj"
    _write(ws, ".delfin/mcp_servers.json",
           {"servers": {"delfin-tools": {"enabled": False},
                        "mine": {"command": "python", "args": ["-m", "mine"]}}})
    offer, = WT.offers(WT.KIND_MCP_SERVERS, ws)
    assert offer.count == 1, 'the switched-off one would not execute'
    assert offer.declares == 2, 'and it is still named'


def test_the_notice_counts_what_would_run_and_not_what_is_switched_off(
        tmp_path):
    """The number in the sentence is the risk, so it is the runnable one."""
    ws = tmp_path / "proj"
    _write(ws, ".delfin/mcp_servers.json",
           {"servers": {"off_one": {"enabled": False},
                        "off_two": {"enabled": False},
                        "mine": {"command": "python"}}})
    decision = WT.gate(WT.KIND_MCP_SERVERS, ws)
    assert decision.withheld_count == 1, decision.notice


def test_a_disabled_builtin_survives_the_gate_and_reaches_the_loader(tmp_path):
    """The whole point of the distinction, asked where the user feels it.

    This is the same claim as ``test_a_project_may_disable_a_builtin`` in
    tests/test_a_folder_cannot_start_a_process.py, kept here as well because
    that one asserts on the layering and this one on the reason the layering
    ever gets the chance to run.
    """
    ws = tmp_path / "proj"
    _write(ws, ".delfin/mcp_servers.json",
           {"servers": {"delfin-tools": {"enabled": False}}})
    WT.trust_workspace(ws, [WT.KIND_MCP_SERVERS], actor=WT.ACTOR_USER)
    assert "delfin-tools" not in M._load_configs(ws)


# ---------------------------------------------------------------------------
# A link out of the workspace counts however empty it is
# ---------------------------------------------------------------------------

def test_an_escaping_link_is_the_finding_whatever_it_points_at(tmp_path):
    """There the link itself is what the user is being told about."""
    outside = tmp_path / "elsewhere" / "settings.json"
    outside.parent.mkdir(parents=True, exist_ok=True)
    outside.write_text(json.dumps({"hooks": {}}), encoding="utf-8")
    ws = tmp_path / "proj"
    (ws / ".delfin").mkdir(parents=True, exist_ok=True)
    try:
        (ws / ".delfin" / "settings.json").symlink_to(outside)
    except (OSError, NotImplementedError):
        pytest.skip('this filesystem does not make symlinks')
    decision = WT.gate(WT.KIND_HOOKS, ws)
    assert decision.state != WT.STATE_NOTHING_OFFERED
    assert decision.withheld and decision.withheld[0].escapes
