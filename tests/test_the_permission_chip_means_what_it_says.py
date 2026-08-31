"""A permission rung must deliver what its label promises.

Reported from the field: Bypass was set and the session still asked for
approval, repeatedly. It was not a fault in the gate — the gate did
exactly what it was told. The dashboard mapped the ``all_free`` profile to
the CLI mode ``"auto"``, which ``_map_kit_permission_mode`` caps at
``acceptEdits``, and in ``acceptEdits`` every shell command that is not on
the auto-allow list is still put to the user.

Four statements described one mechanism and disagreed with it:

    the chip                "Bypass"
    the ladder docstring    "Bypass  asks nothing"
    the warning banner      "unrestricted access to files, shell commands"
    the mapping comment     "user chose all free -> skip permission prompts"

The measurement that made it concrete — a single grep runs silently, the
loop over nine result folders asks:

    grep ENERGY data.out                                    prompts=0
    for d in */ESD; do grep "…" "$d/S1.out"; done           prompts=1

A compound command is auto-allowed only when EVERY segment is, and a
``for`` loop does not decompose into such segments. Reading nine
calculation folders is a loop, so the rung that promised no prompts
produced one per folder.

The cap was defensible on its own terms. What it could not do is call
itself Bypass: a user who sets a control and is asked anyway learns the
control is decorative, and then stops reading labels that do matter.

So the rung now passes through, and these tests pin BOTH halves — that it
asks nothing, and that the three protections the banner still promises
(deny-list, sandbox, read-only archive) are untouched by it. The second
half is the one that matters in a year: "no prompts" must never quietly
become "no limits".
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from delfin.agent.api_client import (
    KitToolPermissions, _doc_executor, _map_kit_permission_mode,
)
from delfin.dashboard.tab_agent import (
    PROFILE_TO_CLI_PERM, _perm_options_for_mode,
)


class _Broker:
    """Stands in for KitConfirmBroker: records, then approves.

    Approving is deliberate — a rung that asks is caught by the RECORD,
    not by a denial, so the same scene measures prompts under every mode
    without changing what the tools do.
    """

    def __init__(self):
        self.asked: list = []
        self.last_timed_out = False

    def callback(self, *args, **kwargs):
        self.asked.append(args[0] if args else kwargs)
        return True


@pytest.fixture
def scene():
    """A workspace, a read-only archive beside it, and a live gate."""
    with tempfile.TemporaryDirectory(prefix="chip-ws-") as ws_dir, \
            tempfile.TemporaryDirectory(prefix="chip-arc-") as arc_dir:
        ws, arc = Path(ws_dir), Path(arc_dir)
        (ws / "data.out").write_text("FINAL SINGLE POINT ENERGY  -1.0\n")
        (arc / "kept.txt").write_text("original\n")

        def build(profile: str):
            broker = _Broker()
            perms = KitToolPermissions(
                mode=_map_kit_permission_mode(PROFILE_TO_CLI_PERM[profile]),
                workspace=str(ws),
                read_only_workspace_dirs=[str(arc)],
                confirm_callback=broker.callback,
            )
            return perms, broker

        yield build, ws, arc


# The command from the report: read-only, compound, and the natural way to
# read nine result folders.
LOOP = 'for d in */ESD; do grep "FINAL SINGLE POINT ENERGY" "$d/S1.out"; done'


def _run(name, args, perms):
    return json.loads(_doc_executor.execute(name, args, perms))


# ---------------------------------------------------------------------------
# The label and the mechanism
# ---------------------------------------------------------------------------

def test_every_offered_rung_has_a_mapping():
    """A rung the user can select must resolve to a mode.

    ``.get(profile, "default")`` at the call site means a missing entry
    silently lands on the most restrictive rung instead of failing — the
    exact shape that let the drift live.
    """
    for label, value in _perm_options_for_mode("code"):
        assert value in PROFILE_TO_CLI_PERM, f"{label!r} maps to nothing"


def test_bypass_asks_nothing(scene):
    """The promise on the chip, measured against the gate."""
    build, _ws, _arc = scene
    perms, broker = build("all_free")
    _run("bash", {"command": LOOP}, perms)
    assert broker.asked == [], (
        "the rung labelled Bypass asked for approval")


def test_the_lower_rungs_still_ask(scene):
    """Guards the fix from becoming a blanket opening.

    If this ever goes green by accident, Bypass has stopped being a
    distinct rung and the ladder has collapsed into one setting.
    """
    build, _ws, _arc = scene
    perms, broker = build("repo_free")
    _run("bash", {"command": LOOP}, perms)
    assert broker.asked, "Accept Edits stopped asking before a shell command"


def test_no_rung_is_downgraded_on_the_way_to_the_gate(scene):
    """The mapper is the other place the setting could quietly change."""
    for profile, cli_perm in PROFILE_TO_CLI_PERM.items():
        assert _map_kit_permission_mode(cli_perm) == cli_perm, profile


# ---------------------------------------------------------------------------
# What "asks nothing" must never come to mean
# ---------------------------------------------------------------------------

def test_the_archive_stays_read_only_under_bypass(scene):
    """The exception the warning banner names, by all three routes.

    Write, edit and a shell redirection are separate code paths; a fix
    that covered two of them would leave the third as the way in.
    """
    build, _ws, arc = scene
    perms, _broker = build("all_free")

    assert "error" in _run(
        "write_file", {"path": str(arc / "new.txt"), "content": "x"}, perms)
    assert "error" in _run(
        "edit_file", {"path": str(arc / "kept.txt"),
                      "old_string": "original", "new_string": "changed"}, perms)
    assert "error" in _run(
        "bash", {"command": f"echo x > {arc}/via_shell.txt"}, perms)

    assert sorted(p.name for p in arc.iterdir()) == ["kept.txt"]
    assert (arc / "kept.txt").read_text() == "original\n"


def test_the_deny_list_still_bites_under_bypass(scene):
    """No confirmation is not the same as no rule."""
    build, _ws, _arc = scene
    perms, broker = build("all_free")
    out = _run("bash", {"command": "rm -rf /tmp/definitely-not-here"}, perms)
    assert "error" in out
    assert broker.asked == [], "a denied command must not be offered instead"


def test_a_write_outside_every_root_is_still_refused(scene):
    """The sandbox is a boundary, not a prompt, so no mode lifts it."""
    build, _ws, _arc = scene
    perms, _broker = build("all_free")
    assert "error" in _run(
        "write_file", {"path": "/etc/delfin-should-never-exist",
                       "content": "x"}, perms)
