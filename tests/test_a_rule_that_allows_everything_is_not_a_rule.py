"""remember_permission would persist an auto-allow matching every command.

`extra_dir` refuses `/` however the confirm dialog is answered — "it is
your home directory or a system root, which would let the agent write
anywhere". That refusal exists because approval in the moment is not
informed consent to an unbounded grant.

`allow_pattern` had the same hole one field over. A rule of `^.*$` skips
the confirmation for EVERY shell command, in this session and in every
later one that reads the settings file, and the model composes both the
pattern and the rationale the user reads on the dialog — so "so the build
commands don't interrupt you" is a plausible caption for a total bypass.
Driven with a confirm callback that answers yes, `^.*$`, `.*` and `.+`
all persisted.

The hazard was already written down in the executor, in the comment above
the scope_locked branch: "allow_pattern can persist a rule that
auto-approves every shell command from now on. Neither is a decision this
session may make for the ones after it." That refusal only fired for a
LOCKED scope, which is not the session most turns run in.

Found by driving a tool that had never been called: remember_permission
is one of 33 advertised tools with zero calls across 2509 recorded
benchmark runs.

How the check works, and why it is not a list of spellings: the question
asked is not "is this regex equivalent to .*", which is not worth
deciding, but "does this pattern separate ANY two commands". A real rule
names something and fails on unrelated input. Every spelling of a
catch-all passes all five probes at once.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest


@pytest.fixture
def session(tmp_path, monkeypatch):
    """A workspace and a user who approves everything — the worst case."""
    from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", lambda: home)

    class _AlwaysYes:
        last_timed_out = False

        def confirm(self, *a, **k):
            return True

    workspace = tmp_path / "ws"
    workspace.mkdir()
    perms = KitToolPermissions(workspace=str(workspace))
    perms.mode = "acceptEdits"
    perms.task_session_id = "unbounded-allow"
    perms.confirm_callback = _AlwaysYes().confirm
    return _DocToolExecutor(), perms


def _remember(session, value, kind="allow_pattern"):
    executor, perms = session
    return json.loads(executor.execute(
        "remember_permission",
        {"kind": kind, "value": value, "scope": "repo",
         "rationale": "keeps the build from interrupting you"}, perms))


@pytest.mark.parametrize("catch_all", [
    "^.*$", ".*", ".+", r"[\s\S]*", "(?s).*", ".", "^", r"(?:.|\n)*",
])
def test_a_pattern_that_matches_everything_is_refused(session, catch_all):
    out = _remember(session, catch_all)
    assert "error" in out, f"{catch_all!r} was persisted: {out}"
    assert "every command" in out["error"]


@pytest.mark.parametrize("real_rule", [
    r"^git (status|diff|log)\b",
    r"^pytest\b",
    "^python -m pytest",
    "^npm (run )?build$",
    r"^\s*ruff\s+(?:check|format)\b",
])
def test_a_rule_that_names_what_it_allows_still_persists(session, real_rule):
    out = _remember(session, real_rule)
    assert out.get("status") == "persisted", out
    assert out["value"] == real_rule


def test_the_refusal_survives_a_user_who_says_yes(session):
    """The point of the floor. The dialog is not the safeguard here,
    because the model writes the sentence the dialog shows."""
    executor, perms = session
    assert perms.confirm_callback("anything") is True, "premise gone"
    assert "error" in _remember(session, "^.*$")


def test_nothing_reached_the_settings_file(session, tmp_path):
    _remember(session, "^.*$")
    written = list((tmp_path / "ws").rglob("settings.json"))
    if written:
        body = json.loads(written[0].read_text(encoding="utf-8"))
        blob = json.dumps(body)
        assert "^.*$" not in blob, f"the refused rule was written: {blob}"


def test_a_deny_rule_is_left_alone(session):
    """A catch-all DENY fails safe — it refuses everything rather than
    allowing it — so the floor does not apply to it."""
    out = _remember(session, "^.*$", kind="deny_pattern")
    assert out.get("status") == "persisted", out


def test_the_detector_answers_the_narrow_question():
    from delfin.agent.api_client import _DocToolExecutor

    check = _DocToolExecutor._pattern_constrains_nothing
    assert check(".*") and check("^.*$") and check(".")
    assert not check(r"^git\b")
    assert not check("^rm -rf /$"), "a specific dangerous rule is not a catch-all"
    # An unparseable regex is not this check's business; the settings
    # loader is what refuses it.
    assert not check("([unclosed")


def test_the_shipped_bundle_carries_no_catch_all():
    """remember_permission_bundle takes a profile name from an enum, not a
    pattern, so a model cannot inject one — but the shipped list should
    hold the same line."""
    from delfin.agent.api_client import _DocToolExecutor

    for profile, rules in _DocToolExecutor._BUNDLE_PROFILES.items():
        for rule in rules:
            pattern = rule if isinstance(rule, str) else (
                rule.get("value") or rule.get("pattern") or "")
            if not isinstance(pattern, str) or not pattern:
                continue
            assert not _DocToolExecutor._pattern_constrains_nothing(pattern), (
                f"{profile} ships a catch-all: {pattern!r}")
