"""One malformed ``timeout`` silently disabled every hook the user had.

``_merge_hooks`` computed the timeout as::

    float(c.get("timeout", 30)) / 1000.0
    if c.get("timeout", 0) and c.get("timeout", 0) > 100 else ...

``"abc" > 100`` raises TypeError, and the exception escapes the whole
merge -- so ONE typo in ANY settings layer took out every hook in every
layer, including the user's own blocking PreToolUse hooks. Every call
site wraps ``load_hooks`` in ``except Exception`` and carries on, so
the only symptom was hooks that stopped firing, which looks exactly
like hooks that had nothing to complain about.

Verified against the pre-change tree: ``"timeout": "abc"`` raised
``TypeError: '>' not supported between instances of 'str' and 'int'``.

A configuration error must cost the entry that has it and nothing else,
and it must be said out loud.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import hooks as H
from delfin.agent import workspace_trust as WT


def _ws(tmp_path, entries) -> "object":
    ws = tmp_path / "project"
    (ws / ".delfin").mkdir(parents=True, exist_ok=True)
    (ws / ".delfin" / "settings.json").write_text(json.dumps({
        "hooks": {"PreToolUse": [{"matcher": "", "hooks": entries}]}
    }), encoding="utf-8")
    WT.trust_workspace(ws, [WT.KIND_HOOKS], actor=WT.ACTOR_USER)
    return ws


@pytest.mark.parametrize("bad", ["abc", [], {}, "30s", "null"])
def test_a_bad_timeout_costs_only_its_own_entry(tmp_path, bad):
    ws = _ws(tmp_path, [
        {"type": "command", "command": "broken", "timeout": bad},
        {"type": "command", "command": "the user's own linter"},
    ])
    cfg = H.load_hooks(ws)
    assert [c.command for c in cfg.for_event("PreToolUse")] == [
        "the user's own linter"], (
        "one malformed entry still disables the hooks around it")


def test_the_dropped_entry_is_reported(tmp_path):
    """A hook that stops running is indistinguishable from a hook that
    found nothing wrong, so the drop has to be said."""
    ws = _ws(tmp_path, [
        {"type": "command", "command": "broken", "timeout": "abc"},
    ])
    cfg = H.load_hooks(ws)
    warning = " ".join(cfg.warnings)
    assert "timeout" in warning
    assert "broken" in warning
    assert "settings.json" in warning


def test_a_negative_timeout_is_dropped_too(tmp_path):
    ws = _ws(tmp_path, [{"type": "command", "command": "x", "timeout": -5}])
    cfg = H.load_hooks(ws)
    assert cfg.for_event("PreToolUse") == []
    assert "positive" in " ".join(cfg.warnings)


def test_the_users_own_blocking_hook_survives_a_repo_typo(tmp_path,
                                                          monkeypatch):
    """The case that makes this a security bug rather than a papercut: a
    settings file the user did not write could switch off the settings
    file they did."""
    home = tmp_path / "home" / ".delfin" / "settings.json"
    home.parent.mkdir(parents=True)
    home.write_text(json.dumps({
        "hooks": {"PreToolUse": [{"matcher": "", "hooks": [
            {"type": "command", "command": "block-secrets"}]}]}
    }), encoding="utf-8")
    monkeypatch.setattr(H, "_user_settings_path", lambda: home)
    ws = _ws(tmp_path, [{"type": "command", "command": "x", "timeout": "abc"}])
    cfg = H.load_hooks(ws)
    assert "block-secrets" in [c.command for c in cfg.for_event("PreToolUse")]


# ---------------------------------------------------------------------------
# The units the schema documents still work
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw,seconds", [
    (5000, 5.0),        # milliseconds, the canonical schema
    (5, 5.0),           # seconds, the tolerated shorthand
    ("5000", 5.0),      # a numeric string is a number
    (None, 30.0),       # absent -> the default
    (0, 30.0),          # falsy -> the default, as before
])
def test_a_well_formed_timeout_is_unchanged(tmp_path, raw, seconds):
    entry = {"type": "command", "command": "x"}
    if raw is not None:
        entry["timeout"] = raw
    ws = _ws(tmp_path, [entry])
    assert H.load_hooks(ws).for_event("PreToolUse")[0].timeout_s == seconds
