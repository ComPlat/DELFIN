"""Behavioral sandbox tests for weak-model output patterns.

We can't run an actual 7B local model in CI, but we can mock the
output patterns weak models routinely emit and assert the framework's
parsing + dispatcher handles them gracefully. This is the "autonomous
training" loop in miniature: feed the dispatcher synthetic weak-model
output, observe what gets extracted, optimise when a real model
exhibits a pattern we don't yet handle.

Coverage:
- Core-tool filter: weak-model schema is trimmed to 15 tools
- ACTION fault-tolerance: 5 variant formats parse to the same command
- Bare-slash false-positive guard: filesystem paths and URLs do NOT
  trigger the dispatcher
- /done sentinel: weak-model bare-slash form still terminates loop
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Core-tool filter (api_client _WEAK_MODEL_CORE_TOOLS)
# ---------------------------------------------------------------------------


def test_weak_model_core_tools_set_has_only_essentials():
    """The core set must cover read/write/edit/bash/grep + planning +
    subagent + ask_user_question — and nothing exotic that confuses
    weak attention (cron / worktree / push_notification / etc.)."""
    from delfin.agent.api_client import _WEAK_MODEL_CORE_TOOLS
    must_have = {
        "read_file", "write_file", "edit_file", "multi_edit",
        "grep_file", "list_files",
        "bash",
        "ask_user_question",
        "task_create", "task_update", "task_list",
        "subagent",
    }
    assert must_have.issubset(_WEAK_MODEL_CORE_TOOLS), (
        f"missing core tools: {must_have - _WEAK_MODEL_CORE_TOOLS}"
    )
    # And exotic ones must NOT be in the weak surface
    must_exclude = {
        "cron_create", "cron_delete", "cron_list",
        "enter_worktree", "exit_worktree",
        "schedule_wakeup", "remote_trigger", "push_notification",
        "notebook_read", "notebook_edit",
        "apply_patch",
    }
    assert must_exclude.isdisjoint(_WEAK_MODEL_CORE_TOOLS), (
        f"weak surface leaked exotic tools: "
        f"{must_exclude & _WEAK_MODEL_CORE_TOOLS}"
    )


def test_weak_model_core_tools_is_at_most_20():
    """Hard cap: even with future additions the weak surface must
    stay small enough for 7B-tier models to route reliably."""
    from delfin.agent.api_client import _WEAK_MODEL_CORE_TOOLS
    assert len(_WEAK_MODEL_CORE_TOOLS) <= 20, (
        f"weak surface grew to {len(_WEAK_MODEL_CORE_TOOLS)} tools — "
        "weak models stop routing reliably past ~15."
    )


def test_core_filter_path_present_in_source():
    """The filter must actually be applied — guard against accidental
    removal in a future cleanup."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "agent" / "api_client.py")
    text = p.read_text(encoding="utf-8")
    idx = text.find("Weak-model core-tool filter")
    assert idx > 0, "core-tool filter block missing"
    snippet = text[idx: idx + 1500]
    assert "PromptLoader" in snippet
    assert "_is_weak_model" in snippet
    assert "_WEAK_MODEL_CORE_TOOLS" in snippet


# ---------------------------------------------------------------------------
# ACTION fault-tolerance — five variants must all parse to /tab calc
# ---------------------------------------------------------------------------


_ACTION_VARIANTS = [
    "ACTION: /tab calc",        # canonical
    "ACTION:/tab calc",          # no space after colon
    "ACTION /tab calc",          # no colon
    "Action: /tab calc",         # lower-case keyword
    "/tab calc",                 # bare slash
]


@pytest.mark.parametrize("variant", _ACTION_VARIANTS)
def test_action_cmd_parses_all_variants(variant):
    """All five variants must yield the same canonical command text."""
    # Mirror the production helper here so the test runs without
    # spinning up the dashboard widget tree.
    import re as _re
    KNOWN = frozenset((
        "/tab", "/control", "/orca", "/submit", "/recalc", "/cancel",
        "/calc", "/analyze", "/jobs", "/ui", "/workspace", "/done",
    ))
    rx = _re.compile(
        r"^\s*(?:ACTION\s*[:.]?\s*)?(/\S+.*?)\s*$", _re.IGNORECASE,
    )

    def _action_cmd(ln: str) -> str:
        s = ln.strip()
        m = rx.match(s)
        if not m:
            return ""
        cmd = m.group(1).strip()
        if s.upper().startswith("ACTION"):
            return cmd
        first = cmd.split()[0].lower()
        return cmd if first in KNOWN else ""

    assert _action_cmd(variant) == "/tab calc"


def test_action_cmd_rejects_filesystem_paths():
    """Critical guard: a line like '/home/user/x.py' must NOT be
    parsed as a slash-command — false positives would route prose
    through the dispatcher."""
    import re as _re
    KNOWN = frozenset((
        "/tab", "/control", "/orca", "/calc", "/done",
    ))
    rx = _re.compile(
        r"^\s*(?:ACTION\s*[:.]?\s*)?(/\S+.*?)\s*$", _re.IGNORECASE,
    )

    def _action_cmd(ln: str) -> str:
        s = ln.strip()
        m = rx.match(s)
        if not m:
            return ""
        cmd = m.group(1).strip()
        if s.upper().startswith("ACTION"):
            return cmd
        first = cmd.split()[0].lower()
        return cmd if first in KNOWN else ""

    for bare in [
        "/home/user/file.py",                       # absolute path
        "/etc/hosts",                                # system path
        "/var/log/syslog",                           # ditto
        "/usr/local/bin/python3",                    # binary path
        "/some/random/path/that/looks/like/cmd",     # generic
    ]:
        assert _action_cmd(bare) == "", f"false positive on {bare!r}"


def test_action_cmd_with_ACTION_prefix_accepts_unknown_slash():
    """When the ACTION: keyword is present, any /-prefixed string is
    accepted — the agent might legitimately want a slash-command we
    haven't pre-registered (e.g., a user-defined /command from
    ~/.delfin/commands/)."""
    import re as _re
    KNOWN = frozenset(("/tab",))
    rx = _re.compile(
        r"^\s*(?:ACTION\s*[:.]?\s*)?(/\S+.*?)\s*$", _re.IGNORECASE,
    )

    def _action_cmd(ln: str) -> str:
        s = ln.strip()
        m = rx.match(s)
        if not m:
            return ""
        cmd = m.group(1).strip()
        if s.upper().startswith("ACTION"):
            return cmd
        first = cmd.split()[0].lower()
        return cmd if first in KNOWN else ""

    assert _action_cmd("ACTION: /my-custom-command foo") == "/my-custom-command foo"


def test_helper_block_present_in_source():
    """Regression guard for the helper itself."""
    p = (Path(__file__).resolve().parent.parent
         / "delfin" / "dashboard" / "tab_agent.py")
    text = p.read_text(encoding="utf-8")
    assert "_ACTION_RE = _re.compile" in text
    assert "def _action_cmd(ln: str)" in text
    assert "_KNOWN_SLASH_PREFIXES" in text
    # /done detection now uses the helper
    assert "_action_cmd(ln).lower().startswith(\"/done\")" in text


# ---------------------------------------------------------------------------
# Multi-step weak-model output: all five variants in one response
# ---------------------------------------------------------------------------


def test_weak_model_mixed_action_formats_in_one_response():
    """Synthetic weak-model output that mixes canonical + bare-slash
    + missing-colon forms — the dispatcher must extract ALL of them."""
    import re as _re
    KNOWN = frozenset((
        "/tab", "/control", "/orca", "/submit", "/recalc", "/cancel",
        "/calc", "/analyze", "/jobs", "/ui", "/workspace", "/done",
    ))
    rx = _re.compile(
        r"^\s*(?:ACTION\s*[:.]?\s*)?(/\S+.*?)\s*$", _re.IGNORECASE,
    )

    def _action_cmd(ln: str) -> str:
        s = ln.strip()
        m = rx.match(s)
        if not m:
            return ""
        cmd = m.group(1).strip()
        if s.upper().startswith("ACTION"):
            return cmd
        first = cmd.split()[0].lower()
        return cmd if first in KNOWN else ""

    output = """\
Sure, here's what I'll do:

ACTION: /control key functional BP86
/tab submit
ACTION /done
"""
    cmds = [_action_cmd(ln) for ln in output.splitlines()]
    cmds = [c for c in cmds if c]
    assert cmds == [
        "/control key functional BP86",
        "/tab submit",
        "/done",
    ]
