"""The containment panel un-labelled its most alarming events.

``security_events`` knew 9 kinds. The gate records 25. An unknown kind
falls back to ``("•", <raw kind string>)``, so the events nobody had added
a label for rendered as an unlabelled bullet — and those were the loud
ones: every escape from a locked scope, an MCP call dispatched without
permissions, a write into a read-only location, and a refusal the agent
tried to walk around. In the panel they read as less serious than an
ordinary denied command, which has a label and a glyph.

This test reads the kinds out of the gate itself, so a new kind added
without a label fails here instead of quietly rendering as a bullet.
"""

from __future__ import annotations

import ast
from pathlib import Path

import _package_ast
from delfin.agent import security_events


def _recorded_kinds() -> set[str]:
    """Every literal kind a module of ``delfin/agent`` records through
    ``_record_security_event``, read from the source. Every module, not
    only ``api_client``: ``file_confirm``, ``cli`` and ``exchange`` record
    kinds of their own, and a scan of one file let three of them render
    as an unlabelled bullet."""
    agent = Path(__file__).resolve().parents[1] / "delfin" / "agent"
    kinds: set[str] = set()
    for src in sorted(agent.rglob("*.py")):
        kinds |= _kinds_in(_package_ast.parse(src))
    return kinds


def _kinds_in(tree) -> set[str]:
    kinds: set[str] = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        name = (func.attr if isinstance(func, ast.Attribute)
                else func.id if isinstance(func, ast.Name) else "")
        if name != "_record_security_event" or not node.args:
            continue
        first = node.args[0]
        if isinstance(first, ast.Constant) and isinstance(first.value, str):
            kinds.add(first.value)
        elif isinstance(first, ast.IfExp):        # "self_mod" if x else "calc_edit"
            for branch in (first.body, first.orelse):
                if isinstance(branch, ast.Constant) and isinstance(branch.value, str):
                    kinds.add(branch.value)
    return kinds


def test_every_recorded_kind_has_a_label():
    recorded = _recorded_kinds()
    assert recorded, "no security events found — the scan is broken, not the code"
    missing = sorted(recorded - security_events.known_kinds())
    assert not missing, f"kinds recorded but not labelled: {missing}"


def test_the_loud_ones_are_labelled():
    known = security_events.known_kinds()
    for kind in ("denied_again", "locked_scope_read", "locked_scope_bash",
                 "locked_scope_symlink", "locked_scope_exec",
                 "locked_scope_widen", "locked_scope_parse",
                 "read_only_write", "calc_edit", "test_tamper",
                 "mcp_bash_no_perms", "mcp_write_no_perms",
                 "role_denied_mcp", "plan_mode_mcp"):
        assert kind in known, kind


def test_an_unknown_kind_still_renders():
    """Fail soft: an unlabelled kind is a test failure, never a crash."""
    security_events.clear()
    security_events.record("something_new", "bash", "detail")
    html = security_events.format_panel_html()
    assert "something_new" in html
    security_events.clear()


def test_a_labelled_kind_renders_its_label():
    security_events.clear()
    security_events.record("locked_scope_symlink", "bash", "cat shared/x")
    html = security_events.format_panel_html()
    assert "Locked scope: link out of the folder" in html
    assert "•" not in html
    security_events.clear()
