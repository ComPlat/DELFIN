"""The protocol pressed through the only tool it was given.

A dashboard-protocol role drives the UI by writing `ACTION: /command`
lines. When a tool surface is advertised alongside, a model sometimes
emits the protocol as a structured call instead, and the repair path
exists for exactly that: recognise it, register the command on the text
channel, and return something constructive rather than an unknown-tool
or role-denial error.

It recognised two shapes, both keyed on the tool NAME — a tool called
`ACTION`, or the slash command used as the tool name. Measured at the
protocol level on kit.glm-5.3, 3 of 3, the third shape:

    read_file(path="/orca show")
    read_file(path="/memories")
    read_file(path="/jobs")

The slash command in the ARGUMENT of a real tool. The same request
without streaming returns the three ACTION lines as text; with a tool
surface and streaming it comes back as three file reads of paths that do
not exist, and `workflow_parallel_checks` failed 0/3 with the whole
answer being "Klar — alle drei auf einen Blick:".

Narrow on purpose. Only a path-shaped key, only a value that is entirely
a slash command with no second separator, and only for the roles that
drive the UI by text — which have no file tools at all, so the
alternative for such a call was always an error.
"""

from __future__ import annotations

import pytest

from delfin.agent import action_protocol as ap


@pytest.mark.parametrize("value", ["/orca show", "/memories", "/jobs",
                                   "/tab calc", "/done"])
def test_a_command_in_a_path_argument_is_the_protocol(value):
    assert ap.is_action_style_call("read_file", {"path": value}), value
    assert ap.extract_slash_command("read_file", {"path": value}) == value


def test_the_namespaced_form_too():
    assert ap.is_action_style_call(
        "mcp__delfin-docs__read_file", {"path": "/orca show"})


@pytest.mark.parametrize("value", [
    "/etc/passwd",
    "/home/user/notes.md",
    "/tmp/x/y",
    "delfin/agent/engine.py",
    "bookmarks.json",
    "",
])
def test_a_real_path_is_still_a_real_path(value):
    assert not ap.is_action_style_call("read_file", {"path": value}), value


def test_a_prose_key_does_not_become_a_command():
    """Only path-shaped keys. A query or a remembered note may quote a
    command without being one."""
    assert not ap.is_action_style_call(
        "search_docs", {"query": "/orca show"})
    assert not ap.is_action_style_call(
        "remember", {"text": "/tab calc wechselt den Tab"})


def test_the_two_older_shapes_still_hold():
    assert ap.is_action_style_call("ACTION", {"command": "/orca show"})
    assert ap.is_action_style_call("/orca show")
    assert not ap.is_action_style_call("read_file", {"path": "x.py"})
