"""A wake-up is not the user speaking, and should not look like it.

A finished job, a scheduled run and a message from another session all
start a turn through the same input box a person types into
(``_send_on_its_own`` fills it and presses send). In the transcript that
made them indistinguishable from something the user wrote — a user read
one back and asked why the dashboard was quoting them.

The message itself has said so since the wake notice was fixed. This is the
other half: the bubble.

What does NOT change is the role. ``/retry`` re-sends the last message
whose role is "user"; the citation check reads the user's own words for
the paths a task names; the turn counter counts them. A different role
would quietly change all three, so the origin travels as a field beside
the role, and only the rendering reads it.
"""

from __future__ import annotations

import ast
import inspect

from delfin.dashboard import tab_agent as T


def _nested_source(name: str) -> str:
    tree = ast.parse(inspect.getsource(T))
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return ast.unparse(node)
    raise AssertionError(f"{name} is not defined in tab_agent")


def _render():
    ns: dict = {"_md_to_html": lambda s: f"<p>{s}</p>",
                "_assistant_turn_is_wordless": lambda m: False}
    import html as _html
    ns["_html"] = _html
    exec(_nested_source("_render_single_msg"), ns)
    return ns["_render_single_msg"]


def test_a_typed_message_is_still_a_user_bubble():
    out = _render()({"role": "user", "content": "do the thing"})
    assert "delfin-chat-user" in out
    assert "do the thing" in out


def test_a_turn_nobody_typed_is_rendered_as_an_event():
    out = _render()({"role": "user", "content": "[watch] a job finished",
                     "origin": "event"})
    assert "delfin-chat-user" not in out, (
        "a wake-up rendered as the user's own words: " + out)
    assert "delfin-chat-system" in out
    assert "a job finished" in out, "and it still says what happened"


def test_the_role_is_never_compared_to_the_origin():
    """The flows that key on the role must not be able to tell the
    difference — that is why the origin is a field, not a role."""
    tree = ast.parse(_nested_source("_render_single_msg"))
    for node in ast.walk(tree):
        if not isinstance(node, ast.Compare):
            continue
        against = {c.value for c in node.comparators
                   if isinstance(c, ast.Constant)}
        if "event" not in against:
            continue
        left = node.left
        assert isinstance(left, ast.Call) and getattr(
            left.func, "attr", "") == "get", ast.unparse(node)
        assert left.args and getattr(left.args[0], "value", "") == "origin", (
            "\"event\" is being used as a role: " + ast.unparse(node))


def test_every_send_site_marks_a_self_started_turn():
    """Four places append what was sent. A wake-up can arrive through any
    of them, so a site that forgets the mark shows a bubble again."""
    src = inspect.getsource(T)
    appends = src.count('_append_chat_message(\n')
    marked = src.count('"origin": "event"')
    assert marked >= 4, (
        f"only {marked} send sites mark a self-started turn")
    assert appends >= marked
