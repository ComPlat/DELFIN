"""An auto-verify round showed the finished answer twice.

The notice went under the answer's bubble without closing it. The next
token found the notice at the end of the chat, opened a new bubble and
filled it with the whole buffer -- the answer again, with the new part
appended (report 20260915-084010).

The buffer is not simply emptied at the notice: twenty places read
``"".join(chunks)`` after the turn, and a notice that ends a turn must not
take the answer from them. It is marked as shown instead, and emptied when
new text arrives.
"""

from __future__ import annotations

import ast
import inspect
import pathlib

from delfin.dashboard import tab_agent as T

_SRC = pathlib.Path(inspect.getfile(T)).read_text(encoding="utf-8")


def _fn(name: str) -> ast.FunctionDef:
    for node in ast.walk(ast.parse(_SRC)):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise AssertionError(f"{name} not found")


def test_a_notice_closes_the_answer_before_it():
    src = ast.unparse(_fn("_on_notice"))
    assert "_update_last_assistant" in src
    assert src.index("_update_last_assistant") < src.index("_append_system_message")
    assert "answer_shown[0] = True" in src


def test_a_notice_does_not_empty_the_buffer():
    assert "chunks.clear()" not in ast.unparse(_fn("_on_notice"))


def test_text_after_a_notice_starts_its_own_bubble():
    src = ast.unparse(_fn("_on_token"))
    i = src.index("if answer_shown[0]:")
    assert "chunks.clear()" in src[i:i + 120]
    assert i < src.index("chunks.append(text)")


def test_nothing_draws_a_shown_answer_a_second_time():
    for name in ("_on_tool_use", "_on_permission_denied"):
        src = ast.unparse(_fn(name))
        assert "if not answer_shown[0]:" in src, name
        assert "answer_shown[0] = False" in src, name
    assert ('if chunks and not answer_shown[0]:\n'
            '                        _update_last_assistant("".join(chunks), '
            'role_label, finalize=True)') in _SRC
