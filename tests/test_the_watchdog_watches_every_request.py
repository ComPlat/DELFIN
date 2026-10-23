"""A request that never started held a solo turn for 934 s.

The first-token budget was judged by a flag set once per TURN: after the
turn's first token, every later request -- the one after each tool result
-- counted as mid-stream, and in solo mode mid-stream has no budget. So the
kill watch returned for good, and a request the endpoint never started
waited until the gateway cut it (report 20260915-084010). The flag is per
request now, and where there is no budget the watch waits for the next
first-token wait instead of ending.
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


def test_the_kill_asks_about_this_request_not_the_turn():
    src = ast.unparse(_fn("_check_kill"))
    assert "waiting_for_first = not state.get('_request_saw_output')" in src


def test_a_tool_result_starts_a_new_first_token_wait():
    body = _fn("_on_tool_result").body
    reset_at = min(i for i, stmt in enumerate(body)
                   if "_request_saw_output" in ast.dump(stmt))
    return_at = min(i for i, stmt in enumerate(body)
                    if isinstance(stmt, ast.If)
                    and any(isinstance(s, ast.Return) for s in stmt.body))
    assert reset_at < return_at


def test_any_output_ends_the_first_token_wait():
    assert "_request_saw_output" in ast.unparse(_fn("_mark_output"))


def test_no_budget_mid_stream_does_not_end_the_watch():
    for node in ast.walk(_fn("_check_kill")):
        if isinstance(node, ast.If) and ast.unparse(node.test) == "budget <= 0":
            body = ast.unparse(node)
            assert "_threading.Timer" in body and ".daemon = True" in body
            return
    raise AssertionError("the budget check is gone from _check_kill")


def test_a_stall_after_work_does_not_say_the_turn_cost_nothing():
    src = ast.unparse(_fn("_check_kill"))
    assert "The rounds completed so far are kept" in src
    assert "No tokens were produced, so this turn cost nothing." in src
