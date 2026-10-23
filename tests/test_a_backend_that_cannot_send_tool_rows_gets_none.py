"""A restored history carries its tool rounds only to a client that can send them."""
from delfin.agent import turn_history as th


def _rows():
    return [{"role": "user", "content": "q"},
            {"role": "assistant", "content": "", "tool_calls": [{"id": "c1"}]},
            {"role": "tool", "tool_call_id": "c1", "content": "r"},
            {"role": "assistant", "content": "note", "tool_calls": [{"id": "c2"}]},
            {"role": "tool", "tool_call_id": "c2", "content": "r2"},
            {"role": "assistant", "content": "answer"}]


def test_tool_rows_and_call_keys_are_stripped():
    out = th.strip_tool_rows(_rows())
    assert [m["role"] for m in out] == ["user", "assistant", "assistant"]
    assert out[1]["content"] == "note" and "tool_calls" not in out[1]
    assert out[2]["content"] == "answer"


def test_restore_strips_for_a_client_without_the_flag():
    import inspect
    from delfin.agent.engine import AgentEngine
    src = inspect.getsource(AgentEngine.restore_state)
    assert "strip_tool_rows" in src
    assert '"KEEPS_TOOL_HISTORY"' in src
