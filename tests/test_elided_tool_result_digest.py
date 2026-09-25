"""Integration contract: an elided tool result leaves a digest behind.

The elision itself lives in api_client (security-owned; the operator
builds the wiring). This test fixes the CONTRACT first: when
_elide_old_tool_results replaces an old tool result, the replacement
must carry the digest of that call (tool name + the facts from
delfin.agent.tool_digest), not only the bare placeholder.

Until the wiring exists, this is xfail(strict=True): it must FAIL on the
current commit and will turn XPASS->failure the moment the wiring lands,
so the integration cannot regress silently.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent.api_client import _elide_old_tool_results


def _api_messages_with_read(path: str, body: str) -> list[dict]:
    """The exact wire shape api_client builds: an assistant message with
    tool_calls, followed by the matching tool result."""
    call_id = "call_1"
    return [
        {"role": "user", "content": "read it"},
        {
            "role": "assistant",
            "content": None,
            "tool_calls": [{
                "id": call_id,
                "type": "function",
                "function": {
                    "name": "read_file",
                    "arguments": json.dumps({"path": path}),
                },
            }],
        },
        {"role": "tool", "tool_call_id": call_id, "content": body},
        # one recent tool result so keep_recent leaves the old one alone
        {"role": "assistant", "content": None, "tool_calls": [{
            "id": "call_2", "type": "function",
            "function": {"name": "list_files",
                         "arguments": json.dumps({"pattern": "*"})},
        }]},
        {"role": "tool", "tool_call_id": "call_2", "content": "recent"},
    ]


@pytest.mark.xfail(strict=True, reason=(
    "elision wiring not built yet: the placeholder must gain the "
    "tool_digest digest of the elided call"))
def test_elided_tool_result_carries_its_digest():
    path = "delfin/somewhere/deep.py"
    body = "\n".join(f"{i + 1}  line {i}" for i in range(400))
    messages = _api_messages_with_read(path, body)
    elided = _elide_old_tool_results(
        messages, char_budget=200, keep_recent=1)
    assert elided >= 1, "the old result should have been elided"
    replaced = str(messages[2]["content"])
    assert path in replaced, (
        "the elided copy must still say WHICH file was read")
    assert "read_file" in replaced
