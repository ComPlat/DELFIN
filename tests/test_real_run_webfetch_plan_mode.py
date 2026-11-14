"""Real run: one row of the admission matrix driven through the real entry
path. The form `chat -p` actually runs for an ollama/API session is
``OpenAIClient.stream_message`` (delfin/agent/api_client.py builds it in
make_client for the ollama backend), so the driver here feeds a fake SDK
endpoint (the tests/test_tc1_promise_endpoint.py shape: only
``client.chat.completions.create`` is replaced; the chunk loop, the
tool-call accumulator and the whole admission machinery are production
code) and drives stream_message with plan-mode KitToolPermissions.

Loopback TCP is not an option in the agent harness: its seccomp cage
denies bind() and connect() outright (the pytest run is inside the same
cage), so the pty/``chat -p`` subprocess form of this driver cannot run
there. The SDK-level driver is the same admission path minus the CLI
argparse shell.

What the real session decides for web_fetch of an outside URL in plan
mode must match the matrix row: deny (plan mode is read-only).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace as NS

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tests"))
sys.path.insert(0, str(REPO))

from delfin.agent.api_client import (  # noqa: E402
    KitToolPermissions,
    OpenAIClient,
)

URL = "https://sink.example/collect"


def text_chunk(t, finish=None):
    return NS(choices=[NS(delta=NS(content=t, tool_calls=None,
                                   reasoning_content=None),
                          finish_reason=finish)], usage=None)


def tool_chunk(name, args, idx=0, finish=None, call_id="call_1"):
    return NS(choices=[NS(delta=NS(
        content=None,
        tool_calls=[NS(index=idx, id=call_id,
                       function=NS(name=name, arguments=args))],
        reasoning_content=None), finish_reason=finish)], usage=None)


def _fake_endpoint(events_by_round):
    """create(**kw) -> per-round generator of chunk objects."""
    rounds = iter(events_by_round)
    sent: list = []

    def create(**kw):
        sent.append(kw.get("messages"))
        def gen():
            for c in next(rounds):
                yield c
        return gen()
    return NS(chat=NS(completions=NS(create=create))), sent


def test_real_session_denies_webfetch_in_plan_mode(tmp_path):
    ws = tmp_path / "ws"
    perms = KitToolPermissions(workspace=ws, mode="plan")
    client = OpenAIClient(api_key="sk-test", model="probe-model",
                          base_url="http://probe.invalid/v1",
                          provider="kit", permissions=perms)
    fake, sent = _fake_endpoint([
        [tool_chunk("web_fetch", json.dumps({"url": URL})),
         text_chunk("", finish="tool_calls")],
        [text_chunk("done", finish="stop")],
    ])
    client.client = fake

    results = []
    notices = []
    for ev in client.stream_message("t", [{"role": "user",
            "content": f"do-webfetch {URL}"}], max_tokens=16):
        if ev.type == "tool_result":
            results.append(ev.tool_output)
        elif ev.type == "notice":
            notices.append(ev.text)

    assert len(sent) >= 2, "the endpoint never saw the follow-up round"
    assert results, "the web_fetch call never produced a tool_result"
    denial = results[0]
    assert json.loads(denial).get("error", "").startswith(
        "plan mode (read-only)"), denial[:200]
    # The refusal names the reason the matrix row records: the tool is not
    # on the plan-mode read-only list, refused before its arguments are read.
    assert "not on the plan-mode read-only tool list" in denial
    assert any("plan mode" in n.lower() and "refused" in n.lower()
               for n in notices), notices
    # The engine corrected itself afterwards: the follow-up round arrived
    # (the model was told the call was refused), and nothing fetched.
    assert any("web_fetch" in json.dumps(rounds) for rounds in sent[:1])
