"""Behavioural red tests: the four promise shapes, the length-cut
write_file and the control, driven through the REAL
``OpenAIClient.stream_message`` stream parser by a fake SDK endpoint.

The fake endpoint replaces only the SDK transport
(``self.client.chat.completions.create``); everything downstream -- the
SSE chunk loop, the tool-call accumulator, the DSML/Harmony recovery and
the whole turn loop with its corrective triggers -- is the production
code.
"""

from __future__ import annotations

from types import SimpleNamespace as NS

import pytest

from delfin.agent import agent_tasks
from delfin.agent.api_client import OpenAIClient, KitToolPermissions

SHAPE1 = ("I am executing the write_file tool call now -- creating "
          "`tools/adapters/orca_scan_opi.py` with its complete content. "
          "This is the actual tool invocation being executed.")
SHAPE2 = ("I am now emitting the actual tool call block. No prose -- "
          "just the invocation.")
SHAPE3 = ("I'll apply the `multi_edit` now with the exact baseline from "
          "lines 355-360 ...")
SHAPE4 = ("I'll run git status/log and the two greps: "
          "<system_warning>THIS IS THE ACTION BLOCK -- executing the "
          "batched tools.</system_warning>")
CONTROL = "All three checks pass; nothing else to do."

CORRECTIVE = ("No tool call arrived. Make the call now; if the file is "
              "large, write it in parts.")


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


def run_stream(client, streams):
    """Drive one real stream_message over a scripted fake endpoint.

    ``streams`` is a list of rounds; each round is a list of chunk
    objects the fake endpoint yields for one ``create`` call.
    """
    rounds = iter([list(r) for r in streams])
    sent: list[list] = []

    def create(**kw):
        sent.append([m for m in (kw.get("messages") or [])])
        def gen():
            try:
                for c in next(rounds):
                    yield c
            except StopIteration:
                return
        return gen()  # a generator has .close(), like the real SDK stream

    fake = NS(chat=NS(completions=NS(create=create)))
    client.client = fake
    text = ""
    stops: list[str] = []
    notices: list[str] = []
    tool_outputs: list[str] = []
    for ev in client.stream_message("t", [{"role": "user", "content": "x"}],
                                    max_tokens=16):
        if ev.type == "text_delta":
            text += ev.text
        elif ev.type == "message_delta":
            stops.append(str(ev.stop_reason))
        elif ev.type == "notice":
            notices.append(str(ev.text))
        elif ev.type == "tool_result":
            tool_outputs.append(str(ev.tool_output))
    return {"text": text, "stops": stops, "notices": notices,
            "tool_outputs": tool_outputs, "sent": sent}


def make_client(workspace, mode="acceptEdits"):
    perms = KitToolPermissions(workspace=workspace, mode=mode)
    return OpenAIClient(api_key="sk-test", model="deepseek-v4-flash",
                        base_url="http://probe.invalid/v1",
                        provider="kit", permissions=perms)


@pytest.fixture
def open_tasks_workspace(tmp_path):
    agent_tasks._STORES.clear()
    store = agent_tasks.get_store(tmp_path)
    store.create("write the adapter module", "step 1")
    yield tmp_path
    agent_tasks._STORES.clear()


class TestPromiseShapesAreCorrected:
    """Each of the four field shapes gets the corrective turn instead of
    ending the turn unrecovered."""

    @pytest.mark.parametrize("shape", [SHAPE1, SHAPE2, SHAPE3, SHAPE4])
    def test_a_promise_gets_the_corrective_turn(self, open_tasks_workspace,
                                                shape):
        client = make_client(open_tasks_workspace)
        # First round: the promise. Second round: the model makes the
        # call it promised (as a real tool call).
        out = run_stream(client, [
            [text_chunk(shape)],
            [text_chunk("done."),
             tool_chunk("write_file", '{"path": "a.txt", "content": "b"}',
                        finish="tool_calls")],
        ])
        assert any(CORRECTIVE in n for n in out["notices"]), out["notices"]
        # The corrective row also went to the MODEL (not just the chat).
        assert any(
            CORRECTIVE in str(m.get("content") or "")
            for m in out["sent"][-1] if m.get("role") == "user"), out["sent"][-1]
        # A tool call ran after the corrective turn.
        assert out["tool_outputs"], out
        assert len(out["stops"]) == 1, out["stops"]

    def test_repeated_promises_are_corrected_at_most_twice_then_reported(
            self, open_tasks_workspace):
        """A model that keeps promising gets at most _PROMISED_CONT_MAX
        corrective rows, then the turn ends with the clear non-done stop
        reason a supervisor resumes."""
        client = make_client(open_tasks_workspace)
        promise = [text_chunk(SHAPE2)]
        out = run_stream(client, [promise] * 6)
        corrective = sum(1 for n in out["notices"] if CORRECTIVE in n)
        assert corrective == 2, corrective
        assert out["stops"][-1] == "promised_no_tool_call", out["stops"]
        # Only the two corrective rounds were spent: 3 rounds served
        # (promise, correction, promise), then the turn ended.
        assert len(out["sent"]) == 3, len(out["sent"])


class TestTruncatedWriteFile:
    def test_a_length_cut_write_file_tells_the_model(self,
                                                     open_tasks_workspace):
        """A write_file whose arguments are cut by max_tokens must come
        back to the model as a cut-call error naming smaller parts --
        not vanish."""
        client = make_client(open_tasks_workspace)
        cut_args = ('{"path": "big.txt", "content": "line one\\nline '
                    'two\\nline three\\nand the string just keeps going')
        out = run_stream(client, [
            [text_chunk("writing the file now."),
             tool_chunk("write_file", cut_args, finish="length")],
            [text_chunk("The file is not written yet.", finish="stop")],
        ])
        # The cut_call continuation row went to the MODEL with the brief's
        # wording (a harness row, not a chat notice).
        rows = [str(m.get("content") or "") for m in out["sent"][-1]
                if m.get("role") == "user"]
        joined = "\n".join(rows + out["notices"])
        assert "the arguments were cut off at the output limit" in joined, \
            out["sent"][-1]
        assert "write the file in smaller parts" in joined
        # The message reached the MODEL (verified above) and the turn did
        # not end before the continuation: two rounds were served.
        assert len(out["sent"]) >= 2, out["sent"]


class TestControl:
    def test_a_plain_answer_ends_normally_first_try(self,
                                                    open_tasks_workspace):
        """No promise: the turn ends on the first stop, no corrective
        round is spent."""
        client = make_client(open_tasks_workspace)
        out = run_stream(client, [[text_chunk(CONTROL, finish="stop")]])
        assert len(out["stops"]) == 1 and \
            out["stops"][0].startswith("end_turn"), out["stops"]
        assert not any(CORRECTIVE in n for n in out["notices"])
