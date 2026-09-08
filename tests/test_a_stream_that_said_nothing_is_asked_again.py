"""A stream that ends having produced nothing is asked again, once.

Measured at the protocol level on 2026-09-07, kit.glm-5.3 with a tool
surface advertised: the streamed request returns an empty content channel
— no content, no tool call, no reasoning, finish_reason "stop" — while the
identical request without ``stream`` returns the answer. Twelve times out
of twelve, and independent of reasoning_effort. The engine reported
"[empty turn]" for a question the model had answered.

The cure already existed for a neighbouring case: a litellm proxy that
cannot stream a request 400s with a message naming __aiter__, and the
client retries it without streaming. Nothing connected that machinery to
"the stream ended empty".

A round that produced reasoning, or a tool call, is NOT this case and must
not be re-sent — it is a model working, and re-asking would double the
work and could double a side effect.
"""

import types

import pytest


class _Chunk:
    def __init__(self):
        self.usage = None
        self.choices = []


class _Msg:
    def __init__(self, content=None, tool_calls=None):
        self.content = content
        self.tool_calls = tool_calls or []


class _Choice:
    def __init__(self, message, finish_reason="stop"):
        self.message = message
        self.finish_reason = finish_reason


class _Resp:
    def __init__(self, message):
        self.usage = types.SimpleNamespace(
            prompt_tokens=11, completion_tokens=7)
        self.choices = [_Choice(message)]


def test_the_helper_reads_a_non_streaming_answer():
    from delfin.agent.api_client import _absorb_non_streaming

    content, calls, finish, (tin, tout, cached) = _absorb_non_streaming(
        _Resp(_Msg(content="ACTION: /tab calc")))
    assert content == "ACTION: /tab calc"
    assert calls == {}
    assert finish == "stop"
    assert (tin, tout) == (11, 7)


def test_the_helper_backfills_a_missing_tool_call_id():
    """The streaming accumulator does the same; a second copy of this
    logic would be a second place to forget it."""
    from delfin.agent.api_client import _absorb_non_streaming

    call = types.SimpleNamespace(
        id="", function=types.SimpleNamespace(
            name="read_file", arguments='{"path": "a.py"}'))
    _, calls, _, _ = _absorb_non_streaming(_Resp(_Msg(tool_calls=[call])))
    assert calls[0]["id"] == "ns_0"
    assert calls[0]["name"] == "read_file"
    assert calls[0]["arguments_parts"] == ['{"path": "a.py"}']


def test_the_helper_survives_a_response_with_nothing_in_it():
    from delfin.agent.api_client import _absorb_non_streaming

    empty = types.SimpleNamespace(usage=None, choices=[])
    assert _absorb_non_streaming(empty) == ("", {}, "", (0, 0, 0))


def test_the_retry_is_conditioned_on_an_empty_round():
    """Reasoning or a tool call means the model was working; re-asking
    would double the work and could double a side effect."""
    import inspect

    from delfin.agent import api_client

    src = inspect.getsource(api_client.OpenAIClient.stream_message)
    i = src.index('if (kwargs.get("stream") and not _text_chunks')
    condition = src[i:i + 200]
    assert "not _tool_calls" in condition
    assert "not _saw_reasoning" in condition


def test_reasoning_is_remembered_not_only_forwarded():
    """_saw_reasoning has to be set where the reasoning delta arrives, or
    the condition above reads False for a model that thinks out loud."""
    import inspect

    from delfin.agent import api_client

    src = inspect.getsource(api_client.OpenAIClient.stream_message)
    assert "_saw_reasoning = False" in src
    assert "_saw_reasoning = True" in src
    assert src.index("_saw_reasoning = False") < src.index("_saw_reasoning = True")


def test_the_two_paths_share_one_absorber():
    """The proxy-cannot-stream path and the empty-stream path must not
    grow separate copies."""
    import inspect

    from delfin.agent import api_client

    src = inspect.getsource(api_client.OpenAIClient.stream_message)
    assert src.count("_absorb_non_streaming(") == 2


def _stub_client(empty_stream: bool = True, reasoning: bool = False,
                 tool_call: bool = False):
    """A client whose stream says nothing (or something), and whose
    non-streaming call answers."""
    import threading
    import types

    from delfin.agent import api_client as ac

    seen = {"stream": 0, "plain": 0}

    class _Stream:
        def __iter__(self):
            if empty_stream:
                return iter(())
            delta = types.SimpleNamespace(
                content=None, tool_calls=None,
                reasoning_content="denke nach" if reasoning else None)
            if tool_call:
                delta.tool_calls = [types.SimpleNamespace(
                    index=0, id="c1", function=types.SimpleNamespace(
                        name="read_file", arguments='{"path":"a.py"}'))]
            chunk = types.SimpleNamespace(
                usage=None,
                choices=[types.SimpleNamespace(delta=delta,
                                               finish_reason=None)])
            return iter((chunk,))

        def close(self):
            pass

    class _Stub:
        class chat:
            class completions:
                @staticmethod
                def create(**kw):
                    if kw.get("stream"):
                        seen["stream"] += 1
                        return _Stream()
                    seen["plain"] += 1
                    msg = types.SimpleNamespace(
                        content="ACTION: /tab calc", tool_calls=[])
                    return types.SimpleNamespace(
                        usage=types.SimpleNamespace(prompt_tokens=5,
                                                    completion_tokens=3),
                        choices=[types.SimpleNamespace(
                            message=msg, finish_reason="stop")])

    c = ac.OpenAIClient.__new__(ac.OpenAIClient)
    c.client = _Stub()
    c.model = "kit.glm-5.3"
    c._provider = "kit"
    c._base_url = "https://ki-toolbox.scc.kit.edu/api/v1"
    c._api_key = "x"
    c.effort = ""
    c._permissions = None
    c.on_model_switched = None
    c._steer_lock = threading.Lock()
    c._steer_queue = []
    c._run_notes = []
    c._stop_flag = False
    return c, seen


def _first_text(client) -> str:
    """Drive the turn until the first answer text, then stop — the stub is
    only complete up to that point."""
    out = []
    try:
        for ev in client.stream_message(
                messages=[{"role": "user", "content": "öffne Calculations"}],
                system="You are a test.", max_tokens=64):
            if getattr(ev, "type", "") == "text_delta":
                out.append(ev.text)
                break
    except Exception:
        pass
    return "".join(out)


def test_an_empty_stream_is_retried_and_the_answer_comes_back():
    client, seen = _stub_client(empty_stream=True)
    assert _first_text(client) == "ACTION: /tab calc"
    assert seen["stream"] == 1
    assert seen["plain"] == 1, "the empty stream was not retried"


def test_a_round_that_reasoned_is_not_re_sent():
    """A model thinking out loud is working. Re-asking would double the
    work and could double a side effect."""
    client, seen = _stub_client(empty_stream=False, reasoning=True)
    _first_text(client)
    assert seen["plain"] == 0, "a reasoning round must not be re-sent"


def test_a_round_that_called_a_tool_is_not_re_sent():
    client, seen = _stub_client(empty_stream=False, tool_call=True)
    _first_text(client)
    assert seen["plain"] == 0, "a tool-calling round must not be re-sent"
