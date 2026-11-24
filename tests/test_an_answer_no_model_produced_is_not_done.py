"""Twenty-two agents ended "done" on an answer no model had produced.

MEASURED on 2026-11-23, wave 40: every agent on Blablador ended within
one to four minutes, each with the same result line:

    "text": "There's an issue with the selected model (alias-glm-5.3-
             flash). It may not exist or you may not have access to it.
             Run --model to pick a different model."
    "input_tokens": 0, "output_tokens": 0, "stop_reason": "done"

The key was valid (an authenticated request listed 76 models, ours among
them), the endpoint was up, and the same call answered "OK" an hour and a
half later. The models were not being served for a while, and the
gateway said so as a TEXT reply rather than as an HTTP error. DELFIN
passed that text on as the model's answer and called the turn done --
so a provider outage read as twenty-two finished agents, and the first
diagnosis went looking for a bad key and bad model names.

THE SIGNAL IS THE TOKEN COUNT, NOT THE WORDING. A model that answers has
read the prompt, so it reports input tokens; the successful Blablador
runs of the same morning reported 10,602,162 of them. Zero input and zero
output tokens on a turn that produced text is an answer no model ran
for. Matching the gateway's sentence would break the day it is
reworded; the count does not.

ONLY WHERE IT IS MEASURED. Every OpenAI-compatible request asks for
usage (``stream_options.include_usage``), but whether a provider sends it
was measured for KIT and Blablador only. An engine of another provider,
or a test double with no provider, keeps today's behaviour exactly.
"""
from __future__ import annotations

from delfin.agent import cli as agent_cli

GATEWAY = ("There's an issue with the selected model (alias-glm-5.3-flash). "
           "It may not exist or you may not have access to it. Run --model "
           "to pick a different model.")


class _Engine:
    session_id = "sess-no-model"
    run_budget_turns = 0

    def __init__(self, provider, text=GATEWAY, *, tokens=(0, 0), tools=0):
        self.provider = provider
        self.text = text
        self._tokens = tokens
        self.tools = tools
        self.token_usage = {"input": 0, "output": 0}

    def stream_response(self, **kw):
        on_tool, on_token = kw.get("on_tool_use"), kw.get("on_token")
        for _ in range(self.tools):
            if on_tool:
                on_tool("bash", '{"command": "ls"}')
        if on_token:
            on_token(self.text)
        self.token_usage = {"input": self._tokens[0],
                            "output": self._tokens[1]}
        return ""


def _run(engine):
    return agent_cli._run_once(engine, "Reply with exactly: OK")


# --- the measured case ---------------------------------------------------

def test_a_zero_token_answer_from_blablador_is_an_error_not_done():
    out = _run(_Engine("blablador"))
    assert out["stop_reason"] != "done", out
    assert out["error"], "the outage must be named, not passed as an answer"
    assert "without running the model" in out["error"]
    # the gateway's own words travel with it, so the reader sees the cause
    assert "issue with the selected model" in out["error"]


def test_the_same_holds_for_kit():
    out = _run(_Engine("kit"))
    assert out["stop_reason"] != "done" and out["error"]


# --- what it must not touch ---------------------------------------------

def test_a_real_answer_is_done():
    out = _run(_Engine("blablador", "OK", tokens=(412, 3)))
    assert out["stop_reason"] == "done" and not out["error"], out


def test_a_provider_whose_usage_is_not_measured_is_left_alone():
    """Ollama is asked for usage like the others; whether it sends it was
    not measured here, so a zero there proves nothing."""
    out = _run(_Engine("ollama", "OK"))
    assert out["stop_reason"] == "done" and not out["error"], out


def test_a_test_double_without_a_provider_is_left_alone():
    class _Stub:
        token_usage = {"input": 0, "output": 0}

        def stream_response(self, **kwargs):
            return "x"

    out = agent_cli._run_once(_Stub(), "hi")
    assert not out["error"], out


def test_a_turn_that_called_tools_is_not_judged_by_this():
    """Tool rounds make the count a different question; only a bare text
    answer with nothing spent is the shape measured."""
    out = _run(_Engine("blablador", tools=1))
    assert "without running the model" not in str(out["error"] or "")
