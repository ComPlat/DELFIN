"""Ask the endpoint one word before the user does.

A model selected in the dashboard was a name in a dropdown until the
first turn, and the first turn on a queued endpoint was a spinner for
ten minutes. Measured on the KIT deployment on 2026-09-11: kit.glm-5.3
took 459 s to its first token for a one-word question with a one-line
system prompt, kit.deepseek-v4-flash 0.6 s. Nothing in the dashboard
said which of the two the user had just picked.

The probe sends one word and stops at the first chunk. Its cost is one
prompt of ~20 tokens and one output token; its answer is a number the
user can act on before they wait.
"""

from __future__ import annotations

import time
from dataclasses import dataclass


@dataclass(frozen=True)
class ProbeResult:
    model: str
    seconds: float
    answered: bool
    error: str = ""


def probe_first_token(sdk_client, model: str, *, provider: str = "",
                      timeout_s: float = 60.0,
                      reasoning_effort: str = "") -> ProbeResult:
    """One word to ``model`` through an OpenAI-compatible ``sdk_client``;
    the time until the first streamed chunk, or the reason there was none.

    ``reasoning_effort`` is sent when given -- a reasoning model asked
    with nothing set spends its whole budget thinking before the first
    visible token, which is the wait being measured, not the queue.
    """
    kwargs = {
        "model": model,
        "messages": [{"role": "user", "content": "hi"}],
        "max_tokens": 8,
        "stream": True,
        "timeout": float(timeout_s),
    }
    if reasoning_effort:
        kwargs["reasoning_effort"] = reasoning_effort
    t0 = time.monotonic()
    try:
        stream = sdk_client.chat.completions.create(**kwargs)
        try:
            for _chunk in stream:
                return ProbeResult(model, time.monotonic() - t0, True)
        finally:
            try:
                stream.close()
            except Exception:
                pass
        return ProbeResult(model, time.monotonic() - t0, False,
                           "the stream ended without a chunk")
    except Exception as exc:  # timeout, 5xx, connection -- all one answer
        return ProbeResult(model, time.monotonic() - t0, False,
                           f"{type(exc).__name__}")


def describe(result: ProbeResult, *, slow_cold_start_s: float = 0.0,
             first_token_budget_s: float = 600.0) -> str:
    """What the dashboard says about the probe, in one line."""
    name = result.model or "the model"
    if result.answered:
        secs = result.seconds
        if secs < 2.0:
            return f"{name} answered a one-word probe in {secs:.1f} s."
        text = f"{name} answered a one-word probe in {secs:.0f} s."
        if secs >= 30.0:
            text += (" The endpoint is queueing; a real turn waits at least "
                     "that long before its first token.")
        return text
    waited = result.seconds
    text = (f"{name}: no answer to a one-word probe after {waited:.0f} s"
            + (f" ({result.error})" if result.error else "") + ".")
    text += (" The endpoint is queueing. A turn is given "
             f"{first_token_budget_s / 60:.0f} min for its first token; "
             "another model answers sooner, or wait and watch the count.")
    return text
