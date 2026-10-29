"""A model KIT stopped serving must switch to another, a gateway-wide
exhaustion must keep retrying the same model, and the switch is bounded.

2026-10-29, measured for about an hour: KIT answered every request to one
model with HTTP 400 and the text "Model 'glm-5.3-flash' is temporarily
unavailable. Please try again later." DELFIN classified that text as gateway
exhaustion (_INFRA_EXHAUSTION_MARKERS, "temporarily unavailable") and retried
the same model; the switch to another model (_kit_fallback_model) only fired
for "Model not found", once per turn. Ten agents stood still for an hour.

A MODEL being unavailable is a property of the model, not of the gateway:
switching to the next family KIT serves fixes it. A gateway that is out of
connections ("too many connections") is not fixed by switching models, and
must keep retrying the same model. No live model: create() is faked.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent import mcp_client as M
from delfin.agent import model_capabilities as mc


_MODEL_UNAVAILABLE_BODY = (
    "Error code: 400 - {'error': 'Model \\'glm-5.3-flash\\' is temporarily "
    "unavailable. Please try again later.'}"
)
_GATEWAY_TOO_MANY_CONNECTIONS = "Error code: 400 - too many connections for role"


class _Status400(Exception):
    def __init__(self, msg):
        super().__init__(msg)
        self.status_code = 400


def _model_unavailable_error():
    return _Status400("Model 'glm-5.3-flash' is temporarily unavailable. "
                      "Please try again later.")


# --- plain-classifier checks ------------------------------------------------


def test_a_model_fulling_out_is_not_that():
    assert A._is_model_unavailable_error(_model_unavailable_error())


def test_unavailable_names_another_model_too():
    assert A._is_model_unavailable_error(
        _Status400("Model 'gpt-oss-120b' is not available"))
    assert A._is_model_unavailable_error(
        _Status400("The model qwen3.5-397b-A17b is currently unavailable"))


def test_gateway_exhaustion_is_not_model_unavailable():
    assert not A._is_model_unavailable_error(
        _Status400(_GATEWAY_TOO_MANY_CONNECTIONS))
    assert not A._is_model_unavailable_error(
        _Status400(_GATEWAY_TOO_MANY_CONNECTIONS))
    assert not A._is_model_unavailable_error(
        _Status400("Error code: 400 - service temporarily unavailable"))
    assert not A._is_model_unavailable_error(
        _Status400("Error code: 400 - connection pool exhausted"))


def test_a_model_unavailable_400_is_not_gateway_exhaustion():
    # It must NOT path into the transient same-model retry.
    assert not A._is_transient_api_error(_model_unavailable_error())


def test_gateway_exhaustion_is_still_transient():
    assert A._is_transient_api_error(
        _Status400(_GATEWAY_TOO_MANY_CONNECTIONS))
    assert A._is_transient_api_error(
        _Status400("Error code: 400 - service temporarily unavailable"))


# --- through the tool loop --------------------------------------------------

class _Delta:
    def __init__(self, content=None):
        self.content = content
        self.tool_calls = None


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta = delta
        self.finish_reason = finish


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices = choices
        self.usage = usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


class _EmptyReg:
    def discover_all(self): return []
    def discover_resources(self): return []
    def discover_prompts(self): return []


@pytest.fixture
def client(monkeypatch, tmp_path):
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    monkeypatch.setattr(mc, "kit_models_served",
                        lambda *a, **k: {"kit.deepseek-v4-flash"})
    cl = A.create_client(backend="api", provider="ollama",
                         model="kit.glm-5.3", cwd=str(tmp_path))
    cl._provider = "kit"
    cl.set_permissions(A.KitToolPermissions(
        workspace=tmp_path, mode="default", task_session_id="sess-switch"))
    return cl


def _notices(events):
    return "".join(getattr(e, "text", "") or "" for e in events
                   if getattr(e, "type", "") == "notice")


def test_a_model_kits_stops_serving_switches_to_the_best_other(client):
    asked = []

    def _create(**kw):
        asked.append(kw.get("model"))
        if kw.get("model") == "kit.glm-5.3":
            raise _Status400("Model 'glm-5.3-flash' is temporarily "
                             "unavailable. Please try again later.")
        return _Stream([_Chunk([_Choice(_Delta("Hallo."), finish="stop")],
                               usage=_Usage())])

    client.client.chat.completions.create = _create
    events = list(client.stream_message(
        "sys", [{"role": "user", "content": "hi"}], max_tokens=50))
    text = "".join(getattr(e, "text", "") or "" for e in events
                   if getattr(e, "type", "") == "text_delta")
    assert text == "Hallo."
    assert asked[0] == "kit.glm-5.3"
    assert asked[-1] == "kit.deepseek-v4-flash"
    assert client.model == "kit.deepseek-v4-flash"
    assert "is temporarily unavailable" in _notices(events)
    assert "continues on `kit.deepseek-v4-flash`" in _notices(events)
    assert getattr(client, "_turn_model_switches", []) == [{
        "before": "kit.glm-5.3", "after": "kit.deepseek-v4-flash",
        "reason": "model unavailable"}]


def test_gateway_exhaustion_retries_the_same_model(client):
    asked = []

    def _create(**kw):
        asked.append(kw.get("model"))
        raise RuntimeError("Error code: 400 - too many connections for role")

    client.client.chat.completions.create = _create
    with pytest.raises(RuntimeError):
        list(client.stream_message(
            "sys", [{"role": "user", "content": "hi"}], max_tokens=50))
    assert all(m == "kit.glm-5.3" for m in asked)
    assert len(asked) >= 2           # retried the SAME model, did not switch
    assert getattr(client, "_turn_model_switches", []) == []
    assert client.model == "kit.glm-5.3"


def test_the_switch_is_bounded(client):
    # Every ranked KIT model answers unavailable: the round must stop after
    # a bounded number of switches, never loop through roster forever.
    asked = []

    def _create(**kw):
        asked.append(kw.get("model"))
        raise _Status400(f"Model '{kw.get('model')}' is temporarily "
                         "unavailable. Please try again later.")

    client.client.chat.completions.create = _create
    with pytest.raises(Exception):
        list(client.stream_message(
            "sys", [{"role": "user", "content": "hi"}], max_tokens=50))
    # At most 3 switches plus the startup request; never an unbounded loop.
    assert len(asked) <= 1 + A._MODEL_FALLBACK_MAX + 1
    switches = getattr(client, "_turn_model_switches", [])
    assert len(switches) <= A._MODEL_FALLBACK_MAX


def test_two_models_down_the_turn_ends_on_the_third(client, monkeypatch):
    # GLM and DeepSeek both report themselves unavailable, Mistral answers.
    # The turn leaves GLM, then DeepSeek, and lands on Mistral -- it never
    # hands the rank slot back to a model it already left (GLM -> DeepSeek
    # -> GLM until the budget is spent).
    monkeypatch.setattr(mc, "kit_models_served", lambda *a, **k: {
        "kit.glm-5.3-flash", "kit.deepseek-v4-flash",
        "kit.mistral-small-4-119b-a8b"})
    client.model = "kit.glm-5.3-flash"
    asked = []

    def _create(**kw):
        asked.append(kw.get("model"))
        if kw.get("model") in ("kit.glm-5.3-flash", "kit.deepseek-v4-flash"):
            raise _Status400(f"Model '{kw.get('model')}' is temporarily "
                             "unavailable. Please try again later.")
        return _Stream([_Chunk([_Choice(_Delta("Fertig."), finish="stop")],
                               usage=_Usage())])

    client.client.chat.completions.create = _create
    events = list(client.stream_message(
        "sys", [{"role": "user", "content": "hi"}], max_tokens=50))
    text = "".join(getattr(e, "text", "") or "" for e in events
                   if getattr(e, "type", "") == "text_delta")
    assert text == "Fertig."
    assert client.model == "kit.mistral-small-4-119b-a8b"
    switches = getattr(client, "_turn_model_switches", [])
    assert [(s["before"], s["after"]) for s in switches] == [
        ("kit.glm-5.3-flash", "kit.deepseek-v4-flash"),
        ("kit.deepseek-v4-flash", "kit.mistral-small-4-119b-a8b")]
    first_left = asked.index("kit.deepseek-v4-flash")
    assert "kit.glm-5.3-flash" not in asked[first_left:]
