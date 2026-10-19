"""KIT runs on its best model, and a model KIT does not serve is explained.

2026-10-19: a request to kit.glm-5.3 failed with "400 - Model not found" and
nothing else, and the dashboard's cheap KIT default was azure.gpt-5-nano.
The user's choice since then: GLM-5.3 is the KIT default, DeepSeek V4 Flash
the next and the cheap tier. When KIT does not serve the model, the session
continues once per turn on the best KIT model it serves, with a notice.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent import mcp_client as M
from delfin.agent import model_capabilities as mc
from delfin.agent import model_routing as MR


class _NotFound(Exception):
    status_code = 400

    def __init__(self):
        super().__init__("Error code: 400 - {'detail': 'Model not found'}")


def test_the_ranking_puts_glm_first_and_deepseek_next():
    assert mc._KIT_PREFERENCE[:2] == ("kit.glm-5.3", "kit.deepseek-v4-flash")
    assert MR._PROVIDER_TIERS["kit"] == {"strong": "kit.glm-5.3",
                                         "cheap": "kit.deepseek-v4-flash"}


def test_model_not_found_is_recognised_and_other_400s_are_not():
    assert A._is_model_not_found_error(_NotFound())
    assert not A._is_model_not_found_error(ValueError("context length exceeded"))
    other = _NotFound()
    other.status_code = 500
    assert not A._is_model_not_found_error(other)


@pytest.mark.parametrize("served,current,expected", [
    ({"kit.deepseek-v4-flash", "azure.gpt-5.4"}, "kit.glm-5.3",
     "kit.deepseek-v4-flash"),
    ({"kit.glm-5.3", "kit.deepseek-v4-flash"}, "azure.gpt-5-nano",
     "kit.glm-5.3"),
    (set(), "kit.glm-5.3", "kit.deepseek-v4-flash"),     # listing unreadable
    ({"azure.gpt-5.4"}, "kit.glm-5.3", ""),              # nothing ranked
])
def test_the_fallback_is_the_best_served(monkeypatch, served, current,
                                         expected):
    monkeypatch.setattr(mc, "kit_models_served", lambda *a, **k: served)
    assert A._kit_fallback_model(current, "kit") == expected


def test_another_provider_gets_no_fallback():
    assert A._kit_fallback_model("gpt-5.4", "openai") == ""


# -- through the tool loop --------------------------------------------------

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
        workspace=tmp_path, mode="default", task_session_id="sess-kit"))
    return cl


def test_a_missing_model_continues_on_the_best_served_one(client):
    asked = []

    def _create(**kw):
        asked.append(kw.get("model"))
        if kw.get("model") == "kit.glm-5.3":
            raise _NotFound()
        return _Stream([_Chunk([_Choice(_Delta("Hallo."), finish="stop")],
                               usage=_Usage())])

    client.client.chat.completions.create = _create
    events = list(client.stream_message(
        "sys", [{"role": "user", "content": "hi"}], max_tokens=50))
    notices = "".join(getattr(e, "text", "") or "" for e in events
                      if getattr(e, "type", "") == "notice")
    assert asked == ["kit.glm-5.3", "kit.deepseek-v4-flash"]
    assert "does not serve `kit.glm-5.3`" in notices
    assert "`kit.deepseek-v4-flash`" in notices
    assert client.model == "kit.deepseek-v4-flash"
