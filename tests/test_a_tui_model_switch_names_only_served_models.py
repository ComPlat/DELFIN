"""/model in the window: what it offers, what it refuses, what it leaves alone.

A bare /model opens the window's picker over session.choices; /model <name>
goes to the session as typed and through the REPL's own command
(repl_commands._model -> client.switch_model), the one path the dashboard
and the terminal use too. These tests pin three things about that path:

* a name the provider's /v1/models does not list is refused with the
  reason and the nearest names, and the model stays as it was -- an
  unknown name used to be taken and failed on the next turn's request;
* when the listing cannot be read, a KIT session's picker offers the
  curated roster (model_capabilities._KIT_PREFERENCE) instead of the
  current model alone, and a typed name is switched to unchecked, saying so;
* a switch changes the model and nothing else: the approval posture is
  the one it was, and the window hears only the model change.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from delfin.agent import model_capabilities as mc
from delfin.agent import repl_commands as rc
from delfin.agent.service import choices
from tui_worker_client import WorkerClient

BASE = "https://kit.example/v1"
SERVED = [{"id": "kit.deepseek-v4-flash"}, {"id": "kit.glm-5.3"}, {"id": "azure.gpt-5.5"},
          {"id": "kit.qwen3-embedding-8b"}]


class _Client:
    def __init__(self, model: str = "kit.deepseek-v4-flash"):
        self.model = model
        self._base_url = BASE
        self._api_key = "k"
        self._provider = "kit"
        self.switched: list[str] = []

    def switch_model(self, name: str) -> None:
        self.switched.append(name)
        self.model = name


@pytest.fixture
def listing(monkeypatch):
    """The endpoint's /v1/models, swappable; nothing cached from before."""
    answer = {"entries": SERVED}

    def fetch(base, key):
        entries = answer["entries"]
        if entries is not None:
            mc._remember_served(base, entries)
        return entries

    monkeypatch.setattr(mc, "_fetch_openai_models", fetch)
    monkeypatch.setattr(mc, "_SERVED_CACHE", {})
    monkeypatch.setattr(choices, "_models_cache", {})
    return answer


def _ctx(client):
    return SimpleNamespace(engine=SimpleNamespace(client=client, provider="kit"))


def test_an_unknown_name_is_refused_with_the_nearest_ones(listing):
    client = _Client()
    out = rc._model(_ctx(client), "kit.glm-5.4").output
    assert "unknown model 'kit.glm-5.4'" in out
    assert "kit.glm-5.3" in out                         # the nearest name
    assert client.switched == [] and client.model == "kit.deepseek-v4-flash"
    # Not a chat model: listed by the endpoint, refused all the same.
    out = rc._model(_ctx(client), "kit.qwen3-embedding-8b").output
    assert "unknown model" in out and client.switched == []


def test_a_served_name_is_switched_to(listing):
    client = _Client()
    assert rc._model(_ctx(client), "kit.glm-5.3").output == "model → kit.glm-5.3"
    assert client.switched == ["kit.glm-5.3"]


def test_an_unreadable_listing_switches_unchecked_and_says_so(listing):
    listing["entries"] = None
    client = _Client()
    out = rc._model(_ctx(client), "kit.something-new").output
    assert out.startswith("model → kit.something-new")
    assert "not checked" in out
    assert client.switched == ["kit.something-new"]


def test_a_kit_picker_offers_the_roster_when_the_listing_cannot_be_read(listing):
    listing["entries"] = None
    engine = SimpleNamespace(client=_Client("kit.glm-5.3"), provider="kit")
    got = choices.choices(engine, "model")
    values = [c["value"] for c in got["choices"]]
    assert got["current"] == "kit.glm-5.3"
    assert set(mc._KIT_PREFERENCE) <= set(values)
    assert values.count("kit.glm-5.3") == 1
    assert all("could not be read" in c["detail"] for c in got["choices"] if c["value"] != got["current"])


def test_a_picker_lists_what_the_endpoint_serves_and_marks_the_current(listing):
    engine = SimpleNamespace(client=_Client(), provider="kit")
    got = choices.choices(engine, "model")
    assert got["current"] == "kit.deepseek-v4-flash"
    assert [c["value"] for c in got["choices"]] == ["kit.deepseek-v4-flash", "kit.glm-5.3", "azure.gpt-5.5"]


def test_a_switch_through_the_session_changes_the_model_and_nothing_else(tmp_path):
    with WorkerClient(tmp_path) as w:
        approval = w.init_result["session"]["approval"]
        mark = len(w.events)
        w.request("session.prompt", {"text": "/model fake-model-2"})
        updated, _ = w.wait_for(lambda e: e.get("type") == "session.updated" and "model" in e, after=mark)
        assert updated["model"] == "fake-model-2"
        assert "approval" not in updated
        status = w.request("session.status")
        assert status["approval"] == approval
        notice, _ = w.wait_for(lambda e: e.get("type") == "notice" and "model" in e.get("text", ""), after=mark)
        assert notice["text"] == "model → fake-model-2"
