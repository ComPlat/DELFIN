"""Per-model tuning for the two KIT flagships, each knob switchable.

Field data (shared bug archive, 33 sessions, 2026-09-16 .. 10-14) and the
benchmark store (950 DeepSeek attempts, 273 GLM attempts):

* kit.deepseek-v4-flash decodes at about 6 s per 100 output tokens
  (request p50, archive) and spends 81% of its turn time generating. Every
  one of its 38 archived tool-call messages carried prose beside the call,
  p50 593 characters ("Let me read ... and then look at ..."), against
  259 of 712 on kit.glm-5.3 at p50 160.

Each knob is a ModelProfile field, so ``agent.model_overrides`` switches
it per model without editing Python, and a benchmark arm that does so is
stamped with the override in its profile name.
"""

from __future__ import annotations

import inspect
import json

import pytest

from delfin.agent import model_profiles as mp
from delfin.agent.prompt_loader import PromptLoader


@pytest.fixture
def settings(tmp_path, monkeypatch):
    """A settings file this test owns; the user's is never read."""
    path = tmp_path / "settings.json"
    path.write_text(json.dumps({"agent": {}}))

    def _write(overrides):
        path.write_text(json.dumps({"agent": {"model_overrides": overrides}}))
        mp._OVERRIDES_CACHE = None

    monkeypatch.setattr("delfin.user_settings.get_settings_path",
                        lambda *a, **k: path)
    mp._OVERRIDES_CACHE = None
    yield _write
    mp._OVERRIDES_CACHE = None


_DEEPSEEK = "kit.deepseek-v4-flash"
_GLM = "kit.glm-5.3-flash"


def _solo_sections(model: str):
    return PromptLoader().compose_sections(
        role_id="solo_agent", mode_id="solo", task_text="hallo",
        session_key=f"s-{model}", model=model)


def _solo_prompt(model: str) -> str:
    return "\n\n".join(s.content for s in _solo_sections(model))


# --- terse tool rounds ------------------------------------------------------

def test_deepseek_is_told_to_keep_tool_rounds_terse(settings):
    assert mp.get_profile(_DEEPSEEK).terse_tool_rounds is True
    note = mp.MODEL_NOTES["terse_tool_rounds"]
    assert note in _solo_prompt(_DEEPSEEK)


def test_glm_is_not(settings):
    """GLM-flash decodes ~7x faster (0.9 s per 100 tokens) and narrates in
    36% of its tool rounds; the note would cost it prompt and buy little."""
    assert mp.get_profile(_GLM).terse_tool_rounds is False
    assert mp.MODEL_NOTES["terse_tool_rounds"] not in _solo_prompt(_GLM)


def test_the_note_can_be_switched_off_per_model(settings):
    settings({"kit.deepseek": {"terse_tool_rounds": False}})
    assert mp.get_profile(_DEEPSEEK).terse_tool_rounds is False
    assert mp.MODEL_NOTES["terse_tool_rounds"] not in _solo_prompt(_DEEPSEEK)


def test_the_note_sits_in_the_stable_head(settings):
    """It is fixed for the model, so it belongs before every per-turn
    section: a note after the volatile material would be re-sent cold
    whenever the git status line moved."""
    sections = _solo_sections(_DEEPSEEK)
    names = [s.name for s in sections]
    assert "model_notes" in names
    note = sections[names.index("model_notes")]
    assert note.layer == PromptLoader.LAYER_STABLE
    layers = [s.layer for s in sections]
    assert layers == sorted(layers)


# --- concise final answers --------------------------------------------------

def test_deepseek_is_told_to_answer_concisely(settings):
    """Benchmark, 950 DeepSeek attempts: within one task, answers longer
    than that task's median passed 0.79 (n=406) and shorter ones 0.75
    (n=399) -- no gain -- at a duration p50 of 55 s against 40 s."""
    assert mp.get_profile(_DEEPSEEK).concise_answers is True
    assert mp.MODEL_NOTES["concise_answers"] in _solo_prompt(_DEEPSEEK)


def test_glm_answers_are_left_as_they_are(settings):
    """GLM's passing answers are already short (p50 710 chars against
    DeepSeek's 2146), and its within-task longer answers passed less
    often, not more -- nothing for the note to buy."""
    assert mp.get_profile(_GLM).concise_answers is False
    assert mp.MODEL_NOTES["concise_answers"] not in _solo_prompt(_GLM)


def test_session_affinity_is_a_glm_knob_and_off_by_default(settings):
    # Off by default: the A/B showed no cache or latency gain (see the profile).
    assert mp.get_profile(_GLM).session_affinity is False
    assert mp.get_profile(_DEEPSEEK).session_affinity is False
    settings({"kit.glm": {"session_affinity": True}})
    assert mp.get_profile(_GLM).session_affinity is True
    assert mp.get_profile("kit.glm-5.3").session_affinity is True


# --- session affinity: the request carries a stable opaque `user` -----------
#
# Archive: 408 of 2066 kit.glm-5.3 requests went cold (cached/input < 0.2)
# with the system prompt and the turn prefix unchanged and under 300 s
# since the previous request -- ~740 min of time to first token, the
# largest single cost in the data. Same bytes, short gap, still cold fits
# replica load-balancing without prefix affinity. The OpenAI `user` field
# is what an OpenAI-compatible router can key a sticky route on; the KIT
# endpoint accepts it (probed 2026-10-14 on kit.glm-5.3-flash and
# kit.deepseek-v4-flash, HTTP 200). Whether it routes on it is what the
# A/B measures.

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
    prompt_tokens_details = None


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


class _FakeRegistry:
    def discover_all(self):
        return []

    def discover_resources(self):
        return []

    def discover_prompts(self):
        return []


@pytest.fixture
def kit_client(monkeypatch, tmp_path, settings):
    from delfin.agent import mcp_client as M
    from delfin.agent import model_capabilities as mc
    import delfin.agent.api_client as A

    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "sk-test-not-a-key")
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _FakeRegistry())
    monkeypatch.setattr(A.time, "sleep", lambda *a, **k: None)

    def _make(model: str, session: str = "trace-1234"):
        caps = mc.ModelCapabilities(model=model, provider="kit",
                                    context_window=131_072,
                                    supports_tools=True)
        monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
        client = A.create_client(backend="api", provider="kit", model=model,
                                 cwd=str(tmp_path))
        client.metrics_session = session
        sent: list[dict] = []

        def _create(**kwargs):
            sent.append(kwargs)
            return _Stream([_Chunk([_Choice(_Delta("ok"), finish="stop")],
                                   usage=_Usage())])

        client.client.chat.completions.create = _create
        return client, sent

    return _make


def _turn(client):
    return [ev for ev in client.stream_message(
        "sys", [{"role": "user", "content": "hallo"}], max_tokens=100)]


def test_glm_requests_carry_a_stable_opaque_user(kit_client, settings):
    settings({"kit.glm": {"session_affinity": True}})
    client, sent = kit_client(_GLM, session="trace-1234")
    _turn(client)
    _turn(client)
    users = [k.get("user") for k in sent]
    assert len(users) == 2 and users[0] and users[0] == users[1]
    assert users[0].startswith("delfin-")
    # Opaque: the session key itself never leaves the machine.
    assert "trace-1234" not in users[0] and "1234" not in users[0]


def test_another_session_gets_another_route_key(kit_client, settings):
    settings({"kit.glm": {"session_affinity": True}})
    a, sent_a = kit_client(_GLM, session="trace-a")
    _turn(a)
    b, sent_b = kit_client(_GLM, session="trace-b")
    _turn(b)
    assert sent_a[0]["user"] != sent_b[0]["user"]


def test_the_route_key_is_not_computed_from_the_session(monkeypatch, settings):
    settings({"kit.glm": {"session_affinity": True}})
    # Drawn at random per key, not derived: the same key in a fresh
    # registry gets another id, so no digest of session data is sent.
    first = mp.session_route_key(_GLM, "trace-same")
    monkeypatch.setattr(mp, "_ROUTE_IDS", {})
    second = mp.session_route_key(_GLM, "trace-same")
    assert first.startswith("delfin-") and second.startswith("delfin-")
    assert first != second
    assert "hashlib" not in inspect.getsource(mp.session_route_key)


def test_deepseek_requests_are_unchanged(kit_client):
    client, sent = kit_client(_DEEPSEEK)
    _turn(client)
    assert "user" not in sent[0]


def test_glm_sends_no_user_by_default(kit_client):
    client, sent = kit_client(_GLM)
    _turn(client)
    assert "user" not in sent[0]


def test_no_session_key_sends_no_user(kit_client):
    client, sent = kit_client(_GLM, session="")
    _turn(client)
    assert "user" not in sent[0]


def test_each_note_switches_alone(settings):
    settings({"kit.deepseek": {"concise_answers": False}})
    text = _solo_prompt(_DEEPSEEK)
    assert mp.MODEL_NOTES["concise_answers"] not in text
    assert mp.MODEL_NOTES["terse_tool_rounds"] in text
