"""What DELFIN offers must be what the endpoint serves.

On 2026-09-07 the KIT roster dropped qwen3.5-397b, gpt-oss-120b,
gemma4-31b and minimax in one go. DELFIN went on naming
``kit.qwen3.5-397b-A17b`` as the best KIT model in the warning it prints
when a user picks a weak one — advice to switch to a model that answers
with a 404. The dropdown's offline fallback offered four of them.

Nothing here reaches the network: the listing is stubbed.
"""

import ast
from pathlib import Path

import pytest

from delfin.agent import model_capabilities as mc
from delfin.agent.model_profiles import get_profile

_REPO = Path(__file__).resolve().parents[1]


@pytest.fixture
def served(monkeypatch):
    """Stub the endpoint listing; returns the mutable set it will answer."""
    ids = set(mc._KIT_PREFERENCE)

    def _fake(base_url="", api_key="", *, fetch=False):
        return set(ids)

    monkeypatch.setattr(mc, "kit_models_served", _fake)
    return ids


# --- the recommendation names something that exists --------------------

def test_the_best_model_is_the_first_one_actually_served(served):
    assert mc.kit_best_model() == mc._KIT_PREFERENCE[0]
    served.discard(mc._KIT_PREFERENCE[0])
    assert mc.kit_best_model() == mc._KIT_PREFERENCE[1]


def test_an_unreachable_listing_falls_back_to_the_constant(served):
    served.clear()
    assert mc.kit_best_model() == mc.KIT_BEST_MODEL


def test_a_retired_model_is_named_as_retired_not_as_weak(served):
    """It was the recommended model until the roster changed. Saying it is
    weak would be both wrong and unhelpful."""
    caps = mc.resolve("kit", "kit.qwen3.5-397b-A17b", "", allow_live=False)
    msg = mc.kit_recommendation("kit.qwen3.5-397b-A17b", caps)
    assert "no longer served" in msg
    assert "weak" not in msg
    assert mc.kit_best_model() in msg


def test_a_retired_model_is_not_waved_through(served):
    """Its static entry still says 128k window and native tools -- true
    while it was served -- so every quality check passes it. Being gone
    has to be checked on its own or the user is told nothing."""
    for name in sorted(mc._KIT_RETIRED):
        caps = mc.resolve("kit", name, "", allow_live=False)
        assert caps.supports_tools and caps.context_window >= mc._KIT_MIN_WINDOW
        assert mc.kit_recommendation(name, caps), name


def test_a_served_model_draws_no_warning(served):
    for name in mc._KIT_PREFERENCE:
        caps = mc.resolve("kit", name, "", allow_live=False)
        assert mc.kit_recommendation(name, caps) == "", name


def test_retired_and_recommended_are_disjoint():
    assert not (mc._KIT_RETIRED & mc._KIT_RECOMMENDED)


# --- a recommended model is a tuned model ------------------------------

@pytest.mark.parametrize("name", mc._KIT_PREFERENCE)
def test_every_recommended_model_is_known_to_both_registries(name):
    caps = mc.resolve("kit", name, "", allow_live=False)
    assert caps.source == "static", (name, caps.source)
    assert caps.context_window >= mc._KIT_MIN_WINDOW
    profile = get_profile(name)
    assert profile.notes
    assert profile is not None


def test_glm_is_marked_a_reasoning_family():
    """It spends the completion budget on hidden reasoning before any
    content: max_tokens=32 returned 32 reasoning tokens and an empty
    message on every probe. api_client floors the budget for a reasoning
    model; without the flag the floor never applies and the turn is
    reported empty."""
    caps = mc.resolve("kit", "kit.glm-5.3", "", allow_live=False)
    assert caps.is_reasoning


def test_glm_is_given_longer_than_its_measured_cold_start():
    """Cold prompt heads measured at 199 / 266 / 268 s on 2026-09-07. A
    120 s stale-kill would cut off a turn that was about to answer, and
    the retry pays the same cold prefill again."""
    assert get_profile("kit.glm-5.3").stale_kill_after_s >= 300.0


def test_a_point_revision_keeps_the_profile():
    assert get_profile("kit.glm-5.9").stale_kill_after_s >= 300.0
    assert get_profile("kit.deepseek-v5-flash").notes


def test_the_qwen_tuning_survived_the_roster_change():
    """Retired from the endpoint, not deleted from the code: the numbers
    were measured, and a restored session may still name the model."""
    p = get_profile("kit.qwen3.5-397b-A17b")
    assert p.max_tool_rounds == 20
    assert p.stale_kill_after_s == 90.0
    assert "Qwen" in p.notes
    caps = mc.resolve("kit", "kit.qwen3.5-397b-A17b", "", allow_live=False)
    assert caps.context_window == 128_000


# --- the offline dropdown ----------------------------------------------

def _fallback_kit_models() -> list[str]:
    """The KIT list from the dashboard's offline fallback.

    It lives inside ``create_tab``'s closure, so it is read from the source
    rather than imported. Any literal named _PROVIDER_MODELS_FALLBACK does.
    """
    tree = ast.parse((_REPO / "delfin" / "dashboard" / "tab_agent.py").read_text())
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        names = [t.id for t in node.targets if isinstance(t, ast.Name)]
        if "_PROVIDER_MODELS_FALLBACK" not in names:
            continue
        table = ast.literal_eval(node.value)
        return [mid for _label, mid in table["kit"]]
    raise AssertionError("_PROVIDER_MODELS_FALLBACK not found")


def test_the_offline_dropdown_offers_no_retired_model():
    offered = _fallback_kit_models()
    assert offered
    assert not (set(offered) & mc._KIT_RETIRED), offered


def test_the_offline_dropdown_offers_the_kit_hosted_models():
    offered = set(_fallback_kit_models())
    assert set(mc._KIT_PREFERENCE) <= offered, sorted(offered)


def test_the_settings_suggestions_offer_no_retired_model():
    src = (_REPO / "delfin" / "dashboard" / "tab_settings.py").read_text()
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        names = [t.id for t in node.targets if isinstance(t, ast.Name)]
        if "_JOBMON_MODEL_SUGGESTIONS" not in names:
            continue
        table = ast.literal_eval(node.value)
        assert not (set(table["kit"]) & mc._KIT_RETIRED), table["kit"]
        return
    raise AssertionError("_JOBMON_MODEL_SUGGESTIONS not found")


def test_the_help_examples_name_a_model_that_exists():
    """A --model example is advice. It must not name a retired model."""
    for rel in ("delfin/agent/cli.py",
                "delfin/agent/AGENT_USER_GUIDE.txt",
                "delfin/agent/benchmark_runner.py"):
        text = (_REPO / rel).read_text()
        for retired in mc._KIT_RETIRED:
            assert retired not in text, (rel, retired)


def test_the_probe_budget_clears_the_measured_response_time():
    """Four consecutive /v1/models fetches on 2026-09-07: 5.02 / 4.77 /
    4.77 / 7.66 s. A timeout inside that band loses the race often, and
    losing it means a 512k model is driven as a 32k one."""
    assert mc._TIMEOUT_MODELS >= 15.0


# --- one fetch answers both questions ----------------------------------

def test_the_capability_probe_records_the_roster(monkeypatch):
    """The roster question and the window question are answered by the
    same HTTP response, and the window question runs on every preflight.
    Asking the wire twice for one answer is what the recording avoids."""
    import json
    from io import BytesIO

    payload = {"data": [
        {"id": "kit.glm-5.3", "max_model_len": 524_288},
        {"id": "kit.whisper-large-v3"},
    ]}

    class _Resp(BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    calls = []

    def _urlopen(req, timeout=0):
        calls.append(getattr(req, "full_url", str(req)))
        return _Resp(json.dumps(payload).encode())

    monkeypatch.setattr(mc.urllib.request, "urlopen", _urlopen)
    mc._SERVED_CACHE.clear()
    mc.clear_cache()

    base = "https://example.invalid/api/v1"
    caps = mc.resolve("kit", "kit.glm-5.3", base, api_key="k")
    assert caps.context_window == 524_288
    assert len(calls) == 1

    # The roster is now known without a second request -- and the speech
    # model is not offered as a chat model.
    assert mc.kit_models_served(base) == {"kit.glm-5.3"}
    assert len(calls) == 1


def test_the_warning_text_does_not_go_to_the_wire(monkeypatch):
    """A warning that blocks for the length of an HTTP timeout is worse
    than a slightly stale one."""
    def _boom(*a, **k):
        raise AssertionError("kit_recommendation reached the network")

    monkeypatch.setattr(mc.urllib.request, "urlopen", _boom)
    mc._SERVED_CACHE.clear()
    caps = mc.resolve("kit", "kit.qwen3.5-397b-A17b", "", allow_live=False)
    msg = mc.kit_recommendation("kit.qwen3.5-397b-A17b", caps)
    assert "no longer served" in msg
    assert mc.KIT_BEST_MODEL in msg


def test_a_failed_probe_is_not_remembered_as_an_empty_roster(monkeypatch):
    def _boom(*a, **k):
        raise OSError("endpoint down")

    monkeypatch.setattr(mc.urllib.request, "urlopen", _boom)
    mc._SERVED_CACHE.clear()
    base = "https://example.invalid/api/v1"
    assert mc.kit_models_served(base, fetch=True) == set()
    assert base not in mc._SERVED_CACHE
