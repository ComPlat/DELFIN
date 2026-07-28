"""Cheap-tier model routing for read-only subagents.

Read-only presets (explore / plan / code-reviewer) do their work fine on
the provider's cheap-tier model; running them on the parent's expensive
model burns budget for nothing. ``run_subagent`` therefore switches the
ISOLATED client copy to ``tier_model(provider, "cheap")`` when:

  - the preset is read-only (writers always keep the parent model),
  - the ``agent.subagents.cheap_tier`` setting is on (default ON) or the
    per-call override ``model="cheap"`` is passed,
  - the candidate resolves, differs from the parent model, is not marked
    broken, and supports tool calling.

Any failure silently keeps the parent model — routing must never break a
subagent run. The model actually used is recorded in the telemetry.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from types import SimpleNamespace

import pytest

from delfin.agent import model_capabilities as MC
from delfin.agent import model_routing as MR
from delfin.agent import subagents as SA


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _isolate_routing_state(monkeypatch, tmp_path):
    """Keep the runtime broken-model list and capability resolution away
    from the real ``~/.delfin`` / network. Capabilities default to a
    tool-capable model; individual tests override."""
    monkeypatch.setattr(MR, "_RUNTIME_BROKEN_PATH",
                        tmp_path / "broken_models.json")
    monkeypatch.setattr(
        MC, "resolve",
        lambda *a, **kw: SimpleNamespace(supports_tools=True),
    )
    yield


def _settings(cheap_tier=None, cheap_model="small-model"):
    """Settings dict: per-provider cheap tier + optional subagents flag."""
    sub: dict = {}
    if cheap_tier is not None:
        sub["cheap_tier"] = cheap_tier
    routing: dict = {}
    if cheap_model:
        routing = {"providers": {"prov": {"cheap_model": cheap_model}}}
    return {"agent": {"subagents": sub, "routing": routing}}


def _patch_settings(monkeypatch, settings: dict):
    monkeypatch.setattr("delfin.user_settings.load_settings",
                        lambda: settings)


@dataclass
class _Ev:
    type: str = ""
    text: str = ""
    tool_name: str = ""
    tool_input: dict | None = None
    tool_output: str = ""
    input_tokens: int = 0
    output_tokens: int = 0


class _FakeClient:
    """Minimal parent client. ``stream_models`` is a shared (shallow-copied)
    list, so the model the COPY streamed with is visible on the parent."""

    def __init__(self, model="big-model", provider="prov"):
        self.model = model
        self._provider = provider
        self._base_url = ""
        self._permissions = None
        self.stream_models: list[str] = []

    def set_permissions(self, p):
        self._permissions = p

    def switch_model(self, m):
        if m and m != self.model:
            self.model = m

    def stream_message(self, messages, system, max_tokens):
        self.stream_models.append(self.model)
        yield _Ev(type="text_delta", text="ok done")
        yield _Ev(type="message_delta", input_tokens=10, output_tokens=5)


class _NoSwitchClient:
    """Client WITHOUT ``switch_model`` — the model attribute is set directly."""

    def __init__(self, model="big-model", provider="prov"):
        self.model = model
        self._provider = provider
        self._base_url = ""
        self._permissions = None
        self.stream_models: list[str] = []

    def set_permissions(self, p):
        self._permissions = p

    def stream_message(self, messages, system, max_tokens):
        self.stream_models.append(self.model)
        yield _Ev(type="text_delta", text="ok done")
        yield _Ev(type="message_delta", input_tokens=10, output_tokens=5)


def _run(parent, subagent_type="explore", **kw):
    return SA.run_subagent(
        subagent_type=subagent_type,
        description="probe",
        prompt="self-contained briefing of more than 20 chars",
        parent_client=parent,
        parent_perms=None,
        **kw,
    )


# ---------------------------------------------------------------------------
# Routing matrix
# ---------------------------------------------------------------------------

def test_readonly_preset_routes_to_cheap_tier(monkeypatch):
    _patch_settings(monkeypatch, _settings())
    parent = _FakeClient()
    res = _run(parent, "explore")
    assert res.error == ""
    assert parent.stream_models == ["small-model"]
    assert res.model == "small-model"
    assert res.model_tier == "cheap"
    # The PARENT client is never mutated — only the isolated copy switched.
    assert parent.model == "big-model"
    assert res.to_payload()["model"] == "small-model"


def test_writer_preset_never_switches(monkeypatch):
    _patch_settings(monkeypatch, _settings())
    parent = _FakeClient()
    res = _run(parent, "general-purpose")
    assert res.error == ""
    assert parent.stream_models == ["big-model"]
    assert res.model == "big-model"
    assert res.model_tier == "parent"


def test_writer_preset_ignores_cheap_override(monkeypatch):
    _patch_settings(monkeypatch, _settings())
    parent = _FakeClient()
    res = _run(parent, "general-purpose", model="cheap")
    assert parent.stream_models == ["big-model"]
    assert res.model_tier == "parent"


def test_cheap_tier_setting_disabled_keeps_parent(monkeypatch):
    _patch_settings(monkeypatch, _settings(cheap_tier=False))
    parent = _FakeClient()
    res = _run(parent, "explore")
    assert parent.stream_models == ["big-model"]
    assert res.model == "big-model"
    assert res.model_tier == "parent"


def test_override_cheap_beats_disabled_setting(monkeypatch):
    _patch_settings(monkeypatch, _settings(cheap_tier=False))
    parent = _FakeClient()
    res = _run(parent, "explore", model="cheap")
    assert parent.stream_models == ["small-model"]
    assert res.model_tier == "cheap"


def test_override_parent_pins_parent_model(monkeypatch):
    _patch_settings(monkeypatch, _settings())
    parent = _FakeClient()
    res = _run(parent, "explore", model="parent")
    assert parent.stream_models == ["big-model"]
    assert res.model == "big-model"
    assert res.model_tier == "parent"


# ---------------------------------------------------------------------------
# Fallbacks — routing must never break a run
# ---------------------------------------------------------------------------

def test_unresolvable_tier_falls_back_to_parent(monkeypatch):
    # No cheap model configured and provider "prov" has no built-in tier.
    _patch_settings(monkeypatch, _settings(cheap_model=""))
    parent = _FakeClient()
    res = _run(parent, "explore")
    assert res.error == ""
    assert parent.stream_models == ["big-model"]
    assert res.model_tier == "parent"


def test_broken_cheap_model_falls_back_to_parent(monkeypatch, tmp_path):
    _patch_settings(monkeypatch, _settings())
    broken = tmp_path / "broken_models.json"
    monkeypatch.setattr(MR, "_RUNTIME_BROKEN_PATH", broken)
    broken.write_text(json.dumps(["small-model"]), encoding="utf-8")
    parent = _FakeClient()
    res = _run(parent, "explore")
    assert res.error == ""
    assert parent.stream_models == ["big-model"]
    assert res.model_tier == "parent"


def test_cheap_model_without_tool_support_falls_back(monkeypatch):
    _patch_settings(monkeypatch, _settings())
    monkeypatch.setattr(
        MC, "resolve",
        lambda *a, **kw: SimpleNamespace(supports_tools=False),
    )
    parent = _FakeClient()
    res = _run(parent, "explore")
    assert res.error == ""
    assert parent.stream_models == ["big-model"]
    assert res.model_tier == "parent"


def test_tier_resolution_raising_never_breaks_run(monkeypatch):
    _patch_settings(monkeypatch, _settings())

    def _boom(*a, **kw):
        raise RuntimeError("routing table on fire")

    monkeypatch.setattr(MR, "tier_model", _boom)
    parent = _FakeClient()
    res = _run(parent, "explore")
    assert res.error == ""
    assert parent.stream_models == ["big-model"]
    assert res.model_tier == "parent"


def test_client_without_provider_keeps_parent(monkeypatch):
    _patch_settings(monkeypatch, _settings())
    parent = _FakeClient(provider="")
    res = _run(parent, "explore")
    assert res.error == ""
    assert parent.stream_models == ["big-model"]
    assert res.model_tier == "parent"


def test_client_without_switch_model_gets_model_attribute(monkeypatch):
    _patch_settings(monkeypatch, _settings())
    parent = _NoSwitchClient()
    assert not hasattr(parent, "switch_model")
    res = _run(parent, "explore")
    assert res.error == ""
    assert parent.stream_models == ["small-model"]
    assert res.model_tier == "cheap"
    assert parent.model == "big-model"


# ---------------------------------------------------------------------------
# Telemetry
# ---------------------------------------------------------------------------

def test_telemetry_records_routed_model(monkeypatch):
    _patch_settings(monkeypatch, _settings())
    parent = _FakeClient()
    _run(parent, "explore")
    records = SA.read_telemetry(SA._TELEMETRY_PATH)
    assert records, "run must append a telemetry record"
    assert records[-1]["model"] == "small-model"
    assert records[-1]["model_tier"] == "cheap"


def test_telemetry_records_parent_model_when_not_routed(monkeypatch):
    _patch_settings(monkeypatch, _settings())
    parent = _FakeClient()
    _run(parent, "general-purpose")
    records = SA.read_telemetry(SA._TELEMETRY_PATH)
    assert records[-1]["model"] == "big-model"
    assert records[-1]["model_tier"] == "parent"


# ---------------------------------------------------------------------------
# Resolver unit checks
# ---------------------------------------------------------------------------

def test_resolver_same_model_as_parent_keeps_parent(monkeypatch):
    _patch_settings(monkeypatch, _settings(cheap_model="big-model"))
    parent = _FakeClient(model="big-model")
    model, tier = SA._resolve_subagent_model(parent, "explore")
    assert (model, tier) == ("big-model", "parent")


def test_resolver_cheap_tier_default_is_enabled(monkeypatch):
    _patch_settings(monkeypatch, _settings())  # no cheap_tier key at all
    parent = _FakeClient()
    model, tier = SA._resolve_subagent_model(parent, "explore")
    assert (model, tier) == ("small-model", "cheap")
