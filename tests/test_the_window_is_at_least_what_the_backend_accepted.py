"""An accepted request is proof of the window, and the floor never lowers it.

The static table understates: 131,072 for kit.deepseek-v4-flash against a
single accepted request of 166,379 tokens, and no entry at all for
kit.glm-5.3-flash, which fell to the 32,768 heuristic while 85,535-token
requests went through (field metrics, 2026-10-08). Every protection of a
long session fires at a share of the window, so a window that is too
small is not conservative -- it throws away context the model would have
held.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import observed_window as OW


@pytest.fixture
def floors(tmp_path, monkeypatch):
    monkeypatch.setattr(OW, "_PATH", tmp_path / "observed_windows.json")
    return tmp_path / "observed_windows.json"


class TestTheFloor:
    def test_nothing_known_is_zero(self, floors):
        assert OW.floor("kit.x") == 0
        assert OW.at_least("kit.x", 131_072) == 131_072

    def test_an_accepted_request_raises_the_window(self, floors):
        OW.note("kit.deepseek-v4-flash", 166_379)
        assert OW.floor("kit.deepseek-v4-flash") == 166_379
        assert OW.at_least("kit.deepseek-v4-flash", 131_072) == 166_379

    def test_it_never_lowers(self, floors):
        OW.note("m", 90_000)
        OW.note("m", 40_000)
        assert OW.floor("m") == 90_000
        assert OW.at_least("m", 131_072) == 131_072

    def test_it_is_per_model(self, floors):
        OW.note("a", 150_000)
        assert OW.floor("b") == 0

    def test_it_persists_owner_only(self, floors):
        OW.note("m", 123_456)
        assert json.loads(floors.read_text())["m"] == 123_456
        assert (floors.stat().st_mode & 0o777) == 0o600

    def test_garbage_is_ignored(self, floors):
        OW.note("", 5); OW.note("m", 0); OW.note("m", "x"); OW.note("m", -3)
        assert OW.floor("m") == 0
        floors.write_text("not json")
        assert OW.floor("m") == 0
        OW.note("m", 7)                       # rewrites the broken file
        assert OW.floor("m") == 7


class TestTheTable:
    def test_the_flash_model_is_no_longer_a_guess(self):
        from delfin.agent import model_capabilities as mc
        caps = mc.resolve("kit", "kit.glm-5.3-flash", "")
        assert caps.source == "static", caps.source
        assert caps.context_window >= 85_535

    def test_an_unknown_model_still_falls_back(self):
        from delfin.agent import model_capabilities as mc
        caps = mc.resolve("kit", "totally-unknown-model-z", "")
        assert caps.source == "heuristic"


class TestTheEngine:
    def test_resolution_applies_the_floor(self, floors, monkeypatch):
        """Through the engine's own resolution path, with the live probe
        answered here: the table says 131k, the backend accepted 166k."""
        from delfin.agent import engine as E
        from delfin.agent import model_capabilities as mc

        class _Caps:
            context_window = 131_072
        monkeypatch.setattr(mc, "resolve", lambda *a, **k: _Caps())
        OW.note("kit.deepseek-v4-flash", 166_379)
        eng = E.AgentEngine.__new__(E.AgentEngine)
        eng.client = type("_C", (), {"_api_key": "",
                                     "model": "kit.deepseek-v4-flash"})()
        eng.provider = "kit"
        eng.context_window_tokens = 100_000
        eng._active_capabilities = None
        eng._refresh_context_window(background=False)
        assert eng.context_window_tokens == 166_379

    def test_without_a_floor_the_table_stands(self, floors, monkeypatch):
        from delfin.agent import engine as E
        from delfin.agent import model_capabilities as mc

        class _Caps:
            context_window = 131_072
        monkeypatch.setattr(mc, "resolve", lambda *a, **k: _Caps())
        eng = E.AgentEngine.__new__(E.AgentEngine)
        eng.client = type("_C", (), {"_api_key": "", "model": "kit.glm-5.3"})()
        eng.provider = "kit"
        eng.context_window_tokens = 100_000
        eng._active_capabilities = None
        eng._refresh_context_window(background=False)
        assert eng.context_window_tokens == 131_072
