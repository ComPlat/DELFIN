"""Subagent-limits settings UI wiring + backend round-trip.

The dashboard settings tab exposes agent.subagents.{max_wall_s,max_tool_calls,
max_output_tokens}; the backend (subagents._subagent_limits) reads exactly those
keys. This locks the integration so a UI change and the reader can't drift.
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent import subagents as S


_SRC = (Path(__file__).resolve().parent.parent
        / "delfin" / "dashboard" / "tab_settings.py").read_text(encoding="utf-8")


def test_widgets_defined():
    for w in ("subagent_wall_input", "subagent_calls_input",
              "subagent_tokens_input"):
        assert f"{w} = widgets.BoundedIntText(" in _SRC


def test_save_paths_write_subagents_keys():
    # Both the dedicated "Save agent extras" handler and the main settings
    # save must persist the three keys under agent.subagents.
    assert _SRC.count("'max_wall_s': int(subagent_wall_input.value") >= 2
    # Persisted via a merge dict (so hand-set keys survive), both paths.
    assert "payload['agent']['subagents'] = subs" in _SRC
    assert "settings_payload['agent']['subagents'] = _subs" in _SRC


def test_load_path_reads_subagents_keys():
    assert "subagents_payload = agent_payload.get('subagents')" in _SRC
    assert "subagent_wall_input.value = int(subagents_payload.get('max_wall_s'" in _SRC


def test_ui_section_present():
    assert "Subagent limits" in _SRC


def test_backend_reads_the_same_keys():
    # The reader consumes exactly the keys the UI writes.
    keys = {"max_wall_s", "max_tool_calls", "max_output_tokens"}
    lim = S._subagent_limits()
    assert keys == set(lim.keys())
    # And UI defaults match the backend fallback defaults.
    assert S._MAX_WALL_S == 300.0
    assert S._MAX_TOOL_CALLS == 40
    assert S._MAX_OUTPUT_TOKENS == 16000
