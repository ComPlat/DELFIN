"""Weak-model scaffolding: per-model budgets, tool-call text recovery,
argument repair, and near-miss tool-name suggestions.

The per-model knobs (max_tool_rounds, tool_result_cap_kb) existed in
model_profiles but were never read by the tool loop; the text-channel
tool-call recovery only knew the harmony grammar; unparseable argument
JSON silently became {} (misleading downstream errors); unknown tool
names got no suggestion.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as A
from delfin.agent import text_sanitize as ts
from delfin.agent.model_profiles import ModelProfile


# ---------------------------------------------------------------------------
# Per-model budgets
# ---------------------------------------------------------------------------

def _no_settings(monkeypatch):
    from delfin import user_settings
    monkeypatch.setattr(user_settings, "load_settings", lambda: {})


def test_profile_max_tool_rounds_used_without_setting(monkeypatch):
    _no_settings(monkeypatch)
    from delfin.agent import model_profiles as mp
    monkeypatch.setattr(
        mp, "get_profile",
        lambda model, caps=None: ModelProfile(max_tool_rounds=10))
    assert A._resolve_max_tool_rounds("some-weak-model") == 10


def test_explicit_setting_beats_profile(monkeypatch):
    from delfin import user_settings
    monkeypatch.setattr(user_settings, "load_settings",
                        lambda: {"agent": {"max_tool_rounds": 123}})
    from delfin.agent import model_profiles as mp
    monkeypatch.setattr(
        mp, "get_profile",
        lambda model, caps=None: ModelProfile(max_tool_rounds=10))
    assert A._resolve_max_tool_rounds("some-weak-model") == 123


def test_zero_setting_disables_cap(monkeypatch):
    from delfin import user_settings
    monkeypatch.setattr(user_settings, "load_settings",
                        lambda: {"agent": {"max_tool_rounds": 0}})
    assert A._resolve_max_tool_rounds("m") == 100_000


def test_no_model_falls_back_to_legacy_default(monkeypatch):
    _no_settings(monkeypatch)
    assert A._resolve_max_tool_rounds("") == 500


def test_tool_result_cap_from_profile(monkeypatch):
    from delfin.agent import model_profiles as mp
    monkeypatch.setattr(
        mp, "get_profile",
        lambda model, caps=None: ModelProfile(tool_result_cap_kb=3))
    assert A._resolve_tool_result_cap("weak") == 3 * 1024
    assert A._resolve_tool_result_cap("") == 5000


# ---------------------------------------------------------------------------
# Argument JSON repair
# ---------------------------------------------------------------------------

def test_repair_trailing_comma():
    assert A._repair_json_args('{"path": "a.py",}') == {"path": "a.py"}


def test_repair_single_quotes():
    assert A._repair_json_args("{'path': 'a.py'}") == {"path": "a.py"}


def test_repair_fenced_block():
    raw = '```json\n{"query": "x"}\n```'
    assert A._repair_json_args(raw) == {"query": "x"}


def test_repair_python_literal():
    assert A._repair_json_args("{'n': 3, 'flag': True}") == {
        "n": 3, "flag": True}


def test_repair_gives_up_on_garbage():
    assert A._repair_json_args('{"path": "unterminated') is None
    assert A._repair_json_args("") is None
    assert A._repair_json_args(None) is None


# ---------------------------------------------------------------------------
# Text-channel tool-call recovery (qwen/gemma grammars)
# ---------------------------------------------------------------------------

def test_parse_qwen_tool_call_markup():
    text = ('I will read the file now.\n<tool_call>\n'
            '{"name": "read_file", "arguments": {"path": "a.py"}}\n'
            '</tool_call>')
    calls = ts.parse_leaked_tool_calls(text)
    assert calls == [{"name": "read_file", "arguments": {"path": "a.py"}}]


def test_parse_qwen_parameters_key():
    text = ('<tool_call>{"name": "bash", "parameters": '
            '{"command": "ls"}}</tool_call>')
    calls = ts.parse_leaked_tool_calls(text)
    assert calls == [{"name": "bash", "arguments": {"command": "ls"}}]


def test_fenced_json_call_object_recovered():
    text = '```json\n{"name": "read_file", "arguments": {"path": "b.py"}}\n```'
    calls = ts.parse_leaked_tool_calls(text)
    assert calls == [{"name": "read_file", "arguments": {"path": "b.py"}}]


def test_ordinary_json_output_is_never_a_call():
    """Printed JSON that merely CONTAINS a name field must not execute."""
    text = ('```json\n{"name": "benzene", "arguments": {"x": 1}, '
            '"extra": true}\n```')
    assert ts.parse_leaked_tool_calls(text) == []
    text2 = '```json\n{"name": "benzene", "smiles": "c1ccccc1"}\n```'
    assert ts.parse_leaked_tool_calls(text2) == []


def test_sanitize_strips_qwen_markup_and_reports_name():
    text = ('done.\n<tool_call>{"name": "bash", "arguments": '
            '{"command": "ls"}}</tool_call>')
    res = ts.sanitize_agent_text(text)
    assert "tool_call" not in res.text
    assert "bash" in res.leaked_tools


# ---------------------------------------------------------------------------
# Unknown-tool near-miss suggestion
# ---------------------------------------------------------------------------

def test_unknown_tool_suggests_near_miss():
    out = json.loads(A._doc_executor.execute("read_fle", {}))
    assert "Unknown tool" in out["error"]
    assert "read_file" in out["error"]


def test_unknown_tool_without_near_miss_has_no_hint():
    out = json.loads(A._doc_executor.execute("zzqqxx", {}))
    assert "Unknown tool" in out["error"]
    assert "Did you mean" not in out["error"]
