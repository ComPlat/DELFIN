"""Tests for harmony-leak / glitch-token sanitisation.

Anchored on the real production failure (azure.gpt-5.4 via KIT): the model
emitted ``to=search_docs {json}`` tool-channel syntax as text, interleaved
with Chinese/Korean/Armenian glitch tokens.
"""

from __future__ import annotations

from delfin.agent.text_sanitize import (
    sanitize_agent_text,
    parse_leaked_tool_calls,
)


# The exact shape observed in the bug report.
_BAD = (
    'to=search_docs  手机天天中彩票{"query":"CASSCF ORCA setup", "top_k": 5}\n'
    'to=read_section ացինjson_schema{"doc_id":"orca","section_id":"6.15"} ыҟоуп출장샵\n'
    'Ja — für **CASSCF in ORCA** ist die relevante Stelle Abschnitt 6.15.'
)


def test_clean_text_unchanged():
    txt = "Im %casscf-Block: nel, norb, mult. Die Antwort ist 42."
    res = sanitize_agent_text(txt)
    assert res.text == txt
    assert res.changed is False
    assert res.leaked_tools == []
    assert res.glitch_chars == 0


def test_leaked_tools_detected_and_stripped():
    res = sanitize_agent_text(_BAD)
    assert "search_docs" in res.leaked_tools
    assert "read_section" in res.leaked_tools
    assert "to=search_docs" not in res.text
    assert "{" not in res.text                       # leaked JSON removed


def test_glitch_tokens_removed():
    res = sanitize_agent_text(_BAD)
    assert res.glitch_chars > 0
    for ch in ("手", "机", "출", "장", "ա", "ы"):
        assert ch not in res.text


def test_real_answer_survives():
    res = sanitize_agent_text(_BAD)
    assert "CASSCF in ORCA" in res.text
    assert "6.15" in res.text
    assert res.changed is True


def test_recover_intended_tool_calls():
    calls = parse_leaked_tool_calls(_BAD)
    names = [c["name"] for c in calls]
    assert names == ["search_docs", "read_section"]
    assert calls[0]["arguments"]["query"].startswith("CASSCF")
    assert calls[1]["arguments"]["section_id"] == "6.15"


def test_german_umlauts_preserved():
    txt = "Für die Berechnung benötigst du größere Aktivräume — schöne Grüße."
    res = sanitize_agent_text(txt)
    assert res.text == txt                           # umlauts are NOT glitch
    assert res.changed is False


def test_empty_text():
    res = sanitize_agent_text("")
    assert res.text == ""
    assert res.changed is False
