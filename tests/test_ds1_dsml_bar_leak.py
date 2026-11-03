"""Task ds1: DeepSeek full-width-bar DSML tool-call markup leak.

deepseek-v4-flash on the KIT endpoint sometimes writes its tool calls as
literal DSML tags with full-width bars (U+FF5C) in the text channel:

    <｜DSML｜tool_calls><｜DSML｜invoke name="bash">
    <｜DSML｜parameter name="command" string="true">git status</｜DSML｜parameter>
    <｜DSML｜/invoke></｜DSML｜tool_calls>

A well-formed block (a tool name plus parseable arguments) must be
RECOVERED as a real tool call; a malformed or partial one must be
STRIPPED from the visible text and must never be executed.  The plain
``<invoke name="...">`` variant keeps its wave-14 semantics: it is a
re-request case, NOT silently redispatched by ``parse_leaked_tool_calls``.
"""

from __future__ import annotations

import json

from delfin.agent.text_sanitize import (
    parse_leaked_tool_calls,
    sanitize_agent_text,
)

_BAR_BASH = (
    '<｜DSML｜tool_calls><｜DSML｜invoke name="bash">\n'
    '<｜DSML｜parameter name="command" string="true">git log --oneline -3</｜DSML｜parameter>\n'
    '<｜DSML｜parameter name="description" string="true">Read the last commits</｜DSML｜parameter>\n'
    '<｜DSML｜/invoke></｜DSML｜tool_calls>'
)

_BAR_JSON_ARGS = (
    '<｜DSML｜invoke name="read_file">'
    '<｜DSML｜arguments>{"path": "delfin/agent/text_sanitize.py"}</｜DSML｜arguments>'
    '<｜DSML｜/invoke>'
)

_BAR_MULTI = (
    '<｜DSML｜tool_calls>'
    '<｜DSML｜invoke name="bash">'
    '<｜DSML｜parameter name="command">git status</｜DSML｜parameter>'
    '<｜DSML｜/invoke>'
    '<｜DSML｜invoke name="read_file">'
    '<｜DSML｜parameter name="path">README.md</｜DSML｜parameter>'
    '<｜DSML｜/invoke>'
    '</｜DSML｜tool_calls>'
)

_ANSWER = "Die CASSCF-Aktivraumwahl steht in Abschnitt 6.15 des Handbuchs."


# --- Recovery: well-formed blocks become real tool calls ---

def test_well_formed_parameter_block_is_recovered():
    calls = parse_leaked_tool_calls(_BAR_BASH)
    assert len(calls) == 1
    assert calls[0]["name"] == "bash"
    assert calls[0]["arguments"] == {
        "command": "git log --oneline -3",
        "description": "Read the last commits",
    }


def test_well_formed_json_arguments_block_is_recovered():
    calls = parse_leaked_tool_calls(_BAR_JSON_ARGS)
    assert calls == [{"name": "read_file",
                      "arguments": {"path": "delfin/agent/text_sanitize.py"}}]


def test_wrapper_with_two_invokes_recovers_both():
    calls = parse_leaked_tool_calls(_BAR_MULTI)
    assert [c["name"] for c in calls] == ["bash", "read_file"]
    assert calls[0]["arguments"] == {"command": "git status"}
    assert calls[1]["arguments"] == {"path": "README.md"}


def test_recovery_survives_answer_interleaved():
    text = _BAR_BASH + "\n\n" + _ANSWER
    calls = parse_leaked_tool_calls(text)
    assert [c["name"] for c in calls] == ["bash"]
    assert calls[0]["arguments"]["command"].startswith("git log")


def test_recovery_with_json_arguments_does_not_read_parameter_pairs():
    # When the block carries a JSON arguments object, that is THE argument
    # payload; a parameter pair with the same name must not double it up.
    text = (
        '<｜DSML｜invoke name="write_file">'
        '<｜DSML｜arguments>{"path": "a.py", "content": "x"}</｜DSML｜arguments>'
        '<｜DSML｜/invoke>'
    )
    calls = parse_leaked_tool_calls(text)
    assert calls == [{"name": "write_file",
                      "arguments": {"path": "a.py", "content": "x"}}]


# --- Never execute malformed or partial markup ---

def test_truncated_block_is_never_recovered():
    # Streaming cut off inside the block: no closing bar tag.
    text = '<｜DSML｜invoke name="bash"><｜DSML｜parameter name="command">git stat'
    assert parse_leaked_tool_calls(text) == []


def test_unclosed_json_arguments_is_never_recovered():
    text = (
        '<｜DSML｜invoke name="read_file">'
        '<｜DSML｜arguments>{"path": "a.py"'
    )
    assert parse_leaked_tool_calls(text) == []


def test_broken_json_arguments_is_never_recovered():
    text = (
        '<｜DSML｜invoke name="read_file">'
        '<｜DSML｜arguments>{path: a.py}</｜DSML｜arguments>'
        '<｜DSML｜/invoke>'
    )
    assert parse_leaked_tool_calls(text) == []


def test_invoke_without_name_is_never_recovered():
    text = (
        '<｜DSML｜invoke><｜DSML｜parameter name="command">git status'
        '</｜DSML｜parameter><｜DSML｜/invoke>'
    )
    assert parse_leaked_tool_calls(text) == []


# --- Stripping: malformed markup leaves the visible text clean ---

def test_well_formed_block_is_stripped_and_answer_survives():
    res = sanitize_agent_text(_BAR_BASH + "\n\n" + _ANSWER)
    assert "｜DSML｜" not in res.text
    assert "DSML" not in res.text
    assert "invoke" not in res.text
    assert _ANSWER in res.text
    assert "bash" in res.leaked_tools
    assert res.dsml_bar_stripped is True


def test_malformed_block_is_stripped_not_shown():
    text = ('Ich prüfe das jetzt.\n'
            '<｜DSML｜invoke name="bash"><｜DSML｜parameter name="command">git stat')
    res = sanitize_agent_text(text)
    assert "｜" not in res.text
    assert "DSML" not in res.text
    assert "Ich prüfe das jetzt." in res.text
    assert res.dsml_bar_stripped is True


def test_wrapper_tags_alone_are_stripped():
    text = '<｜DSML｜tool_calls>\n' + _ANSWER + '\n</｜DSML｜tool_calls>'
    res = sanitize_agent_text(text)
    assert "｜" not in res.text
    assert _ANSWER in res.text


def test_clean_text_untouched_by_bar_path():
    res = sanitize_agent_text(_ANSWER)
    assert res.text == _ANSWER
    assert res.changed is False
    assert res.dsml_bar_stripped is False


# --- Streamed chunk boundaries inside the markup ---

def _chunks_at(text, positions):
    return [text[i:j] for i, j in zip([0] + positions,
                                      positions + [len(text)])]


def test_recovery_identical_for_every_chunk_boundary():
    for cut in range(1, len(_BAR_BASH)):
        chunks = _chunks_at(_BAR_BASH, [cut])
        joined = "".join(chunks)
        assert joined == _BAR_BASH
        calls = parse_leaked_tool_calls(joined)
        assert [c["name"] for c in calls] == ["bash"]
        assert calls[0]["arguments"] == {
            "command": "git log --oneline -3",
            "description": "Read the last commits",
        }


def test_json_recovery_identical_for_every_chunk_boundary():
    for cut in range(1, len(_BAR_JSON_ARGS)):
        joined = "".join(_chunks_at(_BAR_JSON_ARGS, [cut]))
        calls = parse_leaked_tool_calls(joined)
        assert calls == [{"name": "read_file",
                          "arguments": {"path": "delfin/agent/text_sanitize.py"}}]


def test_stream_chunks_produce_one_recovered_call_not_many():
    # A pathological per-character stream still yields exactly ONE call.
    chars = list(_BAR_BASH)
    joined = "".join(chars)
    calls = parse_leaked_tool_calls(joined)
    assert len(calls) == 1
    assert json.dumps(calls[0]["arguments"]) .count("command") == 1


# --- The plain <invoke> variant keeps its wave-14 semantics ---

def test_plain_invoke_variant_is_not_silently_redispatched():
    text = ('<invoke name="bash">'
            '<parameter name="command">git status</parameter>'
            '</invoke>')
    assert parse_leaked_tool_calls(text) == []


def test_harmony_and_qwen_grammars_unaffected():
    harmony = 'to=search_docs {"query": "CASSCF"}'
    qwen = '<tool_call>{"name": "read_file", "arguments": {"path": "x"}}</tool_call>'
    h = parse_leaked_tool_calls(harmony)
    q = parse_leaked_tool_calls(qwen)
    assert h and h[0]["name"] == "search_docs"
    assert q and q[0]["name"] == "read_file"


def test_two_different_calls_to_one_tool_are_both_recovered():
    """Two bash calls with different commands are two calls; only an exact
    repeat collapses."""
    from delfin.agent.text_sanitize import parse_leaked_tool_calls

    def block(cmd):
        return ('<｜DSML｜tool_calls><｜DSML｜invoke name="bash">'
                f'<｜DSML｜parameter name="command" string="true">{cmd}'
                '</｜DSML｜parameter><｜DSML｜/invoke></｜DSML｜tool_calls>')

    calls = parse_leaked_tool_calls(block("git status") + block("git log -1"))
    assert [c["arguments"]["command"] for c in calls] == ["git status", "git log -1"]
    repeat = parse_leaked_tool_calls(block("git status") + block("git status"))
    assert len(repeat) == 1
