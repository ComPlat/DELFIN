"""Shell output reaches the model history without colour codes, and a
bash call that only reads or searches gets a one-line pointer to the
dedicated tool.

Field data (33 sessions, 2026-09-16..10-14): 149 pytest runs through bash
against 66 through run_tests, with about 97 k characters of ANSI escapes
in the archived outputs, every one of them re-sent with the history on
each later round; 391 of 1 179 bash calls (33 %) only read or searched
(grep 207, sed -n 67, ls 52, cat/head/tail 30, find 15). The hint does
not refuse anything: the result is returned as before.
"""

from __future__ import annotations

import json
from unittest.mock import patch

import pytest

import delfin.agent.api_client as A
from delfin.agent import mcp_client as M
from delfin.agent import model_capabilities as mc
from delfin.agent.tool_output_hygiene import shell_tool_hint, strip_ansi


# ---------------------------------------------------------------------------
# strip_ansi: raw text and the JSON a bash result is serialised as
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw, clean", [
    ("\x1b[31mFAILED\x1b[0m t.py", "FAILED t.py"),
    ("\x1b[1m\x1b[32m3 passed\x1b[0m", "3 passed"),
    ("ok\x1b]8;;http://x\x1b\\link\x1b]8;;\x1b\\", "oklink"),
    ("a\x1b[?25lb", "ab"),
    ("plain\ttext\n", "plain\ttext\n"),
])
def test_escape_sequences_are_removed_from_raw_text(raw, clean):
    assert strip_ansi(raw) == clean


def test_escape_sequences_are_removed_inside_a_json_result():
    raw = json.dumps({"exit_code": 1,
                      "stdout": "\x1b[31mFAILED\x1b[0m tests/t.py::x\n"
                                "\x1b[1m1 failed\x1b[0m",
                      "stderr": ""})
    assert "\\u001b" in raw
    out = strip_ansi(raw)
    assert "\\u001b" not in out
    assert json.loads(out)["stdout"] == "FAILED tests/t.py::x\n1 failed"


def test_a_literal_backslash_u_in_json_is_left_alone():
    raw = json.dumps({"stdout": "the text \\u001b[31m is documented"})
    assert strip_ansi(raw) == raw


# ---------------------------------------------------------------------------
# shell_tool_hint: plain reads/searches and pytest, nothing else
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "grep -n foo delfin/agent/x.py",
    "cd /ws && grep -rn 'def foo' delfin | head -20",
    "rg foo",
    "cat README.md",
    "head -50 a.py && tail -20 b.py",
    "ls -la delfin/agent",
    "find . -name '*.py'",
    "sed -n 10,40p delfin/agent/x.py",
])
def test_a_plain_read_or_search_names_the_dedicated_tool(cmd, tmp_path):
    hint = shell_tool_hint(cmd, tmp_path)
    assert "read_file" in hint and "grep_file" in hint
    assert "\n" not in hint.strip()


@pytest.mark.parametrize("cmd", [
    "python -m pytest tests/test_x.py -q",
    "cd /ws && pytest -x tests 2>&1 | tail -30",
    "PYTHONPATH=. python3 -m pytest -q tests/test_y.py",
])
def test_pytest_through_bash_points_to_run_tests(cmd, tmp_path):
    hint = shell_tool_hint(cmd, tmp_path)
    assert "run_tests" in hint
    assert "\n" not in hint.strip()


@pytest.mark.parametrize("cmd", [
    "git status",
    "python script.py",
    "grep foo a.txt > hits.txt",
    "cat a.py && rm b.py",
    "sed -i 's/a/b/' x.py",
    "find . -name '*.pyc' -delete",
    "find . -exec rm {} +",
    "ls $(cat list.txt)",
    "cat /etc/hosts",
    "grep -r token ~/.config",
    "",
])
def test_anything_else_gets_no_hint(cmd, tmp_path):
    assert shell_tool_hint(cmd, tmp_path) == ""


# ---------------------------------------------------------------------------
# Driven through the real client loop: what the NEXT request carries
# ---------------------------------------------------------------------------

class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content = content
        self.tool_calls = tool_calls
        self.reasoning_content = None


class _Choice:
    def __init__(self, delta, finish=None):
        self.index = 0
        self.delta = delta
        self.finish_reason = finish


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


class _ToolCall:
    def __init__(self, name, arguments, cid):
        self.index = 0
        self.id = cid
        self.type = "function"
        self.function = type("F", (), {"name": name,
                                       "arguments": arguments})()


class _FakeRegistry:
    def discover_all(self):
        return []

    def discover_resources(self):
        return []

    def discover_prompts(self):
        return []


def _drive(ws, commands):
    """Run each bash command as one tool round; return the tool rows the
    final request carried."""
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    seen: list = []
    with patch.object(mc, "resolve", lambda *a, **k: caps), \
         patch.object(M, "get_registry", lambda *a, **k: _FakeRegistry()), \
         patch.object(A.time, "sleep", lambda *a, **k: None):
        client = A.create_client(backend="api", provider="ollama",
                                 model="qwen2.5-coder:7b", cwd=str(ws))
        client.set_permissions(A.KitToolPermissions(
            workspace=ws, mode="bypassPermissions"))
        n = {"i": 0}

        def _create(**kwargs):
            seen.append(kwargs.get("messages") or [])
            i = n["i"]
            n["i"] += 1
            if i < len(commands):
                return _Stream([
                    _Chunk([_Choice(_Delta(tool_calls=[_ToolCall(
                        "bash", json.dumps({"command": commands[i]}),
                        f"c{i}")]))]),
                    _Chunk([_Choice(_Delta(), finish="tool_calls")]),
                ])
            return _Stream([_Chunk([_Choice(_Delta(content="done"),
                                            finish="stop")])])

        client.client.chat.completions.create = _create
        list(client.stream_message(
            "sys", [{"role": "user", "content": "go"}], max_tokens=64))
    return [m for m in seen[-1] if isinstance(m, dict)
            and m.get("role") == "tool"]


def test_driven_colour_is_gone_and_the_read_hint_is_attached(tmp_path):
    (tmp_path / "log.txt").write_text(
        "\x1b[31mFAILED\x1b[0m tests/t.py::x\n\x1b[1m1 failed\x1b[0m\n",
        encoding="utf-8")
    rows = _drive(tmp_path, ["grep -n failed log.txt",
                             "grep -n FAILED log.txt"])
    assert len(rows) == 2
    for row in rows:
        content = str(row["content"])
        assert "\x1b" not in content and "\\u001b" not in content, content
        assert "FAILED" in content or "failed" in content
    # One hint per kind and turn: the first read carries it, the second
    # does not pay for it again.
    assert "grep_file" in str(rows[0]["content"])
    assert "grep_file" not in str(rows[1]["content"])
