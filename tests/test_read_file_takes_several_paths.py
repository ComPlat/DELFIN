"""read_file reads several files in one call, under every per-path gate.

Input: a read_file call whose ``paths`` lists up to five files (``path``,
when given, counts as the first). Output: each file's text, or the error a
single read of that file returns, in the same request. Every mechanism a
single read passes through (secret deny-globs, locked scope, the read
baseline that edit_file requires, the reread note) holds for every path.

Driven through the streaming tool loop (the model's entry point) and through
the executor's public entry (Code Mode's entry point).
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import mcp_client as M


@pytest.fixture
def ws(tmp_path):
    w = tmp_path / "ws"
    w.mkdir()
    (w / "a.py").write_text("ALPHA = 1\n")
    (w / "b.py").write_text("BETA = 2\n")
    (w / ".env").write_text("TOKEN=s3cr3t-value\n")
    (tmp_path / "outside.txt").write_text("OUTSIDE CONTENT\n")
    return w


# -- a scripted model ---------------------------------------------------------

class _Fn:
    def __init__(self, name=None, arguments=None):
        self.name = name
        self.arguments = arguments


class _TC:
    def __init__(self, index, id=None, name=None, arguments=None):
        self.index = index
        self.id = id
        self.function = _Fn(name, arguments)


class _Delta:
    def __init__(self, content=None, tool_calls=None):
        self.content = content
        self.tool_calls = tool_calls


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta = delta
        self.finish_reason = finish


class _Usage:
    prompt_tokens = 5
    completion_tokens = 3


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


class _EmptyReg:
    def discover_all(self): return []
    def discover_resources(self): return []
    def discover_prompts(self): return []


def _call(cid, name, args):
    return _Stream([_Chunk([_Choice(_Delta(tool_calls=[_TC(
        0, id=cid, name=name, arguments=json.dumps(args))]),
        finish="tool_calls")], usage=_Usage())])


def _say(text):
    return _Stream([_Chunk([_Choice(_Delta(content=text), finish="stop")],
                           usage=_Usage())])


def _make_client(monkeypatch, ws, **perm_kw):
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    from delfin.agent.api_client import KitToolPermissions, create_client
    cl = create_client(backend="api", provider="ollama",
                       model="qwen2.5-coder:7b", cwd=str(ws))
    perm_kw.setdefault("mode", "acceptEdits")
    cl.set_permissions(KitToolPermissions(
        workspace=ws, task_session_id="sess-paths", **perm_kw))
    return cl


def _drive(client, streams):
    """Run one turn; return the request lists the model was sent."""
    sent = []

    def _create(**kw):
        sent.append(list(kw.get("messages") or []))
        return streams[min(len(sent) - 1, len(streams) - 1)]

    client.client.chat.completions.create = _create
    list(client.stream_message(
        "sys", [{"role": "user", "content": "lies die Dateien"}],
        max_tokens=100))
    return sent


def _rows_after(messages, n_assistant):
    """The tool rows answering the n-th assistant message (1-based)."""
    seen = 0
    rows = []
    for m in messages:
        if m.get("role") == "assistant" and m.get("tool_calls"):
            seen += 1
            continue
        if seen == n_assistant and m.get("role") == "tool":
            rows.append(str(m.get("content")))
        elif seen > n_assistant:
            break
    return rows


def _single_error(ws, perms_kw, path):
    """What a single read_file of *path* answers, byte for byte."""
    from delfin.agent.api_client import KitToolPermissions, _doc_executor
    kw = dict(perms_kw)
    kw.setdefault("mode", "acceptEdits")
    perms = KitToolPermissions(workspace=ws, task_session_id="s-single", **kw)
    out = _doc_executor.execute("read_file", {"path": path}, perms)
    assert out.startswith('{"error"'), out
    return out


# -- through the tool loop ----------------------------------------------------

def test_two_files_are_read_in_one_request(monkeypatch, ws):
    cl = _make_client(monkeypatch, ws)
    sent = _drive(cl, [_call("c1", "read_file", {"paths": ["a.py", "b.py"]}),
                       _say("fertig")])
    assert len(sent) == 2          # one model request for both files
    text = "\n".join(_rows_after(sent[-1], 1))
    assert "ALPHA = 1" in text and "BETA = 2" in text


def test_path_and_paths_combine(monkeypatch, ws):
    cl = _make_client(monkeypatch, ws)
    sent = _drive(cl, [_call("c1", "read_file",
                             {"path": "a.py", "paths": ["b.py"]}),
                       _say("fertig")])
    text = "\n".join(_rows_after(sent[-1], 1))
    assert "ALPHA = 1" in text and "BETA = 2" in text


def test_a_secret_among_them_is_denied_and_the_other_returned(
        monkeypatch, ws):
    cl = _make_client(monkeypatch, ws)
    sent = _drive(cl, [_call("c1", "read_file", {"paths": [".env", "a.py"]}),
                       _say("fertig")])
    text = "\n".join(_rows_after(sent[-1], 1))
    assert "s3cr3t-value" not in text
    assert "ALPHA = 1" in text
    assert _single_error(ws, {}, ".env") in text


def test_an_outside_path_in_a_locked_scope_is_refused(monkeypatch, ws):
    outside = str(ws.parent / "outside.txt")
    cl = _make_client(monkeypatch, ws, lock_workspace=True)
    sent = _drive(cl, [_call("c1", "read_file",
                             {"paths": ["a.py", outside]}),
                       _say("fertig")])
    text = "\n".join(_rows_after(sent[-1], 1))
    assert "OUTSIDE CONTENT" not in text
    assert "ALPHA = 1" in text
    assert _single_error(ws, {"lock_workspace": True}, outside) in text


def test_a_file_read_through_paths_can_be_edited(monkeypatch, ws):
    cl = _make_client(monkeypatch, ws)
    sent = _drive(cl, [
        _call("c1", "read_file", {"paths": ["a.py", "b.py"]}),
        _call("c2", "edit_file", {"path": "b.py", "old_string": "BETA = 2",
                                  "new_string": "BETA = 3"}),
        _say("fertig")])
    edit_rows = _rows_after(sent[-1], 2)
    assert edit_rows and "before editing" not in edit_rows[0], edit_rows
    assert (ws / "b.py").read_text() == "BETA = 3\n"


def test_a_file_not_read_still_cannot_be_edited(monkeypatch, ws):
    """The control: the read baseline is per file, not per call."""
    (ws / "c.py").write_text("GAMMA = 3\n")
    cl = _make_client(monkeypatch, ws)
    sent = _drive(cl, [
        _call("c1", "read_file", {"paths": ["a.py", "b.py"]}),
        _call("c2", "edit_file", {"path": "c.py", "old_string": "GAMMA = 3",
                                  "new_string": "GAMMA = 4"}),
        _say("fertig")])
    edit_rows = _rows_after(sent[-1], 2)
    assert "before editing" in edit_rows[0]
    assert (ws / "c.py").read_text() == "GAMMA = 3\n"


def test_a_second_read_of_the_same_files_gets_the_reread_note(
        monkeypatch, ws):
    cl = _make_client(monkeypatch, ws)
    sent = _drive(cl, [
        _call("c1", "read_file", {"paths": ["a.py", "b.py"]}),
        _call("c2", "read_file", {"paths": ["a.py", "b.py"]}),
        _say("fertig")])
    text = "\n".join(_rows_after(sent[-1], 2))
    assert "unchanged" in text
    assert "ALPHA = 1" not in text and "BETA = 2" not in text


def test_more_than_five_files_are_refused(monkeypatch, ws):
    for i in range(6):
        (ws / f"f{i}.py").write_text(f"F{i} = {i}\n")
    cl = _make_client(monkeypatch, ws)
    sent = _drive(cl, [_call("c1", "read_file",
                             {"paths": [f"f{i}.py" for i in range(6)]}),
                       _say("fertig")])
    text = "\n".join(_rows_after(sent[-1], 1))
    assert "at most 5" in text
    assert "F0 = 0" not in text


def test_a_single_path_reads_as_before(monkeypatch, ws):
    cl = _make_client(monkeypatch, ws)
    sent = _drive(cl, [_call("c1", "read_file", {"path": "a.py"}),
                       _say("fertig")])
    assistant = [m for m in sent[-1] if m.get("role") == "assistant"]
    assert [tc["function"]["name"] for tc in assistant[0]["tool_calls"]] \
        == ["read_file"]
    assert assistant[0]["tool_calls"][0]["id"] == "c1"
    rows = _rows_after(sent[-1], 1)
    assert len(rows) == 1 and "ALPHA = 1" in rows[0]


# -- through the executor's public entry (Code Mode) --------------------------

def _perms(ws, **kw):
    from delfin.agent.api_client import KitToolPermissions
    kw.setdefault("mode", "acceptEdits")
    return KitToolPermissions(workspace=ws, task_session_id="s-exec", **kw)


def test_the_executor_reads_each_path_under_its_gate(ws):
    from delfin.agent.api_client import _doc_executor
    perms = _perms(ws)
    out = _doc_executor.execute(
        "read_file", {"paths": ["a.py", ".env", "b.py"]}, perms)
    assert "ALPHA = 1" in out and "BETA = 2" in out
    assert "s3cr3t-value" not in out
    assert _single_error(ws, {}, ".env") in out
    # each header names its path
    for p in ("a.py", ".env", "b.py"):
        assert f"=== {p} ===" in out
    # the read baseline exists for the files that were read
    assert str((ws / "a.py").resolve()) in perms.read_tracker
    assert str((ws / ".env").resolve()) not in perms.read_tracker


def test_the_executor_refuses_an_outside_path_in_a_locked_scope(ws):
    from delfin.agent.api_client import _doc_executor
    outside = str(ws.parent / "outside.txt")
    out = _doc_executor.execute(
        "read_file", {"paths": [outside, "a.py"]},
        _perms(ws, lock_workspace=True))
    assert "OUTSIDE CONTENT" not in out and "ALPHA = 1" in out
    assert _single_error(ws, {"lock_workspace": True}, outside) in out


def test_the_executor_refuses_more_than_five(ws):
    from delfin.agent.api_client import _doc_executor
    out = _doc_executor.execute(
        "read_file", {"paths": ["a.py"] * 3 + ["b.py"] * 3,
                      "path": "x.py"}, _perms(ws))
    assert json.loads(out)["error"].startswith("read_file reads at most 5")


@pytest.mark.parametrize("bad", ["a.py", [1, 2], [""], []])
def test_a_malformed_paths_argument_is_refused(ws, bad):
    from delfin.agent.api_client import _doc_executor
    out = _doc_executor.execute("read_file", {"paths": bad}, _perms(ws))
    assert "paths must be" in json.loads(out)["error"]


def test_the_schema_offers_paths():
    from delfin.agent.api_client import _DOC_TOOLS_OPENAI
    spec = next(t["function"] for t in _DOC_TOOLS_OPENAI
                if t.get("function", {}).get("name") == "read_file")
    props = spec["parameters"]["properties"]
    assert props["paths"]["type"] == "array"
    assert props["paths"]["items"] == {"type": "string"}
