"""A page, a script or a config the agent edits is checked like code.

Session 78c36f8c (2026-10-19) built a browser game in one HTML file over two
hours. The turn-end check looked only at ``.py`` files, so not one edit was
checked by the harness; the model, which could not open a browser, asked the
user to test each one. artifact_check checks each file by its kind.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import artifact_check as AC
from delfin.agent import mcp_client as M

# A node of the test's own on PATH, so the checks run on every machine
# (tests/test_a_test_is_bound_to_delfin_not_to_a_machine.py: supply it).
# `node --check` replays what node 18.19.1 printed for these two syntax
# errors, recorded 2026-10-19 -- the location line, the source line, the
# message and the stack trace artifact_check must drop -- and passes every
# other file.
_NODE_STUB = """#!{python}
import sys
args = sys.argv[1:]
if args[:1] != ["--check"] or len(args) < 2:
    sys.exit(0)
path = args[1]
for n, line in enumerate(open(path, encoding="utf-8").read().splitlines(), 1):
    if "function loop( {{" in line or line.strip() == "let = ;":
        sys.stderr.write(
            f"{{path}}:{{n}}\\n{{line}}\\n\\n\\n"
            "SyntaxError: Unexpected token '{{'\\n"
            "    at internalCompileFunction (node:internal/vm:73:18)\\n"
            "    at wrapSafe (node:internal/modules/cjs/loader:1274:20)\\n"
            "\\nNode.js v18.19.1\\n")
        sys.exit(1)
sys.exit(0)
"""


@pytest.fixture(autouse=True)
def node_on_path(tmp_path_factory, monkeypatch):
    import os
    import sys
    bin_dir = tmp_path_factory.mktemp("bin")
    node = bin_dir / "node"
    node.write_text(_NODE_STUB.format(python=sys.executable))
    node.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")

GOOD = """<html><head><link rel="stylesheet" href="style.css"></head>
<body>
<script>
const BX = 10;
function loop() { return BX; }
</script>
<img src="https://example.org/x.png">
</body></html>
"""
BROKEN_SCRIPT = GOOD.replace("function loop() {", "function loop( {")


@pytest.fixture
def ws(tmp_path):
    (tmp_path / "style.css").write_text("body{}")
    return tmp_path


def test_a_sound_page_has_no_problem(ws):
    (ws / "index.html").write_text(GOOD)
    assert AC.problems([ws / "index.html"]) == ""


def test_a_broken_inline_script_is_named_at_its_html_line(ws):
    (ws / "index.html").write_text(BROKEN_SCRIPT)
    out = AC.problems([ws / "index.html"])
    assert "inline <script> at line 3" in out
    assert "index.html:5" in out and "SyntaxError" in out
    assert "    at " not in out            # no node stack trace


def test_a_missing_local_file_is_named_and_a_url_is_not(ws):
    (ws / "style.css").unlink()
    (ws / "index.html").write_text(GOOD)
    out = AC.problems([ws / "index.html"])
    assert "references style.css, which does not exist" in out
    assert "example.org" not in out


def test_a_module_script_is_checked_as_a_module(ws):
    (ws / "m.html").write_text(
        '<script type="module">\nimport x from "./x.js";\nexport {};\n'
        '</script>\n<script type="text/x-template">not js {</script>')
    assert AC.problems([ws / "m.html"]) == ""


def test_a_js_file_json_and_shell_are_checked(ws):
    (ws / "a.js").write_text("let = ;")
    (ws / "c.json").write_text('{"a": 1,}')
    (ws / "s.sh").write_text("if then fi\n")
    out = AC.problems([ws / "a.js", ws / "c.json", ws / "s.sh"])
    assert "SyntaxError" in out and "invalid JSON" in out and "s.sh" in out


def test_a_python_file_is_still_compiled(ws):
    (ws / "m.py").write_text("def f(:\n")
    assert "m.py" in AC.problems([ws / "m.py"])


def test_what_cannot_be_checked_is_left_alone(ws):
    (ws / "notes.md").write_text("# hi")
    (ws / "x.css").write_text("{{{")
    assert AC.problems([ws / "notes.md", ws / "x.css", ws / "gone.js"]) == ""
    assert not AC.is_checkable("notes.md") and AC.is_checkable("a/b.HTML")


def test_without_node_a_script_is_not_called_broken(ws, monkeypatch):
    monkeypatch.setattr(AC, "_node", lambda: "")
    (ws / "index.html").write_text(BROKEN_SCRIPT)
    assert AC.problems([ws / "index.html"]) == ""


# -- the turn-end check sends the model back ---------------------------------

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


def _write(content):
    return _Stream([_Chunk([_Choice(_Delta(tool_calls=[_TC(
        0, id="w1", name="write_file",
        arguments=json.dumps({"path": "index.html", "content": content}))]),
        finish="tool_calls")], usage=_Usage())])


def _say(text):
    return _Stream([_Chunk([_Choice(_Delta(content=text), finish="stop")],
                           usage=_Usage())])


@pytest.fixture
def client(monkeypatch, ws):
    from delfin.agent import model_capabilities as mc
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    monkeypatch.setenv("DELFIN_AUTO_VERIFY", "smart")
    from delfin.agent.api_client import KitToolPermissions, create_client
    cl = create_client(backend="api", provider="ollama",
                       model="qwen2.5-coder:7b", cwd=str(ws))
    cl.set_permissions(KitToolPermissions(
        workspace=ws, mode="acceptEdits", task_session_id="sess-ac"))
    return cl


def _run(client, streams):
    sent = []

    def _create(**kw):
        sent.append(list(kw.get("messages") or []))
        return streams[min(len(sent) - 1, len(streams) - 1)]

    client.client.chat.completions.create = _create
    events = list(client.stream_message(
        "sys", [{"role": "user", "content": "bau die Seite"}], max_tokens=100))
    notices = "".join(getattr(e, "text", "") or "" for e in events
                      if getattr(e, "type", "") == "notice")
    return sent, notices


def test_a_broken_page_written_by_the_agent_is_sent_back(client):
    sent, notices = _run(client, [_write(BROKEN_SCRIPT), _say("Fertig."),
                                  _say("Behoben.")])
    assert "Auto-verify" in notices
    last = " ".join(str(m.get("content") or "") for m in sent[-1])
    assert "inline <script> at line 3" in last


def test_a_sound_page_ends_the_turn_without_a_skip_notice(client):
    sent, notices = _run(client, [_write(GOOD), _say("Fertig.")])
    assert len(sent) == 2
    assert "Auto-verify" not in notices and "NOT machine-verified" not in notices


def test_a_module_page_is_not_judged_by_a_file_url(ws, monkeypatch):
    """From file:// every module import is refused by the browser (CORS,
    origin null): that is the check's error, not the page's, and it came
    back on every machine with Playwright and Chromium."""
    called: list = []
    monkeypatch.setattr(AC, "_browser_errors", lambda p: called.append(p) or ["x"])
    (ws / "x.js").write_text("export default 1;\n")
    (ws / "m.html").write_text('<script type="module" src="./x.js"></script>')
    (ws / "plain.html").write_text("<p>hi</p>")
    assert AC.problems([ws / "m.html"]) == ""
    assert AC.problems([ws / "plain.html"]) != ""
    assert called == [ws / "plain.html"]
