"""Filing a bug report from the terminal and from the command line.

Input: a live terminal engine (``/bug``) or a saved session
(``delfin-agent bug report``). Output: the same report directory the
dashboard button writes, in the archive ``resolve_archive_dir`` names, with
a ``surface`` field saying which door it came through.

Every store is redirected: HOME, the archive, the session store, the tool
trace and the credential file all live under ``tmp_path``.
"""

from __future__ import annotations

import hashlib
import io
import json
import string
from contextlib import redirect_stdout
from pathlib import Path

import pytest

from delfin.agent import bug_report as br
from delfin.agent import repl_commands as rc
from delfin.agent import session_store as ss
from delfin.agent import tool_trace as tt


def _token(seed: str, n: int = 40) -> str:
    """High-entropy and prefix-free: the shape checks cannot see it, only
    the exact-value pass over the credentials DELFIN holds can."""
    out: list[str] = []
    h = seed.encode()
    alphabet = string.ascii_letters + string.digits
    while len(out) < n:
        h = hashlib.sha256(h).digest()
        out += [alphabet[b % len(alphabet)] for b in h]
    return "".join(out[:n])


STORED_KEY = _token("stored-credential")
ENV_KEY = _token("env-credential")


@pytest.fixture
def stores(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", lambda: home)
    archive = tmp_path / "archive"
    monkeypatch.setenv("DELFIN_BUG_ARCHIVE", str(archive))
    monkeypatch.setattr(br, "_FALLBACK_DIR", home / ".delfin" / "agent_bugs")
    sessions = tmp_path / "sessions"
    sessions.mkdir()
    monkeypatch.setattr(ss, "_SESSIONS_DIR", sessions)
    monkeypatch.setattr(tt, "_DIR", tmp_path / "tool_traces")
    ws = tmp_path / "ws"
    ws.mkdir()
    # One credential in DELFIN's own store (the conftest points the file
    # store into tmp and disables the OS keyring), one in the environment.
    from delfin.agent import credentials as cred
    cred.set_credential("KIT_TOOLBOX_API_KEY", STORED_KEY, store="file")
    monkeypatch.setenv("OPENAI_API_KEY", ENV_KEY)
    return {"archive": archive, "ws": ws, "home": home}


class _Client:
    model = "fixture-model"
    effort = "high"


class _Perms:
    def __init__(self, ws):
        self.workspace = str(ws)
        self.mode = "default"


class _Engine:
    """The attributes the terminal reads off a live engine."""

    def __init__(self, ws, sid="termsess1"):
        self.session_id = sid
        self.client = _Client()
        self.kit_permissions = _Perms(ws)
        self.repo_dir = str(ws)
        self.last_system_prompt = "system prompt of the last turn"
        self.last_compaction_info = {"compactions": 1}
        self.messages = [
            {"role": "user",
             "content": f"my key is {STORED_KEY} and {ENV_KEY}; read notes.txt"},
            {"role": "assistant", "content": [
                {"type": "text", "text": "Reading it."},
                {"type": "tool_use", "id": "t1", "name": "read_file",
                 "input": {"path": "notes.txt"}}]},
            {"role": "user", "content": [
                {"type": "tool_result", "tool_use_id": "t1",
                 "content": "notes"}]},
            {"role": "assistant", "content": "The file says notes."},
        ]

    def trace_session(self):
        return self.session_id

    def get_status(self):
        return {"mode": "solo", "provider": "kit", "backend": "api",
                "role": "solo", "session_id": self.session_id,
                "input_tokens": 120, "output_tokens": 30, "cost_usd": 0.0}


def _trace(sid):
    tt.record(sid, tool="read_file", tool_input={"path": "notes.txt"},
              output=f"notes with {STORED_KEY}", ok=True)
    tt.record(sid, tool="bash", tool_input={"command": "rm -rf build"},
              output="", ok=False,
              error="Command not on the auto-allow list: rm -rf build")


def _ctx(ws, engine, answers, said):
    asked: list[str] = []

    def ask(prompt):
        asked.append(prompt)
        return answers.pop(0)

    ctx = rc.ReplContext(engine=engine, workspace=ws, session_id=engine.session_id,
                         ask=ask, say=said.append)
    return ctx, asked


def _reports(archive: Path) -> list[Path]:
    if not archive.is_dir():
        return []
    return sorted(d for d in archive.iterdir() if (d / "report.json").is_file())


def _all_text(report_dir: Path) -> str:
    return "".join(p.read_text(encoding="utf-8", errors="replace")
                   for p in report_dir.rglob("*") if p.is_file())


# ---------------------------------------------------------------------------
# /bug in the terminal
# ---------------------------------------------------------------------------

def test_bug_writes_the_full_report_with_every_held_secret_redacted(stores):
    engine = _Engine(stores["ws"])
    _trace(engine.session_id)
    said: list[str] = []
    ctx, asked = _ctx(stores["ws"], engine, ["y"], said)

    res = rc.dispatch("/bug the answer ignored the file", ctx)

    reports = _reports(stores["archive"])
    assert len(reports) == 1, res.output
    d = reports[0]
    assert str(d) in res.output
    rep = json.loads((d / "report.json").read_text(encoding="utf-8"))
    assert rep["surface"] == "terminal"
    assert rep["description"] == "the answer ignored the file"
    assert rep["session_id"] == "termsess1"
    assert rep["model"] == "fixture-model"
    assert rep["provider"] == "kit"
    assert rep["mode"] == "solo"
    assert rep["effort"] == "high"
    assert rep["perms"] == "default"
    assert rep["workspace"] == str(stores["ws"])
    assert rep["input_tokens"] == 120
    assert "git_head" in rep
    assert rep["system_prompt"] == "system prompt of the last turn"
    assert rep["last_compaction_info"] == {"compactions": 1}
    assert len(rep["chat_messages"]) == 3
    assert len(rep["engine_messages"]) == 4
    assert [e["tool"] for e in rep["tool_trace"]] == ["read_file", "bash"]
    assert any("rm -rf build" in c for c in rep["denied_commands"])
    assert "notes.txt" in json.dumps(rep["referenced_files"])
    assert rep["settings"]["agent"] is not None
    # The preview came before the question, and named what goes in.
    preview = "\n".join(said)
    assert "4 engine message(s)" in preview
    assert "2 tool call(s)" in preview
    assert "1 refused command(s)" in preview
    assert str(stores["archive"]) in preview
    assert len(asked) == 1 and "y/N" in asked[0]
    # Neither the stored credential nor the environment key is anywhere
    # in what was written, although the conversation and the trace held
    # both.
    blob = _all_text(d)
    assert STORED_KEY not in blob
    assert ENV_KEY not in blob
    assert (d / "REDACTIONS.txt").is_file()


def test_answering_no_writes_nothing(stores):
    engine = _Engine(stores["ws"])
    said: list[str] = []
    ctx, _ = _ctx(stores["ws"], engine, ["n"], said)
    res = rc.dispatch("/bug it hung", ctx)
    assert _reports(stores["archive"]) == []
    assert "nothing was written" in res.output


def test_an_interrupted_question_writes_nothing(stores):
    engine = _Engine(stores["ws"])

    def ask(_prompt):
        raise EOFError

    ctx = rc.ReplContext(engine=engine, workspace=stores["ws"], ask=ask,
                         say=lambda _s: None)
    res = rc.dispatch("/bug it hung", ctx)
    assert _reports(stores["archive"]) == []
    assert "nothing was written" in res.output


def test_without_a_description_it_asks_for_one(stores):
    engine = _Engine(stores["ws"])
    said: list[str] = []
    ctx, asked = _ctx(stores["ws"], engine, ["it stalled after the tool", "y"],
                      said)
    rc.dispatch("/bug", ctx)
    assert len(asked) == 2
    rep = json.loads((_reports(stores["archive"])[0] / "report.json")
                     .read_text(encoding="utf-8"))
    assert rep["description"] == "it stalled after the tool"


def test_a_real_engine_yields_the_fields_the_report_needs(stores, monkeypatch):
    """The live collector reads attributes off a real AgentEngine, not
    only off the stand-in above."""
    from delfin.agent import bug_filing
    from delfin.agent import model_capabilities as mc
    from delfin.agent.engine import AgentEngine
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: mc.ModelCapabilities(
        model="m", provider="ollama", context_window=32768,
        supports_tools=True))
    eng = AgentEngine(repo_dir=str(stores["ws"]), backend="api",
                      provider="ollama", model="probe-model")
    eng.messages.append({"role": "user", "content": "hello"})
    tt.record(eng.trace_session(), tool="read_file",
              tool_input={"path": "a.txt"}, output="x", ok=True)
    got = bug_filing.collect_live(eng, stores["ws"])
    assert got["model"] == "probe-model"
    assert got["provider"] == "ollama"
    assert got["trace_session"] == eng.trace_session()
    assert got["chat_messages"] == [{"role": "user", "content": "hello"}]
    assert got["referenced_files"] == ["a.txt"]
    assert got["workspace"]


def test_an_empty_description_files_nothing(stores):
    engine = _Engine(stores["ws"])
    ctx, asked = _ctx(stores["ws"], engine, ["   "], [])
    res = rc.dispatch("/bug", ctx)
    assert len(asked) == 1
    assert _reports(stores["archive"]) == []
    assert "nothing was written" in res.output


def test_bug_is_a_listed_builtin():
    assert "/bug" in rc.BUILTINS
    assert any(row[1] == "/bug" for row in rc.palette_rows())


def test_the_terminal_hands_its_reader_to_the_command(stores):
    """Through the real TerminalAgent line handler: the question is read
    with the terminal's own reader and the preview is printed as chrome."""
    from delfin.agent import repl
    engine = _Engine(stores["ws"])
    answers = ["y"]
    prompts: list[str] = []

    def read_line(prompt):
        prompts.append(prompt)
        return answers.pop(0)

    out, err = io.StringIO(), io.StringIO()
    agent = repl.TerminalAgent(engine, repl.ReplOptions(cwd=stores["ws"]),
                               out=out, err=err, read_line=read_line)
    assert agent._handle_line("/bug the terminal path") == ""
    assert prompts == ["File this report? [y/N] "]
    screen = out.getvalue() + err.getvalue()
    assert "The bug report will contain" in screen
    assert "bug report written" in screen
    rep = json.loads((_reports(stores["archive"])[0] / "report.json")
                     .read_text(encoding="utf-8"))
    assert rep["surface"] == "terminal"


# ---------------------------------------------------------------------------
# delfin-agent bug report, for a session that already ended
# ---------------------------------------------------------------------------

def _saved_session(ws, sid="savedsess42"):
    ss.save_session(
        sid, mode="solo", provider="kit", model="saved-model", effort="low",
        perm_profile="acceptEdits", workspace=str(ws),
        chat_messages=[{"role": "user", "content": "plot the spectrum"},
                       {"role": "assistant", "content": "Done."}],
        engine_messages=[{"role": "user", "content": "plot the spectrum"},
                         {"role": "assistant", "content": "Done."}],
        token_usage={"input": 50, "output": 7}, cost_usd=0.01,
        last_compaction_info={"compactions": 0})
    ss.release_session_lock(sid)
    _trace(sid)
    return sid


def _cli(argv):
    from delfin.agent import cli
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc_ = cli.main(argv)
    return rc_, buf.getvalue()


def test_bug_report_files_the_latest_session_of_this_workspace(stores):
    sid = _saved_session(stores["ws"])
    code, out = _cli(["bug", "report", "--session", "latest",
                      "-m", "the plot was empty", "--yes",
                      "--workspace", str(stores["ws"])])
    assert code == 0, out
    reports = _reports(stores["archive"])
    assert len(reports) == 1
    assert str(reports[0]) in out
    rep = json.loads((reports[0] / "report.json").read_text(encoding="utf-8"))
    assert rep["surface"] == "cli"
    assert rep["session_id"] == sid
    assert rep["description"] == "the plot was empty"
    assert rep["model"] == "saved-model"
    assert rep["perms"] == "acceptEdits"
    assert rep["input_tokens"] == 50
    assert len(rep["chat_messages"]) == 2
    assert len(rep["tool_trace"]) == 2
    assert any("rm -rf build" in c for c in rep["denied_commands"])
    blob = _all_text(reports[0])
    assert STORED_KEY not in blob


def test_bug_report_takes_an_id_prefix(stores):
    sid = _saved_session(stores["ws"])
    code, out = _cli(["bug", "report", "--session", sid[:8], "-m", "x",
                      "--yes", "--workspace", str(stores["ws"])])
    assert code == 0, out
    assert len(_reports(stores["archive"])) == 1


def test_bug_report_without_a_session_says_so(stores):
    code, out = _cli(["bug", "report", "--session", "latest", "-m", "x",
                      "--yes", "--workspace", str(stores["ws"])])
    assert code != 0
    assert _reports(stores["archive"]) == []


def test_bug_report_without_yes_and_without_a_terminal_writes_nothing(
        stores, monkeypatch):
    _saved_session(stores["ws"])
    monkeypatch.setattr("sys.stdin", io.StringIO(""))
    code, out = _cli(["bug", "report", "-m", "x",
                      "--workspace", str(stores["ws"])])
    assert code != 0
    assert _reports(stores["archive"]) == []


# ---------------------------------------------------------------------------
# The surface field, and the listing that reads it
# ---------------------------------------------------------------------------

def test_the_dashboard_button_records_its_surface():
    src = (Path(__file__).resolve().parent.parent / "delfin" / "dashboard"
           / "tab_agent.py").read_text(encoding="utf-8")
    i = src.index("def _on_bug_report")
    assert 'surface="dashboard"' in src[i:i + 6000]


def test_a_report_without_the_field_reads_as_dashboard(tmp_path):
    d = br.write_bug_report(chat_messages=[], session_id="old",
                            archive_dir=tmp_path)
    rep = json.loads((d / "report.json").read_text(encoding="utf-8"))
    rep.pop("surface", None)
    (d / "report.json").write_text(json.dumps(rep), encoding="utf-8")
    assert br.report_surface(rep) == "dashboard"
    assert br.list_reports(tmp_path)[0]["surface"] == "dashboard"


def test_bug_ls_lists_a_terminal_report_with_its_surface(stores):
    _saved_session(stores["ws"])
    _cli(["bug", "report", "-m", "from the cli", "--yes",
          "--workspace", str(stores["ws"])])
    code, out = _cli(["bug", "ls", "--archive", str(stores["archive"])])
    assert code == 0
    name = _reports(stores["archive"])[0].name
    assert name in out
    assert "cli" in out
    assert "from the cli" in out


def test_triage_input_names_the_surface(stores):
    from delfin.agent import bug_watcher as bw
    _saved_session(stores["ws"])
    _cli(["bug", "report", "-m", "x", "--yes", "--workspace", str(stores["ws"])])
    d = _reports(stores["archive"])[0]
    body = bw._report_body(d, br.load_report(d))
    assert "surface" in body and "cli" in body


def test_the_writer_itself_removes_a_held_credential_without_a_prefix(stores):
    """Shared by all three surfaces: the redactor in front of every report
    file matched shapes only, and a key with no recognisable prefix (a KIT
    key) passed it verbatim into the shared archive."""
    d = br.write_bug_report(
        chat_messages=[{"role": "user", "content": f"key {STORED_KEY}"}],
        engine_messages=[{"role": "user", "content": f"env {ENV_KEY}"}],
        session_id="direct", archive_dir=stores["archive"])
    blob = _all_text(d)
    assert STORED_KEY not in blob
    assert ENV_KEY not in blob
