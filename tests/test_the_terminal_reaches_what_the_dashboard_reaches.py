"""The terminal reaches what the dashboard reaches.

The browser surface could list the memory store, stage and approve a
diff, restore a file from the undo journal, read the attention inbox and
enumerate the model's tools. The terminal could talk to the model and
little else, so a session started there could be discussed but not
driven.

Three rules earn a test per command, and they are the three ways this
kind of addition goes wrong.

IT IS IN THE TABLE. The table is data; a handler nobody registered is
dead code.

/help LISTS IT. Help is GENERATED from the table, so a command that is
accepted and undocumented means the generation was bypassed.

AN EMPTY STORE IS A SENTENCE. A command that reaches an optional module
must degrade to one line of text when the module is missing or the store
is empty — a traceback out of a command handler ends the line the user
was typing, and the store being empty is the normal case on day one.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from delfin.agent import repl_commands as rc


# Every command this file added to the terminal. Adding one here and
# forgetting the table (or the reverse) fails immediately below.
NEW_COMMANDS = (
    "/tools", "/usage", "/export",
    "/memories", "/forget",
    "/undo", "/undo-file",
    "/pending", "/approve", "/reject",
    "/git", "/agents", "/skills", "/hooks", "/bash", "/attention",
    "/plans", "/commands", "/trace",
)

# The module each command degrades through when it is missing. Used to
# simulate an install that does not ship it.
BACKING_MODULE = {
    "/tools": "delfin.agent.api_client",
    "/memories": "delfin.agent.memory_store",
    "/forget": "delfin.agent.memory_store",
    "/undo-file": "delfin.agent.change_journal",
    "/pending": "delfin.agent.pending_changes",
    "/approve": "delfin.agent.pending_changes",
    "/reject": "delfin.agent.pending_changes",
    "/agents": "delfin.agent.subagents",
    "/skills": "delfin.agent.skills",
    "/hooks": "delfin.agent.hooks_editor",
    "/bash": "delfin.agent.bash_jobs",
    "/attention": "delfin.agent.attention",
    "/plans": "delfin.agent.memory_store",
    "/commands": "delfin.agent.slash_commands",
    "/trace": "delfin.agent.tool_trace",
}


@pytest.fixture(autouse=True)
def _stores_under_tmp(tmp_path, monkeypatch):
    """Every store these commands read is redirected away from the user's.

    They are keyed by SESSION ID, not by workspace, so a tmp_path
    workspace does not isolate them: a test that stages a diff or
    journals an edit would leave it in the store a real session reads,
    and the next real /pending would show it.
    """
    from delfin.agent import (
        attention, audit_log, change_journal, pending_changes, tool_trace,
    )

    root = tmp_path / "state"
    monkeypatch.setattr(pending_changes, "_pending_root",
                        lambda: root / "pending")
    monkeypatch.setattr(change_journal, "_undo_root", lambda: root / "undo")
    monkeypatch.setattr(tool_trace, "_DIR", root / "tool_traces")
    monkeypatch.setattr(attention, "_inbox_path",
                        lambda: root / "attention.json")
    monkeypatch.setattr(audit_log, "_default_log_path",
                        lambda: root / "audit.log")


class _Engine:
    """The least an engine can be and still be asked these questions."""

    def __init__(self, **kw):
        self.messages: list = []
        self.client = type("C", (), {"model": "test-model"})()
        self.session_id = ""
        self.kit_permissions = None
        self.__dict__.update(kw)

    def get_status(self):
        return {"provider": "kit", "role": "solo_agent", "input_tokens": 0,
                "output_tokens": 0, "cached_tokens": 0, "cost_usd": 0.0}

    def trace_session(self):
        return self.session_id


def _ctx(tmp_path, engine=None):
    return rc.ReplContext(engine=engine or _Engine(), workspace=tmp_path)


# ---------------------------------------------------------------------------
# The three rules, per command
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", NEW_COMMANDS)
def test_the_command_is_in_the_table(name):
    cmd = rc.BUILTINS.get(name)
    assert cmd is not None, f"{name} is documented nowhere the router reads"
    assert callable(cmd.handler)
    assert cmd.summary.strip(), "a row with no summary renders a blank help line"


@pytest.mark.parametrize("name", NEW_COMMANDS)
def test_help_lists_the_command(name):
    """Generated, never hand-written: an unlisted command is unfindable."""
    assert name in rc.dispatch("/help", rc.ReplContext()).output


@pytest.mark.parametrize("name", NEW_COMMANDS)
def test_an_empty_store_is_a_sentence_not_a_traceback(tmp_path, name):
    result = rc.BUILTINS[name].handler(_ctx(tmp_path), "")
    assert isinstance(result, rc.CommandResult)
    assert result.output.strip(), f"{name} answered with nothing at all"


@pytest.mark.parametrize(
    "name", sorted(set(NEW_COMMANDS) & set(BACKING_MODULE)))
def test_a_missing_module_degrades_to_one_line(tmp_path, monkeypatch, name):
    """An optional module that is not installed is not a crash.

    ``None`` in sys.modules is exactly what an interpreter reports for a
    module that cannot be imported, so this is the real failure and not a
    mock of it.
    """
    monkeypatch.setitem(sys.modules, BACKING_MODULE[name], None)
    result = rc.BUILTINS[name].handler(_ctx(tmp_path), "")
    assert isinstance(result, rc.CommandResult)
    assert result.output.strip()


@pytest.mark.parametrize("name", NEW_COMMANDS)
def test_no_engine_at_all_is_survivable(tmp_path, name):
    """The first line of a session is typed before any engine has run."""
    result = rc.BUILTINS[name].handler(
        rc.ReplContext(engine=None, workspace=tmp_path), "")
    assert result.output.strip()


# ---------------------------------------------------------------------------
# The table stays data
# ---------------------------------------------------------------------------

def test_every_command_is_reached_through_the_table_not_an_if_chain():
    """dispatch looks the name up; it does not know any of them by name."""
    import inspect

    src = inspect.getsource(rc.dispatch)
    for name in NEW_COMMANDS:
        assert name not in src, (
            f"{name} is special-cased in dispatch; the table is the router")


def test_help_counts_every_row_the_table_holds():
    text = rc.dispatch("/help", rc.ReplContext()).output
    for name in rc.BUILTINS:
        if name == "/quit":
            continue                      # an alias of /exit
        assert name in text, f"{name} is accepted but undocumented"


# ---------------------------------------------------------------------------
# Reading is reading
# ---------------------------------------------------------------------------

def test_tools_enumerates_the_catalogue_the_model_is_offered(tmp_path):
    """From api_client's own catalogue, not from a list kept here.

    A second copy of the tool list is a copy that goes stale the next
    time a tool is added.
    """
    from delfin.agent import api_client

    out = rc.BUILTINS["/tools"].handler(_ctx(tmp_path), "").output
    assert "read_file" in out
    known = {t.get("function", {}).get("name")
             for t in api_client._DOC_TOOLS_OPENAI}
    assert any(n and n in out for n in known)


def test_tools_can_be_filtered(tmp_path):
    out = rc.BUILTINS["/tools"].handler(_ctx(tmp_path), "notebook").output
    assert "notebook" in out.lower()
    assert "read_file " not in out


def test_usage_names_the_rate_it_used(tmp_path):
    """A cost figure with no rate beside it cannot be checked."""
    out = rc.BUILTINS["/usage"].handler(_ctx(tmp_path), "").output
    assert "rate" in out
    assert "test-model" in out


def test_an_unpriced_model_says_so_instead_of_printing_zero(tmp_path):
    engine = _Engine()
    engine.client = type("C", (), {"model": "a-model-nobody-priced"})()
    out = rc.BUILTINS["/usage"].handler(_ctx(tmp_path, engine), "").output
    assert "$0.0000" not in out or "unmeasured" in out or "no per-token" in out


def test_memories_lists_what_the_store_holds(tmp_path):
    from delfin.agent import memory_store

    memory_store.save_typed_memory(
        "the build needs python 3.11", repo_root=tmp_path)
    out = rc.BUILTINS["/memories"].handler(_ctx(tmp_path), "").output
    assert "python 3.11" in out or "build" in out


def test_hooks_is_read_only(tmp_path):
    """Listing a hook must not be able to write one.

    Editing settings the user owns is not something a listing does, and
    the terminal has no confirmation surface to do it behind.
    """
    import inspect

    src = inspect.getsource(rc._hooks)
    for verb in ("add_hook", "remove_hook", "_write_settings"):
        assert verb not in src, f"the hooks listing reaches {verb}"


def test_attention_is_read_only(tmp_path):
    import inspect

    src = inspect.getsource(rc._attention)
    for verb in ("answer_item(", "dismiss_item(", "clear_all(", "resolve("):
        assert verb not in src, f"the attention listing reaches {verb}"


def test_trace_reports_the_calls_of_this_session(tmp_path):
    from delfin.agent import tool_trace

    tool_trace.record("sid-1", tool="read_file", tool_input={"path": "x.py"},
                      output="ok")
    engine = _Engine(session_id="sid-1")
    out = rc.BUILTINS["/trace"].handler(_ctx(tmp_path, engine), "").output
    assert "read_file" in out


# ---------------------------------------------------------------------------
# /export
# ---------------------------------------------------------------------------

def test_export_writes_the_conversation_as_markdown(tmp_path, monkeypatch):
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    engine = _Engine()
    engine.messages = [{"role": "user", "content": "why is it slow"},
                       {"role": "assistant", "content": "the loop re-reads"}]
    out = rc.BUILTINS["/export"].handler(_ctx(tmp_path, engine), "").output
    written = list((tmp_path / ".delfin" / "exports").glob("*.md"))
    assert len(written) == 1, out
    text = written[0].read_text(encoding="utf-8")
    assert "why is it slow" in text and "the loop re-reads" in text


def test_export_reads_content_blocks_not_their_repr(tmp_path, monkeypatch):
    """Some backends store a list of blocks where others store a string.

    Writing the list out verbatim would put JSON in the file where the
    user expects their conversation.
    """
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    engine = _Engine()
    engine.messages = [{"role": "assistant",
                        "content": [{"type": "text", "text": "the answer"}]}]
    rc.BUILTINS["/export"].handler(_ctx(tmp_path, engine), "")
    text = list((tmp_path / ".delfin" / "exports").glob("*.md"))[0].read_text()
    assert "the answer" in text
    assert "'type'" not in text


# ---------------------------------------------------------------------------
# /undo drops context. /undo-file restores files. They are not the same.
# ---------------------------------------------------------------------------

def test_undo_drops_the_whole_last_turn(tmp_path):
    engine = _Engine()
    engine.messages = [
        {"role": "user", "content": "first"},
        {"role": "assistant", "content": "answer"},
        {"role": "user", "content": "second"},
        {"role": "assistant", "content": "tool call"},
        {"role": "tool", "content": "tool result"},
        {"role": "assistant", "content": "final"},
    ]
    out = rc.BUILTINS["/undo"].handler(_ctx(tmp_path, engine), "").output
    assert [m["content"] for m in engine.messages] == ["first", "answer"], (
        "a turn is the user message and everything the agent did after it")
    assert "no file" in out, "the two undos must not be mistaken for each other"


def test_undo_with_no_user_turn_does_not_empty_the_context(tmp_path):
    """Emptying it is /clear's job, and nobody asked for that."""
    engine = _Engine()
    engine.messages = [{"role": "system", "content": "you are an agent"}]
    rc.BUILTINS["/undo"].handler(_ctx(tmp_path, engine), "")
    assert len(engine.messages) == 1


def _record_edit(session_id, target, before, after):
    """Journal one agent edit exactly as the executor does: write, then record.

    record_change hashes the file the write produced, so the order is not
    cosmetic — recording first would store a post-hash that never
    existed and the revert guard would call the undo a conflict.
    """
    from delfin.agent import change_journal

    target.write_text(before)
    target.write_text(after)
    return change_journal.record_change(
        session_id, tool="edit_file", path=target,
        old_text=before, new_text=after)


def test_undo_file_restores_the_pre_image_from_the_journal(tmp_path):
    target = tmp_path / "calc.py"
    _record_edit("sid-undo", target, "original\n", "rewritten\n")

    engine = _Engine(session_id="sid-undo")
    out = rc.BUILTINS["/undo-file"].handler(_ctx(tmp_path, engine), "last").output
    assert target.read_text() == "original\n", out
    assert "restored" in out


def test_undo_file_lists_before_it_touches_anything(tmp_path):
    target = tmp_path / "calc.py"
    _record_edit("sid-list", target, "original\n", "new\n")

    engine = _Engine(session_id="sid-list")
    out = rc.BUILTINS["/undo-file"].handler(_ctx(tmp_path, engine), "").output
    assert "calc.py" in out
    assert target.read_text() == "new\n", "a listing must not restore anything"


def test_undo_file_refuses_a_scope_it_does_not_know(tmp_path):
    engine = _Engine(session_id="sid-x")
    out = rc.BUILTINS["/undo-file"].handler(
        _ctx(tmp_path, engine), "everything").output
    assert "unknown scope" in out


# ---------------------------------------------------------------------------
# Staging: approve and reject act only on an id they found
# ---------------------------------------------------------------------------

def _stage_one(tmp_path, session_id="sid-pending"):
    from delfin.agent import pending_changes

    target = tmp_path / "notes.md"
    target.write_text("before\n")
    return pending_changes.stage(
        session_id, tool="edit_file", path=target,
        old_text="before\n", new_text="after\n")


def test_pending_lists_the_staged_diff(tmp_path):
    rec = _stage_one(tmp_path)
    engine = _Engine(session_id="sid-pending")
    out = rc.BUILTINS["/pending"].handler(_ctx(tmp_path, engine), "").output
    assert "notes.md" in out
    assert str(rec.get("id")) in out


def test_approve_applies_the_named_change(tmp_path):
    _stage_one(tmp_path, "sid-approve")
    engine = _Engine(session_id="sid-approve")
    out = rc.BUILTINS["/approve"].handler(_ctx(tmp_path, engine), "1").output
    assert (tmp_path / "notes.md").read_text() == "after\n", out
    assert "applied" in out


def test_reject_discards_without_writing(tmp_path):
    _stage_one(tmp_path, "sid-reject")
    engine = _Engine(session_id="sid-reject")
    out = rc.BUILTINS["/reject"].handler(_ctx(tmp_path, engine), "1").output
    assert (tmp_path / "notes.md").read_text() == "before\n", out
    assert "reject" in out.lower() or "1" in out


def test_approving_an_id_that_does_not_exist_claims_nothing(tmp_path):
    """A queue that reports work it did not do is worse than one that refuses."""
    _stage_one(tmp_path, "sid-missing")
    engine = _Engine(session_id="sid-missing")
    out = rc.BUILTINS["/approve"].handler(_ctx(tmp_path, engine), "999").output
    assert "applied" not in out
    assert (tmp_path / "notes.md").read_text() == "before\n"


def test_rejecting_an_id_that_does_not_exist_claims_nothing(tmp_path):
    _stage_one(tmp_path, "sid-missing-2")
    engine = _Engine(session_id="sid-missing-2")
    out = rc.BUILTINS["/reject"].handler(_ctx(tmp_path, engine), "999").output
    assert "rejected 1" not in out


def test_forgetting_a_memory_that_is_not_there_deletes_nothing(tmp_path):
    from delfin.agent import memory_store

    memory_store.save_typed_memory("keep this one", repo_root=tmp_path)
    before = sorted(p.name for p in tmp_path.rglob("*.md"))
    out = rc.BUILTINS["/forget"].handler(_ctx(tmp_path), "no-such-memory").output
    assert "no memory named" in out
    assert sorted(p.name for p in tmp_path.rglob("*.md")) == before


def test_forgetting_a_memory_names_what_it_deleted(tmp_path):
    from delfin.agent import memory_store

    _path, slug, _type = memory_store.save_typed_memory(
        "the parser is in code_nav", repo_root=tmp_path)
    out = rc.BUILTINS["/forget"].handler(_ctx(tmp_path), slug).output
    assert "deleted" in out and slug in out


# ---------------------------------------------------------------------------
# /git: through the gate, and read-only
# ---------------------------------------------------------------------------

class _Gated(_Engine):
    def __init__(self, **kw):
        super().__init__(**kw)
        self.ran: list = []
        self.kit_permissions = type(
            "P", (), {"confirm_callback": None,
                      "matches_bash_auto_allow": lambda self, cmd: True})()

    def run_gated_bash(self, command):
        self.ran.append(command)
        return '{"exit_code": 0, "stdout": " M calc.py\\n", "stderr": ""}'


def test_git_goes_through_the_agents_own_gate(tmp_path):
    """Never around it.

    A /git that shelled out directly would be a way to run a command from
    the agent's prompt without the deny-list or the approval prompt, and
    it would look like a convenience.
    """
    engine = _Gated()
    out = rc.BUILTINS["/git"].handler(_ctx(tmp_path, engine), "status").output
    assert engine.ran == ["git status --short"]
    assert "M calc.py" in out
    assert "exit_code" not in out, "the envelope is for the model, not here"


def test_git_never_runs_a_subcommand_that_is_not_read_only(tmp_path):
    engine = _Gated()
    for attempt in ("push", "commit -m x", "reset --hard", "status; rm -rf /"):
        out = rc.BUILTINS["/git"].handler(_ctx(tmp_path, engine), attempt).output
        assert "usage:" in out, attempt
    assert engine.ran == [], "a /git pass-through is a second shell escape"


def test_bash_kill_acts_only_on_a_job_id_it_was_given(tmp_path):
    """Killing "the last one" would reach a job the user is not looking at."""
    out = rc.BUILTINS["/bash"].handler(_ctx(tmp_path), "kill").output
    assert "usage:" in out
    out = rc.BUILTINS["/bash"].handler(_ctx(tmp_path), "kill nope-42").output
    assert "nope-42" in out and "unknown" in out


@pytest.mark.parametrize("name", ["/skills", "/commands", "/plans"])
def test_a_name_that_is_not_there_is_not_invented(tmp_path, name):
    out = rc.BUILTINS[name].handler(_ctx(tmp_path), "no-such-thing").output
    assert "no-such-thing" in out
    assert "no " in out.lower()


def test_git_without_a_gate_says_there_is_none(tmp_path):
    out = rc.BUILTINS["/git"].handler(_ctx(tmp_path), "status").output
    assert "gate" in out


def test_git_does_not_park_a_question_it_cannot_answer(tmp_path):
    """The gate parks approvals for the MAIN thread, which is this one.

    Calling into the prompt from a command handler is the asker waiting
    for the answerer, i.e. for itself. Found by reading the same deadlock
    the shell escape already documents.
    """
    engine = _Gated()
    engine.kit_permissions = type(
        "P", (), {"confirm_callback": lambda *a: True,
                  "matches_bash_auto_allow": lambda self, cmd: False})()
    out = rc.BUILTINS["/git"].handler(_ctx(tmp_path, engine), "status").output
    assert engine.ran == []
    assert "!git status --short" in out


# ---------------------------------------------------------------------------
# /trust names a remedy that exists here
# ---------------------------------------------------------------------------

def test_trust_no_longer_points_at_an_editor_the_terminal_cannot_open(tmp_path):
    """The old text named "the hooks / MCP editors" and nothing here
    opened either, so the only advice it gave was unreachable."""
    ws = tmp_path / "repo"
    (ws / ".delfin").mkdir(parents=True)
    (ws / ".delfin" / "settings.json").write_text(
        '{"hooks": {"PreToolUse": [{"matcher": "Edit",'
        ' "hooks": [{"type": "command", "command": "echo hi"}]}]}}')
    out = rc.BUILTINS["/trust"].handler(rc.ReplContext(workspace=ws), "").output
    assert "editors" not in out
    assert "/trust grant" in out


def test_trust_grants_through_the_same_wrapper_the_dashboard_uses(tmp_path,
                                                                  monkeypatch):
    """Trust has exactly two grant call sites and this adds no third."""
    from delfin.agent import hooks_editor

    seen: list = []
    monkeypatch.setattr(hooks_editor, "trust_this_workspace",
                        lambda ws: seen.append(ws) or {"kinds": {}})
    out = rc.BUILTINS["/trust"].handler(
        rc.ReplContext(workspace=tmp_path), "grant hooks").output
    assert seen == [tmp_path]
    assert "trusted" in out


def test_trust_revokes_through_the_same_wrapper(tmp_path, monkeypatch):
    from delfin.agent import mcp_editor

    seen: list = []
    monkeypatch.setattr(mcp_editor, "untrust_this_workspace",
                        lambda ws: bool(seen.append(ws)) or True)
    out = rc.BUILTINS["/trust"].handler(
        rc.ReplContext(workspace=tmp_path), "revoke mcp").output
    assert seen == [tmp_path]
    assert "withdrawn" in out


def test_a_typo_in_the_trust_subcommand_grants_nothing(tmp_path, monkeypatch):
    from delfin.agent import hooks_editor

    monkeypatch.setattr(hooks_editor, "trust_this_workspace",
                        lambda ws: pytest.fail("a typo granted trust"))
    for typo in ("grnt hooks", "trust hooks", "hooks"):
        out = rc.BUILTINS["/trust"].handler(
            rc.ReplContext(workspace=tmp_path), typo).output
        assert "usage:" in out


def test_the_grant_is_still_only_reachable_through_the_two_wrappers():
    """The chokepoint test reads the source tree; keep this file out of it."""
    import inspect

    src = inspect.getsource(rc)
    assert "trust_workspace(" not in src
    assert "revoke_workspace(" not in src


# ---------------------------------------------------------------------------
# A command that never once did what it says
# ---------------------------------------------------------------------------

def test_mcp_lists_the_servers_that_would_load(tmp_path):
    """Found by running it: `/mcp` had never listed a server.

    `effective_servers` takes a WORKSPACE and was handed a registry
    object, so every call raised a TypeError — which the handler's broad
    except turned into "MCP registry unavailable", a sentence that reads
    like a diagnosis of the environment and was a diagnosis of the line
    above it. The built-in servers are always configured, so an empty
    listing here can only mean the call failed.
    """
    from delfin.agent import repl_commands as rc

    class _Ctx:
        workspace = tmp_path
        engine = None

    out = rc.BUILTINS["/mcp"].handler(_Ctx(), "").output
    assert "unavailable" not in out, out
    assert "delfin-tools" in out, "the built-ins are always configured"
    assert "built-in" in out, "and each row says where it came from"


def test_mcp_renders_rows_rather_than_printing_them(tmp_path):
    """`effective_servers` returns dicts. The old code interpolated each
    one straight into a line, so even with the right argument the output
    would have been a repr per server."""
    from delfin.agent import repl_commands as rc

    class _Ctx:
        workspace = tmp_path
        engine = None

    out = rc.BUILTINS["/mcp"].handler(_Ctx(), "").output
    assert "{'name'" not in out and "'command':" not in out


def test_mcp_still_answers_when_the_registry_cannot_be_read(monkeypatch,
                                                            tmp_path):
    from delfin.agent import mcp_client
    from delfin.agent import repl_commands as rc

    def _boom(_ws):
        raise RuntimeError("config is corrupt")

    monkeypatch.setattr(mcp_client, "effective_servers", _boom)

    class _Ctx:
        workspace = tmp_path
        engine = None

    assert "unavailable" in rc.BUILTINS["/mcp"].handler(_Ctx(), "").output
