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
    '/tools',
    '/usage',
    '/export',
    '/memories',
    '/forget',
    '/agents',
    '/skills',
    '/hooks',
    '/attention',
    '/plans',
    '/commands',
    '/trace',
)

# The module each command degrades through when it is missing. Used to
# simulate an install that does not ship it.
BACKING_MODULE = {
    '/tools': 'delfin.agent.api_client',
    '/memories': 'delfin.agent.memory_store',
    '/forget': 'delfin.agent.memory_store',
    '/agents': 'delfin.agent.subagents',
    '/skills': 'delfin.agent.skills',
    '/hooks': 'delfin.agent.hooks_editor',
    '/attention': 'delfin.agent.attention',
    '/plans': 'delfin.agent.memory_store',
    '/commands': 'delfin.agent.slash_commands',
    '/trace': 'delfin.agent.tool_trace',
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


@pytest.mark.parametrize("name", ["/skills", "/commands", "/plans"])
def test_a_name_that_is_not_there_is_not_invented(tmp_path, name):
    out = rc.BUILTINS[name].handler(_ctx(tmp_path), "no-such-thing").output
    assert "no-such-thing" in out
    assert "no " in out.lower()
