"""What a line starting with punctuation means, and what it must not.

The routing order is the dashboard's, kept identical on purpose: a
`/deploy` that means one thing in the browser and another in a terminal is
worse than not having it at all.

Two rules earn their own tests.

A PASTED PATH IS NOT A COMMAND. `/home/user/notes.md` has to reach the
model as text.

`!cmd` GOES THROUGH THE GATE, not around it. A shell escape that skipped
the deny-list and the approval prompt would be a way to do, from the
agent's own prompt, exactly what the agent may not do — and it would look
like a feature.
"""

from __future__ import annotations

import io

from delfin.agent import repl, repl_commands as rc, repl_render as rr


PLAIN = rr.Theme(enabled=False)


class _Tty(io.StringIO):
    def isatty(self):
        return True


class _Engine:
    def __init__(self):
        self.messages = [{"role": "user", "content": "hi"}]
        self.token_usage = {"input": 0, "output": 0}
        self.kit_permissions = type("P", (), {"mode": "plan"})()
        self.client = type("C", (), {"model": "m"})()
        self.mode = "solo"
        self.session_id = "sid"

    def get_status(self):
        return {"input_tokens": 12, "output_tokens": 3, "cached_tokens": 0,
                "cost_usd": 0.5, "mode": "solo", "provider": "kit",
                "role": "solo_agent", "session_id": "sid"}


def _agent(tmp_path, engine=None):
    engine = engine or _Engine()
    out, err = io.StringIO(), _Tty()
    agent = repl.TerminalAgent(
        engine, repl.ReplOptions(cwd=tmp_path), out=out, err=err)
    agent.transcript.theme = PLAIN
    return agent, engine, err


# ---------------------------------------------------------------------------
# Command or path
# ---------------------------------------------------------------------------

def test_a_pasted_path_is_not_a_command():
    assert rc.looks_like_command("/home/user/notes.md") is False
    assert rc.looks_like_command("/etc/hosts") is False
    assert rc.looks_like_command("/.hidden") is False
    assert rc.looks_like_command("/") is False
    assert rc.looks_like_command("") is False


def test_a_command_is_a_command():
    assert rc.looks_like_command("/help") is True
    assert rc.looks_like_command("/session search foo") is True


def test_the_two_surfaces_agree_on_what_routes_where():
    """A shared name must not mean two different things."""
    from delfin.dashboard import tab_agent

    for name in sorted(rc.BUILTINS):
        ours = rc.looks_like_command(name)
        assert ours is True, f"{name} does not read as a command"
        if tab_agent.slash_command_routes_to_builtin(name):
            continue                      # the dashboard has it too: fine
    assert tab_agent.slash_command_routes_to_builtin("/home/u/x.py") is False


# ---------------------------------------------------------------------------
# The four-step route
# ---------------------------------------------------------------------------

def test_a_builtin_is_answered_here_and_not_sent_to_the_model():
    result = rc.dispatch("/status", rc.ReplContext(engine=_Engine()))
    assert result.handled is True
    assert result.prompt == ""
    assert "tokens" in result.output


def test_a_user_command_is_expanded_and_sent(tmp_path):
    cmds = tmp_path / ".delfin" / "commands"
    cmds.mkdir(parents=True)
    (cmds / "deploy.md").write_text("Deploy to $ARGUMENTS and report back.")

    result = rc.dispatch("/deploy staging", rc.ReplContext(workspace=tmp_path))
    assert result.handled is True
    assert result.prompt == "Deploy to staging and report back."


def test_a_subagent_command_is_expanded(tmp_path):
    result = rc.dispatch("/explore where is the parser",
                         rc.ReplContext(workspace=tmp_path))
    assert result.handled is True
    assert result.prompt, "the builtin subagent commands must still route"


def test_an_unknown_slash_reaches_the_model_as_text(tmp_path):
    result = rc.dispatch("/zzz whatever", rc.ReplContext(workspace=tmp_path))
    assert result.handled is False, (
        "refusing here would make a typo unanswerable")


def test_a_handler_that_raises_does_not_end_the_session():
    def _boom(ctx, args):
        raise RuntimeError("handler is broken")

    original = rc.BUILTINS["/status"]
    rc.BUILTINS["/status"] = rc.ReplCommand(
        "/status", "session", "x", _boom)
    try:
        result = rc.dispatch("/status", rc.ReplContext())
        assert "failed" in result.output
    finally:
        rc.BUILTINS["/status"] = original


# ---------------------------------------------------------------------------
# Help cannot drift
# ---------------------------------------------------------------------------

def test_help_lists_every_command_the_table_accepts():
    text = rc.dispatch("/help", rc.ReplContext()).output
    for name in rc.BUILTINS:
        if name == "/quit":
            continue                      # an alias of /exit
        assert name in text, f"{name} is accepted but undocumented"


def test_help_is_generated_from_the_table_not_written_by_hand():
    import inspect
    src = inspect.getsource(rc._help)
    assert "help_gen" in src and "palette_rows" in src, (
        "a hand-written help page and a command table are one edit away "
        "from disagreeing; help_gen exists because that already happened")


def test_a_command_added_without_a_summary_is_visible_in_help():
    rc.BUILTINS["/temp"] = rc.ReplCommand(
        "/temp", "session", "a temporary thing", lambda c, a: rc.CommandResult())
    try:
        assert "/temp" in rc.dispatch("/help", rc.ReplContext()).output
    finally:
        del rc.BUILTINS["/temp"]


# ---------------------------------------------------------------------------
# @path
# ---------------------------------------------------------------------------

def test_at_paths_are_annotated_and_never_read(tmp_path):
    (tmp_path / "calc.py").write_text("secret contents")
    out = rc.expand_at_references("look at @calc.py please", tmp_path)
    assert "calc.py" in out
    assert "secret contents" not in out, (
        "reading it here would put file contents into the conversation "
        "without going through the gate that decides whether they may be")


def test_a_missing_at_path_says_so(tmp_path):
    out = rc.expand_at_references("check @nope.py", tmp_path)
    assert "not found" in out


def test_completion_stays_inside_the_workspace(tmp_path):
    (tmp_path / "src").mkdir()
    (tmp_path / "src" / "app.py").touch()
    (tmp_path / "readme.md").touch()

    assert "readme.md" in rc.complete_path("@read", tmp_path)
    assert "src/" in rc.complete_path("@", tmp_path)
    assert rc.complete_path("@../", tmp_path) == [], (
        "offering paths outside the workspace invites a request that is "
        "then refused")
    assert rc.complete_path("@/etc/", tmp_path) == []


def test_hidden_files_appear_only_when_asked_for(tmp_path):
    (tmp_path / ".env").touch()
    (tmp_path / "app.py").touch()
    assert ".env" not in rc.complete_path("@", tmp_path)
    assert ".env" in rc.complete_path("@.e", tmp_path)


# ---------------------------------------------------------------------------
# The prefixes, from the loop
# ---------------------------------------------------------------------------

def test_a_plain_line_goes_to_the_model(tmp_path):
    agent, _engine, _err = _agent(tmp_path)
    assert agent._handle_line("fix the failing test") == "fix the failing test"


def test_a_builtin_is_consumed_by_the_loop(tmp_path):
    agent, _engine, err = _agent(tmp_path)
    assert agent._handle_line("/status") == ""
    assert "tokens" in err.getvalue()


def test_clear_starts_a_fresh_conversation(tmp_path):
    agent, engine, _err = _agent(tmp_path)
    agent._handle_line("/clear")
    assert engine.messages == []


def test_exit_ends_the_loop(tmp_path):
    agent, _engine, _err = _agent(tmp_path)
    agent._handle_line("/exit")
    assert agent._quit is True


def test_a_shell_escape_goes_through_the_agents_own_gate(tmp_path):
    """Never around it.

    A `!cmd` that skipped the deny-list and the approval prompt would be a
    way to do, from the agent's prompt, exactly what the agent may not —
    and it would look like a feature rather than a hole.
    """
    ran: list[str] = []

    class _Gated(_Engine):
        def run_gated_bash(self, command):
            ran.append(command)
            return "total 0"

    agent, _engine, err = _agent(tmp_path, _Gated())
    assert agent._handle_line("!ls -la") == ""
    assert ran == ["ls -la"]
    assert "total 0" in err.getvalue()


def test_a_backend_with_no_gate_refuses_the_shell_escape(tmp_path):
    agent, _engine, err = _agent(tmp_path)      # no run_gated_bash
    assert agent._handle_line("!rm -rf /") == ""
    assert "same gate" in err.getvalue()


def test_shell_output_cannot_repaint_the_terminal(tmp_path):
    class _Gated(_Engine):
        def run_gated_bash(self, command):
            return "ok\x1b[2J\x1b]0;pwned\x07"

    agent, _engine, err = _agent(tmp_path, _Gated())
    agent._handle_line("!echo x")
    assert "\x1b" not in err.getvalue()


def test_a_hash_line_writes_a_memory_marked_as_the_users(tmp_path, monkeypatch):
    saved: dict = {}
    monkeypatch.setattr(
        "delfin.agent.memory_store.save_typed_memory",
        lambda **kw: saved.update(kw))
    agent, _engine, err = _agent(tmp_path)
    assert agent._handle_line("#the build needs python 3.11") == ""
    assert saved.get("text") == "the build needs python 3.11"
    assert saved.get("author") == "user", (
        "model-written and user-written memory must stay distinguishable")
    assert "remembered" in err.getvalue()


def test_the_shell_escape_reaches_the_real_gate(tmp_path, monkeypatch):
    """Not a fake: the engine method has to hand the command to the gate.

    Everything else in this file uses a stub, which would happily pass if
    the entry point did not exist at all — and it did not, until this
    test asked for it.
    """
    from unittest.mock import MagicMock, patch
    from delfin.agent import engine as E
    from delfin.agent.api_client import KitToolPermissions

    seen: list = []

    class _Executor:
        def execute(self, name, arguments, permissions):
            seen.append((name, arguments, permissions))
            return "gate saw it"

    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = E.AgentEngine(repo_dir=tmp_path, backend="api", provider="kit",
                            model="m", mode="solo")
    perms = KitToolPermissions(workspace=tmp_path)
    monkeypatch.setattr(type(eng), "kit_permissions",
                        property(lambda self: perms))
    monkeypatch.setattr("delfin.agent.api_client._doc_executor", _Executor())

    assert eng.run_gated_bash("ls -la") == "gate saw it"
    assert seen == [("bash", {"command": "ls -la"}, perms)], (
        "the command must arrive as a normal bash tool call, gate and all")


def test_a_backend_without_permissions_has_no_shell_either(tmp_path):
    from unittest.mock import MagicMock, patch
    from delfin.agent import engine as E
    import pytest

    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        eng = E.AgentEngine(repo_dir=tmp_path, backend="api", provider="claude",
                            model="m", mode="solo")
    if eng.kit_permissions is None:
        with pytest.raises(RuntimeError):
            eng.run_gated_bash("ls")


def test_a_shell_escape_that_needs_approval_does_not_deadlock(tmp_path):
    """The asker must not be the answerer.

    run_gated_bash goes through the permission gate, and the gate parks
    its question for THE MAIN THREAD. Running it directly from the main
    thread is therefore a deadlock by construction: it waits for itself.
    Found by watching a live session stop responding, not by reading.
    """
    import threading
    from delfin.agent import terminal_confirm as tc

    broker = tc.TerminalConfirmBroker(timeout_s=5)
    answered = threading.Event()

    class _Gated(_Engine):
        def run_gated_bash(self, command):
            # Exactly what the gate does: ask, and block until answered.
            allowed = broker.callback("bash", {"command": command},
                                      f"$ {command}")
            answered.set()
            return "ran it" if allowed else "refused"

    agent, _engine, err = _agent(tmp_path, _Gated())
    agent.broker = broker
    agent._raw_for_prompt = lambda: _KeyFeed(["y"])

    done = threading.Event()
    threading.Thread(
        target=lambda: (agent._handle_line("!ls"), done.set()),
        daemon=True).start()

    assert done.wait(timeout=10), "the shell escape deadlocked"
    assert answered.is_set()
    assert "ran it" in err.getvalue()


class _KeyFeed:
    active = True

    def __init__(self, keys):
        self._keys = list(keys)

    def read_ready(self, timeout):
        return self._keys.pop(0) if self._keys else "n"


def test_the_shell_escape_shows_output_not_the_envelope(tmp_path):
    """The tool returns JSON because a MODEL needs it.

    A person who typed `!git status` wants the output of git status.
    """
    class _Gated(_Engine):
        def run_gated_bash(self, command):
            return ('{"exit_code": 0, "stdout": " M calc.py\\n", '
                    '"stderr": "", "command": "git status"}')

    agent, _engine, err = _agent(tmp_path, _Gated())
    agent._handle_line("!git status")
    text = err.getvalue()
    assert "M calc.py" in text
    assert "exit_code" not in text, "the envelope is for the model, not here"


def test_a_refused_shell_escape_says_refused(tmp_path):
    class _Gated(_Engine):
        def run_gated_bash(self, command):
            return '{"error": "command is on the deny-list"}'

    agent, _engine, err = _agent(tmp_path, _Gated())
    agent._handle_line("!rm -rf /")
    assert "refused" in err.getvalue()
    assert "deny-list" in err.getvalue()


def test_a_failing_command_reports_its_exit_code(tmp_path):
    class _Gated(_Engine):
        def run_gated_bash(self, command):
            return ('{"exit_code": 2, "stdout": "", "stderr": "no such file"}')

    agent, _engine, err = _agent(tmp_path, _Gated())
    agent._handle_line("!cat nope")
    assert "no such file" in err.getvalue()
    assert "exit 2" in err.getvalue()
