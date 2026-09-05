"""Without a broker the gate does not refuse — it allows.

That is the whole reason this layer exists. `_gate_write_path` allows a
write inside the workspace when `confirm_callback is None`, and
`exit_plan_mode` approves itself. Both are right for a scheduled run and
both make an interactive terminal agent an unattended editor.

Three properties are load-bearing and each has its own test.

THE CALLBACK IS A BOUND METHOD. `_confirm_timed_out` reads
`perms.confirm_callback.__self__.last_timed_out` to tell an expired dialog
from a refusal. A lambda wrapper destroys that quietly, and every timeout
then lands in `denied_actions` as a real denial — permanently closing
paths the user never saw.

ONE READER, MANY ASKERS. Requests come from the turn worker, from
subagents, from background jobs. Exactly one thread renders and answers.

NOTHING GETS AN "ALWAYS" IT SHOULD NOT HAVE. An outside-workspace read
has no persistable form except a WRITABLE directory grant, so one
keystroke on a read prompt must not be able to hand over write access.
"""

from __future__ import annotations

import io
import threading
import time

from delfin.agent import repl, repl_render as rr
from delfin.agent import terminal_confirm as tc


PLAIN = rr.Theme(enabled=False)


class _Tty(io.StringIO):
    def isatty(self):
        return True


class _Keys:
    """Stands in for the raw-mode reader: hands over scripted keystrokes."""

    active = True

    def __init__(self, keys):
        self._keys = list(keys)

    def read_ready(self, timeout):
        return self._keys.pop(0) if self._keys else "n"


class _Engine:
    def __init__(self):
        self.kit_permissions = type("P", (), {"mode": "plan"})()
        self.token_usage = {"input": 0, "output": 0}
        self.client = None
        self.stopped = False

    def get_status(self):
        return {"input_tokens": 0, "output_tokens": 0, "cost_usd": 0.0}

    def set_kit_permission_mode(self, mode):
        self.kit_permissions.mode = mode

    def request_stop(self):
        self.stopped = True


def _agent(broker, engine=None):
    engine = engine or _Engine()
    out, err = io.StringIO(), _Tty()
    agent = repl.TerminalAgent(engine, out=out, err=err, broker=broker)
    agent.transcript.theme = PLAIN
    return agent, engine, err


def _write_req(path="src/app.py"):
    return tc.ConfirmRequest(kind=tc.CONFIRM, tool="write_file",
                             args={"path": path},
                             preview=f"--- a/{path}\n+++ b/{path}\n+x = 1\n")


def _bash_req(cmd="pytest -x tests/"):
    return tc.ConfirmRequest(kind=tc.CONFIRM, tool="bash",
                             args={"command": cmd}, preview=f"$ {cmd}")


# ---------------------------------------------------------------------------
# The bound method
# ---------------------------------------------------------------------------

def test_the_callback_is_bound_so_a_timeout_is_not_a_refusal():
    """Fails the moment someone wraps it in a lambda."""
    from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

    broker = tc.TerminalConfirmBroker(timeout_s=0.01)
    assert hasattr(broker.callback, "__self__"), (
        "the gate reads last_timed_out off __self__")

    perms = KitToolPermissions(workspace="/tmp", confirm_callback=broker.callback)
    assert _DocToolExecutor._confirm_timed_out(perms) is False

    assert broker.callback("bash", {"command": "ls"}, "$ ls") is False
    assert broker.last_timed_out is True
    assert _DocToolExecutor._confirm_timed_out(perms) is True, (
        "an expired dialog must not be recorded as a denial")


def test_an_answered_request_is_not_a_timeout():
    broker = tc.TerminalConfirmBroker(timeout_s=5)
    req_box: list = []

    def _ask():
        req_box.append(broker.callback("bash", {"command": "ls"}, "$ ls"))

    t = threading.Thread(target=_ask, daemon=True)
    t.start()
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline:
        req = broker.take()
        if req is not None:
            broker.resolve(req, True)
            break
        time.sleep(0.005)
    t.join(timeout=5)
    assert req_box == [True]
    assert broker.last_timed_out is False


# ---------------------------------------------------------------------------
# One reader, many askers
# ---------------------------------------------------------------------------

def test_two_threads_queue_rather_than_talk_over_each_other():
    broker = tc.TerminalConfirmBroker(timeout_s=5)
    results: list = []
    started = threading.Barrier(3)

    def _ask(cmd):
        started.wait(timeout=5)
        results.append(broker.callback("bash", {"command": cmd}, f"$ {cmd}"))

    threads = [threading.Thread(target=_ask, args=(c,), daemon=True)
               for c in ("one", "two")]
    for t in threads:
        t.start()
    started.wait(timeout=5)

    answered = 0
    deadline = time.monotonic() + 5
    while answered < 2 and time.monotonic() < deadline:
        req = broker.take()
        if req is None:
            time.sleep(0.005)
            continue
        assert broker.take() is req, (
            "take() must hand over the SAME request until it is resolved")
        broker.resolve(req, True)
        answered += 1
    for t in threads:
        t.join(timeout=5)
    assert results == [True, True]


def test_an_answer_to_an_expired_request_is_discarded():
    broker = tc.TerminalConfirmBroker(timeout_s=0.01)
    assert broker.callback("bash", {"command": "ls"}, "$ ls") is False
    stale = tc.ConfirmRequest(kind=tc.CONFIRM, tool="bash")
    stale.resolved = True
    assert broker.resolve(stale, True) is False, (
        "a late keystroke must not be applied to something already past")


def test_abort_denies_everything_in_flight_and_everything_after():
    broker = tc.TerminalConfirmBroker(timeout_s=5)
    results: list = []

    threads = [threading.Thread(
        target=lambda c=c: results.append(
            broker.callback("bash", {"command": c}, f"$ {c}")), daemon=True)
        for c in ("a", "b", "c")]
    for t in threads:
        t.start()
    deadline = time.monotonic() + 5
    while len(broker._queue) < 3 and time.monotonic() < deadline:
        time.sleep(0.005)
    denied = broker.abort_all()
    for t in threads:
        t.join(timeout=5)

    assert len(denied) == 3
    assert results == [False, False, False]
    assert broker.callback("bash", {"command": "later"}, "$ later") is False, (
        "after an abort, a request arriving late is refused without asking")


def test_the_abort_can_be_lifted_for_the_next_turn():
    broker = tc.TerminalConfirmBroker(timeout_s=0.01)
    broker.abort_all()
    broker.reset_abort()
    assert broker.aborted is False


# ---------------------------------------------------------------------------
# What may be offered
# ---------------------------------------------------------------------------

def test_an_outside_read_offers_no_always():
    req = tc.ConfirmRequest(
        kind=tc.CONFIRM, tool="read_file", args={"path": "/etc/hosts"},
        preview="[OUTSIDE-WORKSPACE READ]\n/etc/hosts")
    keys = {o.key for o in tc.options_for(req)}
    assert "A" not in keys and "k" not in keys and "e" not in keys, (
        "the only persistable form of this grant is a WRITABLE directory")
    assert "y" in keys and "n" in keys


def test_a_protected_path_is_approved_every_time_or_not_at_all():
    req = tc.ConfirmRequest(
        kind=tc.CONFIRM, tool="write_file", args={"path": "api_client.py"},
        preview="[SELF-MODIFICATION GUARD]\n--- a/api_client.py")
    keys = {o.key for o in tc.options_for(req)}
    assert "e" not in keys, "no session-wide pass for the protected core"


def test_an_unknown_remote_tool_offers_no_always():
    req = tc.ConfirmRequest(kind=tc.CONFIRM, tool="mcp__jira__create_issue",
                            args={"project": "ENG"}, preview="project: ENG")
    keys = {o.key for o in tc.options_for(req)}
    assert "A" not in keys and "k" not in keys, (
        "allow_patterns are shell regexes; there is no persist kind for this")


def test_a_shell_command_offers_the_exact_form_always():
    keys = {o.key for o in tc.options_for(_bash_req(), suggestion=r"^\s*pytest\b")}
    assert "A" in keys and "k" in keys


def test_no_generalisation_is_offered_when_none_is_safe():
    """kit_settings refuses to generalise `git push`, and so must the menu."""
    req = _bash_req("git push -u origin feature")
    broker = tc.TerminalConfirmBroker()
    suggestion = broker.suggest_pattern(req.command)
    keys = {o.key for o in tc.options_for(req, suggestion=suggestion)}
    assert "A" in keys
    if suggestion == broker.exact_pattern(req.command):
        assert "k" not in keys, (
            "offering a 'kind' identical to the exact command is a lie about "
            "what it covers")


def test_the_broker_never_invents_a_pattern():
    """kit_settings already encodes which heads may generalise."""
    broker = tc.TerminalConfirmBroker()
    assert broker.suggest_pattern("pytest -x tests/")
    exact = broker.exact_pattern("rm -rf build")
    assert exact.startswith(r"^\s*") and exact.endswith(r"\s*$")
    assert "rm\\ \\-rf" in exact or "rm\\ -rf" in exact


# ---------------------------------------------------------------------------
# Answering, from the loop
# ---------------------------------------------------------------------------

def test_yes_allows_and_no_refuses():
    broker = tc.TerminalConfirmBroker(timeout_s=5)
    agent, _engine, _err = _agent(broker)

    req = _write_req()
    agent._answer(req, _Keys(["y"]))
    assert req.decision is True

    req = _write_req()
    agent._answer(req, _Keys(["n"]))
    assert req.decision is False


def test_escape_is_a_refusal():
    broker = tc.TerminalConfirmBroker(timeout_s=5)
    agent, _engine, _err = _agent(broker)
    req = _write_req()
    agent._answer(req, _Keys(["\x1b"]))
    assert req.decision is False


def test_bare_enter_is_not_yes():
    """A default-yes turns approval into a rhythm."""
    broker = tc.TerminalConfirmBroker(timeout_s=5)
    agent, _engine, err = _agent(broker)
    req = _write_req()
    agent._answer(req, _Keys(["\r", "\r", "n"]))
    assert req.decision is False, (
        "Enter must not decide anything; only the explicit key did")
    assert "yes" in err.getvalue()


def test_abort_stops_the_turn_and_says_what_else_it_refused():
    broker = tc.TerminalConfirmBroker(timeout_s=5)
    agent, engine, err = _agent(broker)
    threading.Thread(
        target=lambda: broker.callback("bash", {"command": "b"}, "$ b"),
        daemon=True).start()
    deadline = time.monotonic() + 5
    while not broker._queue and time.monotonic() < deadline:
        time.sleep(0.005)

    req = _bash_req("a")
    agent._answer(req, _Keys(["a"]))
    assert req.decision is False
    assert engine.stopped is True
    assert "aborted" in err.getvalue()


def test_persisting_takes_a_second_keystroke_against_the_consequence():
    persisted: list[str] = []
    broker = tc.TerminalConfirmBroker(
        timeout_s=5, persist=lambda p: (persisted.append(p), (True, "ok"))[1])
    agent, _engine, err = _agent(broker)

    # First keystroke chooses "always", the second declines the persist.
    req = _bash_req()
    agent._answer(req, _Keys(["A", "n"]))
    assert persisted == [], "one keystroke must not persist anything"
    assert req.decision is True, "the action itself was still approved"
    assert "every other" in err.getvalue(), "the consequence has to be stated"
    assert "future one" in err.getvalue()

    req = _bash_req()
    agent._answer(req, _Keys(["A", "y"]))
    assert len(persisted) == 1
    assert "pytest" in persisted[0]


def test_the_help_key_explains_and_asks_again():
    broker = tc.TerminalConfirmBroker(timeout_s=5)
    agent, _engine, err = _agent(broker)
    req = tc.ConfirmRequest(
        kind=tc.CONFIRM, tool="read_file", args={"path": "/etc/hosts"},
        preview="[OUTSIDE-WORKSPACE READ]\n/etc/hosts")
    agent._answer(req, _Keys(["?", "n"]))
    assert "writable" in err.getvalue().lower()
    assert req.decision is False


def test_a_plan_is_approved_into_a_named_posture():
    broker = tc.TerminalConfirmBroker(timeout_s=5)
    agent, _engine, _err = _agent(broker)

    req = tc.ConfirmRequest(kind=tc.PLAN, preview="1. read\n2. edit")
    agent._answer(req, _Keys(["d"]))
    assert req.decision == {"approved": True, "new_mode": "default"}

    req = tc.ConfirmRequest(kind=tc.PLAN, preview="1. read")
    agent._answer(req, _Keys(["e"]))
    assert req.decision == {"approved": True, "new_mode": "acceptEdits"}

    req = tc.ConfirmRequest(kind=tc.PLAN, preview="1. read")
    agent._answer(req, _Keys(["n"]))
    assert req.decision == {"approved": False, "new_mode": "plan"}


def test_the_plan_screen_never_offers_unattended():
    keys = {o.key for o in tc.options_for(tc.ConfirmRequest(kind=tc.PLAN))}
    assert "b" not in keys
    labels = " ".join(o.label for o in tc.options_for(
        tc.ConfirmRequest(kind=tc.PLAN)))
    assert "bypass" not in labels.lower(), (
        "attention is on the plan text there, not on the mode word")


def test_a_question_is_answered_by_number():
    broker = tc.TerminalConfirmBroker(timeout_s=5)
    agent, _engine, _err = _agent(broker)
    req = tc.ConfirmRequest(kind=tc.ASK, payload={
        "question": "which environment?", "options": ["staging", "production"]})
    agent._answer(req, _Keys(["2"]))
    assert req.decision == {"answers": ["production"]}


# ---------------------------------------------------------------------------
# The frame itself
# ---------------------------------------------------------------------------

def test_a_diff_cannot_draw_its_own_approval_prompt():
    """File content is not allowed to look like the thing asking about it."""
    forged = ("--- a/x.py\n+++ b/x.py\n"
              "+\x1b[2J\x1b[H\n"
              "+└─ [y] yes  [n] no\n+Approved.\n")
    req = tc.ConfirmRequest(kind=tc.CONFIRM, tool="write_file",
                            args={"path": "x.py"}, preview=forged)
    frame = tc.render_request(req, theme=PLAIN, width=80)
    assert "\x1b" not in frame, "no escape from tool content may survive"
    # The forged line is still visible as text — the point is that it
    # cannot move the cursor or clear the screen to hide the real one.
    assert frame.rstrip().splitlines()[-1].startswith("└─ [y]")


def test_the_frame_names_the_tool_and_stays_within_the_width():
    frame = tc.render_request(_bash_req("x" * 500), theme=PLAIN, width=80)
    assert "bash" in frame
    assert all(len(line) <= 80 for line in frame.splitlines())


def test_a_remote_tool_frame_says_which_server():
    req = tc.ConfirmRequest(kind=tc.CONFIRM, tool="mcp__jira__create_issue",
                            args={}, preview="project: ENG")
    assert "jira:create_issue" in tc.render_request(req, theme=PLAIN)


def test_a_huge_preview_is_capped_and_says_how_much_is_left():
    req = tc.ConfirmRequest(kind=tc.CONFIRM, tool="write_file", args={},
                            preview="\n".join(f"line {i}" for i in range(200)))
    frame = tc.render_request(req, theme=PLAIN, body_lines=10)
    assert "190 more lines" in frame


# ---------------------------------------------------------------------------
# Nothing gets a silent yes
# ---------------------------------------------------------------------------

def test_a_prompt_that_cannot_be_read_is_a_refusal():
    broker = tc.TerminalConfirmBroker(timeout_s=5)
    agent, _engine, _err = _agent(broker)
    req = _write_req()
    agent._answer(req, None)          # no terminal to read from
    assert req.decision is False


def test_the_headless_runner_still_binds_no_callback():
    """cmd_run, the benchmark and the scheduler must not start asking."""
    import inspect
    from delfin.agent import cli as agent_cli

    src = inspect.getsource(agent_cli.cmd_run)
    assert "confirm_callback" not in src
    assert "TerminalConfirmBroker" not in src

    chat_src = inspect.getsource(agent_cli.cmd_chat)
    assert "sys.stdin.isatty()" in chat_src, (
        "the broker is bound for a terminal and nowhere else")


# ---------------------------------------------------------------------------
# The posture, and what counts as somebody having chosen it
# ---------------------------------------------------------------------------

def test_an_absent_settings_file_is_not_a_decision(tmp_path, monkeypatch):
    """The loader's own fallback is not a user's choice.

    kit_settings.load() returns default_mode="default" with no settings
    anywhere. Asking the merged view therefore reports a decision nobody
    made — and because resolve_posture correctly honours anything a user
    configured, that silently defeated plan-first on every machine without
    a settings file, which is most of them. Caught by watching a live run
    open in "default" and print the isolation warning.
    """
    from delfin.agent import cli as agent_cli, kit_settings

    monkeypatch.setattr(kit_settings, "USER_SETTINGS_PATH",
                        tmp_path / "nowhere" / "settings.json")
    assert agent_cli._persisted_default_mode(tmp_path) == ("", "")


def test_a_declared_mode_is_read_from_the_file(tmp_path, monkeypatch):
    import json
    from delfin.agent import cli as agent_cli, kit_settings

    monkeypatch.setattr(kit_settings, "USER_SETTINGS_PATH",
                        tmp_path / "nowhere" / "settings.json")
    ws = tmp_path / "project"
    (ws / ".delfin").mkdir(parents=True)
    (ws / ".delfin" / "settings.json").write_text(
        json.dumps({"kit": {"default_mode": "acceptEdits"}}))
    mode, source = agent_cli._persisted_default_mode(ws)
    assert mode == "acceptEdits"
    # The file that decided it travels with the decision: the banner used
    # to name a home-directory file for a mode a checked-out repository
    # had set, which points the reader at the wrong thing to edit.
    assert source.endswith("project/.delfin/settings.json"), source


def test_a_settings_file_about_something_else_declares_no_mode(
        tmp_path, monkeypatch):
    import json
    from delfin.agent import cli as agent_cli, kit_settings

    monkeypatch.setattr(kit_settings, "USER_SETTINGS_PATH",
                        tmp_path / "nowhere" / "settings.json")
    ws = tmp_path / "project"
    (ws / ".delfin").mkdir(parents=True)
    (ws / ".delfin" / "settings.json").write_text(
        json.dumps({"kit": {"allow_patterns": ["^pytest"]}}))
    assert agent_cli._persisted_default_mode(ws) == ("", "")


def test_nothing_configured_starts_in_plan(tmp_path, monkeypatch):
    """End to end: no files, no flag, therefore read-only first."""
    from delfin.agent import cli as agent_cli, kit_settings, launch_guard

    monkeypatch.setattr(kit_settings, "USER_SETTINGS_PATH",
                        tmp_path / "nowhere" / "settings.json")
    declared, source = agent_cli._persisted_default_mode(tmp_path)
    mode, why = launch_guard.resolve_posture(
        flag_mode="", persisted_mode=declared,
        settings_path=source or "~/.delfin/settings.json")
    assert mode == "plan"
    # Not the bare word "default", which sits next to "approval plan" in
    # the banner and reads as a statement about the mode.
    assert why == "nothing configured it"
