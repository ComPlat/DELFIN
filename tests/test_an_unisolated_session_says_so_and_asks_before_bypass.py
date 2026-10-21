"""Without filesystem isolation the posture is stated, and bypass is asked.

Where no layer (bubblewrap, Landlock, Seatbelt) holds a shell command, or
the user switched isolation off (``--no-isolate``, ``agent.bash_isolation
= "off"``), what stops a command is a reading of its text, and a command
that assembles its words at run time walks past any reading. Pinned here,
on a host simulated without any layer:

* the status (isolated or not, by which layer, and why not) is one
  function, shown at session start and by /status and the status line;
* bypassPermissions: an attended start asks once per session; with nobody
  to ask, or a "no", the session starts in ``default`` with a notice and
  keeps running; acceptEdits, default and plan start as asked;
* a command that cannot be confined in a session that promises
  confinement (locked scope, forced isolation) is asked in every mode,
  and refused with nobody to ask, naming what works instead;
* an unattended ``-p`` run on such a host still completes a normal task:
  read, write in the workspace, run the tests.
"""

from __future__ import annotations

import io
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from delfin.agent import api_client as A
from delfin.agent import cli as agent_cli
from delfin.agent import isolation_posture as P


# -- the simulated host ------------------------------------------------------

def _host(monkeypatch, *, bwrap=False, landlock=False, seatbelt=False,
          setting="auto"):
    monkeypatch.setattr(A.shutil, "which",
                        lambda name, *a, **k: "/usr/bin/bwrap" if (
                            name == "bwrap" and bwrap) else None)
    monkeypatch.setattr(A, "_bwrap_functional", lambda: bwrap)
    monkeypatch.setattr(A, "_landlock_functional", lambda: landlock)
    monkeypatch.setattr(A, "_seatbelt_functional", lambda: seatbelt)
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: False)
    import delfin.user_settings as US
    monkeypatch.setattr(US, "load_settings",
                        lambda *a, **k: {"agent": {"bash_isolation": setting}})


@pytest.fixture(autouse=True)
def _fresh(monkeypatch):
    A.set_bash_isolation_override("")
    monkeypatch.setattr(P, "_BYPASS_CONFIRMED", False)
    yield
    A.set_bash_isolation_override("")


class _Engine:
    """The parts of an engine the start check and cmd_run touch."""

    session_id = "sess-iso-1"
    token_usage = {"input": 1, "output": 1}

    def __init__(self, workspace: Path, mode: str):
        self.kit_permissions = A.KitToolPermissions(workspace=workspace,
                                                    mode=mode)
        self.messages = [{"role": "user", "content": "x"}]

    def set_kit_permission_mode(self, mode):
        self.kit_permissions.mode = mode
        return True

    def record_cycle_outcome(self, *a, **k):
        return None

    def export_state(self):
        return {"engine_state": {}}


def _ws(tmp_path) -> Path:
    ws = tmp_path / "ws"
    ws.mkdir(exist_ok=True)
    return ws


# -- the status --------------------------------------------------------------

def test_no_layer_reads_as_not_isolated_and_says_why(monkeypatch):
    _host(monkeypatch)
    st = P.isolation_status()
    assert not st.isolated and not st.refuses and st.layer == ""
    assert "NOT isolated" in st.line()
    assert "bubblewrap" in st.line() and "Landlock" in st.line()


@pytest.mark.parametrize("layer,name", [("bwrap", "bwrap"),
                                        ("landlock", "Landlock"),
                                        ("seatbelt", "Seatbelt")])
def test_a_working_layer_is_named(monkeypatch, layer, name):
    _host(monkeypatch, **{layer: True})
    st = P.isolation_status()
    assert st.isolated and st.layer == name
    assert name in st.line()


def test_no_isolate_is_not_isolated_on_a_host_that_could(monkeypatch):
    _host(monkeypatch, bwrap=True)
    A.set_bash_isolation_override("off")
    st = P.isolation_status()
    assert not st.isolated
    assert "--no-isolate" in st.line()


def test_the_setting_off_is_not_isolated_and_named(monkeypatch):
    _host(monkeypatch, bwrap=True, setting="off")
    st = P.isolation_status()
    assert not st.isolated
    assert "agent.bash_isolation" in st.line()


def test_a_locked_scope_without_a_layer_refuses_rather_than_runs(
        monkeypatch, tmp_path):
    _host(monkeypatch)
    # The office role is one of the carriers of a locked scope.
    perms = A.KitToolPermissions(workspace=_ws(tmp_path), mode="default",
                                 agent_role="office_agent")
    assert perms.scope_locked
    st = P.isolation_status(perms)
    assert st.refuses and not st.isolated


# -- the start check, every mode ---------------------------------------------

def _start(monkeypatch, tmp_path, mode, *, can_ask, answer=None):
    asked: list[str] = []

    def _ask(question):
        asked.append(question)
        return bool(answer)

    eng = _Engine(_ws(tmp_path), mode)
    notes = P.apply_at_start(eng, can_ask=can_ask,
                             ask=_ask if can_ask else None)
    return eng.kit_permissions.mode, notes, asked


@pytest.mark.parametrize("mode", ["plan", "default"])
@pytest.mark.parametrize("can_ask", [True, False])
def test_plan_and_default_start_as_asked_without_a_question(
        monkeypatch, tmp_path, mode, can_ask):
    _host(monkeypatch)
    got, _notes, asked = _start(monkeypatch, tmp_path, mode, can_ask=can_ask)
    assert got == mode and asked == []


@pytest.mark.parametrize("can_ask", [True, False])
def test_accept_edits_starts_with_a_notice_and_no_question(
        monkeypatch, tmp_path, can_ask):
    _host(monkeypatch)
    got, notes, asked = _start(monkeypatch, tmp_path, "acceptEdits",
                               can_ask=can_ask)
    assert got == "acceptEdits" and asked == []
    assert any("NOT isolated" in n for n in notes), notes


def test_unattended_bypass_falls_back_to_default_and_says_so(
        monkeypatch, tmp_path):
    _host(monkeypatch)
    got, notes, asked = _start(monkeypatch, tmp_path, "bypassPermissions",
                               can_ask=False)
    assert got == "default" and asked == []
    text = " ".join(notes)
    assert "bypassPermissions" in text and "default" in text
    assert "NOT isolated" in text


def test_attended_bypass_asks_once_naming_the_consequence(
        monkeypatch, tmp_path):
    _host(monkeypatch)
    got, notes, asked = _start(monkeypatch, tmp_path, "bypassPermissions",
                               can_ask=True, answer=True)
    assert got == "bypassPermissions" and len(asked) == 1
    assert "outside the workspace" in asked[0]
    # One dialog per session: a second start check in the same process
    # does not ask again.
    got2, _n, asked2 = _start(monkeypatch, tmp_path, "bypassPermissions",
                              can_ask=True, answer=False)
    assert got2 == "bypassPermissions" and asked2 == []


def test_attended_bypass_refused_starts_in_default(monkeypatch, tmp_path):
    _host(monkeypatch)
    got, notes, asked = _start(monkeypatch, tmp_path, "bypassPermissions",
                               can_ask=True, answer=False)
    assert got == "default" and len(asked) == 1
    assert any("default" in n for n in notes)


def test_no_isolate_keeps_bypass_after_the_confirmation(monkeypatch, tmp_path):
    _host(monkeypatch, bwrap=True)
    A.set_bash_isolation_override("off")
    got, _notes, asked = _start(monkeypatch, tmp_path, "bypassPermissions",
                                can_ask=True, answer=True)
    assert got == "bypassPermissions" and len(asked) == 1


@pytest.mark.parametrize("mode", ["plan", "default", "acceptEdits",
                                  "bypassPermissions"])
@pytest.mark.parametrize("can_ask", [True, False])
def test_an_isolated_host_asks_nothing_and_changes_nothing(
        monkeypatch, tmp_path, mode, can_ask):
    _host(monkeypatch, landlock=True)
    got, notes, asked = _start(monkeypatch, tmp_path, mode, can_ask=can_ask)
    assert got == mode and asked == [] and notes == []


# -- the -p path -------------------------------------------------------------

def _one_shot(monkeypatch, tmp_path, run_once, *, mode="bypassPermissions"):
    ws = _ws(tmp_path)
    holder = {}

    def _build(args):
        holder["engine"] = _Engine(ws, args.permission_mode or "default")
        return holder["engine"]

    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    monkeypatch.setattr(agent_cli, "_build_engine", _build)
    monkeypatch.setattr(agent_cli, "_run_once", run_once)
    monkeypatch.setattr(agent_cli, "_save_session", lambda *a, **k: "sess-iso-1")
    monkeypatch.setattr(agent_cli, "_resume_or_create", lambda *a, **k: None)
    monkeypatch.setattr("sys.stdin", io.TextIOWrapper(io.BytesIO(b"")))
    # cmd_run moves into --cwd and stays there.
    monkeypatch.chdir(Path.cwd())
    # --unattended: a headless bypass run needs it (_headless_posture).
    rc = agent_cli.main(["-p", "do it", "--permission-mode", mode,
                         "--unattended", "--cwd", str(ws)])
    return rc, holder["engine"]


def _answer(text="done"):
    return {"text": text, "tool_calls": [], "input_tokens": 1,
            "output_tokens": 1, "error": "", "session_id": "sess-iso-1"}


def test_unattended_p_bypass_runs_in_default_with_a_notice(
        monkeypatch, tmp_path, capsys):
    _host(monkeypatch)
    seen = {}

    def _run_once(engine, prompt, **kw):
        seen["mode"] = engine.kit_permissions.mode
        return _answer()

    rc, _eng = _one_shot(monkeypatch, tmp_path, _run_once)
    err = capsys.readouterr().err
    assert rc == 0
    assert seen["mode"] == "default"
    assert "NOT isolated" in err and "bypassPermissions" in err


def test_unattended_p_bypass_keeps_bypass_where_a_layer_holds(
        monkeypatch, tmp_path, capsys):
    _host(monkeypatch, landlock=True)
    seen = {}

    def _run_once(engine, prompt, **kw):
        seen["mode"] = engine.kit_permissions.mode
        return _answer()

    rc, _eng = _one_shot(monkeypatch, tmp_path, _run_once)
    assert rc == 0 and seen["mode"] == "bypassPermissions"
    assert "NOT isolated" not in capsys.readouterr().err


def test_an_unattended_p_run_without_isolation_completes_a_normal_task(
        monkeypatch, tmp_path, capsys):
    """Read, write in the workspace, run the tests: through the real
    executor and gate, after the fallback from bypass to default."""
    _host(monkeypatch)
    results = {}

    def _run_once(engine, prompt, **kw):
        perms = engine.kit_permissions
        ex = A._DocToolExecutor()
        results["write_src"] = ex.execute("write_file", {
            "path": "calc.py",
            "content": "def add(a, b):\n    return a + b\n"}, perms)
        results["write_test"] = ex.execute("write_file", {
            "path": "test_calc.py",
            "content": "from calc import add\n\n\n"
                       "def test_add():\n    assert add(2, 3) == 5\n"},
            perms)
        results["read"] = ex.execute("read_file", {"path": "calc.py"}, perms)
        results["tests"] = ex.execute("run_tests",
                                      {"target": "test_calc.py"}, perms)
        return _answer()

    rc, eng = _one_shot(monkeypatch, tmp_path, _run_once)
    assert rc == 0
    ws = eng.kit_permissions.workspace
    assert (Path(ws) / "calc.py").read_text().startswith("def add")
    assert (Path(ws) / "test_calc.py").exists()
    assert "return a + b" in results["read"], results["read"]
    rt = results["tests"]
    assert "error" not in json.loads(rt) or not json.loads(rt).get("error"), rt
    assert "1 passed" in rt or '"passed": 1' in rt, rt


# -- a command that cannot be confined ---------------------------------------

def _locked(tmp_path, mode, callback):
    perms = A.KitToolPermissions(workspace=_ws(tmp_path), mode=mode,
                                 agent_role="office_agent")
    assert perms.scope_locked
    perms.confirm_callback = callback
    return perms


@pytest.mark.parametrize("mode", ["default", "acceptEdits",
                                  "bypassPermissions"])
def test_an_unconfinable_command_is_asked_in_every_mode(
        monkeypatch, tmp_path, mode):
    _host(monkeypatch)
    asked = []

    def _cb(tool, args, preview):
        asked.append(preview)
        return True

    perms = _locked(tmp_path, mode, _cb)
    cmd = "ls"
    assert A._DocToolExecutor()._run_permission_gate(
        "bash", {"command": cmd}, perms) is None
    assert any("[NO ISOLATION]" in p for p in asked), asked
    # Approved once: that command runs, outside a cage nothing here builds.
    argv = A._bash_isolation_argv(cmd, perms.workspace, perms)
    assert argv[-1] == cmd and "refused" not in argv
    # Only once: the next run of it is refused again unless asked again.
    argv2 = A._bash_isolation_argv(cmd, perms.workspace, perms)
    assert "refused" in argv2


def test_an_unconfinable_command_denied_is_not_run(monkeypatch, tmp_path):
    _host(monkeypatch)
    perms = _locked(tmp_path, "default", lambda *a: False)
    err = A._DocToolExecutor()._run_permission_gate(
        "bash", {"command": "ls"}, perms)
    assert err and "denied" in err


def test_an_unconfinable_command_with_nobody_to_ask_names_what_works(
        monkeypatch, tmp_path):
    _host(monkeypatch)
    perms = _locked(tmp_path, "bypassPermissions", None)
    err = A._DocToolExecutor()._run_permission_gate(
        "bash", {"command": "ls"}, perms)
    assert err and err.startswith("blocked")
    assert "read_file" in err and "write_file" in err


@pytest.mark.parametrize("deny", [None, "[SCRIPT IN AN UNKNOWN DIRECTORY]",
                                  "[NO ISOLATION]"])
def test_an_unknown_script_dir_and_no_isolation_are_both_asked(
        monkeypatch, tmp_path, deny):
    """Two consequences, two dialogs: approving the unknown-directory
    script does not skip the [NO ISOLATION] question, and a no to either
    refuses the command."""
    _host(monkeypatch)
    asked = []

    def _cb(tool, args, preview):
        asked.append(preview)
        return not (deny and preview.startswith(deny))

    perms = _locked(tmp_path, "default", _cb)
    verdict = A._DocToolExecutor()._run_permission_gate(
        "bash", {"command": "cd $D && python3 y.py"}, perms)
    heads = [p.split("\n", 1)[0] for p in asked]
    if deny is None:
        assert heads[:2] == ["[SCRIPT IN AN UNKNOWN DIRECTORY]",
                             "[NO ISOLATION]"], heads
    else:
        assert verdict and "denied" in verdict, verdict
        assert heads[-1] == deny, heads


def test_a_confinable_command_asks_nothing_extra(monkeypatch, tmp_path):
    _host(monkeypatch, landlock=True)
    perms = _locked(tmp_path, "bypassPermissions", None)
    assert A._DocToolExecutor()._run_permission_gate(
        "bash", {"command": "ls"}, perms) is None


# -- /status and the status line ---------------------------------------------

def test_status_names_the_isolation(monkeypatch, tmp_path):
    _host(monkeypatch)
    from delfin.agent import repl_commands as RC

    class _Ctx:
        engine = _Engine(_ws(tmp_path), "default")
        workspace = tmp_path

    _Ctx.engine.get_status = lambda: {}
    out = RC._status(_Ctx, "").output
    assert "isolation" in out and "NOT isolated" in out


def test_the_status_line_marks_an_unisolated_session(monkeypatch, tmp_path):
    _host(monkeypatch)
    from delfin.agent import repl as R
    agent = R.TerminalAgent.__new__(R.TerminalAgent)
    eng = _Engine(_ws(tmp_path), "default")
    eng.get_status = lambda: {}
    agent.engine = eng
    agent._spin = 0
    agent._turn_t0 = 0.0
    agent._turn_base = (0, 0, 0.0)
    agent._streamed_chars = 0

    class _Theme:
        @staticmethod
        def dim(s):
            return s

    class _Tr:
        theme = _Theme()
        width = 200

    agent.transcript = _Tr()
    assert "unisolated" in agent._status_line()


def test_the_interactive_banner_names_a_setting_that_switched_it_off(
        monkeypatch, tmp_path):
    _host(monkeypatch, bwrap=True, setting="off")
    eng = _Engine(_ws(tmp_path), "default")

    monkeypatch.setattr(agent_cli, "_wordmark", lambda w: [])
    banner = agent_cli._startup_banner(eng, _Report(), _ws(tmp_path))
    assert "agent.bash_isolation" in banner
    assert "isolation  on" not in banner


# -- the interactive path ----------------------------------------------------

class _Report:
    refused = False
    questions = ()
    granted_dirs = ()
    read_dirs = ()
    git = SimpleNamespace(is_repo=False, branch="", dirty=())

    def render(self):
        return ""


class _TtyIn(io.StringIO):
    def isatty(self):
        return True


def _interactive(monkeypatch, tmp_path, argv, answer):
    from delfin.agent import launch_guard, repl
    ws = _ws(tmp_path)
    holder = {}

    def _build(args):
        holder["engine"] = _Engine(ws, args.permission_mode or "plan")
        return holder["engine"]

    class _Agent:
        def __init__(self, engine, broker=None, opts=None):
            holder["banner"] = opts.banner
            holder["mode_at_run"] = engine.kit_permissions.mode

        def run(self, positional):
            return 0

        def leave_held(self):
            import contextlib
            return contextlib.nullcontext(self)

        def begin_leaving(self, code, *, terminal_gone=False):
            pass

        def leave_stages(self):
            return None         # session-end stages inline

    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    monkeypatch.setattr(launch_guard, "inspect_launch_dir",
                        lambda *a, **k: _Report())
    monkeypatch.setattr(agent_cli, "_build_engine", _build)
    monkeypatch.setattr(agent_cli, "_open_session", lambda *a, **k: True)
    monkeypatch.setattr(agent_cli, "_persisted_default_mode",
                        lambda ws: ("", ""))
    monkeypatch.setattr(agent_cli, "_save_session", lambda *a, **k: "s")
    monkeypatch.setattr(agent_cli, "_wordmark", lambda w: [])
    monkeypatch.setattr(repl, "TerminalAgent", _Agent)
    monkeypatch.setattr("sys.stdin", _TtyIn(answer))
    cwd = Path.cwd()
    try:
        rc = agent_cli.main([*argv, "--cwd", str(ws)])
    finally:
        import os
        os.chdir(cwd)
    return rc, holder


def test_interactive_bypass_without_isolation_asks_and_keeps_it_on_yes(
        monkeypatch, tmp_path, capsys):
    _host(monkeypatch)
    rc, h = _interactive(monkeypatch, tmp_path,
                         ["--permission-mode", "bypassPermissions",
                          "--unattended"], "y\n")
    err = capsys.readouterr().err
    assert rc == 0 and h["mode_at_run"] == "bypassPermissions"
    assert "outside the workspace" in err
    assert "NOT isolated" in h["banner"]


def test_interactive_bypass_without_isolation_on_no_starts_in_default(
        monkeypatch, tmp_path):
    _host(monkeypatch)
    rc, h = _interactive(monkeypatch, tmp_path,
                         ["--permission-mode", "bypassPermissions",
                          "--unattended"], "n\n")
    assert rc == 0 and h["mode_at_run"] == "default"
    assert "NOT isolated" in h["banner"]


@pytest.mark.parametrize("mode", ["plan", "default", "acceptEdits"])
def test_interactive_other_modes_start_without_a_question(
        monkeypatch, tmp_path, capsys, mode):
    _host(monkeypatch)
    rc, h = _interactive(monkeypatch, tmp_path,
                         ["--permission-mode", mode], "")
    assert rc == 0 and h["mode_at_run"] == mode
    assert "[y/N]" not in capsys.readouterr().out
    assert "NOT isolated" in h["banner"]


def test_interactive_no_isolate_is_confirmed_then_kept(
        monkeypatch, tmp_path):
    _host(monkeypatch, bwrap=True)
    rc, h = _interactive(monkeypatch, tmp_path,
                         ["--permission-mode", "bypassPermissions",
                          "--unattended", "--no-isolate"], "y\n")
    assert rc == 0 and h["mode_at_run"] == "bypassPermissions"
    assert "--no-isolate" in h["banner"]


def test_interactive_isolated_bypass_asks_nothing(
        monkeypatch, tmp_path, capsys):
    _host(monkeypatch, bwrap=True)
    rc, h = _interactive(monkeypatch, tmp_path,
                         ["--permission-mode", "bypassPermissions",
                          "--unattended"], "")
    assert rc == 0 and h["mode_at_run"] == "bypassPermissions"
    assert "[y/N]" not in capsys.readouterr().out


# -- the unattended fallback: off-list commands are refused ----------------

#: Harmless and off the auto-allow list (asserted in the test).
_OFF_LIST = "expr 2 + 3"


def _fallback_perms(monkeypatch, tmp_path, runner):
    _host(monkeypatch)
    monkeypatch.setenv("DELFIN_CLASSIFIER", "on")
    eng = _Engine(_ws(tmp_path), "bypassPermissions")
    P.apply_at_start(eng, can_ask=False)
    perms = eng.kit_permissions
    assert perms.mode == "default" and perms.unisolated_fallback
    assert perms.confirm_callback is None
    assert not perms.matches_bash_auto_allow(_OFF_LIST)
    perms.classifier_runner = runner
    return perms


def _run(perms, cmd=_OFF_LIST):
    return A._DocToolExecutor().execute("bash", {"command": cmd}, perms)


def _down(prompt, model):
    raise ConnectionError("endpoint down")


@pytest.mark.parametrize("answer", [
    "DECISION: allow\nREASON: arithmetic, harmless.",
    "DECISION: refuse\nREASON: not part of the task.",
    "DECISION: ask\nREASON: unclear.",
    "raise", None])
def test_the_fallback_refuses_an_off_list_command_whatever_the_classifier_says(
        monkeypatch, tmp_path, answer):
    """No layer holds the command and nobody reads it: refused, naming
    what works, without consulting the classifier
    (test_an_unattended_run_without_isolation_runs_only_the_auto_allow_
    list.py)."""
    seen = []
    if answer is None:
        runner = None
    elif answer == "raise":
        runner = _down
    else:
        def runner(prompt, model, _a=answer):
            seen.append(prompt)
            return _a
    out = _run(_fallback_perms(monkeypatch, tmp_path, runner))
    assert "blocked" in out and "auto-allow" in out, out
    assert "read_file" in out and "run_tests" in out, out
    assert seen == []


def test_the_fallback_keeps_the_hard_refusals(monkeypatch, tmp_path):
    def _allow(prompt, model):
        return "DECISION: allow\nREASON: fine."

    perms = _fallback_perms(monkeypatch, tmp_path, _allow)
    out = _run(perms, "cat ~/.ssh/id_rsa")
    # Refused by the secret-deny tier, before the classifier is reached.
    assert "secret-deny" in out and "auto-allow list" not in out, out


def test_without_the_fallback_mark_nothing_changes(monkeypatch, tmp_path):
    _host(monkeypatch)
    monkeypatch.setenv("DELFIN_CLASSIFIER", "on")
    perms = A.KitToolPermissions(workspace=_ws(tmp_path), mode="default")
    perms.classifier_runner = lambda p, m: "DECISION: allow\nREASON: x."
    out = _run(perms)
    assert "auto-allow" in out or "allow_pattern" in out, out


# -- the dashboard start: asked in the approval panel ------------------------

def _wait(predicate, seconds=5.0):
    import time
    end = time.monotonic() + seconds
    while time.monotonic() < end:
        if predicate():
            return True
        time.sleep(0.02)
    return False


@pytest.mark.parametrize("answer", [True, False])
def test_the_dashboard_start_asks_in_the_approval_panel(
        monkeypatch, tmp_path, answer):
    from delfin.agent.kit_confirm import KitConfirmBroker
    _host(monkeypatch)
    broker = KitConfirmBroker(default_timeout_s=10)
    eng = _Engine(_ws(tmp_path), "bypassPermissions")
    confirmed = []

    def _ask(question):
        return broker.callback("bypassPermissions", {"isolation": "none"},
                               "[NO ISOLATION]\n" + question)

    notes = P.apply_at_start(eng, can_ask=True, ask=_ask, background=True,
                             on_confirmed=lambda: confirmed.append(1))
    # Not held: the session runs in default while the panel waits.
    assert eng.kit_permissions.mode == "default"
    assert any("approval panel" in n for n in notes), notes
    assert _wait(lambda: broker.pending_requests())
    req = broker.pending_requests()[0]
    assert "[NO ISOLATION]" in req.preview
    assert "outside the workspace" in req.preview
    broker.resolve(req, answer)
    if answer:
        assert _wait(lambda: confirmed)
        assert eng.kit_permissions.mode == "bypassPermissions"
    else:
        assert not _wait(lambda: confirmed, 0.3)
        assert eng.kit_permissions.mode == "default"


def test_the_dashboard_calls_the_start_check_where_it_builds_and_switches():
    import inspect
    from delfin.dashboard import tab_agent
    src = inspect.getsource(tab_agent)
    # Engine construction, the chip cycle, a session restore and the
    # permissions dropdown.
    assert src.count("_settle_isolation_posture(") >= 5
    built = src.index("ctx.agent_engine = engine")
    assert "_settle_isolation_posture(engine)" in src[built:built + 200]
