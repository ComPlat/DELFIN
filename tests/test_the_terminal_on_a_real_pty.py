"""`delfin-agent chat` driven through a real pty, against a scripted model.

Each test starts the real front door (``python -m delfin.agent.cli chat
--provider ollama``) in a pexpect pty, with the model replaced by the
local endpoint in ``tests/agent_terminal_probes/fake_model.py`` and HOME,
the scratch state and the workspace under ``tmp_path``. What is asserted
is the screen as the wrapping model terminal in ``tests/vt.py`` draws it,
or the requests the endpoint received -- never the raw byte stream, which
carries every intermediate frame of the bottom rows.

The items are the ones a scored drive of the terminal found broken:

  * a credential inside a command the model runs was masked in the
    approval dialog and printed whole on the tool line above it;
  * the approval key row was cut in the middle at 100 columns, so the
    key and the scope of "always commands like this" were not on screen;
  * the launch posture was announced to the first turn as a mid-turn
    switch by the user, which cost one extra model round;
  * a shell result over the 2000-character display slice was drawn as a
    cut JSON envelope instead of as its output;
  * a mistyped slash command was sent to the model as a prompt.
"""

from __future__ import annotations

import json

import pytest

pytest.importorskip("pexpect")
pytest.importorskip("termios")

from agent_terminal_probes import fake_model            # noqa: E402
from agent_terminal_probes.pty_chat import Chat, close_all  # noqa: E402

DEFAULT = ("--permission-mode", "default")
BYPASS = ("--permission-mode", "bypassPermissions", "--unattended")


@pytest.fixture
def chat(tmp_path):
    opened: list[Chat] = []

    def _open(*, args=DEFAULT, cols=100, rows=30, agent_settings=None) -> Chat:
        root = tmp_path / f"s{len(opened)}"
        if agent_settings:
            # The child's own settings file, in the HOME Chat gives it.
            (root / "home").mkdir(parents=True, exist_ok=True)
            (root / "home" / ".delfin_settings.json").write_text(
                json.dumps({"agent": dict(agent_settings)}))
        session = Chat(root, cols=cols, rows=rows, args=tuple(args))
        opened.append(session)
        assert session.wait_for("idle", 60) == "idle", (
            "no prompt within 60 s:\n" + "\n".join(session.visible()))
        return session

    yield _open
    close_all(opened)


def _requests(session: Chat) -> list[dict]:
    if not session.model_log.exists():
        return []
    return [json.loads(line)
            for line in session.model_log.read_text().splitlines()]


def _dialog(session: Chat) -> list[str]:
    rows = session.visible()
    start = max(i for i, row in enumerate(rows) if row.startswith("┌─"))
    return rows[start:]


# ---------------------------------------------------------------------------
# Approval
# ---------------------------------------------------------------------------

def test_a_credential_in_a_command_is_on_no_row_of_the_screen(chat):
    session = chat()
    session.send("do-bash-secret\r", 0.5)
    assert session.wait_for("approval", 30) == "approval", (
        "\n".join(session.visible()))
    screen = session.text()
    assert fake_model.FAKE_SECRET not in screen, (
        "the token the model put in the command is on screen:\n"
        + "\n".join(r for r in screen.splitlines()
                    if fake_model.FAKE_SECRET in r))
    session.send("n", 0.5)
    assert session.wait_for("idle", 30) == "idle"
    assert fake_model.FAKE_SECRET not in session.text()


@pytest.mark.parametrize("cols", [100, 60])
def test_every_approval_key_and_what_always_saves_is_on_screen(chat, cols):
    # The caged rule off for this child: where bubblewrap holds, `sleep 0`
    # is admitted unasked (#212), and the dialog this test reads existed
    # only on hosts without the cage. A command that is always asked
    # (`rm`) offers no "commands like this" pattern, so it cannot stand in.
    session = chat(cols=cols,
                   agent_settings={"caged_commands_run_unasked": False})
    session.send("do-bash sleep 0\r", 0.5)
    assert session.wait_for("approval", 30) == "approval", (
        "\n".join(session.visible()))
    dialog = _dialog(session)
    text = "\n".join(dialog)
    for key, label in (("y", "yes"), ("n", "no"),
                       ("A", "always this exact command"),
                       ("k", "always commands like this"),
                       ("a", "abort")):
        assert f"[{key}] {label}" in text.replace("\n   ", " "), (
            f"[{key}] {label} is not readable in the dialog:\n{text}")
    # What "always" covers, without pressing [?] first.
    assert "every future session" in text, text
    keys = [r for r in dialog if "[y] yes" in r or r.startswith("   [")]
    assert keys and not any("…" in r for r in keys), (
        "a key row was cut:\n" + text)
    # Enter is not an answer: the dialog is still open after it.
    session.send("\r", 0.8)
    assert session.state() == "approval"
    session.send("n", 0.5)
    assert session.wait_for("idle", 30) == "idle"


# ---------------------------------------------------------------------------
# The first turn
# ---------------------------------------------------------------------------

def test_the_launch_posture_is_not_announced_as_a_switch(chat):
    """No --permission-mode and nothing configured: the session starts in
    plan. That is the starting rule, stated in the banner and the prompt,
    not something the user changed while a turn ran."""
    session = chat(args=())
    assert session.turn("hello there") == "idle"
    assert "switched this session" not in session.text(), (
        "\n".join(session.visible()))
    asked = _requests(session)
    assert len(asked) == 1, (
        f"one typed line cost {len(asked)} model requests: "
        + json.dumps([a["last_user"][:80] for a in asked]))


# ---------------------------------------------------------------------------
# Long output
# ---------------------------------------------------------------------------

def _result_rows(session: Chat) -> list[str]:
    rows = session.screen().text()
    start = max(i for i, row in enumerate(rows) if "⏺ " in row)
    return rows[start:]


def test_a_long_shell_result_shows_its_last_lines(chat):
    """seq 1 400: stdout (1.9 kB escaped) fits the 2000-character display
    slice, the envelope fields after it do not -- so the slice is no longer
    valid JSON although the output in it is complete."""
    session = chat(args=BYPASS)
    assert session.turn("do-bash seq 1 400") == "idle"
    rows = _result_rows(session)
    text = "\n".join(rows)
    assert '"exit_code"' not in text, text
    assert [r.strip() for r in rows if r.strip() in ("398", "399", "400")] \
        == ["398", "399", "400"], text


def test_a_result_cut_for_display_says_so_and_shows_its_head(chat):
    session = chat(args=BYPASS)
    assert session.turn("do-bash seq 1 2000") == "idle"
    rows = _result_rows(session)
    text = "\n".join(rows)
    assert '"exit_code"' not in text, text
    stripped = [r.strip() for r in rows]
    assert "1" in stripped and "2" in stripped, text
    assert "more" in text and "ctrl+o" in text.lower(), text
    assert "100" not in stripped, "the head is three lines, not the slice"
    # The line points at Ctrl+O, so Ctrl+O at the idle prompt must show
    # what was kept -- as output lines, not as the envelope.
    session.send("draft", 0.3)
    session.send("\x0f", 1.0)
    assert session.wait_for("idle", 10) == "idle"
    after = [r.strip() for r in session.screen().text()]
    assert "100" in after, "\n".join(session.visible())
    assert '"exit_code"' not in "\n".join(after)
    assert session.visible()[-3].startswith("> draft"), (
        "the typed draft did not survive Ctrl+O:\n"
        + "\n".join(session.visible()[-4:]))


# ---------------------------------------------------------------------------
# Slash commands
# ---------------------------------------------------------------------------

def test_a_mistyped_command_names_the_nearest_and_sends_nothing(chat):
    session = chat()
    session.send("/hlep\r", 1.5)
    session.wait_until(lambda s: "unknown command /hlep" in s.text())
    assert session.wait_for("idle", 15) == "idle"
    text = session.text()
    assert "/help" in text and "unknown command /hlep" in text, (
        "\n".join(session.visible()))
    assert _requests(session) == [], "the typo went to the model"
    session.send("\x1b[A", 0.6)
    assert session.wait_until(
        lambda s: s.visible()[-3].startswith("> /hlep")), (
        "↑ did not bring the line back:\n" + "\n".join(session.visible()[-4:]))


@pytest.mark.parametrize("line,near", [
    ("/hlep", "/help"), ("/stauts", "/status"), ("/comapct now", "/compact"),
    ("/tmp", ""), ("/scratch what is here", ""), ("/data", ""),
    ("/home", ""), ("/help", ""), ("/home/u/x.py", ""),
])
def test_only_a_near_miss_of_a_builtin_is_held_back(line, near):
    from delfin.agent import repl_commands as rc
    assert rc.nearest_builtin(line) == near


# ---------------------------------------------------------------------------
# Teardown
# ---------------------------------------------------------------------------

def test_a_session_with_a_draft_in_the_box_still_leaves_on_ctrl_d(chat):
    """Ctrl+D on a non-empty line deletes forward (readline), not leaves.

    A close() that sends only Ctrl+D leaves such a session to pexpect's
    forced close, whose SIGKILL is given 0.1 s to be reaped.
    """
    session = chat()
    session.send("a draft that is never sent", 0.5)
    assert session.close() == 0
    assert session.child.signalstatus is None


_IGNORES = ("import signal, sys, time\n"
            "for s in {sigs}: signal.signal(s, signal.SIG_IGN)\n"
            "print('ready', flush=True)\n"
            "time.sleep(120)\n")


def _stub(sigs: str):
    import sys
    import pexpect
    child = pexpect.spawn(sys.executable, ["-c", _IGNORES.format(sigs=sigs)],
                          encoding="utf-8", timeout=10)
    child.expect("ready")
    return child


def test_a_child_that_only_sigterm_ends_fails_its_teardown():
    """Ending on SIGTERM is not leaving: the teardown says which step."""
    from agent_terminal_probes.pty_chat import end_child
    child = _stub("(signal.SIGHUP, signal.SIGINT)")
    with pytest.raises(AssertionError, match=r"Ctrl\+D.*SIGTERM ended it"):
        end_child(child, leave_s=1.0)
    assert not child.isalive() and child.signalstatus == 15


def test_a_child_that_ignores_every_leave_is_killed_and_reported():
    from agent_terminal_probes.pty_chat import end_child
    child = _stub("(signal.SIGHUP, signal.SIGINT, signal.SIGTERM)")
    with pytest.raises(AssertionError, match=r"after SIGTERM.*SIGKILL ended"):
        end_child(child, leave_s=1.0)
    assert not child.isalive() and child.signalstatus == 9
    assert child.closed


def test_a_path_is_still_a_message(chat):
    """The suggestion must not swallow text that only starts with a slash."""
    session = chat()
    assert session.turn("/tmp/notes.txt what is in it") == "idle"
    asked = _requests(session)
    assert len(asked) == 1
    assert asked[0]["last_user"].startswith("/tmp/notes.txt")
