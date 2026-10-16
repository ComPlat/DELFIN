"""The terminal's status row says what the running turn waits on.

Input: ``delfin-agent chat`` on a real pty against the scripted endpoint
(``tests/agent_terminal_probes``), with a prompt that makes the endpoint
connect and stay silent (``do-delay``), run a shell command
(``do-bash``), or ask for an approval. Output: the status row of the
turn, read off the screen model at two moments.

What is asserted: the row names the phase -- waiting for the model,
running a tool, waiting for the person -- with a clock for THAT phase
that advances while nothing arrives; a phase past its threshold adds one
line saying what is slow and how to stop it; the row is gone once the
turn ends, also when it failed.

The defect: the row read ``⠴ 12s fake-model default ↑0 ↓… esc to
interrupt`` in every phase. A connected endpoint that had sent no token
(KIT GLM: 77-136 s to the first token) and a test suite running under
``bash`` looked the same, and the only phase-specific line, the
first-byte tick, exists only until the response headers arrive.
"""

from __future__ import annotations

import re
import time

import pytest

pytest.importorskip("pexpect")
pytest.importorskip("termios")

from agent_terminal_probes.pty_chat import Chat         # noqa: E402

DEFAULT = ("--permission-mode", "default")
BYPASS = ("--permission-mode", "bypassPermissions", "--unattended")


@pytest.fixture
def chat(tmp_path):
    opened: list[Chat] = []

    def _open(*, args=DEFAULT, env=None) -> Chat:
        session = Chat(tmp_path / f"s{len(opened)}", cols=120, rows=30,
                       args=tuple(args), env=env)
        opened.append(session)
        assert session.wait_for("idle", 60) == "idle", (
            "no prompt within 60 s:\n" + "\n".join(session.visible()))
        return session

    yield _open
    for session in opened:
        session.close()


def _status(session: Chat) -> str:
    rows = [r for r in session.visible()[-8:] if "esc to interrupt" in r]
    return rows[-1] if rows else ""


def _seconds(row: str, phrase: str) -> int | None:
    m = re.search(re.escape(phrase) + r"[^…]*… (\d+) s", row)
    return int(m.group(1)) if m else None


def _start(session: Chat, prompt: str) -> None:
    session.send(prompt, 0.2)
    session.send("\r", 0.3)


def _sample_until(session: Chat, at: float, t0: float) -> str:
    session.pump(max(0.0, t0 + at - time.monotonic()))
    return _status(session)


def test_a_connected_silent_endpoint_is_waited_on_with_a_clock(chat):
    session = chat()
    t0 = time.monotonic()
    _start(session, "do-delay 6")
    early = _sample_until(session, 2.0, t0)
    late = _sample_until(session, 4.5, t0)
    a = _seconds(early, "waiting for fake-model")
    b = _seconds(late, "waiting for fake-model")
    assert a is not None and b is not None, (
        "the row does not say it waits for the model:\n"
        f"  at 2 s:   {early!r}\n  at 4.5 s: {late!r}")
    assert b > a, f"the clock did not advance: {early!r} -> {late!r}"
    assert session.wait_for("idle", 30) == "idle"
    assert _status(session) == "", "the row outlived the turn"


def test_past_the_threshold_one_line_says_what_is_slow_and_how_to_stop(chat):
    session = chat(env={"DELFIN_SLOW_FIRST_TOKEN_S": "2"})
    _start(session, "do-delay 5")
    session.pump(3.5)
    text = session.text()
    hints = [r for r in text.splitlines() if "slow:" in r]
    assert hints, "no slow line past the threshold:\n" + "\n".join(
        session.visible())
    assert "fake-model" in hints[0] and "Esc" in hints[0], hints[0]
    assert session.wait_for("idle", 30) == "idle"
    assert len([r for r in session.text().splitlines() if "slow:" in r]) == 1, (
        "the slow line is said once per phase, not per repaint")


def test_a_running_tool_is_named_with_its_argument_and_a_clock(chat):
    # The action watcher is off: it may ask about any command it cannot
    # judge, and an approval is the other test's phase, not this one's.
    session = chat(args=BYPASS, env={"DELFIN_CLASSIFIER": "off"})
    t0 = time.monotonic()
    _start(session, "do-bash sleep 5")
    # Sampled every half second for as long as the call can take: the
    # first request alone has taken 3.5 s on a loaded host, so fixed
    # moments early in the turn may all fall before the call starts.
    rows = []
    at = 0.5
    while at <= 15.0:
        rows.append(_sample_until(session, at, t0))
        if rows[-1] and "running bash" not in rows[-1] and any(
                "running bash" in r for r in rows):
            break
        at += 0.5
    running = [r for r in rows if "running bash: sleep 5" in r]
    assert len(running) >= 2, "no running-tool row:\n" + "\n".join(rows)
    secs = [_seconds(r, "running bash: sleep 5") for r in running]
    assert secs[-1] > secs[0], running
    assert session.wait_for("idle", 30) == "idle"


def test_an_approval_says_the_turn_waits_for_you(chat):
    session = chat()
    _start(session, "do-bash rm -f nothing.txt")
    assert session.wait_for("approval", 30) == "approval", (
        "\n".join(session.visible()))
    assert "waiting for your approval" in session.text(), (
        "\n".join(session.visible()))
    session.send("n", 0.5)
    assert session.wait_for("idle", 30) == "idle"


def test_a_failed_turn_leaves_no_indicator_behind(chat):
    session = chat()
    _start(session, "do-fail401")
    assert session.wait_for("idle", 60) == "idle", "\n".join(session.visible())
    bottom = "\n".join(session.visible()[-6:])
    assert "waiting for fake-model" not in bottom, bottom
    assert "esc to interrupt" not in bottom, bottom
