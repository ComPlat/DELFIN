"""Per-turn tool-executor state belongs to ONE session's turn, not to the
process.

``_doc_executor`` is a module-level singleton (``delfin/agent/api_client.py``
``_doc_executor = _DocToolExecutor()``) shared by every session in a process:
the dashboard's several sessions, and a parent agent with a background
sub-agent. Three pieces of state that describe a single turn live on it:

* ``_turn_change_seqs`` / ``_turn_seqs_by_session`` — the undo journal seqs
  recorded this turn, read by ``_turn_seqs_for(session_id)`` and by
  ``undo_changes(scope="turn")``.
* ``_pending_view_images`` — images a ``view_image`` call opened, drained by
  the next model round that finishes.
* ``live_messages`` — the conversation ``history_search`` resolves ``live:``
  refs against.

The turn boundary and the per-session key meet in ``_note_change_seq``
(api_client.py:16551-16563): it keys the write by
``perms.task_session_id``, but it detects "a new turn" by
``len(_turn_change_seqs) < _turn_seqs_seen`` and then clears the WHOLE
``_turn_seqs_by_session`` dict. The client marks the new turn by emptying the
flat list at the top of ``stream_message`` (api_client.py:33336). So a second
session starting a turn empties the first session's window and, on that
session's next write, wipes the first session's per-session entry: the first
session's undo window then falls back to the shared flat list, which holds
the OTHER session's seqs.

Each test below states the behaviour a user would expect and marks the ones
today's code does not provide ``xfail(strict=True)`` with the shared field and
its line. Where today's code already behaves (one session, or the
single-session control) the test is a plain test.

Determinism: the state tests drive the real ``_DocToolExecutor`` object with
the real ``_note_change_seq`` / ``_turn_seqs_for`` / ``_execute_view_image``
methods and the exact ``_turn_change_seqs = []`` reset the client performs at
turn start. No sleep, no threads, no network — the interleave is stated by
calling the transitions in the order the client would. The real-client test at
the bottom drives two ``OpenAIClient`` turns against the scripted endpoint
(``tests/agent_terminal_probes/fake_model.py``, ``serve(port=0)``); it is
gated on the host permitting a loopback connect to the endpoint, because some
sandboxes hold every loopback connect (the same condition
``tests/test_the_first_remember_on_a_real_pty.py`` guards) and a test of the
sessions cannot say anything on a host where no session can reach the
endpoint at all.
"""
from __future__ import annotations

import json
import threading
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions

# The executor attribute names under test, saved and restored around every
# test because the object is a module singleton.
_SAVED = ("_turn_change_seqs", "_turn_seqs_by_session", "_turn_seqs_seen",
          "_pending_view_images", "live_messages")


@pytest.fixture(autouse=True)
def _pristine_executor():
    """Start and end each test with no per-turn state, and restore it after,
    so this file never leaks a turn into another test (or another file)."""
    ex = A._doc_executor
    had = {k: ex.__dict__[k] for k in _SAVED if k in ex.__dict__}
    for k in _SAVED:
        ex.__dict__.pop(k, None)
    try:
        yield ex
    finally:
        for k in _SAVED:
            ex.__dict__.pop(k, None)
        ex.__dict__.update(had)


def _perms(ws: Path, sid: str) -> KitToolPermissions:
    p = KitToolPermissions(workspace=ws, mode="bypassPermissions")
    p.task_session_id = sid
    return p


def _start_turn() -> None:
    """What the client does at the top of ``stream_message``
    (api_client.py:33336): empty the flat turn window."""
    A._doc_executor._turn_change_seqs = []


def _write(ex, perms, seq: int) -> None:
    """One file write's journal entry, through the real method every write
    path calls (``_capture_change`` / ``_capture_binary_change`` ->
    ``_note_change_seq``)."""
    ex._note_change_seq({"seq": seq}, perms)


def _drain_pending_images(ex) -> list:
    """What the round loop does with opened images (api_client.py:
    35091-35093): take the whole list and clear it."""
    imgs = getattr(ex, "_pending_view_images", None)
    if not imgs:
        return []
    ex._pending_view_images = []
    return list(imgs)


def _png(path: Path) -> None:
    Image = pytest.importorskip("PIL.Image")
    Image.new("RGB", (8, 8), "red").save(path)


# ---------------------------------------------------------------------------
# Single-session control: the machinery itself works.
# ---------------------------------------------------------------------------


def test_control_one_session_undo_window_is_its_own(tmp_path, _pristine_executor):
    """One session, same steps as the two-session test: its turn window is
    exactly the seqs it wrote, and its undo reads them. This shows the test
    machinery works, so a failing two-session test is about sharing."""
    ex = _pristine_executor
    s = _perms(tmp_path, "solo")
    _start_turn()
    _write(ex, s, 11)
    _write(ex, s, 12)
    assert ex._turn_seqs_for("solo") == [11, 12]


def test_control_one_session_image_reaches_its_own_round(tmp_path,
                                                        _pristine_executor):
    """One session: an image it opens in its turn is delivered to that turn's
    next round."""
    ex = _pristine_executor
    _png(tmp_path / "shot.png")
    s = _perms(tmp_path, "solo")
    out = json.loads(ex._execute_view_image({"path": "shot.png"}, s))
    assert out["status"] == "ok"
    delivered = _drain_pending_images(ex)
    assert [img.source_path.name for img in delivered] == ["shot.png"]


# ---------------------------------------------------------------------------
# Two sessions, undo window: A's undo must see only A's writes.
# ---------------------------------------------------------------------------


def test_two_sessions_a_window_is_correct_before_b_writes(
        tmp_path, _pristine_executor):
    """Session A writes inside its own turn; its undo window holds A's seqs.
    No other session has started, so today's code does this right — it pins
    the healthy half and the machinery."""
    ex = _pristine_executor
    a = _perms(tmp_path / "wsA", "A")
    (tmp_path / "wsA").mkdir()

    _start_turn()
    _write(ex, a, 11)
    _write(ex, a, 12)
    assert ex._turn_seqs_for("A") == [11, 12]


@pytest.mark.xfail(
    strict=True,
    reason="_turn_change_seqs/_turn_seqs_by_session are one flat field pair on "
           "the process-wide _doc_executor: B's turn start empties "
           "_turn_change_seqs (api_client.py:33336) and B's next "
           "_note_change_seq clears all of _turn_seqs_by_session "
           "(api_client.py:16558-16559), so A's window falls back to B's seq",
)
def test_two_sessions_a_undo_window_keeps_only_a_after_b_turn(
        tmp_path, _pristine_executor):
    """Session A writes; session B then runs a turn while A's turn is still
    open and writes its own file. A's undo window (what
    ``undo_changes(scope="turn")`` reads via ``_turn_seqs_for("A")``) must
    hold A's seqs, not B's.

    Today's code does not: B's turn start empties the shared flat
    ``_turn_change_seqs``, and B's next write sees
    ``len(_turn_change_seqs) < _turn_seqs_seen`` and clears the whole
    ``_turn_seqs_by_session`` dict, so A's window falls back to the flat
    list, which holds B's seq. Measured on the base commit: after B writes,
    ``A._turn_seqs_for("A")`` returns ``[21]`` (B's seq), not ``[11, 12]``."""
    ex = _pristine_executor
    a = _perms(tmp_path / "wsA", "A")
    b = _perms(tmp_path / "wsB", "B")
    (tmp_path / "wsA").mkdir()
    (tmp_path / "wsB").mkdir()

    # A's turn: two writes.
    _start_turn()
    _write(ex, a, 11)
    _write(ex, a, 12)

    # B starts a turn while A's turn is still open, then writes.
    _start_turn()
    _write(ex, b, 21)

    a_window_after_b = ex._turn_seqs_for("A")
    # What a user expects: A's undo reverts A's two writes and nothing of B's.
    assert a_window_after_b == [11, 12], (
        "session A's turn-undo window after session B ran a turn; it must "
        f"still hold A's own seqs [11, 12], not B's; got {a_window_after_b}"
    )


def test_two_sessions_b_undo_window_is_its_own(tmp_path, _pristine_executor):
    """The other side of the same coin, which today's code DOES get right:
    after B's write clears the dict, B's own window holds B's seq. This pins
    that the clearing is symmetric — it is not that A is privileged."""
    ex = _pristine_executor
    a = _perms(tmp_path / "wsA", "A")
    b = _perms(tmp_path / "wsB", "B")
    (tmp_path / "wsA").mkdir()
    (tmp_path / "wsB").mkdir()

    _start_turn()
    _write(ex, a, 11)
    _start_turn()
    _write(ex, b, 21)

    assert ex._turn_seqs_for("B") == [21]


# ---------------------------------------------------------------------------
# Two sessions, view_image: A's image must not be delivered to B's round.
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    reason="_pending_view_images is one flat list on the process-wide "
           "_doc_executor (appended at api_client.py:19301-19304, drained at "
           "35091-35093 with no session key): an image session A opened is "
           "handed to whichever session's round drains the list first",
)
def test_two_sessions_an_image_reaches_the_session_that_opened_it(
        tmp_path, _pristine_executor):
    """Session A opens an image in its turn. Session B's round then finishes
    and drains the pending images. The image must go to A, which opened it,
    and not to B."""
    ex = _pristine_executor
    a = _perms(tmp_path / "wsA", "A")
    (tmp_path / "wsA").mkdir()
    _png(tmp_path / "wsA" / "a.png")

    # A opens an image.
    out = json.loads(ex._execute_view_image({"path": "a.png"}, a))
    assert out["status"] == "ok"

    # B's round drains the shared list.
    delivered_to_b = _drain_pending_images(ex)

    assert delivered_to_b == [], (
        "session B's round was handed the image session A opened: "
        f"{[img.source_path.name for img in delivered_to_b]}"
    )


# ---------------------------------------------------------------------------
# Real client run: two OpenAIClient turns in one process against the
# scripted OpenAI-compatible endpoint (tests/agent_terminal_probes/
# fake_model.py). Gated on the host permitting a loopback connect to the
# endpoint; some sandboxes hold every loopback connect, and on such a host
# no session can reach the endpoint (tests/test_the_first_remember_on_a_
# real_pty.py guards the same condition).
# ---------------------------------------------------------------------------


def _loopback_reachable(port: int) -> bool:
    import socket
    try:
        with socket.create_connection(("127.0.0.1", port), timeout=3):
            return True
    except OSError:
        return False


def _fake_caps(monkeypatch):
    """Resolve model capabilities to a fixed record, so a turn does not probe
    the endpoint's /v1/models over the network."""
    from delfin.agent import model_capabilities as mc
    from tests.agent_terminal_probes import fake_model
    caps = mc.ModelCapabilities(
        model=fake_model.MODEL, provider="ollama",
        context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)


def _drive_turn(client, prompt: str, sid: str,
                started: threading.Event | None = None) -> list:
    """One real turn: drain stream_message, record every event."""
    if started is not None:
        started.set()
    events = []
    for ev in client.stream_message(
            system="you are a session in a two-session test",
            messages=[{"role": "user", "content": prompt}],
            max_tokens=8192, session_id=sid):
        events.append(ev)
    return events


@pytest.mark.xfail(
    strict=True,
    reason="same shared _turn_change_seqs/_turn_seqs_by_session on "
           "_doc_executor as the state tests: with A and B running turns in "
           "one process, B's turn start empties the shared flat list "
           "(api_client.py:33336) and B's next write clears all of "
           "_turn_seqs_by_session (api_client.py:16558-16559), so A's window "
           "falls back to the flat list that holds B's seq and the two "
           "windows overlap",
)
def test_real_client_two_sessions_keep_their_own_turn_windows(
        tmp_path, monkeypatch, _pristine_executor):
    """Real run: two real OpenAIClient sessions in one process, each driving
    a ``do-write`` turn against the scripted endpoint. Session A writes, then
    session B starts a turn while A's turn is still open (released by a
    threading.Event in the endpoint, never by a sleep), and A's undo window
    must still hold A's write.

    Skipped where the host holds loopback connects: there the endpoint cannot
    be reached at all (see the module docstring)."""
    fake_model = pytest.importorskip("tests.agent_terminal_probes.fake_model")

    server, port = fake_model.serve(port=0,
                                    log_path=str(tmp_path / "model.log"))
    try:
        if not _loopback_reachable(port):
            pytest.skip("this host's sandbox holds loopback connects, so a "
                        "session cannot reach the scripted endpoint")
        # The endpoint is a loopback URL the client's httpx must reach
        # directly: the harness exports its egress proxy in HTTP(S)_PROXY,
        # which refuses a non-80/443 port.
        monkeypatch.setenv("NO_PROXY", "127.0.0.1,localhost,::1")
        monkeypatch.setenv("no_proxy", "127.0.0.1,localhost,::1")
        monkeypatch.setenv("OLLAMA_HOST", f"http://127.0.0.1:{port}")
        _fake_caps(monkeypatch)

        ws_a = tmp_path / "wsA"
        ws_b = tmp_path / "wsB"
        ws_a.mkdir()
        ws_b.mkdir()
        cli_a = A.create_client(backend="api", provider="ollama",
                                model=fake_model.MODEL, cwd=str(ws_a))
        cli_b = A.create_client(backend="api", provider="ollama",
                                model=fake_model.MODEL, cwd=str(ws_b))
        cli_a._permissions.task_session_id = "A"
        cli_b._permissions.task_session_id = "B"

        # A's turn runs while B's turn starts: A signals when it has been
        # handed to the endpoint, B starts, then A finishes. The Event is set
        # by A's thread the instant its request begins, so the order is
        # deterministic without a sleep.
        a_started = threading.Event()
        errors: list = []

        def run_a():
            try:
                _drive_turn(cli_a, "do-write a.txt from-a", "A",
                            started=a_started)
            except BaseException as exc:        # noqa: BLE001 - reported
                errors.append(("A", exc))

        t = threading.Thread(target=run_a, daemon=True)
        t.start()
        assert a_started.wait(30), "session A's turn never started"
        _drive_turn(cli_b, "do-write b.txt from-b", "B")
        t.join(60)
        assert not t.is_alive(), "session A's turn did not finish"
        assert not errors, f"a session's turn raised: {errors!r}"

        # A wrote a.txt and B wrote b.txt, in their own workspaces.
        assert (ws_a / "a.txt").exists(), "session A's write did not land"
        assert (ws_b / "b.txt").exists(), "session B's write did not land"

        # Each session's undo window holds that session's own write seq and
        # not the other session's. The two windows share no seq.
        a_window = A._doc_executor._turn_seqs_for("A")
        b_window = A._doc_executor._turn_seqs_for("B")
        assert a_window, "session A's undo window is empty"
        assert b_window, "session B's undo window is empty"
        assert set(a_window).isdisjoint(set(b_window)), (
            "session A and session B undo windows overlap: A="
            f"{a_window} B={b_window}"
        )
    finally:
        server.shutdown()
        server.server_close()
