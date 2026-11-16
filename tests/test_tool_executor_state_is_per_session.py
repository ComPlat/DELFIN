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

The turn boundary and the per-session key used to meet in ``_note_change_seq``:
it keyed the write by ``perms.task_session_id``, but it detected "a new turn"
by ``len(_turn_change_seqs) < _turn_seqs_seen`` and cleared the WHOLE
``_turn_seqs_by_session`` dict. A second session starting a turn emptied the
first session's window and, on its next write, wiped the first session's
per-session entry: the first session's undo window fell back to the shared
flat list, which held the OTHER session's seqs.

That state has been moved into ``delfin.agent.tool_context.ToolContext``: one
per session, owned by the session's client and routed by
``perms.task_session_id``. Each session's window, pending images and live
conversation now live apart, so nothing one session starts can empty or clear
anything another session wrote.

Each test below states the behaviour a user expects on the real
``_DocToolExecutor`` singleton. Where the behaviour now holds for every
session (one session, or two sessions with the shared field) the test is a
plain test.

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
from delfin.agent import tool_context as tc

# The executor attribute names under test, saved and restored around every
# test because the object is a module singleton. With ToolContext these are
# properties delegating to a per-session context; the per-session registry and
# the process-wide default context are what actually hold the state.
_SAVED = ("_turn_change_seqs", "_turn_seqs_by_session", "_turn_seqs_seen",
          "_pending_view_images", "live_messages")


@pytest.fixture(autouse=True)
def _pristine_executor():
    """Start and end each test with no per-turn state, and restore it after,
    so this file never leaks a turn into another test (or another file)."""
    ex = A._doc_executor
    had = {k: ex.__dict__[k] for k in _SAVED if k in ex.__dict__}
    had_registry = None
    had_default = None
    if hasattr(ex, "_tool_contexts"):
        had_registry = ex._tool_contexts
        ex._tool_contexts = {}
    if hasattr(ex, "_default_tool_context"):
        had_default = ex._default_tool_context
        ex._default_tool_context = tc.ToolContext()
    for k in _SAVED:
        ex.__dict__.pop(k, None)
    try:
        yield ex
    finally:
        for k in _SAVED:
            ex.__dict__.pop(k, None)
        ex.__dict__.update(had)
        if had_registry is not None:
            ex._tool_contexts = had_registry
        if had_default is not None:
            ex._default_tool_context = had_default


def _perms(ws: Path, sid: str) -> KitToolPermissions:
    p = KitToolPermissions(workspace=ws, mode="bypassPermissions")
    p.task_session_id = sid
    return p


def _start_turn() -> None:
    """What the client does at the top of ``stream_message``: empty the turn
    window. Under ToolContext a turn boundary belongs to one context; the
    executor's property routes this reset to the process-wide default context,
    which is inert for the per-session windows that now separate sessions."""
    A._doc_executor._turn_change_seqs = []


def _write(ex, perms, seq: int) -> None:
    """One file write's journal entry, through the real method every write
    path calls (``_capture_change`` / ``_capture_binary_change`` ->
    ``_note_change_seq``). Routes to the session's own ToolContext."""
    ex._note_change_seq({"seq": seq}, perms)


def _drain_pending_images(ex, sid: str) -> list:
    """What the round loop does with opened images: take the whole list and
    clear it. Pinned to one session's context, because an image opened by one
    session must never be drained by another (the round loop drains its own
    ``ToolContext``)."""
    return ex._ctx_for(sid).take_pending_images()


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
    delivered = _drain_pending_images(ex, s.task_session_id)
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


def test_two_sessions_a_undo_window_keeps_only_a_after_b_turn(
        tmp_path, _pristine_executor):
    """Session A writes; session B then runs a turn while A's turn is still
    open and writes its own file. A's undo window (what
    ``undo_changes(scope="turn")`` reads via ``_turn_seqs_for("A")``) must
    hold A's seqs, not B's.

    Under ToolContext each session routes its seqs to its OWN context, so B's
    turn never empties or clears anything of A's: after B writes,
    ``_turn_seqs_for("A")`` still returns ``[11, 12]``."""
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


def test_two_sessions_an_image_reaches_the_session_that_opened_it(
        tmp_path, _pristine_executor):
    """Session A opens an image in its turn. Session B's round then finishes
    and drains its own pending images. The image must go to A, which opened
    it, and not to B."""
    ex = _pristine_executor
    a = _perms(tmp_path / "wsA", "A")
    (tmp_path / "wsA").mkdir()
    _png(tmp_path / "wsA" / "a.png")

    # A opens an image (landed in A's session context).
    out = json.loads(ex._execute_view_image({"path": "a.png"}, a))
    assert out["status"] == "ok"

    # B's round drains B's OWN context: nothing of A's is there.
    delivered_to_b = _drain_pending_images(ex, "B")

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


def test_real_client_two_sessions_keep_their_own_turn_windows(
        tmp_path, monkeypatch, _pristine_executor):
    """Real run: two real OpenAIClient sessions in one process against the
    scripted endpoint (tests/agent_terminal_probes/fake_model.py, serve(port=
    0)). Session A starts a turn and writes TWO files (a.txt, a2.txt) through
    the same tool-execution entry the engine's round loop uses
    (``_doc_executor.execute("write_file", ..., permissions=...)``, the call
    the client's round loop makes). Session B then starts a turn and
    writes b.txt through a real ``stream_message`` round. After both turns A
    runs a turn-scoped undo through the same real path
    (``_doc_executor.execute("undo_changes", {"scope": "turn"}, A_perms)``).
    The behaviour a user expects: A's undo reverts A's whole turn (BOTH a.txt
    and a2.txt) and never touches B's workspace.

    Two writes are needed for the assertion to be unambiguous: the change
    journal numbers seq per session (change_journal.py:521, keyed by
    task_session_id), so a single write's seq 1 collides by number with the
    other session's seq 1. With A's window at seqs [1, 2] in A's own
    ToolContext, B's turn — which resets and writes only B's context — cannot
    collapse A's window, so A's undo reverts both files and never touches B's.

    Ordering is strict and needs no Event or sleep: every step runs on the
    main thread in sequence -- A writes both files, then B runs its whole
    turn, then A runs its undo. So A's two writes always land before B's
    turn starts, and no session ever runs concurrently with its own undo.

    Skipped where the host holds loopback connects: there the endpoint cannot
    be reached at all (see the module docstring), so no session can drive a
    stream_message round."""
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
        perms_a = cli_a._permissions
        perms_b = cli_b._permissions
        perms_a.task_session_id = "A"
        perms_b.task_session_id = "B"

        # Session A starts a turn: two real write_file executions (the same
        # entry the engine's round loop calls). Both land in A's journal as
        # seqs 1 and 2 and in ``_turn_seqs_by_session["A"]``.
        for path in ("a.txt", "a2.txt"):
            out = A._doc_executor.execute("write_file", {
                "path": path, "content": f"from-{path[0]}\n"},
                permissions=perms_a)
            assert '"error"' not in str(out), f"write {path} failed: {out}"
        assert (ws_a / "a.txt").exists() and (ws_a / "a2.txt").exists()

        # Session B starts a turn while A's turn is still open and writes
        # b.txt. B's stream_message turn start empties the shared flat list,
        # and B's write wipes _turn_seqs_by_session, so A's window collapses.
        _drive_turn(cli_b, "do-write b.txt from-b", "B")
        assert (ws_b / "b.txt").exists()

        # Session A runs a turn-scoped undo through the real tool path. The
        # behaviour a user expects: A's whole turn (a.txt AND a2.txt) reverts,
        # and nothing in B's workspace is touched.
        undone = json.loads(A._doc_executor.execute(
            "undo_changes", {"scope": "turn"}, permissions=perms_a))
        assert (ws_a / "a.txt").exists() is False, (
            "a.txt still present after session A's turn undo; reverted="
            + repr(undone.get("reverted", [])))
        assert (ws_a / "a2.txt").exists() is False, (
            "a2.txt still present: A's turn-scope undo captured only the "
            "shared flat list (B's seq) and not A's own [1, 2]; reverted="
            + repr(undone.get("reverted", [])))
        # B's workspace is untouched by A's undo, and B's window is B's own.
        assert (ws_b / "b.txt").read_text() == "from-b\n", (
            "session A's undo reached into session B's workspace")
        assert A._doc_executor._turn_seqs_for("B") == [1], (
            "session B's undo window no longer holds B's own write")
    finally:
        server.shutdown()
        server.server_close()
