"""The crash checkpoint stored pre-redaction text and replayed it.

The checkpoint is written from inside the tool loop -- every ten tool
events or every sixty seconds -- carrying the turn's user message and the
partial answer so far. The output guard runs at the END of the turn, so
everything the checkpoint holds reaches disk before anything has looked at
it.

The file itself is 0600, so this is not primarily a permissions problem.
It is a bypass and a replay: on the next start ``consume_crash_recovery
note`` appends the stored partial straight into the message list, where it
re-enters the model context and is written again by the next full session
save -- which is the file a bug report bundles and ships.

Redacting in the writer guards both at once, which is why the fix belongs
there rather than at either use site.
"""

from __future__ import annotations

import hashlib
import string

import pytest

from delfin.agent import session_store


def _token(n: int, seed: str) -> str:
    out: list[str] = []
    h = seed.encode()
    alphabet = string.ascii_letters + string.digits
    while len(out) < n:
        h = hashlib.sha256(h).digest()
        out += [alphabet[b % len(alphabet)] for b in h]
    return "".join(out[:n])


@pytest.fixture
def sessions(tmp_path, monkeypatch):
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", tmp_path / "sessions")
    return tmp_path / "sessions"


# ---------------------------------------------------------------------------
# The disk
# ---------------------------------------------------------------------------

def test_a_credential_in_the_partial_answer_never_reaches_the_file(sessions):
    leaked = _token(48, "partial")
    p = session_store.save_turn_checkpoint("s1", {
        "user_message": "why is the request rejected?",
        "partial_response": f"the call sends OPENAI_API_KEY={leaked} and 401s",
        "tool_calls": 4,
    })
    assert p is not None
    text = p.read_text(encoding="utf-8")
    assert leaked not in text
    assert "redacted" in text


def test_a_credential_in_the_user_message_never_reaches_the_file(sessions):
    """The user pasting their own key into the question is the ordinary
    way this happens, not the exotic one."""
    leaked = "sk-ant-" + _token(40, "user")
    p = session_store.save_turn_checkpoint("s2", {
        "user_message": f"this key does not work: {leaked}",
        "partial_response": "checking",
    })
    assert leaked not in p.read_text(encoding="utf-8")


def test_a_credential_nested_in_the_payload_is_found(sessions):
    """The payload is a dict the engine builds; a later field carrying a
    structure must not walk past the guard."""
    leaked = _token(48, "nested")
    p = session_store.save_turn_checkpoint("s3", {
        "user_message": "x",
        "extra": {"headers": [f"Authorization: Bearer {leaked}"]},
    })
    assert leaked not in p.read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# ...and the replay
# ---------------------------------------------------------------------------

def test_the_recovery_note_does_not_re_inject_the_credential(sessions):
    """The note is appended to the message list on the next start, so an
    unredacted partial would re-enter the model context and be written
    again by the next full session save."""
    leaked = _token(48, "replay")
    session_store.save_turn_checkpoint("s4", {
        "user_message": "debug the auth failure",
        "partial_response": f"found it: api_key={leaked}",
        "tool_calls": 7,
        "ts": 9_999_999_999.0,        # newer than any saved session
    })
    note = session_store.consume_crash_recovery_note("s4")
    assert note, "no recovery note was produced"
    assert leaked not in note
    assert "[recovered]" in note


def test_the_recovery_note_still_says_what_was_interrupted(sessions):
    """A note stripped of its content would be a different bug: its whole
    job is telling the next turn what was in flight."""
    session_store.save_turn_checkpoint("s5", {
        "user_message": "rebuild the ligand table",
        "partial_response": "I had written three rows",
        "tool_calls": 3,
        "ts": 9_999_999_999.0,
    })
    note = session_store.consume_crash_recovery_note("s5")
    assert "rebuild the ligand table" in note
    assert "3 tool calls" in note
    assert "I had written three rows" in note


def test_ordinary_progress_text_is_left_alone(sessions):
    """A redactor that mangles the partial answer costs more than it
    saves -- the note exists to be read."""
    text = ("Ran xtb on 3f2a1b9c4d5e6f7a8b9c0d1e2f3a4b5c6d7e8f90, "
            "got -12.4 kcal/mol for CC(=O)Oc1ccccc1C(=O)O")
    p = session_store.save_turn_checkpoint("s6", {
        "user_message": "optimise aspirin",
        "partial_response": text,
    })
    assert text in p.read_text(encoding="utf-8")


def test_the_checkpoint_survives_a_broken_redactor(sessions, monkeypatch):
    """A checkpoint that fails to save costs a crashed turn's progress, so
    the guard stays best-effort."""
    import delfin.agent.output_guard as og

    def boom(*a, **kw):
        raise RuntimeError("redactor exploded")

    monkeypatch.setattr(og, "_redact_secrets", boom)
    p = session_store.save_turn_checkpoint("s7", {
        "user_message": "keep going", "partial_response": "half done"})
    assert p is not None and p.is_file()
    assert "half done" in p.read_text(encoding="utf-8")


def test_the_checkpoint_file_is_owner_only(sessions):
    import os
    import stat

    p = session_store.save_turn_checkpoint("s8", {"user_message": "x"})
    assert stat.S_IMODE(os.stat(p).st_mode) == 0o600
