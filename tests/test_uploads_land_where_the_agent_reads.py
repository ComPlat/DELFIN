"""An attached file must land where the agent is allowed to read it.

The destination was chosen when the file was dropped. Before the first
message of a session there is no engine, so it fell back to ctx.repo_dir
-- the DELFIN checkout -- and in office mode that folder is not even
readable. Field case 20260805-105640: the UI reported

    📎 1 file(s) saved ... /pfs/.../software/delfin/.delfin/uploads/x.pdf

and the session's single tool call came back

    read denied: ... is outside /pfs/.../office ... no confirmation can
    grant it

Every upload made before the first message failed that way, and the agent
spent a turn discovering it. Buffering moves the decision to the moment
the workspace is known, and the target is then verified against the roots
the permissions actually grant rather than assumed to be reachable.

The fix deliberately does NOT widen the office root or unlock the scope.
That would trade a correct containment property for a UI convenience.
"""

from __future__ import annotations

import pathlib
import re

_SOURCE = pathlib.Path(
    __import__("delfin.dashboard.tab_agent", fromlist=["x"]).__file__
).read_text(encoding="utf-8")


def _block(start: str, end: str) -> str:
    return _SOURCE[_SOURCE.index(start):_SOURCE.index(end)]


def _code_only(text: str) -> str:
    """The block with its docstring and comments removed.

    A test that greps raw source matches the prose explaining a bug as
    readily as the bug itself -- this one did, on the very comment that
    records why the old destination was wrong.
    """
    quote_markers = ('"' * 3, "'" * 3)
    out: list[str] = []
    in_doc = False
    for line in text.splitlines():
        stripped = line.strip()
        opener = next((q for q in quote_markers if stripped.startswith(q)), None)
        if opener is not None:
            one_liner = len(stripped) > 3 and stripped.endswith(opener)
            if not one_liner:
                in_doc = not in_doc
            continue
        if in_doc or stripped.startswith("#"):
            continue
        out.append(line.split("#", 1)[0])
    return "\n".join(out)


# ---------------------------------------------------------------------------
# The decision moved
# ---------------------------------------------------------------------------

def test_the_drop_handler_no_longer_picks_a_folder():
    """It cannot: at drop time the workspace is not known yet."""
    handler = _code_only(
        _block("def _on_image_upload", "def _materialise_uploads"))
    assert "ctx.repo_dir" not in handler, (
        "the drop handler guesses a destination again")
    assert "write_bytes" not in handler, (
        "the drop handler writes again, before the workspace is known")


def test_the_bytes_are_buffered_at_drop_time():
    handler = _block("def _on_image_upload", "def _materialise_uploads")
    assert 'state["_pending_uploads"] = merged' in handler, (
        "buffered at drop time -- merged with what is already waiting")


def test_writing_happens_at_send_time():
    assert "_materialise_uploads(engine)" in _SOURCE, (
        "uploads are no longer written when the message is sent, so the "
        "workspace is guessed again")


# ---------------------------------------------------------------------------
# The destination is verified, not assumed
# ---------------------------------------------------------------------------

def test_the_target_is_checked_against_the_permitted_roots():
    writer = _block("def _materialise_uploads", "image_upload.observe")
    assert "find_readable_root_for" in writer, (
        "the upload target is assumed reachable again")


def test_an_unreachable_target_is_reported_instead_of_promised():
    """An upload that cannot be read is worse than one that was refused:
    the agent spends a turn discovering it."""
    writer = _block("def _materialise_uploads", "image_upload.observe")
    assert "was NOT attached" in writer
    assert "continue" in writer, "an unreachable file must be skipped"


def test_the_fix_does_not_widen_the_office_root():
    """The tempting shortcut -- add the uploads dir as an extra root, or
    drop the lock -- would trade a containment property for convenience."""
    writer = _block("def _materialise_uploads", "image_upload.observe")
    for shortcut in ("add_extra_dir", "lock_workspace", "scope_locked = False"):
        assert shortcut not in writer, shortcut


# ---------------------------------------------------------------------------
# Bounds
# ---------------------------------------------------------------------------

def test_the_whole_batch_is_capped_not_just_each_file():
    """Buffering means the bytes sit in memory until the next send; the
    per-file cap bounded one drop, not fifty."""
    assert "_UPLOAD_BUFFER_CAP" in _SOURCE
    handler = _block("def _on_image_upload", "def _materialise_uploads")
    assert "_UPLOAD_BUFFER_CAP" in handler


def test_the_per_file_cap_survives():
    handler = _block("def _on_image_upload", "def _materialise_uploads")
    assert "_UPLOAD_SIZE_CAP" in handler


def test_the_buffer_is_cleared_after_writing():
    """A second send must not re-attach the same files."""
    writer = _block("def _materialise_uploads", "image_upload.observe")
    assert re.search(r'state\["_pending_uploads"\] = \[\]', writer)


def test_no_stale_state_key_remains():
    """The old key held Paths, the new one holds bytes. A leftover reader
    would silently see the wrong shape."""
    assert "_pending_images" not in _SOURCE



# ---------------------------------------------------------------------------
# the folder is one a listing shows
# ---------------------------------------------------------------------------

def test_the_visible_folder_comes_first_and_the_hidden_one_stays_as_fallback(tmp_path):
    """A GLM session was handed the absolute path of an attachment under
    .delfin/uploads, explored the workspace with a listing instead, and
    found nothing: hidden folders are not listed (2026-09-11)."""
    from delfin.dashboard.tab_agent import _upload_dir_candidates
    cands = _upload_dir_candidates(tmp_path / "agent_workspace", tmp_path / "ws")
    assert [str(c.relative_to(tmp_path)) for c in cands] == [
        "agent_workspace/uploads", "ws/.delfin/uploads"]
    assert _upload_dir_candidates("", tmp_path / "ws") == [tmp_path / "ws" / ".delfin" / "uploads"]


def test_the_writer_takes_the_first_candidate_the_session_may_read():
    writer = _block("def _materialise_uploads", "image_upload.observe")
    assert "_upload_dir_candidates(" in writer
    assert "upload_dir = candidates[-1]" in writer, "the hidden folder is the fallback"
    assert "kp.find_readable_root_for(probe) is not None" in writer
    assert 'Path(ws) / ".delfin" / "uploads"' not in writer



# ---------------------------------------------------------------------------
# every attachment survives until it is sent
# ---------------------------------------------------------------------------

def test_a_second_drop_joins_the_first_instead_of_replacing_it():
    """Field report 2026-09-11: 07a attached, then 07 attached, one file
    written -- the second. The widget's value is the latest selection,
    and the buffer was replaced by it."""
    from delfin.dashboard.tab_agent import _merge_uploads
    first = [("07a_Auftrag.md", b"a")]
    second = [("07_Hintergrund.md", b"b")]
    assert _merge_uploads(first, second) == [("07a_Auftrag.md", b"a"), ("07_Hintergrund.md", b"b")]
    # the same name dropped again is the newer bytes, once
    assert _merge_uploads(first, [("07a_Auftrag.md", b"a2")]) == [("07a_Auftrag.md", b"a2")]
    assert _merge_uploads([], second) == second and _merge_uploads(first, []) == first


def test_the_drop_handler_merges_and_caps_the_whole_buffer():
    handler = _block("def _on_image_upload", "def _materialise_uploads")
    assert "_merge_uploads(state.get(\"_pending_uploads\") or [], buffered)" in handler
    assert 'state["_pending_uploads"] = merged' in handler
    assert 'state["_pending_uploads"] = buffered' not in handler
    assert "total = sum(len(c) for _n, c in merged)" in handler


def test_a_file_attached_mid_run_is_written_and_named_to_the_running_agent():
    """The note said the file would be written when the message is sent,
    and a mid-run message is one; the buffer was left for the next full
    send and the running agent looked for a file that was not there."""
    i = _SOURCE.index("Mid-loop steering: for the API/KIT/Ollama engine")
    body = _SOURCE[i:i + 2200]
    assert "_materialise_uploads(_seng)" in body
    assert "The user attached these files" in body
    assert "_seng.steer(_steer_text)" in body



def test_a_missing_attachment_is_to_be_said_not_searched():
    """A turn asked about an attachment that was not there spent 239k
    tokens and six and a half minutes on five tool rounds of ls and find
    (measured 2026-09-11). The note now says what to do instead."""
    assert _SOURCE.count("If one of these paths does not exist, say which") == 2, (
        "both the send-time note and the mid-run note carry the rule")
