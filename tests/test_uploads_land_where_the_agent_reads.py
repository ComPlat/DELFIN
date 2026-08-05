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
    assert 'state["_pending_uploads"] = buffered' in handler


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
