"""Controls for append_write and the write_file(mode="append") wiring.

Red on the previous commit: ``delfin.agent.append_write`` does not
exist (Phase 3 of the LB assignment), and ``write_file`` has no
``mode`` argument, so an append call overwrites the file with only the
tail. The xfail(strict=True) integration test pins the handler
contract: append must go through the same journal / read-baseline /
diff path as a plain write, because undo and the changes report must
see a full before/after pair.
"""

import pytest

from delfin.agent import append_write
from delfin.agent.api_client import _DocToolExecutor, KitToolPermissions


def _perms(ws):
    perms = KitToolPermissions(workspace=str(ws))
    perms.mode = "acceptEdits"
    perms.task_session_id = "append-write-test"
    return perms


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "ws"
    d.mkdir()
    (d / "log.txt").write_text("alpha\nbeta", encoding="utf-8")
    return d


# ---------------------------------------------------------------------------
# build_new_text: the pure append transformation
# ---------------------------------------------------------------------------

def test_append_after_final_newline_is_plain_concatenation():
    assert append_write.build_new_text("a\n", "b\n") == "a\nb\n"


def test_append_inserts_the_missing_final_newline():
    # The fixture's shape: a file without a trailing newline. Appending
    # must not glue the tail onto the last line ("betaTAIL").
    assert append_write.build_new_text("alpha\nbeta", "gamma\n") \
        == "alpha\nbeta\ngamma\n"


def test_append_to_empty_text_is_the_tail():
    assert append_write.build_new_text("", "x\n") == "x\n"


def test_append_empty_tail_is_rejected():
    # An append that adds nothing is a no-op write; refusing it keeps
    # the journal free of changes that did not change anything.
    with pytest.raises(append_write.AppendRejected):
        append_write.build_new_text("a\n", "")


def test_append_tail_that_is_only_whitespace_is_rejected():
    with pytest.raises(append_write.AppendRejected):
        append_write.build_new_text("a\n", "   \n")


def test_new_file_append_is_allowed_and_is_just_the_tail():
    # mode="append" on a file that does not exist yet is the same as
    # creating it -- the heredoc-in-pieces pattern starts with file
    # creation, and refusing it would force a mode switch mid-stream.
    assert append_write.build_new_text(None, "first\n") == "first\n"


# ---------------------------------------------------------------------------
# Integration through the public call path (xfail until the handler is
# wired by the operator).
# ---------------------------------------------------------------------------

@pytest.mark.xfail(
    strict=True,
    reason="write_file does not accept mode=append yet; wiring lives in "
           "api_client.py (operator's area). See .gate/SCHEMA-APPEND.md",
)
def test_write_file_append_goes_through_the_write_path(ws):
    ex = _DocToolExecutor()
    perms = _perms(ws)
    # The read baseline an append needs: same rule as an overwrite.
    ex._execute_read_file({"path": "log.txt"}, perms)
    out = ex._execute_write_file(
        {"path": "log.txt", "content": "gamma\n", "mode": "append"}, perms)
    assert "overwritten" not in out
    assert (ws / "log.txt").read_text(encoding="utf-8") \
        == "alpha\nbeta\ngamma\n"
    # The journal record must be a full write record like any other:
    # a pre-image hash (of "alpha\nbeta"), a post hash (of the joined
    # text), created=False. This is what makes undo and the changes
    # report treat an append exactly like a plain write.
    from delfin.agent import change_journal as cj
    recs = cj._read_journal("append-write-test")
    recs = [r for r in recs if r.get("tool") == "write_file"]
    assert recs, "no write_file record journalled"
    last = recs[-1]
    assert last["created"] is False
    assert last["pre_hash"] == cj.sha256_text("alpha\nbeta")
    assert last["post_hash"] == cj.sha256_text("alpha\nbeta\ngamma\n")
