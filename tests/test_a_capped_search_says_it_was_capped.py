"""A search that stopped at its cap said nothing about having stopped.

Measured through the executor's own dispatch (the path a model's call
takes), on a temporary workspace, before the fix:

    grep_file(pattern="hit", path="many.txt")   -> 30 lines, last
                                                  "many.txt:30: hit 30"
    list_files(pattern="*.txt")                 -> 50 entries, last
                                                  "f48.txt"
    read_file(path="five.txt", offset=6)        -> ""

The file held 59 hits, the folder 60 files, the file 5 lines. In every
case the answer was byte-indistinguishable from a search that had found
everything there was, a folder holding nothing else, and an empty file.

That is not a cosmetic silence. The engine records a result as cut only
when the text says so (``api_client.py`` sets ``_tool_cut`` from
"truncated," / "[truncated:" / ``"truncated": true``, and
``engine.py:4275`` repeats the test), so a count taken from a capped
``grep_file`` reached ``verify_guard.record_tool_numbers`` as a
MEASUREMENT. ``read_file`` has stated its cut since it was written
("... (N lines total, showing A-B)"); the two search tools never did.

An empty result was the same defect in its cheapest form: a read past the
end of the file and an empty file both returned "". One is a wrong offset
the model must fix, the other is the file's content.

``assert_the_note_is_not_a_hit`` below is the assertion that could have
failed the wrong way: the note is appended to a result whose every line is
parsed as ``path:line:`` by the observed-files ledger, so a note shaped
like a hit would invent a file the agent never read.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client
from delfin.agent.api_client import KitToolPermissions, _paths_in_grep_output

GREP_CAP = 30
LIST_CAP = 50


@pytest.fixture
def box(tmp_path):
    (tmp_path / "many.txt").write_text(
        "".join(f"hit {i}\n" for i in range(1, 60)), encoding="utf-8")
    (tmp_path / "few.txt").write_text("hit a\nhit b\nhit c\n",
                                      encoding="utf-8")
    (tmp_path / "five.txt").write_text("l1\nl2\nl3\nl4\nl5\n",
                                       encoding="utf-8")
    (tmp_path / "empty.txt").write_text("", encoding="utf-8")
    return tmp_path


@pytest.fixture
def call(box):
    perms = KitToolPermissions(mode="bypassPermissions", workspace=str(box))

    def _call(tool, args):
        return str(api_client._doc_executor.execute(tool, args, perms))

    return _call


# ---------------------------------------------------------------------------
# grep_file: the cap it never mentioned
# ---------------------------------------------------------------------------

def test_a_grep_that_stopped_at_the_cap_says_so(call):
    out = call("grep_file", {"pattern": "hit", "path": "many.txt"})
    hits = [ln for ln in out.splitlines() if ln.startswith("many.txt:")]
    assert len(hits) == GREP_CAP, "the cap itself is not what this pins"
    assert f"stopped at max_results={GREP_CAP}" in out, out
    assert "more may exist" in out, out


def test_a_grep_that_found_everything_says_nothing_extra(call):
    """The control: a result that fits must not grow a note, or the note
    stops meaning 'you were shown less than there is'."""
    out = call("grep_file", {"pattern": "hit", "path": "few.txt"})
    assert len(out.splitlines()) == 3, out
    assert "stopped at" not in out, out


def test_no_matches_is_still_not_a_cut(call):
    out = call("grep_file", {"pattern": "zzz-not-there", "path": "."})
    assert out == "No matches found.", out


def test_the_note_is_not_a_hit(call):
    """Every line of a grep result is parsed as ``path:line:`` to decide
    which files the agent demonstrably saw. A note that looks like a hit
    would put a file in the ledger that was never shown."""
    out = call("grep_file", {"pattern": "hit", "path": "many.txt"})
    assert _paths_in_grep_output(out) == {"many.txt"}, out


# ---------------------------------------------------------------------------
# list_files: the same cap, silently
# ---------------------------------------------------------------------------

def test_a_listing_that_stopped_at_the_cap_says_so(tmp_path, call):
    for i in range(60):
        (tmp_path / f"f{i:02d}.txt").write_text("x\n", encoding="utf-8")
    out = call("list_files", {"pattern": "*.txt"})
    entries = [ln for ln in out.splitlines() if ln.endswith(".txt")]
    assert len(entries) == LIST_CAP, "the cap itself is not what this pins"
    assert f"stopped at {LIST_CAP} entries" in out, out


def test_a_listing_that_fits_says_nothing_extra(box, call):
    out = call("list_files", {"pattern": "*.txt"})
    assert "stopped at" not in out, out


# ---------------------------------------------------------------------------
# read_file: past the end of the file
# ---------------------------------------------------------------------------

def test_a_read_past_the_end_names_the_line_count(call):
    out = call("read_file", {"path": "five.txt", "offset": 6})
    assert out.lstrip().startswith('{"error"'), out
    assert "5 lines" in out, out


def test_a_read_inside_the_file_is_unchanged(call):
    out = call("read_file", {"path": "five.txt", "offset": 5})
    assert out == "5  l5", out


def test_the_last_line_of_a_file_is_readable(call):
    """A window that REACHES the end is not past it."""
    out = call("read_file", {"path": "five.txt", "offset": 4, "limit": 2})
    assert "4  l4" in out and "5  l5" in out, out
    assert "past the end" not in out, out


def test_an_empty_file_is_not_a_wrong_offset(call):
    """Reading a 0-byte file is a legitimate read; it must not be turned
    into an error by the past-the-end repair."""
    out = call("read_file", {"path": "empty.txt"})
    assert not out.lstrip().startswith('{"error"'), out


def test_a_past_the_end_read_grounds_nothing(call, box):
    """The observed-files ledger asks whether a read SHOWED something. A
    wrong offset shows nothing, whatever the tool now says about it."""
    seen: set[str] = set()
    api_client._observe_read_files(
        seen, "read_file", {"path": "five.txt"},
        call("read_file", {"path": "five.txt", "offset": 6}),
        workspace=box)
    assert seen == set()
