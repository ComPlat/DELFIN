"""grep_file and read_file disagreed about what line 35 is.

grep emitted ``i + 1`` (1-based); read_file emitted ``i + offset``
(0-based). Reproduced on this repository before the fix:

    grep:                 delfin/agent/code_nav.py:35: _MAX_MATCHES = 50
    read_file(offset=35):  35  _MAX_FILES_GREP = 5000

Both print "35" and they mean different lines.

The standard loop is grep -> read -> edit. The agent greps a symbol, gets
`:35`, reads at 35, lands one line off, and copies the wrong block into
`edit_file(old_string=...)`. It then gets "old_string not found" -- or
worse, the fuzzy fallback finds a different unique match and edits THAT.
These are the two most-used tools in the recorded traces: read_file 396
calls, grep_file 76.

read_file is now 1-based, which is the half that had to move: grep's
numbering is what the rest of the world uses, and it is what the model
sees in every error message and stack trace it reads.
"""

from __future__ import annotations

import pathlib

import pytest

from delfin.agent import api_client as A


@pytest.fixture
def ex(tmp_path):
    f = tmp_path / "sample.py"
    f.write_text("\n".join(f"line{i}" for i in range(1, 21)) + "\n",
                 encoding="utf-8")
    perms = A.KitToolPermissions(workspace=tmp_path)
    return A._DocToolExecutor.__new__(A._DocToolExecutor), perms


# ---------------------------------------------------------------------------
# They agree
# ---------------------------------------------------------------------------

def test_grep_and_read_name_the_same_line(ex):
    executor, perms = ex
    hit = executor._execute_grep_file(
        {"pattern": "line7", "path": "sample.py"}, perms).strip()
    assert ":7:" in hit, hit
    got = executor._execute_read_file(
        {"path": "sample.py", "offset": 7, "limit": 1}, perms)
    assert got.splitlines()[0].startswith("7  line7"), got


def test_the_first_line_is_line_one(ex):
    executor, perms = ex
    got = executor._execute_read_file(
        {"path": "sample.py", "offset": 1, "limit": 1}, perms)
    assert got.splitlines()[0].startswith("1  line1"), got


def test_offset_zero_is_taken_as_the_beginning(ex):
    """A model that read the schema as 0-based gets the first line, not
    an error -- the numbers still line up afterwards."""
    executor, perms = ex
    got = executor._execute_read_file(
        {"path": "sample.py", "offset": 0, "limit": 2}, perms)
    assert got.splitlines()[0].startswith("1  line1"), got


def test_a_window_is_numbered_from_its_own_start(ex):
    executor, perms = ex
    got = executor._execute_read_file(
        {"path": "sample.py", "offset": 5, "limit": 3}, perms)
    nums = [ln.split("  ")[0] for ln in got.splitlines() if ln[:1].isdigit()]
    assert nums[:3] == ["5", "6", "7"], got


# ---------------------------------------------------------------------------
# The footer says what it did
# ---------------------------------------------------------------------------

def test_the_truncation_footer_names_the_real_range(ex):
    executor, perms = ex
    got = executor._execute_read_file(
        {"path": "sample.py", "offset": 3, "limit": 4}, perms)
    assert "showing 3-6" in got, got


def test_a_window_reaching_the_end_has_no_footer(ex):
    executor, perms = ex
    got = executor._execute_read_file(
        {"path": "sample.py", "offset": 18, "limit": 50}, perms)
    assert "lines total" not in got, got
    assert got.splitlines()[-1].startswith("20  line20"), got


# ---------------------------------------------------------------------------
# The schema says which convention it is
# ---------------------------------------------------------------------------

def test_the_schema_declares_the_convention():
    for tool in A._DOC_TOOLS_OPENAI:
        fn = tool.get("function", {})
        if fn.get("name") != "read_file":
            continue
        desc = str((fn.get("parameters", {}).get("properties", {})
                    .get("offset", {})).get("description", ""))
        assert "1-based" in desc, desc
        return
    pytest.fail("read_file is not in the catalogue")
