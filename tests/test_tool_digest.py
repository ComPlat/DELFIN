"""Controls for delfin.agent.tool_digest — the wanted-poster that replaces an
elided tool result.

Every case here encodes the REAL shape of a tool result as api_client emits
it (numbered lines for read_file, "path:line: text" for grep_file, a JSON
payload for bash), not a guess. digest() must be deterministic, never leak
secrets, and stay under ~300 chars.
"""

from __future__ import annotations

from delfin.agent.tool_digest import digest


READ_BODY = "\n".join([
    "1  from __future__ import annotations",
    "2  ",
    "3  ",
    "4  def first(arg):",
    "5      return arg",
    "6  ",
    "7  ",
    "8  class Thing:",
    "9      pass",
    "10  ",
    "... (120 lines total, showing 1-10)",
])

GREP_BODY = "\n".join([
    "delfin/api.py:12: def alpha():",
    "delfin/api.py:34: class Beta:",
    "tests/test_a.py:7: alpha",
    "tests/test_b.py:99: alpha",
])


def test_read_file_names_path_range_and_defs():
    d = digest("read_file", {"path": "delfin/mod.py", "offset": 1, "limit": 10},
               READ_BODY)
    assert "delfin/mod.py" in d
    assert "lines 1-10" in d          # the slice actually shown
    assert "120 lines" in d           # total
    assert "def first" in d
    assert "class Thing" in d
    assert len(d) <= 300


def test_read_file_range_from_tail_marker_only():
    d = digest("read_file", {"path": "x.py"}, READ_BODY)
    # Without offset/limit args the tail marker is the evidence.
    assert "lines 1-10" in d and "120 lines" in d


def test_grep_names_pattern_count_and_first_hits():
    d = digest("grep_file", {"pattern": "alpha", "path": "."}, GREP_BODY)
    assert "alpha" in d
    assert "4 hits" in d
    assert "delfin/api.py:12" in d
    assert "delfin/api.py:34" in d
    assert "tests/test_a.py:7" in d
    assert "tests/test_b.py:99" not in d   # only the first 3


def test_grep_no_matches_says_so():
    d = digest("grep_file", {"pattern": "zzz"}, "No matches found.")
    assert "zzz" in d
    assert "0 hits" in d or "no matches" in d.lower()


def test_bash_json_gives_command_exit_and_last_line():
    content = (
        '{"exit_code": 1, "elapsed_s": 2.0, "stdout": "line one\\nline two\\n'
        '3 failed, 12 passed in 4.5s", "stderr": "", "command": "pytest -q", '
        '"description": "", "cwd": "."}'
    )
    d = digest("bash", {"command": "pytest -q"}, content)
    assert "pytest -q" in d
    assert "exit 1" in d
    assert "3 failed, 12 passed" in d


def test_bash_test_run_lists_failed_names():
    content = (
        '{"exit_code": 1, "elapsed_s": 9.0, "stdout": "FAILED '
        'tests/test_x.py::test_a - assert 0\\nFAILED '
        'tests/test_x.py::test_b - assert 1\\nFAILED '
        'tests/test_x.py::test_c - assert 2\\nFAILED '
        'tests/test_x.py::test_d - assert 3\\n2 passed, 4 failed in 3s", '
        '"stderr": "", "command": "gate tests/test_x.py", '
        '"description": "", "cwd": "."}'
    )
    d = digest("run_tests", {"target": "tests/test_x.py"}, content)
    assert "test_a" in d and "test_b" in d and "test_c" in d
    assert "test_d" not in d          # only the first FAILED names
    assert "2 passed, 4 failed" in d


def test_bash_plain_text_falls_back_to_last_line():
    d = digest("bash", {"command": "make build"}, "step 1\nstep 2\nall done")
    assert "make build" in d
    assert "all done" in d


def test_fallback_first_line_and_size():
    d = digest("search_docs", {"query": "freq"}, "first result line\nsecond\n")
    assert "search_docs" in d
    assert "first result line" in d


def test_digest_is_deterministic():
    a = digest("read_file", {"path": "x.py"}, READ_BODY)
    b = digest("read_file", {"path": "x.py"}, READ_BODY)
    assert a == b


def test_secrets_are_scrubbed():
    d = digest(
        "bash", {"command": "cat .env"},
        '{"exit_code": 0, "stdout": "api_key = AKIA1234567890EXAMPLE\\n", '
        '"stderr": "", "command": "cat .env", "cwd": "."}')
    assert "AKIA1234567890EXAMPLE" not in d
