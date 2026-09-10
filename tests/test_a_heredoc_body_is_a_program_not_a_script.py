"""A here-document body was read as a chain of shell commands.

`python3 - <<'EOF' … EOF` is `python -c` with different punctuation: the
program is a literal in the command line. It is the third-largest group
in the denial log — 29 of 437 auto-allow refusals — and nothing looked
at it.

_split_shell_segments split on newlines, so the body became one
pseudo-segment per line of the program: `import os`, `x = 1`,
`os.system("rm -rf /tmp/x")`. The auto-allow check requires every
segment to be allowed, and those are not commands, so it refused —
correctly, by accident, and only while the lines fail to READ as
commands. With a user grant of `^\\s*python3\\b`, which is exactly what
the refusal tells the model to ask for, a body of `ls` / `cat` was
allowed through without one line of it being read as Python.

The body now stays with its command as one opaque unit, and the write
gate is shown the payload: a body of `open('/etc/evil','w')` produced no
write target at all before, and with the same grant it ran.

Unquoted `<<EOF` is expanded by the shell first, so the text analysed is
not necessarily the text that runs. The write gate scans it anyway —
that path REFUSES on what it finds, so scanning an expandable body can
only block more. The argument does not hold where a decision GRANTS, so
`_inline_payload_is_readable` still takes `-c` only.
"""

from __future__ import annotations

import json
import tempfile

import pytest

import delfin.agent.api_client as A
from delfin.agent.inline_payload import analyze_payload, extract_stdin_payloads

_n = [0]


@pytest.fixture
def ws():
    with tempfile.TemporaryDirectory() as tmp:
        yield tmp


def _perms(ws, grant=False):
    _n[0] += 1
    p = A.KitToolPermissions(mode="default", workspace=ws)
    p.task_session_id = f"hd-{_n[0]}"
    if grant:
        p.bash_auto_allow_patterns = tuple(p.bash_auto_allow_patterns) + (
            r"^\s*python3\b",)
    return p


def _run(ws, cmd, grant=False):
    out = A._doc_executor.execute(
        "bash", {"command": cmd, "description": "d"}, _perms(ws, grant))
    try:
        return json.loads(out)
    except Exception:
        return {"raw": out}


# ---------------------------------------------------------------------------
# The body stays with its command
# ---------------------------------------------------------------------------

def test_a_heredoc_is_one_segment():
    segs = A._split_shell_segments(
        "python3 - <<'EOF'\nimport os\nx = 1; y = 2\nEOF")
    assert len(segs) == 1, segs
    assert "import os" in segs[0]


def test_a_command_after_the_terminator_is_its_own_segment():
    """The newline after the terminator separates commands; swallowing it
    hid `echo done` inside the heredoc."""
    segs = A._split_shell_segments("cat <<'X'\na;b\nX\necho done")
    assert len(segs) == 2, segs
    assert segs[1] == "echo done"


def test_a_here_string_is_not_a_here_document():
    """`<<<` needs no terminator. Stepping over only two of its three
    characters made the remaining `<<` look like a heredoc named 'hi'."""
    segs = A._split_shell_segments('cat <<<"hi"; echo x')
    assert len(segs) == 2, segs


def test_a_tab_stripped_heredoc_ends_at_its_indented_terminator():
    segs = A._split_shell_segments("python3 - <<-'T'\n\tprint(1)\n\tT\necho after")
    assert len(segs) == 2, segs
    assert segs[1] == "echo after"


def test_two_heredocs_are_two_segments():
    segs = A._split_shell_segments("cat <<'A'\nx\nA\ncat <<'B'\ny\nB")
    assert len(segs) == 2, segs


def test_an_unterminated_body_runs_to_the_end():
    """What the shell would read. Stopping early would hand the caller a
    segment the shell never sees."""
    segs = A._split_shell_segments("python3 - <<'EOF'\nprint(1)")
    assert len(segs) == 1, segs


@pytest.mark.parametrize("cmd,want", [
    ("echo a; echo b", 2),
    ("ls | wc -l", 2),
    ('grep "a|b" f', 1),
    ("cat a.txt", 1),
    ("ls && echo ok", 2),
    ("cmd 2>&1", 1),
])
def test_ordinary_splitting_is_unchanged(cmd, want):
    assert len(A._split_shell_segments(cmd)) == want, cmd


# ---------------------------------------------------------------------------
# The write gate sees the payload
# ---------------------------------------------------------------------------

def test_a_write_outside_the_workspace_is_blocked_even_with_a_grant(ws):
    cmd = "python3 - <<'EOF'\nopen('/etc/evil','w').write('x')\nEOF"
    assert A._bash_write_targets(cmd), "the write is invisible to the gate"
    out = _run(ws, cmd, grant=True)
    assert "error" in out, out
    assert "/etc/evil" in out["error"]


def test_an_unquoted_body_is_scanned_too(ws):
    """It may be expanded before the interpreter sees it, so what is read
    here is a lower bound on what it does. This path refuses, so a lower
    bound is the safe direction."""
    cmd = "python3 - <<EOF\nopen('/etc/evil2','w').write('x')\nEOF"
    out = _run(ws, cmd, grant=True)
    assert "error" in out and "/etc/evil2" in out["error"]


def test_a_write_inside_the_workspace_is_what_the_grant_allows(ws):
    cmd = "python3 - <<'EOF'\nopen('out.txt','w').write('x')\nEOF"
    assert A._bash_write_targets(cmd) == ["out.txt"]
    assert "error" not in _run(ws, cmd, grant=True)


def test_a_computation_writes_nothing_and_runs_under_a_grant(ws):
    cmd = "python3 - <<'EOF'\nprint(1+1)\nEOF"
    assert not A._bash_write_targets(cmd)
    out = _run(ws, cmd, grant=True)
    assert "error" not in out, out
    assert "2" in (out.get("stdout") or "")


def test_the_deny_list_still_reads_the_body(ws):
    """It scans the whole command string, so a body calling out to `rm
    -rf` was already caught — and must stay caught now the body is one
    segment."""
    out = _run(ws, "python3 - <<'EOF'\nimport os\nos.system('rm -rf /tmp/zz')\nEOF",
               grant=True)
    assert "error" in out and "deny-pattern" in out["error"]


# ---------------------------------------------------------------------------
# Extraction: quoted is literal, unquoted is not
# ---------------------------------------------------------------------------

def test_only_a_quoted_delimiter_is_taken_by_default():
    assert extract_stdin_payloads(
        "python3 - <<'EOF'\nprint(1)\nEOF") == ["print(1)"]
    assert extract_stdin_payloads("python3 - <<EOF\nprint(1)\nEOF") == []


def test_an_unquoted_body_is_available_to_a_refusing_caller():
    assert extract_stdin_payloads(
        "python3 - <<EOF\nprint(1)\nEOF", include_expanded=True) == ["print(1)"]


def test_a_heredoc_into_something_other_than_python_is_not_a_payload():
    assert extract_stdin_payloads("cat <<'X'\nnot python\nX") == []


def test_the_tabs_a_dash_heredoc_strips_are_stripped():
    """`<<-` strips leading tabs before the interpreter sees the text.
    Leaving them in hands ast.parse an IndentationError and the payload
    comes back opaque for a reason that is not in the program."""
    body = extract_stdin_payloads("python3 - <<-'T'\n\tprint(1)\n\tT")
    assert body == ["print(1)"]
    assert not analyze_payload(body[0], "/tmp").opaque


def test_each_of_two_heredocs_is_extracted():
    assert extract_stdin_payloads(
        "python3 - <<'A'\nprint(1)\nA\npython3 - <<'B'\nprint(2)\nB"
    ) == ["print(1)", "print(2)"]


def test_a_command_with_no_heredoc_yields_nothing():
    assert extract_stdin_payloads("python3 script.py") == []
    assert extract_stdin_payloads("") == []


def test_the_granting_path_still_takes_c_payloads_only(ws):
    """Being wrong in the allow direction runs something unread. A quoted
    heredoc is literal and could be read the same way, but that is a
    separate decision and this is not it."""
    assert not A._inline_payload_is_readable(
        "python3 - <<'EOF'\nprint(1)\nEOF", ws)
    assert A._inline_payload_is_readable('python3 -c "print(1)"', ws)
