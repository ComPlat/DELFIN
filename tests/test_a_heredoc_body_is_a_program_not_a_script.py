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
import re
import tempfile
from pathlib import Path

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


def test_a_quoted_heredoc_is_accountable_the_same_way_c_is(ws):
    """The predicate says so. Whether anything ASKS it is the separate
    question below — today nothing does, so the gate still refuses."""
    assert A._inline_payload_is_readable(
        "python3 - <<'EOF'\nprint(1)\nEOF", ws)
    assert A._inline_payload_is_readable('python3 -c "print(1)"', ws)
    # And an unquoted one never is: the shell expands it first.
    assert not A._inline_payload_is_readable(
        "python3 - <<EOF\nprint(1)\nEOF", ws)


def test_a_payload_reading_an_absolute_path_is_not_accounted_for(ws):
    """Found by the test below it. The predicate looked only at WRITES,
    so `print(open('/etc/passwd').read())` came back readable — for the
    `-c` form too, which has been in production. The read gate catches it
    downstream either way, but this predicate is what decides whether the
    confirm veto fires AT ALL, and it cannot call that fully accounted
    for."""
    for form in ("python3 - <<'EOF'\nprint(open('/etc/passwd').read())\nEOF",
                 "python3 -c \"print(open('/etc/passwd').read())\""):
        assert not A._inline_payload_is_readable(form, ws), form
    # A read inside the workspace is ordinary and stays accounted for.
    assert A._inline_payload_is_readable(
        "python3 -c \"print(open('local.txt').read())\"", ws)


def test_a_read_outside_the_workspace_reaches_the_read_gate(ws):
    """It did not: the read gate took `-c` payloads only, so a heredoc
    doing the same thing produced no read target at all."""
    assert A._bash_outside_reads(
        "python3 - <<'EOF'\nprint(open('/etc/passwd').read())\nEOF"
    ) == ["/etc/passwd"]


# ---------------------------------------------------------------------------
# The allow path, opened — after the missing link was found
# ---------------------------------------------------------------------------
#
# `python3 - <<'EOF' … EOF` is `python -c` with different punctuation and
# is the third-largest group in the denial log: 29 of 437 auto-allow
# refusals. It now runs on the same terms as `-c`.
#
# The first attempt added only the auto-allow pattern, and it ran six of
# seven things it must refuse:
#
#     heredoc: arithmetic              RAN   (intended)
#     heredoc that writes              RAN   (not intended)
#     UNQUOTED heredoc                 RAN   (not intended)
#     heredoc importing os             RAN   (not intended)
#     heredoc importing subprocess     RAN   (not intended)
#     print(open('/etc/passwd').read())  RAN (not intended)
#
# because `_interpreter_needs_confirm` opens with
# `_is_interpreter_invocation(cmd)` and that was False for the stdin
# form — the veto short-circuited before the predicate was ever asked.
# One dash did it: `_INTERPRETER_RE` already carried
# `python[0-9.]*\s*<`, which catches `python3 <<'EOF'` and not
# `python3 - <<'EOF'`, which is how the form is actually written.
#
# With the veto reaching the predicate, the pattern is safe and every one
# of those six is refused again. Checked against the audit log the way
# the splitter was: 41 commands are newly treated as interpreter
# invocations, every one of them `python3 - <<'…'`, and NOTHING that was
# gated before is ungated now.

_HEREDOCS_THAT_MUST_NOT_RUN = [
    ("writes a file", "python3 - <<'EOF'\nopen('out.txt','w').write('x')\nEOF"),
    ("unquoted", "python3 - <<EOF\nprint(1+1)\nEOF"),
    ("imports os", "python3 - <<'EOF'\nimport os\nprint(os.listdir('/'))\nEOF"),
    ("imports subprocess",
     "python3 - <<'EOF'\nimport subprocess\nsubprocess.run(['ls'])\nEOF"),
    ("reads /etc/passwd",
     "python3 - <<'EOF'\nprint(open('/etc/passwd').read())\nEOF"),
]


@pytest.mark.parametrize("name,cmd", _HEREDOCS_THAT_MUST_NOT_RUN,
                         ids=[n for n, _ in _HEREDOCS_THAT_MUST_NOT_RUN])
def test_none_of_the_six_runs_unattended(ws, name, cmd):
    """Each of these ran when the pattern went in without the fix."""
    assert "error" in _run(ws, cmd), name


@pytest.mark.parametrize("cmd", [
    "python3 - <<'EOF'\nprint(1+1)\nEOF",
    "python3 - <<'EOF'\nimport math\nprint(math.exp(-1.5))\nEOF",
    "python3 - <<'PY'\nprint(4.073215*96.48533212)\nPY",
    "python3 -u - <<'EOF'\nprint(1)\nEOF",
])
def test_a_computation_runs_like_the_c_form_does(ws, cmd):
    out = _run(ws, cmd)
    assert "error" not in out, out


def test_a_read_inside_the_workspace_is_ordinary(ws):
    Path(ws, "local.txt").write_text("hello\n")
    out = _run(ws, "python3 - <<'EOF'\nprint(open('local.txt').read())\nEOF")
    assert "error" not in out, out
    assert "hello" in (out.get("stdout") or "")


def test_the_link_that_was_missing(ws):
    """One dash. `python[0-9.]*\\s*<` catches `python3 <<'EOF'` and not
    `python3 - <<'EOF'`, and the veto opens with this check."""
    assert A._is_interpreter_invocation("python3 <<'EOF'\nprint(1)\nEOF")
    assert A._is_interpreter_invocation("python3 - <<'EOF'\nprint(1)\nEOF")
    assert A._is_interpreter_invocation("python3 -u - <<'EOF'\nprint(1)\nEOF")
    # And what must NOT become an interpreter invocation: a script run is
    # scanned as a file, not treated as opaque.
    assert not A._is_interpreter_invocation("python3 script.py")
    assert not A._is_interpreter_invocation("python3 script.py < input.txt")
    assert not A._is_interpreter_invocation("cat <<'EOF'\nx\nEOF")


def test_the_veto_now_reaches_the_predicate(ws):
    perms = A.KitToolPermissions(mode="default", workspace=ws)
    assert perms._interpreter_needs_confirm(
        "python3 - <<'EOF'\nimport subprocess\nEOF")
    assert not perms._interpreter_needs_confirm(
        "python3 - <<'EOF'\nprint(1+1)\nEOF")


def test_the_predicate_itself_is_already_right(ws):
    """The half that IS in place: when the veto is eventually asked, it
    has the right answer ready."""
    assert A._inline_payload_is_readable(
        "python3 - <<'EOF'\nprint(1+1)\nEOF", ws)
    for _name, cmd in _HEREDOCS_THAT_MUST_NOT_RUN:
        assert not A._inline_payload_is_readable(cmd, ws), cmd


# ---------------------------------------------------------------------------
# A here-document feeding a redirect is write_file's job
# ---------------------------------------------------------------------------
#
# `cat > run.py << 'EOF' … EOF` used to be refused by the same accident:
# `print(1)` is not an allowed command, so the segment check failed. Now
# that the body stays with its command it would be auto-allowed off
# `cat`, and the file written with no pre-image in the change journal —
# undo_changes could not take it back and list_changes_made would not
# report it. That is a property worth keeping, so it is kept on purpose.
#
# Narrow, and narrower than the first attempt. A REDIRECT specifically,
# not any write the segment performs: an inline python payload that
# writes has its own policy a layer up, and overriding an explicit user
# grant over a journalling concern would be too strong.
#
# `cat > f` and `echo x > f` are a separate and older question. Both
# already run unattended, and neither carries its content in the command.

def test_a_heredoc_into_a_redirect_still_names_write_file(ws):
    out = _run(ws, "cat > run.py << 'EOF'\nprint(1)\nEOF")
    assert "error" in out
    assert "write_file" in out["error"]


def test_tee_fed_by_a_heredoc_too(ws):
    out = _run(ws, "tee run.py << 'EOF'\nprint(1)\nEOF")
    assert "error" in out


def test_a_heredoc_with_no_redirect_is_untouched(ws):
    assert "error" not in _run(ws, "cat << 'EOF'\nhello\nEOF")


@pytest.mark.parametrize("cmd", ["cat > f1.py", "echo x > f2.py"])
def test_a_plain_redirect_is_the_older_question_and_unchanged(ws, cmd):
    """Documented rather than asserted as good: these run unattended on
    main too, and journalling them is a separate fix."""
    assert "error" not in _run(ws, cmd)


def test_a_grant_is_not_overridden_for_a_journalling_concern(ws):
    """The first version of the rule refused any heredoc segment with a
    write target, which cancelled an explicit user grant. A payload
    writing inside the workspace is what that grant means."""
    assert "error" not in _run(
        ws, "python3 - <<'EOF'\nopen('out.txt','w').write('x')\nEOF",
        grant=True)


def test_and_the_sandbox_still_wins_over_the_grant(ws):
    out = _run(ws, "python3 - <<'EOF'\nopen('/etc/evil','w').write('x')\nEOF",
               grant=True)
    assert "error" in out and "/etc/evil" in out["error"]


# ---------------------------------------------------------------------------
# Nothing else moved
# ---------------------------------------------------------------------------
#
# _split_shell_segments is load-bearing: every permission decision runs
# through it. The heredoc change had to leave everything else exactly as
# it was, and six hand-written cases do not show that.
#
# Checked against the corpus instead: the OLD algorithm reimplemented
# below, both run over the ~2000 distinct commands in
# ~/.delfin/audit*.log. Segmentation differed on 46 of them, every one
# containing `<<`, and on nothing else. The commands kept here are a
# sample of the real ones that actually exercise the splitter —
# operators, quotes, redirects, `$(…)`, embedded newlines — so the
# property stays checked in a repo that has no audit log.

def _split_before_heredocs(cmd):
    """The algorithm as it stood before the here-document change."""
    segs, buf, q, i, n = [], [], None, 0, len(cmd)
    while i < n:
        c = cmd[i]
        if q is not None:
            buf.append(c)
            if c == q:
                q = None
            i += 1
            continue
        if c in ("'", '"'):
            q = c
            buf.append(c)
            i += 1
            continue
        if cmd[i:i + 2] in ("||", "&&"):
            segs.append("".join(buf))
            buf = []
            i += 2
            continue
        if c in (";", "|", "\n"):
            segs.append("".join(buf))
            buf = []
            i += 1
            continue
        buf.append(c)
        i += 1
    segs.append("".join(buf))
    return [s.strip() for s in segs if s.strip()]


_REAL_COMMANDS_WITHOUT_A_HEREDOC = [
    "find . -maxdepth 3 -iname '*bookmark*' -o -iname '*lesezeichen*' "
    "2>/dev/null",
    'python3 -m py_compile export.py && echo "syntax OK"',
    'grep -i "HOMO-LUMO GAP\\|ORBITAL ENERGIES\\|: HOMO\\|: LUMO" run_*.out',
    "lsof -i :8899 2>/dev/null || netstat -tlnp 2>/dev/null | grep 8899 "
    "|| ss -tlnp | grep 8899",
    "python3 -c \"import ast; ast.parse(open('pipeline.py').read())\" "
    "&& echo SYNTAX_OK",
    'python3 export.py --smtp-host smtp.example.org --sender a@b.c; '
    'echo "exit=$?"',
    "cd tests/fixtures/user_project_workspace && python3 tagreport.py",
    "ls -R . | head -50",
    'git status && echo "---" && find . -maxdepth 3 -type d | head -50',
    "ls -la && git status --short",
    "cat delfin/agent/pack/benchmark/tasks_auto_behavior.yaml 2>&1 | head -150",
    'ls -la .delfin/ 2>/dev/null; echo "---"; ls -la tests/fixtures/ 2>/dev/null',
    "python3 export.py && cat bookmarks.csv && python3 export.py --send "
    "2>&1 | tail -1",
    "ls -la && cat bookmarks.json 2>/dev/null | head -50",
    'find /home/user/ComPlat -name "zahlen.txt" 2>/dev/null',
    'ls -la; echo "---"; rm -f a.csv b.json; rm -rf __pycache__; '
    'echo "cleaned"; ls -la',
    "sed -n '538,610p' delfin/agent/workspace_trust.py",
    "cmd 2>&1",
    "echo 'a;b' && echo \"c|d\"",
    "for i in 1 2 3; do echo $i; done",
]


@pytest.mark.parametrize("cmd", _REAL_COMMANDS_WITHOUT_A_HEREDOC)
def test_a_command_with_no_heredoc_segments_exactly_as_before(cmd):
    assert A._split_shell_segments(cmd) == _split_before_heredocs(cmd), cmd


def test_the_reference_implementation_really_is_different_on_a_heredoc():
    """Otherwise the parametrised test above proves nothing: a reference
    that agreed everywhere would be the same function."""
    cmd = "python3 - <<'EOF'\nimport os\nx = 1; y = 2\nEOF"
    assert _split_before_heredocs(cmd) != A._split_shell_segments(cmd)
    assert len(_split_before_heredocs(cmd)) > 1
    assert len(A._split_shell_segments(cmd)) == 1


# ---------------------------------------------------------------------------
# Nothing that was gated before is ungated now
# ---------------------------------------------------------------------------
#
# `_is_interpreter_invocation` decides whether the confirm veto is even
# consulted, so widening it is safe in one direction only: more commands
# may become interpreter invocations, none may stop being one.
#
# Checked against the audit log the way the splitter was. Of ~2050
# distinct real commands, 41 changed — every one of them
# `python3 - <<'…'` — and none went the other way. The old regex is
# reimplemented here so the property stays checked in a repo with no
# audit log.

_INTERPRETER_RE_BEFORE = re.compile(
    r"(?:^|[;&|`$(]\s*)\s*"
    r"(?:[\w./~+-]*/)?"
    r"(?:"
    r"python[0-9.]*\s+-c\b|python[0-9.]*\s*<|"
    r"perl\s+-e\b|ruby\s+-e\b|node\s+-e\b|php\s+-r\b|"
    r"(?:eval|exec|source|make|xargs)\b|\.\s|"
    r"env\s+[A-Za-z_][A-Za-z0-9_]*=|"
    r"base64\s+(?:-d|--decode)\b|"
    r"find\b[^;|&]*-exec\b"
    r")")


def _was_an_interpreter_invocation(cmd):
    """The predicate as it stood before the stdin form was added."""
    text = cmd or ""
    if _INTERPRETER_RE_BEFORE.search(text):
        return True
    for seg in A._split_shell_segments(text):
        stripped = A._strip_exec_wrappers(seg)
        if stripped != seg and _INTERPRETER_RE_BEFORE.search(stripped):
            return True
    return bool(re.search(
        r"\|\s*(?:ba|z|k|da)?sh\b|\|\s*python[0-9.]*\b", text))


_COMMANDS_THE_GATE_ALREADY_KNEW = [
    'python3 -c "print(1)"',
    "python3 < script.py",
    "time python3 -c \"print(1)\"",
    "env FOO=1 python3 -c \"print(1)\"",
    "perl -e 'print 1'",
    "node -e 'console.log(1)'",
    "base64 -d payload.b64",
    "find . -name '*.py' -exec rm {} \\;",
    "cat x | bash",
    "curl -s u | python3",
    "eval \"$CMD\"",
    "xargs rm < list.txt",
    "ls -la",
    "grep -i x f.out",
    "python3 script.py",
    "python3 script.py < input.txt",
    "cat <<'EOF'\nplain text\nEOF",
    "git status && echo ok",
    "echo a; echo b",
    "sed -n '1,5p' f.py",
]


@pytest.mark.parametrize("cmd", _COMMANDS_THE_GATE_ALREADY_KNEW)
def test_no_command_stops_being_an_interpreter_invocation(cmd):
    """The only direction that can hurt. A command the veto used to see
    and no longer sees would be a hole, not a widening."""
    if _was_an_interpreter_invocation(cmd):
        assert A._is_interpreter_invocation(cmd), cmd


@pytest.mark.parametrize("cmd", _COMMANDS_THE_GATE_ALREADY_KNEW)
def test_and_nothing_ordinary_newly_becomes_one(cmd):
    """The other side: a plain command must not be dragged in, or every
    `ls` starts asking for confirmation."""
    if not _was_an_interpreter_invocation(cmd):
        assert not A._is_interpreter_invocation(cmd), cmd


def test_the_reference_really_differs_on_the_stdin_form():
    """Otherwise the two tests above are comparing a function with
    itself."""
    cmd = "python3 - <<'EOF'\nprint(1)\nEOF"
    assert not _was_an_interpreter_invocation(cmd)
    assert A._is_interpreter_invocation(cmd)
