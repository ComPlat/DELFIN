"""An auto-allowed interpreter turned every write gate into an optional one.

In `default` mode the shipped auto-allow table contained
``^\\s*python(?:3(?:\\.\\d+)?)?\\s+-c\\s+``, and the guard that keeps an
interpreter off that table fired only under a locked scope. Filesystem
containment is off in those modes too, ``_bash_write_targets`` recognises
cp/mv/tee/mkdir/rm/sed -i/dd and not python, and the payload scanner reads
referenced script FILES — of which ``-c`` has none.

So: write_file on a protected path is refused, the user clicks Deny, and
the next call is

    bash python3 -c "open('delfin/agent/api_client.py','a').write(...)"

with no dialog, no isolation and no record. One route past the
read-before-write contract, the protected-path globs, the
self-modification guard, the workspace boundary and the calc confirm at
the same time.

A second, quieter hole in the same guard: the interpreter may be named
through a path. ``.venv/bin/python -c "…"`` read as an ordinary word to
the guard's regex and was auto-allowed even under a locked scope, where
the whole promise is that an opaque command reaches the user.

Cost of the fix, measured over 2458 recorded interactive bash calls
(``~/.delfin/audit-*.log``, bypassPermissions excluded because it never
consults the auto-allow table): 107 calls that used to run unattended now
raise one confirm — 6.1% of them, 9.0% of the unique commands, and 89 of
the 90 unique commands are ``python -c``. The routine surface (ls, git
status, grep, pytest, ``python script.py``, ``python -m pytest``) is
untouched, which is what the second half of this file pins down.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from delfin.agent import api_client as A


def _perms(tmp_path, **kw):
    return A.KitToolPermissions(workspace=tmp_path, **kw)


def _locked(tmp_path, **kw):
    return A.KitToolPermissions(
        workspace=tmp_path, agent_role="office_agent", **kw)


# ---------------------------------------------------------------------------
# The hole
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    'python3 -c "open(\'delfin/agent/api_client.py\',\'a\').write(1)"',
    'python -c "import shutil; shutil.rmtree(\'/home/u/data\')"',
    "python3 -c 'print(1)'",
    "make build",
    "cat list.txt | xargs cat",
    "find . -name x -exec cat {} +",
    'eval "$(echo ls)"',
    "echo Y2F0 | base64 -d | bash",
])
@pytest.mark.parametrize("mode", ["default", "acceptEdits"])
def test_an_interpreter_is_not_auto_allowed_in_any_mode(tmp_path, cmd, mode):
    perms = _perms(tmp_path, mode=mode)
    assert perms.matches_bash_auto_allow(cmd) is False, cmd


@pytest.mark.parametrize("cmd", [
    '.venv/bin/python -c "open(\'x\',\'w\')"',
    '/usr/bin/python3 -c "print(1)"',
    'app/.venv-proj/bin/python3.11 -c "print(1)"',
])
def test_a_path_named_interpreter_is_still_an_interpreter(tmp_path, cmd):
    """The venv fallback auto-allowed `.venv/bin/python -c` even under a
    lock, because the guard's regex only knew the bare command name."""
    assert A._is_interpreter_invocation(cmd) is True, cmd
    assert _perms(tmp_path).matches_bash_auto_allow(cmd) is False, cmd
    assert _locked(tmp_path).matches_bash_auto_allow(cmd) is False, cmd


def test_the_denied_write_has_no_shell_route_left(tmp_path):
    """The incident end to end: a write the user refused, retried through
    `python -c`. The command no longer runs unattended, and the refusal
    ledger then refuses it outright rather than asking a second time."""
    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    target = tmp_path / "keep.txt"
    target.write_text("original", encoding="utf-8")
    asked: list[str] = []

    def _deny(name, args, preview):
        asked.append(name)
        return False

    perms = _perms(tmp_path, mode="default", confirm_callback=_deny)
    first = ex._run_permission_gate(
        "write_file", {"path": "keep.txt", "content": "x"}, perms)
    assert first is not None and "denied" in first
    assert asked == ["write_file"]

    cmd = f'python3 -c "open(\'{target}\',\'a\').write(1)"'
    # 1. it no longer runs unattended...
    assert perms.matches_bash_auto_allow(cmd) is False
    # 2. ...and the path gate _execute_bash runs next refuses it outright,
    #    because the file itself is what the user declined.
    blocked = ex._gate_bash_write_targets(cmd, {}, perms)
    assert blocked is not None and "keep.txt" in blocked
    assert asked == ["write_file"]                 # nothing re-asked
    assert target.read_text(encoding="utf-8") == "original"


# ---------------------------------------------------------------------------
# ...without making the gate one nobody leaves on
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "ls -la",
    "pwd",
    "git status",
    "git diff --stat",
    "grep -rn foo .",
    "rg pattern src/",
    "cat notes.md",
    "pytest -q tests/",
    "python3 --version",
    "python -m pytest tests/ -v",
    "python -m molkit 'Fe2(SO4)3'",
    "python3 scripts/report.py --out x",
    "ruff check .",
    "mypy delfin/",
    "cd /tmp && ls",
    'find . -name "*.py"',
    "python3.10 --version || which python3.10 || echo nf",
])
def test_the_routine_case_still_runs_unattended(tmp_path, cmd):
    """A gate that asks about everything gets switched off. Running a
    script and running a module are equally arbitrary code and stay
    auto-allowed: what changed is only the commands whose target cannot
    be read at all."""
    assert _perms(tmp_path).matches_bash_auto_allow(cmd) is True, cmd


def test_a_rule_the_user_wrote_still_applies(tmp_path):
    """The block message tells the model to ask the user for an
    allow_pattern (which always confirms before it is stored). If the
    stored rule then did nothing, the advice would send the agent into a
    loop it cannot leave."""
    base = _perms(tmp_path)
    with_rule = _perms(
        tmp_path,
        bash_auto_allow_patterns=tuple(base.bash_auto_allow_patterns)
        + (r"^\s*python3?\s+-c\s+",),
    )
    assert with_rule.matches_bash_auto_allow('python3 -c "print(1)"') is True
    # ...and it is the USER's rule that does it, not a shipped default.
    assert base.matches_bash_auto_allow('python3 -c "print(1)"') is False


def test_no_rule_reopens_a_locked_scope(tmp_path):
    """Where the folder IS the promise, the answer was decided in advance."""
    locked = _locked(
        tmp_path,
        bash_auto_allow_patterns=tuple(_perms(tmp_path).bash_auto_allow_patterns)
        + (r"^\s*python3?\s+-c\s+",),
    )
    assert locked.matches_bash_auto_allow('python3 -c "print(1)"') is False


def test_a_default_pattern_cannot_pass_for_a_user_rule():
    """The built-in table is exactly what the rule overrides."""
    ws = Path(tempfile.mkdtemp()).resolve()
    perms = A.KitToolPermissions(workspace=ws)
    for pat in A._DEFAULT_BASH_AUTO_ALLOW:
        assert pat in A._BUILTIN_BASH_AUTO_ALLOW
    assert perms._custom_allow_matches("python3 -c 'x'") is False
