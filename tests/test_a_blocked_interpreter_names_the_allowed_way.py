"""Blocking `python -c` without naming `py_compile` costs the agent its
only way to check its own work.

Found by the new everyday-task probes on their first outing. Asked for a
Tetris, the agent did every single thing right and still ended up flying
blind:

    write_file           ok     tetris.py written
    bash                 BLOCKED  python3 -c "import ast; ast.parse(…)"
    edit_file            ok     )
    multi_edit x3        ok     ) three more passes, none of them verified
    remember_permission  BLOCKED  asked for ^\\s*python3\\s+-c\\b

Both blocks are correct. Inline interpreter code carries its own program
text and so reaches past every write gate — that hole was closed once
already — and granting an allow-pattern for it would re-open exactly that.

But the job the agent wanted is legitimate and HAS a sanctioned spelling:
``python3 -m py_compile <file>`` is on the auto-allow list, and so is
``python3 <file>``. Both act on a file the write gate has already seen,
which is the part ``-c`` skips. The refusal never said so, and its "do
not work around this block with alternative commands" actively pushed
against the one alternative that is not a workaround.

Same shape as the plan-mode refusal fixed the same day: a gate that says
no without saying what yes looks like.

**Update.** The premise above -- "inline interpreter code carries its own
program text and so reaches past every write gate" -- was half right. It
carries its own program text, and that is the reason it CAN be read: the
program is a literal in the command line, unlike `xargs`, `make` or a
base64 payload. Both checks the blanket ban stood in for run on it now
(delfin/agent/inline_payload.py), so the ast.parse in the incident above
runs, and the ban narrowed to payloads the analysis cannot account for.

The hint still matters and is still tested here, because that narrower
set is where the agent lands next: `python3 -c "import os; …"` is opaque
and `python3 -m py_compile` is still the sanctioned spelling.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def ws():
    with tempfile.TemporaryDirectory(prefix="interp-") as tmp:
        d = Path(tmp)
        (d / "t.py").write_text("def f():\n    return 1\n")
        yield d


def run(cmd: str, ws: Path) -> dict:
    perms = KitToolPermissions(mode="default", workspace=str(ws))
    return json.loads(_doc_executor.execute("bash", {"command": cmd}, perms))


OPAQUE = 'python3 -c "import os; print(os.getcwd())"'


@pytest.mark.parametrize("cmd", [
    OPAQUE,
    'python -c "exec(open(\'p\').read())"',
    'python3.11 -c "import subprocess; subprocess.run([\'ls\'])"',
])
def test_inline_interpreter_code_is_still_refused(cmd, ws):
    """The hint must not soften the block it is attached to."""
    out = run(cmd, ws)
    assert "not on the auto-allow list" in out.get("error", "")


def test_the_refusal_names_the_spelling_that_works(ws):
    err = run(OPAQUE, ws).get("error", "")
    assert "py_compile" in err


def test_it_says_why_that_is_not_a_workaround(ws):
    """Without this the hint reads as a loophole, and the next reader
    deletes it."""
    err = run(OPAQUE, ws).get("error", "")
    assert "write gate already saw" in err


def test_the_job_the_agent_wanted_now_runs(ws):
    """The incident's own command, verbatim. It was the whole reason this
    file exists: four edit passes, none of them verified, because the one
    way to check the work was refused."""
    out = run(
        'python3 -c "import ast; ast.parse(open(\'t.py\').read());'
        ' print(\'syntax OK\')"', ws)
    assert out.get("exit_code") == 0, out
    assert "syntax OK" in str(out.get("stdout", ""))


def test_the_named_alternative_actually_runs(ws):
    """A hint pointing at a command that is ALSO blocked would be worse
    than no hint."""
    out = run("python3 -m py_compile t.py", ws)
    assert out.get("exit_code") == 0, out


def test_running_the_file_itself_is_allowed_too(ws):
    out = run("python3 t.py", ws)
    assert out.get("exit_code") == 0, out


def test_a_plain_blocked_command_gets_no_interpreter_hint(ws):
    """The hint is for one situation. Attached to everything it becomes
    noise, and noise is what people learn to skip."""
    err = run("curl https://example.com", ws).get("error", "")
    assert "not on the auto-allow list" in err
    assert "py_compile" not in err


def test_the_cd_hint_still_fires(ws):
    """It shares the branch; a careless edit would shadow it.

    The command has to be one the gate actually refuses -- `cd /tmp && ls`
    is on the auto-allow list and simply runs, so a hint test built on it
    would be green whatever the branch did.
    """
    err = run("cd /tmp && curl https://example.com", ws).get("error", "")
    assert "not on the auto-allow list" in err
    assert "cwd" in err


def test_asking_for_a_c_allow_pattern_is_still_refused(ws):
    """The escape the agent reached for. Granting it would undo the fix
    that took `-c` off the list in the first place."""
    perms = KitToolPermissions(mode="default", workspace=str(ws))
    out = _doc_executor.execute(
        "remember_permission",
        {"kind": "allow_pattern", "value": r"^\s*python3\s+-c\b",
         "rationale": "syntax checks", "scope": "repo"}, perms)
    assert "error" in str(out)
