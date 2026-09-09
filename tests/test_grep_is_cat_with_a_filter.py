"""Four other spellings of `cat` walked past the gate that blocks `cat`.

``_gate_bash_read_paths`` exists for one observed incident: three
read_file calls were explicitly denied by the user, and the agent
obtained the same three files seconds later with `cat`. Its docstring
states the principle -- "a refusal that one tool honours and the next
ignores is not a refusal" -- and then the set of commands it inspects was
written narrowly enough to honour it for `cat` alone.

Measured against main on 2026-09-09, workspace elsewhere:

    cat /etc/passwd                 blocked
    grep root /etc/passwd           printed it
    awk '{print}' /etc/passwd       printed it
    sed -n 1,2p /etc/passwd         printed it
    cut -d: -f1 /etc/passwd         printed it

Every one of those four is on the auto-allow list, so an agent whose
read_file was refused could reach the same bytes with the next line it
typed -- and `grep .` is `cat`.

Found while checking whether a change of my own had opened a read path:
`python3 -c "print(open('/etc/passwd').read())"` became auto-allowed when
inline payloads started being read, and the control run against main
showed both that the payload route was new (mine to close, and closed by
collecting the payload's reads alongside its writes) and that four older
routes were already open. The control answered two questions for the
price of one, which is the argument for running it.

Only ABSOLUTE and `~` arguments are inspected, so nothing about ordinary
work changes: `grep -rn foo .`, `sed -n 1,50p src/x.py`, and every
relative path behave exactly as before.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

import delfin.agent.api_client as A


@pytest.fixture
def ws():
    with tempfile.TemporaryDirectory(prefix="rd-") as tmp:
        d = Path(tmp)
        (d / "local.txt").write_text("hello local\n", encoding="utf-8")
        (d / "sub").mkdir()
        (d / "sub" / "a.py").write_text("def f():\n    pass\n", encoding="utf-8")
        yield d


def _run(cmd: str, ws: Path) -> dict:
    perms = A.KitToolPermissions(mode="default", workspace=str(ws))
    return json.loads(A._doc_executor.execute(
        "bash", {"command": cmd, "description": "d"}, perms))


@pytest.mark.parametrize("cmd", [
    "cat /etc/passwd",
    "grep root /etc/passwd",
    "awk '{print}' /etc/passwd",
    "sed -n 1,2p /etc/passwd",
    "cut -d: -f1 /etc/passwd",
    "sort /etc/passwd",
    "rg root /etc/passwd",
    "jq . /etc/passwd",
    "python3 -c \"print(open('/etc/passwd').read())\"",
])
def test_every_way_of_printing_a_file_outside_is_gated(cmd, ws):
    out = _run(cmd, ws)
    assert "error" in out, f"{cmd} -> {str(out)[:120]}"
    assert "/etc/passwd" in out["error"]


@pytest.mark.parametrize("cmd", [
    "grep -rn def .",
    "sed -n 1,5p sub/a.py",
    "cat local.txt",
    "awk '{print}' local.txt",
    "cut -d: -f1 local.txt",
    "sort local.txt",
])
def test_relative_work_is_untouched(cmd, ws):
    """The gate looks at absolute arguments only. If this ever starts
    failing, the fix has become friction: `grep -rn` is the single
    commonest command an agent runs."""
    out = _run(cmd, ws)
    assert out.get("exit_code") == 0, out


def test_an_absolute_path_inside_the_workspace_still_reads(ws):
    out = _run(f"grep hello {ws}/local.txt", ws)
    assert out.get("exit_code") == 0, out
    assert "hello local" in out.get("stdout", "")


@pytest.mark.parametrize("cmd", [
    "wc -l /etc/passwd",
    "file /etc/passwd",
    "stat /etc/passwd",
])
def test_metadata_is_not_content(cmd, ws):
    """A size, a type or a digest is not the file. Gating a question about
    whether a path exists would be friction with nothing behind it."""
    out = _run(cmd, ws)
    assert "error" not in out, out


def test_the_secret_deny_list_still_answers_first(ws):
    """A path on the secret list gets its own message, not the generic
    one -- the more specific reason is the more useful one."""
    out = _run("grep -r . ~/.ssh/id_rsa", ws)
    assert "secret-deny path" in out.get("error", "")


# ---------------------------------------------------------------------------
# list_files: a required argument that meant "everything"
# ---------------------------------------------------------------------------

def test_a_directory_argument_is_honoured_not_ignored(ws):
    """`path` was accepted and silently dropped, so `list_files(path="src")`
    answered with every file in the workspace -- a listing of everything
    presented as the answer to a question about one folder. Callers pass
    it: this repository's own tool-surface test does."""
    perms = A.KitToolPermissions(mode="default", workspace=str(ws))
    out = A._doc_executor.execute("list_files", {"path": "sub"}, perms)
    assert "a.py" in out
    assert "local.txt" not in out


def test_the_default_pattern_is_the_documented_one(ws):
    """`pattern` was marked required while the code has always defaulted it
    to '*'. The schema described a contract the executor did not have."""
    entry = next(t for t in A._DOC_TOOLS_OPENAI
                 if t["function"]["name"] == "list_files")
    assert "required" not in entry["function"]["parameters"]
    perms = A.KitToolPermissions(mode="default", workspace=str(ws))
    out = A._doc_executor.execute("list_files", {}, perms)
    assert "local.txt" in out


def test_a_directory_outside_the_workspace_is_refused(ws):
    perms = A.KitToolPermissions(mode="default", workspace=str(ws))
    out = A._doc_executor.execute("list_files", {"path": "../.."}, perms)
    assert "outside the allowed workspace roots" in out


def test_a_path_that_is_not_a_directory_says_so(ws):
    perms = A.KitToolPermissions(mode="default", workspace=str(ws))
    out = A._doc_executor.execute("list_files", {"path": "nope"}, perms)
    assert "not a directory" in out
