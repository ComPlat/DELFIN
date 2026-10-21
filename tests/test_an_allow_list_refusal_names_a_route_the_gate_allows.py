"""An auto-allow refusal names one route, and the gate allows that route.

Input: a bash command refused by the auto-allow tier of
``_DocToolExecutor._run_permission_gate`` in the two modes that reach that
tier with nobody to ask (``default`` and ``acceptEdits`` without a confirm
callback; ``plan`` refuses bash earlier, ``bypassPermissions`` returns
before the tier).

Output under test: the refusal text. For each cause it must name a
concrete route in one line, and the route must be what the same gate, in
the same mode, actually runs. A suggestion the gate refuses again is a
second refusal with extra words, and the measured habit after a refusal
is another spelling of the same command.

Causes, by share of the auto-allow refusals in the tool-trace denial
report (``report_denials`` category ``allowlist``, 67 % of all refusals):

* inline interpreter code -- ``python -c`` (about half) and a
  here-document fed to python -- route: the program as a file through
  ``write_file``, run as ``python3 <file>``;
* a chain or pipe with one refused part -- route: the refused part by
  name, the allowed parts as calls of their own;
* a program not on the list -- the spellings of it that are refused the
  same way (a path, ``env``, ``nohup``, ``time``, ``sh -c``), and the user
  as the route;
* a shell write outside the roots -- route: a workspace-relative path.

Every case drives ``_doc_executor.execute`` for real: the refused
command, then the route, then the refused command again (the route text
must not soften the block it is attached to).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

import delfin.agent.api_client as A

MODES = ("default", "acceptEdits")

OPAQUE_C = 'python3 -c "import os; print(os.getcwd())"'
HEREDOC_DASH = "python3 - <<'EOF'\nimport os\nprint(os.getcwd())\nEOF"
HEREDOC_BARE = "python3 <<EOF\nimport os\nprint(os.getcwd())\nEOF"
PROGRAM = "import os\nprint(os.getcwd())\n"

REFUSED = "not on the auto-allow list"
_BASE_END = "Then STOP and wait."


@pytest.fixture(autouse=True)
def _no_cage_admission(monkeypatch):
    """The refusal text under test is the one given where the caged rule
    does not apply (no bubblewrap cage on the host, or the rule switched
    off). Where it applies, the command runs inside the cage instead:
    tests/test_a_caged_command_runs_without_asking.py."""
    monkeypatch.setattr(A, "_caged_command_may_run_unasked",
                        lambda *a, **k: False)


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    d = tmp_path / "ws"
    d.mkdir()
    (d / "a.txt").write_text("x\n")
    return d


_n = [0]


def _perms(ws: Path, mode: str) -> A.KitToolPermissions:
    _n[0] += 1
    p = A.KitToolPermissions(mode=mode, workspace=str(ws))
    p.task_session_id = f"route-{_n[0]}"
    return p


def _run(perms, tool: str, args: dict) -> dict:
    out = A._doc_executor.execute(tool, args, perms)
    try:
        return json.loads(out)
    except Exception:
        return {"raw": out}


def _error(result: dict) -> str:
    return str(result.get("error") or result.get("raw") or "")


def _bash(perms, cmd: str) -> dict:
    return _run(perms, "bash", {"command": cmd})


def _route_part(err: str) -> str:
    """The text after the fixed refusal sentence: the hints."""
    assert _BASE_END in err, err
    return err.split(_BASE_END, 1)[1]


# ---------------------------------------------------------------------------
# Inline interpreter code: the program as a file
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd", [OPAQUE_C, HEREDOC_DASH, HEREDOC_BARE],
                         ids=["dash-c", "heredoc-dash", "heredoc-bare"])
def test_inline_code_names_a_file_route_the_gate_runs(ws, mode, cmd):
    perms = _perms(ws, mode)
    err = _error(_bash(perms, cmd))
    assert REFUSED in err, err
    hints = _route_part(err)
    path = re.search(r'write_file\(path="([^"]+)"', hints)
    run = re.search(r"run `(python3 [^`]+)`", hints)
    assert path and run, hints

    wrote = _run(perms, "write_file",
                 {"path": path.group(1), "content": PROGRAM})
    assert "error" not in wrote, wrote
    ran = _bash(perms, run.group(1))
    assert ran.get("exit_code") == 0, ran
    assert str(ws) in str(ran.get("stdout", "")), ran

    # The cleanup the route names takes the file back.
    assert "undo_changes(scope=\"last\")" in hints
    undone = _run(perms, "undo_changes", {"scope": "last"})
    assert "error" not in undone, undone
    assert not (ws / path.group(1)).exists()

    # The block itself is unchanged.
    assert REFUSED in _error(_bash(_perms(ws, mode), cmd))


@pytest.mark.parametrize("mode", MODES)
def test_the_file_route_avoids_a_file_that_exists(ws, mode):
    """write_file refuses to overwrite an unread file; naming one would
    be a route the gate refuses."""
    (ws / "_inline_check.py").write_text("print('mine')\n")
    hints = _route_part(_error(_bash(_perms(ws, mode), OPAQUE_C)))
    path = re.search(r'write_file\(path="([^"]+)"', hints).group(1)
    assert not (ws / path).exists(), path


def test_the_file_route_keeps_the_interpreter_the_model_used(ws):
    cmd = 'python3.11 -c "import os; print(os.getcwd())"'
    hints = _route_part(_error(_bash(_perms(ws, "default"), cmd)))
    assert "run `python3.11 _inline_check.py`" in hints, hints
    assert _perms(ws, "default").matches_bash_auto_allow(
        "python3.11 _inline_check.py")


def test_a_cwd_call_gets_no_relative_file_route(ws):
    """A relative path in the route would mean another directory for
    write_file (workspace-relative) than for the command (cwd-relative)."""
    (ws / "sub").mkdir()
    out = _run(_perms(ws, "default"), "bash",
               {"command": OPAQUE_C, "cwd": "sub"})
    err = _error(out)
    assert REFUSED in err, err
    assert "write_file(path=" not in err


# ---------------------------------------------------------------------------
# Chains and pipes: the refused part by name
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd,refused,allowed", [
    ("ls && frobnicate --x", "frobnicate --x", ["ls"]),
    ("cat a.txt; frobnicate; wc -l a.txt", "frobnicate",
     ["cat a.txt", "wc -l a.txt"]),
    ("cat a.txt | xargs echo", "xargs echo", ["cat a.txt"]),
])
def test_a_chain_names_the_part_the_list_refuses(ws, mode, cmd, refused,
                                                 allowed):
    err = _error(_bash(_perms(ws, mode), cmd))
    assert REFUSED in err, err
    hints = _route_part(err)
    assert f"`{refused}` is the part not on the list" in hints, hints
    for part in allowed:
        assert f"`{part}`" in hints, hints
        ran = _bash(_perms(ws, mode), part)
        assert ran.get("exit_code") == 0, (part, ran)
    assert REFUSED in _error(_bash(_perms(ws, mode), refused))


def test_a_chain_with_inline_code_carries_the_file_route(ws):
    cmd = "ls && " + OPAQUE_C
    hints = _route_part(_error(_bash(_perms(ws, "default"), cmd)))
    assert f"`{OPAQUE_C}` is the part not on the list" in hints, hints
    assert 'write_file(path="_inline_check.py"' in hints, hints


def test_an_all_refused_chain_names_no_allowed_part(ws):
    """Nothing to run separately: the line would invite a retry."""
    err = _error(_bash(_perms(ws, "default"), "frobnicate && xargs echo"))
    assert REFUSED in err, err
    assert "is the part not on the list" not in err


# ---------------------------------------------------------------------------
# A program not on the list: its other spellings are refused too
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("cmd", ["frobnicate --x", "xargs echo"])
def test_a_refused_program_lists_spellings_the_gate_also_refuses(
        ws, mode, cmd):
    err = _error(_bash(_perms(ws, mode), cmd))
    assert REFUSED in err, err
    hints = _route_part(err)
    line = re.search(r"Running it as (.+?) is refused the same way", hints)
    assert line, hints
    spellings = re.findall(r"`([^`]+)`", line.group(1))
    assert len(spellings) >= 3, spellings
    for spelling in spellings:
        assert REFUSED in _error(_bash(_perms(ws, mode), spelling)), spelling
    assert "say in your answer which command you need" in hints


def test_with_a_question_dialog_the_route_is_ask_user_question(ws):
    perms = _perms(ws, "default")
    perms.ask_user_callback = lambda q: {"answers": [q["options"][0]["label"]]}
    hints = _route_part(_error(_bash(perms, "frobnicate --x")))
    assert "ask_user_question" in hints, hints
    asked = _run(perms, "ask_user_question", {
        "question": "Run frobnicate --x?",
        "options": [{"label": "yes"}, {"label": "no"}]})
    assert "error" not in asked, asked


def test_without_a_question_dialog_ask_user_question_is_not_named(ws):
    hints = _route_part(_error(_bash(_perms(ws, "default"), "frobnicate")))
    assert "ask_user_question" not in hints


# ---------------------------------------------------------------------------
# A shell write outside the roots: a workspace-relative path
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mode", MODES)
def test_an_outside_write_names_a_relative_path_the_gate_runs(ws, mode):
    err = _error(_bash(_perms(ws, mode), "echo x > ../outside.txt"))
    assert "outside what you may modify" in err, err
    assert "relative paths" in err
    ran = _bash(_perms(ws, mode), "echo x > outside.txt")
    assert ran.get("exit_code") == 0, ran
    assert (ws / "outside.txt").read_text() == "x\n"


# ---------------------------------------------------------------------------
# No route weakens a protection
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    OPAQUE_C, HEREDOC_DASH, "ls && frobnicate --x",
    "cat a.txt | xargs echo", "frobnicate --x",
])
def test_no_route_names_a_switch_past_a_protection(ws, cmd):
    hints = _route_part(_error(_bash(_perms(ws, "default"), cmd)))
    for word in ("bypass", "--no-isolate", "allow_pattern", "dangerously",
                 "remember_permission"):
        assert word not in hints, (word, hints)


@pytest.mark.parametrize("cmd", [OPAQUE_C, HEREDOC_DASH,
                                 "ls && frobnicate --x", "frobnicate --x"])
def test_each_route_is_one_line(ws, cmd):
    hints = _route_part(_error(_bash(_perms(ws, "default"), cmd)))
    routes = [m for m in re.findall(r"ROUTE: [^\n]*", hints)]
    assert routes, hints
    assert "\n" not in hints, hints
