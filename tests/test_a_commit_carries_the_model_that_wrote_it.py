"""The co-author trailer is mounted: a turn sets it, a commit carries it.

The trailer mechanism was complete and unreachable. ``_commit_trailer_value``
builds the line, ``_rewrite_segments_for_trailer`` inserts ``--trailer`` into
every ``git commit`` in a command string -- quote-aware, heredoc-aware,
deduplicating on ``--amend`` -- and ``agent.commit_coauthor`` defaults to on.
But ``set_commit_trailer_model`` was called by NOTHING in the product: a grep
over delfin/ found the definition and no caller, so the context variable was
always empty and ``_commit_trailer_enabled`` was always false. The REWRITE
therefore never fired.

What did fire was the prompt. The system prompt tells the model to end its
commits with the line (tests/test_a_commit_names_the_agent_and_its_model.py),
and main holds 1,858 commits that carry it -- 1,208 naming
kit.glm-5.3-flash, 499 kit.deepseek-v4-flash, 228 kit.glm-5.3, the last of
them 2026-11-19. Three carry the template unexpanded, as
``DELFIN-Agent (<model>)`` and ``DELFIN-Agent (model)``, which is what an
instruction produces and a rewrite cannot: the model copied the shape of the
line instead of filling it in.

So the line was never guaranteed, only asked for. A model that did not
write it left no record, and an orchestrator that told its agents to write
no trailers switched the whole thing off -- which is what happened to the
swarm waves of 2026-11-20: their commits carry nothing.

Its own tests passed throughout, because they call the setter themselves
(tests/test_commit_trailer.py) and then check the rewrite. That is the
shape: a mechanism whose parts are each tested, whose guarantee is wired to
nothing, and whose effect therefore depends on the model remembering.

So this file tests the wiring, not the parts:

* the turn entry point sets the variable from the live client, and sets it
  per turn rather than once, because a successor switch changes the model
  mid-run and the trailer must name the model that produced the diff;
* with it set, a real ``git commit`` command string comes out of the
  isolation layer carrying the trailer.

Input: a client with a model, a git commit command, ``DELFIN_AGENT_LABEL``.
Output: the trailer line, and the rewritten command.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

from delfin.agent import api_client as A

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(autouse=True)
def _no_label(monkeypatch):
    """The label is an environment variable; the suite must not read the
    operator's own."""
    monkeypatch.delenv("DELFIN_AGENT_LABEL", raising=False)
    yield
    A.set_commit_trailer_model("")


# ------------------------------------------------------------- the wiring

def test_the_turn_entry_point_sets_the_model():
    """The mounting itself, read from the source of the method that runs a
    turn.

    A behavioural check would have to drive a whole turn against a live
    backend. What can be checked without one is that the call is in the
    method every turn goes through, and that is exactly what was missing:
    the function existed for months with no caller anywhere.
    """
    src = (ROOT / "delfin" / "agent" / "engine.py").read_text()
    tree = ast.parse(src)
    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "stream_response":
            for call in ast.walk(node):
                if isinstance(call, ast.Call):
                    fn = call.func
                    name = getattr(fn, "id", None) or getattr(fn, "attr", None)
                    if name == "set_commit_trailer_model":
                        found.append(node.name)
    assert found, (
        "nothing in stream_response calls set_commit_trailer_model, so the "
        "whole trailer path is unreachable again and no commit will carry "
        "the model")


def test_the_product_calls_the_setter_somewhere():
    """A tree-wide check, because the defect was the absence of a caller.

    An importer grep over one module would not have found this: the caller
    has to exist SOMEWHERE in delfin/, and that is the invariant.
    """
    callers = []
    for path in (ROOT / "delfin").rglob("*.py"):
        text = path.read_text(errors="replace")
        if "set_commit_trailer_model" not in text:
            continue
        for n, line in enumerate(text.splitlines(), 1):
            if "set_commit_trailer_model" in line and "def " not in line:
                callers.append(f"{path.relative_to(ROOT)}:{n}")
    assert callers, "no caller of set_commit_trailer_model in delfin/"


# ------------------------------------------------------- what the line says

def test_the_trailer_names_the_model():
    A.set_commit_trailer_model("deepseek-v4.1-flash-at-hzdr")
    line = A._commit_trailer_value()
    assert line.startswith("Co-Authored-By: DELFIN-Agent "
                           "(deepseek-v4.1-flash-at-hzdr) <")
    assert "\n" not in line, "no label set, so the trailer is one line"


def test_an_empty_model_turns_it_off():
    A.set_commit_trailer_model("")
    assert A._commit_trailer_enabled() is False


def test_the_label_names_the_run_that_produced_it(monkeypatch):
    """Which agent of which wave is the orchestrator's knowledge, not
    DELFIN's, and the orchestrator is not in this repository -- so the
    environment is the seam."""
    monkeypatch.setenv("DELFIN_AGENT_LABEL", "night3/n3b via blablador")
    A.set_commit_trailer_model("alias-huge")
    line = A._commit_trailer_value()
    assert line.endswith("\nAgent: night3/n3b via blablador")
    assert line.count("\n") == 1


@pytest.mark.parametrize("raw,expect", [
    ("evil\nCo-Authored-By: someone <x@y>", "evil"),
    ("night3/n3b via blablador", "night3/n3b via blablador"),
    ("  padded  ", "padded"),
    ("semi;colon && rm -rf /", "semicolon  rm -rf /"),
    ("$(whoami)", "whoami"),
    ("`id`", "id"),
    ("x" * 300, "x" * 120),
    ("", ""),
])
def test_the_label_cannot_write_a_trailer_of_its_own(monkeypatch, raw, expect):
    """The value lands in a commit message, which other tools parse. A
    newline would end the trailer block and let the rest be read as a
    trailer; shell metacharacters have no business there either."""
    monkeypatch.setenv("DELFIN_AGENT_LABEL", raw)
    assert A._agent_label() == expect


# --------------------------------------------- and the commit carries it

def test_a_git_commit_command_comes_out_carrying_the_trailer(tmp_path, monkeypatch):
    """End to end through the layer every bash command passes, with the
    model set the way a turn sets it."""
    monkeypatch.setenv("DELFIN_AGENT_LABEL", "night3/n3b via blablador")
    A.set_commit_trailer_model("alias-code")
    perms = A.KitToolPermissions(workspace=tmp_path, mode="default",
                                 confirm_callback=None)
    argv = A._bash_isolation_argv("git add f.py && git commit -m 'fix it'",
                                  str(tmp_path), perms)
    whole = " ".join(argv)
    assert "--trailer" in whole
    assert "DELFIN-Agent (alias-code)" in whole
    assert "Agent: night3/n3b via blablador" in whole


def test_a_command_that_is_no_commit_is_untouched(tmp_path):
    A.set_commit_trailer_model("alias-code")
    perms = A.KitToolPermissions(workspace=tmp_path, mode="default",
                                 confirm_callback=None)
    argv = A._bash_isolation_argv("git status --porcelain", str(tmp_path), perms)
    assert "--trailer" not in " ".join(argv)


# ------------------------------------------- and it always comes back out

@pytest.mark.timeout(20)
@pytest.mark.parametrize("cmd,trailer", [
    # The command that hung CI shard 1 on 2026-11-20: no commit in it, a
    # lone `&` in `2>&1`, and once a turn armed the trailer the walk over
    # it never ended (tests/test_a_push_over_https_uses_the_gh_login_
    # outside_the_cage.py::test_no_caged_command_sees_the_token, 300 s).
    ("printf 'protocol=https\\nhost=github.com\\n\\n' | "
     "git credential fill 2>&1; env", False),
    ("git commit -m x 2>&1", True),
    ("git add . && git commit -m 'x' >/dev/null 2>&1", True),
    ("ls &> /tmp/x", False),
    ("sleep 1 &", False),
    ("echo hi \\", False),          # a backslash with nothing to escape
    ("&", False),
    ("\\", False),
])
def test_a_lone_ampersand_does_not_stop_the_walk(cmd, trailer):
    """Input: a command whose unquoted ``&`` is a redirect or a background
    operator, not the ``&&`` the segment branch takes. Output: the command,
    rewritten only where a ``git … commit`` stands.

    The walk consumed characters word by word and broke each word on
    ``&``, while no branch above consumed a lone one -- so the offset
    stood still and the loop ran forever, inside the one layer every bash
    command of every turn passes. The marker bounds the test because the
    regression is a hang, not a wrong answer.
    """
    out = A._rewrite_segments_for_trailer(cmd, "T")
    assert ("--trailer 'T'" in out) is trailer
    assert out.replace(" --trailer 'T'", "") == cmd
