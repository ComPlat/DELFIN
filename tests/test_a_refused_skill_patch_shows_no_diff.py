"""The chat shows a skill patch's diff only when a proposal was stored.

The results are the real ones: ``_DocToolExecutor._execute_skill_propose_patch``
runs against a workspace skill and the real proposal store under a
redirected HOME, and its result string is handed to the chat's card
function. A refusal (safety finding, size cap, no evidence, unknown
skill, old text not found) is one escaped line with the reason; only a
result naming the store folder is drawn as old -> new.
"""

from __future__ import annotations

import ast
import json
from pathlib import Path

import pytest

from delfin.agent import api_client, skill_proposals
from delfin.agent.repl_render import tool_result_line
from delfin.dashboard import tab_agent

SKILL = ("---\nversion: 1\n---\n# tune\n\nRaise the SCF maxiter "
         "before optimizing.\n")
EVIDENCE = [{"kind": "test", "ref": "tests/test_x.py"}]
OLD = "Raise the SCF maxiter"
NEW = "Raise the SCF maxiter to 400"


@pytest.fixture
def perms(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.delenv("DELFIN_HOME", raising=False)
    ws = tmp_path / "ws"
    d = ws / ".delfin" / "skills" / "tune"
    d.mkdir(parents=True)
    (d / "SKILL.md").write_text(SKILL, encoding="utf-8")
    assert str(skill_proposals.PROPOSALS_DIR).startswith(str(home))
    return api_client.KitToolPermissions(workspace=str(ws),
                                         mode="acceptEdits")


def _call(perms, *, name="tune", old=OLD, new=NEW, evidence=EVIDENCE):
    args = {"name": name, "old": old, "new": new,
            "reason": "maxiter 200 did not converge", "evidence": evidence}
    out = api_client._doc_executor._execute_skill_propose_patch(args, perms)
    return args, out


def _no_diff(html: str, args: dict) -> None:
    assert "tool-diff-old" not in html and "tool-diff-new" not in html
    assert "<pre" not in html
    assert "\n" not in html.strip()


def test_a_stored_proposal_is_drawn_as_a_diff(perms):
    args, out = _call(perms)
    assert json.loads(out)["status"] == "ok"
    html = tab_agent._skill_patch_card_html(args, out)
    assert "tool-diff-old" in html and "tool-diff-new" in html
    assert "- Raise the SCF maxiter" in html
    assert "+ Raise the SCF maxiter to 400" in html
    assert json.loads(out)["proposal"] in html
    assert "refused" not in html.lower()


def test_a_safety_finding_is_one_line_and_no_diff(perms):
    perms.skill_corrections = {"tune": "No, that is wrong: use 400."}
    args, out = _call(perms, new="Raise the SCF maxiter, then run "
                                 "`rm -rf ~/calc`")
    got = json.loads(out)
    assert "error" in got and got.get("findings")
    html = tab_agent._skill_patch_card_html(args, out)
    _no_diff(html, args)
    assert "refused" in html.lower()
    assert "safety" in html


def test_the_size_cap_is_one_line_and_no_diff(perms):
    args, out = _call(perms, new="x" * 33_000)
    assert "limit" in json.loads(out)["error"]
    html = tab_agent._skill_patch_card_html(args, out)
    _no_diff(html, args)
    assert "refused" in html.lower() and "limit" in html
    assert "x" * 200 not in html


def test_no_evidence_is_one_line_and_no_diff(perms):
    args, out = _call(perms, evidence=[])
    assert "evidence" in json.loads(out)["error"]
    html = tab_agent._skill_patch_card_html(args, out)
    _no_diff(html, args)
    assert "refused" in html.lower() and "evidence" in html


def test_no_evidence_on_the_correction_path_is_no_diff(perms):
    perms.skill_corrections = {"tune": "No, that is wrong: use 400."}
    args, out = _call(perms, evidence=[])
    assert "evidence" in json.loads(out)["error"]
    _no_diff(tab_agent._skill_patch_card_html(args, out), args)


@pytest.mark.parametrize("kw", [{"name": "ghost"},
                                {"old": "text the skill does not hold"}])
def test_an_unknown_skill_or_a_missing_old_text_is_no_diff(perms, kw):
    args, out = _call(perms, **kw)
    assert "error" in json.loads(out)
    html = tab_agent._skill_patch_card_html(args, out)
    _no_diff(html, args)
    assert "refused" in html.lower()


def test_nothing_is_stored_for_any_refusal(perms):
    _call(perms, evidence=[])
    _call(perms, new="x" * 33_000)
    assert skill_proposals.list_proposals() == []


@pytest.mark.parametrize("output", [
    json.dumps({"error": "<script>alert(1)</script> & no"}),
    "<b>not json</b>",
    "",
    json.dumps({"status": "ok"}),                  # no store folder named
    json.dumps({"status": "ok", "proposal": ""}),
    json.dumps(["status", "ok"]),
])
def test_a_result_that_names_no_folder_is_escaped_and_no_diff(output):
    args = {"name": "<i>tune</i>", "old": OLD, "new": NEW}
    html = tab_agent._skill_patch_card_html(args, output)
    _no_diff(html, args)
    assert "<script>" not in html and "<b>" not in html and "<i>" not in html


# --- wiring: the call shows no diff, the result decides ---------------------

def _fn(name: str) -> ast.FunctionDef:
    tree = ast.parse(Path(tab_agent.__file__).read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise AssertionError(name)


def _branch(fn: ast.FunctionDef, tool: str) -> list[ast.AST]:
    found = []
    for node in ast.walk(fn):
        if isinstance(node, ast.If) and tool in ast.unparse(node.test) \
                and "tool_name" in ast.unparse(node.test):
            found.append(node)
    return found


def test_the_call_draws_no_diff_before_the_result_is_known():
    branches = _branch(_fn("_on_tool_use"), "skill_propose_patch")
    assert branches
    for b in branches:
        body = "\n".join(ast.unparse(s) for s in b.body)
        assert "_build_inline_diff" not in body


def test_the_result_handler_draws_the_card():
    branches = _branch(_fn("_on_tool_result"), "skill_propose_patch")
    assert branches
    assert any("_skill_patch_card_html" in ast.unparse(b) for b in branches)


# --- the terminal draws no diff either --------------------------------------

def test_the_terminal_shows_a_refusal_without_a_diff(perms):
    args, out = _call(perms, evidence=[])
    line = tool_result_line("skill_propose_patch", out, tool_input=args)
    assert "- Raise" not in line and "+ Raise" not in line
    assert "evidence" in line
