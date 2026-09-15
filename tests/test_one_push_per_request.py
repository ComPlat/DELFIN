"""One request, one push.

Report 20260915-085107: asked to commit and push, the agent did -- and in
the same turn, after an auto-verify round, chained
``git add && git commit && git push origin main`` onto a test fix nobody had
asked for. The session ran all_free, which asks nothing, and the prompt's
"never push unprompted" was a sentence. A push is now granted by the user's
message, spent by the push that goes through, and refused without a grant
where nobody can be asked.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as A
from delfin.agent import job_monitor as jm
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

_PUSH_OUTPUT = ("To github.com:ComPlat/DELFIN.git\n"
                "   032b8093..d507039a  main -> main\n")


@pytest.fixture(autouse=True)
def _watch_index(tmp_path, monkeypatch):
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")


def _perms(tmp_path, mode="bypassPermissions"):
    return KitToolPermissions(workspace=tmp_path, mode=mode)


def _gate(perms, cmd):
    return _DocToolExecutor()._run_permission_gate(
        "bash", {"command": cmd}, perms)


def _pushed(perms, exit_code=0, output=_PUSH_OUTPUT):
    return A._after_push(
        {"command": "git push origin main"},
        json.dumps({"exit_code": exit_code, "stdout": "", "stderr": output}),
        perms)


@pytest.mark.parametrize("text", [
    "commiten und pushen mit Co-Authored-By: DELFIN-Agent GLM-5.3",
    "please push it",
    "puschen bitte",
])
def test_a_request_to_push_is_recognised(tmp_path, text):
    perms = _perms(tmp_path)
    A._grant_push_from(perms, text, new_request=True)
    assert perms.push_grants.get("push") == 1


@pytest.mark.parametrize("text", [
    "richte die CI wieder",
    "bitte nicht pushen",
    "don't push yet",
    "commit it",
])
def test_other_requests_grant_no_push(tmp_path, text):
    perms = _perms(tmp_path)
    A._grant_push_from(perms, text, new_request=True)
    assert not perms.push_grants.get("push")


def test_without_a_request_the_unattended_profile_refuses_a_push(tmp_path):
    perms = _perms(tmp_path)
    A._grant_push_from(perms, "fix the collector test", new_request=True)
    err = _gate(perms, "git add t.py && git commit -m fix && git push origin main")
    assert err and "has not asked for a push" in err


def test_the_request_lets_one_push_through_and_not_a_second(tmp_path):
    perms = _perms(tmp_path)
    A._grant_push_from(perms, "commiten und pushen", new_request=True)

    assert _gate(perms, "git push origin main") is None
    assert "ci:ComPlat/DELFIN@d507039a" in _pushed(perms)

    err = _gate(perms, "git push origin main")
    assert err and "has not asked for a push" in err


def test_a_push_that_failed_keeps_the_grant(tmp_path):
    perms = _perms(tmp_path)
    A._grant_push_from(perms, "push it", new_request=True)
    assert _pushed(perms, exit_code=1) == ""
    assert _gate(perms, "git push origin main") is None


def test_a_steered_request_adds_a_grant_and_a_steered_remark_keeps_it(tmp_path):
    perms = _perms(tmp_path)
    A._grant_push_from(perms, "fix the viewer", new_request=True)
    assert not perms.push_grants.get("push")
    A._grant_push_from(perms, "commiten und pushen", new_request=False)
    assert perms.push_grants.get("push") == 1
    A._grant_push_from(perms, "use a shorter commit message", new_request=False)
    assert perms.push_grants.get("push") == 1


def test_a_new_request_takes_an_unspent_grant_back(tmp_path):
    perms = _perms(tmp_path)
    A._grant_push_from(perms, "push it", new_request=True)
    A._grant_push_from(perms, "look at the CI log", new_request=True)
    assert not perms.push_grants.get("push")


def test_only_a_push_is_a_push():
    assert A._is_git_push("git push origin main")
    assert A._is_git_push("cd repo && git -C sub push --force-with-lease")
    assert not A._is_git_push("git status && echo push")
    assert not A._is_git_push("grep -n 'git push' notes.txt")


def test_a_human_is_asked_instead_of_refused(tmp_path):
    perms = _perms(tmp_path, mode="default")
    asked: list[str] = []
    perms.confirm_callback = lambda name, args, preview: asked.append(preview) or True
    assert _gate(perms, "git push origin main") is None
    assert asked and "not asked for a push" in asked[0]


def test_an_answer_picked_in_the_dialog_grants_the_push(tmp_path):
    """Report 20260915-132613: the user picked "Ja, committen und pushen"
    in ask_user_question, and the gate refused the push twice."""
    perms = _perms(tmp_path)
    A._grant_push_from(perms, "mach die Tests grün", new_request=True)
    A._grant_push_from_answer(perms, json.dumps(
        {"answers": ["Ja, committen und pushen"], "multiSelect": False}))
    assert _gate(perms, "git push origin main") is None


def test_an_answer_that_declines_grants_nothing(tmp_path):
    perms = _perms(tmp_path)
    A._grant_push_from(perms, "mach die Tests grün", new_request=True)
    A._grant_push_from_answer(perms, json.dumps(
        {"answers": ["Nein, noch nicht pushen"]}))
    assert "has not asked for a push" in _gate(perms, "git push origin main")


def test_a_question_to_the_user_ends_the_turn():
    """After "Soll ich pushen?" auto-continue sent the agent back in, and it
    tried the push it had just asked about."""
    import inspect
    assert A._ends_with_a_question(
        "Soll ich `git push origin main` jetzt ausführen?")
    assert A._ends_with_a_question("Welche Variante willst du?**")
    assert not A._ends_with_a_question("Push ist durch. CI läuft.")
    src = inspect.getsource(A.OpenAIClient.stream_message)
    i = src.index("_did_tools_since_cont and _auto_cont_count < _AUTO_CONT_CAP")
    assert "not _ends_with_a_question(" in src[i:i + 400]


def test_the_push_arms_a_watch_on_its_ci(tmp_path):
    perms = _perms(tmp_path)
    A._grant_push_from(perms, "push it", new_request=True)
    _pushed(perms)
    jobs = jm.load_watched(tmp_path / ".delfin" / "agent_watched_jobs.json")["jobs"]
    entry = jobs["ci:ComPlat/DELFIN@d507039a"]
    assert entry["kind"] == "ci"
    assert entry["branch"] == "main"
