"""The push gate sees every spelling of git, and grants only on a request.

Review 2026-09-16 of the one-push-per-request gate: `/usr/bin/git push`
passed both the contributor block and the grant check, because the
pattern wanted `git` at a word boundary; and any message with "push" in it
granted one -- "der Push ist gestern fehlgeschlagen" included -- which in
the bypass profile is the only gate there is.
"""
import pytest

from delfin.agent import api_client as A


@pytest.mark.parametrize("cmd", [
    "/usr/bin/git push origin main",
    "cd repo && /opt/homebrew/bin/git -C sub push --force-with-lease",
    "git -c core.sshCommand=ssh push origin HEAD:main",
    "sh -c 'git push'",
    'bash -lc "cd repo && git push origin feature"',
])
def test_every_spelling_of_git_push_is_a_push(cmd):
    assert A._is_git_push(cmd)


@pytest.mark.parametrize("cmd", [
    "git status && echo push",
    "grep -n 'git push' notes.txt",
    "gitpush",
    "digit push",
    "echo 'git push' > notes.txt",
])
def test_what_is_not_a_push_stays_not_one(cmd):
    assert not A._is_git_push(cmd)


@pytest.mark.parametrize("text", [
    "commiten und pushen",
    "please push it",
    "puschen bitte",
    "kannst du das auf main pushen?",
    "commit and push",
    "pushe den Branch",
])
def test_a_request_grants(text):
    assert A._asks_for_push(text) is True


@pytest.mark.parametrize("text", [
    "der Push ist gestern fehlgeschlagen",
    "warum hat der letzte Push nicht funktioniert?",
    "the push failed yesterday, look at the log",
    "ich habe schon gepusht",
    "I pushed it an hour ago, CI is red",
    "bitte nicht pushen",
    "don't push yet",
    "commit it",
    "",
])
def test_a_report_or_a_refusal_grants_nothing(text):
    assert A._asks_for_push(text) is False


def test_a_report_and_a_request_in_one_message_is_a_request():
    assert A._asks_for_push("der Push ist fehlgeschlagen, fix es und push nochmal") is True



def test_a_push_inside_a_shell_string_names_its_targets(tmp_path):
    assert A._push_targets("sh -c 'git push origin HEAD:main'", tmp_path) == {"main"}
