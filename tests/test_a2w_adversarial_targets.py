"""Package A2W reviewer (nacht-s12): adversarial tests for
delfin.agent.workspace_selfservice.classify.

Pins:
  * a -r/--requirement value that escapes the workspace in ANY position of
    a multi--r command is refused (the s25 test
    tests/test_a2_adversarial_targets.py::test_pip_multiple_r_flags_last_wins_is_none
    constructs its "outside" file INSIDE the workspace, so it is red against
    correct code -- see REIBUNG_s12.md; these tests pin the real escapes);
  * spellings that are not the pinned word form are refused, never allowed;
  * the git-switch shape refuses when a redirect points outside.
"""

import delfin.agent.workspace_selfservice as ws


def _venv(tmp_path):
    venv = tmp_path / ".venv"
    (venv / "bin").mkdir(parents=True)
    (tmp_path / "requirements.txt").write_text("numpy==1.26.4\n")
    return venv


# --- true multi--r escapes ---------------------------------------------------

def test_double_r_second_outside_is_none(tmp_path):
    # The workspace is tmp_path/ws; the second -r names a file one level
    # ABOVE it -- a real escape, unlike the s25 test whose "outside" file
    # sits inside the workspace.
    ws_root = tmp_path / "ws"
    _venv(ws_root)
    (tmp_path / "evil.txt").write_text("numpy==1.26.4\n")
    cmd = f".venv/bin/pip install -r requirements.txt -r {tmp_path / 'evil.txt'} --require-hashes"
    assert ws.classify(cmd, str(ws_root)) is None


def test_double_r_first_outside_is_none(tmp_path):
    ws_root = tmp_path / "ws"
    _venv(ws_root)
    (tmp_path / "evil.txt").write_text("numpy==1.26.4\n")
    cmd = f".venv/bin/pip install -r {tmp_path / 'evil.txt'} -r requirements.txt --require-hashes"
    assert ws.classify(cmd, str(ws_root)) is None


def test_double_r_both_inside_is_auto(tmp_path):
    # Both -r files inside the workspace: pip processes every -r it is
    # given, so the shape's risk is unchanged; containment is
    # workspace-scoped. Pins the correct reading that the s25 test's
    # assertion (refuse when the second file is merely inside) contradicts.
    _venv(tmp_path)
    (tmp_path / "more.txt").write_text("numpy==1.26.4\n")
    cmd = ".venv/bin/pip install -r requirements.txt -r more.txt --require-hashes"
    assert ws.classify(cmd, str(tmp_path)) == "auto"


# --- spellings that are NOT the pinned form ----------------------------------

def test_attached_short_r_value_is_none(tmp_path):
    # pip's optparse accepts -rreqs as -r reqs; the matcher's word list
    # cannot bind the attached value to the containment check, so the
    # pinned form cannot be verified and must refuse (None), never auto.
    _venv(tmp_path)
    cmd = ".venv/bin/pip install -rrequirements.txt --require-hashes"
    assert ws.classify(cmd, str(tmp_path)) is None


def test_requirement_equals_value_is_none(tmp_path):
    # --requirement=<path> is not a spelling the matcher parses; refusal
    # (None) is the safe direction.
    _venv(tmp_path)
    cmd = ".venv/bin/pip install --requirement=requirements.txt --require-hashes"
    assert ws.classify(cmd, str(tmp_path)) is None


def test_pip3_in_venv_is_none(tmp_path):
    # Only the exact <venv>/bin/pip and <venv>/bin/python -m pip forms are
    # pinned; pip3 is a different binary the shape does not rate.
    _venv(tmp_path)
    (tmp_path / ".venv" / "bin" / "pip3").write_text("#!/bin/sh\n")
    cmd = ".venv/bin/pip3 install -r requirements.txt --require-hashes"
    assert ws.classify(cmd, str(tmp_path)) is None


# --- git switch: redirects and branch-name attacks ---------------------------

def test_git_switch_redirect_outside_is_none(tmp_path):
    cmd = "git switch feature > /etc/a2w_escape_marker"
    assert ws.classify(cmd, str(tmp_path)) is None


def test_git_switch_branch_with_dotdot_middle_is_none(tmp_path):
    cmd = "git switch fe..ature"
    assert ws.classify(cmd, str(tmp_path)) is None


def test_git_switch_branch_unicode_is_none(tmp_path):
    # The safe refname pattern is ASCII-only; a look-alike branch name is
    # refused rather than guessed about.
    cmd = "git switch featur\u00e9"
    assert ws.classify(cmd, str(tmp_path)) is None


def test_git_checkout_b_extra_word_is_none(tmp_path):
    cmd = "git checkout -b feature extra"
    assert ws.classify(cmd, str(tmp_path)) is None


# --- tool --version: no smuggled words ---------------------------------------

def test_version_flag_with_extra_word_is_none(tmp_path):
    assert ws.classify("git --version --exec-path=/etc", str(tmp_path)) is None


def test_version_of_unknown_tool_is_none(tmp_path):
    assert ws.classify("curl --version", str(tmp_path)) is None


def test_version_of_gh_bare_word_is_none(tmp_path):
    # gh has no `version` subcommand (s25 F4); only `git version` is bare.
    assert ws.classify("gh version", str(tmp_path)) is None
