""""Allow for good" covers the part of the command that was asked about.

Field report 20261020-140059 (dashboard, KIT, acceptEdits). Two commands
in a row were confirmed with "Allow + Permanent":

    python3 - <<'EOF' ... EOF            (asked: an interpreter hides its target)
    node --check /tmp/t.js && echo OK

    python3 - <<'EOF' ... EOF
    node --check .t.js && echo OK; rm -f .t.js     (asked: rm -f .t.js)

Both offered ^\\s*python3\\b -- the first word of the whole command. For
the second that rule does not cover the segment that was asked about, so
the click persisted nothing new, the agent was told twice that python3 is
now allowed, and the same shape of command would be asked about again.

The suggestion is now built from the segments the gate refuses, and a
rule that is already in force is not announced a second time.
"""

from __future__ import annotations

import re

from delfin.agent import api_client as A
from delfin.agent.kit_settings import suggest_pattern_for_command

FIRST = ("python3 - <<'EOF'\nsrc=open('3d-tetris/index.html').read()\n"
         "open('/tmp/t.js','w').write(src)\nEOF\n"
         "node --check /tmp/t.js && echo \"JS SYNTAX OK\"")
SECOND = ("python3 - <<'EOF'\nsrc=open('3d-tetris/index.html').read()\n"
          "open('.t.js','w').write(src)\nEOF\n"
          "node --check .t.js && echo \"JS SYNTAX OK\"; rm -f .t.js")


def _perms(tmp_path, *extra):
    p = A.KitToolPermissions(workspace=str(tmp_path), mode="acceptEdits")
    p.bash_auto_allow_patterns = p.bash_auto_allow_patterns + tuple(extra)
    return p


def test_the_suggestion_lets_the_asked_command_through(tmp_path):
    for cmd in (FIRST, SECOND):
        assert not _perms(tmp_path).matches_bash_auto_allow(cmd)
        pattern = suggest_pattern_for_command(cmd)
        re.compile(pattern)
        assert _perms(tmp_path, pattern).matches_bash_auto_allow(cmd), (cmd, pattern)


def test_after_the_first_grant_the_second_suggestion_is_about_rm(tmp_path):
    granted = _perms(tmp_path, suggest_pattern_for_command(FIRST))
    pattern = suggest_pattern_for_command(SECOND, perms=granted)
    # The rm it was asked about, exactly -- never a bare rm.
    assert re.search(pattern, "rm -f .t.js")
    assert not re.search(pattern, "rm -rf /home/user")
    assert not re.search(pattern, "rm -f other.js")
    granted.bash_auto_allow_patterns += (pattern,)
    assert granted.matches_bash_auto_allow(SECOND)


def test_a_single_command_keeps_its_old_suggestion():
    assert suggest_pattern_for_command("pytest -xvs tests/foo.py") == r"^\s*pytest\b"
    assert suggest_pattern_for_command("git status") == r"^\s*git\s+status\b"
    assert suggest_pattern_for_command("rm build/tmp.txt") == r"^\s*rm\ build/tmp\.txt\s*$"


def test_a_push_inside_a_compound_stays_exact(tmp_path):
    cmd = "git add -A && git push -u origin feature"
    pattern = suggest_pattern_for_command(cmd)
    p = _perms(tmp_path, pattern)
    assert p.matches_bash_auto_allow(cmd)
    assert not p.matches_bash_auto_allow("git push --force origin main")


def test_a_rule_already_in_force_is_not_announced_again(tmp_path, monkeypatch):
    from delfin.agent import kit_settings
    from delfin.agent.engine import AgentEngine

    told = []
    import types

    eng = AgentEngine.__new__(AgentEngine)
    eng.client = types.SimpleNamespace(_permissions=_perms(tmp_path))
    eng.repo_dir = str(tmp_path)
    eng._tell_the_running_turn = told.append
    monkeypatch.setattr(kit_settings, "persist_pattern", lambda *a, **k: None)
    pattern = r"^\s*python3\b"
    assert AgentEngine.persist_kit_pattern(eng, pattern)[0]
    assert AgentEngine.persist_kit_pattern(eng, pattern)[0]
    assert len(told) == 1
    assert eng.kit_permissions.bash_auto_allow_patterns.count(pattern) == 1
