"""Characterization for the shell-word / security-event move.

Pins today's behaviour BEFORE the move: every name moved to a leaf module
(delfin/agent/shell_words.py, delfin/agent/security_events.py) is the SAME
object when imported from the old place (delfin.agent.api_client) and the
new place, and importing the leaf alone never pulls the 39k-line api_client
into sys.modules.
"""
import os
import subprocess
import sys
from pathlib import Path

# The helpers other modules import from api_client that this PR moves to
# delfin.agent.shell_words.py (a leaf importing only stdlib).
SHELL_NAMES = [
    "_split_shell_segments",
    "_shell_segments",
    "_with_shell_bodies",
    "_command_words",
]
CONSTANTS = [
    "_consume_heredoc",
    "_strip_heredoc_bodies",
    "_lines_from",
    "_HEREDOC_START_RE",
    "_SHELL_BODY_RE",
    "_GH_SEGMENT_SPLIT_RE",
    "_COMMAND_WRAPPERS",
    "_GH_OPTS_WITH_VALUE",
]


def test_record_security_event_is_the_same_object_from_old_and_new():
    import delfin.agent.api_client as A
    import delfin.agent.security_events as S
    assert A._record_security_event is S._record_security_event


def test_shell_word_helpers_are_the_same_object_from_old_and_new():
    import delfin.agent.api_client as A
    import delfin.agent.shell_words as W
    for name in SHELL_NAMES + CONSTANTS:
        assert getattr(A, name) is getattr(W, name), name


def test_importing_the_leaf_does_not_pull_api_client():
    # A fresh interpreter: this test process may already hold api_client.
    code = (
        "import sys;"
        "import delfin.agent.shell_words as W;"
        "assert W._split_shell_segments('a; b') == ['a', 'b'];"
        "assert 'delfin.agent.api_client' not in sys.modules, "
        "'api_client imported by the leaf';"
        "import delfin.agent.security_events as S;"
        "S.clear();"
        "S._record_security_event('x-test', 'bash', 'd', blocked=True);"
        "assert any(e.kind == 'x-test' for e in S.recent());"
        "assert 'delfin.agent.api_client' not in sys.modules, "
        "'api_client imported via security_events'"
    )
    root = Path(__file__).resolve().parents[1]
    env = dict(os.environ)
    env["PYTHONPATH"] = str(root) + (
        os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    env["PYTHONNOUSERSITE"] = "1"
    p = subprocess.run([sys.executable, "-c", code], capture_output=True,
                       text=True, cwd=str(root), env=env)
    assert p.returncode == 0, p.stderr
