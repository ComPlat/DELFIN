"""Real-SIGKILL crash test for Target 1 (crash-resume) — reviewer s2.

QS lens: the kill must be REAL (os.kill SIGKILL on a live child), never
a mock. The child process writes a mid-turn checkpoint with the REAL
writer, then blocks forever on a pipe read (deterministic, no timing
dependence); the parent SIGKILLs it and verifies what a restart reads.

On main 04bf9703b the checkpoint part works (green) — the note content
contract is covered by test_f1_adversarial_crash_resume.py (red). This
battery adds the real-kill proof of the durable artefacts: the
checkpoint survives the kill, and the step journal's open step survives
it too.
"""

from __future__ import annotations

import os
import signal
import subprocess
import sys
import textwrap
import time
from pathlib import Path

import pytest

from delfin.agent import session_store
from delfin.agent import step_journal
from delfin.agent.state_paths import RedirectedUserState

_CHILD = textwrap.dedent("""
    import sys, time
    from pathlib import Path
    sys.path.insert(0, {repo!r})
    from delfin.agent import session_store, step_journal
    from delfin.agent.state_paths import RedirectedUserState
    with RedirectedUserState({home!r}):
        session_store.save_turn_checkpoint("sess-kill", {{
            "user_message": "the task that was running",
            "partial_response": "half-written",
            "tool_calls": 3,
            "ts": time.time(),
        }})
        step_journal.start_step("sess-kill", turn=4, step=2,
                                tool="git_commit", arguments={{"m": "fix"}})
        ready = Path({ready!r})
        ready.write_text("ok", encoding="utf-8")
        # Block forever on stdin: the parent kills us mid-"turn".
        sys.stdin.buffer.read()
""")
REPO = str(Path.cwd())
TIMEOUT_S = 60


def _spawn_child(home: Path, ready: Path) -> subprocess.Popen:
    code = _CHILD.format(repo=REPO, home=str(home), ready=str(ready))
    return subprocess.Popen([sys.executable, "-c", code],
                            stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)


class TestRealSigkill:
    def test_checkpoint_and_open_step_survive_a_real_sigkill(self, tmp_path):
        ready = tmp_path / "child_ready"
        with RedirectedUserState(tmp_path):
            proc = _spawn_child(tmp_path, ready)
            try:
                # Wait for the child to signal its artefacts are written.
                deadline = time.monotonic() + TIMEOUT_S
                while not ready.exists():
                    if proc.poll() is not None:
                        out, err = proc.communicate(timeout=5)
                        pytest.fail(
                            "child died before writing: "
                            f"{err.decode(errors='replace')[-800:]}")
                    assert time.monotonic() < deadline, (
                        "child never signalled readiness")
                    time.sleep(0.05)
                # THE REAL KILL: SIGKILL mid-turn, no cleanup, no mock.
                proc.send_signal(signal.SIGKILL)
                proc.wait(timeout=10)
                assert proc.returncode == -signal.SIGKILL
                # Everything the restart reads is still on disk:
                assert session_store.load_turn_checkpoint("sess-kill") \
                    is not None, "checkpoint must survive SIGKILL"
                open_out, notes = step_journal.open_steps("sess-kill")
                assert [(e["turn"], e["step"]) for e in open_out] == [(4, 2)]
                assert open_out[0]["side_effect"] is True
                assert open_out[0]["tool"] == "git_commit"
                # And the recovery path the restart takes reads them:
                note = session_store.consume_crash_recovery_note("sess-kill")
                assert note is not None
                assert "sess-kill" not in note  # note is about the task
                assert "git_commit" in note
                assert "step 2" in note
                assert "CHECK BEFORE REPEATING" in note
            finally:
                if proc.poll() is None:
                    proc.kill()
                    proc.wait(timeout=10)
