"""Red test for the commits column: _default_commits must read the
member's own workspace (the record's ``workspace`` field, what the
presence record announces), not the calling process's cwd -- the service
runs with cwd= its runtime dir, so the cwd route answers 0 commits for
every member whose repo is elsewhere. Run from the worktree root."""
from pathlib import Path
import subprocess

from delfin.agent.swarm import observe as _observe


def test_default_commits_read_the_members_workspace(tmp_path):
    ws = tmp_path / "ws"
    ws.mkdir()
    env = {"GIT_AUTHOR_NAME": "rr4", "GIT_AUTHOR_EMAIL": "r@e.x",
           "GIT_COMMITTER_NAME": "rr4", "GIT_COMMITTER_EMAIL": "r@e.x",
           "PATH": "/usr/bin:/bin"}
    def git(*args, cwd):
        subprocess.run(["git", *args], cwd=str(cwd), env=env,
                       capture_output=True, check=True)
    git("init", "-q", "-b", "main", ".", cwd=ws)
    (ws / "a.txt").write_text("base\n")
    git("add", "a.txt", cwd=ws)
    git("commit", "-qm", "base", cwd=ws)
    git("checkout", "-qb", "swarm/rr4-pr", cwd=ws)
    (ws / "fix.txt").write_text("fix\n")
    git("add", "fix.txt", cwd=ws)
    git("commit", "-qm", "rr4: fix", cwd=ws)
    git("checkout", "-q", "main", cwd=ws)

    records = [{"name": "rr4-pr", "session_id": "rr4-pr-sess1",
                "workspace": str(ws)}]
    # The overview's default path: the record's workspace is where git runs.
    rows = _observe.swarm_overview({"rr4-pr": {"package": "rr4", "role": "pr"}},
                                   records=records,
                                   steps_fn=lambda *a, **k: ([], 0))
    assert rows[0]["commits_on_branch"] == 1
    assert rows[0]["last_commit_subject"] == "rr4: fix"
