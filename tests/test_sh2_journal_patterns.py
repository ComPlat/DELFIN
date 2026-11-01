"""Pin: shell_changes journaling catches the week's shell mistakes.

The week's real mistakes (wave 26 brief): a `> 1e-6:` inside a python
one-liner read as a redirection; inline heredocs; a `git checkout HEAD
-- file` that threw away 291 uncommitted lines. shell_changes.before/
after must see each of these as a journalled change with an undoable
pre-image, so undo_changes can take it back.

These are pins of EXISTING behaviour (all four were green on
f45828f79 before any change), not a fix: the audit found the patterns
covered. Documented boundaries this does not contradict: a change that
leaves size, mtime, inode and ctime all equal is not seen
(shell_changes.py docstring, "Not seen" list) -- a coarse-timestamp
filesystem showed that once in the audit probe, and it stays a stated
limitation rather than a defect (the module's "Rejected alternatives"
paragraph rules out the re-read that would fix it).
"""
from __future__ import annotations

import subprocess
from pathlib import Path

from delfin.agent import shell_changes as sc


def _run(cmd: str, cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run(["sh", "-c", cmd], cwd=cwd, capture_output=True,
                          text=True)


def test_git_checkout_discard_is_journalled_with_pre_image(tmp_path):
    # Repo first with the base content, THEN the uncommitted WIP lines,
    # so `git checkout --` discards real uncommitted work.
    proc = _run("git init -q . && printf 'def a():\\n    return 1\\n' "
                "> analysis.py && git add analysis.py && "
                "git -c user.name=t -c user.email=t@e.invalid "
                "commit -q -m x", tmp_path)
    assert proc.returncode == 0, proc.stderr
    file = tmp_path / "analysis.py"
    file.write_text("def a():\n    return 1\n" + "# WIP\n" * 10,
                    encoding="utf-8")
    snap = sc.before([tmp_path])
    proc = _run("git checkout -- analysis.py", tmp_path)
    assert proc.returncode == 0, proc.stderr
    changes = sc.after(snap)
    kinds = {Path(c.path).name: (c.kind, c.pre is not None, c.reason)
             for c in changes}
    assert kinds["analysis.py"][0] == "modified", kinds
    # The pre-image is held, so the discarded WIP lines are undoable.
    assert kinds["analysis.py"][1] is True, kinds
    assert kinds["analysis.py"][2] == "", kinds


def test_stray_redirection_file_is_journalled_as_created(tmp_path):
    snap = sc.before([tmp_path])
    # The week's mistake: `> 1e-6:` inside a one-liner -- the shell
    # creates a file literally named '1e-6:'.
    proc = _run("printf 'x\\n' > 1e-6:", tmp_path)
    assert proc.returncode == 0, proc.stderr
    changes = sc.after(snap)
    created = [c for c in changes if Path(c.path).name == "1e-6:"]
    assert created, changes
    assert created[0].kind == "created"
    assert created[0].reason == ""


def test_heredoc_written_file_is_journalled_as_created(tmp_path):
    snap = sc.before([tmp_path])
    proc = _run("cat > out.txt <<'EOF'\nhello\nEOF", tmp_path)
    assert proc.returncode == 0, proc.stderr
    changes = sc.after(snap)
    created = [c for c in changes if Path(c.path).name == "out.txt"]
    assert created, changes
    assert created[0].kind == "created"


def test_pipe_to_file_overwrite_is_journalled_with_pre_image(tmp_path):
    file = tmp_path / "log.txt"
    file.write_text("original\n", encoding="utf-8")
    snap = sc.before([tmp_path])
    proc = _run("echo new > log.txt", tmp_path)
    assert proc.returncode == 0, proc.stderr
    changes = sc.after(snap)
    mod = [c for c in changes if Path(c.path).name == "log.txt"]
    assert mod, changes
    assert mod[0].kind == "modified"
    assert mod[0].pre == b"original\n"
