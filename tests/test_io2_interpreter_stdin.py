"""IO2 round 10, QS finding 2: a shell fed its script on stdin runs it.

`bash <<'EOF' … EOF` is the SAME program as `bash -c '…'`, and
``_strip_wrappers`` already unwraps the ``-c`` spelling, so the gate
refused that one and allowed its here-document twin: a submit or a tree
walk placed in a body handed to a shell reached sbatch. The module
already holds the opposite premise elsewhere — ``check()`` (the rule set
IO1 wired into ``cluster_submit``) finds a ``du -sh /pfs/…`` in such a
body, and ``sbatch_gate`` returned allow for the same text (QS lead
nacht-s8, round 10).

The rule this file pins: an INTERPRETER THAT READS ITS STDIN AS SHELL
COMMANDS (bash, sh, dash, zsh, ksh, also as an absolute path, also via
``env``, also with ``-s``) is judged by its body — the body's own
submits are gated and the body's text goes through the metadata rules.
Every other reader keeps its old behaviour, because a body is DATA to it:
``python3 - <<'PY'`` and a harmless shell body must still run, and the
shell_words note (delfin/agent/shell_words.py:241) records why a python
body must never be judged as shell.

Both directions are pinned: 6 shapes refused, 5 that must still run.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin import hpc_io_check

LUSTRE_MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/home/x"}
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"
CLEAN = "#!/bin/bash\nexport ORCA_SCRATCH=/tmp/scr\n"


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read",
                        lambda self: LUSTRE_MOUNTS)


def _script(tmp_path: Path, text: str, name: str = "walk.sh") -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _gate(cmd: str):
    return hpc_io_check.sbatch_gate(cmd, env=dict(ENV))


# --- refused: the body of a shell IS that shell's script ----------------

def test_bash_heredoc_submitting_a_walking_script_is_refused(lustre, tmp_path):
    script = _script(tmp_path, WALK)
    cmd = f"bash <<'EOF'\nsbatch {script}\nEOF"
    verdict, text, may_run = _gate(cmd)
    assert verdict == "refuse" and not may_run
    assert "tree_walk_on_lustre" in text


def test_sh_s_heredoc_submitting_a_walking_script_is_refused(lustre, tmp_path):
    script = _script(tmp_path, WALK)
    cmd = f"sh -s <<'EOF'\nsbatch {script}\nEOF"
    verdict, text, may_run = _gate(cmd)
    assert verdict == "refuse" and not may_run
    assert "tree_walk_on_lustre" in text


def test_absolute_path_and_env_forms_are_the_same_shell(lustre, tmp_path):
    script = _script(tmp_path, WALK)
    for head in ("/bin/bash", "/bin/sh", "env bash", "env -i bash", "dash",
                 "zsh", "ksh"):
        cmd = f"{head} <<'EOF'\nsbatch {script}\nEOF"
        verdict, text, may_run = _gate(cmd)
        assert verdict == "refuse" and not may_run, head
        assert "tree_walk_on_lustre" in text, head


def test_bash_heredoc_walk_is_refused_like_check_finds_it(lustre):
    """The row where the module contradicted itself: check() refused the
    body while the gate said allow for the same text."""
    body = "du -sh /pfs/data6/workspace/\n"
    findings = hpc_io_check.check(
        body, env=dict(ENV),
        mounts_path=hpc_io_check.MountsSource(text=LUSTRE_MOUNTS))
    assert any(f.rule == "tree_walk_on_lustre" for f in findings)
    verdict, text, may_run = _gate(f"bash <<'EOF'\n{body}EOF")
    assert verdict == "refuse" and not may_run
    assert "tree_walk_on_lustre" in text


def test_bash_heredoc_submit_relative_after_cd_in_the_body(lustre, tmp_path):
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "walk.sh").write_text(WALK, encoding="utf-8")
    cmd = f"bash <<'EOF'\ncd {sub}\nsbatch walk.sh\nEOF"
    verdict, text, may_run = _gate(cmd)
    assert verdict == "refuse" and not may_run
    assert "tree_walk_on_lustre" in text


def test_nested_interpreter_heredoc_is_judged(lustre, tmp_path):
    """`bash <<A` whose body starts `sh <<B` still runs commands."""
    script = _script(tmp_path, WALK)
    cmd = f"bash <<'A'\nsh <<'B'\nsbatch {script}\nB\nA"
    verdict, text, may_run = _gate(cmd)
    assert verdict == "refuse" and not may_run
    assert "tree_walk_on_lustre" in text


# --- allowed: data stays data, clean work keeps running ----------------

def test_python_heredoc_still_runs(lustre):
    """A body is DATA to python -- shell_words.py:241 records the bug that
    came from judging it as shell."""
    cmd = "python3 - <<'PY'\nprint('du -sh /pfs/data6/workspace/')\nPY"
    verdict, text, may_run = _gate(cmd)
    assert verdict == "allow" and may_run, text


def test_harmless_shell_body_still_runs(lustre):
    cmd = "bash <<'EOF'\necho hello\nexport X=/tmp/y\nEOF"
    verdict, text, may_run = _gate(cmd)
    assert verdict == "allow" and may_run, text


def test_shell_body_with_a_clean_submit_still_runs(lustre, tmp_path):
    script = _script(tmp_path, CLEAN, name="clean.sh")
    cmd = f"bash <<'EOF'\nsbatch {script}\nEOF"
    verdict, text, may_run = _gate(cmd)
    assert verdict == "allow" and may_run, text


def test_bash_dash_c_spelling_unchanged(lustre, tmp_path):
    script = _script(tmp_path, WALK)
    verdict, text, may_run = _gate(f"bash -c 'sbatch {script}'")
    assert verdict == "refuse" and not may_run


def test_shell_heredoc_of_data_for_sbatch_is_still_allowed(lustre, tmp_path):
    """`sbatch <file> <<EOF` submits the FILE; a shell head is not present,
    so the body stays the command's stdin and is not the script."""
    script = _script(tmp_path, CLEAN, name="clean.sh")
    cmd = f"sbatch {script} <<'EOF'\nsbatch /pfs/data6/should-not-matter.sh\nEOF"
    verdict, text, may_run = _gate(cmd)
    assert verdict == "allow" and may_run, text
