"""IO2 round 10: a newline and a here-document may not hide a submit.

Found by the QS lead (nacht-s8) on the draft PR: the bash gate splits a
command line on ``;`` ``&&`` ``|`` ``||`` ``&`` at quote depth 0, but not
on a newline, so ``echo start`` + newline + ``sbatch <walk.sh>`` reached
sbatch with the gate reporting nothing. And a here-document handed to
sbatch (``sbatch <<'EOF' ... EOF``) IS that job's script -- judged as
script text, exactly like ``--wrap``, never split into fake segments.

Both directions are pinned here: the shapes that must be refused and the
multi-line / here-document shapes that are ordinary work and must still
run (a gate that refuses those would be a defect of its own).

Every case drives the public decision surface -- ``sbatch_gate`` -- which
is the function the protected ``api_client.py`` hook calls.
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


# --- refused: the submit is hidden behind a newline or a heredoc ---------

def test_newline_before_sbatch_is_refused(lustre, tmp_path):
    """The QS case verbatim: `echo start` newline `sbatch <abs>/walk.sh`."""
    script = _script(tmp_path, WALK)
    verdict, text, may_run = _gate(f"echo start\nsbatch {script}")
    assert verdict == "refuse"
    assert not may_run, "a refused sbatch command must not run"
    assert "tree_walk_on_lustre" in text


def test_indented_second_line_submit_is_refused(lustre, tmp_path):
    """A continuation-formatted command is still two commands."""
    script = _script(tmp_path, WALK)
    verdict, text, may_run = _gate(f"echo start\n  sbatch {script}")
    assert verdict == "refuse" and not may_run
    assert "tree_walk_on_lustre" in text


def test_submit_on_a_later_line_of_three_is_refused(lustre, tmp_path):
    script = _script(tmp_path, WALK)
    cmd = f"mkdir -p out\ncd {tmp_path}\nsbatch {script.name}"
    verdict, text, may_run = _gate(cmd)
    assert verdict == "refuse" and not may_run
    assert "tree_walk_on_lustre" in text


def test_heredoc_script_is_checked_as_text(lustre, tmp_path):
    """`sbatch <<EOF` submits the BODY: it is script text, not segments."""
    cmd = ("sbatch <<'EOF'\n"
           "#!/bin/bash\n"
           "du -sh /pfs/data6/workspace/\n"
           "EOF")
    verdict, text, may_run = _gate(cmd)
    assert verdict == "refuse" and not may_run
    assert "tree_walk_on_lustre" in text


def test_heredoc_script_with_options_is_checked_as_text(lustre, tmp_path):
    """Options on the heredoc's own line do not turn the body into data."""
    cmd = ("sbatch -J io2 <<EOF\n"
           "#!/bin/bash\n"
           "find /pfs/data6/workspace/ -type f\n"
           "EOF")
    verdict, text, may_run = _gate(cmd)
    assert verdict == "refuse" and not may_run
    assert "tree_walk_on_lustre" in text


# --- allowed: ordinary multi-line work and non-script heredocs ----------

def test_multiline_clean_command_is_allowed(lustre):
    verdict, text, may_run = _gate("echo a\necho b")
    assert verdict == "allow" and may_run
    assert text == ""


def test_multiline_clean_submit_is_allowed(lustre, tmp_path):
    script = _script(tmp_path, CLEAN, name="clean.sh")
    verdict, text, may_run = _gate(f"echo start\nsbatch {script}")
    assert verdict == "allow" and may_run, text


def test_backslash_continuation_stays_one_command(lustre, tmp_path):
    """A `sbatch \\` + newline + script line is ONE command; refused."""
    script = _script(tmp_path, WALK)
    verdict, text, may_run = _gate(f"sbatch \\\n  {script}")
    assert verdict == "refuse" and not may_run
    assert "tree_walk_on_lustre" in text


def test_clean_heredoc_script_is_allowed(lustre):
    cmd = ("sbatch <<'EOF'\n"
           "#!/bin/bash\n"
           "export ORCA_SCRATCH=/tmp/scr\n"
           "EOF")
    verdict, text, may_run = _gate(cmd)
    assert verdict == "allow" and may_run, text


def test_heredoc_data_for_a_reader_is_not_a_script(lustre):
    """`python - <<'PY'` feeds an interpreter, not sbatch: no findings."""
    cmd = ("python3 - <<'PY'\n"
           "print('du -sh /pfs/data6/workspace/')\n"
           "PY")
    verdict, text, may_run = _gate(cmd)
    assert verdict == "allow" and may_run, text


def test_space_form_wrap_text_is_judged(lustre):
    """`sbatch --wrap '<walk>'` (space form) hides the same walk.

    Pre-existing defect found by my own adversarial pass in round 10: the
    wrap regex captured the surrounding quotes, so the text arrived as
    "'du -sh /pfs/...'" and the rule never saw the command.
    """
    cmd = "sbatch --wrap 'du -sh /pfs/data6/workspace/'"
    verdict, text, may_run = _gate(cmd)
    assert verdict == "refuse" and not may_run
    assert "tree_walk_on_lustre" in text


def test_space_form_clean_wrap_is_allowed(lustre):
    cmd = "sbatch --wrap 'echo hello'"
    verdict, text, may_run = _gate(cmd)
    assert verdict == "allow" and may_run, text


def test_heredoc_then_a_clean_submit_is_allowed(lustre, tmp_path):
    script = _script(tmp_path, CLEAN, name="clean.sh")
    cmd = ("cat <<'NOTE' > note.txt\n"
           "du -sh /pfs/data6/workspace/\n"
           "NOTE\n"
           f"sbatch {script}")
    verdict, text, may_run = _gate(cmd)
    assert verdict == "allow" and may_run, text
