"""IO2 round 10, QS finding 3: an option's VALUE is not a script argument.

`_shell_command_head` exists so that `bash <<'EOF' ... EOF` is judged as
the program it is. It returned None on the first token that is not an
option, calling that token a script argument — but two token kinds are
NOT script arguments while the shell is still reading stdin:

- the VALUE of an option that takes one: ``-o posix``, ``+o posix``,
  ``-O extglob``, ``+O extglob``;
- the ``$0`` NAME after ``-s``/``--stdin``: ``bash -s myname``.

So ``bash -o posix <<'EOF' … sbatch walk.sh … EOF`` — one line anyone can
type — reached sbatch unchecked, while bare ``bash`` was refused. QS lead
nacht-s8 found it (their probe drove REAL bash with the body on a pipe
and a marker file as the witness), and the six rows below were red on
10bc22c19 before this file's fix.

The counter-row is why the rule cannot be "any token after an option is
an option": ``bash -eu pipefail`` runs a SCRIPT called ``pipefail``,
never stdin, and must stay allowed.

Both directions: 8 shapes refused, 5 that must still run.
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


# --- the six rows that were red: an option VALUE still means stdin -------

@pytest.mark.parametrize("head", [
    "bash -s myname",     # $0 name after -s
    "bash -o posix",      # -o takes a value
    "bash +o posix",      # +o takes a value
    "bash -O extglob",    # -O takes a value
    "bash +O extglob",    # +O takes a value
    "bash --stdin name",  # the long spelling of -s, then $0
])
def test_option_value_does_not_cancel_stdin(lustre, tmp_path, head):
    script = _script(tmp_path, WALK)
    cmd = f"{head} <<'EOF'\nsbatch {script}\nEOF"
    verdict, text, may_run = _gate(cmd)
    assert verdict == "refuse" and not may_run, head
    assert "tree_walk_on_lustre" in text, head


def test_option_value_forms_allow_a_bare_walk(lustre):
    """No submit at all, same three rows: the walk is still the finding."""
    for head in ("bash -o posix", "bash +o posix", "bash -s myname"):
        verdict, text, may_run = _gate(f"{head} <<'EOF'\ndu -sh /pfs/data6\nEOF")
        assert verdict == "refuse" and not may_run, head
        assert "tree_walk_on_lustre" in text, head


# --- what must keep working ---------------------------------------------

def test_a_real_script_argument_still_wins(lustre, tmp_path):
    """`bash <script> <<EOF` runs the FILE: the body is stdin data."""
    clean = _script(tmp_path, CLEAN, name="clean.sh")
    cmd = f"bash {clean} <<'EOF'\nsbatch /pfs/data6/nope.sh\nEOF"
    verdict, text, may_run = _gate(cmd)
    assert verdict == "allow" and may_run, text


def test_script_argument_after_an_option_still_wins(lustre, tmp_path):
    clean = _script(tmp_path, CLEAN, name="clean.sh")
    cmd = f"bash -x {clean} <<'EOF'\nsbatch /pfs/data6/nope.sh\nEOF"
    verdict, text, may_run = _gate(cmd)
    assert verdict == "allow" and may_run, text


def test_pipefail_is_a_script_name_not_an_option_value(lustre):
    """The counter-row: `bash -eu pipefail` runs a script called pipefail.

    This is why the rule cannot be "every token after an option is an
    option value" — that would make a script argument look like stdin.
    """
    cmd = "bash -eu pipefail <<'EOF'\nsbatch /pfs/data6/nope.sh\nEOF"
    verdict, text, may_run = _gate(cmd)
    assert verdict == "allow" and may_run, text


def test_python_and_cat_guards_unchanged(lustre):
    assert _gate("python3 - <<'PY'\nprint(1)\nPY")[0] == "allow"
    assert _gate("cat <<'NOTE'\ndu -sh /pfs/data6\nNOTE")[0] == "allow"


def test_clean_shell_body_still_runs(lustre, tmp_path):
    clean = _script(tmp_path, CLEAN, name="clean.sh")
    verdict, text, may_run = _gate(f"bash -o posix <<'EOF'\nsbatch {clean}\nEOF")
    assert verdict == "allow" and may_run, text
