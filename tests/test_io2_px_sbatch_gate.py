"""IO2 phase 3 (px): the agent's own bash route checks `sbatch` scripts.

Route (b) of the package's reason for existing: an `sbatch <script>`
the agent types into its bash tool submitted with no Lustre check.
The hook lives in `api_client.py` (PROTECTED — the .gate/px pair is
built by the operator); this red test pins the CONTRACT the px must
satisfy, against a small reference implementation of the hook's
decision logic (`_io2_sbatch_gate` in delfin/hpc_io_check.py — my
file) that the px will call.

Red on the previous commit: the reference hook does not exist yet, so
every case errors; after the helper hook exists but BEFORE the px is
applied, the api_client-side tests fail differently (they drive the
real gate decision function directly — the operator builds the px from
the copies, and its own wiring test is theirs).

Cases (reviewer nacht-s2's measured shapes):
- first word `sbatch <script>` -> refused when the script walks Lustre
- `sbatch` after `;`, `&&`, `|`, `env X=1`, `exec`, `nice` -> caught
- script argument attached: `-o file`, `-ofile`
- `--wrap=<text>` -> the wrapped text is checked as text
- refuse -> the command is NOT run and the finding names it
- warn -> the command runs and the warning is in the result
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from delfin import hpc_io_check


LUSTRE_MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"
SCRATCH = "#!/bin/bash\nexport ORCA_SCRATCH=/pfs/data6/scr/\n"


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read",
                        lambda self: LUSTRE_MOUNTS)


def _script(tmp_path: Path, text: str, name: str = "run.sh") -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _gate(cmd: str, tmp_path: Path, run: bool = True):
    """The contract surface: the hook the px calls. Returns
    ``(verdict, text, ran)`` where verdict is "allow" | "refuse" |
    "warn", text is the finding text ("" when allow) and ran says
    whether the command may execute."""
    from delfin import hpc_io_check
    decision = hpc_io_check.sbatch_gate(cmd, env=dict(_ENV))
    return decision


_ENV = {"HOME": "/home/x"}


def test_first_word_sbatch_script_refused(lustre, tmp_path):
    script = _script(tmp_path, WALK)
    verdict, text, ran = _gate(f"sbatch {script}", tmp_path)
    assert verdict == "refuse"
    assert not ran, "a refused sbatch command must not run"
    assert "tree_walk_on_lustre" in text


def test_sbatch_after_separators_and_wrappers_is_caught(lustre, tmp_path):
    script = _script(tmp_path, WALK)
    for prefix in ("cd /tmp && ", "true; ", "echo x | ", "env X=1 ",
                   "exec ", "nice "):
        verdict, _, ran = _gate(f"{prefix} sbatch {script}", tmp_path)
        assert verdict == "refuse", prefix
        assert not ran, prefix


def test_attached_script_argument(lustre, tmp_path):
    script = _script(tmp_path, WALK)
    for spelling in (f"sbatch -o out.txt {script}",
                     f"sbatch -o{script}"):
        if spelling.endswith(str(script)) or "-o" + str(script) in spelling:
            verdict, text, ran = _gate(spelling, tmp_path)
            # -o names the OUTPUT file, not the script; the script is
            # the last positional argument or the token after -o.
            # Both spellings here name the walking script as the job
            # script, so both must refuse.
            assert verdict == "refuse", spelling


def test_wrap_text_is_checked_as_text(lustre):
    verdict, text, ran = _gate(
        "sbatch --wrap='du -sh /pfs/data6/workspace/'", None)
    assert verdict == "refuse"
    assert "tree_walk_on_lustre" in text


def test_warn_runs_and_says_so(lustre, tmp_path):
    script = _script(tmp_path, SCRATCH)
    verdict, text, ran = _gate(f"sbatch {script}", tmp_path)
    assert verdict == "warn"
    assert ran, "a warning must not stop the submit"
    assert "scratch_on_lustre" in text


def test_clean_script_runs(lustre, tmp_path):
    script = _script(tmp_path, "#!/bin/bash\necho hello\n")
    verdict, text, ran = _gate(f"sbatch {script}", tmp_path)
    assert verdict == "allow"
    assert ran


def test_non_sbatch_commands_pass_through(lustre, tmp_path):
    verdict, _, ran = _gate("ls -la /tmp", tmp_path)
    assert verdict == "allow" and ran
