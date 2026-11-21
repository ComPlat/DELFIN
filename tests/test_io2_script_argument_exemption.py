"""IO2 round 10: a token exempts the body only when a FILE is behind it.

Three rounds produced three leaks in the same function (`-s NAME`,
`-o/+o VALUE`, `--rcfile VALUE`), each because the gate enumerated what to
SKIP -- a hand-written approximation of bash's argument grammar, which
leaks by design. QS lead nacht-s8 ruled the shape instead of the spelling:

1. `-c`/`--command` present -> the body is an ARGUMENT (data).
2. Consume the values real bash consumes: `-o +o -O +O --rcfile
   --init-file` take one word, `--long=value` carries its own, other
   `--long` is a switch.
3. Any remaining non-option token is a SCRIPT argument -- and therefore
   makes the body that script's stdin, data -- ONLY IF IT NAMES AN
   EXISTING FILE. Otherwise the body is the PROGRAM and is judged.
4. No token at all -> the body is the program.

Point 3 is what ends the class: an unrecognised token is no longer a
licence to skip the check. It is the same fail-closed principle
`check_submit` already applies to a path that is not a file.

The deviation this buys, stated so nobody rediscovers it: `bash -eu
pipefail` with NO file named `pipefail` does not run stdin in real bash,
yet the gate now judges the body anyway. That is fail-closed on purpose
(an unknown token is not trusted), and the cost is bounded -- a body with
no finding still passes.

Both directions: 7 refused, 9 that must still run.
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


def _file(tmp_path: Path, text: str, name: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _gate(cmd: str):
    return hpc_io_check.sbatch_gate(cmd, env=dict(ENV))


# --- refused: no existing file behind the token, so the body is judged ---

def test_rcfile_value_is_consumed_so_the_body_is_judged(lustre, tmp_path):
    """Finding 5: `--rcfile FILE` takes a value; real bash then runs stdin."""
    rc = _file(tmp_path, "# rc\n", "rc")
    walk = _file(tmp_path, WALK, "walk.sh")
    verdict, text, may_run = _gate(f"bash --rcfile {rc} <<'EOF'\nsbatch {walk}\nEOF")
    assert verdict == "refuse" and not may_run
    assert "tree_walk_on_lustre" in text


def test_init_file_value_is_consumed_so_the_body_is_judged(lustre, tmp_path):
    rc = _file(tmp_path, "# rc\n", "rc")
    walk = _file(tmp_path, WALK, "walk.sh")
    verdict, text, may_run = _gate(
        f"bash --init-file {rc} <<'EOF'\nsbatch {walk}\nEOF")
    assert verdict == "refuse" and not may_run
    assert "tree_walk_on_lustre" in text


def test_unknown_token_without_a_file_fails_closed(lustre, tmp_path):
    """`bash nosuchfile.sh` -- no file, so the body is judged, not trusted.

    Real bash would fail on the missing script and run nothing; the gate
    cannot know that, so it does not treat an unknown word as a licence.
    """
    walk = _file(tmp_path, WALK, "walk.sh")
    verdict, text, may_run = _gate(
        f"bash nosuchfile.sh <<'EOF'\nsbatch {walk}\nEOF")
    assert verdict == "refuse" and not may_run
    assert "tree_walk_on_lustre" in text


def test_pipefail_without_a_file_is_judged(lustre, tmp_path):
    """The documented deviation: `-eu pipefail` + no such file -> judged."""
    walk = _file(tmp_path, WALK, "walk.sh")
    verdict, text, may_run = _gate(
        f"bash -eu pipefail <<'EOF'\nsbatch {walk}\nEOF")
    assert verdict == "refuse" and not may_run
    assert "tree_walk_on_lustre" in text


def test_bare_walk_in_an_rcfile_body_is_refused(lustre, tmp_path):
    rc = _file(tmp_path, "# rc\n", "rc")
    verdict, text, may_run = _gate(f"bash --rcfile {rc} <<'EOF'\ndu -sh /pfs/data6\nEOF")
    assert verdict == "refuse" and not may_run
    assert "tree_walk_on_lustre" in text


def test_both_rcfile_and_a_missing_script_is_judged(lustre, tmp_path):
    rc = _file(tmp_path, "# rc\n", "rc")
    walk = _file(tmp_path, WALK, "walk.sh")
    verdict, text, may_run = _gate(
        f"bash --norc --rcfile {rc} missing.sh <<'EOF'\nsbatch {walk}\nEOF")
    assert verdict == "refuse" and not may_run


def test_relative_missing_script_after_a_cd_is_judged(lustre, tmp_path):
    sub = tmp_path / "sub"
    sub.mkdir()
    walk = _file(sub, WALK, "walk.sh")
    cmd = f"cd {sub} && bash nothere.sh <<'EOF'\nsbatch {walk.name}\nEOF"
    verdict, text, may_run = _gate(cmd)
    assert verdict == "refuse" and not may_run


# --- allowed: an EXISTING file is what makes the body data --------------

def test_existing_script_argument_makes_the_body_data(lustre, tmp_path):
    clean = _file(tmp_path, CLEAN, "clean.sh")
    verdict, text, may_run = _gate(
        f"bash {clean} <<'EOF'\nsbatch /pfs/data6/nope.sh\nEOF")
    assert verdict == "allow" and may_run, text


def test_existing_script_after_an_option_makes_the_body_data(lustre,
                                                            tmp_path):
    clean = _file(tmp_path, CLEAN, "clean.sh")
    verdict, text, may_run = _gate(
        f"bash -x {clean} <<'EOF'\nsbatch /pfs/data6/nope.sh\nEOF")
    assert verdict == "allow" and may_run, text


def test_existing_pipefail_file_makes_the_body_data(lustre, tmp_path):
    """With the file present, `-eu pipefail` really does run that script."""
    pf = _file(tmp_path, CLEAN, "pipefail")
    verdict, text, may_run = _gate(
        f"bash -eu {pf} <<'EOF'\nsbatch /pfs/data6/nope.sh\nEOF")
    assert verdict == "allow" and may_run, text


def test_existing_script_after_rcfile_value_is_data(lustre, tmp_path):
    rc = _file(tmp_path, "# rc\n", "rc")
    clean = _file(tmp_path, CLEAN, "clean.sh")
    verdict, text, may_run = _gate(
        f"bash --rcfile {rc} {clean} <<'EOF'\nsbatch /pfs/data6/nope.sh\nEOF")
    assert verdict == "allow" and may_run, text


def test_existing_script_after_double_dash_is_data(lustre, tmp_path):
    clean = _file(tmp_path, CLEAN, "clean.sh")
    verdict, text, may_run = _gate(
        f"bash -- {clean} <<'EOF'\nsbatch /pfs/data6/nope.sh\nEOF")
    assert verdict == "allow" and may_run, text


def test_rcfile_body_with_a_clean_submit_still_runs(lustre, tmp_path):
    rc = _file(tmp_path, "# rc\n", "rc")
    clean = _file(tmp_path, CLEAN, "clean.sh")
    verdict, text, may_run = _gate(
        f"bash --rcfile {rc} <<'EOF'\nsbatch {clean}\nEOF")
    assert verdict == "allow" and may_run, text


def test_command_argument_is_still_data(lustre, tmp_path):
    walk = _file(tmp_path, WALK, "walk.sh")
    verdict, text, may_run = _gate(f"bash -c 'sbatch {walk}' <<'EOF'\ndata\nEOF")
    assert verdict == "refuse" and not may_run


def test_readers_and_plain_shell_bodies_unchanged(lustre):
    assert _gate("python3 - <<'PY'\nprint(1)\nPY")[0] == "allow"
    assert _gate("cat <<'NOTE'\ndu -sh /pfs/data6\nNOTE")[0] == "allow"
    assert _gate("bash <<'EOF'\necho hello\nEOF")[0] == "allow"
