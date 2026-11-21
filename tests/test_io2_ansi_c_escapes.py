r"""IO2 round 16: numeric escapes in ANSI-C quoting hide a script file.

Measured at de7a21c8d: `sbatch $'\\\\x2ftmp/.../walk.sh'` is ALLOWED while
the same file by its plain path refuses. bash expands `\\x2f` to `/`
(confirmed: `printf '%s' $'\\x2f'` prints `/`), so sbatch receives a real,
walking file -- and the gate's name test sees the literal `$'\x2ftmp/...'`,
which is no file, so the deliberate missing-name pass allows.

This is the same fail-open as round 15 in a different escape, which is why
the fix is not "add \x to the list": the DETECTABLE fact is whether the
text the shell would hand sbatch names a real file. Two directions are
pinned so the fail-closed side cannot swallow legitimate submits.
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check

MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/home/x"}
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"
CLEAN = "#!/bin/bash\necho hello\n"


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read", lambda self: MOUNTS)


def _v(cmd):
    return hpc_io_check.sbatch_gate(cmd, env=dict(ENV))


def _hex_quoted(path: str) -> str:
    """The path as an ANSI-C word with its leading slash as \\x2f."""
    return "$'\\x2f" + path.lstrip("/") + "'"


def test_a_hex_escaped_walking_file_is_refused(lustre, tmp_path):
    walk = tmp_path / "walk.sh"
    walk.write_text(WALK, encoding="utf-8")
    cmd = f"sbatch {_hex_quoted(str(walk))}"
    verdict, text, may_run = _v(cmd)
    assert verdict == "refuse", (
        f"FAIL-OPEN: bash expands \\x2f, so sbatch runs the walking file, "
        f"but the gate allowed it (verdict={verdict!r})")
    assert may_run is False


def test_a_hex_escaped_walking_file_with_options_is_refused(lustre, tmp_path):
    walk = tmp_path / "walk.sh"
    walk.write_text(WALK, encoding="utf-8")
    cmd = f"sbatch --test-only -p dev_cpu {_hex_quoted(str(walk))}"
    assert _v(cmd)[0] == "refuse"


def test_a_hex_escaped_clean_file_still_runs(lustre, tmp_path):
    clean = tmp_path / "clean.sh"
    clean.write_text(CLEAN, encoding="utf-8")
    cmd = f"sbatch {_hex_quoted(str(clean))}"
    verdict, text, may_run = _v(cmd)
    assert verdict == "allow", (
        f"OVER-REFUSAL: a clean file named by a \\x2f escape was refused: "
        f"{text[:160]}")
    assert may_run is True


def test_a_hex_escape_that_names_nothing_is_still_allowed(lustre, tmp_path):
    """The deliberate missing-name pass survives: sbatch fails on it."""
    cmd = f"sbatch {_hex_quoted(str(tmp_path / 'nosuch.sh'))}"
    verdict, text, may_run = _v(cmd)
    assert verdict == "allow", (
        f"a name that exists nowhere should go to sbatch, which fails; "
        f"got {verdict!r}")
    assert may_run is True


def test_plain_paths_are_unchanged(lustre, tmp_path):
    walk = tmp_path / "walk.sh"
    walk.write_text(WALK, encoding="utf-8")
    clean = tmp_path / "clean.sh"
    clean.write_text(CLEAN, encoding="utf-8")
    assert _v(f"sbatch {walk}")[0] == "refuse"
    assert _v(f"sbatch {clean}")[0] == "allow"


# --------------------------------------------------------------------------
# the expansion must COMPOSE with a cd target, not bypass it
# --------------------------------------------------------------------------

@pytest.mark.parametrize("tail,expected", [
    ("sbatch walk.sh", "refuse"),
    ("sbatch $'walk.sh'", "refuse"),
    ("sbatch $'\\x77alk.sh'", "refuse"),
    ("sbatch clean.sh", "allow"),
    ("sbatch $'clean.sh'", "allow"),
    ("sbatch $'\\x63lean.sh'", "allow"),
], ids=["plain_walk", "ansi_c_walk", "hex_walk",
        "plain_clean", "ansi_c_clean", "hex_clean"])
def test_a_cd_target_composes_with_the_expansion(lustre, tmp_path, tail,
                                                 expected):
    """`cd X && sbatch $'<relative>'` resolves the EXPANDED name in X.

    A refuse is only meaningful if it names the CONTENT: the resolution
    refusal ("not found from the cd target") would also refuse a walk that
    was never resolved, so the reason is asserted too (round-11 lesson).
    """
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "walk.sh").write_text(WALK, encoding="utf-8")
    (sub / "clean.sh").write_text(CLEAN, encoding="utf-8")
    verdict, text, _may_run = _v(f"cd {sub} && {tail}")
    assert verdict == expected, (
        f"{tail!r} from the cd target: expected {expected}, got {verdict!r} "
        f"({text[:120]})")
    if expected == "refuse":
        assert "tree_walk_on_lustre" in text, (
            f"the refusal must come from the CONTENT, not from a failed "
            f"resolution: {text[:160]}")
        assert "not found from the cd target" not in text
