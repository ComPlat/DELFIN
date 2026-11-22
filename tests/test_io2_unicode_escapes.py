r"""IO2 round 17: unicode ANSI-C escapes hide a script file too.

Measured at b0202f6ef: `sbatch $'\u002ftmp/.../walk.sh'` and the `\U`
form are both ALLOWED while the plain path refuses. bash expands them to
`/` (verified: `printf '%s' $'\u002f'` prints `/`), so sbatch receives a
real walking file; the gate's name test sees the literal escape, finds no
file, and the deliberate missing-name pass allows.

This is the THIRD escape family in the same fail-open (round 15: the
`$'...'` quoting itself; round 16: `\xHH`/`\NNN`; now `\uHHHH`/`\UHHHHHHHH`).
Each earlier fix was correct for its family and the next family was still
open, which is the clearest evidence in this package that decoding must be
COMPLETE for bash's escape grammar rather than added to case by case.

Two directions pinned, so the fix cannot buy safety with over-refusal.
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check

MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/home/x"}
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"
CLEAN = "#!/bin/bash\necho hello\n"
ESC = "0x2f"


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read", lambda self: MOUNTS)


def _v(cmd):
    return hpc_io_check.sbatch_gate(cmd, env=dict(ENV))


def _u_quoted(path: str, form: str = "u") -> str:
    body = path.lstrip("/")
    if form == "u":
        return "$'\\u002f" + body + "'"
    return "$'\\U0000002f" + body + "'"


@pytest.mark.parametrize("form", ["u", "U"], ids=["u4", "U8"])
def test_a_unicode_escaped_walking_file_is_refused(lustre, tmp_path, form):
    walk = tmp_path / "walk.sh"
    walk.write_text(WALK, encoding="utf-8")
    cmd = f"sbatch {_u_quoted(str(walk), form)}"
    verdict, text, may_run = _v(cmd)
    assert verdict == "refuse", (
        f"FAIL-OPEN: bash expands the unicode escape, so sbatch runs the "
        f"walking file, but the gate allowed it (verdict={verdict!r})")
    assert may_run is False


@pytest.mark.parametrize("form", ["u", "U"], ids=["u4", "U8"])
def test_a_unicode_escaped_clean_file_still_runs(lustre, tmp_path, form):
    clean = tmp_path / "clean.sh"
    clean.write_text(CLEAN, encoding="utf-8")
    cmd = f"sbatch {_u_quoted(str(clean), form)}"
    verdict, text, may_run = _v(cmd)
    assert verdict == "allow", (
        f"OVER-REFUSAL: a clean file named by a unicode escape was "
        f"refused: {text[:160]}")
    assert may_run is True


def test_a_unicode_escape_that_names_nothing_is_still_allowed(lustre,
                                                              tmp_path):
    """The deliberate missing-name pass survives."""
    cmd = f"sbatch {_u_quoted(str(tmp_path / 'nosuch.sh'))}"
    verdict, text, may_run = _v(cmd)
    assert verdict == "allow", (
        f"a name that exists nowhere should go to sbatch, which fails; "
        f"got {verdict!r}")
    assert may_run is True


def test_dollar_double_quote_is_ordinary_quoting(lustre, tmp_path):
    r"""`$"..."` is locale translation, not an escape: bash passes the path.

    Measured: bash yields the real path, so the gate must refuse a walking
    file named that way. Pinned because it is the one `$`-prefixed form my
    expansion must NOT mangle.
    """
    walk = tmp_path / "walk.sh"
    walk.write_text(WALK, encoding="utf-8")
    cmd = 'sbatch $"' + str(walk) + '"'
    assert _v(cmd)[0] == "refuse"


def test_dotdot_traversal_from_a_cd_target_is_refused(lustre, tmp_path):
    outer = tmp_path / "outer"
    inner = outer / "inner"
    inner.mkdir(parents=True)
    (outer / "walk.sh").write_text(WALK, encoding="utf-8")
    cmd = f"cd {inner} && sbatch ../walk.sh"
    verdict, text, _may_run = _v(cmd)
    assert verdict == "refuse", (
        f"`sbatch ../walk.sh` from a cd target was allowed: {verdict!r}")
    assert "tree_walk_on_lustre" in text, (
        f"the refusal must come from the CONTENT: {text[:160]}")
