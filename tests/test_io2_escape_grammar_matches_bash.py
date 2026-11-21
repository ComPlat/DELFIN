"""The gate's decoding of a script word must agree with what bash produces.

Why a property test and not another list of spellings: three escape families
were found by inspection in three consecutive rounds (`$'...'` itself, then
`\\xHH`/`\\NNN`, then `\\uHHHH`/`\\UHHHHHHHH`), each discovered because someone
thought of it. A list of families is a list of one's imagination. This file
enumerates a space instead, asks a REAL bash what each word denotes, and
requires the gate to agree.

The property, stated once:

    for every spelling W of an existing walking script's path,
    if `bash -c "printf %s W"` yields that path,
    then `sbatch W` must be REFUSED.

The converse is deliberately NOT asserted: when bash yields something other
than the walking file, other rules decide the verdict, and pinning them here
would turn this into a spelling list again.

A row whose bash expansion does not reproduce the path is INVALID and fails
the test: an invalid row measures nothing, which is the mistake this file
exists to avoid.
"""
from __future__ import annotations

import subprocess

import pytest

from delfin import hpc_io_check

MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/home/x"}
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"

# Every way bash lets one write the path separator inside $'...'. The list is
# a space to enumerate, not a list of defects: a form the decoder misses shows
# up here as a red row, which is the point.
ESCAPED_SEPARATOR = [
    r"\x2f",          # hex, lower
    r"\x2F",          # hex, upper
    r"\057",          # octal
    r"\u002f",        # unicode, 4 hex digits
    r"\U0000002f",    # unicode, 8 hex digits
]

# MEASURED AND EXCLUDED: `\/` is NOT a separator escape inside $'...'.
# bash keeps the backslash (`$'a\/b'` -> `a\/b`), so a row using it does not
# produce the walking path and measures nothing. The invalid-row assertion
# below is what caught this: it is the check that stops this file decaying
# into a list of spellings that "look like they should work".


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read", lambda self: MOUNTS)


def _bash_says(word: str) -> str:
    """What bash makes of this word, as text. Real binary, no re-implementation."""
    proc = subprocess.run(["bash", "-c", f"printf '%s' {word}"],
                          capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, (
        f"bash could not evaluate {word!r}: rc={proc.returncode} "
        f"stderr={proc.stderr[:200]!r}")
    return proc.stdout


def _rows(path: str) -> list[tuple[str, str]]:
    """(label, word) for the plain path and every escaped separator form."""
    rows = [("plain", path)]
    for form in ESCAPED_SEPARATOR:
        rows.append((f"$'{form}'", "$'" + path.replace("/", form) + "'"))
    return rows


@pytest.mark.parametrize("form", ["plain"] + ESCAPED_SEPARATOR,
                         ids=["plain"] + [f.replace("\\", "") for f
                                          in ESCAPED_SEPARATOR])
def test_a_word_bash_resolves_to_the_walking_file_is_refused(lustre, tmp_path,
                                                             form):
    script = tmp_path / "walk.sh"
    script.write_text(WALK, encoding="utf-8")
    path = str(script)

    word = path if form == "plain" else "$'" + path.replace("/", form) + "'"

    # The row is only a measurement if bash really produces the walking file.
    assert _bash_says(word) == path, (
        f"INVALID ROW [{form}]: bash yields {_bash_says(word)!r}, not the "
        f"existing walking script — this row measures nothing and must be "
        f"fixed or removed, not asserted")

    verdict, text, may_run = hpc_io_check.sbatch_gate(f"sbatch {word}",
                                                      env=dict(ENV))
    assert verdict == "refuse", (
        f"FAIL-OPEN [{form}]: bash resolves this word to the existing walking "
        f"script {path!r}, but the gate allows it (verdict={verdict!r}, "
        f"text={text[:120]!r})")
    assert may_run is False, f"[{form}] refused but may_run={may_run!r}"


def test_the_row_set_is_not_empty_and_covers_the_plain_path(tmp_path):
    """Guard on the generator itself: an empty or degenerate set measures nothing."""
    script = tmp_path / "walk.sh"
    path = str(script)
    rows = _rows(path)
    assert len(rows) == 1 + len(ESCAPED_SEPARATOR)
    assert rows[0] == ("plain", path)
    labels = [label for label, _ in rows]
    assert len(set(labels)) == len(labels), "duplicate labels hide a row"
