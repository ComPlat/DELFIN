"""IO2 round 14: the QS's here-string matrix, including ANSI-C quoting.

QS s8 measured these six shapes at f13ec9f91 and four were ALLOW:
`$'...'` (ANSI-C quoting) and options AFTER the redirect both defeated the
reader. ANSI-C quoting is how anyone writes a multi-line script inline, so
it is exactly the shape the stdin channel has to catch.

Under ruling (a) every row refuses, whatever the quoting and wherever the
options sit, because the channel is judged rather than the text.
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check

MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/home/x"}
WALK = "du -sh /pfs/data6/workspace/"
WALK_SCRIPT = "#!/bin/bash\\ndu -sh /pfs/data6/workspace/\\n"


@pytest.fixture
def lustre(monkeypatch):
    monkeypatch.setattr(hpc_io_check.MountsSource, "read", lambda self: MOUNTS)


def _v(cmd):
    return hpc_io_check.sbatch_gate(cmd, env=dict(ENV))


MATRIX = [
    ("plain_double_quoted", f'sbatch <<< "{WALK}"'),
    ("ansi_c_no_options", f"sbatch <<< $'{WALK}'"),
    ("options_before_double", f'sbatch --test-only -p dev_cpu -t 5 <<< "{WALK}"'),
    ("options_before_ansi_c",
     f"sbatch --test-only -p dev_cpu -t 5 <<< $'{WALK}'"),
    ("ansi_c_multiline_script",
     f"sbatch --test-only -p dev_cpu -t 5 <<< $'{WALK_SCRIPT}'"),
    ("options_after", f'sbatch <<< "{WALK}" --test-only -p dev_cpu'),
]


@pytest.mark.parametrize("name,cmd", MATRIX, ids=[c[0] for c in MATRIX])
def test_every_here_string_shape_is_refused(lustre, name, cmd):
    verdict, text, may_run = _v(cmd)
    assert verdict == "refuse", (
        f"FAIL-OPEN [{name}]: the here-string channel was allowed "
        f"(verdict={verdict!r})")
    assert may_run is False
    assert "job.sh" in text, (
        f"the refusal must name the file route: {text[:160]}")


@pytest.mark.parametrize("name,cmd", MATRIX, ids=[c[0] for c in MATRIX])
def test_the_same_shapes_with_clean_text_are_also_refused(lustre, name, cmd):
    """The channel is refused, so the text does not rescue it."""
    clean = cmd.replace(WALK, "echo hi").replace(
        WALK_SCRIPT, "#!/bin/bash\\\\necho hi\\\\n")
    verdict, text, may_run = _v(clean)
    assert verdict == "refuse", (
        f"[{name}] a clean here-string must still refuse under ruling (a); "
        f"got {verdict!r}")
    assert may_run is False
