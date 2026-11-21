"""A `cd` inside a SUBSHELL reaches a relative script argument.

The cd-propagation rule resolves a relative `sbatch x.sh` against the target
of a `cd` that precedes it in the same command. A subshell -- `( cd X;
sbatch y.sh )` -- separates the two by parentheses, which is a different
shape from `cd X && sbatch y.sh` and `cd X; sbatch y.sh`.

Both directions are pinned:
  * a WALKING script inside a subshell must be REFUSED, otherwise the
    exemption would be a way to hide a walk behind parentheses;
  * a CLEAN script inside the same shape must still be ALLOWED, because a
    gate that refuses every subshell is as wrong as one that allows every
    subshell.

The absolute-path and here-document rows are here because they are the two
ways the same paragraph could leak: an absolute script must be judged
regardless of where the `cd` points, and when the body IS the program (a
`bash` reading its here-document) the `cd` must propagate into that body.
"""
from __future__ import annotations

import pytest

from delfin import hpc_io_check

MOUNTS = "mgmt@o2ib:fs /pfs lustre rw 0 0\n"
ENV = {"HOME": "/home/x"}
WALK = "#!/bin/bash\ndu -sh /pfs/data6/workspace/\n"
CLEAN = "#!/bin/bash\nsrun hostname\n"


@pytest.fixture
def lustre(monkeypatch):
    """A Lustre mount in the fixture, so the rule fires on content."""
    monkeypatch.setattr(hpc_io_check.MountsSource, "read", lambda self: MOUNTS)


def _script(directory, body, name="run.sh"):
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text(body, encoding="utf-8")
    return path


def test_subshell_cd_relative_walk_script_is_refused(lustre, tmp_path):
    """The refusal must be about the CONTENT, not about a missing script.

    Asserts the reason, not just the verdict: the gate also refuses a
    relative script that does not resolve under the cd target, with "not
    found from the cd target". A bare `verdict == "refuse"` would pass even
    if the cd never reached the file and its content was never read -- a
    false green. The walk reason and the resolved file's line number are
    required, and the resolution reason is explicitly excluded.
    """
    sub = tmp_path / "sub"
    _script(sub, WALK)
    verdict, text, may_run = hpc_io_check.sbatch_gate(
        f"( cd {sub}; sbatch run.sh )", env=dict(ENV))
    assert verdict == "refuse", (verdict, text)
    assert may_run is False
    assert "tree_walk_on_lustre" in text, (
        "refused for the wrong reason -- content was not read")
    assert "not found from the cd target" not in text, (
        "the cd did not reach the file; this refusal proves nothing about "
        "content")
    # The walking line of the RESOLVED file, so the content of sub/run.sh
    # demonstrably came from the cd target and not from the process cwd.
    assert "run.sh:2" in text, text


def test_subshell_cd_missing_script_refuses_for_resolution(lustre, tmp_path):
    """The control that keeps the test above honest.

    A subshell cd naming a script that does not exist must give the
    RESOLUTION reason. If this ever carried the walk reason, the content
    assertion above would be measuring the wrong thing.
    """
    sub = tmp_path / "sub"
    sub.mkdir()
    verdict, text, _may_run = hpc_io_check.sbatch_gate(
        f"( cd {sub}; sbatch nosuch.sh )", env=dict(ENV))
    assert verdict == "refuse"
    assert "not found from the cd target" in text, text
    assert "tree_walk_on_lustre" not in text, text


def test_subshell_cd_relative_clean_script_still_runs(lustre, tmp_path):
    sub = tmp_path / "sub"
    _script(sub, CLEAN)
    verdict, text, may_run = hpc_io_check.sbatch_gate(
        f"( cd {sub}; sbatch run.sh )", env=dict(ENV))
    assert may_run is True, (verdict, text)


def test_subshell_cd_absolute_walk_script_is_refused(lustre, tmp_path):
    sub = tmp_path / "sub"
    path = _script(sub, WALK)
    verdict, _text, may_run = hpc_io_check.sbatch_gate(
        f"( cd {sub}; sbatch {path} )", env=dict(ENV))
    assert verdict == "refuse"
    assert may_run is False


def test_subshell_cd_reaches_a_here_document_body(lustre, tmp_path):
    sub = tmp_path / "sub"
    _script(sub, WALK)
    verdict, _text, may_run = hpc_io_check.sbatch_gate(
        f"( cd {sub}; bash ) <<'EOF'\nsbatch run.sh\nEOF\n", env=dict(ENV))
    assert verdict == "refuse"
    assert may_run is False
