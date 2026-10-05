"""Adversarial review, finding F4 (newline/CR/NUL/tab smuggling into the job script).

slurm_tests.build_job passes repo, ref, partition, mem, job_name and output_dir
through shlex.quote, which does NOT neutralise newlines/CR/NUL. A payload like
    ref='main\\ntouch /tmp/X'
becomes, after shlex.quote, the single-quoted string 'main\\ntouch /tmp/X' --
bash single quotes keep the newline, so 'touch /tmp/X' lands as its OWN command
line in the rendered script (code execution on the compute node). Likewise
    partition='cpu\\n#SBATCH --export=ALL,EVIL=1'
injects the forbidden '--export=ALL,VAR' SBATCH form (the exact thing the
cluster-incident rule forbids).

These tests encode the fix the builder must land: strict allow-lists BEFORE
rendering, each parameter refused with ValueError when it carries a control
character / bad form, and an invariant that no rendered line beyond the fixed
#SBATCH set starts with '#SBATCH'.

ALL of these must be RED on the current (unfixed) code: today the bad values
are accepted and rendered into the script, so the ValueError assertions fail.
Reviewer: nacht-s24, package T6.
"""
import pytest

from delfin.agent.slurm_tests import build_job

#: the fixed #SBATCH lines build_job may emit (only these, and only their
#: exact prefix). Any other line starting with '#SBATCH' is an injection.
_FIXED_SBATCH_STARTS = (
    "#SBATCH --job-name=",
    "#SBATCH --nodes=",
    "#SBATCH --ntasks=",
    "#SBATCH --cpus-per-task=",
    "#SBATCH --mem=",
    "#SBATCH --time=",
    "#SBATCH --partition=",
    "#SBATCH --export=ALL",
    "#SBATCH --output=",
    "#SBATCH --error=",
)


def _good() -> dict:
    return dict(
        repo="delfin",
        ref="main",
        test_paths=["tests/a.py"],
        partition="cpu",
        minutes=10,
        mem="8G",
        job_name="delfin-agent-tests",
        output_dir="logs",
    )


def _payloads():
    """(name, payload) control characters a caller could smuggle in."""
    return [
        ("newline", "\n"),
        ("CR", "\r"),
        ("NUL", "\x00"),
        ("tab", "\t"),
    ]


@pytest.mark.parametrize("name,ctl", _payloads())
@pytest.mark.parametrize(
    "param",
    ["ref", "partition", "mem", "job_name", "output_dir", "repo"],
)
def test_control_char_refused(param, ctl, name):
    """Every injectable parameter must refuse a control character (ValueError)."""
    good = _good()
    good[param] = good[param] + ctl + "bogus"
    with pytest.raises(ValueError):
        build_job(**good)


@pytest.mark.parametrize("param", ["ref", "partition", "mem", "job_name", "output_dir", "repo"])
def test_leading_dash_refused(param):
    """A leading '-' must be refused (option-swallowing on the #SBATCH line)."""
    good = _good()
    good[param] = "-watch"
    with pytest.raises(ValueError):
        build_job(**good)


def test_ref_newline_cannot_execute_on_its_own_line():
    """A newline in ref must not reach the script: build_job refuses it (ValueError)."""
    good = _good()
    good["ref"] = "main\ntouch MARKER_NEWLINE"
    with pytest.raises(ValueError):
        build_job(**good)


def test_partition_cannot_inject_sbatch_line():
    """A newline in partition must not smuggle a '#SBATCH --export=ALL,VAR' line: refused."""
    good = _good()
    good["partition"] = "cpu\n#SBATCH --export=ALL,EVIL=1"
    with pytest.raises(ValueError):
        build_job(**good)


def test_no_stray_sbatch_lines():
    """A control char in any #SBATCH-valued param is refused (no stray directive)."""
    for ctl in ("\n", "\r", "\t", "\x00"):
        good = _good()
        good["partition"] = "cpu" + ctl
        with pytest.raises(ValueError):
            build_job(**good)


def test_good_job_still_renders_the_fixed_header():
    """The happy path still carries exactly the fixed #SBATCH header."""
    script = build_job(**_good())
    body = script.splitlines()
    evil = [L for L in body if L.startswith("#SBATCH") and not any(L.startswith(p) for p in _FIXED_SBATCH_STARTS)]
    assert not evil
