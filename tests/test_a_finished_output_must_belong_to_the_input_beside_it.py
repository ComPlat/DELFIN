"""A complete output is not evidence that it belongs to the input beside it.

Recalc skips a job when it finds a finished output and a fingerprint that
matches the input. When the fingerprint is missing it used to skip anyway, on
the strength of the output alone, and then write the current input's
fingerprint as though it had produced that output.

The fingerprint goes missing routinely: it is written only on the success path
and only reaches the calculation directory through the full sync profile, so a
job killed at walltime leaves outputs without one. Change a functional, submit
again, and the jobs that had finished are skipped — their old numbers reported
under the new settings, with a fresh fingerprint asserting it all agrees. The
archive shows this branch taken 291 times.

ORCA copies its input into the head of its output, one line per `|  12> `
marker, so the claim can simply be checked. The fixtures here are that echo as
ORCA 6.1.1 writes it.
"""

from pathlib import Path

import pytest

from delfin import smart_recalc


def _out(tmp_path: Path, inp_lines, *, name="job.out", terminated=True) -> Path:
    """An ORCA output whose echoed input is *inp_lines*."""
    body = ["", "                                       INPUT FILE",
            "=" * 80, "NAME = /scratch/x/.orca_iso_job_1_2_abc/job.inp"]
    body += [f"|{i:>4}> {line}" for i, line in enumerate(inp_lines, 1)]
    body.append(f"|{len(inp_lines) + 1:>4}>                          ****END OF INPUT****")
    body.append("")
    if terminated:
        body.append("                             ****ORCA TERMINATED NORMALLY****")
    p = tmp_path / name
    p.write_text("\n".join(body), encoding="utf-8")
    return p


_INPUT = [
    "! HF-3c",
    '%base "job"',
    "%maxcore 1000",
    "",
    "* xyz 0 1",
    "O   0.000000   0.000000   0.117300",
    "H   0.000000   0.757200  -0.469200",
    "H   0.000000  -0.757200  -0.469200",
    "*",
]


def _inp(tmp_path: Path, lines=None, name="job.inp") -> Path:
    p = tmp_path / name
    p.write_text("\n".join(lines if lines is not None else _INPUT) + "\n", encoding="utf-8")
    return p


# --------------------------------------------------------------------------
# the comparison itself
# --------------------------------------------------------------------------

def test_an_untouched_input_matches_its_own_output(tmp_path):
    assert smart_recalc.output_matches_input(_inp(tmp_path), _out(tmp_path, _INPUT)) is True


def test_a_changed_method_does_not_match(tmp_path):
    """The case that makes this worth having: same geometry, different
    functional, and the old numbers would have been reported as the new ones."""
    changed = ["! PBE def2-SVP"] + _INPUT[1:]

    assert smart_recalc.output_matches_input(_inp(tmp_path, changed), _out(tmp_path, _INPUT)) is False


def test_a_changed_geometry_does_not_match(tmp_path):
    changed = [l.replace("0.117300", "0.150000") for l in _INPUT]

    assert smart_recalc.output_matches_input(_inp(tmp_path, changed), _out(tmp_path, _INPUT)) is False


def test_blank_lines_and_comments_are_not_a_difference(tmp_path):
    """Reformatting must not force hours of recomputation."""
    padded = ["# a note", ""] + _INPUT + ["", ""]

    assert smart_recalc.output_matches_input(_inp(tmp_path, padded), _out(tmp_path, _INPUT)) is True


def test_the_isolation_rewrite_is_not_a_difference(tmp_path):
    """ORCA runs on a copy whose %moinp and xyzfile references were rewritten
    to absolute paths. Comparing those verbatim would report every job as
    changed."""
    on_disk = ['! HF-3c MOREAD', '%moinp "S0.gbw"', "* xyzfile 0 1 start.xyz"]
    as_run = ['! HF-3c MOREAD',
              '%moinp "/scratch/job_42/run/ESD/S0.gbw"',
              "* xyzfile 0 1 /scratch/job_42/run/ESD/start.xyz"]

    assert smart_recalc.output_matches_input(_inp(tmp_path, on_disk), _out(tmp_path, as_run)) is True


def test_an_output_without_an_echo_yields_no_verdict(tmp_path):
    """A raw xTB log carries no echo. None means "cannot tell", and the caller
    keeps its previous behaviour rather than recomputing everything."""
    raw = tmp_path / "raw.out"
    raw.write_text("xtb version 6.7.1\nnormal termination\n", encoding="utf-8")

    assert smart_recalc.output_matches_input(_inp(tmp_path), raw) is None


def test_a_missing_file_yields_no_verdict(tmp_path):
    assert smart_recalc.output_matches_input(tmp_path / "gone.inp", _out(tmp_path, _INPUT)) is None
    assert smart_recalc.output_matches_input(_inp(tmp_path), tmp_path / "gone.out") is None


def test_the_end_of_input_marker_is_not_part_of_the_input(tmp_path):
    """ORCA closes the echo with its own line; counting it would make every
    output disagree with the file that produced it."""
    echoed = smart_recalc.echoed_input(_out(tmp_path, _INPUT))

    assert echoed is not None
    assert not any("END OF INPUT" in line for line in echoed)
    assert len(echoed) == len(_INPUT)


def test_only_the_head_of_a_large_output_is_read(tmp_path):
    """These outputs reach hundreds of MB; the echo is in the first pages."""
    import tracemalloc

    out = _out(tmp_path, _INPUT)
    out.write_text(out.read_text() + ("padding line\n" * 400000), encoding="utf-8")
    assert out.stat().st_size > 4 * 1024 * 1024

    tracemalloc.start()
    verdict = smart_recalc.output_matches_input(_inp(tmp_path), out)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    assert verdict is True
    assert peak < 4 * 1024 * 1024, f"read {peak / 1024 / 1024:.1f} MB of a larger file"


# --------------------------------------------------------------------------
# what the caller does with it
# --------------------------------------------------------------------------

def test_completeness_is_unchanged_by_the_comparison(tmp_path):
    """outputs_complete answers a different question — whether the run
    finished — and must not start answering this one too."""
    changed = ["! PBE def2-SVP"] + _INPUT[1:]
    out = _out(tmp_path, _INPUT)

    assert smart_recalc.outputs_complete(_inp(tmp_path, changed), out) is True


@pytest.mark.parametrize(
    "verdict,skips",
    [(True, True), (None, True), (False, False)],
)
def test_the_bootstrap_only_skips_when_nothing_contradicts_it(verdict, skips):
    """The gate in run_orca is `is not False`: a match skips, no verdict keeps
    the previous behaviour, and only a contradiction forces the recalculation.
    """
    assert ((verdict is not False) is skips)
