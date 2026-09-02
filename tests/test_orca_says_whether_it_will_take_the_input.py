"""An input that ORCA will not start on is worth knowing before the queue.

A job submitted with a keyword ORCA does not recognise waits for a node, takes
it, and stops on the first line it reads.  That is a round trip through the
scheduler to be told something the input said all along.

So the input is offered to ORCA here, on one core and a small memory ceiling,
and the only question asked is whether it gets going.  What a calculation does
after that is what a real run is for -- an SCF that will not converge in hour
three cannot be found this way and is not claimed to be.
"""

from __future__ import annotations

import pathlib
import time

import pytest

from delfin.dashboard import tab_orca_builder as builder
from delfin.dashboard.context import DashboardContext


# ---------------------------------------------------------------------------
# cut down to something a login node can run
# ---------------------------------------------------------------------------
def test_the_check_runs_on_one_core_and_a_small_ceiling():
    """The point is to run it where the dashboard is.  Sixteen cores and six
    gigabytes is not that."""
    text = builder.input_for_check(
        "! B3LYP def2-SVP\n%pal\n  nprocs 16\nend\n%maxcore 6000\n"
        "* xyz 0 1\nH 0 0 0\n*\n")

    assert "nprocs 1\n" in text
    assert "%maxcore 1000" in text
    assert "nprocs 16" not in text and "6000" not in text


def test_the_keyword_form_of_pal_is_taken_out_too():
    """``! PAL8`` sends it just as wide as the block does."""
    text = builder.input_for_check("! B3LYP PAL8 def2-SVP\n* xyz 0 1\nH 0 0 0\n*\n")

    assert "PAL8" not in text
    assert "B3LYP" in text and "def2-SVP" in text


def test_nothing_else_about_the_input_is_touched():
    """A check that quietly rewrote the input would be checking a different
    input."""
    body = ("! PBE0 OPT def2-TZVP D4\n%scf\n  maxiter 300\nend\n"
            "* xyzfile 0 1 thing.xyz\n")

    text = builder.input_for_check(body)

    assert "PBE0 OPT def2-TZVP D4" in text
    assert "maxiter 300" in text
    assert "* xyzfile 0 1 thing.xyz" in text


# ---------------------------------------------------------------------------
# reading what ORCA said
# ---------------------------------------------------------------------------
def test_an_unrecognised_keyword_is_reported_with_the_keyword_in_it():
    """Which is the whole use of it.  Measured against ORCA 6.1.1, the lines
    that matter are the banner, the sentence and the keyword under it."""
    said = (
        "                     * O   R   C   A *\n"
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
        "                INPUT ERROR\n"
        "  UNRECOGNIZED OR DUPLICATED KEYWORD(S) IN SIMPLE INPUT LINE\n"
        "    FOOBARNONSENSE\n"
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    )

    ok, headline, detail = builder.orca_startup_report(said, 1)

    assert ok is False
    assert "will not take" in headline
    assert "FOOBARNONSENSE" in detail, "the offending keyword is what to show"


def test_getting_into_the_calculation_is_the_answer_yes():
    ok, headline, _ = builder.orca_startup_report(
        "* O R C A *\nORCA GTO INTEGRAL CALCULATION\nSCF SETTINGS\n", 0)

    assert ok is True and "read the input" in headline


def test_still_running_is_also_the_answer_yes():
    """It got past the input, which is the only thing being asked."""
    ok, headline, _ = builder.orca_startup_report(
        "* O R C A *\nSCF SETTINGS\n", None, still_running=True)

    assert ok is True and "still running" in headline


def test_stopping_without_saying_why_is_still_reported_as_stopping():
    ok, headline, detail = builder.orca_startup_report(
        "* O R C A *\nsomething went wrong somewhere\n", 127)

    assert ok is False
    assert "127" in headline
    assert "something went wrong" in detail, "the last of what it said"


# ---------------------------------------------------------------------------
# the button, over a stand-in ORCA
# ---------------------------------------------------------------------------
def _an_orca(room, name, body):
    """A stand-in that behaves the way ORCA does in one particular way."""
    where = room / name
    where.write_text("#!/bin/sh\n" + body)
    where.chmod(0o755)
    return where


@pytest.fixture
def builder_tab(tmp_path, monkeypatch):
    pytest.importorskip("ipywidgets")

    (tmp_path / "calc").mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "calc",
                           office_dir=tmp_path / "calc")
    ctx.run_js = lambda _script: None
    _widget, refs = builder.create_tab(ctx)
    refs["orca_coords"].value = (
        "3\nwater\nO 0.0 0.0 0.0\nH 0.0 0.0 0.96\nH 0.93 0.0 -0.24\n")
    refs["update_orca_preview"]()
    refs["tmp_path"] = tmp_path
    return refs


def _pressed(refs, capsys, seconds=60):
    """Press it and wait for the answer.

    Read off stdout rather than out of the Output widget: without a kernel
    around it, an Output captures nothing and the text goes where print always
    sends it.
    """
    refs["orca_check_btn"].click()
    for _ in range(int(seconds * 4)):
        time.sleep(0.25)
        if not refs["orca_check_btn"].disabled:
            break
    return capsys.readouterr().out


def test_the_button_says_what_orca_said_about_the_keyword(builder_tab,
                                                          monkeypatch, capsys):
    fake = _an_orca(builder_tab["tmp_path"], "orca_bad", """
cat <<'EOT'
                INPUT ERROR
  UNRECOGNIZED OR DUPLICATED KEYWORD(S) IN SIMPLE INPUT LINE
    FOOBARNONSENSE
EOT
exit 1
""")
    monkeypatch.setattr("delfin.dashboard.saddle.find_orca", lambda: str(fake))

    said = _pressed(builder_tab, capsys)

    assert "STOPPED" in said
    assert "FOOBARNONSENSE" in said


def test_the_button_says_so_when_orca_gets_going(builder_tab, monkeypatch,
                                                capsys):
    fake = _an_orca(builder_tab["tmp_path"], "orca_ok",
                    "echo 'ORCA GTO INTEGRAL CALCULATION'\nexit 0\n")
    monkeypatch.setattr("delfin.dashboard.saddle.find_orca", lambda: str(fake))

    said = _pressed(builder_tab, capsys)

    assert "OK" in said
    assert "PAL 1" in said, "and at what it was checked"


def test_one_that_keeps_running_is_stopped_and_counted_as_started(
        builder_tab, monkeypatch, capsys):
    """A check that waited for the calculation would be the calculation."""
    fake = _an_orca(builder_tab["tmp_path"], "orca_slow",
                    "echo 'SCF SETTINGS'\nsleep 300\n")
    monkeypatch.setattr("delfin.dashboard.saddle.find_orca", lambda: str(fake))
    monkeypatch.setattr(builder, "CHECK_SECONDS", 2.0)

    began = time.monotonic()
    said = _pressed(builder_tab, capsys)
    took = time.monotonic() - began

    assert "OK" in said and "still running" in said
    assert took < 45, f"the check waited {took:.0f}s on a job that never ends"


def test_no_orca_to_check_with_is_said_rather_than_guessed(builder_tab,
                                                           monkeypatch, capsys):
    monkeypatch.setattr("delfin.dashboard.saddle.find_orca", lambda: None)

    said = _pressed(builder_tab, capsys)

    assert "No ORCA" in said


def test_the_check_writes_nothing_into_the_calculation_directory(builder_tab,
                                                                 monkeypatch,
                                                                 capsys):
    """A check is not a job.  SUBMIT saves; this does not."""
    fake = _an_orca(builder_tab["tmp_path"], "orca_ok2",
                    "echo 'SCF SETTINGS'\nexit 0\n")
    monkeypatch.setattr("delfin.dashboard.saddle.find_orca", lambda: str(fake))

    _pressed(builder_tab, capsys)

    assert list((builder_tab["tmp_path"] / "calc").iterdir()) == []


# ---------------------------------------------------------------------------
# where it stands
# ---------------------------------------------------------------------------
def test_it_stands_before_the_press_it_is_there_to_come_before():
    """The order on screen is the order of the work: ask whether ORCA will
    take it, then send it to the queue."""
    source = pathlib.Path(builder.__file__).read_text(encoding="utf-8")
    row = source.split("orca_save_btn, orca_slurm_time")[1].split("]")[0]

    assert row.index("orca_check_btn") < row.index("orca_submit_btn")
