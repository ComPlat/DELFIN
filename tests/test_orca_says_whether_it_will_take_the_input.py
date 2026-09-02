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
        "* O R C A *\nSCF SETTINGS\n", None, still_running=True, waited=0.3)

    assert ok is True and "still running" in headline
    assert "0.3 s" in headline, (
        "how long it actually took, not how long it was allowed to"
    )


def test_the_answer_is_given_as_soon_as_orca_has_given_it():
    """Waiting out the whole window to repeat what is already on the screen is
    a minute nobody has a reason to spend.  ORCA names an unrecognised keyword
    in the first second and reaches the integrals in a few more."""
    assert builder._orca_has_spoken("... SCF SETTINGS ...") is True
    assert builder._orca_has_spoken(
        "UNRECOGNIZED OR DUPLICATED KEYWORD(S)") is True
    assert builder._orca_has_spoken("* O R C A *\nreading input\n") is False


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


def _use(monkeypatch, orca):
    """Put *orca* where the check looks, which is where a real run looks.

    Both resolvers: the check asks DELFIN's own first and falls back to the
    editor's, so a stand-in has to stand in both places or the real ORCA on
    the machine answers instead.
    """
    monkeypatch.setattr("delfin.orca.find_orca_executable", lambda: orca)
    monkeypatch.setattr("delfin.dashboard.saddle.find_orca", lambda: orca)


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
    _use(monkeypatch, str(fake))

    said = _pressed(builder_tab, capsys)

    assert "STOPPED" in said
    assert "FOOBARNONSENSE" in said


def test_the_button_says_so_when_orca_gets_going(builder_tab, monkeypatch,
                                                capsys):
    fake = _an_orca(builder_tab["tmp_path"], "orca_ok",
                    "echo 'ORCA GTO INTEGRAL CALCULATION'\nexit 0\n")
    _use(monkeypatch, str(fake))

    said = _pressed(builder_tab, capsys)

    assert "OK" in said
    assert "PAL 1" in said, "and at what it was checked"


def test_one_that_keeps_running_is_stopped_and_counted_as_started(
        builder_tab, monkeypatch, capsys):
    """A check that waited for the calculation would be the calculation."""
    fake = _an_orca(builder_tab["tmp_path"], "orca_slow",
                    "echo 'SCF SETTINGS'\nsleep 300\n")
    _use(monkeypatch, str(fake))
    monkeypatch.setattr(builder, "CHECK_SECONDS", 2.0)

    began = time.monotonic()
    said = _pressed(builder_tab, capsys)
    took = time.monotonic() - began

    assert "OK" in said and "still running" in said
    assert took < 45, f"the check waited {took:.0f}s on a job that never ends"


def test_a_good_input_is_answered_in_seconds_not_at_the_window(
        builder_tab, monkeypatch, capsys):
    """The window is the fallback, not the wait.  It said "Starting ORCA..."
    for a minute before this, on an input that was fine."""
    fake = _an_orca(builder_tab["tmp_path"], "orca_quick",
                    "echo 'SCF SETTINGS'\nsleep 300\n")
    _use(monkeypatch, str(fake))
    monkeypatch.setattr(builder, "CHECK_SECONDS", 45.0)

    began = time.monotonic()
    said = _pressed(builder_tab, capsys)
    took = time.monotonic() - began

    assert "OK" in said
    assert took < 15, (
        f"answered after {took:.0f}s, though ORCA had said so at once"
    )


def test_no_orca_to_check_with_is_said_rather_than_guessed(builder_tab,
                                                           monkeypatch, capsys):
    _use(monkeypatch, None)

    said = _pressed(builder_tab, capsys)

    assert "No ORCA" in said


def test_the_check_writes_nothing_into_the_calculation_directory(builder_tab,
                                                                 monkeypatch,
                                                                 capsys):
    """A check is not a job.  SUBMIT saves; this does not."""
    fake = _an_orca(builder_tab["tmp_path"], "orca_ok2",
                    "echo 'SCF SETTINGS'\nexit 0\n")
    _use(monkeypatch, str(fake))

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


# ---------------------------------------------------------------------------
# a check that cannot finish is not a check
# ---------------------------------------------------------------------------
def test_orca_is_given_no_standing_input_to_wait_on():
    """A wrapper that reads stdin waits for something nobody is going to type,
    and the dashboard sits on "Starting ORCA..." for ever."""
    source = pathlib.Path(builder.__file__).read_text(encoding="utf-8")
    started = source.split("running = subprocess.Popen(")[1].split("\n        )")[0]

    assert "stdin=subprocess.DEVNULL" in started


def test_one_that_says_nothing_and_waits_is_still_answered(builder_tab,
                                                           monkeypatch, capsys):
    """Silent and blocked on stdin is the shape that hung.  It has to come
    back with an answer like everything else."""
    fake = _an_orca(builder_tab["tmp_path"], "orca_mute",
                    "cat > /dev/null\nsleep 300\n")
    _use(monkeypatch, str(fake))
    monkeypatch.setattr(builder, "CHECK_SECONDS", 2.0)

    began = time.monotonic()
    said = _pressed(builder_tab, capsys)
    took = time.monotonic() - began

    assert "OK" in said or "STOPPED" in said, "an answer either way"
    assert took < 30, f"still waiting after {took:.0f}s"


def test_it_says_which_orca_it_found_before_it_waits(builder_tab, monkeypatch,
                                                     capsys):
    """So that a wait is never a blank one, and so that a report of it hanging
    says where it was standing."""
    fake = _an_orca(builder_tab["tmp_path"], "orca_named",
                    "echo 'SCF SETTINGS'\nexit 0\n")
    _use(monkeypatch, str(fake))

    said = _pressed(builder_tab, capsys)

    assert "Looking for ORCA" in said
    assert "orca_named" in said, "the one it found, by name"


def test_no_answer_at_all_is_said_rather_than_left_standing():
    """A button that stays greyed with one line under it is not an answer."""
    source = pathlib.Path(builder.__file__).read_text(encoding="utf-8")
    guard = source.split("def nothing_came_back")[1].split("\n    watchdog")[0]

    assert "orca_check_btn.disabled = False" in guard, "the button comes back"
    assert "No answer after" in guard
    assert "watchdog.cancel()" in source, "and it is called off when one comes"


# ---------------------------------------------------------------------------
# started the way a real run starts it
# ---------------------------------------------------------------------------
def test_it_starts_orca_the_way_delfin_starts_orca():
    """A check that started ORCA its own way would be checking its own way of
    starting ORCA.  The scratch directory, the library path and the MPI
    settings are all in that environment, and without them ORCA can fail for a
    reason that has nothing to do with the input."""
    source = pathlib.Path(builder.__file__).read_text(encoding="utf-8")

    assert "find_orca_executable" in source, "the resolver a real run uses"
    assert "_prepare_orca_environment" in source
    started = source.split("running = subprocess.Popen(")[1].split("\n        )")[0]
    assert "env=environment" in started
    assert "start_new_session=True" in started, (
        "its own process group, the way a real run gets one"
    )


def test_stopping_it_takes_the_children_with_it():
    """ORCA starts several, and one left behind holds the pipe open -- a read
    on it never returns, which is what "Starting ORCA..." and then nothing
    was."""
    source = pathlib.Path(builder.__file__).read_text(encoding="utf-8")
    stopping = source.split("def _stop(running)")[1].split("\n    def ")[0]

    assert "killpg" in stopping
    assert "SIGKILL" in stopping, "and it does not ask twice for ever"
