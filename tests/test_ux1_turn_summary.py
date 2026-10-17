"""UX1 phase-2: the turn summary line, told like Claude Code tells it.

``format_turn_summary`` is pure: no I/O, no clock, no registry. The
caller (CLI footer, dashboard hook) supplies duration, end time,
outcome and the count of background jobs still running.

The line the person sees:

    Surfaced after 31s · done 11:11 · 1 shell still running
"""

import datetime

import pytest

from delfin.agent.turn_summary import format_turn_summary


# --- duration formatting -------------------------------------------------

def test_under_a_minute_is_plain_seconds():
    line = format_turn_summary(duration_s=31, ended_at=None, outcome="ok")
    assert line.startswith("Surfaced after 31s")


def test_seconds_are_rounded_not_floored_with_decimals():
    line = format_turn_summary(duration_s=31.6, ended_at=None, outcome="ok")
    assert line.startswith("Surfaced after 32s")


def test_over_a_minute_is_m_and_s():
    line = format_turn_summary(duration_s=133, ended_at=None, outcome="ok")
    assert line.startswith("Surfaced after 2m 13s")


def test_exactly_one_minute_has_no_zero_seconds_padding():
    line = format_turn_summary(duration_s=60.0, ended_at=None, outcome="ok")
    assert line.startswith("Surfaced after 1m 0s")


def test_over_an_hour_is_h_and_m():
    line = format_turn_summary(duration_s=3840, ended_at=None, outcome="ok")
    assert line.startswith("Surfaced after 1h 4m")


def test_missing_duration_omits_the_part():
    line = format_turn_summary(duration_s=None, ended_at=None, outcome="ok")
    assert "after" not in line


def test_negative_duration_is_treated_as_zero():
    line = format_turn_summary(duration_s=-5, ended_at=None, outcome="ok")
    assert line.startswith("Surfaced after 0s")


# --- the verb ------------------------------------------------------------

def test_default_verb_is_surfaced():
    line = format_turn_summary(duration_s=31, ended_at=None, outcome="ok")
    assert line.startswith("Surfaced after 31s")


def test_rotate_mode_rotates_the_verb():
    # "rotate" is caller-driven: the caller passes ``pick`` (e.g. its
    # turn counter) so this stays a pure function. Five picks must
    # cover the five verbs, six draws must not repeat all the time.
    verbs = [
        format_turn_summary(duration_s=31, ended_at=None, outcome="ok",
                            verb_mode="rotate", pick=n)
        for n in range(6)
    ]
    first_words = [line.split()[0] for line in verbs]
    assert len(set(first_words)) == 5
    assert first_words[0] == first_words[5]  # wraps around


def test_stopped_turn_never_gets_a_success_verb():
    line = format_turn_summary(duration_s=12, ended_at=None, outcome="stopped")
    assert line.startswith("Stopped after 12s")
    assert "Surfaced" not in line


def test_failed_turn_never_gets_a_success_verb():
    line = format_turn_summary(duration_s=12, ended_at=None, outcome="error")
    assert line.startswith("Stopped after 12s")
    assert "Surfaced" not in line


def test_unknown_outcome_is_treated_as_a_stopped_turn():
    line = format_turn_summary(duration_s=12, ended_at=None, outcome="weird")
    assert line.startswith("Stopped after 12s")


def test_rotate_mode_still_never_succeeds_a_failed_turn():
    line = format_turn_summary(duration_s=12, ended_at=None, outcome="error",
                               verb_mode="rotate")
    assert line.startswith("Stopped after 12s")


# --- done time -----------------------------------------------------------

def test_done_time_in_local_hhmm():
    ended = datetime.datetime(2026, 10, 18, 11, 11)
    line = format_turn_summary(duration_s=31, ended_at=ended, outcome="ok")
    assert " · done 11:11" in line


def test_done_time_is_omitted_when_not_given():
    line = format_turn_summary(duration_s=31, ended_at=None, outcome="ok")
    assert "done" not in line


# --- background jobs -----------------------------------------------------

def test_one_running_job_is_counted():
    line = format_turn_summary(duration_s=31, ended_at=None, outcome="ok",
                               running_jobs=1)
    assert " · 1 shell still running" in line


def test_two_running_jobs_pluralise():
    line = format_turn_summary(duration_s=31, ended_at=None, outcome="ok",
                               running_jobs=2)
    assert " · 2 jobs still running" in line


def test_zero_running_jobs_omits_the_part():
    line = format_turn_summary(duration_s=31, ended_at=None, outcome="ok",
                               running_jobs=0)
    assert "running" not in line
    assert line == "Surfaced after 31s"


def test_missing_job_count_omits_the_part():
    line = format_turn_summary(duration_s=31, ended_at=None, outcome="ok")
    assert "running" not in line


# --- the full example line -----------------------------------------------

def test_the_full_line_matches_the_target_shape():
    ended = datetime.datetime(2026, 10, 18, 11, 11)
    line = format_turn_summary(duration_s=31, ended_at=ended, outcome="ok",
                               running_jobs=1)
    assert line == "Surfaced after 31s · done 11:11 · 1 shell still running"


# --- hostile input -------------------------------------------------------

def test_job_name_like_text_cannot_inject_markup():
    # The function only ever counts jobs; it must not echo a job label.
    line = format_turn_summary(
        duration_s=31, ended_at=None, outcome="ok",
        running_jobs="<img src=x onerror=alert(1)>")
    assert "<img" not in line
    assert "onerror" not in line


def test_job_name_string_is_treated_as_a_count_of_one():
    # A caller passing the label by mistake must not render a fake count.
    line = format_turn_summary(duration_s=31, ended_at=None, outcome="ok",
                               running_jobs="abc")
    assert " · 1 shell still running" in line


def test_hostile_outcome_cannot_inject_text():
    line = format_turn_summary(duration_s=31, ended_at=None,
                               outcome="<script>alert(1)</script>")
    assert "<script" not in line
    assert line.startswith("Stopped after 31s")


def test_no_html_ever_in_the_output():
    line = format_turn_summary(duration_s=31, ended_at=None, outcome="ok",
                               running_jobs=3)
    assert "<" not in line and ">" not in line and "&" not in line
