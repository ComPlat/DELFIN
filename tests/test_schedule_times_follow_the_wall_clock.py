"""When a user schedule fires: the written forms, cron, and the wall clock.

A schedule is stored as the text the user wrote (``täglich 07:30``,
``mo-fr 18:00``, ``alle 2h``, ``30 7 * * 1-5``) and its time zone; the next
run is computed from both in that zone, so a run at 07:30 stays at 07:30
across a change between summer and winter time. Times are fixed instants
here, never the clock of the machine running the test.
"""

from __future__ import annotations

from datetime import datetime
from zoneinfo import ZoneInfo

import pytest

from delfin.agent import schedules

BERLIN = "Europe/Berlin"


def _ts(text: str, tz: str = BERLIN) -> float:
    return datetime.fromisoformat(text).replace(tzinfo=ZoneInfo(tz)).timestamp()


def _local(ts: float, tz: str = BERLIN) -> str:
    return datetime.fromtimestamp(ts, ZoneInfo(tz)).strftime("%a %Y-%m-%d %H:%M %Z")


def _next(when: str, after: str, tz: str = BERLIN) -> str:
    spec, rest = schedules.split_when(when)
    assert rest == ""
    return _local(schedules.next_fire(spec, _ts(after, tz), tz), tz)


@pytest.mark.parametrize("text, after, expected", [
    ("täglich 07:30", "2026-10-27T08:00", "Wed 2026-10-28 07:30 CET"),
    ("daily 07:30", "2026-10-27T07:00", "Tue 2026-10-27 07:30 CET"),
    ("jeden Tag um 7:30 Uhr", "2026-10-27T07:30", "Wed 2026-10-28 07:30 CET"),
    ("every day at 7:30", "2026-10-27T06:00", "Tue 2026-10-27 07:30 CET"),
    ("mo-fr 18:00", "2026-10-30T19:00", "Mon 2026-11-02 18:00 CET"),
    ("werktags 18:00", "2026-10-31T10:00", "Mon 2026-11-02 18:00 CET"),
    ("weekdays 18:00", "2026-10-28T10:00", "Wed 2026-10-28 18:00 CET"),
    ("wochenende 09:00", "2026-10-27T10:00", "Sat 2026-10-31 09:00 CET"),
    ("montags 8 uhr", "2026-10-27T10:00", "Mon 2026-11-02 08:00 CET"),
    ("mo,mi,fr 07:00", "2026-10-27T10:00", "Wed 2026-10-28 07:00 CET"),
    ("sun 22:15", "2026-10-27T10:00", "Sun 2026-11-01 22:15 CET"),
    ("täglich 07:30, 18:00", "2026-10-27T08:00", "Tue 2026-10-27 18:00 CET"),
    ("täglich 07:30 und 18:00", "2026-10-27T19:00", "Wed 2026-10-28 07:30 CET"),
    ("30 7 * * 1-5", "2026-10-31T10:00", "Mon 2026-11-02 07:30 CET"),
    ("0 */6 * * *", "2026-10-27T07:00", "Tue 2026-10-27 12:00 CET"),
    ("15 8 1 * *", "2026-10-27T07:00", "Sun 2026-11-01 08:15 CET"),
    ("0 9 * * 0", "2026-10-27T07:00", "Sun 2026-11-01 09:00 CET"),
    ("0 9 * * 7", "2026-10-27T07:00", "Sun 2026-11-01 09:00 CET"),
])
def test_the_written_forms_name_the_next_run(text, after, expected):
    assert _next(text, after) == expected


def test_cron_day_of_month_and_weekday_widen_each_other():
    # As cron does it: both restricted means either one matches.
    assert _next("0 8 13 * 5", "2026-11-01T00:00") == "Fri 2026-11-06 08:00 CET"
    assert _next("0 8 13 * 5", "2026-11-12T00:00") == "Fri 2026-11-13 08:00 CET"


def test_a_run_keeps_its_wall_clock_time_across_the_change_to_winter_time():
    # 2026-10-25 03:00 CEST -> 02:00 CET.
    first = _next("täglich 07:30", "2026-10-24T08:00")
    second = _next("täglich 07:30", "2026-10-25T08:00")
    assert first == "Sun 2026-10-25 07:30 CET"
    assert second == "Mon 2026-10-26 07:30 CET"


def test_a_time_inside_the_hour_that_repeats_runs_once():
    spec, _ = schedules.split_when("täglich 02:30")
    one = schedules.next_fire(spec, _ts("2026-10-25T00:00"), BERLIN)
    two = schedules.next_fire(spec, one, BERLIN)
    assert _local(one) == "Sun 2026-10-25 02:30 CEST"
    assert _local(two) == "Mon 2026-10-26 02:30 CET"


def test_a_time_inside_the_hour_that_is_skipped_still_runs_that_day():
    # 2027-03-28 02:00 CET -> 03:00 CEST: 02:30 does not exist.
    spec, _ = schedules.split_when("täglich 02:30")
    one = schedules.next_fire(spec, _ts("2027-03-28T00:00"), BERLIN)
    assert _local(one).startswith("Sun 2027-03-28 03:30")
    two = schedules.next_fire(spec, one, BERLIN)
    assert _local(two) == "Mon 2027-03-29 02:30 CEST"


def test_the_time_zone_is_the_schedules_own():
    assert _next("täglich 07:30", "2026-10-27T08:00", "America/New_York") \
        == "Wed 2026-10-28 07:30 EDT"
    spec, rest = schedules.split_when("täglich 07:30 UTC Bericht schreiben")
    assert rest == "Bericht schreiben" and spec["tz"] == "UTC"
    spec, rest = schedules.split_when("daily 07:30 tz=Asia/Tokyo report")
    assert rest == "report" and spec["tz"] == "Asia/Tokyo"


def test_an_interval_counts_from_the_last_run():
    spec, rest = schedules.split_when("alle 2h Jobs prüfen")
    assert rest == "Jobs prüfen"
    assert schedules.next_fire(spec, 1000.0, BERLIN) == 1000.0 + 7200
    for text, seconds in (("every 30m x", 1800), ("alle 15 min x", 900),
                          ("stündlich x", 3600), ("hourly x", 3600),
                          ("alle 3 Stunden x", 3 * 3600), ("every 1d x", 86400)):
        spec, rest = schedules.split_when(text)
        assert (spec["every"], rest) == (seconds, "x"), text


@pytest.mark.parametrize("text, words", [
    ("alle 1m x", "5 minutes"),
    ("täglich x", "time"),
    ("mo-fr", "time"),
    ("Bitte die Jobs prüfen", "when"),
    ("täglich 25:00 x", "when"),
    ("61 * * * * x", "when"),
    ("0 8 31 2 * x", "never"),
    ("täglich 07:30 tz=Mars/Base x", "time zone"),
])
def test_what_cannot_be_scheduled_is_refused_with_the_reason(text, words):
    with pytest.raises(schedules.ScheduleError) as info:
        spec, _rest = schedules.split_when(text)
        schedules.next_fire(spec, _ts("2026-10-27T08:00"), spec.get("tz") or BERLIN)
    assert words in str(info.value)


def test_the_task_after_the_time_is_kept_as_written():
    spec, rest = schedules.split_when("mo-fr 18:00: Fasse die Slurm-Jobs zusammen, kurz.")
    assert rest == "Fasse die Slurm-Jobs zusammen, kurz."
    spec, rest = schedules.split_when("täglich um 7 Uhr 3 Berichte lesen")
    assert spec["times"] == [[7, 0]] and rest == "3 Berichte lesen"


def test_the_local_zone_is_named_not_an_offset(monkeypatch):
    monkeypatch.setenv("TZ", "Europe/Berlin")
    assert schedules.local_zone() == "Europe/Berlin"
    monkeypatch.setenv("TZ", "Not/AZone")
    assert schedules.local_zone() in schedules.available_zone_names() or \
        schedules.local_zone() == "UTC"
