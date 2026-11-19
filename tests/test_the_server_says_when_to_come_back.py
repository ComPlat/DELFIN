"""A Retry-After is the one number the service actually knows.

The transient-error ladder (1.5, 3, 6, 30, 60, 120, 240 s with jitter) is
a guess made in this repository. It rides out a proxy outage well, which
is what it was measured for. But when the server itself answers with
``Retry-After``, guessing is the wrong thing to do in both directions: a
guess that is too long wastes the person's turn, and a guess that is too
short adds load to a service that has just said it is struggling.

That matters more now than it did. Blablador is a shared research service
at Forschungszentrum Jülich, and KIT's gateway is shared across the
institute; neither is a dedicated endpoint whose capacity is ours.

The header was read nowhere before this. A quota 429 is a different and
longer decision and keeps its own path (``_quota_429_reset``), which
parks the run until the hourly reset rather than retrying into it.
"""

from __future__ import annotations

import datetime as dt

import pytest

from delfin.agent import api_client as ac


class _Headers(dict):
    """Just enough of a header mapping; real ones are case-insensitive."""


class _Response:
    def __init__(self, headers):
        self.headers = headers


class _Error(Exception):
    def __init__(self, headers=None):
        super().__init__("too many requests")
        if headers is not None:
            self.response = _Response(_Headers(headers))


def test_a_plain_number_of_seconds_is_honoured():
    assert ac._retry_after_seconds(_Error({"retry-after": "12"})) == 12.0


def test_the_header_is_found_whichever_way_it_is_spelled():
    assert ac._retry_after_seconds(_Error({"Retry-After": "7"})) == 7.0


def test_an_http_date_is_resolved_against_now():
    now = 1_700_000_000.0
    when = dt.datetime.fromtimestamp(now + 45, dt.timezone.utc)
    stamp = when.strftime("%a, %d %b %Y %H:%M:%S GMT")

    seconds = ac._retry_after_seconds(_Error({"retry-after": stamp}), now=now)

    assert seconds is not None and 43 <= seconds <= 47, seconds


def test_a_date_in_the_past_waits_the_floor_not_a_negative_time():
    now = 1_700_000_000.0
    when = dt.datetime.fromtimestamp(now - 600, dt.timezone.utc)
    stamp = when.strftime("%a, %d %b %Y %H:%M:%S GMT")

    assert ac._retry_after_seconds(_Error({"retry-after": stamp}),
                                   now=now) == ac._RETRY_AFTER_MIN_S


def test_zero_means_soon_not_instantly():
    """Honouring a literal 0 would retry inside the same millisecond and
    turn a rate limit into a tight loop."""
    assert ac._retry_after_seconds(_Error({"retry-after": "0"})) == \
        ac._RETRY_AFTER_MIN_S


def test_an_absurd_wait_is_refused_so_the_ladder_decides():
    """A turn that sleeps for half an hour on a header nobody checked is
    indistinguishable from a hang. Beyond the cap we do not trust it."""
    assert ac._retry_after_seconds(
        _Error({"retry-after": str(int(ac._RETRY_AFTER_MAX_S) + 60)})) is None


def test_the_cap_itself_is_still_honoured():
    assert ac._retry_after_seconds(
        _Error({"retry-after": str(ac._RETRY_AFTER_MAX_S)})) == \
        ac._RETRY_AFTER_MAX_S


@pytest.mark.parametrize("headers", [
    None,                                   # no response at all
    {},                                     # a response, no header
    {"retry-after": ""},                    # present and empty
    {"retry-after": "soon"},                # not a number, not a date
    {"retry-after": "nan"},                 # a float that is not one
])
def test_anything_unusable_answers_none_rather_than_guessing(headers):
    assert ac._retry_after_seconds(_Error(headers)) is None


def test_an_exception_with_no_response_does_not_raise():
    """Most transient errors are timeouts and connection resets that
    carry no response at all. The reader must never be the thing that
    breaks a retry."""
    assert ac._retry_after_seconds(Exception("connection reset")) is None
    assert ac._retry_after_seconds(None) is None


def test_the_retry_path_prefers_the_header_over_the_ladder():
    """The wiring: the ladder is the fallback, not the first choice."""
    import inspect

    src = inspect.getsource(ac)
    i = src.index("_delay = _retry_after_seconds(")
    block = src[i:i + 260]

    assert "_transient_retry_delay" in block, (
        "the ladder must still be there for when the header is absent")
    assert block.index("_retry_after_seconds") < block.index(
        "_transient_retry_delay"), "the header is asked first"


def test_the_quota_reset_keeps_its_own_longer_path():
    """A quota 429 says the key is spent until the hourly reset. That is a
    different decision from "come back in N seconds" and must not be
    shortened into one."""
    assert hasattr(ac, "_quota_429_reset")
    assert ac._QUOTA_RESET_MAX_S > ac._RETRY_AFTER_MAX_S, (
        "a quota reset may legitimately park a run far longer than any "
        "Retry-After we would accept")


def _real_sdk_error(status: int, headers: dict | None):
    """A genuine SDK exception, not a stand-in with a dict for headers.

    The tests above use a plain dict, which is a `get`-shaped object and
    nothing more. What a retry actually receives is an openai exception
    whose `.response.headers` is the transport's own case-insensitive
    mapping -- a different object with different lookup rules. A reader
    that works on a dict and not on that one would pass every test above
    and read nothing in production.
    """
    openai = pytest.importorskip("openai")
    try:
        import httpx2 as transport          # what openai 3.x ships with
    except ImportError:                      # pragma: no cover - older SDK
        transport = pytest.importorskip("httpx")
    request = transport.Request(
        "POST", "https://api.blablador.fz-juelich.de/v1/chat/completions")
    response = transport.Response(status, request=request,
                                  **({"headers": headers} if headers else {}))
    kind = (openai.RateLimitError if status == 429
            else openai.APIStatusError)
    return kind("rate limited", response=response, body=None)


@pytest.mark.parametrize("spelling", ["Retry-After", "retry-after",
                                      "RETRY-AFTER"])
def test_a_real_sdk_error_is_read_whatever_the_case(spelling):
    exc = _real_sdk_error(429, {spelling: "23"})

    assert ac._retry_after_seconds(exc) == 23.0
    assert ac._is_transient_api_error(exc), (
        "a 429 the reader understands must also be one the retry path "
        "takes; reading the header on an error nobody retries is useless")


def test_a_real_sdk_error_without_the_header_falls_back_to_the_ladder():
    assert ac._retry_after_seconds(_real_sdk_error(429, None)) is None


def test_a_real_sdk_error_with_a_far_future_date_is_not_trusted():
    exc = _real_sdk_error(503, {"retry-after": "Wed, 21 Oct 2099 07:28:00 GMT"})

    assert ac._retry_after_seconds(exc) is None


# ---------------------------------------------------------------------------
# A wait longer than we will sleep must END the attempt, not shorten it
# ---------------------------------------------------------------------------
#
# Measured from the operator's own documentation (AI@JSC, Blablador API
# access, read 2026-11-19): exceeding a model's per-minute request limit
# returns 429 with Retry-After, blocks that model for 60 s, and EVERY
# further violation doubles the block up to 3600 s. External callers get
# 60 requests/min on alias-fast, 40 on alias-code and alias-large, 30 on
# alias-huge, and 20 on any other model name.
#
# The ladder's longest step is 240 s. So before this, a stated wait above
# the 300 s ceiling fell back to the ladder and the round asked again
# inside the block -- a fresh violation each time, which doubles it. Four
# retries turn a 20-minute block into an hour. The ladder is built for a
# proxy hiccup, where retrying costs nothing; against a counter that
# punishes retries it is the wrong instrument.


def test_a_wait_we_will_sit_out_is_not_treated_as_too_long():
    assert ac._retry_after_outlasts_patience(_Error({"retry-after": "30"})) is None


def test_the_ceiling_itself_is_still_a_wait_we_sit_out():
    assert ac._retry_after_outlasts_patience(
        _Error({"retry-after": str(ac._RETRY_AFTER_MAX_S)})) is None


def test_a_wait_beyond_the_ceiling_is_reported_with_its_number():
    """Reported, not just refused: a person told "the service asked for 20
    minutes" can pick another model, while an unexplained failure after
    four silent retries tells them nothing."""
    assert ac._retry_after_outlasts_patience(
        _Error({"retry-after": "1200"})) == 1200.0


def test_blabladors_longest_documented_block_is_recognised():
    assert ac._retry_after_outlasts_patience(
        _Error({"retry-after": "3600"})) == 3600.0


@pytest.mark.parametrize("headers", [None, {}, {"retry-after": "soon"}])
def test_no_usable_header_is_not_a_long_wait(headers):
    """Absent is not long. These must keep falling back to the ladder,
    which is the right instrument for an error nobody stated a wait for."""
    assert ac._retry_after_outlasts_patience(_Error(headers)) is None


def test_the_two_readers_cannot_disagree_about_the_header():
    """One parse, two decisions. A second parser is the shape that drifts:
    exactly one of the pair would then have to be wrong about what the
    header said."""
    exc = _Error({"retry-after": "120"})

    assert ac._retry_after_seconds(exc) == 120.0
    assert ac._retry_after_outlasts_patience(exc) is None

    exc = _Error({"retry-after": "900"})

    assert ac._retry_after_seconds(exc) is None, "too long to sleep"
    assert ac._retry_after_outlasts_patience(exc) == 900.0, "and that is why"


def test_a_real_sdk_error_with_a_long_block_outlasts_patience():
    """On the transport's own case-insensitive header mapping, not a dict."""
    exc = _real_sdk_error(429, {"Retry-After": "1800"})

    assert ac._retry_after_outlasts_patience(exc) == 1800.0
    assert ac._retry_after_seconds(exc) is None


def test_the_retry_path_stops_rather_than_asking_again_sooner():
    """The wiring. A guard that reads the header but still enters the
    retry branch would pass every test above and keep the defect."""
    import inspect

    src = inspect.getsource(ac)
    i = src.index("_stated_wait = _retry_after_outlasts_patience(")
    block = src[i:i + 900]

    assert "_is_transient_api_error" in block, (
        "the guard sits at the retry decision, not somewhere else")
    assert "_stated_wait is None" in block, (
        "a stated long wait must exclude the round from the retry branch")
    assert block.index("_stated_wait = ") < block.index(
        "_is_transient_api_error"), "it is asked before the decision"


def test_the_person_is_told_why_the_round_stopped():
    """Silence here is the failure mode: the round ends, and without a
    notice nobody can tell a rate limit from a broken endpoint."""
    import inspect

    src = inspect.getsource(ac)
    i = src.index("_stated_wait = _retry_after_outlasts_patience(")
    block = src[i:i + 900]

    assert "StreamEvent(type=\"notice\"" in block, "it must be visible"
    assert "lengthens the block" in block, (
        "and it must say why asking again is not the answer")
