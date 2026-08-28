"""Three fixes from one morning's bug reports, and one shape between them:
a rule that lived only in prose while the mechanism allowed its opposite.

Reports 20260828-131034 and -131603 (user ka_gf0106, dashboard_agent):

    14 calls,  7 failed  -> all seven web_search
    79 calls, 31 failed  -> 31 of 32 web_search
    turn 2: 68 tool rounds, 1 419 632 input tokens

Every failure was HTTP 202 from DuckDuckGo. The returned error already said
"tell the user the search could not run" -- and the model reworded the query
and tried again, 32 times. Rewording changed the error TEXT, so the stream
loop's identical-error abort kept resetting; the turn ran to the round cap.

Measured 2026-08-28, same host, same minute, same query:

    Mozilla/5.0 (DELFIN-KIT-agent; +https://kit.edu)   HTTP 202,  0 results
    Mozilla/5.0 (X11 …) Chrome/124 … DELFIN/1.0        HTTP 200, 10 results
    Mozilla/5.0 (X11 …) Chrome/124 …                   HTTP 200, 10 results
    Mozilla/5.0 (X11 …) Firefox/127.0                  HTTP 202,  0 results

So the 202 was our own User-Agent, not the network and not the endpoint --
and the comment in the code claiming otherwise was true but misattributed.
A browser-shaped agent that still names DELFIN gets results, so the tool is
repaired rather than retired. 202 still appears under rapid repeat queries,
which is what the second fix is for.

Nothing here touches the network.
"""

from __future__ import annotations

from delfin.agent import web_tools
from delfin.agent.api_client import _error_signature


# ---------------------------------------------------------------------------
# The cause: our own User-Agent
# ---------------------------------------------------------------------------

def test_the_agent_is_browser_shaped():
    """The measured discriminator. A Firefox-shaped string was refused too,
    so this is not "any browser will do" -- it is the shape that was tested."""
    ua = web_tools._USER_AGENT
    assert "Chrome/" in ua and "AppleWebKit" in ua and "Safari/" in ua


def test_the_agent_still_says_who_is_calling():
    """A search tool that hides its identity would be the wrong trade for a
    fix that does not need it: the DELFIN token still measured HTTP 200."""
    ua = web_tools._USER_AGENT
    assert "DELFIN" in ua and "kit.edu" in ua


def test_the_old_agent_is_gone():
    assert "DELFIN-KIT-agent" not in web_tools._USER_AGENT


def test_the_error_no_longer_invites_a_reworded_retry():
    """The returned text used to end at "ask them for a URL", which left
    rewording as an obvious next move. It now says plainly that the reason
    does not depend on the wording."""
    import inspect

    src = inspect.getsource(web_tools.web_search)
    assert "does not depend on the wording" in src


# ---------------------------------------------------------------------------
# The mechanism: rewording no longer resets the loop guard
# ---------------------------------------------------------------------------

def _refusal(query: str) -> str:
    return ('{"error": "the search backend refused this request (HTTP 202: '
            'an anti-bot challenge, not an empty result set)", '
            f'"query": "{query}", "results": [], "result_count": 0}}')


def test_the_same_reason_is_one_failure_however_it_is_worded():
    """The finding: 31 refusals, one reason, 31 different query strings --
    so the identical-error abort saw 31 different errors and never fired."""
    a = _error_signature(_refusal("ORCA imaginary frequencies"))
    b = _error_signature(_refusal("SCF convergence troubleshooting guide"))
    assert a == b


def test_a_genuinely_different_failure_is_still_different():
    """The abort must not fire on a model that is making progress through
    different problems -- that would end turns that were going fine."""
    same = _error_signature(_refusal("x"))
    other = _error_signature('{"error": "file not found", "query": "x"}')
    assert same != other


def test_every_echoed_argument_is_dropped_not_only_the_query():
    """bash echoes ``command``, the file tools echo ``path`` -- the same
    reset would happen there the moment one of them starts failing."""
    for field in ("command", "path", "file_path", "url", "pattern"):
        one = f'{{"error": "denied", "{field}": "first thing"}}'
        two = f'{{"error": "denied", "{field}": "a completely other thing"}}'
        assert _error_signature(one) == _error_signature(two), field


def test_the_heads_up_advice_still_does_not_mask_a_loop():
    """Why the old normalisation existed; it must survive the new one."""
    base = '{"error": "denied", "heads_up": "try something else"}'
    plain = '{"error": "denied"}'
    assert _error_signature(base) == _error_signature(plain)


def test_the_signature_is_what_the_loop_actually_compares():
    """A helper nobody calls would be a fix in name only."""
    import inspect

    from delfin.agent import api_client

    src = inspect.getsource(api_client)
    assert "_error_signature(r) for r in _round_results" in src


# ---------------------------------------------------------------------------
# Where the reports were filed
# ---------------------------------------------------------------------------

def test_a_remote_path_that_already_names_the_archive_is_not_doubled():
    """Two reports from a second user sat in AGENT_BUGS/AGENT_BUGS/ for
    hours. find_unsolved iterates the immediate subdirs and skips anything
    without a report.json, so nobody's triage list ever showed them.

    People configure the remote path both ways, and appending blindly is
    only correct for one of them."""
    from delfin.agent import bug_report

    seen = {}

    def _fake_run(cmd, **kw):
        seen.setdefault("cmds", []).append(cmd)

        class _R:
            returncode = 0
            stderr = ""
        return _R()

    import unittest.mock as _mock
    with _mock.patch.object(bug_report.subprocess, "run", _fake_run), \
            _mock.patch.object(bug_report.Path, "is_dir", lambda self: True):
        ok, where = bug_report.push_report_to_remote(
            "/tmp/report-dir", host="h", user="u",
            remote_path="/archive/AGENT_BUGS")
    assert ok
    assert where.count("AGENT_BUGS") == 1, where


def test_a_parent_remote_path_still_gets_the_archive_appended():
    """The other half of the same setting, which must keep working."""
    from delfin.agent import bug_report
    import unittest.mock as _mock

    def _fake_run(cmd, **kw):
        class _R:
            returncode = 0
            stderr = ""
        return _R()

    with _mock.patch.object(bug_report.subprocess, "run", _fake_run), \
            _mock.patch.object(bug_report.Path, "is_dir", lambda self: True):
        ok, where = bug_report.push_report_to_remote(
            "/tmp/report-dir", host="h", user="u", remote_path="/archive")
    assert ok
    assert "/archive/AGENT_BUGS/" in where, where


# ---------------------------------------------------------------------------
# The round cap that let 68 rounds happen
# ---------------------------------------------------------------------------

def test_unset_means_unset_and_the_profile_decides():
    """The setting default is None for a documented reason. The settings tab
    showed 500 for an unset value and wrote it back on save, so opening the
    tab for any other field pinned 500 and every per-model profile died --
    measured on two installations, neither of which chose 500."""
    from delfin import user_settings

    assert user_settings.DEFAULT_SETTINGS["agent"]["max_tool_rounds"] is None


def test_the_settings_control_can_express_unset():
    """A plain number field cannot say "unset"; -1 is the value that can."""
    import inspect

    from delfin.dashboard import tab_settings

    src = inspect.getsource(tab_settings)
    assert "value=-1, min=-1" in src
    assert "None if _rounds_val < 0 else _rounds_val" in src
    assert "-1 if _rounds is None else int(_rounds)" in src
