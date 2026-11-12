"""URL grounding: a URL in the final answer must be backed by a tool result
of this session (a web_search hit's ``url`` field, the exact URL a web_fetch
requested, or a URL that appears in a file the model read). An unsupported
URL is removed or marked, with the model told once.

Live failure that shipped this feature (kit.mistral-small): asked "Kannst du
googeln?", the answer linked https://developers.google.com/custom-search/v1/
overview -- a URL no tool returned; DELFIN's web_search is not Google. That
URL is the canonical unsupported case here.
"""

from __future__ import annotations

import json
import textwrap

import pytest

from delfin.agent import api_client as A
from delfin.agent import verify_guard as vg

LIVE_INVENTED = "https://developers.google.com/custom-search/v1/overview"


# ---------------------------------------------------------------------------
# Scanner -- the live example
# ---------------------------------------------------------------------------


def test_live_invented_url_flagged_when_unsupported():
    flags = vg.scan_for_ungrounded_url_claims(
        f"Google it here: {LIVE_INVENTED}", observed_urls=set())
    assert [f.url for f in flags] == [LIVE_INVENTED]
    assert flags[0].kind == "unobserved"
    assert "tool result" in flags[0].message()


def test_live_invented_url_not_flagged_when_a_search_returned_it():
    flags = vg.scan_for_ungrounded_url_claims(
        f"Google it here: {LIVE_INVENTED}", observed_urls={LIVE_INVENTED})
    assert flags == []


def test_url_flagged_but_normalized_match_clears_it():
    # Trailing slash / trailing dot / fragment differences are the same URL.
    flags = vg.scan_for_ungrounded_url_claims(
        f"see {LIVE_INVENTED}/. for details", observed_urls={LIVE_INVENTED})
    assert flags == []


def test_url_without_ledger_flagging_available():
    # Mirror of the location guard: "ledger kept but empty" enforces;
    # "no ledger at all" (backends that keep none) stays silent.
    assert vg.scan_for_ungrounded_url_claims(
        "see http://example.com/x", ledger_available=False) == []
    assert vg.scan_for_ungrounded_url_claims(
        "see http://example.com/x", observed_urls=set()) != []


def test_fenced_and_blockquoted_urls_skipped():
    text = (
        "```\n"
        f"curl {LIVE_INVENTED}\n"
        "```\n"
        f"> quoted: {LIVE_INVENTED}\n"
        "outside here\n"
    )
    flags = vg.scan_for_ungrounded_url_claims(
        text, observed_urls={"https://outside.example/"})
    # The fenced copy is scrubbed; the blockquote line is blanked. The plain
    # line "outside here" carries no URL, so nothing survives.
    assert flags == []


def test_url_cap_and_dedupe():
    urls = [f"https://unseen{i}.example/p" for i in range(10)]
    flags = vg.scan_for_ungrounded_url_claims(
        " ".join(urls), observed_urls=set())
    assert len(flags) == 6  # default cap
    # de-duplicated: same URL repeated counts once
    flags2 = vg.scan_for_ungrounded_url_claims(
        f"a {urls[0]} b {urls[0]} c", observed_urls=set())
    assert len(flags2) == 1


def test_hedged_url_is_still_flagged():
    # A URL is a link, not a hedged claim: saying "maybe wrong" does not make
    # an invented link sourced. Unlike location claims, hedging does not excuse.
    flags = vg.scan_for_ungrounded_url_claims(
        f"Unsure but likely: {LIVE_INVENTED}", observed_urls=set())
    assert [f.url for f in flags] == [LIVE_INVENTED]


def test_web_fetched_url_observed_not_flagged():
    flags = vg.scan_for_ungrounded_url_claims(
        "From the page I fetched:", observed_urls={"https://docs.python.org/3/"},
        ledger_available=True)
    assert flags == []


def test_color_markdown_link_target_flagged_when_unsupported():
    flags = vg.scan_for_ungrounded_url_claims(
        "[Custom Search](https://developers.google.com/custom-search/v1/overview)",
        observed_urls=set())
    assert len(flags) == 1
    assert flags[0].url == LIVE_INVENTED


# ---------------------------------------------------------------------------
# Ledger -- what a tool result contributes as observed URLs
# ---------------------------------------------------------------------------


def test_web_search_hits_are_observed():
    obs: set = set()
    result = json.dumps({
        "query": "orca", "source": "ddg",
        "results": [{"title": "a", "url": "https://a.example/", "snippet": "x"},
                     {"title": "b", "url": "https://b.example/p", "snippet": "y"}],
    })
    A._observe_read_urls(obs, "web_search", {"query": "orca"}, result)
    assert obs == {"https://a.example/", "https://b.example/p"}


def test_web_search_error_no_urls():
    obs: set = set()
    A._observe_read_urls(obs, "web_search", {"query": "x"},
                         '{"error": "failed"}')
    assert obs == set()


def test_web_fetch_url_is_observed():
    obs: set = set()
    result = json.dumps({"url": "https://docs.python.org/3/",
                         "status": 200, "text": "body", "truncated": False})
    A._observe_read_urls(obs, "web_fetch",
                         {"url": "https://docs.python.org/3/"}, result)
    assert obs == {"https://docs.python.org/3/"}


def test_web_fetch_error_still_observes_the_fetched_url():
    # A refused/failed fetch names the URL that was actually requested; the
    # model can still truthfully cite "the page at X (which I could not read)".
    obs: set = set()
    result = json.dumps({"error": "HTTP 403", "url": "https://denied.example/"})
    A._observe_read_urls(obs, "web_fetch", {"url": "https://denied.example/"}, result)
    assert obs == {"https://denied.example/"}


def test_urls_in_a_read_file_are_observed():
    obs: set = set()
    content = textwrap.dedent(f"""\
        # README
        See {LIVE_INVENTED} here and https://a.example/ there.
        """)
    A._observe_read_urls(obs, "read_file", {"path": "README.md"}, content)
    assert obs == {LIVE_INVENTED, "https://a.example/"}


def test_urls_in_a_grep_result_are_observed():
    obs: set = set()
    hits = "README.md:1: See https://grep.example/ in the docs\n"
    A._observe_read_urls(obs, "grep_file", {"pattern": "https"}, hits)
    assert "https://grep.example/" in obs


def test_invented_url_not_added_by_shell_read_with_no_urls():
    obs: set = set()
    A._observe_read_urls(obs, "bash", {"command": "ls"}, "a.md\nb.py\n")
    assert obs == set()


# ---------------------------------------------------------------------------
# Feedback + caveat
# ---------------------------------------------------------------------------


def test_url_claim_feedback_names_the_url():
    flags = vg.scan_for_ungrounded_url_claims(
        f"here {LIVE_INVENTED}", observed_urls=set())
    fb = vg.url_claim_feedback(flags)
    assert LIVE_INVENTED in fb
    assert "Automatic check" in fb
    assert "unverified" in fb.lower() or "remove" in fb.lower()


def test_grounding_caveat_includes_url_flags():
    flags = vg.scan_for_ungrounded_url_claims(
        f"here {LIVE_INVENTED}", observed_urls=set())
    caveat = vg.grounding_caveat([], [], [], url_flags=flags)
    assert caveat
    assert LIVE_INVENTED in caveat
    assert vg.grounding_caveat([], [], [], url_flags=[]) == ""


def test_url_claim_flag_is_frozen():
    flags = vg.scan_for_ungrounded_url_claims(
        f"here {LIVE_INVENTED}", observed_urls=set())
    assert flags
    with pytest.raises(Exception):
        flags[0].kind = "observed"  # type: ignore[misc]
