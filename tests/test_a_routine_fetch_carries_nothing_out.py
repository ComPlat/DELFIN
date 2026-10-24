"""A web fetch is routine only when its URL cannot carry data out.

Input: routine_actions._url_routine, reached through action_classifier.screen
(the FR1 routine path, after the user's "never" rules and the injection
guard). Before: any http(s) GET was an allow without a dialog, with any
query string, fragment or path length -- the shape a hostile page uses to
have the model send workspace content to its own host
(``https://host/?d=<file contents>``). Secrets in a URL are refused
elsewhere; ordinary file contents were not.

Now the routine shape is a GET whose URL has no query and no fragment and
a path of at most 200 characters: a DOI, an article, a manual page. A URL
outside that is judged as before (screen returns None), not refused.
"""
from __future__ import annotations

import pytest

from delfin.agent import action_classifier as ac


def _verdict(url):
    return ac.screen({"tool": "web_fetch", "url": url},
                     {"task": "read the literature"})


@pytest.mark.parametrize("url", [
    "https://doi.org/10.1021/acs.jctc.9b00143",
    "https://arxiv.org/abs/2101.00001",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC123/",
])
def test_a_plain_article_fetch_stays_routine(url):
    v = _verdict(url)
    assert v is not None and v.decision == "allow"


@pytest.mark.parametrize("url", [
    "https://collector.example/?d=c2VjcmV0IGNvbnRlbnRz",
    "https://collector.example/page#c2VjcmV0",
    "https://collector.example/" + "a" * 300,
])
def test_a_url_that_can_carry_data_is_judged_not_waved_through(url):
    assert _verdict(url) is None
