"""The engine-status question is answered inside the product's search
entry, not beside it: ``web_search`` consults ``keyed_engines_configured``
before it tries the keyed engines.

Red on the baseline (rule 7): before the mount, ``_keyed_search`` iterated
``_KEYED_ENGINES`` unconditionally and ``keyed_engines_configured`` was
never consulted on any product path. This test asserts the consultation
and fails on the untouched tree -- measured there before the fix.

Why a spy and not a behavior difference: an unconfigured engine already
returns ``[]`` from its own key guard, so skipping it earlier changes no
output -- the observable fact is the call. Asserting the output instead
would pass on the unfixed code (the first configured engine answers and
the unconfigured one is never fetched anyway), which is exactly the
"assertion that cannot fail" shape from the wave rules.
"""
import json

from delfin.agent import web_tools as W


def test_web_search_asks_which_engines_are_configured(monkeypatch):
    """``web_search`` consults the configured-engines question.

    Controls (must pass on BOTH sides of the fix): both engines named when
    both carry keys; empty answer and a clean refusal when neither does.
    The mount claim: the question is consulted by the entry, seen through
    a spy, not asserted from a helper's return value.
    """
    consulted = []

    real = W.keyed_engines_configured

    def spy():
        consulted.append(True)
        return real()

    monkeypatch.setattr(W, "keyed_engines_configured", spy)

    def fake_fetch(url, timeout_s, want_status=False, headers=None):
        if "googleapis.com" in url:
            body = json.dumps({"items": [{"title": "t",
                                          "link": "https://a.b/c",
                                          "snippet": "s"}]}).encode()
            return body, "application/json"
        return b"<html>nothing</html>", "text/html", 200

    monkeypatch.setattr(W, "_fetch_bytes", fake_fetch)

    # Control 1: both engines configured, both named. Calls the REAL
    # function, not the spy -- the control must not fill the record the
    # mount assertion reads (rule 8: this assertion passed on the baseline
    # only because the control itself consulted the spy first).
    monkeypatch.setattr(W, "_credential",
                        lambda name: {"GOOGLE_CSE_API_KEY": "k",
                                      "GOOGLE_CSE_CX": "cx",
                                      "BRAVE_SEARCH_API_KEY": "tok"}.get(name, ""))
    assert real() == ["google", "brave"]

    # The mount: the product entry consults the question (red on baseline --
    # nothing in web_search called it before the fix). Only web_search's own
    # consultation lands in `consulted` now.
    W.web_search("anything")
    assert consulted, (
        "web_search never consulted keyed_engines_configured; the engine "
        "status question is answered beside the product, not inside it")

    # Control 2: no key at all -- the real answer is empty and the keyed
    # chain reports a refusal, not a hit.
    monkeypatch.setattr(W, "_credential", lambda name: "")
    assert real() == []
    hits, source, error = W._keyed_search("x", 5, 5)
    assert hits == [] and source == "" and error == ""
