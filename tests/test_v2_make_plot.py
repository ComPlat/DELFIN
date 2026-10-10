"""make_plot(spec): the shared DELFIN_CARD: tool result (package V2, phase 3).

make_plot renders a declarative spec (via chat_plots.plot), wraps the figure
as the sandboxed iframe card (chat_plots.to_card) and returns the shared
V1+V2 tool-result string: ``DELFIN_CARD:`` followed by
``{"escape":"html","html":<sandboxed iframe>,"text":<plain path + caption>}``.

The dashboard's role="tool" hook inlines the ``html`` card; a terminal shows
only ``text``. The model never receives the html (operator data-is-never-
instruction rule): text is deliberately plain.

These tests compare structure (the DELFIN_CARD: shape, the escaped sandboxed
iframe, the plain text field, workspace confinement), not pixels.
"""

from __future__ import annotations

import json

import pytest

import delfin.agent.chat_plots as cp
from delfin.agent.chat_plots import SpecError


def _spec() -> dict:
    return {
        "kind": "line",
        "title": "Energies per method",
        "units": "eV",
        "xlabel": "method",
        "ylabel": "energy",
        "x": ["a", "b", "c"],
        "y": [1.0, 2.0, 3.0],
    }


def _payload(out: str) -> dict:
    assert out.startswith("DELFIN_CARD:"), f"no card marker in {out[:60]!r}"
    return json.loads(out[len("DELFIN_CARD:"):])


def test_make_plot_returns_a_shared_card_result(tmp_path):
    out = cp.make_plot(_spec(), out_dir=str(tmp_path))
    assert out.startswith("DELFIN_CARD:")
    obj = _payload(out)
    assert obj["escape"] == "html"


def test_make_plot_html_is_the_sandboxed_iframe_card(tmp_path):
    obj = _payload(cp.make_plot(_spec(), out_dir=str(tmp_path)))
    html = obj["html"]
    assert html.startswith("<iframe")
    # single, isolated iframe: sandbox=allow-scripts, no allow-same-origin.
    assert 'sandbox="allow-scripts"' in html
    assert "allow-same-origin" not in html
    # no event-handler anywhere in the card.
    assert "onerror" not in html and "onclick" not in html and "onload" not in html
    assert "javascript:" not in html
    # the IFrame's OWN attributes (before srcdoc) carry no src=/on*: the figure
    # only travels as a data: URI in the srcdoc value, never loaded by the frame.
    iframe_open = html.split(">", 1)[0]
    own_atts = iframe_open.split('srcdoc="', 1)[0]
    for banned in (" src=", "onerror", "onclick", "onload"):
        assert banned not in own_atts, (banned, own_atts)
    assert " srcdoc=" in iframe_open


def test_make_plot_text_is_plain_and_holds_the_path(tmp_path):
    obj = _payload(cp.make_plot(_spec(), out_dir=str(tmp_path)))
    assert "<" not in obj["text"], "text must be plain, no html"
    assert "Energies per method" in obj["text"]
    assert ".svg" in obj["text"], obj["text"]


def test_make_plot_writes_the_figure_inside_the_workspace(tmp_path):
    cp.make_plot(_spec(), out_dir=str(tmp_path))
    figure = next(tmp_path.rglob("*.svg"))
    assert str(figure).startswith(str(tmp_path))


def test_make_plot_bad_spec_is_refused_clearly(tmp_path):
    before = set(tmp_path.rglob("*"))
    with pytest.raises(SpecError):
        cp.make_plot({"kind": "pie"}, out_dir=str(tmp_path))
    assert set(tmp_path.rglob("*")) == before, "bad spec must not write files"


def test_make_plot_filename_cannot_escape_the_workspace(tmp_path):
    spec = _spec()
    spec["filename"] = "../escape.svg"
    with pytest.raises(SpecError):
        cp.make_plot(spec, out_dir=str(tmp_path))
    assert not (tmp_path.parent / "escape.svg").exists()


def test_make_plot_prefers_the_shared_v1_chat_viewer_emit_point(
        tmp_path, monkeypatch):
    """When the shared emit point (chat_viewer.tool_result) is importable,
    make_plot delegates to it — ONE emit path for V1+V2, not a second one."""
    import types

    import delfin.dashboard

    calls = {}
    fake = types.ModuleType("delfin.dashboard.chat_viewer")

    def tool_result(html: str, text: str) -> str:
        calls["html"], calls["text"] = html, text
        return "DELFIN_CARD:" + json.dumps(
            {"escape": "html", "html": html, "text": text}
        )

    fake.tool_result = tool_result
    # Once V1 is merged, the real chat_viewer is a BOUND ATTRIBUTE of the
    # delfin.dashboard package, so _emit_card's `from delfin.dashboard import
    # chat_viewer` resolves to that attribute and never reads sys.modules.
    # Patch the package attribute itself (raising=False also covers the
    # pre-merge case where it is not yet bound) so the delegation is actually
    # observed on the merged head.
    monkeypatch.setattr(delfin.dashboard, "chat_viewer", fake, raising=False)
    out = cp.make_plot(_spec(), out_dir=str(tmp_path))
    assert calls, "make_plot did not delegate to chat_viewer.tool_result"
    assert calls["text"].startswith("Wrote "), calls["text"]
    assert "<" not in calls["text"]
    assert calls["html"].startswith("<iframe")
    assert out.startswith("DELFIN_CARD:")
    # byte-identity of delegate vs local fallback: the operator-tightened
    # parse_card_result must see exactly the same shape either way.
    fallback = cp._emit_card(calls["html"], calls["text"])
    expected = "DELFIN_CARD:" + json.dumps(
        {"escape": "html", "html": calls["html"], "text": calls["text"]}
    )
    assert fallback == expected
    assert out == expected  # delegated string == documented contract string


def test_fallback_emit_is_byte_identical_to_the_shared_contract(tmp_path):
    """The interim local emit produces exactly `CARD_MARKER + json.dumps(
    {"escape":"html","html":…,"text":…})` with default separators — compact
    ", " and ": ", no sort_keys — so it matches chat_viewer.tool_result byte
    for byte once V1 merges, and the strict parse_card_result is unaffected
    by which path produced the string."""
    cards = {f"<iframe>{i}</iframe>": f"Wrote plot_{i}.svg — line ({i} pts)"
             for i in range(2)}
    for html, text in cards.items():
        got = cp._emit_card(html, text)
        want = "DELFIN_CARD:" + json.dumps(
            {"escape": "html", "html": html, "text": text}
        )
        assert got == want
        assert got.startswith("DELFIN_CARD:")
        # the string stays compact (default separators, no whitespace padding)
        assert ", " in got and ": " in got
