"""A streaming answer is formatted as it grows, not swapped in at the end.

Measured in a real page (2026-10-22): while the reply streamed, the chat
showed it as raw markdown (``**``, ``| --- |``) in a dark box -- the
transparent streaming style lost to the more specific code-block rule
``.delfin-chat-agent pre`` -- and flipped to the formatted table only when
the turn ended. Rendering markdown costs about 0.5 ms per KB here, so the
finished paragraphs are rendered while the last one is still arriving.
"""

from __future__ import annotations

from delfin.dashboard import tab_agent


def test_finished_paragraphs_are_rendered_and_the_growing_one_is_text():
    text = ("**anti** is more stable.\n\n"
            "| Quantity | Value |\n|---|---:|\n| gap | 0.60 |\n\n"
            "GFN2-xTB underestim")
    html = tab_agent._streaming_html(text)
    assert "<strong>anti</strong>" in html
    assert "<table" in html
    tail = html.split("delfin-streaming-pre", 1)[1]
    assert "GFN2-xTB underestim" in tail and "<table" not in tail


def test_a_table_shows_its_finished_rows_while_the_next_one_arrives():
    text = ("Result:\n\n| Quantity | Value |\n|---|---:|\n| gap | 0.60 |\n"
            "| barrier | 2.6 |\n| syn | 4.9")
    html = tab_agent._streaming_html(text)
    head, tail = html.split("delfin-streaming-pre", 1)
    assert "<table" in head and "barrier" in head
    assert "| syn | 4.9" in tail and "<table" not in tail


def test_an_open_code_fence_stays_raw_until_it_closes():
    text = "Here is the input:\n\n```\n! r2SCAN-3c\n* xyzfile 0 1"
    html = tab_agent._streaming_html(text)
    assert "! r2SCAN-3c" in html.split("delfin-streaming-pre", 1)[1]
    assert "<code" not in html.split("delfin-streaming-pre", 1)[0]


def test_markup_in_the_growing_part_is_escaped():
    html = tab_agent._streaming_html("Done.\n\n<img src=x onerror=alert(1)>")
    assert "<img src=x" not in html and "&lt;img" in html


def test_a_long_answer_keeps_its_finished_formatting_across_20000_characters():
    paragraph = "## Result\n\n**Stable** paragraph with `code`.\n\n"
    before = paragraph * (19900 // len(paragraph))
    after = before + paragraph * 10
    assert len(before) < 20000 < len(after)
    prefix = tab_agent._streaming_html(before + 'tail').split('<pre')[0]
    html = tab_agent._streaming_html(after + 'tail')
    assert html.startswith(prefix)
    assert '<strong>Stable</strong>' in html


def test_tokens_in_the_tail_reuse_completed_markdown(monkeypatch):
    tab_agent._streaming_markdown_head.cache_clear()
    calls = []
    monkeypatch.setattr(tab_agent, '_md_to_html', lambda text: calls.append(text) or '<p>head</p>')
    try:
        for tail in ['g', 'gr', 'growing']:
            assert '<p>head</p>' in tab_agent._streaming_html('head\n\n' + tail)
        assert calls == ['head']
    finally:
        tab_agent._streaming_markdown_head.cache_clear()


def test_the_streaming_text_is_not_styled_as_a_code_block():
    css = " ".join(tab_agent._AGENT_CSS.split())
    assert ".delfin-chat-agent pre.delfin-streaming-pre {" in css
    rule = css.split(".delfin-chat-agent pre.delfin-streaming-pre {", 1)[1]
    rule = rule.split("}", 1)[0]
    assert "background: transparent" in rule and "font-family: inherit" in rule
