"""What a fetched page runs must not arrive as what the page says.

``web_fetch`` reduces HTML to plain text and hands that text to the model, so
the boundary between a page's prose and a page's code is where a stranger's
instructions would enter the conversation. Script and style blocks are dropped
for that reason, not for tidiness.

The end tag was matched literally as ``</script>``. HTML allows whitespace
before the ``>``, so a page closing with ``</script >`` kept its block: the
generic tag rule then removed the two tags around it and left the JavaScript
standing in the middle of the extracted text, indistinguishable from a
paragraph. Whether that is done by a careless templating engine or on purpose,
the outcome is the same.

An unclosed ``<script>`` had the same effect and is treated the same way here
as a browser treats it — everything after it belongs to the script.
"""

import pytest

from delfin.agent.web_tools import _html_to_text


def _text(html: str) -> str:
    return _html_to_text(html.encode("utf-8"))


# --------------------------------------------------------------------------
# the end tag, in the spellings HTML permits
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "close",
    ["</script>", "</script >", "</script  >", "</script\t>", "</script\n>",
     "</SCRIPT>", "</Script >"],
)
def test_a_script_block_is_dropped_however_its_end_tag_is_spelled(close):
    text = _text(f"<p>before</p><script>alert('INJECTED'){close}<p>after</p>")

    assert "INJECTED" not in text
    assert "before" in text and "after" in text


@pytest.mark.parametrize("close", ["</style>", "</style >", "</style\n>"])
def test_a_style_block_is_dropped_however_its_end_tag_is_spelled(close):
    text = _text(f"<p>before</p><style>.a{{content:'INJECTED'}}{close}<p>after</p>")

    assert "INJECTED" not in text
    assert "before" in text and "after" in text


def test_the_prose_around_a_dropped_block_survives():
    text = _text(
        "<h1>Title</h1><script >var x = 1;</script ><p>Body text.</p>"
    )

    assert "Title" in text and "Body text." in text
    assert "var x" not in text


# --------------------------------------------------------------------------
# a block that never closes
# --------------------------------------------------------------------------

def test_an_unclosed_script_runs_to_the_end_of_the_document():
    text = _text("<p>before</p><script>alert('INJECTED'); more(); ")

    assert "INJECTED" not in text
    assert "before" in text


def test_an_unclosed_style_runs_to_the_end_of_the_document():
    text = _text("<p>before</p><style>.a{content:'INJECTED'}")

    assert "INJECTED" not in text
    assert "before" in text


# --------------------------------------------------------------------------
# the block, whatever it holds
# --------------------------------------------------------------------------

def test_a_script_carrying_its_own_end_tag_in_a_string_still_ends():
    """A page can write "</scr" + "ipt>" to embed one. The first real end tag
    is still the end of the block, which is what a browser does too."""
    text = _text(
        "<p>before</p><script>var s = '<' + '/script>';</script >"
        "<p>after</p>"
    )

    assert "after" in text


def test_a_script_with_attributes_is_still_recognised():
    text = _text(
        "<script type='text/javascript' async>alert('INJECTED')</script >"
        "<p>after</p>"
    )

    assert "INJECTED" not in text
    assert "after" in text


def test_several_blocks_in_one_page_are_all_dropped():
    text = _text(
        "<script>a('ONE')</script >"
        "<p>keep</p>"
        "<script>b('TWO')</script>"
        "<p>this</p>"
        "<style>.c{d:'THREE'}</style >"
    )

    assert "ONE" not in text and "TWO" not in text and "THREE" not in text
    assert "keep" in text and "this" in text
