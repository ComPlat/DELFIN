"""A model's Markdown, shown on the phone in SimpleX's own markup.

Seen on the phone on 2026-11-12: answers arrived with **bold** and
[text](url) shown literally. SimpleX bolds *one-star* text, italicises
_underscored_ text and has no headings or link syntax.
"""
from __future__ import annotations

from delfin.messenger.render import simplex_text


def test_strong_and_heading_become_simplex_bold():
    assert simplex_text("Hier ist **der Link**.") == "Hier ist *der Link*."
    assert simplex_text("__fett__ und *auch*") == "*fett* und *auch*"
    assert simplex_text("## Ergebnis\nfertig") == "*Ergebnis*\nfertig"


def test_links_become_text_and_url():
    assert simplex_text("[Repo](https://github.com/x/y)") == "Repo (https://github.com/x/y)"
    assert (simplex_text("👉 **[https://github.com/ComPlat/DELFIN](https://github.com/ComPlat/DELFIN)**")
            == "👉 *https://github.com/ComPlat/DELFIN*")


def test_code_and_bullets():
    text = "* eins\n* zwei\n```\nx = **2**\n```\nund `a*b*c`"
    assert simplex_text(text) == "• eins\n• zwei\n```\nx = **2**\n```\nund `a*b*c`"


def test_plain_text_is_unchanged():
    for text in ("Energie -5.0704 Eh, 3 Schritte.", "2 * 3 = 6", "a_b_c", ""):
        assert simplex_text(text) == text
