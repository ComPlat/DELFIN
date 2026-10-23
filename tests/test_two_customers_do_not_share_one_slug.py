"""Umlauts were deleted from slugs, so two customers collided.

The slug rule stripped anything outside ``[a-z0-9\\s-]``. That DELETES an
umlaut instead of spelling it, and German distinguishes words by exactly
that letter:

    Müller GmbH -> mller-gmbh      Möller GmbH -> mller-gmbh
    Größe       -> gre             Grüße       -> gre

The filename got a numeric suffix, so nothing crashed. But the slug is
what the frontmatter carries, what a merge and a delete match on, and
what the user reads as the link text — so one customer's note answered
for the other's.

There was also a SECOND slugifier in another module, which returned a
third answer for the same title (``m-ller-gmbh``). One title, two stores,
three keys.
"""

from __future__ import annotations

import pytest

from delfin.agent import german as G
from delfin.agent.deliverables import _slugify as report_slug
from delfin.agent.memory_store import _slugify as memory_slug


_COLLIDING = (
    ("Müller GmbH", "Möller GmbH"),
    ("Größe", "Grüße"),
    ("Schön", "Schon"),
    ("Bäcker", "Becker"),
    ("Müll", "Mull"),
)


@pytest.mark.parametrize("left,right", _COLLIDING)
@pytest.mark.parametrize("slugify", [memory_slug, report_slug, G.slugify],
                         ids=["memory", "report", "shared"])
def test_two_names_that_differ_only_by_an_umlaut_get_two_slugs(
        slugify, left, right):
    assert slugify(left) != slugify(right)


@pytest.mark.parametrize("title,expected", [
    ("Müller GmbH", "mueller-gmbh"),
    ("Möller GmbH", "moeller-gmbh"),
    ("Größe", "groesse"),
    ("Grüße", "gruesse"),
    ("Straße 5", "strasse-5"),
    ("Kostenstellenübersicht", "kostenstellenuebersicht"),
    ("Öffnungszeiten", "oeffnungszeiten"),
])
def test_an_umlaut_is_spelled_out_and_not_dropped(title, expected):
    assert G.slugify(title) == expected


@pytest.mark.parametrize("title", [
    "Müller GmbH", "Größe", "Rechnung 2026/07", "Kostenstelle 4711",
    "Ä Ö Ü", "Straße 5",
])
def test_the_two_stores_agree_on_every_title(title):
    """One rule, one answer. The two used to disagree with each other and
    with the shared one. Only the empty fallback differs, on purpose: a
    nameless memory and a nameless report are different things."""
    assert memory_slug(title, 48) == report_slug(title, 48)


def test_an_empty_title_still_gets_a_usable_name():
    assert memory_slug("") == "memory"
    assert report_slug("") == "report"
    assert G.slugify("///") == "eintrag"


def test_a_slug_stays_within_its_length_and_has_no_trailing_dash():
    long = "Übersicht über alle Kostenstellen des Instituts im Jahr 2026"
    for slugify, limit in ((memory_slug, 60), (report_slug, 48)):
        out = slugify(long)
        assert len(out) <= limit
        assert not out.endswith("-")
        assert out == out.lower()


def test_the_transliteration_covers_the_sharp_s():
    assert G.transliterate("Straße") == "Strasse"
    assert G.transliterate("Müller") == "Mueller"
    assert G.transliterate("ASCII only") == "ASCII only"
