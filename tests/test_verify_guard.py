"""Tests for the runtime verify-enforcement guard.

Covers the two detectors in delfin.agent.verify_guard plus the
dashboard-side _build_verify_hint trigger.  No network / RDKit needed.
"""

from __future__ import annotations

from delfin.agent import verify_guard as vg


# ---------------------------------------------------------------------------
# Ground-truth loading
# ---------------------------------------------------------------------------

def test_namespace_loads_and_is_large():
    ns = vg.load_orca_namespace()
    assert len(ns) > 1000          # the committed manual extract is ~1600
    assert "nel" in ns             # real %casscf keyword
    assert "maxiter" in ns         # real %scf/%geom keyword


def test_fakes_are_absent_from_namespace():
    ns = vg.load_orca_namespace()
    fakes = vg.known_fake_keywords()
    assert "nactel" in fakes
    assert "nactorb" in fakes
    # known_fake_keywords must never contain a real manual keyword
    assert not (fakes & ns)


# ---------------------------------------------------------------------------
# Fake-keyword detector (high precision)
# ---------------------------------------------------------------------------

def test_fake_keyword_flagged_in_prose():
    flags = vg.scan_for_unverified_keywords(
        "Im %casscf-Block nutzt du Nactel für aktive Elektronen."
    )
    assert any(f.keyword == "nactel" and f.reason == "fake" for f in flags)


def test_real_keywords_not_flagged():
    flags = vg.scan_for_unverified_keywords(
        "Im %casscf-Block: nel, norb, mult und nroots."
    )
    assert flags == []


def test_ambiguous_word_not_flagged_in_prose():
    # "multiplicity" is a forbid entry but also an ordinary word — it must
    # NOT fire when used as plain prose.
    flags = vg.scan_for_unverified_keywords(
        "Die multiplicity is 3 and there are excited states in ORCA."
    )
    assert flags == []


def test_ambiguous_word_flagged_as_keyword():
    # ...but it MUST fire when presented as a keyword assignment.
    flags = vg.scan_for_unverified_keywords(
        "Setze `multiplicity = 3` im %casscf-Block."
    )
    assert any(f.keyword == "multiplicity" for f in flags)


# ---------------------------------------------------------------------------
# Unknown-keyword detector (conservative)
# ---------------------------------------------------------------------------

def test_unknown_keyword_flagged_in_orca_context():
    flags = vg.scan_for_unverified_keywords(
        "In ORCA setzt du `foobarbaz = 5` im %scf-Block."
    )
    assert any(f.keyword == "foobarbaz" and f.reason == "unknown"
               for f in flags)


def test_unknown_token_not_flagged_outside_orca_context():
    # No ORCA / % context → the unknown detector stays silent.
    flags = vg.scan_for_unverified_keywords(
        "Die `schema.json` Datei ist groß, width = 500."
    )
    assert flags == []


def test_detect_unknown_false_restricts_to_fakes():
    text = "In ORCA: `foobarbaz = 1` und Nactel im %casscf."
    only_fakes = vg.scan_for_unverified_keywords(text, detect_unknown=False)
    assert all(f.reason == "fake" for f in only_fakes)
    assert any(f.keyword == "nactel" for f in only_fakes)


def test_empty_text_returns_no_flags():
    assert vg.scan_for_unverified_keywords("") == []
    assert vg.scan_for_unverified_keywords("   ") == []


# ---------------------------------------------------------------------------
# Correction feedback + flag message
# ---------------------------------------------------------------------------

def test_correction_feedback_lists_keywords():
    flags = [vg.VerifyFlag("nactel", "fake", "ist kein echtes Keyword")]
    fb = vg.correction_feedback(flags)
    assert "nactel" in fb
    assert "search_docs" in fb


def test_flag_message_starts_with_warning():
    flag = vg.VerifyFlag("nactel", "fake", "ist kein echtes Keyword.")
    assert flag.message().startswith("⚠️ Verify:")
    assert "nactel" in flag.message()


# ---------------------------------------------------------------------------
# Dashboard verify-hint trigger
# ---------------------------------------------------------------------------

def test_verify_hint_fires_on_keyword_question():
    from delfin.dashboard.tab_agent import _build_verify_hint
    assert _build_verify_hint(
        "welche EXAKTEN keyword-namen nutzt ORCA im %casscf-block?"
    )
    assert _build_verify_hint("wie heißt das keyword für aktive Orbitale?")


def test_verify_hint_silent_on_navigation():
    from delfin.dashboard.tab_agent import _build_verify_hint
    assert _build_verify_hint("wechsel zu Submit") == ""
    assert _build_verify_hint("erkläre mir kurz wie DFT funktioniert") == ""
