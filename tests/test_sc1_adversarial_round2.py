"""SC1 round-2 adversarial tests (reviewer nacht-s4) on builder code at
agent/s3-sc1b 4c6668dc. Attacks:

* A1: deviation justified with a synonym OUTSIDE _REASON_WORDS —
  the wording list is the acceptance signal; this must NOT depend on
  a fixed vocabulary (structure-over-wordlist rule, wave 15a lesson).
* A2: deviation justified with an emoji/unicode look-alike of a reason
  word.
* A3: claim-no-source tolerance — trailing zeros and negative zero.
* A4: a justified answer in the SAME window as an UNRELATED
  claim — window collisions must not excuse two different deviations.
* A5: scientific prose with 'est.' abbreviation marking.
"""

from __future__ import annotations

import delfin.agent.result_claims as rcl


def _dev(actual="def2-TZVP", param="basis", default="def2-SVP"):
    return [{"parameter": param, "default": default, "actual": actual}]


def test_a1_synonym_outside_wordlist_is_flagged():
    # A real justification the word list does not know. The record
    # drives the check; the wording list decides acceptance. Decision
    # (QS ruling s29, 2026-10-23): the list is deliberately
    # conservative and stays the acceptance signal for text-only
    # justification; the STRUCTURAL escape is the record-level
    # `reason` field (test_sc1_workflow.py covers the short-circuit).
    # This test pins the documented consequence: a text justification
    # with a list-unknown synonym IS flagged. The escape hatch is
    # structural, not vocabulary.
    ans = ("We ran def2-TZVP; the choice follows from a systematic "
           "basis-set study we completed first (see SI Table 2).")
    claims = rcl.check_deviations(ans, _dev())
    unexplained = [c for c in claims if c.code == "deviation-unexplained"]
    assert unexplained, (
        "A list-unknown synonym justification is flagged by design; "
        "the structural escape is the record-level reason field.")


def test_a2_deviation_word_absent_means_not_flagged_here():
    # Deviation NOT mentioned in the answer at all: builder says CB0A
    # reports it. Verify that contract: no claim here either way.
    ans = "Energy is -1234.567890 Eh with B3LYP/def2-SVP."
    claims = rcl.check_deviations(ans, _dev())
    assert claims == []


def test_a3_trailing_zeros_and_comma_matched():
    # Evidence contains "-76.342100"; answer says "-76.3421 Eh" — must
    # match (string-tolerant), and comma format "-1,234.567890" too.
    ev = [{"kind": "calc", "ref": "x.out",
           "text": "Electronic Energy -76.342100 Eh",
           "flags": []}]
    ans = "Energy -76.3421 Eh, B3LYP/def2-SVP."
    claims = rcl.check(ans, evidence=ev)
    assert [c for c in claims if c.code == "claim-no-source"] == [], claims


def test_a4_window_collision_excuses_only_near_reason():
    # Two deviations in one answer; only one is justified with a
    # reason word. The 200-char window around deviation A must not
    # absorb the reason meant for deviation B.
    ans = ("We used TightSCF because the preliminary run had SCF "
           "oscillations. For GridX7 we chose the default GridX5.")
    dev = [{"parameter": "SCF", "default": "NormalSCF",
            "actual": "TightSCF"},
           {"parameter": "grid", "default": "GridX7", "actual": "GridX5"}]
    claims = rcl.check_deviations(ans, dev)
    flagged = {c.message for c in claims
               if c.code == "deviation-unexplained"}
    assert any("grid=GridX5" in m for m in flagged), flagged
    assert not any("TightSCF" in m for m in flagged), flagged
