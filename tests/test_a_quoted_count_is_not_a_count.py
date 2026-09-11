"""An answer cited the prompt's own rule -- "never >200 lines at once" --
and the truncation caveat told the reader that "200 lines" was a number
whose only source was a cut-short read_file (operator interview,
2026-09-11). A figure inside quotes or code formatting is repeated, not
counted from anything the turn read.
"""

from __future__ import annotations

from delfin.agent import verify_guard as vg


def test_a_count_inside_quotes_is_not_a_claim():
    text = 'Die Regel sagt "Grep before Read — never >200 lines at once" (Provider Profile).'
    assert vg.scan_for_counts_over_truncated_output(text, ["read_file"]) == []
    text2 = "Der Prompt verlangt `never >200 lines`, ich habe 31 Dateien gelesen."
    claims = vg.scan_for_counts_over_truncated_output(text2, ["read_file"])
    assert any("31 Dateien" in c for c in claims), claims
    assert not any("200" in c for c in claims), claims


def test_german_quotation_marks_count_as_quotes():
    text = "Er schreibt \u201eich habe 29 PDF-Dateien verifiziert\u201c, aber der Bericht sagt nichts."
    assert vg.scan_for_counts_over_truncated_output(text, ["fill_series"]) == []


def test_an_unquoted_count_over_truncated_output_still_warns():
    text = "Ich habe 31 PDF-Dateien verifiziert."
    claims = vg.scan_for_counts_over_truncated_output(text, ["fill_series"])
    assert any("31 PDF-Dateien" in c for c in claims), claims


def test_the_quote_state_is_per_line():
    text = 'Zeile eins mit einem " Anfang\nIch habe 12 Ordner gelesen.'
    claims = vg.scan_for_counts_over_truncated_output(text, ["list_files"])
    assert any("12 Ordner" in c for c in claims), claims
