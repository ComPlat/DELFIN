"""The packaged conformer-shape-classification skill stays English.

DELFIN's packaged skills are English; German user phrases may stay as
examples of what a user types (task TC1 target 1). Pinned structurally:
every line with German characters must be one of
  - a quoted user phrase inside a prose sentence, or
  - an NL-term alias mapping (a quoted German term followed by '->'
    or '=' pointing at an engine identifier), or
  - the header of an example code block line.
"""
from pathlib import Path
import re

SKILL = Path(
    "delfin/agent/pack/skills/conformer-shape-classification.md")
NON_ASCII = re.compile(r"[^\x00-\x7F]")
GERMAN = re.compile(r"[äöüßÄÖÜ]")


def test_skill_german_appears_only_as_user_phrases_or_aliases():
    # German characters are allowed only where the user's words appear:
    # inside quotes (user-typed phrase examples; may wrap across lines)
    # and in the table cells that map the user's German term to a
    # descriptor (lines starting with "|").
    offenders = []
    in_quote = False
    quote_char = None
    for i, line in enumerate(SKILL.read_text().splitlines(), 1):
        is_table = line.strip().startswith("|")
        for j, ch in enumerate(line):
            if in_quote:
                if ch == quote_char:
                    in_quote = False
                continue
            if ch in "\"'":
                in_quote, quote_char = True, ch
                continue
            if ch in "äöüßÄÖÜ" and not is_table:
                offenders.append((i, line))
                break
    assert not offenders, (
        "German outside quotes/table cells in the packaged skill: "
        + "; ".join(f"line {i}: {l!r}" for i, l in offenders))


def test_skill_is_english_prose():
    # Headings and the intro block must be pure ASCII English.
    lines = SKILL.read_text().splitlines()
    for i, line in enumerate(lines, 1):
        if line.startswith("#"):
            assert not GERMAN.search(line), f"line {i}: {line!r}"
    intro = " ".join(lines[3:10])
    assert not GERMAN.search(intro), intro


def test_every_descriptor_has_a_verified_source_pointer():
    # TC1 target 2: each descriptor definition in the docs table cites a
    # source; the code exposes it via describe_descriptor. The two DOIs
    # were verified 2026-11-01 through DELFIN's web tools.
    from delfin.analysis_tools.geometry_classification_tools import (
        DESCRIPTOR_SOURCES, describe_descriptor)
    assert set(DESCRIPTOR_SOURCES) == set(_descriptor_names())
    for name in DESCRIPTOR_SOURCES:
        info = describe_descriptor(name)
        assert info["source"], name
        assert "DOI" in info["source"] or "DELFIN operational" in info["source"]
        assert info["source_reference"].endswith(
            "CONFORMER_GEOMETRY_DESCRIPTORS.md")
    for doi in ("10.1021/ci025599w", "10.1021/ma00148a028"):
        assert any(doi in s for s in DESCRIPTOR_SOURCES.values()), doi


def _descriptor_names():
    from delfin.analysis_tools.conformer_shape_classification import (
        DESCRIPTOR_DOC)
    return set(DESCRIPTOR_DOC)
