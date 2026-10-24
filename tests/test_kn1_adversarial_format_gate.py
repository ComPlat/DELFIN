"""Adversarial format gate for pack/knowledge/operating/*.md (package KN1).

Loads EVERY operating knowledge file and checks the format the existing
files establish:
- front matter with title, topic, scope, provenance, review_by
- provenance cites a real artifact (doc:, PR, commit, test id or run)
- body has a Situation, a Lesson and a Why paragraph
- no duplicate of an existing lesson (title must be unique)

Reviewer nacht-s17, wave 22 stage 3a. Builder files that violate any of
these must fail this test.
"""
from __future__ import annotations

import re
from pathlib import Path

import yaml

import pytest

OPERATING = Path(__file__).resolve().parents[1] / "delfin" / "agent" / "pack" / "knowledge" / "operating"

REQUIRED_FIELDS = ("title", "topic", "scope", "provenance", "review_by")
REQUIRED_SECTIONS = ("Situation:", "Lesson:", "Why:")

# evidence markers a provenance line must carry
EVIDENCE_RE = re.compile(
    r"(doc:[\w/.]+|pr:\d+|PR ?#\d+|#\d+|commit[ :][0-9a-f]{7,}|run[ :]|DOI \S+)", re.IGNORECASE
)


def _operating_files() -> list[Path]:
    files = sorted(OPERATING.glob("*.md"))
    assert files, "no operating knowledge files found at all"
    return files


def test_every_operating_file_has_required_front_matter():
    for path in _operating_files():
        text = path.read_text(encoding="utf-8")
        assert text.startswith("---"), f"{path.name}: no front matter"
        end = text.find("\n---", 3)
        assert end > 0, f"{path.name}: front matter not closed"
        meta = yaml.safe_load(text[3:end])
        for field in REQUIRED_FIELDS:
            assert field in meta, f"{path.name}: missing front-matter field {field!r}"
            assert str(meta[field]).strip(), f"{path.name}: empty front-matter field {field!r}"


def test_every_operating_file_cites_evidence_in_provenance():
    for path in _operating_files():
        text = path.read_text(encoding="utf-8")
        end = text.find("\n---", 3)
        meta = yaml.safe_load(text[3:end])
        provenance = str(meta.get("provenance", ""))
        assert EVIDENCE_RE.search(provenance), (
            f"{path.name}: provenance {provenance!r} cites no evidence "
            "(doc:, PR#, commit or run)"
        )


def test_every_operating_file_has_situation_lesson_why():
    for path in _operating_files():
        text = path.read_text(encoding="utf-8")
        for section in REQUIRED_SECTIONS:
            assert section in text, f"{path.name}: missing {section!r} paragraph"


def test_every_operating_file_mentions_a_delfin_code_or_pr_anchor():
    # A lesson for DELFIN agents must point at where it held: a delfin path,
    # a test file, a PR number or a CONTROL key.
    anchor = re.compile(
        r"(delfin/|tests/|docs/|CONTROL|PR ?#\d+|#\d+|\.py|\.md|pr:\d+|DOI \S+|GFN2)", re.IGNORECASE
    )
    for path in _operating_files():
        body = path.read_text(encoding="utf-8")
        assert anchor.search(body), (
            f"{path.name}: no DELFIN code/test/PR anchor anywhere in the file"
        )


def test_no_two_operating_files_share_a_title():
    titles = []
    for path in _operating_files():
        text = path.read_text(encoding="utf-8")
        end = text.find("\n---", 3)
        meta = yaml.safe_load(text[3:end])
        titles.append(str(meta["title"]).strip().lower())
    assert len(titles) == len(set(titles)), "duplicate lesson titles: extend the existing file instead"


def test_operating_lessons_are_plain_english_no_secret_shapes():
    forbidden = re.compile(
        r"(\b(api[_-]?key|token|password|secret)\b\s*[:=])|(\bAKIA[0-9A-Z]{16}\b)", re.IGNORECASE
    )
    for path in _operating_files():
        text = path.read_text(encoding="utf-8")
        hit = forbidden.search(text)
        assert hit is None, f"{path.name}: looks like a secret shape near {hit.group(0)!r}"


def test_review_by_dates_are_future_iso_dates():
    for path in _operating_files():
        text = path.read_text(encoding="utf-8")
        end = text.find("\n---", 3)
        meta = yaml.safe_load(text[3:end])
        review_by = str(meta["review_by"])
        assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", review_by), f"{path.name}: review_by not ISO date"
        assert review_by >= "2026-10-26", f"{path.name}: review_by already in the past"
