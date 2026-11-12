"""KN1 phase 1 red test: every shipped operating-knowledge file
(delfin/agent/pack/knowledge/operating/*.md) passes the entry template
(front matter + Situation/Lesson/Why) and cites evidence — a lesson
without evidence never ships (task brief KN1 item 4). Red until the
wave-22 lessons are distilled and shipped.
"""
from pathlib import Path

import delfin.agent.knowledge as kn

PACK = kn.global_dir()
OPERATING = PACK / "operating"


def _operating_entries():
    assert OPERATING.is_dir(), "operating/ missing"
    slugs = [p.stem for p in sorted(OPERATING.glob("*.md"))]
    assert slugs, "no operating entries shipped"
    return slugs


def test_every_operating_entry_passes_the_template():
    topics = kn.read_topics(PACK)
    bad = []
    for slug in _operating_entries():
        e = kn.parse_entry(
            (OPERATING / f"{slug}.md").read_text(encoding="utf-8"),
            slug=slug, layer=kn.GLOBAL)
        if kn.template_findings(e, topics=topics, layer=kn.GLOBAL):
            bad.append(slug)
    assert bad == [], bad


def test_every_operating_entry_carries_evidence():
    for slug in _operating_entries():
        e = kn.parse_entry(
            (OPERATING / f"{slug}.md").read_text(encoding="utf-8"),
            slug=slug, layer=kn.GLOBAL)
        assert e.provenance, f"{slug}: no evidence"
        for p in e.provenance:
            assert p.evidence.split(":", 1)[0] in kn.PUBLIC_EVIDENCE_KINDS, \
                (slug, p.evidence)


def test_every_operating_entry_has_situation_lesson_why():
    for slug in _operating_entries():
        text = (OPERATING / f"{slug}.md").read_text(encoding="utf-8")
        for head in ("Situation:", "Lesson:", "Why:"):
            assert head in text, (slug, head)


def test_every_cited_path_exists():
    # A cited "doc:<path>" evidence pointer points at a file in this
    # repository; a path that does not exist makes the evidence
    # unverifiable (operator, wave 22).
    import re
    missing = []
    for slug in _operating_entries():
        text = (OPERATING / f"{slug}.md").read_text(encoding="utf-8")
        for path in re.findall(r"doc:(\S+?\.(?:py|md|out|json|txt))", text):
            if not Path(path).is_file():
                missing.append((slug, path))
    assert missing == [], missing
