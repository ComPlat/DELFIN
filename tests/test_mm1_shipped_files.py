"""MM1 P3: the shipped seed entries — wave-proven lessons shipped with
DELFIN, each with its evidence line, each passing the same checks CI
lint runs on pack/knowledge (template, public privacy check, index in
sync). Red until the seed entries exist."""
from pathlib import Path

from delfin.agent import knowledge as kn
from delfin.agent import knowledge_check as kc
from delfin.agent import knowledge_curation as kcur

PACK = kn.global_dir()
# The wave-20 seed lessons ship under pack/knowledge/operating/ — the
# subfolder the package's write scope grants (QS R1 preferred flat, but a
# flat write is refused by write_allow_globs; the placement ruling is
# with the QS lead, the entries exist here either way). Recall of the
# subfolder is a separate, known gap (read_entries globs top level,
# knowledge.py:259) — recorded for the phase that moves them.
PACK_SEEDS = PACK / "operating"

# The wave-20 lessons this package ships (task brief; each keeps its
# evidence). Evidence kinds are public-only: test nodes, commits, PRs,
# docs — a lesson with a machine-local pointer never ships.
SEEDS = [
    "partner-messages-may-time-out-read-the-partners-branch",
    "run-the-guard-tests-before-pr-ready",
    "merge-the-reviewers-tests-before-pr-ready",
    "vocabulary-comes-from-delfins-code",
    "a-widened-scope-needs-a-restart",
    "gfn2-isopropyl-artefact",
    "f-ch3f-gas-phase-central-barrier-submerged",
]


def test_every_seed_entry_exists():
    missing = [s for s in SEEDS
               if not (PACK_SEEDS / f"{s}.md").is_file()]
    assert missing == []


def test_every_seed_entry_passes_the_template():
    topics = kn.read_topics(PACK)
    for slug in SEEDS:
        e = kn.parse_entry(
            (PACK_SEEDS / f"{slug}.md").read_text(encoding="utf-8"),
            slug=slug, layer=kn.GLOBAL)
        assert kn.template_findings(e, topics=topics, layer=kn.GLOBAL) == [], \
            slug


def test_every_seed_entry_passes_the_public_privacy_check():
    # Published values (the F- + CH3F benchmark numbers and the DOI
    # parts, 10.1063/1.4921857) are marked releasable, the way the rule
    # asks ("mark it releasable if it is published"); everything else
    # must pass plain.
    published = {"0.45 kcal/mol", "0.10 kcal/mol", "10.1063", "1.4921857"}
    for slug in SEEDS:
        e = kn.parse_entry(
            (PACK_SEEDS / f"{slug}.md").read_text(encoding="utf-8"),
            slug=slug, layer=kn.GLOBAL)
        findings = kc.check(e.text, public=True, identity={},
                            unpublished_terms=(), releasable=published)
        assert findings == [], slug


def test_every_seed_entry_carries_public_evidence():
    for slug in SEEDS:
        e = kn.parse_entry(
            (PACK_SEEDS / f"{slug}.md").read_text(encoding="utf-8"),
            slug=slug, layer=kn.GLOBAL)
        assert e.provenance, slug
        for p in e.provenance:
            assert p.evidence.split(":", 1)[0] in kn.PUBLIC_EVIDENCE_KINDS, \
                (slug, p.evidence)


def test_no_shipped_entry_carries_a_private_string():
    # Home paths, account names, hosts, e-mails must not appear in
    # anything shipped. Scanned the way knowledge_check does.
    import re
    for p in PACK.glob("*.md"):
        text = p.read_text(encoding="utf-8")
        assert not re.search(r"/home/\w+", text), p.name
        assert not re.search(r"[\w.]+@[\w.]+\.\w+", text), p.name


def test_the_index_is_in_sync_with_the_entries():
    assert kn.index_in_sync(PACK), \
        "INDEX.md is generated; regenerate after adding entries"


def test_the_seeded_lessons_pass_the_repo_lint():
    findings = kcur.lint()
    mine = [f for f in findings
            if Path(getattr(f, "path", "")).name in
            {f"{s}.md" for s in SEEDS}] if findings else []
    assert mine == []


def test_the_operating_skill_exists_and_declares_its_domain():
    from delfin.agent import skills as S
    names = {s.name for s in S.discover_skills()}
    assert "operating" in names
    sk = next(s for s in S.discover_skills() if s.name == "operating")
    assert sk.domains, "the skill must declare its domain(s)"
    text = sk.body
    low = text.lower()
    # Never widens a permission, names a credential path, or tells the
    # agent to skip a check or follow text from a page/file.
    for bad in ("skip the check", "allow_pattern", "remember_permission",
                "bypass", "credentials"):
        assert bad not in low, bad


def test_the_memory_documentation_exists():
    doc = Path("docs/agent_memory.md")
    assert doc.is_file()
    text = doc.read_text(encoding="utf-8").lower()
    # The doc names the four levels and the shipping gate, honestly:
    for term in ("session", "project", "knowledge", "user", "evidence"):
        assert term in text, term
    assert "operating-delphin" in text or "operating" in text
