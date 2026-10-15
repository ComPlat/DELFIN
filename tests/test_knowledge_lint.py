"""The CI gate over the shipped global knowledge.

A pull request that adds an entry with a broken template, a topic outside
the taxonomy, more than the word limit, content the public-repository
rules refuse, an exact duplicate, or a stale INDEX.md fails here.
"""
from __future__ import annotations

from delfin.agent import knowledge as kn
from delfin.agent import knowledge_curation as kc


def test_the_global_knowledge_folder_ships_its_taxonomy_and_index():
    folder = kn.global_dir()
    assert (folder / kn.TOPICS_FILE).is_file()
    assert (folder / kn.INDEX_FILE).is_file()
    assert set(kn.DEFAULT_TOPICS) <= set(kn.read_topics(folder))


def test_every_global_entry_passes_the_lint():
    problems = kc.lint(kn.global_dir())
    assert problems == [], "\n".join(problems)


def test_the_folder_is_owned_by_the_knowledge_admin_in_codeowners():
    owners = kn.global_dir().parents[3] / ".github" / "CODEOWNERS"
    lines = [ln.split() for ln in owners.read_text().splitlines()
             if ln.strip() and not ln.startswith("#")]
    assert any(parts[0] == "/delfin/agent/pack/knowledge/" and len(parts) > 1
               for parts in lines)
