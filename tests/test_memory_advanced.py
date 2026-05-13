"""Tests for the second-wave memory features:

- ``[[name]]`` cross-link resolution at prompt build time
- Auto-classification when ``/remember`` is called without a type prefix
- Stale-reference verification for ``/memories verify``
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import memory_store as ms


# ---------------------------------------------------------------------------
# Auto-classification (parse_memory_type when no explicit prefix)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "text, expected",
    [
        ("Don't mock the DB", "feedback"),
        ("Never amend a pushed commit", "feedback"),
        ("Always run pytest after edits", "feedback"),
        ("Prefer one bundled PR for refactors", "feedback"),
        ("Stop summarising every diff", "feedback"),
        ("Deadline 2026-05-20 for the milestone", "project"),
        ("Mobile release branch cut on Friday", "project"),
        ("Sprint freeze starts next week", "project"),
        ("Pipeline bugs tracked in Linear INGEST", "reference"),
        ("Dashboard at https://grafana.internal/d/api", "reference"),
        ("See pr 8423 on github.com/foo/bar", "reference"),
        ("Prefers terse answers in German", "user"),  # fallback
        ("Senior data scientist", "user"),
    ],
)
def test_parse_memory_type_auto_classifier(text, expected):
    t, body = ms.parse_memory_type(text)
    assert t == expected, (text, t)
    # Body must not be the type prefix or empty
    assert body == text


def test_parse_memory_type_explicit_prefix_beats_heuristic():
    """An explicit prefix wins even when the heuristic would pick something
    else — the user's intent is sacred."""
    t, body = ms.parse_memory_type("user: don't mock the DB")
    assert t == "user"
    assert body == "don't mock the DB"


# ---------------------------------------------------------------------------
# Wiki-style [[name]] cross-link resolution
# ---------------------------------------------------------------------------

def _write_memory(memory_dir: Path, kind: str, name: str, desc: str, body: str = "body") -> Path:
    memory_dir.mkdir(parents=True, exist_ok=True)
    p = memory_dir / f"{kind}_{name}.md"
    p.write_text(
        f"---\nname: {name}\ndescription: {desc}\nmetadata:\n  type: {kind}\n---\n\n{body}\n",
        encoding="utf-8",
    )
    return p


def test_resolve_wikilinks_resolved(tmp_path):
    _write_memory(tmp_path, "feedback", "no-mocks", "Integration tests must hit a real DB")
    text = "We follow [[no-mocks]] when writing tests."
    out = ms.resolve_wikilinks(text, tmp_path)
    assert "[no-mocks](feedback_no-mocks.md)" in out
    assert "Integration tests must hit a real DB" in out


def test_resolve_wikilinks_unresolved_marked(tmp_path):
    out = ms.resolve_wikilinks("Track [[missing-memory]] later.", tmp_path)
    assert "[[missing-memory]] (not yet written)" in out


def test_resolve_wikilinks_no_links_passthrough(tmp_path):
    text = "Just prose, no links here."
    assert ms.resolve_wikilinks(text, tmp_path) == text


def test_resolve_wikilinks_multiple_in_one_line(tmp_path):
    _write_memory(tmp_path, "user", "prefers-german", "Terse, exact, numeric")
    _write_memory(tmp_path, "project", "ml-refiner", "ML refiner Q2 work")
    out = ms.resolve_wikilinks(
        "Linked: [[prefers-german]] and [[ml-refiner]] and [[unknown]].",
        tmp_path,
    )
    assert "[prefers-german]" in out
    assert "[ml-refiner]" in out
    assert "[[unknown]] (not yet written)" in out


def test_resolve_wikilinks_priority_feedback_first(tmp_path):
    """If a name accidentally exists under multiple type prefixes,
    feedback wins (it's the most-edited / most-correctional layer)."""
    _write_memory(tmp_path, "user", "shared-name", "user version")
    _write_memory(tmp_path, "feedback", "shared-name", "feedback version")
    out = ms.resolve_wikilinks("[[shared-name]]", tmp_path)
    assert "feedback_shared-name.md" in out
    assert "feedback version" in out


# ---------------------------------------------------------------------------
# Stale-reference verification (/memories verify backend)
# ---------------------------------------------------------------------------

def test_find_stale_references_detects_missing_paths(tmp_path):
    real = tmp_path / "subdir" / "existing.py"
    real.parent.mkdir(parents=True)
    real.write_text("# real file\n", encoding="utf-8")

    text = (
        "Look at `subdir/existing.py:42` for the impl, and the broken "
        "ref `subdir/gone.py` should be flagged."
    )
    stale = ms.find_stale_references(text, tmp_path)
    assert any("gone.py" in s for s in stale), stale
    assert not any("existing.py" in s for s in stale), stale


def test_find_stale_references_ignores_urls(tmp_path):
    text = "Docs at https://example.com/docs/api.html and grafana.internal/d/x"
    stale = ms.find_stale_references(text, tmp_path)
    # URLs and dotted hostnames must not be flagged as paths
    assert not any("example.com" in s for s in stale), stale


def test_find_stale_references_ignores_prose(tmp_path):
    text = "Just words about doing, working, helping. No paths here."
    assert ms.find_stale_references(text, tmp_path) == []


def test_verify_typed_memories_round_trip(tmp_path, monkeypatch):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    repo = tmp_path / "myrepo"
    repo.mkdir()
    # File that DOES exist in repo
    (repo / "real_module.py").write_text("# real", encoding="utf-8")
    # Memory referencing one existing + one missing path
    ms.save_typed_memory(
        "feedback: see real_module.py:10 and gone_module.py for fix",
        repo_root=repo,
    )
    issues = ms.verify_typed_memories(repo)
    assert issues
    assert any("gone_module.py" in r for r in issues[0]["stale_refs"])
    assert not any("real_module.py" in r for r in issues[0]["stale_refs"])
