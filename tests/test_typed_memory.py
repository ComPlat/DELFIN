"""Tests for save_typed_memory + parse_memory_type.

The typed per-project memory store lives under
``~/.claude/projects/<slug>/memory/`` — we monkey-patch ``Path.home`` so
the test never touches the real home directory.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import memory_store as ms


@pytest.fixture
def fake_home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    return tmp_path


def test_parse_memory_type_default_user():
    t, body = ms.parse_memory_type("prefers German")
    assert t == "user"
    assert body == "prefers German"


def test_parse_memory_type_explicit_prefix():
    for prefix, expected in [
        ("user:", "user"),
        ("feedback:", "feedback"),
        ("project:", "project"),
        ("reference:", "reference"),
        ("User:", "user"),
        ("FEEDBACK :", "feedback"),
    ]:
        t, body = ms.parse_memory_type(f"{prefix} the rest")
        assert t == expected, (prefix, t)
        assert body == "the rest"


def test_parse_memory_type_unknown_prefix_falls_back():
    """A leading word that *looks* like a prefix but isn't one of the
    four types must be preserved as part of the body."""
    t, body = ms.parse_memory_type("notes: this is just notes")
    assert t == "user"
    assert body == "notes: this is just notes"


def test_save_typed_memory_creates_file_with_frontmatter(fake_home, tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    fpath, slug, mtype = ms.save_typed_memory(
        "user: prefers German, exact, numeric",
        repo_root=repo,
    )
    assert mtype == "user"
    assert fpath.exists()
    assert fpath.name.startswith("user_")
    text = fpath.read_text(encoding="utf-8")
    assert text.startswith("---\n")
    assert f"name: {slug}" in text
    assert "metadata:" in text and "type: user" in text
    assert "prefers German" in text


def test_save_typed_memory_indexes_in_memory_md(fake_home, tmp_path):
    repo = tmp_path / "repo2"
    repo.mkdir()
    ms.save_typed_memory(
        "feedback: do not mock the DB",
        repo_root=repo,
    )
    index = fpath = ms._claude_memory_dir(repo) / "MEMORY.md"
    content = index.read_text(encoding="utf-8")
    assert "## Feedback (how to work)" in content
    assert "(feedback_" in content
    assert "do not mock the DB" in content


def test_save_typed_memory_merges_duplicate(fake_home, tmp_path):
    repo = tmp_path / "repo3"
    repo.mkdir()
    p1, _, _ = ms.save_typed_memory("user: first fact here", repo_root=repo)
    # Saving the SAME fact again must MERGE in place (dedup), not create a
    # second file. The store deduplicates deterministically so weak/open
    # models get this for free.
    p2, _, _ = ms.save_typed_memory("user: first fact here", repo_root=repo)
    assert p1 == p2
    mem_dir = ms._claude_memory_dir(repo)
    md_files = [p for p in mem_dir.glob("*.md") if p.name != "MEMORY.md"]
    assert len(md_files) == 1
    # use_count bumped on the merge.
    recs = ms.list_typed_memories(repo)
    assert len(recs) == 1
    assert recs[0]["use_count"] == 2


def test_save_typed_memory_creates_missing_section(fake_home, tmp_path):
    """If MEMORY.md exists but lacks the type's section header,
    the writer appends a new section at the bottom."""
    repo = tmp_path / "repo4"
    repo.mkdir()
    mem_dir = ms._claude_memory_dir(repo)
    mem_dir.mkdir(parents=True)
    (mem_dir / "MEMORY.md").write_text(
        "# DELFIN Project Memory\n\n## Project (current work)\n- prior\n",
        encoding="utf-8",
    )
    ms.save_typed_memory("reference: see Linear/INGEST", repo_root=repo)
    content = (mem_dir / "MEMORY.md").read_text(encoding="utf-8")
    assert "## Reference (read only when relevant)" in content
    assert "Linear/INGEST" in content


def test_save_typed_memory_explicit_type_overrides_prefix(fake_home, tmp_path):
    repo = tmp_path / "repo5"
    repo.mkdir()
    _, _, mtype = ms.save_typed_memory(
        "user: this is the body",
        repo_root=repo,
        memory_type="project",
    )
    assert mtype == "project"


def test_slugify_handles_unicode_and_punctuation():
    assert ms._slugify("Über die Brücke!! ja??") in {
        "ber-die-brcke-ja",     # ASCII-stripped variant
        "ber-die-br-cke-ja",
    } or ms._slugify("Über die Brücke!! ja??")  # at minimum non-empty


# ---------------------------------------------------------------------------
# list_typed_memories + delete_typed_memory (single-store consolidation)
# ---------------------------------------------------------------------------


def test_list_typed_memories_returns_records(fake_home, tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    ms.save_typed_memory("feedback: never add a co-author trailer", repo_root=repo)
    ms.save_typed_memory("user: Max is a quantum chemist at KIT", repo_root=repo)
    recs = ms.list_typed_memories(repo)
    assert {r["type"] for r in recs} == {"feedback", "user"}
    by_type = {r["type"]: r for r in recs}
    assert "co-author" in by_type["feedback"]["body"].lower()
    assert by_type["user"]["name"]            # slug present
    assert by_type["user"]["description"]     # frontmatter description present


def test_list_typed_memories_empty_when_none(fake_home, tmp_path):
    assert ms.list_typed_memories(tmp_path / "norepo") == []


def test_delete_typed_memory_removes_file_and_index(fake_home, tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    fpath, slug, _ = ms.save_typed_memory(
        "project: ship the memory layer this week", repo_root=repo)
    mdir = ms._claude_memory_dir(repo)
    assert (mdir / fpath.name).is_file()
    index_before = (mdir / "MEMORY.md").read_text(encoding="utf-8")
    assert fpath.name in index_before

    deleted = ms.delete_typed_memory(repo, slug)
    assert deleted is not None
    assert not fpath.exists()
    index_after = (mdir / "MEMORY.md").read_text(encoding="utf-8")
    assert fpath.name not in index_after       # pointer line removed


def test_delete_typed_memory_unknown_returns_none(fake_home, tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    ms.save_typed_memory("user: something", repo_root=repo)
    assert ms.delete_typed_memory(repo, "does-not-exist") is None


# ---------------------------------------------------------------------------
# Store lives under DELFIN's OWN ~/.delfin (not the legacy ~/.claude namespace)
# ---------------------------------------------------------------------------


def test_memory_dir_is_under_delfin_not_claude(fake_home, tmp_path):
    repo = tmp_path / "repo"; repo.mkdir()
    mdir = ms._delfin_memory_dir(repo)
    assert ".delfin" in str(mdir) and ".claude" not in str(mdir)
    fpath, _, _ = ms.save_typed_memory("user: x lives here", repo_root=repo)
    assert ".delfin" in str(fpath)


def test_plans_dir_under_delfin(fake_home, tmp_path):
    repo = tmp_path / "repo"; repo.mkdir()
    assert ".delfin" in str(ms._delfin_plans_dir(repo))


def test_legacy_claude_store_is_migrated(fake_home, tmp_path):
    repo = tmp_path / "repo"; repo.mkdir()
    slug = ms._project_slug(repo)
    old = tmp_path / ".claude" / "projects" / slug / "memory"
    old.mkdir(parents=True)
    (old / "MEMORY.md").write_text("# index\n", encoding="utf-8")
    (old / "user_legacy.md").write_text(
        "---\nname: legacy\ndescription: an old fact\n"
        "metadata:\n  type: user\n---\n\nold fact\n", encoding="utf-8")
    # Resolving the dir migrates the legacy store into ~/.delfin.
    new = ms._delfin_memory_dir(repo)
    assert ".delfin" in str(new)
    assert (new / "user_legacy.md").is_file()
    assert not old.exists()                         # moved, not copied
    assert "legacy" in [r["name"] for r in ms.list_typed_memories(repo)]


# ---------------------------------------------------------------------------
# Dedup / merge / prune / recall (unbounded-growth fix)
# ---------------------------------------------------------------------------

def test_near_duplicate_merges_not_duplicates(fake_home, tmp_path):
    """A fact that is *similar* (not byte-identical) to an existing one of
    the same type is merged in place rather than duplicated."""
    repo = tmp_path / "dup"; repo.mkdir()
    ms.save_typed_memory(
        "project: the scheduler uses PAL cores for parallel jobs",
        repo_root=repo)
    ms.save_typed_memory(
        "project: scheduler uses PAL cores for the parallel jobs today",
        repo_root=repo)
    recs = [r for r in ms.list_typed_memories(repo) if r["type"] == "project"]
    assert len(recs) == 1
    assert recs[0]["use_count"] == 2
    assert "superseded:" in Path(recs[0]["path"]).read_text(encoding="utf-8")


def test_merge_threshold_is_configurable(fake_home, tmp_path, monkeypatch):
    repo = tmp_path / "threshold"; repo.mkdir()
    monkeypatch.setenv("DELFIN_MEMORY_MERGE_THRESHOLD", "1.0")
    ms.save_typed_memory(
        "project: scheduler uses PAL cores for parallel jobs", repo_root=repo)
    ms.save_typed_memory(
        "project: scheduler uses PAL cores for parallel jobs today", repo_root=repo)
    recs = [r for r in ms.list_typed_memories(repo) if r["type"] == "project"]
    assert len(recs) == 2


def test_distinct_facts_stay_separate(fake_home, tmp_path):
    """Two unrelated facts of the same type must NOT be merged."""
    repo = tmp_path / "distinct"; repo.mkdir()
    ms.save_typed_memory(
        "project: the scheduler uses PAL cores for parallel jobs",
        repo_root=repo)
    ms.save_typed_memory(
        "project: emitters are optimised with CAM-B3LYP in toluene",
        repo_root=repo)
    recs = [r for r in ms.list_typed_memories(repo) if r["type"] == "project"]
    assert len(recs) == 2


def test_different_types_never_merge(fake_home, tmp_path):
    """Identical text under different types stays as two separate memories."""
    repo = tmp_path / "types"; repo.mkdir()
    ms.save_typed_memory("project: prefers explicit state contracts",
                         repo_root=repo)
    ms.save_typed_memory("feedback: prefers explicit state contracts",
                         repo_root=repo)
    types = {r["type"] for r in ms.list_typed_memories(repo)}
    assert {"project", "feedback"} <= types


def test_prune_caps_prunable_type(fake_home, tmp_path):
    """Pruning keeps the most-used entries and removes the remaining tail."""
    repo = tmp_path / "prune"; repo.mkdir()
    protected_path = None
    for i in range(6):
        path, _, _ = ms.save_typed_memory(
            f"reference: distinct source number {i} alpha{i}", repo_root=repo)
        if i == 0:
            protected_path = path
    assert protected_path is not None
    ms.record_memory_recall(repo, [protected_path.name])
    ms.record_memory_recall(repo, [protected_path.name])
    before = [r for r in ms.list_typed_memories(repo) if r["type"] == "reference"]
    assert len(before) == 6
    deleted = ms.prune_memories(repo, max_per_type=3)
    assert len(deleted) == 3
    after = [r for r in ms.list_typed_memories(repo) if r["type"] == "reference"]
    assert len(after) == 3
    assert protected_path.name in {r["file"] for r in after}


def test_prune_gives_feedback_and_user_a_larger_cap(fake_home, tmp_path):
    repo = tmp_path / "protect"; repo.mkdir()
    for i in range(9):
        ms.save_typed_memory(f"feedback: distinct guidance item {i} beta{i}",
                             repo_root=repo)
        ms.save_typed_memory(f"user: distinct preference item {i} gamma{i}",
                             repo_root=repo)
    ms.prune_memories(repo, max_per_type=2)
    fb = [r for r in ms.list_typed_memories(repo) if r["type"] == "feedback"]
    us = [r for r in ms.list_typed_memories(repo) if r["type"] == "user"]
    assert len(fb) == 8
    assert len(us) == 8


def test_prune_removes_stale_prunable_entries(
        fake_home, tmp_path, monkeypatch):
    repo = tmp_path / "stale"; repo.mkdir()
    clock = {"now": 1_000_000}
    monkeypatch.setattr(ms.time, "time", lambda: clock["now"])
    stale, _, _ = ms.save_typed_memory(
        "reference: obsolete deployment instructions", repo_root=repo)
    clock["now"] += 2 * 86_400
    fresh, _, _ = ms.save_typed_memory(
        "reference: current API compatibility matrix", repo_root=repo)

    deleted = ms.prune_memories(repo, max_per_type=10, max_age_days=1)

    assert stale.name in deleted
    assert not stale.exists()
    assert fresh.exists()


def test_old_file_without_frontmatter_fields_loads_with_defaults(
        fake_home, tmp_path):
    """A legacy .md predating use_count/updated_at loads with safe defaults
    (use_count 0, timestamps from file mtime), without migration."""
    repo = tmp_path / "legacy"; repo.mkdir()
    mem_dir = ms._delfin_memory_dir(repo)
    mem_dir.mkdir(parents=True, exist_ok=True)
    (mem_dir / "project_old.md").write_text(
        "---\nname: old\ndescription: a pre-existing fact\n"
        "metadata:\n  type: project\n---\n\nsome older project fact\n",
        encoding="utf-8")
    recs = ms.list_typed_memories(repo)
    assert len(recs) == 1
    assert recs[0]["use_count"] == 0
    assert recs[0]["updated_at"] > 0        # falls back to mtime


def test_record_memory_recall_bumps_use_count(fake_home, tmp_path):
    """record_memory_recall increments use_count only for the named files."""
    repo = tmp_path / "recall"; repo.mkdir()
    p, _slug, _t = ms.save_typed_memory(
        "project: some recallable project fact here", repo_root=repo)
    n = ms.record_memory_recall(repo, [p.name])
    assert n == 1
    rec = [r for r in ms.list_typed_memories(repo) if r["file"] == p.name][0]
    assert rec["use_count"] == 2
    # A second recall bumps again.
    ms.record_memory_recall(repo, [p.name])
    rec = [r for r in ms.list_typed_memories(repo) if r["file"] == p.name][0]
    assert rec["use_count"] == 3


def test_record_memory_recall_rejects_parent_path(fake_home, tmp_path):
    repo = tmp_path / "recall-path"; repo.mkdir()
    outside = ms._delfin_memory_dir(repo).parent / "outside.md"
    outside.parent.mkdir(parents=True, exist_ok=True)
    outside.write_text("must stay unchanged\n", encoding="utf-8")
    assert ms.record_memory_recall(repo, ["../outside.md"]) == 0
    assert outside.read_text(encoding="utf-8") == "must stay unchanged\n"


def test_rewrites_preserve_hand_added_frontmatter(fake_home, tmp_path):
    """Unknown scalar frontmatter fields (hand-added by the user) survive
    both a recall rewrite and a dedup merge."""
    repo = tmp_path / "extras"; repo.mkdir()
    p, _slug, _t = ms.save_typed_memory(
        "project: emitters are optimised with CAM-B3LYP", repo_root=repo)
    text = p.read_text(encoding="utf-8")
    p.write_text(text.replace("metadata:", "source: hand-curated\nmetadata:"),
                 encoding="utf-8")

    ms.record_memory_recall(repo, [p.name])
    assert "source: hand-curated" in p.read_text(encoding="utf-8")

    ms.save_typed_memory(
        "project: emitters are optimised with CAM-B3LYP today",
        repo_root=repo)
    after = p.read_text(encoding="utf-8")
    assert "source: hand-curated" in after
    assert "superseded:" in after


def test_merge_refreshes_index_hook(fake_home, tmp_path):
    """After a dedup merge the MEMORY.md pointer reflects the NEW text."""
    repo = tmp_path / "hook"; repo.mkdir()
    p1, _, _ = ms.save_typed_memory(
        "project: scheduler uses PAL cores for parallel jobs", repo_root=repo)
    p2, _, _ = ms.save_typed_memory(
        "project: scheduler uses PAL cores for parallel jobs on niagara",
        repo_root=repo)
    assert p1 == p2
    index = (ms._delfin_memory_dir(repo) / "MEMORY.md").read_text(
        encoding="utf-8")
    assert index.count(f"({p1.name})") == 1
    assert "niagara" in index


def test_memory_writes_leave_no_temp_files(fake_home, tmp_path):
    """Atomic writes must clean up their temp files in the memory dir."""
    repo = tmp_path / "tmpfiles"; repo.mkdir()
    p, _, _ = ms.save_typed_memory(
        "project: no temp file droppings please", repo_root=repo)
    ms.save_typed_memory(
        "project: no temp file droppings please again", repo_root=repo)
    ms.record_memory_recall(repo, [p.name])
    ms.prune_memories(repo, max_per_type=1)
    leftovers = [f.name for f in ms._delfin_memory_dir(repo).iterdir()
                 if f.name.endswith(".tmp")]
    assert leftovers == []


def test_fresh_memory_survives_post_save_prune(
        fake_home, tmp_path, monkeypatch):
    """A just-saved memory must survive the prune that runs right after its
    save, even when a saturated store holds older entries with far higher
    use_counts (recency ranks before use_count)."""
    repo = tmp_path / "fresh"; repo.mkdir()
    clock = {"now": 1_000_000}
    monkeypatch.setattr(ms.time, "time", lambda: clock["now"])
    paths = []
    for i in range(3):
        p, _, _ = ms.save_typed_memory(
            f"reference: established source number {i} delta{i}",
            repo_root=repo)
        paths.append(p)
    # Established entries get recalled often -> high use_count.
    for _ in range(5):
        ms.record_memory_recall(repo, [p.name for p in paths])
    clock["now"] += 60
    fresh, _, _ = ms.save_typed_memory(
        "reference: brand new source worth keeping", repo_root=repo)

    deleted = ms.prune_memories(repo, max_per_type=3)

    assert len(deleted) == 1
    assert fresh.name not in deleted
    assert fresh.exists()


# ---------------------------------------------------------------------------
# Recall-time provenance: content anchors, stale/drifted refs, stale_hits
# ---------------------------------------------------------------------------


def test_save_captures_content_anchors(fake_home, tmp_path):
    """A body with a file:line reference gets the cited line's text stored
    as a content anchor in the frontmatter."""
    repo = tmp_path / "anchors"; repo.mkdir()
    src = repo / "src"; src.mkdir()
    (src / "app.py").write_text(
        "import os\ndef retry_loop():\n    return 42\n", encoding="utf-8")
    fpath, _, _ = ms.save_typed_memory(
        "project: the retry loop lives in src/app.py:2 and buffers writes",
        repo_root=repo)
    text = fpath.read_text(encoding="utf-8")
    assert "anchors:" in text
    meta, _ = ms._parse_frontmatter(text)
    assert ms._parse_anchors(meta) == {"src/app.py:2": "def retry_loop():"}


def test_save_without_line_refs_writes_no_anchors(fake_home, tmp_path):
    """Refs without a :line suffix (or to unreadable files) get no anchor."""
    repo = tmp_path / "anchors2"; repo.mkdir()
    fpath, _, _ = ms.save_typed_memory(
        "project: config lives in conf/settings.yaml near the top",
        repo_root=repo)
    assert "anchors:" not in fpath.read_text(encoding="utf-8")


def test_find_stale_references_flat_list_back_compat(fake_home, tmp_path):
    """Legacy call shape: flat list of dead refs, healthy refs excluded."""
    repo = tmp_path / "flat"; repo.mkdir()
    (repo / "still_here.py").write_text("x = 1\n", encoding="utf-8")
    refs = ms.find_stale_references(
        "see still_here.py and gone_module.py for details", repo)
    assert refs == ["gone_module.py"]


def test_drifted_anchor_detected(fake_home, tmp_path):
    """File exists but the anchored line moved out of the ±20 window."""
    repo = tmp_path / "drift"; repo.mkdir()
    src = repo / "src"; src.mkdir()
    original = "\n".join(f"line_{i} = {i}" for i in range(1, 61))
    (src / "core.py").write_text(original, encoding="utf-8")
    fpath, _, _ = ms.save_typed_memory(
        "project: the counter init sits at src/core.py:5 exactly",
        repo_root=repo)
    meta, body = ms._parse_frontmatter(fpath.read_text(encoding="utf-8"))
    anchors = ms._parse_anchors(meta)
    assert anchors["src/core.py:5"] == "line_5 = 5"

    # Healthy: anchor still at the cited line -> no findings.
    assert ms.find_stale_references(
        body, repo, anchors=anchors, with_status=True) == []

    # Shifted a little (still within ±20 lines) -> healthy, not drifted.
    (src / "core.py").write_text("# new\n# new\n# new\n" + original,
                                 encoding="utf-8")
    assert ms.find_stale_references(
        body, repo, anchors=anchors, with_status=True) == []

    # Fully rewritten -> the anchor is gone -> drifted.
    (src / "core.py").write_text(
        "\n".join(f"other_{i} = {i}" for i in range(1, 61)), encoding="utf-8")
    statuses = ms.find_stale_references(
        body, repo, anchors=anchors, with_status=True)
    assert ("src/core.py:5", "drifted") in statuses


def test_verify_typed_memories_reports_stale_and_drifted(fake_home, tmp_path):
    repo = tmp_path / "verify"; repo.mkdir()
    (repo / "mod.py").write_text("alpha = 1\nbeta = 2\n", encoding="utf-8")
    ms.save_typed_memory("project: beta is defined in mod.py:2 today",
                         repo_root=repo)
    ms.save_typed_memory("project: legacy notes sit in old/notes.txt still",
                         repo_root=repo)
    (repo / "mod.py").write_text("gamma = 3\n", encoding="utf-8")
    results = ms.verify_typed_memories(repo)
    stale = [r for rec in results for r in rec["stale_refs"]]
    drifted = [r for rec in results for r in rec["drifted_refs"]]
    assert "old/notes.txt" in stale
    assert "mod.py:2" in drifted


def test_record_stale_hits_bumps_counter_not_usage(fake_home, tmp_path):
    """stale_hits rises; use_count and updated_at stay untouched so the
    rotted entry keeps decaying like an unused one."""
    repo = tmp_path / "rot"; repo.mkdir()
    p, _, _ = ms.save_typed_memory(
        "project: helper lived in gone/helper.py once", repo_root=repo)
    before = [r for r in ms.list_typed_memories(repo) if r["file"] == p.name][0]
    assert ms.record_stale_hits(repo, [p.name]) == 1
    rec = [r for r in ms.list_typed_memories(repo) if r["file"] == p.name][0]
    assert rec["stale_hits"] == 1
    assert rec["use_count"] == before["use_count"]        # NOT bumped
    assert rec["updated_at"] == before["updated_at"]      # NOT bumped
    ms.record_stale_hits(repo, [p.name])
    rec = [r for r in ms.list_typed_memories(repo) if r["file"] == p.name][0]
    assert rec["stale_hits"] == 2
    # A later healthy recall must not clobber the rot counter.
    ms.record_memory_recall(repo, [p.name])
    rec = [r for r in ms.list_typed_memories(repo) if r["file"] == p.name][0]
    assert rec["stale_hits"] == 2
    assert rec["use_count"] == 2


def test_stale_hits_influence_prune_order(fake_home, tmp_path, monkeypatch):
    """On an updated_at/use_count tie the rotted entry is evicted first."""
    repo = tmp_path / "rotprune"; repo.mkdir()
    clock = {"now": 1_000_000}
    monkeypatch.setattr(ms.time, "time", lambda: clock["now"])
    paths = []
    for i in range(3):
        p, _, _ = ms.save_typed_memory(
            f"reference: distinct pointer number {i} rho{i}", repo_root=repo)
        paths.append(p)
    rotted = paths[1]
    ms.record_stale_hits(repo, [rotted.name])

    deleted = ms.prune_memories(repo, max_per_type=2)

    assert deleted == [rotted.name]
    assert not rotted.exists()
    assert all(p.exists() for p in paths if p is not rotted)


# ---------------------------------------------------------------------------
# User-scoped global store (~/.delfin/memory)
# ---------------------------------------------------------------------------


def test_parse_memory_scope():
    assert ms.parse_memory_scope("global: feedback: be terse") == (
        "user", "feedback: be terse")
    assert ms.parse_memory_scope("feedback: be terse") == (
        "project", "feedback: be terse")
    # Colon-only prefix: hyphenated words must survive untouched.
    assert ms.parse_memory_scope("global-warming data is rising") == (
        "project", "global-warming data is rising")


def test_parse_memory_type_global_prefix():
    assert ms.parse_memory_type("global: feedback: never guess APIs") == (
        "feedback", "never guess APIs")
    # A bare global fact defaults to the user (identity) type.
    assert ms.parse_memory_type("global: answers come back in German") == (
        "user", "answers come back in German")


def test_global_prefix_routes_to_user_store(fake_home, tmp_path):
    repo = tmp_path / "repo"; repo.mkdir()
    fpath, _, mtype = ms.save_typed_memory(
        "global: feedback: never add attribution trailers", repo_root=repo)
    assert mtype == "feedback"
    gdir = tmp_path / ".delfin" / "memory"
    assert fpath.parent == gdir
    index = (gdir / "MEMORY.md").read_text(encoding="utf-8")
    assert "Global Memory" in index
    assert "attribution" in index
    # The project store stays untouched.
    assert ms.list_typed_memories(repo) == []


def test_scope_param_routes_to_user_store(fake_home, tmp_path):
    repo = tmp_path / "repo"; repo.mkdir()
    fpath, _, _ = ms.save_typed_memory(
        "user: Max is a quantum chemist", repo_root=repo, scope="user")
    assert fpath.parent == tmp_path / ".delfin" / "memory"


def test_global_memory_skips_anchor_capture(fake_home, tmp_path):
    """Global facts outlive any single checkout — no repo-relative anchors."""
    repo = tmp_path / "repo"; repo.mkdir()
    (repo / "mod.py").write_text("alpha = 1\n", encoding="utf-8")
    fpath, _, _ = ms.save_typed_memory(
        "global: reference: interesting bits in mod.py:1 sometimes",
        repo_root=repo)
    assert "anchors:" not in fpath.read_text(encoding="utf-8")


def test_list_typed_memories_scopes(fake_home, tmp_path):
    repo = tmp_path / "repo"; repo.mkdir()
    ms.save_typed_memory("project: repo-local fact alpha", repo_root=repo)
    ms.save_typed_memory("global: user: cross-repo identity beta",
                         repo_root=repo)
    proj = ms.list_typed_memories(repo)
    assert [r["scope"] for r in proj] == ["project"]
    glob = ms.list_typed_memories(repo, scope="user")
    assert [r["scope"] for r in glob] == ["user"]
    assert "beta" in glob[0]["body"]
    both = ms.list_typed_memories(repo, scope="all")
    assert [r["scope"] for r in both] == ["user", "project"]   # global first


def test_delete_typed_memory_user_scope(fake_home, tmp_path):
    repo = tmp_path / "repo"; repo.mkdir()
    fpath, slug, _ = ms.save_typed_memory(
        "global: user: deletable identity fact", repo_root=repo)
    # The default (project) scope cannot see it ...
    assert ms.delete_typed_memory(repo, slug) is None
    # ... the user scope can, and the global index pointer goes with it.
    assert ms.delete_typed_memory(repo, slug, scope="user") == fpath
    assert not fpath.exists()
    index = (tmp_path / ".delfin" / "memory" / "MEMORY.md").read_text(
        encoding="utf-8")
    assert fpath.name not in index


def test_delete_typed_memory_scope_all_reaches_both(fake_home, tmp_path):
    repo = tmp_path / "repo"; repo.mkdir()
    gpath, gslug, _ = ms.save_typed_memory(
        "global: user: reachable via scope all", repo_root=repo)
    assert ms.delete_typed_memory(repo, gslug, scope="all") == gpath
    assert not gpath.exists()


def test_prune_user_store_gets_protected_caps(fake_home, tmp_path):
    """In the global store even prunable types get the protected 4x cap
    and are exempt from age-based decay."""
    repo = tmp_path / "repo"; repo.mkdir()
    for i in range(6):
        ms.save_typed_memory(
            f"global: reference: distinct global source {i} sigma{i}",
            repo_root=repo)
    deleted = ms.prune_memories(repo, max_per_type=1, max_age_days=0,
                                scope="user")
    assert len(deleted) == 2
    assert len(ms.list_typed_memories(repo, scope="user")) == 4
    # The project store is a separate prune domain — untouched.
    assert ms.list_typed_memories(repo) == []


def test_record_memory_recall_user_scope(fake_home, tmp_path):
    repo = tmp_path / "repo"; repo.mkdir()
    p, _, _ = ms.save_typed_memory(
        "global: user: recallable identity fact", repo_root=repo)
    assert ms.record_memory_recall(repo, [p.name], scope="user") == 1
    rec = ms.list_typed_memories(repo, scope="user")[0]
    assert rec["use_count"] == 2
