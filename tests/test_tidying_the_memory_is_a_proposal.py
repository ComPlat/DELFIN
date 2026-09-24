"""Tidying the memory store proposes; it does not act until asked.

Input: a typed-memory store. Output: a proposal — pairs to merge, entries
to retire — with counts before and after. Nothing on disk changes until
``apply`` is called, and ``apply`` never deletes.

Why a proposal rather than a sweep. A store is what the agent learned; an
unsupervised pass over it can drop a fact and nobody notices for weeks.
The opposite failure was measured on 2026-09-24: 358 of 371 notes had
been invisible for months because their store was keyed to paths that no
longer existed, and that went unnoticed too. Both directions are silent,
so the change is made visible instead of automatic.

What it reuses, and does not reimplement:
  - ``_jaccard`` and ``_merge_similarity_threshold`` decide near-duplicates
  - ``list_typed_memories`` supplies records with usage and decay metadata
  - ``_AGENT_MEMORY_MAX_AGE_DAYS`` defines disuse

What it adds: merging applies today only when a memory is WRITTEN, so
look-alikes already in the store stay there; and ``prune_memories``
deletes on its own schedule without showing anyone what it took.

Retirement moves the file to ``<store>/retired/`` with the date. The
store shrinks, the text stays readable, and a wrong call costs a move
rather than a fact.
"""

from __future__ import annotations

import time

import pytest

from delfin.agent import memory_tidy as T


NOW = int(time.time())
DAY = 86400


def _write(store, name, mtype, body, *, age_days=0, source="agent",
           use_count=3):
    store.mkdir(parents=True, exist_ok=True)
    stamp = NOW - age_days * DAY
    path = store / f"{mtype}_{name}.md"
    path.write_text(
        "---\n"
        f"name: {name}\n"
        f"description: {name} description\n"
        "metadata:\n"
        f"  type: {mtype}\n"
        f"  source: {source}\n"
        f"  created_at: {stamp}\n"
        f"  updated_at: {stamp}\n"
        f"  use_count: {use_count}\n"
        "---\n\n"
        f"{body}\n",
        encoding="utf-8")
    return path


@pytest.fixture()
def store(tmp_path):
    return tmp_path / "memory"


# -- merging look-alikes ----------------------------------------------------

def test_two_near_duplicates_are_proposed_for_merging(store):
    _write(store, "gate-recipe", "feedback",
           "Run the full gate with TMPDIR set before pushing a branch.")
    _write(store, "how-to-gate", "feedback",
           "Before pushing a branch run the full gate, TMPDIR must be set.")
    p = T.propose(store)
    assert len(p.merges) == 1, p.merges
    pair = p.merges[0]
    assert {pair.keep.name, pair.drop.name} == {"gate-recipe", "how-to-gate"}
    assert pair.similarity >= 0.5


def test_unrelated_notes_are_left_alone(store):
    _write(store, "gate-recipe", "feedback", "Run the gate before pushing.")
    _write(store, "kit-latency", "reference",
           "The endpoint answers a one-word probe in about six seconds.")
    assert T.propose(store).merges == []


def test_notes_of_different_types_are_never_merged(store):
    """Same words, different kind of fact: the type is the user's
    classification and merging across it loses that."""
    text = "The gate must run with TMPDIR set before a push."
    _write(store, "a", "feedback", text)
    _write(store, "b", "project", text)
    assert T.propose(store).merges == []


def test_the_older_note_is_kept_and_the_newer_folded_into_it(store):
    """The name that other memories link to is the one that stays."""
    _write(store, "older", "feedback", "Run the gate before pushing.",
           age_days=30)
    _write(store, "newer", "feedback", "Before pushing, run the gate.",
           age_days=1)
    pair = T.propose(store).merges[0]
    assert pair.keep.name == "older"
    assert pair.drop.name == "newer"


# -- retiring what nobody recalls -------------------------------------------

def test_a_model_written_note_nobody_recalled_is_proposed_for_retirement(store):
    _write(store, "stale", "project", "Something learned long ago.",
           age_days=200, source="agent")
    p = T.propose(store)
    assert [r.name for r in p.retire] == ["stale"]


def test_the_users_own_notes_are_never_retired(store):
    """Disuse is not a reason to drop what the user wrote down."""
    _write(store, "users-rule", "feedback", "Always commit as the repo id.",
           age_days=500, source="user")
    assert T.propose(store).retire == []


def test_a_recalled_note_survives_however_old(store):
    """Recall bumps updated_at, so this is disuse and not age."""
    _write(store, "fresh-recall", "project", "Still in use.",
           age_days=0, source="agent", use_count=40)
    assert T.propose(store).retire == []


# -- the report -------------------------------------------------------------

def test_the_proposal_counts_before_and_after(store):
    _write(store, "a", "feedback", "Run the gate before pushing.")
    _write(store, "b", "feedback", "Before pushing, run the gate.")
    _write(store, "old", "project", "Ancient.", age_days=200)
    p = T.propose(store)
    assert p.before == 3
    assert p.after == 1          # one merged away, one retired
    assert "3" in p.render() and "1" in p.render()


def test_an_empty_store_proposes_nothing(store):
    store.mkdir(parents=True)
    p = T.propose(store)
    assert p.merges == [] and p.retire == [] and p.before == 0


def test_a_missing_store_is_not_an_error(tmp_path):
    p = T.propose(tmp_path / "nothing-here")
    assert p.before == 0 and p.merges == [] and p.retire == []


# -- nothing happens until it is asked for ----------------------------------

def test_proposing_changes_nothing_on_disk(store):
    a = _write(store, "a", "feedback", "Run the gate before pushing.")
    b = _write(store, "b", "feedback", "Before pushing, run the gate.")
    before = {p: p.read_text(encoding="utf-8") for p in (a, b)}
    T.propose(store)
    for path, text in before.items():
        assert path.exists() and path.read_text(encoding="utf-8") == text


def test_applying_folds_the_pair_into_one_file(store):
    _write(store, "older", "feedback", "Run the gate before pushing.",
           age_days=30)
    _write(store, "newer", "feedback", "Before pushing, run the gate.",
           age_days=1)
    p = T.propose(store)
    T.apply(p)
    left = sorted(f.name for f in store.glob("*.md"))
    assert left == ["feedback_older.md"], left
    kept = (store / "feedback_older.md").read_text(encoding="utf-8")
    assert "Run the gate before pushing." in kept
    assert "Before pushing, run the gate." in kept, (
        "the folded text has to survive, or merging loses a fact")


def test_retiring_moves_and_never_deletes(store):
    _write(store, "stale", "project", "Something learned long ago.",
           age_days=200, source="agent")
    p = T.propose(store)
    T.apply(p)
    assert not (store / "project_stale.md").exists()
    moved = list((store / "retired").glob("*.md"))
    assert len(moved) == 1, moved
    assert "Something learned long ago." in moved[0].read_text(
        encoding="utf-8")


def test_applying_an_empty_proposal_is_quiet(store):
    store.mkdir(parents=True)
    T.apply(T.propose(store))


def test_it_never_raises_on_a_damaged_file(store):
    store.mkdir(parents=True)
    (store / "feedback_broken.md").write_text("not frontmatter at all",
                                              encoding="utf-8")
    _write(store, "fine", "feedback", "Run the gate before pushing.")
    T.apply(T.propose(store))


# -- reachable from both surfaces -------------------------------------------

def test_the_cli_offers_it():
    import inspect

    from delfin.agent import cli

    src = inspect.getsource(cli.build_parser)
    assert '"tidy"' in src, "delfin-agent memory tidy has to exist"
    assert '"--apply"' in src


def test_the_dashboard_offers_it_and_calls_the_same_code():
    import inspect

    from delfin.dashboard import tab_agent as T2

    src = inspect.getsource(T2)
    assert "/tidy" in src, "the dashboard needs a way in"
    assert "memory_tidy" in src, (
        "it must call the same proposal code, not a second one")


# -- choosing what to carry out ---------------------------------------------

def test_a_proposal_can_be_narrowed_to_chosen_names(store):
    """All-or-nothing takes every decision away at once."""
    _write(store, "a", "feedback", "Run the gate before pushing.")
    _write(store, "b", "feedback", "Before pushing, run the gate.")
    _write(store, "old", "project", "Ancient.", age_days=200)
    full = T.propose(store)
    assert len(full.merges) == 1 and len(full.retire) == 1

    only_merge = full.select(["a"])
    assert len(only_merge.merges) == 1 and only_merge.retire == []

    only_retire = full.select(["old"])
    assert only_retire.merges == [] and len(only_retire.retire) == 1


def test_a_merge_is_chosen_by_either_of_its_two_names(store):
    _write(store, "older", "feedback", "Run the gate before pushing.",
           age_days=30)
    _write(store, "newer", "feedback", "Before pushing, run the gate.",
           age_days=1)
    full = T.propose(store)
    assert len(full.select(["older"]).merges) == 1
    assert len(full.select(["newer"]).merges) == 1


def test_choosing_nothing_carries_out_nothing(store):
    _write(store, "a", "feedback", "Run the gate before pushing.")
    _write(store, "b", "feedback", "Before pushing, run the gate.")
    chosen = T.propose(store).select([])
    T.apply(chosen)
    assert len(list(store.glob("*.md"))) == 2


def test_the_counts_follow_the_choice(store):
    _write(store, "a", "feedback", "Run the gate before pushing.")
    _write(store, "b", "feedback", "Before pushing, run the gate.")
    _write(store, "old", "project", "Ancient.", age_days=200)
    chosen = T.propose(store).select(["old"])
    assert chosen.before == 3 and chosen.after == 2


def test_an_unknown_name_selects_nothing_rather_than_everything(store):
    _write(store, "a", "feedback", "Run the gate before pushing.")
    _write(store, "b", "feedback", "Before pushing, run the gate.")
    chosen = T.propose(store).select(["no-such-note"])
    assert chosen.merges == [] and chosen.retire == []


def test_no_model_is_called(store):
    _write(store, 'a', 'feedback', 'Run the gate before pushing.')
    _write(store, 'b', 'feedback', 'Before pushing, run the gate.')
    """Deterministic by decision: the same store proposes the same thing
    twice, and a proposal that changes your memory can be recomputed
    rather than believed."""
    import inspect

    from delfin.agent import memory_tidy as M

    import ast

    tree = ast.parse(inspect.getsource(M))
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported.add(node.module or "")
    # Naming the word "model" in prose is fine; reaching one is not.
    forbidden = ("api_client", "engine", "openai", "anthropic", "httpx",
                 "requests", "memory_distill")
    for mod in imported:
        for bad in forbidden:
            assert bad not in mod, f"{mod} can reach a model"

    # And twice over the same store says the same thing.
    first, second = T.propose(store), T.propose(store)
    assert ([ (m.keep.name, m.drop.name) for m in first.merges ]
            == [ (m.keep.name, m.drop.name) for m in second.merges ])


# -- saying so, once, when there is something to say ------------------------

def test_the_hint_names_what_is_waiting(store):
    _write(store, "a", "feedback", "Run the gate before pushing.")
    _write(store, "b", "feedback", "Before pushing, run the gate.")
    line = T.hint(store)
    assert line and "/tidy" in line
    assert "1" in line


def test_a_tidy_store_says_nothing(store):
    _write(store, "only-one", "feedback", "Run the gate before pushing.")
    assert T.hint(store) == ""


def test_a_missing_store_says_nothing(tmp_path):
    assert T.hint(tmp_path / "gone") == ""


def test_the_hint_counts_retirements_too(store):
    _write(store, "old", "project", "Ancient.", age_days=200)
    line = T.hint(store)
    assert line and "1" in line


def test_the_hint_never_raises(store):
    store.mkdir(parents=True)
    (store / "feedback_broken.md").write_text("garbage", encoding="utf-8")
    T.hint(store)


def test_the_banner_asks_for_it():
    """A button nobody knows about is a button nobody presses."""
    import inspect

    from delfin.agent import cli

    src = inspect.getsource(cli._startup_banner)
    assert "memory_tidy" in src or "tidy_hint" in src, (
        "the startup banner has to mention a waiting tidy")
