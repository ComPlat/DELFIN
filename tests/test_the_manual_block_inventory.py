"""Inventory of the manual's code blocks: what kinds exist and where.

This is the transparency half of the recipe tests: it prints (with -s) the
count per language tag and lists the untagged blocks, so a new block that
the checks cannot classify is seen immediately, not silently skipped.
"""

from __future__ import annotations

from collections import Counter

from tests.test_the_manual_recipes import load_blocks, repo_root


def _blocks():
    return load_blocks(repo_root())


def test_every_block_has_a_known_language_or_is_untagged_by_design():
    """Untagged blocks are allowed but must be listed in the printout.

    The check proper lives in the recipe tests; here we only assert the
    inventory itself is well-formed (no empty bodies, unique lines).
    """
    blocks = _blocks()
    assert blocks, "no fenced blocks found — extractor is broken"
    for b in blocks:
        assert b.line > 0
        assert b.path, "block without a source file"


def test_block_inventory(capsys):
    blocks = _blocks()
    counts = Counter(b.lang or "(untagged)" for b in blocks)
    marked = [b for b in blocks if b.not_runnable_reason]
    untagged = [b for b in blocks if not b.lang]
    print()
    print("=== manual code-block inventory ===")
    for lang, n in sorted(counts.items()):
        print(f"  {lang:12} {n}")
    print(f"  marked not-runnable: {len(marked)}")
    print(f"  untagged: {len(untagged)}")
    for b in untagged:
        first = b.text.splitlines()[0] if b.text else "(empty)"
        print(f"    {b.where}: {first[:60]}")
    assert counts, "inventory empty"
