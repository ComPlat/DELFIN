"""FS2: an unchanged text is not a clean text.

Input:  ``output_guard.scrub_secrets_certified`` and the two bug_watcher
        callers whose own docstrings promise to fail closed, each driven
        with one pass of the scrubber made to raise.
Output: whether the caller can tell "no secret present" from "the scan
        never finished", and whether a credential reaches a printed line.
Semantics: ``scrub_secrets`` promises never to raise and to return the
        INPUT when scrubbing fails. That contract is right for a caller
        that merely passes content on -- blanking a tool result on a
        transient keyring failure would cost the run its output -- and
        wrong for a caller whose verdict is a security decision, because
        "unchanged" then means both things at once. Nothing here reads or
        writes the real credential stores: the lookup is replaced.

Run: ~/micromamba/envs/delfin-ci/bin/python -m pytest \
     tests/test_fs2_a_scrub_that_did_not_finish.py -q
"""

import pytest

from delfin.agent import bug_watcher, output_guard

# 32 lowercase hex. Nothing in output_guard._SECRET_CHECKS matches this
# shape, which is the point: a KIT key carries no prefix and no separator,
# so the held-credential lookup is the ONLY pass that can remove it. When
# that pass does not run, no later pass compensates.
KEY = "a3f29c1b7e4d8a6f0b2c5d9e1a4f7c3b"
TEXT = f"RuntimeError: auth failed for token {KEY}"


@pytest.fixture
def lookup_raises(monkeypatch):
    """The held-credential lookup raises, as it does when the keyring is
    locked or the credential file is unreadable."""
    import delfin.agent.web_tools as web_tools

    def boom():
        raise RuntimeError("keyring is locked")

    monkeypatch.setattr(web_tools, "_known_secret_values", boom)


def test_the_shape_passes_cannot_catch_a_key_with_no_shape():
    """The premise of every test below, asserted rather than assumed: if a
    pattern matched this key, a failed exact pass would be harmless."""
    assert KEY in output_guard.scrub_secrets(TEXT, exact_values=())
    scrubbed, certified = output_guard.scrub_secrets_certified(
        TEXT, exact_values=())
    assert KEY in scrubbed
    # An explicitly empty tuple is a caller SAYING "no exact pass", not a
    # failure, so it is still certified -- the degraded-but-known case.
    assert certified is True


def test_a_failed_lookup_is_reported_as_not_certified(lookup_raises):
    value, certified = output_guard.scrub_secrets_certified(TEXT)
    assert certified is False
    assert KEY not in value, ("the certified form never hands back the text "
                              "it could not scrub: a flag beside a raw "
                              "value only helps a caller that reads it")
    assert value == output_guard.WITHHELD


def test_the_certified_form_never_returns_the_text_it_could_not_scrub():
    """The FS1 invariant, pinned at this function rather than only in the
    repo-wide AST guard: whichever pass fails, the value returned is the
    fixed placeholder. Driven over all three failure paths at once."""
    import delfin.agent.web_tools as web_tools

    class Hostile(str):
        def replace(self, *a, **k):
            raise RuntimeError("replace failed part way")

    def boom():
        raise RuntimeError("keyring is locked")

    def boom_pattern(*a, **k):
        raise RuntimeError("catastrophic backtracking")

    cases = []
    import pytest as _pytest
    mp = _pytest.MonkeyPatch()
    try:
        mp.setattr(web_tools, "_known_secret_values", boom)
        cases.append(output_guard.scrub_secrets_certified(TEXT))
        mp.undo()
        cases.append(output_guard.scrub_secrets_certified(
            Hostile(TEXT), exact_values=[KEY]))
        mp.setattr(output_guard, "_redact", boom_pattern)
        cases.append(output_guard.scrub_secrets_certified(
            TEXT, exact_values=[]))
    finally:
        mp.undo()
    assert cases, "three paths were driven"
    for value, certified in cases:
        assert certified is False
        assert value == output_guard.WITHHELD
        assert KEY not in value


def test_a_failed_replace_loop_is_reported(monkeypatch):
    """Half a scrub is not a scrub."""

    class Hostile(str):
        def replace(self, *a, **k):
            raise RuntimeError("replace failed part way")

    _, certified = output_guard.scrub_secrets_certified(
        Hostile(TEXT), exact_values=[KEY])
    assert certified is False


def test_a_failed_shape_pass_is_reported(monkeypatch):
    """_redact_secrets skips a pattern that raises. Skipping silently makes
    a whole class of secret unexamined while the result looks complete."""
    def boom(*a, **k):
        raise RuntimeError("catastrophic backtracking")

    monkeypatch.setattr(output_guard, "_redact", boom)
    failed: list[str] = []
    output_guard._redact_secrets("anything", [], failed=failed)
    assert failed, "every skipped pass is named"
    _, certified = output_guard.scrub_secrets_certified("anything",
                                                        exact_values=())
    assert certified is False


def test_the_detector_fails_closed_as_its_docstring_promises(lookup_raises):
    """_scrubber_leaves says "Fails closed". Before FS2 it asked
    scrub_secrets, which does not raise -- it returns the input -- so the
    except clause could never fire and an unfinished scan read as clean."""
    assert bug_watcher._scrubber_leaves(TEXT) is False


def test_a_credential_is_not_printed_when_the_scan_did_not_finish(lookup_raises):
    """The leak this package exists for: printable_error returns the text
    verbatim for a printed line AND for a webhook, and _scrubber_leaves is
    its only defence. Measured at origin/main with the same fixture, it
    returned 'RuntimeError: auth failed for token <the key>'."""
    out = bug_watcher.printable_error(TEXT)
    assert KEY not in out
    assert "not shown" in out
    assert out.startswith("RuntimeError:"), ("the exception type is still "
                                             "useful and is not a secret")


def test_a_name_field_is_defended_twice(lookup_raises):
    """_shown_name does not depend on the verdict alone: a key-shaped run
    is refused by _KEY_RUN first. Recorded so a later simplification that
    drops one of the two is seen to remove a layer."""
    assert bug_watcher._shown_name(KEY) == bug_watcher._NAME_NOT_SHOWN
    assert bug_watcher._KEY_RUN.search(KEY), "the independent layer"


def test_the_one_value_function_still_never_raises_and_returns_the_input(
        lookup_raises):
    """The old contract is unchanged, because callers that pass content on
    depend on it. Breaking it would blank a tool result whenever the
    keyring is briefly unavailable."""
    assert output_guard.scrub_secrets(TEXT) == TEXT


def test_a_held_credential_is_still_removed_on_the_normal_path(monkeypatch):
    """The fix must not cost the function it protects."""
    import delfin.agent.web_tools as web_tools
    monkeypatch.setattr(web_tools, "_known_secret_values", lambda: (KEY,))
    scrubbed, certified = output_guard.scrub_secrets_certified(TEXT)
    assert KEY not in scrubbed and "[redacted:credential]" in scrubbed
    assert certified is True
    assert bug_watcher.printable_error(TEXT).startswith("RuntimeError: (mess")
