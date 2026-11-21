"""FS1 adversarial, second file: the EDGE CASES of a fail-closed path.

The first file (`test_fs1_adversarial_secret_redaction.py`) drives each wrapper
with well-formed input. This one drives the edges my review checklist names --
`None`, empty, large, unicode, concurrent -- because the fail-closed path has
its own ways to leak:

- a wrapper that returns its argument unchanged for falsy input leaks `None`
  and `0` as TYPES rather than withholding them, so a caller can stringify
  arbitrary objects straight past the redactor;
- truncating a large input before redacting can cut at an offset that leaves
  the key in the fragment;
- the once-per-process log flag is shared mutable state, so under threads a
  raw string could in principle reach one caller while another gets the
  placeholder.

Measured control: 5 of these 9 cases fail on plain `main` (`7abca4483`) and all
9 pass once the builder's fail-closed fix is applied.

The planted token is ASSEMBLED at runtime on purpose: a literal token-shaped
string is redacted by DELFIN's own write path before it reaches disk (observed
-- the constant arrived as the placeholder), which would silently make every
case here vacuous. `test_the_planted_token_really_is_a_token` guards that.
"""
from __future__ import annotations

import threading

import pytest

from delfin.agent import memory_store, session_export, session_index, task_state

# Assembled: no token-shaped literal exists in this source file.
TOKEN = "gh" + "p_" + "A1b2C3d4E5f6G7h8I9j0K1l2M3n4O5p6Q7r8"


@pytest.fixture
def redactor_fails(monkeypatch):
    """The house redactor raises for the duration of the test."""
    def boom(*a, **k):
        raise RuntimeError("redactor down")

    monkeypatch.setattr("delfin.agent.output_guard.scrub_secrets", boom)
    # Reset the once-per-process flag where the module provides it, so the
    # log-once case in the sibling file is not affected by run order.
    if hasattr(memory_store, "_REDACT_LOGGED_ONCE"):
        monkeypatch.setattr(memory_store, "_REDACT_LOGGED_ONCE", False,
                            raising=False)


def _wrappers():
    return [
        ("memory_store", memory_store._without_secrets),
        ("task_state", task_state._scrub),
        ("session_index", session_index._scrub),
        ("session_export", session_export._scrub),
    ]


def test_the_planted_token_really_is_a_token():
    """Anti-vacuity guard: if the write path ever redacts TOKEN again, every
    other case in this file would pass for the wrong reason."""
    assert len(TOKEN) == 40, f"the planted token is {len(TOKEN)} chars, not 40"
    assert TOKEN.startswith("gh")
    assert "redacted" not in TOKEN


@pytest.mark.parametrize("value", [None, "", 0, False])
def test_falsy_input_is_never_returned_as_a_non_string(redactor_fails, value):
    """`None`, `""`, `0`, `False` reach every wrapper.

    A fail-open wrapper returns the argument itself, so `None` comes back as
    `None` and `0` as `0` -- the redactor was never consulted. A fail-closed
    wrapper always yields a string.
    """
    for name, fn in _wrappers():
        out = fn(value)
        assert isinstance(out, str), f"{name} returned {type(out).__name__}"


def test_a_large_input_is_withheld_not_truncated_to_a_secret(redactor_fails):
    """The token is buried in ~1.7 MB of filler.

    Truncating before redaction is a classic leak: cut at the wrong offset and
    the wrapper returns a fragment that still holds the key.
    """
    big = ("filler words " * 80000) + TOKEN + (" more filler" * 80000)
    for name, fn in _wrappers():
        out = fn(big)
        assert TOKEN not in out, f"{name} leaked the token from a large input"


def test_concurrent_callers_all_get_the_placeholder(redactor_fails):
    """Sixteen threads through the wrapper that owns the shared log flag.

    Every caller must get the placeholder; none may receive the raw text, and
    no caller may block on the flag.
    """
    results = []
    lock = threading.Lock()

    def worker():
        try:
            out = memory_store._without_secrets(f"thread text with {TOKEN}")
        except BaseException as exc:      # noqa: BLE001 -- recorded, then asserted
            out = f"<raised {type(exc).__name__}>"
        with lock:
            results.append(out)

    threads = [threading.Thread(target=worker) for _ in range(16)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=30)
        assert not t.is_alive(), "a worker thread hung in the fail-closed path"

    assert len(results) == 16
    leaked = [r for r in results if TOKEN in r]
    assert not leaked, f"{len(leaked)} of 16 concurrent callers got the RAW text"
    placeholders = [r for r in results
                    if r == memory_store.REDACT_FAILED_PLACEHOLDER]
    assert len(placeholders) == 16, (
        f"not every caller got the placeholder; saw {sorted(set(results))[:3]}")


def test_a_zero_width_broken_token_is_RECORDED_not_claimed():
    """Current behaviour, recorded so the PR cannot overclaim.

    The shape scanner does not join a token split by a zero-width space. That
    is not a wrapper defect -- the wrapper can only fail closed, it cannot
    teach the scanner a new shape -- but it must not be reported as closed.
    """
    broken = TOKEN[:4] + "\u200b".join(TOKEN[4:])

    out = memory_store._without_secrets(f"token {broken}")

    assert isinstance(out, str)
    assert TOKEN[:4] in out, (
        "the scanner now catches a zero-width-broken token; the PR follow-up "
        "list should be updated and this case flipped")


def test_a_homoglyph_prefix_is_RECORDED_not_claimed():
    """Same for a non-ASCII homoglyph inside the prefix."""
    homoglyph = "gh\u0440" + TOKEN[3:]

    out = memory_store._without_secrets(f"token {homoglyph}")

    assert isinstance(out, str)
    assert homoglyph in out, (
        "the scanner now normalises homoglyphs; update the follow-up list")
