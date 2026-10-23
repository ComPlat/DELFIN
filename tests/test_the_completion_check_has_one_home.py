"""The completion check has ONE home: delfin/agent/task_evidence.py.

Since the swap (night run 2026-09-25), api_client.check_completion_claim
is a thin delegate to task_evidence. The helpers of the old check were
left in api_client so their tests could move first; they have moved
(night run 2026-09-25, Q), and the dead block went with the next commit.
This test keeps a second implementation from growing back next to the
live one: two copies of the same heuristics drift apart.
"""

from __future__ import annotations

from pathlib import Path

import pytest

_SOURCE = Path(__file__).resolve().parents[1] / "delfin" / "agent" / (
    "api_client.py")

# The dead helpers of the old completion check, by their definition
# names. Everything the old check used except the ones the live call
# site still needs (_journal_changes, _journal_ts_epoch, _task_ts_epoch,
# check_completion_claim itself -- the delegate).
_DEAD_HELPERS = (
    "_ARTIFACT_PROMISES",
    "_ARTIFACT_WORD_RES",
    "_artifact_word",
    "_path_suffixes",
    "_unmet_artifact",
    "_TASK_PATH_EXTS",
    "_TASK_PATH_TOKEN_RE",
    "_WRITE_VERB_RE",
    "_READ_VERB_RE",
    "_TEST_TASK_RE",
    "_paths_in_text",
    "_path_matches",
    "_verdict",
)


def test_api_client_has_no_second_completion_check():
    """No definition of an old-check helper remains in api_client.py.

    Reads the source text (not an import) so the check works whether or
    not the names are importable, and counts only DEFINITIONS, so a
    docstring or comment naming a helper does not trip it.
    """
    text = _SOURCE.read_text(encoding="utf-8")
    survivors = [
        name for name in _DEAD_HELPERS
        if _defines(text, name)
    ]
    assert not survivors, (
        "api_client.py still defines the old completion check's helpers: "
        + ", ".join(survivors)
        + " -- delete them; the check lives in task_evidence.py")


def _defines(source_text: str, name: str) -> bool:
    for line in source_text.splitlines():
        stripped = line.strip()
        if (stripped.startswith(name)
                and len(stripped) > len(name)
                and stripped[len(name)] in ":(= "):
            # A definition: "_name(...)...", "_name: ... =", "_name = ..."
            # (attribute-style "obj._name" is excluded by startswith on
            # the stripped line unless it is at column 0).
            return True
    return False
