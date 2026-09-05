"""The staleness checker read chemistry method strings as file paths.

`_PATH_REF_RE` accepted any token containing a slash, so in a
quantum-chemistry agent -- where `functional/basis/dispersion` is simply
how a method is named -- every memory mentioning one was checked for a
file that never existed and annotated

    [stale: BP86/def2-TZVP/D3BJ no longer exists — verify against the
     current code before relying on this]

Measured on the live store: one memory carries stale_hits 2441 against
use_count 2. It has been injected with that warning attached roughly two
and a half thousand times.

The second-order damage is worse than the noise. A non-empty notes list
routes the file to the "rotted" branch, so `record_memory_recall` is
skipped and `record_stale_hits` runs instead -- which means `updated_at`
freezes. `prune_memories` sorts survivors by `updated_at`. So the newer,
CORRECTIVE memory is the one queued for eviction, while the superseded
one it was written to override is refreshed on every turn and is
effectively immortal.

A path reference has to look like a path: a dotted final segment, or a
prefix that exists on disk.
"""

from __future__ import annotations

import pytest

from delfin.agent import memory_store as M


def _refs(text):
    return [r[1] for r in M._extract_path_refs(text)]


# ---------------------------------------------------------------------------
# What is not a path
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("text", [
    "BP86/def2-TZVP/D3BJ als aktuelle Einstellung gewählt",
    "use CAM-B3LYP/def2-SVP for the excited states",
    "the S1/S2 gap",
    "a signal/noise ratio of 3",
    "wB97X-D/def2-TZVPP",
])
def test_a_method_string_is_not_a_path(text):
    assert _refs(text) == [], (text, _refs(text))


def test_a_bare_two_word_slash_is_not_a_path():
    assert _refs("the input/output behaviour") == []


# ---------------------------------------------------------------------------
# What still is
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("text,want", [
    ("see delfin/agent/engine.py for the loop", "delfin/agent/engine.py"),
    ("`delfin/agent/api_client.py:1654` has it", "delfin/agent/api_client.py"),
    ("the file engine.py holds it", "engine.py"),
    ("check tests/test_agent_engine.py first", "tests/test_agent_engine.py"),
    ("config.yaml is read at start", "config.yaml"),
])
def test_a_real_path_is_still_found(text, want):
    assert want in _refs(text), (text, _refs(text))


def test_a_line_number_is_still_captured():
    refs = M._extract_path_refs("delfin/agent/engine.py:1654 is the call")
    assert refs and refs[0][2] == 1654


def test_an_existing_directory_is_still_a_path(tmp_path, monkeypatch):
    """A directory has no dotted segment, so it needs the other rule."""
    (tmp_path / "delfin" / "agent").mkdir(parents=True)
    monkeypatch.chdir(tmp_path)
    assert "delfin/agent/" in _refs("see delfin/agent/ for the module")


# ---------------------------------------------------------------------------
# The consequence the noise had
# ---------------------------------------------------------------------------

def test_a_method_memory_gets_no_staleness_note(tmp_path):
    body = ("User hat BP86/def2-TZVP/D3BJ als aktuelle Einstellung "
            "gewählt; das überschreibt temporär die gespeicherte "
            "def2-svp-Präferenz.")
    notes = M.recall_reference_notes(body, repo_root=tmp_path)
    assert notes == [], notes


def test_a_memory_naming_a_deleted_file_still_warns(tmp_path):
    notes = M.recall_reference_notes(
        "the retry loop lives in delfin/agent/nope_invented.py",
        repo_root=tmp_path)
    assert notes and "nope_invented.py" in notes[0]
