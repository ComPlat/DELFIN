"""The memory-layer reference says what the store does, and is checked.

Input: docs/MEMORY_LAYERS.md and the modules it describes. Output: nothing
— this file fails when the two disagree.

A reference an agent reads to decide where a fact belongs is only worth
its accuracy. The document arrived describing the layers correctly and
omitting five fields that decide whether a given note reaches a prompt at
all: a reader would have known when facts are retrieved, and not that a
note saved out of an office turn is dropped from a code turn, or that a
model-written note expires while the user's identical one does not.

Each claim here names the thing it is about, so a rename breaks this file
rather than quietly making the document wrong.
"""

from __future__ import annotations

import pathlib

import pytest

_ROOT = pathlib.Path(__file__).resolve().parents[1]
_DOC = _ROOT / "docs" / "MEMORY_LAYERS.md"


def _doc() -> str:
    return _DOC.read_text(encoding="utf-8")


def test_the_reference_exists():
    assert _DOC.is_file(), "the memory-layer reference is gone"


# -- the fields that decide whether a note is recalled ---------------------

@pytest.mark.parametrize("field, module, symbol", [
    ("domain", "memory_store", "memory_text_domain"),
    ("source", "memory_store", "SOURCE_AGENT"),
    ("learned_at", "memory_store", "_head_ref"),
    ("use_count", "memory_store", "use_count"),
    ("stale_hits", "memory_store", "stale_hits"),
])
def test_a_documented_field_exists_in_the_store(field, module, symbol):
    from delfin.agent import memory_store

    assert field in _doc(), f"{field} is not in the reference"
    src = pathlib.Path(memory_store.__file__).read_text(encoding="utf-8")
    assert symbol in src, f"{field} is documented but {symbol} is gone"


def test_the_expiry_rule_is_stated_with_its_number():
    from delfin.agent.memory_store import _AGENT_MEMORY_MAX_AGE_DAYS

    assert str(_AGENT_MEMORY_MAX_AGE_DAYS) in _doc(), (
        "the reference names no expiry, or names a different one than the "
        "store applies")


def test_the_store_key_and_its_default_are_stated():
    from delfin.user_settings import DEFAULT_SETTINGS

    doc = _doc()
    assert "agent.memory_key" in doc
    shipped = (DEFAULT_SETTINGS.get("agent") or {}).get("memory_key")
    assert f'`"{shipped}"`' in doc, (
        f"the reference does not name the shipped default ({shipped!r})")
    assert "-repo-" in doc, "the repository-keyed slug shape is not shown"


def test_retiring_is_documented_as_moving_not_deleting():
    from delfin.agent import memory_tidy

    doc = _doc()
    assert "retired/" in doc, "the reference does not say where a note goes"
    src = pathlib.Path(memory_tidy.__file__).read_text(encoding="utf-8")
    assert "retired" in src and "_retire_file" in src


def test_the_unapplied_report_is_documented_as_unapplied():
    """The one list tidy shows and never acts on. A reader who thinks it
    is applied would not check it."""
    doc = _doc()
    assert "squash" in doc, (
        "the reference does not say why the unlanded-work list is only a "
        "report")
    from delfin.agent import memory_tidy
    src = pathlib.Path(memory_tidy.__file__).read_text(encoding="utf-8")
    assert "unlanded" in src


def test_every_layer_in_the_table_names_a_module_that_exists():
    """The overview table is the part a reader trusts first."""
    import importlib

    doc = _doc()
    for name in ("project_memory", "memory_store", "session_index",
                 "memory_tidy", "memory_nudge"):
        assert f"`{name}`" in doc, f"{name} is missing from the reference"
        importlib.import_module(f"delfin.agent.{name}")
