"""The ``subagent`` tool's schema follows the preset registry.

Two seams, both from the catalogue being a module-level literal whose
``subagent`` entry was fixed at import:

* the ``subagent_type`` enum was resolved once, so a preset dropped into
  ``~/.delfin/subagents/`` during a session was not offered to the model
  until the process restarted -- which defeats "a specialist is one file";
* the description named the writers in prose: "All types are read-only
  except 'general-purpose' and 'verifier'". That is a sentence about a
  configuration. Measured 2026-11-20 with the two chemistry departments
  switched on, the registry holds FIVE types that can change something --
  general-purpose, metal-smiles-builder, reactivity-department, verifier,
  xyz-editor -- so the sentence was already false and told the model that
  two writers were read-only. A sentence that can be false is worse than
  none, because the model believes it.

``tool_catalog.refresh_subagent_tool`` rebuilds the entry from the live
registry, and ``subagents.reload_subagent_presets`` calls it.

Input: the registry. Output: the ``subagent`` entry's enum and
description.

Semantics pinned here: a preset added and reloaded is offered; one
removed is not; the description names exactly the types that can change
something; and the refresh never raises, because a catalogue that cannot
be refreshed must keep what it had rather than break a turn.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import subagents as S
from delfin.agent import tool_catalog as TC


def _subagent_entry() -> dict:
    for tool in TC._DOC_TOOLS_OPENAI:
        fn = tool.get("function") or {}
        if fn.get("name") == "subagent":
            return fn
    raise AssertionError("the catalogue has no 'subagent' tool")


def _enum() -> list[str]:
    props = (_subagent_entry().get("parameters") or {}).get("properties") or {}
    return list((props.get("subagent_type") or {}).get("enum") or [])


@pytest.fixture
def user_dir(tmp_path, monkeypatch):
    """A stand-in for ``~/.delfin/subagents/``; the operator's own is never
    read or written here."""
    d = tmp_path / "subagents"
    d.mkdir()
    pack = Path(S.__file__).resolve().parent / "pack" / "agents"
    monkeypatch.setattr(S, "_md_preset_search_dirs", lambda: [pack, d])
    yield d
    monkeypatch.undo()
    S.reload_subagent_presets()


def test_a_preset_added_and_reloaded_is_offered(user_dir):
    """The seam: without the refresh this preset is in the registry and
    absent from the schema, so the model cannot name it."""
    before = _enum()
    assert "late-arrival" not in before
    (user_dir / "late_arrival_subagent.md").write_text(
        "---\nname: late-arrival\ndescription: dropped in mid-session\n"
        "mode: plan\ntools: [read_file]\n---\n\nLook and report.\n")
    S.reload_subagent_presets()
    assert "late-arrival" in S.SUBAGENT_PRESETS
    assert "late-arrival" in _enum(), (
        "the registry has it and the tool schema does not, so the model is "
        "never offered the type")


def test_a_preset_removed_is_no_longer_offered(user_dir):
    p = user_dir / "short_lived_subagent.md"
    p.write_text("---\nname: short-lived\ndescription: here and gone\n"
                 "mode: plan\ntools: [read_file]\n---\n\nbody\n")
    S.reload_subagent_presets()
    assert "short-lived" in _enum()
    p.unlink()
    S.reload_subagent_presets()
    assert "short-lived" not in _enum(), (
        "a type the registry no longer has is still offered")


def test_the_description_names_exactly_the_writers():
    """Derived, not stated. A writer is a type whose mode is not 'plan' or
    which carries a write scope."""
    writers = TC._subagent_writers()
    assert writers, "no writer at all; general-purpose should be one"
    descr = _subagent_entry()["description"]
    for w in writers:
        assert repr(w) in descr or w in descr, (w, descr)
    read_only = [n for n, p in S.SUBAGENT_PRESETS.items()
                 if n not in writers]
    for n in read_only:
        assert f"'{n}'" not in descr, (
            f"{n} is read-only and the description lists it as a writer")


def test_a_preset_that_gains_a_write_scope_appears_as_a_writer(user_dir):
    """The case the hard-coded sentence could not follow."""
    (user_dir / "scoped_writer_subagent.md").write_text(
        "---\nname: scoped-writer\ndescription: writes inside one folder\n"
        "mode: default\ntools: [read_file, write_file]\n"
        "write_globs: ['scratch/**']\n---\n\nbody\n")
    S.reload_subagent_presets()
    assert "scoped-writer" in TC._subagent_writers()
    assert "scoped-writer" in _subagent_entry()["description"]


def test_the_refresh_never_raises(monkeypatch):
    """A catalogue that cannot be refreshed keeps what it had. Breaking a
    turn over a schema refresh would be a worse failure than a stale enum.
    """
    def _boom():
        raise RuntimeError("simulated")

    monkeypatch.setattr(S, "subagent_type_names", _boom)
    before = _enum()
    TC.refresh_subagent_tool()          # must not raise
    assert _enum() == before


def test_the_enum_is_not_empty():
    """A refresh that emptied the enum would make every delegation invalid
    while every assertion above about membership still passed."""
    assert len(_enum()) >= 4
