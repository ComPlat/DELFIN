"""A pasted brief that starts with a heading is a brief, not a note.

``#`` at the prompt writes a memory. The check was
``line.startswith("#")`` against the WHOLE submitted text, and a brief
pasted into a terminal arrives as one submission. So a task beginning

    # Gemeinsamer Teil

went into _remember() instead of to the model. On 2026-09-20 that was
done to five sessions at once: every one of them printed a message about
remembering and then sat at an empty prompt, and the model never saw a
word of its assignment. The whole first run was lost to it.

Underneath sat a second defect. _remember() called
``save_typed_memory(text=…, kind=…, workspace=…, author=…)`` while the
store takes ``memory_type=``, ``repo_root=`` and ``scope=``, so the note
feature had never worked on this path at all -- the error message about
remembering was the only thing it ever produced.

  one line beginning with #    still a note; that is the feature
  several lines               a brief, and goes to the model
  a heading with nothing else  still a note -- "#" plus one word is not
                               a document, and guessing otherwise would
                               take the feature away
"""

from __future__ import annotations

import pytest


class _Theme:
    def dim(self, t): return t
    def red(self, t): return t
    def bold(self, t): return t


class _Transcript:
    def __init__(self): self.lines = []
    theme = _Theme()
    def chrome(self, line): self.lines.append(line)


def _agent(tmp_path):
    from delfin.agent.repl import TerminalAgent, ReplOptions
    a = object.__new__(TerminalAgent)
    a.opts = ReplOptions(cwd=tmp_path)
    a.transcript = _Transcript()
    a._quit = False
    # Enough engine for _ctx(): a brief that is NOT a note walks on to
    # the command dispatch, which is the path this case is about.
    a.engine = type("E", (), {"session_id": "s-1", "messages": []})()
    return a


class TestWhatCountsAsANote:
    def test_one_line_is_a_note(self):
        from delfin.agent.repl import _is_a_note
        assert _is_a_note("#remember the gate is wide here")

    def test_several_lines_are_not(self):
        from delfin.agent.repl import _is_a_note
        assert not _is_a_note("# Gemeinsamer Teil\n\nDu arbeitest an DELFIN.")

    def test_a_bare_heading_is_still_a_note(self):
        from delfin.agent.repl import _is_a_note
        assert _is_a_note("# one thought")

    def test_text_that_does_not_start_with_it_is_not(self):
        from delfin.agent.repl import _is_a_note
        assert not _is_a_note("the # is in the middle")

    def test_a_lone_hash_is_not_a_note(self):
        from delfin.agent.repl import _is_a_note
        assert not _is_a_note("#")

    def test_trailing_blank_lines_do_not_make_it_a_brief(self):
        from delfin.agent.repl import _is_a_note
        assert _is_a_note("#remember this\n\n")


class TestTheBriefReachesTheModel:
    def test_a_multi_line_brief_is_passed_through(self, tmp_path):
        brief = "# Gemeinsamer Teil\n\nDu arbeitest an DELFINs Agenten."
        agent = _agent(tmp_path)
        assert agent._handle_line(brief) == brief

    def test_nothing_was_remembered(self, tmp_path):
        agent = _agent(tmp_path)
        agent._handle_line("# Auftrag 1\n\nBau die Liste aus.")
        assert not any("remember" in line for line in agent.transcript.lines)


class TestRememberingWorksAtAll:
    def test_a_note_is_stored(self, tmp_path, monkeypatch):
        from delfin.agent import memory_store
        seen = {}

        def _save(text, **kw):
            seen["text"] = text
            seen["kw"] = kw
            return (tmp_path / "m.md", "name", "user")

        monkeypatch.setattr(memory_store, "save_typed_memory", _save)
        agent = _agent(tmp_path)
        agent._handle_line("#the gate is wider in a DELFIN checkout")
        assert seen["text"] == "the gate is wider in a DELFIN checkout"
        assert "remembered" in " ".join(agent.transcript.lines)

    def test_it_passes_the_names_the_store_actually_takes(self, tmp_path,
                                                          monkeypatch):
        import inspect

        from delfin.agent import memory_store
        taken = set(inspect.signature(
            memory_store._save_typed_memory_unlocked).parameters)
        seen = {}
        monkeypatch.setattr(memory_store, "save_typed_memory",
                            lambda text, **kw: seen.update(kw) or
                            (tmp_path / "m.md", "n", "user"))
        _agent(tmp_path)._handle_line("#a thought")
        unknown = sorted(set(seen) - taken)
        assert not unknown, f"the store does not take: {unknown}"
