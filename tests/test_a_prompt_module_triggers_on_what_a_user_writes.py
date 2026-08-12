"""The lazy-module triggers matched identifiers, not what a user writes.

26 % of the largest role prompt (~11 kB of solo_agent.md) sits behind
``<!-- module:X -->`` markers whose triggers are matched against the task
text. The trigger sets were the identifiers a MODEL emits — tool names
(``web_search``, ``notebook_read``), file extensions (``.ipynb``), config
filenames (``pyproject.toml``). Nothing a user types looks like that.

Measured on the shipped sets, all of these activated NOTHING:

    "search the internet for the paper"          (web missing)
    "open the notebook and run the cell"         (notebook missing)
    "install the dependencies"                   (project_dev missing)
    "Suche im Internet nach dem Paper"
    "Öffne das Notizbuch und führe die Zelle aus"
    "Optimiere die Geometrie des Komplexes"      (chemistry missing)
    "lass das im Hintergrund laufen"             (only the English form hit)

So the guidance was absent from exactly the turn that needed it, and
11 kB of prompt was unreachable in practice. Only the documents module had
German terms — in a framework whose users write German.

The corpus below is the contract: every module must activate on at least
three natural English AND three natural German user phrasings. The German
entries are matchers for what the USER writes, so they stay German.
"""

from __future__ import annotations

import pytest


@pytest.fixture
def agent_tree(tmp_path):
    (tmp_path / "pack" / "shared").mkdir(parents=True)
    (tmp_path / "pack" / "agents").mkdir()
    return tmp_path


def _active(agent_tree, text: str, **kw) -> set[str]:
    from delfin.agent.prompt_loader import PromptLoader
    return PromptLoader(agent_tree)._detect_active_modules(text, "solo", **kw)


# --- the corpus ------------------------------------------------------------
# Written as a user would type it: no tool names, no file extensions unless
# a user would naturally name a file.

_ENGLISH: dict[str, tuple[str, ...]] = {
    "chemistry": (
        "run a geometry optimization on the ligand",
        "what does the ORCA output say about the SCF",
        "check the imaginary frequencies of the complex",
    ),
    "web": (
        "search the internet for the paper",
        "can you look it up online",
        "check the release notes on the project's web page",
    ),
    "notebook": (
        "open the notebook and run the cell",
        "add a markdown cell to the analysis notebook",
        "the jupyter kernel keeps dying",
    ),
    "documents": (
        "read the spreadsheet and total column C",
        "fill in the PDF form for me",
        "summarise budget.xlsx",
    ),
    "project_dev": (
        "install the dependencies",
        "set up the project in a virtual environment",
        "the dependencies in requirements.txt are out of date",
    ),
    "kit": (
        "permanently allow pytest",
        "add /home/max/proj to the allowed directories",
        "the sandbox blocked my command again",
    ),
    "bash_bg": (
        "run the full suite in the background",
        "this is a long job, don't block on it",
        "kick it off and watch progress",
    ),
}

_GERMAN: dict[str, tuple[str, ...]] = {
    "chemistry": (
        "Optimiere die Geometrie des Komplexes",
        "Berechne die Schwingungsfrequenzen des Moleküls",
        "Welches Lösungsmittel hast du im letzten Rechenlauf genommen?",
    ),
    "web": (
        "Suche im Internet nach dem Paper",
        "Kannst du das online nachschlagen?",
        "Schau mal auf der Webseite nach der neuen Version",
    ),
    "notebook": (
        "Öffne das Notizbuch und führe die Zelle aus",
        "Kannst du dir das Jupyter Notebook anschauen?",
        "Füge eine neue Code-Zelle in analysis.ipynb ein",
    ),
    "documents": (
        "werte bitte die Tabelle aus",
        "die Rechnung prüfen",
        "Erstelle ein Anschreiben aus der Vorlage",
    ),
    "project_dev": (
        "Installier bitte die Abhängigkeiten",
        "Das Projekt in einer virtuellen Umgebung aufsetzen",
        "Die Abhängigkeit fehlt noch im Projekt",
    ),
    "kit": (
        "merk dir pytest immer erlauben",
        "Kannst du das dauerhaft erlauben?",
        "Füge /home/max/proj den erlaubten Verzeichnissen hinzu",
    ),
    "bash_bg": (
        "lass das im Hintergrund laufen",
        "Der Job läuft lange, starte ihn nebenher",
        "Starte den Testlauf im Hintergrund und mach weiter",
    ),
}


def test_the_corpus_covers_every_module():
    from delfin.agent.prompt_loader import PromptLoader
    modules = set(PromptLoader._MODULE_TRIGGERS)
    assert set(_ENGLISH) == modules
    assert set(_GERMAN) == modules
    for corpus in (_ENGLISH, _GERMAN):
        for module, phrases in corpus.items():
            assert len(phrases) >= 3, module


@pytest.mark.parametrize(
    "module,phrase",
    [(m, p) for m, ps in _ENGLISH.items() for p in ps],
)
def test_a_module_activates_on_natural_english(agent_tree, module, phrase):
    assert module in _active(agent_tree, phrase)


@pytest.mark.parametrize(
    "module,phrase",
    [(m, p) for m, ps in _GERMAN.items() for p in ps],
)
def test_a_module_activates_on_natural_german(agent_tree, module, phrase):
    assert module in _active(agent_tree, phrase)


# --- the triggers stay selective ------------------------------------------

@pytest.mark.parametrize("task", [
    "fix the failing test in foo.py",
    "format the output of the parser",
    "give me more information about the run",
    "rename the variable in engine.py",
])
def test_ordinary_coding_work_still_activates_nothing(agent_tree, task):
    """Triggering everything would defeat the point of the split."""
    assert _active(agent_tree, task) == set()


# --- what the user said two messages ago still counts ----------------------

def test_the_subject_is_read_from_the_conversation_not_only_the_last_line(
        agent_tree):
    """A task line is one message. The subject is often two messages back
    ("Öffne das Notizbuch" → "ja, mach weiter"), and a resumed session
    starts on a follow-up with the subject only in the history."""
    assert _active(agent_tree, "ja, mach weiter") == set()
    assert "notebook" in _active(
        agent_tree, "ja, mach weiter",
        conversation_text="Öffne das Notizbuch und führe die Zelle aus")


def test_a_module_stripped_from_the_prompt_comes_back_when_it_is_needed(
        agent_tree):
    """End to end through the stripper the prompt actually goes through."""
    from delfin.agent.prompt_loader import PromptLoader

    text = ("## Intro\nkept\n\n<!-- module:web -->\n## Web\n"
            "web module body\n\n## Tail\ntail\n")
    loader = PromptLoader(agent_tree)
    off = loader._strip_lazy_modules(
        text, task_text="rename the variable", mode_id="solo")
    on = loader._strip_lazy_modules(
        text, task_text="Suche im Internet nach dem Paper", mode_id="solo")
    assert "web module body" not in off
    assert "web module body" in on
