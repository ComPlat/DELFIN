"""The lazy-module triggers, measured against what a user actually types.

26 % of the largest role prompt sits behind ``<!-- module:X -->`` markers
whose triggers are matched against the conversation text. If a trigger
set misses the way a user phrases something, the guidance is absent from
exactly the turn that needed it.

**Why this file was rewritten.** Its previous version had 49 passing
assertions and measured nothing: every phrase in its corpus contained a
literal trigger substring, because the corpus and the trigger list were
written together. It asserted that the shipped list contains the keywords
its author put in it. Natural variants of the same intent failed for four
of seven modules —

    Richte bitte die Pakete ein            (project_dev missing)
    Kannst du dir die Freigabe merken?     (kit missing)
    Der Job braucht Stunden, starte ihn parallel   (bash_bg missing)
    Optimiere die Struktur der Verbindung  (chemistry missing)
    werte bitte die Belege aus             (documents missing)

— and the test was green throughout.

**What binds it now.** Two things, and neither works without the other:

* ``_INTENT`` is written from the INTENT. Each line is what a user types
  when they want that module, phrased without looking at the trigger
  list. It is the check that the list covers real speech.
* ``_WITNESS`` gives every trigger the one line that justifies it, and
  the mutation test deletes each trigger in turn and requires its witness
  to go dark. A trigger nothing can show a use for is dead weight that
  reads as coverage — twenty such triggers were in the shipped sets,
  each of them a superstring of another in the same module.

Both corpora are German where a German user would write German: they are
matchers for what the USER writes, so they stay German.
"""

from __future__ import annotations

import pytest

from delfin.agent.prompt_loader import PromptLoader


@pytest.fixture
def agent_tree(tmp_path):
    (tmp_path / "pack" / "shared").mkdir(parents=True)
    (tmp_path / "pack" / "agents").mkdir()
    return tmp_path


def _active(agent_tree, text: str, **kw) -> set[str]:
    return PromptLoader(agent_tree)._detect_active_modules(text, "solo", **kw)


def _active_without(agent_tree, module: str, trigger: str, text: str) -> bool:
    """Whether *module* still activates for *text* with *trigger* removed."""
    loader = PromptLoader(agent_tree)
    loader._MODULE_TRIGGERS = {
        name: tuple(t for t in triggers if t != trigger)
        for name, triggers in PromptLoader._MODULE_TRIGGERS.items()
    }
    return module in loader._detect_active_modules(text, "solo")


# ===========================================================================
# 1. The intent corpus — written from what the user wants, not from the list
# ===========================================================================

_INTENT: dict[str, tuple[str, ...]] = {
    "chemistry": (
        # English
        "run a geometry optimization on the ligand",
        "what does the ORCA output say about the SCF",
        "check the imaginary frequencies of the complex",
        "which basis set did we use for the single point",
        "how much energy does that cost",
        # German
        "Optimiere die Geometrie des Komplexes",
        "Optimiere die Struktur der Verbindung",
        "Berechne die Schwingungsfrequenzen des Moleküls",
        "Welches Lösungsmittel hast du im letzten Rechenlauf genommen?",
        "Wie ist die Bindungslänge im Komplex?",
        "Zeig mir die Ladungsverteilung nach Mulliken",
    ),
    "web": (
        "search the internet for the paper",
        "can you look it up online",
        "check the release notes on the project's web page",
        "is there anything on stackoverflow about this",
        "browse to the documentation",
        "Suche im Internet nach dem Paper",
        "Kannst du das online nachschlagen?",
        "Schau mal auf der Webseite nach der neuen Version",
        "Gibt es dazu eine Veröffentlichung?",
        "Recherchier das bitte kurz",
    ),
    "notebook": (
        "open the notebook and run the cell",
        "add a markdown cell to the analysis notebook",
        "the jupyter kernel keeps dying",
        "run analysis.ipynb from the top",
        "restart the kernel and run everything",
        "Öffne das Notizbuch und führe die Zelle aus",
        "Kannst du dir das Jupyter Notebook anschauen?",
        "Füge eine neue Code-Zelle in analysis.ipynb ein",
        "Bitte alle Zellen ausführen",
        "Mach daraus eine Markdown-Zelle",
    ),
    "documents": (
        "read the spreadsheet and total column C",
        "fill in the PDF form for me",
        "summarise budget.xlsx",
        "make a letter for every row of the table",
        "which worksheet holds the numbers",
        "werte bitte die Belege aus",
        "werte bitte die Tabelle aus",
        "die Rechnung prüfen",
        "Erstelle ein Anschreiben aus der Vorlage",
        "Auf welche Kostenstelle geht die Buchung?",
        "Trag die Umsätze aus dem Kontoauszug ein",
    ),
    "project_dev": (
        "install the dependencies",
        "set up the project in a virtual environment",
        "the dependencies in requirements.txt are out of date",
        "the lockfile is stale",
        "scaffold a new module",
        "Installier bitte die Abhängigkeiten",
        "Richte bitte die Pakete ein",
        "Das Projekt in einer virtuellen Umgebung aufsetzen",
        "Die Abhängigkeit fehlt noch im Projekt",
        "Kannst du die Umgebung einrichten?",
    ),
    "kit": (
        "permanently allow pytest",
        "add /home/max/proj to the allowed directories",
        "the sandbox blocked my command again",
        "just allow the command from now on",
        "please remember the permission",
        "merk dir pytest immer erlauben",
        "Kannst du das dauerhaft erlauben?",
        "Kannst du dir die Freigabe merken?",
        "Füge /home/max/proj den erlaubten Verzeichnissen hinzu",
        "Die Berechtigung fehlt noch",
    ),
    "bash_bg": (
        "run the full suite in the background",
        "this is a long job, don't block on it",
        "kick it off and watch progress",
        "start it and do something else while it runs",
        "that takes a while, start it now",
        "lass das im Hintergrund laufen",
        "Der Job läuft lange, starte ihn nebenher",
        "Der Job braucht Stunden, starte ihn parallel",
        "Starte den Testlauf im Hintergrund und mach weiter",
        "Ich will nebenbei weiterarbeiten",
    ),
}


def test_the_intent_corpus_covers_every_module():
    modules = set(PromptLoader._MODULE_TRIGGERS)
    assert set(_INTENT) == modules
    for module, phrases in _INTENT.items():
        assert len(phrases) >= 10, module


@pytest.mark.parametrize(
    "module,phrase",
    [(m, p) for m, ps in _INTENT.items() for p in ps],
)
def test_a_module_activates_on_what_a_user_writes(agent_tree, module, phrase):
    assert module in _active(agent_tree, phrase)


# ===========================================================================
# 2. The witness table — one justifying line per trigger
# ===========================================================================
#
# A witness contains exactly ONE trigger of its module, so deleting that
# trigger must make the line go dark. Identifier-shaped triggers (tool
# names, config filenames, MCP servers) get transcript-shaped witnesses:
# they are matched against the conversation, and the conversation holds
# what the MODEL emitted as well as what the user typed.

_WITNESS: dict[str, str] = {
    # --- chemistry ---------------------------------------------------------
    "orca": "Was steht in der ORCA-Ausgabe?",
    "dft": "Nimm dafür DFT.",
    "xtb": "Nimm dafür xtb.",
    "calc/": "Der Ordner calc/ ist voll.",
    "archive/": "Verschiebe das nach archive/.",
    ".out": "Öffne bitte job.out.",
    "frequencies": "Show me the frequencies.",
    "orbital": "Which orbital is that?",
    "imag": "Are there imaginary modes?",
    "homo": "Where is the HOMO?",
    "lumo": "Where is the LUMO?",
    "uv/vis": "Do the UV/Vis part.",
    "dipole": "How large is the dipole?",
    "scf": "The SCF failed to converge.",
    "mulliken": "Give me the Mulliken charges.",
    "loewdin": "Give me the Loewdin charges.",
    "extract_": "Call extract_charges.",
    "search_calcs": "Use search_calcs.",
    "search_docs": "Use search_docs.",
    "explain_delfin_feature": "Use explain_delfin_feature.",
    "thermochem": "Add the thermochem block.",
    "vibrational": "Show the vibrational modes.",
    "molecule": "Which molecule is it?",
    "complex": "Is it a complex?",
    "ligand": "Which ligand binds?",
    "ml potential": "Try the ML potential.",
    "geometry optimi": "Run a geometry optimization.",
    "transition state": "Find the transition state.",
    "basis set": "Which basis set do we use?",
    "single point": "Just a single point please.",
    "solvent": "Which solvent was that?",
    "energ": "How much energy is that?",
    "quantum chem": "This is quantum chemistry.",
    "geometrie": "Ist die Geometrie sauber?",
    "komplex": "Ist der Komplex stabil?",
    "molekül": "Wie groß ist das Molekül?",
    "molekul": "Wie gross ist das Molekul?",
    "schwingung": "Zeig mir die Schwingung.",
    "frequenzrechnung": "Mach eine Frequenzrechnung.",
    "übergangszustand": "Such den Übergangszustand.",
    "basissatz": "Welcher Basissatz?",
    "funktional": "Welches Funktional nehmen wir?",
    "lösungsmittel": "Welches Lösungsmittel?",
    "rechenlauf": "Der Rechenlauf ist durch.",
    "quantenchem": "Das ist Quantenchemie.",
    "verbindung": "Ist die Verbindung stabil?",
    "struktur der": "Optimiere die Struktur der Probe.",
    "bindungslänge": "Wie ist die Bindungslänge?",
    "ladungsverteilung": "Zeig die Ladungsverteilung.",
    "anregung": "Welche Anregung ist das?",

    # --- web ---------------------------------------------------------------
    "http://": "Hol mir http://beispiel.de",
    "https://": "Öffne https://beispiel.de",
    "web_search": "Nutze web_search.",
    "web_fetch": "Nutze web_fetch.",
    "duckduckgo": "Frag DuckDuckGo.",
    "google": "Frag Google.",
    "stackoverflow": "Steht auf StackOverflow.",
    "the internet": "Ask the internet.",
    "the web": "Ask the web.",
    "online": "Is it online?",
    "web search": "Do a web search.",
    "search for the paper": "Please search for the paper.",
    "look it up": "Just look it up.",
    "release notes": "Check the release notes.",
    "browse to": "Browse to that page.",
    "im internet": "Schau im Internet.",
    "internet nach": "Durchsuche das Internet nach dem Paper.",
    "im netz": "Steht im Netz.",
    "im web": "Steht im Web.",
    "webseite": "Auf der Webseite steht es.",
    "web-seite": "Auf der Web-Seite steht es.",
    "nachschlagen": "Kannst du das nachschlagen?",
    "recherchier": "Recherchier das bitte.",
    "veröffentlichung": "Such die Veröffentlichung.",
    "publikation": "Such die Publikation.",

    # --- notebook ----------------------------------------------------------
    ".ipynb": "Öffne analysis.ipynb.",
    "jupyter": "Starte Jupyter.",
    "notebook": "Open the notebook.",
    "kernel": "The kernel died again.",
    "notizbuch": "Öffne das Notizbuch.",
    "code-zelle": "Füge eine Code-Zelle ein.",
    "codezelle": "Füge eine Codezelle ein.",
    "zellen ausführen": "Alle Zellen ausführen bitte.",
    "markdown-zelle": "Mach daraus eine Markdown-Zelle.",

    # --- documents ---------------------------------------------------------
    ".xlsx": "Öffne budget.xlsx.",
    ".xlsm": "Öffne makros.xlsm.",
    ".csv": "Öffne daten.csv.",
    ".pdf": "Öffne antrag.pdf.",
    ".docx": "Öffne schreiben.docx.",
    "spreadsheet": "Open the spreadsheet.",
    "excel": "Open it in Excel.",
    "tabelle": "Werte die Tabelle aus.",
    "formular": "Fülle das Formular.",
    "pdf form": "Fill the PDF form.",
    "form field": "Which form field is empty?",
    "acroform": "Is it an AcroForm?",
    "invoice": "Send the invoice.",
    "rechnung": "Prüf die Rechnung.",
    "vorlage": "Nimm die Vorlage.",
    "template": "Use the template.",
    "serienbrief": "Mach einen Serienbrief.",
    "anschreiben": "Schreib das Anschreiben.",
    "read_document": "Nutze read_document.",
    "edit_sheet": "Nutze edit_sheet.",
    "fill_pdf_form": "Nutze fill_pdf_form.",
    "create_docx": "Nutze create_docx.",
    "word file": "Open the Word file.",
    "word document": "Open the Word document.",
    "fill in the": "Please fill in the amounts.",
    "fill out the": "Please fill out the rest.",
    "column": "Total that column.",
    "worksheet": "Which worksheet?",
    "the table": "Look at the table.",
    "a letter": "Write a letter.",
    "mail merge": "Do a mail merge.",
    "arbeitsmappe": "Öffne die Arbeitsmappe.",
    "word-datei": "Öffne die Word-Datei.",
    "pdf-datei": "Öffne die PDF-Datei.",
    "csv-datei": "Öffne die CSV-Datei.",
    "ausfüllen": "Kannst du das ausfüllen?",
    "spalte": "Rechne die Spalte zusammen.",
    "buchung": "Prüf die Buchung.",
    "beleg": "Sortiere die Belege.",
    "kostenstelle": "Welche Kostenstelle?",
    "umsatz": "Wie hoch ist der Umsatz?",
    "journal": "Trag das ins Journal.",
    "kontoauszug": "Lies den Kontoauszug.",
    "quittung": "Wo ist die Quittung?",
    "mahnung": "Schreib eine Mahnung.",
    "gutschrift": "Buch die Gutschrift.",
    "lieferschein": "Wo ist der Lieferschein?",
    "zeile eintragen": "Kannst du die Zeile eintragen?",

    # --- project_dev -------------------------------------------------------
    "pyproject.toml": "Schau in pyproject.toml.",
    "package.json": "Schau in package.json.",
    "cargo.toml": "Schau in cargo.toml.",
    "go.mod": "Schau in go.mod.",
    "venv": "Aktivier das venv.",
    "pip install": "Mach pip install requests.",
    "npm ": "Führe npm run build aus.",
    "pnpm": "Wir nehmen pnpm.",
    "yarn": "Wir nehmen yarn.",
    " cargo ": "Bitte cargo build starten.",
    "requirements.txt": "Schau in requirements.txt.",
    "build script": "Das build script ist kaputt.",
    "dependenc": "The dependencies are stale.",
    "virtual environment": "Use a virtual environment.",
    "install the": "Please install the tools.",
    "set up the project": "Can you set up the project?",
    "lockfile": "Update the lockfile.",
    "scaffold": "Scaffold a new module.",
    "abhängigkeit": "Die Abhängigkeit fehlt.",
    "abhaengigkeit": "Die Abhaengigkeit fehlt.",
    "installier": "Installier das bitte.",
    "projekt aufsetzen": "Kannst du das Projekt aufsetzen?",
    "einrichten": "Kannst du das einrichten?",
    "paket": "Lade die Pakete.",
    "umgebung": "Die Umgebung fehlt.",

    # --- kit ---------------------------------------------------------------
    "kit-toolbox": "Nutze die kit-toolbox.",
    "kit_coding": "Nutze kit_coding.",
    "mcp__kit-coding__": "Der Server mcp__kit-coding__bash antwortet nicht.",
    "remember_permission": "Ruf remember_permission auf.",
    "extra_dir": "Setz extra_dir.",
    "kit mode": "Wir sind im kit mode.",
    "permanently allow": "Can you permanently allow pytest?",
    "allow the command": "Just allow the command.",
    "allowed director": "Add it to the allowed directories.",
    "remember the permission": "Please remember the permission.",
    "sandbox": "The sandbox blocked it.",
    "auto-allow": "Turn on auto-allow.",
    "immer erlauben": "pytest immer erlauben bitte.",
    "erlaubte verzeichnis": "Erlaubte Verzeichnisse anzeigen bitte.",
    "erlaubten verzeichnis": "Füge das den erlaubten Verzeichnissen hinzu.",
    "berechtigung": "Die Berechtigung fehlt.",
    "dauerhaft erlaub": "Kannst du das dauerhaft erlauben?",
    "freigabe": "Merk dir die Freigabe.",
    "freigeben": "Kannst du den Ordner freigeben?",
    "zugriff erlauben": "Bitte den Zugriff erlauben.",
    "nicht mehr fragen": "Bitte nicht mehr fragen.",

    # --- bash_bg -----------------------------------------------------------
    "bash_background": "Nutze bash_background.",
    "long running": "It is a long running task.",
    "long-running": "It is a long-running task.",
    "in the background": "Run it in the background.",
    "background job": "Start a background job.",
    "watch progress": "Let me watch progress.",
    "long job": "That is a long job.",
    "takes a while": "This takes a while.",
    "while it runs": "Do something else while it runs.",
    "kick it off": "Just kick it off.",
    "don't block": "Please don't block on it.",
    "im hintergrund": "Starte das im Hintergrund.",
    "läuft lange": "Der Job läuft lange.",
    "laeuft lange": "Der Job laeuft lange.",
    "lange laufen": "Das wird lange laufen.",
    "nebenher": "Mach das nebenher.",
    "dauert lange": "Das dauert lange.",
    "parallel": "Starte es parallel.",
    "stunden": "Das braucht Stunden.",
    "nebenbei": "Mach das nebenbei.",
    "weiterarbeiten": "Ich will weiterarbeiten.",
}

_ALL_TRIGGERS = [(module, trigger)
                 for module, triggers in PromptLoader._MODULE_TRIGGERS.items()
                 for trigger in triggers]


def test_every_trigger_has_a_witness():
    """A trigger nobody can write a line for has no reason to exist."""
    missing = sorted({t for _m, t in _ALL_TRIGGERS if t not in _WITNESS})
    assert not missing, f"triggers with nothing to show for them: {missing}"


@pytest.mark.parametrize("module,trigger", _ALL_TRIGGERS,
                         ids=[f"{m}:{t}" for m, t in _ALL_TRIGGERS])
def test_deleting_any_single_trigger_breaks_its_witness(
        agent_tree, module, trigger):
    """The mutation check. This is what makes the corpus honest: without
    it a trigger list can grow indefinitely and every addition looks
    covered."""
    line = _WITNESS[trigger]
    assert module in _active(agent_tree, line), (trigger, line)
    assert not _active_without(agent_tree, module, trigger, line), (
        f"{trigger!r} is dead weight: {line!r} still activates {module} "
        f"without it")


def test_no_trigger_is_a_superstring_of_another_in_its_module():
    """The structural form of the same rule, so a new trigger cannot be
    added dead. Twenty triggers were in this state when it was written."""
    dead = []
    for module, triggers in PromptLoader._MODULE_TRIGGERS.items():
        for trigger in triggers:
            for other in triggers:
                if other != trigger and other in trigger:
                    dead.append(f"{module}: {trigger!r} ⊇ {other!r}")
    assert not dead, dead


# ===========================================================================
# 3. The triggers stay selective
# ===========================================================================

@pytest.mark.parametrize("task", [
    "fix the failing test in foo.py",
    "format the output of the parser",
    "give me more information about the run",
    "rename the variable in engine.py",
    "warum ist der Test rot?",
    "Benenn die Variable um",
])
def test_ordinary_coding_work_still_activates_nothing(agent_tree, task):
    """Triggering everything would defeat the point of the split."""
    assert _active(agent_tree, task) == set()


# ===========================================================================
# 4. Where the subject comes from
# ===========================================================================

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
    text = ("## Intro\nkept\n\n<!-- module:web -->\n## Web\n"
            "web module body\n\n## Tail\ntail\n")
    loader = PromptLoader(agent_tree)
    off = loader._strip_lazy_modules(
        text, task_text="rename the variable", mode_id="solo")
    on = loader._strip_lazy_modules(
        text, task_text="Suche im Internet nach dem Paper", mode_id="solo")
    assert "web module body" not in off
    assert "web module body" in on
