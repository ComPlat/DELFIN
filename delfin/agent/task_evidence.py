"""What a completed task can show for its work -- judged on the CLAIM,
not on the words.

Drop-in replacement for ``api_client.check_completion_claim`` with the
same arguments and the same return shape (``{verdict, kind, detail,
note}``, verdict one of ``verified`` / ``unmet`` / ``unchecked``).

WHY THIS MODULE EXISTS. The supervised run of 2026-09-22 (s4) showed the
old check keying on WORDING rather than work: "eine Tabelle der
Testfaelle" in a task that wrote a markdown section was reported
``unmet`` for lack of a .xlsx; a file only mentioned as a reference
("wie in foo.py beschrieben") made the task ``unmet`` until foo.py
itself was written. The pattern behind both: a word or a filename is
treated as a claim whenever it OCCURS in the text, whatever its role.

The distinction this module draws instead: a word or a path is a claim
only when it is the OBJECT of the task. Concretely,

* a named path accuses only when it comes from the SUBJECT, in no
  mention context ("wie in X", "analog zu X", "kam aus X", "as in X",
  "based on X", ...), and no write of it is recorded. A path from the
  DESCRIPTION, or one that is only referenced, is context: it never
  turns a finished task ``unmet``.
* an artefact noun ("Tabelle", "Bericht") promises a file TYPE only
  with an explicit format word (PDF, Excel, CSV, Word, pptx, ...)
  attached: "PDF-Bericht erstellen", "Erstelle den Bericht als PDF",
  "Excel-Tabelle anlegen". Without one, the word describes content, not
  a container -- the change journal decides.
* the write/read intent is read from the SUBJECT, not from subject +
  description: a note instruction in the description ("notiere
  Stichpunkte in der Antwort") cannot turn a read task into a write
  task.

And the bias throughout: where the text cannot decide, the answer is
``unchecked``. A false ``unmet`` makes the agent doubt or redo finished
work; ``unchecked`` is an honest unknown. Every heuristic here that
judged correctly yesterday keeps its verdict -- pinned in
``tests/test_task_evidence_characterization.py``.

The verb tables and the path extraction are the ones the old check
measured against real office runs (``delfin/agent/german.py`` carries
their history); they are reused, not reinvented.
"""

from __future__ import annotations

import re
from pathlib import Path

from .german import GERMAN_READ_VERB_SOURCE, GERMAN_WRITE_VERB_SOURCE
from . import evidence_freshness as evidence

__all__ = ["check_completion_claim"]


# ---------------------------------------------------------------------------
# Reused vocabulary -- measured against real office runs, see german.py
# ---------------------------------------------------------------------------

# Verbs that promise a CHANGE.
_WRITE_VERB_RE = re.compile(
    r"(?i)\b(?:add|create|write|implement|build|generate|produce|export|"
    r"save|fix|repair|patch|refactor|rename|move|update|extend|port|"
    r"migrate|wire|integrate|remove|delete|split|merge|fill\s+in|"
    r"fill\s+out|enter|record|book|sort|"
    r"erstell\w*|schreib\w*|füg\w*|hinzufüg\w*|implementier\w*|bau\w*|"
    r"erzeug\w*|generier\w*|exportier\w*|speicher\w*|beheb\w*|"
    r"reparier\w*|korrigier\w*|refaktor\w*|umbenenn\w*|verschieb\w*|"
    r"aktualisier\w*|erweiter\w*|entfern\w*|lösch\w*|anpass\w*|änder\w*|"
    r"einbau\w*|einbind\w*|ergänz\w*|umstell\w*|überarbeit\w*|"
    r"integrier\w*|umbau\w*|aufräum\w*|bereinig\w*)\b"
    r"|" + GERMAN_WRITE_VERB_SOURCE
)

# Verbs that promise only LOOKING.
_READ_VERB_RE = re.compile(
    r"(?i)\b(?:read|review|analyse|analyze|inspect|examine|check|"
    r"understand|summarise|summarize|compare|explore|investigate|study|"
    r"audit|trace|total|count|"
    r"lies|lese\w*|les\w*|prüf\w*|überprüf\w*|analysier\w*|untersuch\w*|"
    r"sicht\w*|versteh\w*|vergleich\w*|durchsuch\w*|betracht\w*|"
    r"anschau\w*|ansehen|recherchier\w*|bewert\w*)\b"
    r"|" + GERMAN_READ_VERB_SOURCE
)

_TEST_TASK_RE = re.compile(
    r"(?i)(?:\b(?:tests?|testing|testsuite|test-suite|pytest|unittest|"
    r"regression|verify|verification|validate|"
    r"teste\w*|testen|verifizier\w*|validier\w*)\b"
    r"|\b\w+tests\b|\b\w+testsuite\b)"
)

# Extensions a subject can name unambiguously enough to key a check on.
_TASK_PATH_EXTS: frozenset[str] = frozenset({
    "py", "pyi", "ipynb", "js", "jsx", "ts", "tsx", "json", "yaml", "yml",
    "toml", "cfg", "ini", "md", "rst", "txt", "csv", "tsv", "sh", "bash",
    "c", "h", "cc", "cpp", "hpp", "rs", "go", "java", "kt", "rb", "php",
    "jl", "sql", "html", "htm", "css", "scss", "xml", "tex", "pdf",
    "docx", "doc", "xlsx", "xls", "pptx", "png", "svg", "jpg", "jpeg",
    "log", "inp", "out", "xyz", "mol", "sdf", "cif", "dat", "gjf",
})

_TASK_PATH_TOKEN_RE = re.compile(r"[A-Za-z0-9_./~+-]{3,}")


# ---------------------------------------------------------------------------
# A path that is MENTIONED, not targeted
# ---------------------------------------------------------------------------

# Constructions that make a following path a REFERENCE, not the object of
# the task: "wie in foo.py", "analog zu check_base.py", "kam aus foo.py",
# "as described in foo.py", "based on foo.py", "siehe foo.py".
# The path names where the knowledge came from, not where the work goes.
_MENTION_BEFORE_RE = re.compile(
    r"(?i)\b(?:"
    r"wie\s+in|analog\s+zu|analog\w*|vergleichbar\s+mit|entsprechend|"
    r"nach\s+Art\s+von|siehe|vgl\.?|vergleiche\s+mit|"
    r"based\s+on|as\s+described\s+in|as\s+in|as\s+per|following|"
    r"see|cf\.|compared\s+to|"
    # Where something came FROM, when it is not the thing being changed:
    # "der Fehler kam aus foo.py". Only in the passive/past framing --
    # "aus X erstellen" is a write and must not match.
    r"stammt\s+aus|kam\s+aus|kommt\s+aus|kommt\s+von|stammt\s+von|"
    r"originates\s+from|came\s+from|comes\s+from|"
    r"auf\s+Basis\s+von|anhand\s+von"
    r")\s*$"
)


def _is_mentioned(path: str, text: str) -> bool:
    """Whether *path* occurs in *text* only in a mention construction.

    True when every occurrence is preceded (within its sentence fragment)
    by a reference construction: "wie in", "analog zu", "kam aus", "as
    described in", ... A path that also stands on its own -- "Behebe den
    Fehler in bar.py, wie in foo.py" -- is both, and being an object
    anywhere wins.
    """
    try:
        for m in re.finditer(re.escape(path), text, re.IGNORECASE):
            before = text[max(0, m.start() - 40):m.start()]
            # Stop at a clause boundary: a construction only governs the
            # path directly after it.
            before = re.split(r"[.;:!?\n]", before)[-1]
            if _MENTION_BEFORE_RE.search(before):
                continue
            return False
        return bool(re.search(re.escape(path), text, re.IGNORECASE))
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Artefact nouns and the formats they can promise
# ---------------------------------------------------------------------------

# A format word makes a noun promise a FILE TYPE. Without one, "Tabelle"
# and "Bericht" describe content (a table in a markdown section is a
# table) and the change journal decides.
_FORMAT_PROMISES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("pdf", (".pdf",)),
    ("docx", (".docx",)),
    ("word", (".docx", ".doc")),
    ("excel", (".xlsx", ".xls")),
    ("xlsx", (".xlsx",)),
    ("csv", (".csv",)),
    ("pptx", (".pptx",)),
    ("presentation", (".pptx",)),
)

# Format words as whole words (with German inflection), and compounds
# joined by hyphen or directly: "PDF-Bericht", "Excel-Tabelle",
# "Word-Dokument". "bericht als pdf" (format after the noun) is covered
# by the second check below.
_FORMAT_WORD_RE = re.compile(
    # "word" only as a compound ("Word-Dokument"): bare English "word"
    # ("explain the word Tabelle") is not a format promise.
    r"(?i)(?:(?:\b|\w+-)(?:pdf|docx|xlsx|pptx|excel|csv)\b"
    r"|\bword-|\bword\s+datei\b|\bword\s+dokument\w*\b)"
)


def _format_word(text: str) -> str:
    """The file-format word a subject carries, "" when none does."""
    try:
        m = _FORMAT_WORD_RE.search(text or "")
        return m.group(0).lower().lstrip("-") if m else ""
    except Exception:
        return ""


def _unmet_format(subject: str, produced) -> str:
    """The format a subject explicitly promises but no written file has.

    "" when the subject carries no format word, or when a written path
    carries the promised extension. *produced* must be write evidence
    only -- reads prove nothing (see the old check's history).
    """
    try:
        word = _format_word(subject)
        if not word:
            return ""
        # Map the matched word back onto its promise table entry.
        for name, extensions in _FORMAT_PROMISES:
            if name in word:
                suffixes = {Path(str(p).replace("\\", "/")).suffix.lower()
                            for p in (produced or ())}
                if any(e in suffixes for e in extensions):
                    return ""
                return name
        return ""
    except Exception:
        return ""


# ---------------------------------------------------------------------------
# Path extraction and matching
# ---------------------------------------------------------------------------

def _paths_in_text(text) -> list[str]:
    """File paths a text names, in order (same rules as before: a closed
    extension set, so "version 1.2" and "z.B." are not paths)."""
    out: list[str] = []
    try:
        for token in _TASK_PATH_TOKEN_RE.findall(str(text or "")):
            token = token.strip(".,;:/")
            if "." not in token:
                continue
            stem, _, ext = token.rpartition(".")
            if ext.lower() not in _TASK_PATH_EXTS:
                continue
            if not stem or not any(c.isalnum() for c in stem):
                continue
            if token not in out:
                out.append(token)
    except Exception:
        return out
    return out


def _path_matches(candidate: str, recorded: str) -> bool:
    """Suffix match, so "mylib/opt/wrapper.py" matches the absolute
    "/home/u/proj/mylib/opt/wrapper.py" and nothing shorter."""
    try:
        c = str(candidate).replace("\\", "/").strip("/").lower()
        r = str(recorded).replace("\\", "/").strip("/").lower()
        return bool(c) and bool(r) and (r == c or r.endswith("/" + c))
    except Exception:
        return False


def _verdict(kind: str, verdict: str, detail: str = "", note: str = "") -> dict:
    return {"verdict": verdict, "kind": kind, "detail": detail, "note": note}


# ---------------------------------------------------------------------------
# The check
# ---------------------------------------------------------------------------

def check_completion_claim(
    subject: str,
    description: str = "",
    *,
    changes=(),
    observed=None,
    tests=None,
    window_start: float = 0.0,
    current_fingerprint: dict | None = None,
) -> dict:
    """What the session can show for a task about to be marked completed.

    Same contract as ``api_client.check_completion_claim``: returns
    ``{"verdict", "kind", "detail", "note"}`` with verdict
    ``verified`` / ``unmet`` / ``unchecked``; pure over its arguments
    (*changes* the write ledger ``[{path, ts, created}]``, *observed*
    the read ledger, *tests* the test-evidence ledger, *window_start*
    the epoch the task went in_progress).

    The difference is WHAT counts as a claim. Only the object of the
    task accuses: a path from the subject that is not merely mentioned
    (``_is_mentioned``), or a file format the subject explicitly names
    (``_unmet_format``). A path that only occurs in the description, or
    inside a reference construction, is context. Where the text cannot
    decide, the answer is ``unchecked`` -- a false ``unmet`` makes the
    agent doubt or redo finished work, and an honest unknown never
    does.
    """
    try:
        subject_text = str(subject or "")
        changed = [c for c in (changes or ()) if isinstance(c, dict)]
        written = [str(c.get("path", "")) for c in changed if c.get("path")]
        read_paths = None if observed is None else [
            str(p) for p in (observed or ())]

        # Intent comes from the SUBJECT. A description can say "notiere
        # Stichpunkte" without turning a read task into a write task.
        wants_write = bool(_WRITE_VERB_RE.search(subject_text))
        reads_only = (not wants_write
                      and bool(_READ_VERB_RE.search(subject_text)))

        # 1. A path the SUBJECT names as its object. The mention guard
        #    applies only where the path must be WRITTEN: for a read
        #    task the file is the object being looked at, however the
        #    sentence carries it.
        subject_paths = [p for p in _paths_in_text(subject_text)
                         if reads_only or not _is_mentioned(p, subject_text)]
        for cand in subject_paths:
            if any(_path_matches(cand, p) for p in written):
                return _verdict("path_write", "verified", cand)
        unmatched = [c for c in subject_paths
                     if not any(_path_matches(c, p) for p in written)]
        if unmatched:
            cand = unmatched[0]
            seen = (None if read_paths is None
                    else any(_path_matches(cand, p) for p in read_paths))
            if reads_only:
                if seen:
                    return _verdict("path_read", "verified", cand)
                if seen is None:
                    return _verdict("path_read", "unchecked", cand)
                return _verdict(
                    "path_untouched", "unmet", cand,
                    f"names {cand} and this session neither read nor wrote "
                    f"it.")
            if wants_write:
                return _verdict(
                    "path_unwritten", "unmet", cand,
                    f"names {cand} and no write of it is recorded"
                    + (" (it was only read)" if seen else "") + ".")
            if seen:
                return _verdict("path_read", "verified", cand)
            if seen is None:
                return _verdict("path_read", "unchecked", cand)
            return _verdict(
                "path_untouched", "unmet", cand,
                f"names {cand} and this session neither read nor wrote it.")

        # 2. An explicitly promised file format. Bare "Tabelle" /
        #    "Bericht" describe content, not a container; the journal
        #    decides those below.
        if not reads_only:
            missing = _unmet_format(subject_text, written)
            if missing:
                return _verdict(
                    "artifact", "unmet", missing,
                    f"names a {missing} and no {missing} file was written "
                    f"in this session.")
            if _format_word(subject_text):
                return _verdict("artifact", "verified",
                                _format_word(subject_text))

        # 3. A test / verification task needs a green run in the window --
        #    and that run must still describe the CURRENT tree state
        #    (night run 2026-09-26, assignment Y: a green count was quoted
        #    after the code beneath it had changed; the "2 passed" belonged
        #    to an older tree). When a *current_fingerprint* is supplied,
        #    a green run whose state fingerprint is stale does not verify.
        if _TEST_TASK_RE.search(subject_text):
            if tests is None:
                return _verdict("tests", "unchecked", "no test ledger")
            entries = [e for e in tests if isinstance(e, dict)]
            in_window = [
                e for e in entries
                if float(e.get("ts", 0) or 0) >= float(window_start or 0)
            ]
            if current_fingerprint is not None:
                # Entries from ledgers that predate stamping carry no
                # fingerprint; is_stale reports them stale, the safe
                # direction for an unknown state.
                fresh_green = [
                    e for e in in_window
                    if int(e.get("failed", 0) or 0) == 0
                    and str(e.get("status", "")) not in ("failed", "error",
                                                         "gave_up")
                    and not evidence.is_stale(e, current_fingerprint)]
            else:
                fresh_green = [
                    e for e in in_window
                    if int(e.get("failed", 0) or 0) == 0
                    and str(e.get("status", "")) not in ("failed", "error",
                                                         "gave_up")]
            if fresh_green:
                return _verdict("tests", "verified",
                                f"{len(fresh_green)} green run(s)")
            if current_fingerprint is not None and in_window:
                # Green runs exist but every in-window run is stale for
                # the current state: the numbers belong to an older
                # tree -- name that, do not report it as "N failing".
                stale_in_window = [
                    e for e in in_window
                    if evidence.is_stale(e, current_fingerprint)]
                if stale_in_window and len(stale_in_window) == len(in_window):
                    why = evidence.note(evidence.is_stale(
                        stale_in_window[-1], current_fingerprint))
                    return _verdict(
                        "tests_stale_state", "unmet", "state moved on",
                        f"is a test task and every recorded run predates the "
                        f"current tree state -- {why}")
            if in_window:
                worst = max(int(e.get("failed", 0) or 0) for e in in_window)
                return _verdict(
                    "tests_red", "unmet", f"{worst} failing",
                    f"is a test task and the runs recorded since it started "
                    f"were not green ({worst} failure(s) in the last one).")
            if entries:
                return _verdict(
                    "tests_stale", "unmet", "before the task started",
                    "is a test task and every recorded test run predates "
                    "it — nothing was run since the work began.")
            return _verdict(
                "tests_none", "unmet", "no run recorded",
                "is a test task and this session recorded no test run at "
                "all.")

        # 4. An edit / refactor task with no path named: the journal has
        #    to show a mutation. A change earlier in the session still
        #    counts -- editing first and flipping the status afterwards
        #    is doing the work, not faking it.
        if wants_write:
            if any(c.get("ts", 0) >= float(window_start or 0)
                   for c in changed):
                return _verdict("journal_window", "verified",
                                f"{len(changed)} change(s)")
            if changed:
                return _verdict("journal_session", "verified",
                                "changed before this task started")
            return _verdict(
                "no_change", "unmet", "nothing written",
                "promises a change and this session changed no file at "
                "all.")

        return _verdict("", "unchecked", "nothing checkable in the subject")
    except Exception:
        return _verdict("", "unchecked", "check failed")




