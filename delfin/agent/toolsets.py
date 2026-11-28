"""Optional tool groups and the workspace facts that decide them.

Two groups of advertised tools apply to a minority of workspaces:

* ``chemistry`` -- every tool of the two built-in MCP servers
  (``mcp__delfin-ops__*``: ORCA/xTB output parsers, calc status, runtime
  checks; ``mcp__delfin-tools__*``: the DELFIN tools framework) plus
  ``show_molecule``. Measured on a TypeScript repository: 43 schemas,
  6,661 of 16,456 tool-schema tokens per request.
* ``office`` -- the twelve document tools (``read_document``,
  ``edit_sheet``, ``create_docx`` ...): 2,028 tokens.

Input: the session's workspace, whether it is a DELFIN tree, the
directories the session carries (a calc folder is a confirm-write dir),
the ``agent.toolsets`` setting and the role. Output: per group, offered
or not, with the reason. A group is offered when

* the setting says ``on``, or
* the setting says ``auto`` (the default) and the workspace applies:
  chemistry for a DELFIN tree, a session with a calc/archive dir, or
  chemistry files within three directory levels (``.xyz``, ``.gbw``,
  ``.hess``, ``CONTROL.txt``, an ORCA/xTB/Gaussian ``.out``/``.inp``/
  ``.log``); office for a DELFIN tree, the office role, or document
  files (``.xlsx``, ``.docx``, ``.pdf`` ...), or
* the user enabled it (``/tools enable <group>``), or
* the session already called one of its tools (a resumed history
  included), or
* the person's request NAMES the group (``request_groups``: "make me a
  docx", "the LUMO of this molecule"). Read once, from the first turn's
  own words, so a later turn cannot move a head the prefix cache has
  already served.

Semantics, chosen for the prefix cache: the workspace facts are read
ONCE, at the first query (the first turn, after any benchmark setup has
seeded the workspace), and never re-read. A group only ever goes from
hidden to offered, never back, and the caller applies a change at the
next turn boundary (the tool list is assembled once per turn). A file
the agent writes mid-session therefore changes nothing; enabling costs
one cold prefix, once.

Hiding is advertising only. The executor's authorisation is unchanged:
a hidden tool is exactly as allowed as before, so a call that reaches it
runs, and enables its group from the next turn on.

Rejected alternatives: re-detecting per turn (the head would change
whenever the agent writes an ``.xyz``); a GENERAL text classifier on the
task -- one that reads the whole request for topical hints -- because
``energy``, ``scan`` and ``charge`` fire on ordinary code tasks; not
spawning the hidden servers (the tool index of the ops server is
still what ``list_tools`` reads once the group is enabled, and the spawn
cost is paid at startup either way).

Not in that rejected list, and added deliberately: a vocabulary of terms
that NAME a group, matched at a word start (``request_groups``). It is
narrower than a classifier by construction -- every term names a document
format, a document kind or a chemistry method, and the tests pin both
directions (every term has a request that needs it; a request that needs
neither group stays without one). It decides only the first turn, from the
person's own words, so the head cannot move under the prefix cache.
"""

from __future__ import annotations

import os
import re
import threading
from pathlib import Path
from typing import Any, Iterable, Optional

CHEMISTRY = "chemistry"
OFFICE = "office"
GROUPS: tuple[str, ...] = (CHEMISTRY, OFFICE)

_CHEMISTRY_MCP_PREFIXES: tuple[str, ...] = (
    "mcp__delfin-ops__", "mcp__delfin-tools__",
)
_CHEMISTRY_NATIVE: frozenset[str] = frozenset({"show_molecule"})

# What the one-line note names, per group. Short on purpose: the note is
# the price a foreign workspace still pays for the hidden groups.
_GROUP_LABEL: dict[str, str] = {
    CHEMISTRY: ("chemistry tools (ORCA/xTB parsers, calc status, DELFIN "
                "tools framework)"),
    OFFICE: "office tools (read_document, edit_sheet, create_docx ...)",
}

_SETTING_VALUES = ("auto", "on", "off")

# Workspace scan bounds. Depth 3 reaches chem/<molecule>/start.xyz; the
# entry cap keeps a session started in a large tree from paying for it.
_SCAN_MAX_DEPTH = 3
_SCAN_MAX_ENTRIES = 3000
_SNIFF_MAX_FILES = 12
_SNIFF_BYTES = 8192
_SKIP_DIRS = frozenset({
    "node_modules", "venv", "env", "target", "build", "dist",
    "__pycache__", "site-packages", "vendor", "bin", "obj",
})
_CHEMISTRY_DIRS = frozenset({"calc", "archive", "remote_archive"})
_CHEMISTRY_EXTS = frozenset({
    ".xyz", ".gbw", ".hess", ".engrad", ".molden", ".cube", ".cub",
    ".mol2", ".sdf", ".gjf",
})
_CHEMISTRY_NAMES = frozenset({"CONTROL.txt", "control.txt", "coord"})
_SNIFF_EXTS = frozenset({".out", ".inp", ".log"})
_CHEMISTRY_SNIFF_RE = re.compile(
    r"O\s{2,}R\s{2,}C\s{2,}A|\bx T B\b|(?i:\bxtb\b)|Entering Gaussian"
    r"|TURBOMOLE|FINAL SINGLE POINT ENERGY|HOMO-LUMO GAP"
    r"|^\s*\*\s*xyzfile?\b|^\s*%(?:pal|maxcore|scf)\b",
    re.MULTILINE)
_OFFICE_EXTS = frozenset({
    ".xlsx", ".xlsm", ".xls", ".ods", ".docx", ".doc", ".odt", ".pdf",
})


def _office_tool_names() -> frozenset[str]:
    from .api_client import _OFFICE_TOOL_NAMES
    return _OFFICE_TOOL_NAMES


def tool_group(name: str) -> Optional[str]:
    """The optional group *name* belongs to, or None for a core tool."""
    name = str(name or "")
    if name.startswith(_CHEMISTRY_MCP_PREFIXES) or name in _CHEMISTRY_NATIVE:
        return CHEMISTRY
    if name in _office_tool_names():
        return OFFICE
    return None


def _setting(group: str) -> str:
    try:
        from delfin.user_settings import load_settings
        cfg = ((load_settings() or {}).get("agent") or {}).get("toolsets")
        value = str((cfg or {}).get(group, "auto") or "auto").strip().lower()
    except Exception:
        value = "auto"
    return value if value in _SETTING_VALUES else "auto"


def scan_workspace(workspace: Optional[Path]) -> dict[str, str]:
    """``{group: evidence}`` for the groups the workspace's files call for.

    Bounded walk: depth <= 3, at most 3000 entries, hidden directories and
    build/dependency directories skipped, at most twelve ``.out``/``.inp``/
    ``.log`` files sniffed (first 8 KiB each). Evidence is the relative
    path that decided it.
    """
    found: dict[str, str] = {}
    if workspace is None:
        return found
    try:
        root = Path(workspace)
        if not root.is_dir():
            return found
    except Exception:
        return found
    seen = 0
    sniffed = 0
    stack: list[tuple[Path, int]] = [(root, 0)]
    while stack and seen < _SCAN_MAX_ENTRIES and len(found) < len(GROUPS):
        here, depth = stack.pop()
        try:
            entries = sorted(os.scandir(here), key=lambda e: e.name)
        except OSError:
            continue
        for entry in entries:
            seen += 1
            if seen > _SCAN_MAX_ENTRIES:
                break
            name = entry.name
            if name.startswith("."):
                continue
            try:
                is_dir = entry.is_dir(follow_symlinks=False)
            except OSError:
                continue
            rel = os.path.relpath(entry.path, root)
            if is_dir:
                if depth == 0 and name in _CHEMISTRY_DIRS:
                    found.setdefault(CHEMISTRY, rel + "/")
                if name not in _SKIP_DIRS and depth + 1 < _SCAN_MAX_DEPTH:
                    stack.append((Path(entry.path), depth + 1))
                continue
            ext = os.path.splitext(name)[1].lower()
            if CHEMISTRY not in found:
                if ext in _CHEMISTRY_EXTS or name in _CHEMISTRY_NAMES:
                    found[CHEMISTRY] = rel
                elif ext in _SNIFF_EXTS and sniffed < _SNIFF_MAX_FILES:
                    sniffed += 1
                    try:
                        with open(entry.path, "rb") as fh:
                            head = fh.read(_SNIFF_BYTES).decode(
                                "utf-8", "replace")
                    except OSError:
                        head = ""
                    if _CHEMISTRY_SNIFF_RE.search(head):
                        found[CHEMISTRY] = rel
            if OFFICE not in found and ext in _OFFICE_EXTS:
                found[OFFICE] = rel
    return found


#: What a person's own words ask for, per group.
#:
#: The workspace scan above answers "does this FOLDER hold document or
#: chemistry files". It cannot answer "what did the person just ask for",
#: and the measured consequence was that the plain request "make me a
#: two-page docx summary of this text" in a folder holding no document
#: got no document tool at all -- 65 tools offered, none of them able to
#: write a .docx, and a system prompt saying
#: ``Not loaded here: office tools (...). If the task needs them, ask the
#: user for /tools enable office``. The model's best move on an ordinary
#: request was to ask for a switch instead of doing the work. Measured
#: over 11 (role x workspace x request) rows in scratch/offer_matrix.md.
#:
#: This vocabulary is deliberately TIGHTER than the prompt-module
#: triggers in ``prompt_loader._MODULE_TRIGGERS``, which match the same
#: subject matter. Those only decide which prose survives stripping; this
#: decides whether 12 (office) or 43 (chemistry) tool schemas ride on
#: every request of the session, so it carries the false positives this
#: module's own docstring already refuses: ``energy``, ``scan``,
#: ``charge`` fire on ordinary code tasks, and so does ``calc/`` -- CSS
#: ``calc(100% - 20px)`` is one of the commonest strings in a frontend
#: repository. ``excel`` is out because ``excellent`` contains it;
#: ``complex``/``komplex`` because they are programming words;
#: ``homo`` because ``homogeneous`` contains it (``lumo`` has no such
#: trap); ``rechnung`` and ``formular`` because German uses them for
#: "calculation" and "HTML form". What is left names a document format,
#: a document kind, or a calculation method -- words that mean the group.
#: Tests pin both directions: every term has a request that needs it, and
#: every request that does not need a group stays without one.
_REQUEST_TERMS: dict[str, tuple[str, ...]] = {
    CHEMISTRY: (
        # Programs and methods, spelled as only a chemistry task spells them
        "orca", "xtb", "turbomole", "gaussian", "dft", "ccsd", "casscf",
        "nevpt2", "b3lyp", "pbe0", "def2", "molden", ".xyz", ".hess",
        "hessian",
        # Results a chemist asks for
        "lumo", "orbital", "mulliken", "loewdin", "thermochem",
        "vibrational", "gibbs", "enthalp", "hartree", "conformer",
        # What they ask the program to do
        "geometry optimi", "transition state", "basis set", "solvent",
        "quantum chem", "molecule",
        # German
        "geometrie", "schwingung", "frequenz", "übergangszustand",
        "basissatz", "lösungsmittel", "quantenchem", "molekül",
        "frequenzrechnung", "einzelpunkt",
    ),
    OFFICE: (
        "docx", "doc-datei", "xlsx", "xlsm", "ods", "odt", "pdf",
        "spreadsheet", "tabellenkalkulation", "arbeitsmappe",
        "word doc", "text document", "pdf-datei",
        "invoice", "letterhead", "cover letter", "mail merge",
        "serienbrief", "anschreiben",
        # The German worksheet and bookkeeping nouns. A German-speaking
        # office user names the THING -- a user's own words, "ein
        # Tabellenblatt anlegen" and "die Belege sortieren" -- and never
        # names a file extension, so without these the group could only be
        # offered by a folder that already held a document, never on the
        # request that creates one. Measured on the column-split request
        # below: the prompt module fired on it, the tool group did not, so
        # a solo session was offered 0 of 12 document tools and no
        # checking rules while its own prompt said the module was active.
        #
        # Each is a COMPOUND or a bookkeeping word, which is what keeps
        # this tighter than the prompt-module triggers while covering
        # the same speech. The bare ``tabelle``/``spalte`` stay OUT: in
        # German a database table and its columns are "Tabelle" and
        # "Spalte" too, so they fire on a SQL request, and ``journal``
        # is a log file. Find the compound instead -- it is the word
        # that can only mean the document.
        "tabellenblatt", "arbeitsblatt", "arbeitsmappe",
        "beleg", "buchung", "kostenstelle", "kontoauszug", "quittung",
        "mahnung", "gutschrift", "lieferschein", "umsatz",
        "zeile eintragen",
        # English, and unambiguous: a worksheet is a document whether it
        # is punched into a spreadsheet or handed out in a classroom, and
        # no code request means one.
        "worksheet",
    ),
}


def _term_in(term: str, s: str) -> bool:
    """Whether *term* occurs in *s* as a word, not inside another one.

    A bare substring match is right for most of the vocabulary and wrong for
    a short format token: ``"ods"`` is inside ``"methods"`` and ``"periods"``,
    so the plainest code request ("refactor these methods") was offered the
    office group for the whole session -- exactly the trap this module's
    vocabulary refuses elsewhere (``excel`` for ``excellent``, ``homo`` for
    ``homogeneous``).

    Only the LEADING edge is bounded. Half the vocabulary is a deliberate
    prefix of what a person types -- ``enthalp`` for enthalpy/enthalpie,
    ``geometry optimi`` for optimise/optimisation, ``thermochem`` for
    thermochemistry -- and requiring a word end would stop every one of them
    from firing. What a leading boundary forbids is a term that begins
    *inside* a longer word, which is the only way a 3-letter format can
    appear by accident. ``.xyz``/``.hess`` need no boundary: the dot
    already supplies one.
    """
    if term[:1] == ".":
        return term in s
    return re.search(rf"(?<![0-9a-z]){re.escape(term)}", s) is not None


def request_groups(text: Any) -> dict[str, str]:
    """``{group: the term that decided it}`` for what *text* asks for.

    Input: the person's own words (the turn's request, or that plus the
    recent user turns). Output: at most one entry per group, with the term
    that fired, so the reason can be shown by ``/tools groups``. Empty for
    empty or non-text input, and for a request that names neither group --
    a code task pays nothing.
    """
    s = str(text or "").lower()
    if not s.strip():
        return {}
    out: dict[str, str] = {}
    for group, terms in _REQUEST_TERMS.items():
        for term in terms:
            if _term_in(term, s):
                out[group] = term
                break
    return out


def tool_names_in_messages(messages: Iterable[Any]) -> set[str]:
    """Names of every tool call recorded in a chat history.

    Reads the OpenAI shape (``tool_calls[].function.name``) and the
    content-block shape (``{"type": "tool_use", "name": ...}``).
    """
    names: set[str] = set()
    for msg in messages or ():
        if not isinstance(msg, dict):
            continue
        for call in msg.get("tool_calls") or ():
            try:
                fn = call.get("function") if isinstance(call, dict) else None
                name = (fn or {}).get("name") if isinstance(fn, dict) else None
                if name:
                    names.add(str(name))
            except Exception:
                continue
        content = msg.get("content")
        if isinstance(content, list):
            for block in content:
                if (isinstance(block, dict)
                        and block.get("type") == "tool_use"
                        and block.get("name")):
                    names.add(str(block["name"]))
    return names


class WorkspaceToolsets:
    """The per-session decision which optional groups are offered.

    Shared by reference between a session's permissions, its sub-agents
    (``dataclasses.replace`` copies the reference) and its prompt loader,
    so the tool list and the prompt text that describes those tools
    cannot disagree.
    """

    def __init__(self, workspace: Optional[Path] = None, *,
                 is_delfin_workspace: bool = False,
                 calc_dirs: Iterable[Path] = (),
                 session_dirs: Iterable[Path] = (),
                 settings: Optional[dict[str, str]] = None) -> None:
        self.workspace = Path(workspace) if workspace else None
        self.is_delfin_workspace = bool(is_delfin_workspace)
        # Confirm-write dirs: the mechanism exists for calc folders (the
        # dashboard hands its calc dir in this way and nothing else).
        self.calc_dirs = tuple(Path(p) for p in calc_dirs or ())
        self.session_dirs = tuple(Path(p) for p in session_dirs or ())
        self._settings_override = dict(settings) if settings else None
        self._lock = threading.Lock()
        self._decided: Optional[dict[str, str]] = None   # group -> reason
        self._enabled_later: dict[str, str] = {}
        self._from_request: dict[str, str] = {}
        self._stacks: Optional[tuple[str, ...]] = None

    # -- decision ---------------------------------------------------------

    def _setting(self, group: str) -> str:
        if self._settings_override is not None:
            value = str(self._settings_override.get(group, "auto")).lower()
            return value if value in _SETTING_VALUES else "auto"
        return _setting(group)

    def _decide(self) -> dict[str, str]:
        offered: dict[str, str] = {}
        settings = {g: self._setting(g) for g in GROUPS}
        scan: Optional[dict[str, str]] = None
        for group in GROUPS:
            mode = settings[group]
            if mode == "on":
                offered[group] = "setting agent.toolsets." + group + "=on"
                continue
            if mode == "off":
                continue
            if self.is_delfin_workspace:
                offered[group] = "DELFIN workspace"
                continue
            if group == CHEMISTRY:
                if self.calc_dirs:
                    offered[group] = (
                        f"calc directory {self.calc_dirs[0].name}/")
                    continue
                named = [p for p in self.session_dirs
                         if p.name in _CHEMISTRY_DIRS]
                if named:
                    offered[group] = f"session directory {named[0].name}/"
                    continue
            if scan is None:
                scan = scan_workspace(self.workspace)
            if group in scan:
                offered[group] = f"workspace file {scan[group]}"
        return offered

    def decided(self) -> dict[str, str]:
        """``{group: reason}`` for every offered group (workspace decides once).

        Workspace facts are decided once per session; groups enabled since
        (by a tool the history already used, by the request's own words, or
        by ``/tools enable``) are merged on top, never withdrawn.
        """
        with self._lock:
            if self._decided is None:
                self._decided = self._decide()
            return {**self._decided, **self._from_request,
                    **self._enabled_later}

    def offered(self, group: str, *, role: str = "") -> bool:
        if group == OFFICE and role == "office_agent":
            return True
        return group in self.decided()

    def hidden_groups(self, *, role: str = "") -> frozenset[str]:
        return frozenset(g for g in GROUPS if not self.offered(g, role=role))

    def enable(self, group: str, reason: str = "enabled by the user") -> bool:
        """Offer *group* from the next turn on. True when that is a change."""
        if group not in GROUPS:
            raise ValueError(f"unknown tool group {group!r}")
        current = self.decided()
        if group in current:
            return False
        with self._lock:
            self._enabled_later[group] = reason
        return True

    def note_used(self, names: Iterable[str]) -> list[str]:
        """Enable every group one of *names* belongs to. Returns the new."""
        added: list[str] = []
        current = None
        for name in names or ():
            group = tool_group(name)
            if group is None:
                continue
            if current is None:
                current = self.decided()
            if group not in current and group not in added:
                with self._lock:
                    self._enabled_later[group] = f"tool {name} already used"
                added.append(group)
        return added

    def note_request(self, text: Any) -> list[str]:
        """Offer every group the person's own words ask for. New groups.

        Input: the turn's request text (with the recent user turns, as the
        prompt assembly holds them). Output: the groups this added, which
        were not offered before.

        TAKEN AT THE FIRST TURN, like every other fact the offering is built
        from, and for the same measured reason. The offering decides stable
        prompt rows -- the lazy module that survives stripping and the
        "not loaded here" line -- so a request on a LATER turn that named a
        new group would rewrite those rows mid-session. That is exactly the
        cold prefix the module's "Semantics" section refuses: an earlier
        reading of this method folded the request in at every turn, and
        tests/test_the_head_does_not_move_with_the_workspace.py caught it
        (role_prompt row 6 moved; the prefix cache cannot reach it).

        The cost is bounded and stated: a session whose FIRST request names
        neither group keeps the `/tools enable <group>` route, which is the
        turn-boundary switch this module already has. What the package asked
        for -- an office request in an empty folder being answered with the
        document tools instead of a request for a switch -- is a first-turn
        case, and is offered.

        The one thing that outranks a request is the person's own switch:
        ``agent.toolsets.<group>: off`` is honoured and the request is
        ignored, because that setting is somebody having already weighed
        exactly this cost. A group already offered for the workspace
        changes nothing here, so the usual case costs one substring scan.
        """
        with self._lock:
            if self._decided is not None:
                # The session's head is composed; a group added now would
                # move rows a provider prefix cache has already served.
                return []
        wanted = request_groups(text)
        if not wanted:
            return []
        added: list[str] = []
        current = self.decided()
        for group, term in wanted.items():
            if group in current or group in self._from_request:
                continue
            if self._setting(group) == "off":
                continue
            with self._lock:
                self._from_request[group] = f"the request says {term!r}"
            added.append(group)
        return added

    # -- stacks (prompt advice) --------------------------------------------

    def stacks(self) -> tuple[str, ...]:
        """The workspace's language stacks, detected once."""
        with self._lock:
            if self._stacks is None:
                try:
                    from .project_init import detect_stacks
                    self._stacks = tuple(detect_stacks(self.workspace))
                except Exception:
                    self._stacks = ()
            return self._stacks

    # -- text --------------------------------------------------------------

    def hidden_note(self, *, role: str = "") -> str:
        """One line naming the hidden groups and how they are enabled."""
        hidden = [g for g in GROUPS if g in self.hidden_groups(role=role)]
        if not hidden:
            return ""
        labels = "; ".join(_GROUP_LABEL[g] for g in hidden)
        cmds = " / ".join(f"`/tools enable {g}`" for g in hidden)
        return (f"Not loaded here: {labels}. If the task needs them, ask "
                f"the user for {cmds} (offered from the next turn).")


def drop_hidden(tools: list[dict[str, Any]],
                hidden: frozenset[str]) -> list[dict[str, Any]]:
    """*tools* without the entries of a hidden group (order kept)."""
    if not hidden:
        return tools
    return [t for t in tools
            if tool_group(t.get("function", {}).get("name", "")) not in hidden]
