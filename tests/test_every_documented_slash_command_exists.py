"""Every slash command the help or the docs name is a real command.

Assignment L: the terminal's ``/help`` text is generated from the register
(``BUILTINS`` in ``delfin/agent/repl_commands.py``), so it cannot drift —
but README.md and docs/ name commands in prose, and message texts in
``delfin/agent/*.py`` recommend them ("use /approve ...").  Those three
surfaces can drift apart.  This test pins them to the two real registers:

* terminal: ``BUILTINS`` (imported — the module is PTY-free by design),
* dashboard: ``_SLASH_COMMANDS`` and the routing set
  ``_BUILTIN_SLASH_PREFIXES`` in ``delfin/dashboard/tab_agent.py``,
  extracted textually in the style of ``tests/manual_blocks.py`` (that
  module imports half the dashboard, which a doc test must not).

If a doc names a command neither register has, the DOC is wrong (fix the
doc).  If the code recommends a command that is not registered, the code
is wrong — that is reported to the operator, not fixed here.
"""

from __future__ import annotations

import re
from pathlib import Path

from delfin.agent.repl_commands import BUILTINS
from tests.manual_blocks import repo_root

# ---------------------------------------------------------------------------
# The two real registers
# ---------------------------------------------------------------------------

_TOKEN_RE = re.compile(r"/[a-z][a-z0-9_-]{1,20}")


def _block_after(path: Path, name: str) -> str:
    """The source text of a module-level tuple/set, extractor-style.

    ``name`` is anchored to a line start so it cannot match inside a
    LONGER identifier: ``_SLASH_COMMANDS`` also occurs as a suffix of
    ``_USER_ONLY_SLASH_COMMANDS`` (tab_agent.py), whose comment block
    mentions a ``/token`` that is no command at all.
    """
    text = path.read_text(encoding="utf-8")
    match = re.search(rf"^{re.escape(name)}[^=\n]*= ", text, re.MULTILINE)
    assert match is not None, f"{name} not defined at module level in {path}"
    start = match.start()
    end = text.index("\n)", start)
    return text[start:end]


def dashboard_commands() -> set[str]:
    """Slash command names from the dashboard's palette table."""
    src = _block_after(
        repo_root() / "delfin" / "dashboard" / "tab_agent.py",
        "_SLASH_COMMANDS",
    )
    return set(_TOKEN_RE.findall(src))


def dashboard_routed() -> set[str]:
    """Slash prefixes the dashboard send path actually dispatches."""
    src = _block_after(
        repo_root() / "delfin" / "dashboard" / "tab_agent.py",
        "_BUILTIN_SLASH_PREFIXES",
    )
    return set(_TOKEN_RE.findall(src))


def all_commands() -> set[str]:
    return set(BUILTINS) | dashboard_commands()


# ---------------------------------------------------------------------------
# 1. The registers themselves are sane
# ---------------------------------------------------------------------------

def test_every_builtin_has_a_help_text():
    empty = [name for name, cmd in BUILTINS.items() if not cmd.summary.strip()]
    assert empty == [], f"builtins without a summary: {empty}"


def test_every_dashboard_palette_command_is_routed():
    """A palette entry the send path does not route goes to the model as chat
    (the /undo-file bug tab_agent.py documents at its routing set)."""
    unrouted = {c for c in dashboard_commands()
                if c.split(" ")[0] not in dashboard_routed()}
    assert unrouted == set(), f"in palette but never routed: {sorted(unrouted)}"


# ---------------------------------------------------------------------------
# 2. Documentation: every named command exists
# ---------------------------------------------------------------------------

def _doc_tokens(path: Path, only_backticked: bool) -> set[str]:
    """Slash tokens a doc names in prose, outside every fenced block.

    README separates the terminal builtins paragraph from a fenced,
    explicitly not-runnable dashboard block; only backticked tokens count
    as command references, and fenced blocks are skipped wholesale so
    shell paths inside them cannot read as commands.
    """
    from tests.manual_blocks import blocks_of

    lines = path.read_text(encoding="utf-8").splitlines()
    fenced = set()
    for block in blocks_of("\n".join(lines), str(path)):
        for offset in range(len(block.lines())):
            # block.line is 1-based (the opening fence); body starts after
            fenced.add(block.line + offset)
    tokens: set[str] = set()
    for idx, line in enumerate(lines, start=1):
        if idx in fenced:
            continue
        tokens.update(_tokens_in(line, only_backticked))
    return tokens


# things that look like a slash command but are not one: absolute
# paths (/opt, /usr), file extensions, module names. Applied both to
# doc text and to delfin/agent message strings.
NOT_COMMANDS = re.compile(
    r"/(api|dev|tmp|home|proc|usr|etc|var|opt|lib|run|sys|bin|sbin|root|"
    r"delfin|tests|docs|agent|dashboard|pack|\.py|\.md|\.json|\.txt|"
    r"smiles|converter|stdout|stderr|stdin)\b")


def _tokens_in(text: str, only_backticked: bool) -> set[str]:
    if only_backticked:
        return {t for t in
                re.findall(r"`(/[a-z][a-z0-9_-]{1,20})`", text)
                if not NOT_COMMANDS.search(t)}
    return {t for t in _TOKEN_RE.findall(text)
            if not NOT_COMMANDS.search(t)}


DOC_FILES = [
    "README.md",
    "docs/USER_MANUAL.md",
    "docs/SETTINGS_AND_SETUP.md",
    "docs/RETRY_LOGIC.md",
    "docs/methodology.md",
]


def test_every_backticked_command_in_the_docs_exists():
    known = all_commands() | dashboard_routed()
    missing = []
    for rel in DOC_FILES:
        path = repo_root() / rel
        if not path.exists():
            continue
        for token in sorted(_doc_tokens(path, only_backticked=True)):
            if token not in known:
                missing.append(f"{rel}: {token}")
    assert missing == [], f"documented but not a command: {missing}"


# ---------------------------------------------------------------------------
# 3. Message texts: every recommended command exists
# ---------------------------------------------------------------------------

def _agent_message_recommendations() -> list[str]:
    """``/token`` a delfin/agent/*.py message text recommends.

    A recommendation is a slash token a message tells the user to run:
    it follows a leading verb (``use /try /run /type /enter /call /
    send`` …) and starts after a word boundary. That boundary is
    decisive — a slash directly after a letter is the second half of
    word/word pairs that are no commands (`curl/wget`, `image/jpeg`,
    `energies/properties`, `squeue/sacct`) or a path segment
    (`feedback:/project:`, `git -C /elsewhere` only matches without
    the verb rule). The verb anchor keeps ``ACTION: /command`` and
    ``sed -n '/class A/…'`` docstring examples out too.
    """
    rec = re.compile(
        r"\b(?:use|try|run|type|enter|call|send|write|hit|press|invoke|"
        r"execute|start|switch|choose|pick|see|check)\b[^.\n]{0,60}?[ ]"
        r"(?<![\w/*~>.<-])/([a-z][a-z0-9_-]{1,20})"
        r"(?=[ ,.:;)`'\"]|$)")
    root = repo_root() / "delfin" / "agent"
    found: set[str] = set()
    for path in sorted(root.glob("*.py")):
        text = path.read_text(encoding="utf-8")
        for match in rec.finditer(text):
            token = "/" + match.group(1)
            if NOT_COMMANDS.search(token):
                continue
            found.add(token)
    return sorted(found)


def test_every_recommended_command_in_agent_messages_exists(
        capsys, request):
    recs = _agent_message_recommendations()
    # -s visibility: what the extractor found, so an empty or surprising
    # list is never mistaken for "nothing to check".
    print("\nrecommended commands found:", recs or "(none)")
    assert recs, ("extractor found nothing — the verb-anchored regex "
                  "no longer matches real recommendations; loosen it, "
                  "do not delete this guard")
    known = all_commands() | dashboard_routed()
    missing = [t for t in recs if t not in known]
    assert missing == [], f"recommended in code but not a command: {missing}"
