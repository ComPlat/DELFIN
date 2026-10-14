"""Every name docs/SECURITY_MODEL.md puts in backticks exists in the code.

Input: the backticked spans of ``docs/SECURITY_MODEL.md`` outside fenced
blocks. Output: one assertion per category, listing every name the code
does not have. Semantics, per shape of the span:

* ``--flag``                 -> the quoted string ``"--flag"`` occurs in
  ``delfin/`` (the argparse definition, or the bwrap option the argv uses);
* ``DELFIN_NAME``            -> the name occurs in ``delfin/``;
* ``agent.a.b`` / ``kit.a``  -> every key after the namespace occurs as a
  quoted string in ``delfin/`` (the ``.get("a")`` that reads it);
* ``module.symbol``          -> ``delfin/agent/module.py`` defines
  ``symbol`` (``def``, ``class`` or a module-level assignment);
* ``module``                 -> named as a module where the span is a
  bare identifier and ``delfin/agent/module.py`` exists, it must exist;
* ``delfin-agent sub [sub]`` -> ``add_parser("sub"`` in the CLI;
* a bare identifier in the section "Where a decision is recorded" -> a
  kind in ``security_events._KINDS``;
* ``#NNN``                   -> only inside "Changes in review": an open
  pull request is never described as part of ``main``.

Slash commands in the document are checked by
``tests/test_every_documented_slash_command_exists.py`` (it lists the file).

Rejected: checking pull-request numbers against GitHub -- a test must pass
offline and on every machine; the section rule is what keeps a claim about
``main`` from resting on an unmerged change. Rejected: a whitelist of
names -- it is a second copy of the document that drifts with it.
Cost: one read of ``delfin/**/*.py`` (about 30 MB of text), well under a
second.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "SECURITY_MODEL.md"
AGENT = ROOT / "delfin" / "agent"


def _doc_text() -> str:
    """The document without fenced blocks."""
    out, fenced = [], False
    for line in DOC.read_text(encoding="utf-8").splitlines():
        if line.lstrip().startswith("```"):
            fenced = not fenced
            continue
        if not fenced:
            out.append(line)
    return "\n".join(out)


def _section(title: str) -> str:
    text = _doc_text()
    m = re.search(rf"^## {re.escape(title)}\n(.*?)(?=^## |\Z)", text,
                  re.MULTILINE | re.DOTALL)
    assert m, f"section '## {title}' missing from {DOC.name}"
    return m.group(1)


def _spans(text: str) -> list[str]:
    return re.findall(r"`([^`\n]+)`", text)


_SOURCE: str = ""


def _source() -> str:
    global _SOURCE
    if not _SOURCE:
        _SOURCE = "\n".join(p.read_text(encoding="utf-8", errors="replace")
                            for p in sorted((ROOT / "delfin").rglob("*.py")))
    return _SOURCE


def _defines(module: str, symbol: str) -> bool:
    path = AGENT / f"{module}.py"
    if not path.is_file():
        return False
    text = path.read_text(encoding="utf-8")
    return re.search(
        rf"^\s*(?:def|class)\s+{re.escape(symbol)}\b"
        rf"|^{re.escape(symbol)}\s*(?::[^=\n]*)?=",
        text, re.MULTILINE) is not None


_FLAG = re.compile(r"^--[a-z][a-z0-9-]+$")
_ENV = re.compile(r"^DELFIN_[A-Z0-9_]+$")
_SETTING = re.compile(r"^(agent|kit)\.([a-z_]+(?:\.[a-z_]+)*)$")
_SYMBOL = re.compile(r"^([a-z_]+)\.([A-Za-z_][A-Za-z0-9_]*)$")
_IDENT = re.compile(r"^[a-z][a-z0-9]*(?:_[a-z0-9]+)*$")
#: Dotted names that are not ``module.symbol``: the two settings
#: namespaces, and git's own configuration section (``core.hooksPath``).
_NOT_MODULES = frozenset({"agent", "kit", "core"})
_COMMAND = re.compile(r"^delfin-agent ([a-z-]+)(?: ([a-z-]+))?")


def test_the_document_exists_and_names_things():
    spans = _spans(_doc_text())
    assert len(spans) > 50, "the extractor found almost nothing to check"


def test_every_flag_exists():
    flags = sorted({s for s in _spans(_doc_text()) if _FLAG.match(s)})
    assert len(flags) >= 8, flags
    missing = [f for f in flags if f'"{f}"' not in _source()]
    assert missing == [], f"flags named in {DOC.name} but not in the code: {missing}"


def test_every_environment_variable_exists():
    names = sorted({s for s in _spans(_doc_text()) if _ENV.match(s)})
    assert len(names) >= 4, names
    missing = [n for n in names if n not in _source()]
    assert missing == [], f"variables named but not in the code: {missing}"


def test_every_setting_key_is_read_by_the_code():
    keys = sorted({s for s in _spans(_doc_text()) if _SETTING.match(s)})
    assert len(keys) >= 6, keys
    missing = []
    for key in keys:
        for part in _SETTING.match(key).group(2).split("."):
            if f'"{part}"' not in _source():
                missing.append(f"{key} ({part})")
    assert missing == [], f"settings named but never read: {missing}"


def test_every_module_and_symbol_exists():
    modules = {p.stem for p in AGENT.glob("*.py")}
    checked, missing = 0, []
    for span in sorted(set(_spans(_doc_text()))):
        m = _SYMBOL.match(span)
        if m and m.group(1) in modules:
            checked += 1
            if not _defines(m.group(1), m.group(2)):
                missing.append(span)
        elif m and m.group(1) not in _NOT_MODULES:
            missing.append(f"{span} (no module delfin/agent/{m.group(1)}.py)")
    assert checked >= 15, checked
    assert missing == [], f"symbols named but not defined: {missing}"


def test_every_bare_module_name_exists():
    """A bare identifier is a module where the document says so: the names
    it introduces as the implementation of a layer."""
    named = ("landlock_exec", "seatbelt", "socket_guard", "process_guard",
             "egress_proxy", "approval_guard", "security_events",
             "mcp_isolation", "web_tools", "audit_log")
    spans = set(_spans(_doc_text()))
    for module in named:
        assert (module in spans
                or any(s.startswith(module + ".") for s in spans)), \
            f"{module} is no longer named; update this list"
        assert (AGENT / f"{module}.py").is_file(), module


def test_every_cli_command_exists():
    cli = (AGENT / "cli.py").read_text(encoding="utf-8")
    commands = set()
    for span in _spans(_doc_text()):
        m = _COMMAND.match(span)
        if m:
            commands.update(c for c in m.groups() if c)
    assert {"capabilities", "gate", "explain"} <= commands, commands
    missing = [c for c in sorted(commands)
               if not re.search(rf'add_parser\(\s*"{re.escape(c)}"', cli)]
    assert missing == [], f"commands named but not registered: {missing}"


def test_every_recorded_kind_is_a_security_event():
    from delfin.agent.security_events import _KINDS
    kinds = [s for s in _spans(_section("Where a decision is recorded"))
             if _IDENT.match(s) and s not in ("security_events", "audit_log")]
    assert len(kinds) >= 10, kinds
    missing = [k for k in kinds if k not in _KINDS]
    assert missing == [], f"event kinds named but never recorded: {missing}"
    # Kinds named elsewhere in the text, too.
    elsewhere = {"isolation_missing", "approval_from_outside"}
    assert elsewhere <= set(_KINDS)


@pytest.mark.parametrize("title", ["Changes in review"])
def test_open_pull_requests_are_named_only_as_changes_in_review(title):
    text = _doc_text()
    review = _section(title)
    outside = text.replace(review, "")
    stray = re.findall(r"#\d{2,}", outside)
    assert stray == [], f"pull requests named outside '{title}': {stray}"
    assert re.findall(r"#\d{2,}", review), "the review section names no PR"
    assert "not on `main`" in review
