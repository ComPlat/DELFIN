"""Diagnostics after an edit (LSP1).

``check(path, workspace)`` returns the language-level findings for ONE
just-changed file, so the agent sees an undefined name or an unused
import the moment the write lands — not one PR cycle later.

Python: ``ruff check --output-format json`` (the ruff on PATH; the CI
pin is in .github/ci/constraints-ci.txt). ``pyright --outputjson`` is
only consulted when it is installed. Nothing is installed or downloaded
here. Every other language yields an empty list — said once, not per
call.

Bounded by contract: one file, at most 5 s wall per tool, at most 20
diagnostics, never a network call, never raises. A missing tool or a
timeout collapses into a one-line unavailability note on the result
object, not into an exception that could break the write path that
calls this.

Only NEW findings are reported: findings present in the pre-edit text
(``old_text``, journalled by the write path) are filtered out, so an
edit that touches one line does not re-report the file's legacy debt.
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

MAX_SECONDS_PER_TOOL = 5.0
MAX_DIAGNOSTICS = 20
#: 10 MB plus slack so a real ruff run never trips the cap, but a
#: hostile flood is cut before json parsing (and before memory grows).
_MAX_OUTPUT_BYTES = 10 * 1024 * 1024 + 65536

PYTHON_EXTS = {".py", ".pyi"}

#: Ruff rules worth surfacing after an edit — the wave 13-18 findings
#: (undefined names, wrong/unused imports) plus the obvious syntax and
#: hygiene breaks. Kept explicit so a config change in the user's repo
#: cannot flood the edit result.
_SELECT = "F,E999"


@dataclass
class Diagnostic:
    line: int
    column: int
    code: str
    message: str
    severity: str  # "error" | "warning" | "info"

    def to_dict(self) -> dict:
        return {
            "line": self.line, "column": self.column, "code": self.code,
            "message": self.message, "severity": self.severity,
        }


@dataclass
class CheckResult:
    diagnostics: list
    note: str  # "" when available; else "diagnostics unavailable: <why>"

    def to_dict(self) -> dict:
        return {"diagnostics": [d.to_dict() for d in self.diagnostics],
                "note": self.note}


#: ANSI escape sequences never belong in a diagnostic message.
_ANSI_RE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")

#: How many plausible JSON start positions the extractor tries before it
#: gives up. Real ruff/pyright output starts with the value itself (or a
#: banner line before it), so a handful of candidates is generous; a
#: binary splat fails fast here instead of an O(n^2) tail-walk.
_JSON_START_CANDIDATES = 8


def _strip_ansi(text: str) -> str:
    return _ANSI_RE.sub("", text)


def _extract_json_array(text: str):
    """Parse the FIRST JSON value in possibly noisy stdout (hostile or
    chatty tool output); raises ValueError when none parses.

    Bounded, deliberately: only the first _JSON_START_CANDIDATES
    plausible start positions are tried, and raw_decode never rescans a
    failing position's whole tail — a 10 MB binary splat is rejected in
    milliseconds instead of the O(n^2) full-tail walk a naive loop does
    (reviewer s15, 2026-10-18: real timeout found by adversarial run)."""
    text = _strip_ansi(text)
    dec = json.JSONDecoder()
    tried = 0
    for start in range(len(text)):
        if tried >= _JSON_START_CANDIDATES:
            break
        ch = text[start]
        if ch not in "[{":
            continue
        tried += 1
        try:
            value, _end = dec.raw_decode(text, start)
        except ValueError:
            continue
        return value
    raise ValueError("no JSON value found")


def _sev(code: str) -> str:
    return "error" if code.startswith("E999") or code.startswith("F821") else "warning"


def _is_inside(path: Path, workspace: Path) -> bool:
    try:
        path.resolve().relative_to(Path(workspace).resolve())
        return True
    except (ValueError, OSError):
        return False


def _run(cmd: list, cwd: Path) -> tuple[Optional[str], str]:
    """Run one bounded tool. Returns (stdout_text_or_None, note)."""
    if shutil.which(cmd[0]) is None:
        return None, f"diagnostics unavailable: {cmd[0]} not found"
    try:
        proc = subprocess.run(
            cmd, cwd=str(cwd), capture_output=True, timeout=MAX_SECONDS_PER_TOOL)
    except subprocess.TimeoutExpired:
        return None, f"diagnostics unavailable: {cmd[0]} timed out after {MAX_SECONDS_PER_TOOL:g}s"
    except OSError as exc:
        return None, f"diagnostics unavailable: {cmd[0]} failed: {exc}"
    if len(proc.stdout) > _MAX_OUTPUT_BYTES:
        return None, f"diagnostics unavailable: {cmd[0]} output exceeded 10 MB"
    return proc.stdout.decode("utf-8", errors="replace")[:_MAX_OUTPUT_BYTES], ""


def _ruff_findings(path: Path, workspace: Path) -> tuple[list, str]:
    stdout, note = _run(
        ["ruff", "check", "--select", _SELECT, "--output-format", "json",
         str(path)],
        workspace)
    if stdout is None:
        return [], note
    try:
        items = _extract_json_array(stdout or "[]")
        if not isinstance(items, list):
            raise ValueError("ruff output is not a JSON array")
    except (ValueError, TypeError):
        return [], "diagnostics unavailable: ruff output was not valid JSON"
    out = []
    for item in items:
        try:
            loc = item.get("location", {})
            code = str(item.get("code", "") or "")
            out.append(Diagnostic(
                line=int(loc.get("row", 0)), column=int(loc.get("column", 0)),
                code=code, message=_strip_ansi(str(item.get("message", ""))),
                severity=_sev(code),
            ))
        except (AttributeError, TypeError, ValueError):
            continue
        if len(out) >= MAX_DIAGNOSTICS:
            break
    return out, ""


def _pyright_findings(path: Path, workspace: Path) -> tuple[list, str]:
    # Only when installed; a missing pyright is silently skipped (ruff
    # already reported, and the note must not blame a tool nobody has).
    if shutil.which("pyright") is None:
        return [], ""
    stdout, note = _run(
        ["pyright", "--outputjson", str(path)], workspace)
    if stdout is None:
        return [], note
    try:
        data = _extract_json_array(stdout or "{}")
        if not isinstance(data, dict):
            raise ValueError("pyright output shape unexpected")
        items = data.get("generalDiagnostics", [])
        if not isinstance(items, list):
            raise ValueError("pyright output shape unexpected")
    except (ValueError, TypeError):
        return [], "diagnostics unavailable: pyright output was not valid JSON"
    sevmap = {"error": "error", "warning": "warning", "information": "info"}
    out = []
    for item in items:
        try:
            rng = (item.get("range", {}).get("start", {})
                   if isinstance(item.get("range"), dict) else {})
            code = str(item.get("rule", "") or "")
            out.append(Diagnostic(
                line=int(rng.get("line", 0)) + 1,
                column=int(rng.get("character", 0)) + 1,
                code=code, message=_strip_ansi(str(item.get("message", ""))),
                severity=sevmap.get(str(item.get("severity", "warning")),
                                    "warning"),
            ))
        except (AttributeError, TypeError, ValueError):
            continue
        if len(out) >= MAX_DIAGNOSTICS:
            break
    return out, ""


def _key_of(d: Diagnostic) -> tuple:
    return (d.code, d.message)


def check(path, workspace, old_text: Optional[str] = None) -> list:
    """Findings for one just-changed file, minus anything the pre-edit
    text already carried. Returns a plain list of Diagnostic (the brief's
    contract); unavailability and internal errors collapse into an empty
    list — the wiring uses ``check_detailed`` when it needs the one-line
    "diagnostics unavailable: <why>" note. Never raises."""
    return check_detailed(path, workspace, old_text).diagnostics


#: Session-wide flag for the operator ruling (2026-10-18): a missing
#: tool or a timeout is shown as the one-line "diagnostics unavailable:
#: <why>" note ONCE per session, then the block stays quiet — an edit
#: stream must not repeat the same bad-news line after every write.
_unavailability_shown = False


def reset_unavailability_note() -> None:
    """Test/turn-boundary hook: show the note again."""
    global _unavailability_shown
    _unavailability_shown = False


def _take_unavailability(note: str) -> str:
    """Return the note the first time it occurs this session, else ''."""
    global _unavailability_shown
    if not note:
        return ""
    if _unavailability_shown:
        return ""
    _unavailability_shown = True
    return note


def check_detailed(path, workspace,
                   old_text: Optional[str] = None) -> CheckResult:
    """Like :func:`check`, but carries the unavailability note (shown
    once per session, see _take_unavailability)."""
    if workspace is None:
        # No workspace known (e.g. perms without one): never check
        # anything, never raise — the edit result stays untouched.
        return CheckResult([], "")
    try:
        path = Path(path)
        workspace = Path(workspace)
        if path.suffix.lower() not in PYTHON_EXTS:
            return CheckResult([], "")
        if not path.is_file():
            return CheckResult([], "")
        if not _is_inside(path, workspace):
            return CheckResult(
                [], _take_unavailability(
                    "diagnostics unavailable: file outside workspace"))
        findings, note = _ruff_findings(path, workspace)
        pyright_items, pyright_note = _pyright_findings(path, workspace)
        seen = {}
        for d in findings + pyright_items:
            seen.setdefault(_key_of(d), d)
        if old_text:
            old_codes = _old_findings(path, workspace, old_text)
            findings = [
                d for d in seen.values() if _key_of(d) not in old_codes]
        else:
            findings = list(seen.values())
        final_note = _take_unavailability(note or pyright_note)
        if len(findings) > MAX_DIAGNOSTICS:
            findings = findings[:MAX_DIAGNOSTICS]
        return CheckResult(findings, final_note)
    except Exception as exc:  # never break the write path
        return CheckResult([], f"diagnostics unavailable: internal error {type(exc).__name__}")


def _old_findings(path: Path, workspace: Path, old_text: str) -> set:
    """Codes+messages the PRE-EDIT text already carried. The pre-image
    goes into a throwaway copy next to the real file so ruff sees
    identical project context; we clean it up ourselves. The copy is
    dot-prefixed (hidden from lint/glob) and a write failure (read-only
    dir, unencodable text) degrades to "nothing pre-existing" — the
    edit result then shows the file's full findings rather than none.

    Filtering key: (code, message), WITHOUT the line number — a finding
    that merely moved lines counts as pre-existing. Deliberate: line
    numbers shift on every insert above the finding, and re-reporting
    the same defect because it slid down two lines is noise, not
    signal."""
    import os

    tmp = path.with_name("." + path.name + ".lsp1-old.py")
    try:
        tmp.write_text(old_text, encoding="utf-8")
    except Exception:
        return set()
    findings, note = _ruff_findings(tmp, workspace)
    try:
        os.unlink(tmp)
    except OSError:
        pass
    return {_key_of(d) for d in findings}
