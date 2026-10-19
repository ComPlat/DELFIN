"""Check a file the agent just edited, by what kind of file it is.

Input: the paths edited in a turn. Output: :func:`problems` -- a short
summary of what is broken ("" when nothing is), for the turn-end
auto-verify, which sends the model back to fix it before it may finish.

What each kind gets, fast and without running the project:

* ``.py`` -- ``py_compile``;
* ``.js`` / ``.mjs`` / ``.cjs`` -- ``node --check`` (when node is installed);
* ``.html`` / ``.htm`` -- every inline ``<script>`` through ``node --check``,
  reported at its line in the HTML; every local ``src``/``href`` of a
  script, stylesheet or image that does not exist; and, when Playwright is
  importable, the page loaded headless with its uncaught errors and
  console errors;
* ``.json`` -- parsed;
* ``.sh`` / ``.bash`` -- ``bash -n``.

A missing checker (no node, no Playwright) is no problem and no finding:
the check is skipped for that kind, never reported as broken.

Measured 2026-10-19 (session 78c36f8c): a browser game in one HTML file,
built over two hours; the turn-end check looked only at ``.py`` files, so
not one of its edits was ever checked by the harness, and the model -- which
could not open a browser -- asked the user to test each one.

Never raises.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import tempfile
from html.parser import HTMLParser
from pathlib import Path
from typing import Optional

CHECKABLE_SUFFIXES = frozenset({
    ".py", ".js", ".mjs", ".cjs", ".html", ".htm", ".json", ".sh", ".bash",
})
_TIMEOUT_S = 20
_MAX_REPORT = 1500


def is_checkable(path: str) -> bool:
    return Path(str(path or "")).suffix.lower() in CHECKABLE_SUFFIXES


def _run(argv: list, cwd: Optional[str] = None) -> tuple[int, str]:
    try:
        proc = subprocess.run(argv, cwd=cwd, capture_output=True, text=True,
                              timeout=_TIMEOUT_S)
        return proc.returncode, (proc.stderr or proc.stdout or "").strip()
    except Exception as exc:
        return 0, f"(not checked: {exc})"


def _node_message(out: str) -> str:
    """node --check's report without its stack trace and version line."""
    keep = []
    for ln in out.splitlines():
        if ln.startswith("    at ") or ln.startswith("Node.js v"):
            continue
        keep.append(ln)
    return re.sub(r"\n{2,}", "\n", "\n".join(keep)).strip()


def _node() -> str:
    return shutil.which("node") or ""


def _check_py(path: Path) -> str:
    import py_compile
    try:
        py_compile.compile(str(path), doraise=True)
    except py_compile.PyCompileError as exc:
        return str(exc).strip()[:400]
    except Exception:
        return ""
    return ""


def _check_js(path: Path) -> str:
    node = _node()
    if not node:
        return ""
    code, out = _run([node, "--check", str(path)])
    return _node_message(out)[:600] if code else ""


def _check_json(path: Path) -> str:
    try:
        json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        return f"{path.name}: invalid JSON: {exc}"
    except Exception:
        return ""
    return ""


def _check_sh(path: Path) -> str:
    bash = shutil.which("bash")
    if not bash:
        return ""
    code, out = _run([bash, "-n", str(path)])
    return out[:600] if code else ""


class _Scan(HTMLParser):
    """Inline scripts with their starting line, and local references."""

    _JS_TYPES = ("", "text/javascript", "application/javascript", "module")

    def __init__(self) -> None:
        super().__init__(convert_charrefs=False)
        self.scripts: list[tuple[int, str, bool]] = []
        self.refs: list[tuple[int, str]] = []
        self._in: Optional[tuple[int, bool]] = None
        self._buf: list[str] = []

    def handle_starttag(self, tag, attrs):
        a = {k.lower(): (v or "") for k, v in attrs}
        line = self.getpos()[0]
        if tag == "script":
            if a.get("src"):
                self.refs.append((line, a["src"]))
            elif a.get("type", "").strip().lower() in self._JS_TYPES:
                self._in = (line, a.get("type", "").strip().lower() == "module")
                self._buf = []
        elif tag == "link" and "stylesheet" in a.get("rel", "").lower() \
                and a.get("href"):
            self.refs.append((line, a["href"]))
        elif tag == "img" and a.get("src"):
            self.refs.append((line, a["src"]))

    def handle_data(self, data):
        if self._in is not None:
            self._buf.append(data)

    def handle_endtag(self, tag):
        if tag == "script" and self._in is not None:
            line, module = self._in
            self.scripts.append((line, "".join(self._buf), module))
            self._in = None


_LINE_RE = re.compile(r":(\d+)\b")


def _check_html(path: Path) -> str:
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return ""
    scan = _Scan()
    try:
        scan.feed(text)
    except Exception:
        return ""
    probs: list[str] = []
    node = _node()
    if node:
        with tempfile.TemporaryDirectory(prefix="delfin-check-") as tmp:
            for n, (line, body, module) in enumerate(scan.scripts):
                if not body.strip():
                    continue
                f = Path(tmp) / f"s{n}{'.mjs' if module else '.js'}"
                # The script's first line sits on the <script> tag's line:
                # padding keeps node's line numbers those of the HTML.
                f.write_text("\n" * (line - 1) + body, encoding="utf-8")
                code, out = _run([node, "--check", str(f)])
                if code:
                    out = _node_message(out.replace(str(f), path.name))
                    probs.append(f"{path.name}: inline <script> at line "
                                 f"{line}: {out[:500]}")
    for line, ref in scan.refs:
        r = ref.strip()
        if (not r or re.match(r"^[a-z][a-z0-9+.-]*:", r, re.I)
                or r.startswith(("//", "#", "{", "$"))):
            continue
        local = r.split("#", 1)[0].split("?", 1)[0]
        target = (path.parent / local) if not local.startswith("/") else None
        if target is not None and local and not target.exists():
            probs.append(f"{path.name}:{line}: references {r}, which does "
                         "not exist")
    if not probs:
        probs.extend(_browser_errors(path))
    return "\n".join(probs)


def _browser_errors(path: Path) -> list[str]:
    """Uncaught errors and console errors of the page loaded headless, when
    Playwright and its browser are installed; [] otherwise."""
    try:
        from playwright.sync_api import sync_playwright  # type: ignore
    except Exception:
        return []
    errors: list[str] = []
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(headless=True)
            try:
                page = browser.new_page()
                page.on("pageerror", lambda e: errors.append(f"uncaught: {e}"))
                page.on("console", lambda m: errors.append(
                    f"console.error: {m.text}") if m.type == "error" else None)
                page.goto(path.resolve().as_uri(), wait_until="load",
                          timeout=_TIMEOUT_S * 1000)
                page.wait_for_timeout(1500)
            finally:
                browser.close()
    except Exception:
        return []
    return [f"{path.name} in a headless browser: {e}"[:400]
            for e in errors[:5]]


_BY_SUFFIX = {
    ".py": _check_py, ".js": _check_js, ".mjs": _check_js, ".cjs": _check_js,
    ".html": _check_html, ".htm": _check_html, ".json": _check_json,
    ".sh": _check_sh, ".bash": _check_sh,
}


def problems(paths) -> str:
    """What is broken in *paths*, one block per file; "" when nothing is."""
    out: list[str] = []
    for p in paths or ():
        try:
            path = Path(str(p))
            fn = _BY_SUFFIX.get(path.suffix.lower())
            if fn is None or not path.is_file():
                continue
            found = fn(path)
            if found:
                out.append(found)
        except Exception:
            continue
    text = "\n".join(out)
    return text if len(text) <= _MAX_REPORT else text[:_MAX_REPORT] + " …"


__all__ = ["CHECKABLE_SUFFIXES", "is_checkable", "problems"]
