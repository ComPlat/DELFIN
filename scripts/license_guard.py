#!/usr/bin/env python3
# license-guard: allow  (this file defines the forbidden patterns; it is self-exempt)
"""CCDC/CSD License Guard — mechanical, fail-closed enforcement.

Purpose
-------
Make it *structurally impossible* to ship CCDC/CSD-derived material in this
public repository. This is the single source of truth for the forbidden
patterns, consumed by three independent layers (defense in depth):

  * tests/test_license_guard.py   -> pytest / CI scans the whole tracked tree
  * .git/hooks/pre-commit         -> blocks CCDC paths + imports at commit time
  * CLI (this file)               -> manual audit

Policy (see LICENSE_GUARD.md)
-----------------------------
CCDC / CSD data is licensed and may be used ONLY INTERNALLY (validation,
calibration, metrics) and lives outside this repo (agent_workspace/). The
shipped, public artifact must be CCDC-FREE:

  * no CCDC/CSD data files committed (*.mol2, *.npz, ccdc_*fragment_index*, ...)
  * no shipped code importing the `ccdc` / `csd` runtime API
  * no shipped code importing CCDC-derived modules (grip*, mogul*)
  * no shipped code loading CCDC fragment libraries (grip_lib*, mogul_bounds, .npz)

ALLOWED open sources stay shippable: COD (cod_ideals, cod_*), Pyykkö/Cordero
covalent radii, VSEPR angles, ideal coordination polyhedra. These never match
the patterns below.

Escape hatch: a deliberate, documented reference may carry an inline
``# license-guard: allow`` comment on the same line. Use sparingly.

Usage
-----
  python scripts/license_guard.py --package       # scan all git-tracked .py
  python scripts/license_guard.py --staged        # scan staged changes (hook)
  python scripts/license_guard.py PATH [PATH...]   # scan specific files

Exit code 0 = clean, 1 = violation(s) found, 2 = usage error.
"""
from __future__ import annotations

import argparse
import ast
import os
import re
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Forbidden patterns (SINGLE SOURCE OF TRUTH)
# ---------------------------------------------------------------------------

# import roots that pull in the licensed CCDC/CSD runtime API
FORBIDDEN_IMPORT_ROOTS = {"ccdc", "csd"}

# any imported module whose dotted path has a component containing these
# substrings is a CCDC-derived module (grip = CCDC-fragment refinement,
# mogul = CCDC Mogul statistics).
FORBIDDEN_IMPORT_SUBSTR = ("grip", "mogul")

# CCDC data-file / library references inside a string literal.  Deliberately
# specific so the open COD assets (cod_ideals, cod_fragment_index) never match.
FORBIDDEN_CONTENT_RE = re.compile(
    r"(grip_lib|mogul_bounds|ccdc[_\-]?fragment|ccdc_ref|ccdc_validation|\.npz\b)",
    re.IGNORECASE,
)

# committed-file PATH patterns that must never enter the public repo.
FORBIDDEN_PATH_RE = re.compile(
    r"("
    r"\.mol2$|\.npz$|"
    r"grip_lib|mogul_bounds|"
    r"ccdc|CCDC_all|ccdc_ref_xyz|ccdc_exploit/|ccdc_validation\.py|"
    r"fragment_index_v|"
    r"fragment_index.*\.(json|sqlite|db)$|"
    r"paper_data/"
    r")",
    re.IGNORECASE,
)

ALLOW_COMMENT = "license-guard: allow"

# files that legitimately contain the forbidden tokens (as patterns/docs) and
# are therefore exempt from the content scan.
SELF_EXEMPT = {
    "scripts/license_guard.py",
    "tests/test_license_guard.py",
    "LICENSE_GUARD.md",
}


# ---------------------------------------------------------------------------
# Core scanning
# ---------------------------------------------------------------------------

def _line_allows(lines: list[str], lineno: int) -> bool:
    if 1 <= lineno <= len(lines):
        return ALLOW_COMMENT in lines[lineno - 1]
    return False


def scan_py_source(relpath: str, text: str) -> list[tuple[int, str]]:
    """Return [(lineno, message)] violations for one Python source string."""
    if relpath in SELF_EXEMPT:
        return []
    findings: list[tuple[int, str]] = []
    lines = text.splitlines()

    try:
        tree = ast.parse(text, filename=relpath)
    except SyntaxError:
        tree = None

    if tree is not None:
        for node in ast.walk(tree):
            mods: list[str] = []
            if isinstance(node, ast.Import):
                mods = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                mods = [node.module]
            lineno = getattr(node, "lineno", 0)
            for m in mods:
                comps = m.split(".")
                if comps[0].lower() in FORBIDDEN_IMPORT_ROOTS and not _line_allows(lines, lineno):
                    findings.append((lineno, f"imports CCDC/CSD runtime API: '{m}'"))
                elif any(sub in c.lower() for c in comps for sub in FORBIDDEN_IMPORT_SUBSTR) \
                        and not _line_allows(lines, lineno):
                    findings.append((lineno, f"imports CCDC-derived module: '{m}'"))

            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                hit = FORBIDDEN_CONTENT_RE.search(node.value)
                lineno = getattr(node, "lineno", 0)
                if hit and not _line_allows(lines, lineno):
                    findings.append((lineno, f"references CCDC data ('{hit.group(0)}') in a string literal"))
    else:
        # unparseable: regex fallback so nothing slips through silently
        for i, ln in enumerate(lines, 1):
            if ALLOW_COMMENT in ln:
                continue
            hit = FORBIDDEN_CONTENT_RE.search(ln)
            if hit:
                findings.append((i, f"references CCDC data ('{hit.group(0)}') (unparseable file, regex match)"))
    return findings


def check_path(relpath: str) -> str | None:
    """Return a message if the committed PATH itself is forbidden, else None."""
    if relpath in SELF_EXEMPT:
        return None
    if FORBIDDEN_PATH_RE.search(relpath):
        return f"forbidden CCDC/CSD file path: '{relpath}'"
    return None


# ---------------------------------------------------------------------------
# Git helpers
# ---------------------------------------------------------------------------

def _git(*args: str) -> str:
    return subprocess.run(["git", *args], capture_output=True, text=True).stdout


def tracked_py_files() -> list[str]:
    out = _git("ls-files", "*.py")
    return [p for p in out.splitlines() if p]


def staged_files() -> list[str]:
    out = _git("diff", "--cached", "--name-only", "--diff-filter=ACM")
    return [p for p in out.splitlines() if p]


def staged_blob(relpath: str) -> str:
    return _git("show", f":{relpath}")


# ---------------------------------------------------------------------------
# Modes
# ---------------------------------------------------------------------------

def _report(violations: list[str]) -> int:
    if violations:
        print("🚫 LICENSE GUARD: CCDC/CSD violation(s) — BLOCKED")
        for v in violations:
            print(f"   {v}")
        print()
        print("CCDC/CSD material is INTERNAL only (agent_workspace/), never shipped.")
        print("Ship open sources (COD, Pyykkö/Cordero radii, VSEPR) instead.")
        print("Deliberate exception: add '# license-guard: allow' on that line.")
        return 1
    print("✅ LICENSE GUARD: clean (no CCDC/CSD material).")
    return 0


def run_package() -> int:
    violations: list[str] = []
    for rel in tracked_py_files():
        if rel in SELF_EXEMPT:
            continue
        try:
            text = Path(rel).read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for lineno, msg in scan_py_source(rel, text):
            violations.append(f"{rel}:{lineno}: {msg}")
    return _report(violations)


def run_staged() -> int:
    violations: list[str] = []
    for rel in staged_files():
        path_msg = check_path(rel)
        if path_msg:
            violations.append(path_msg)
            continue  # don't also content-scan a file that's already path-forbidden
        if rel.endswith(".py") and rel not in SELF_EXEMPT:
            try:
                text = staged_blob(rel)
            except Exception:
                text = ""
            for lineno, msg in scan_py_source(rel, text):
                violations.append(f"{rel}:{lineno}: {msg}")
    return _report(violations)


def run_paths(paths: list[str]) -> int:
    violations: list[str] = []
    for p in paths:
        rel = os.path.relpath(p)
        path_msg = check_path(rel)
        if path_msg:
            violations.append(path_msg)
        if p.endswith(".py") and rel not in SELF_EXEMPT and Path(p).is_file():
            text = Path(p).read_text(encoding="utf-8", errors="replace")
            for lineno, msg in scan_py_source(rel, text):
                violations.append(f"{rel}:{lineno}: {msg}")
    return _report(violations)


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description="CCDC/CSD license guard (fail-closed).")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--package", action="store_true", help="scan all git-tracked .py files")
    g.add_argument("--staged", action="store_true", help="scan staged changes (pre-commit)")
    ap.add_argument("paths", nargs="*", help="explicit files to scan")
    ns = ap.parse_args(argv)

    if ns.staged:
        return run_staged()
    if ns.package or not ns.paths:
        return run_package()
    return run_paths(ns.paths)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
