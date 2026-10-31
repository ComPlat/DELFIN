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

# ---------------------------------------------------------------------------
# Licensed QM program output (LC1): signatures of DATA produced by licensed
# QM codes that may not ship in this public LGPL repo.
#   ORCA EULA 2025 section 3: DATA may be shared only for ACADEMIC PURPOSES
#   (c), to databases only non-commercially with an EULA notice (d), software
#   using DATA must be under the EULA (e), ML on DATA only academic (g).
#   See LICENSE_GUARD.md. Deliberately line-oriented: one signature per line
#   so the scan stays cheap over a whole tracked tree.
# ---------------------------------------------------------------------------

# The one documented exception: a synthetic fixture declares itself on its
# FIRST line with this marker.
QM_SYNTHETIC_MARKER = "# SYNTHETIC test fixture written for DELFIN -- not ORCA output"

QM_OUTPUT_SIGNATURES = (
    r"^\s*\* O\s+R\s+C\s+A\s+\*$",                    # ORCA banner
    r"ORCA\s+TERMINATED\s+NORMALLY",                   # normal-termination block
    r"^\s*ORCA\s+terminated\s+normally",               # lowercase variant
    r"Program\s+Version\s+6\.\d+\.\d+\s+-\s+RELEASE",  # ORCA version line
    r"^\s*\* Host name:",                              # ORCA run header
    r"^\s*\* Working dir\.:",                          # ORCA run header
    r"^\|\s*\d+>",                                     # ORCA INPUT FILE echo block
    r"THERMOCHEMISTRY\s+AT\s+\d+\.\d+K",               # ORCA thermochem block
)

QM_OUTPUT_RE = re.compile(
    "(" + "|".join(QM_OUTPUT_SIGNATURES) + ")",
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
    r"paper_data/|"
    r"\.(gbw|hess|engrad)$"          # ORCA binary/vector artifacts (EULA DATA)
    r")",
    re.IGNORECASE,
)


# ORCA-generated GOAT ensembles must never be added to DELFIN. Matched by the
# file names GOAT itself writes -- a broader "goat|orca in the path" rule
# would also hit DELFIN's own knowledge cards (delfin/orca_knowledge/**.json).
# Other ORCA output is caught by its content (output signature), not its name.
ORCA_GENERATED_PATH_RE = re.compile(
    r"\.(?:finalensemble|globalminimum)(?:\.\d+)?\.xyz$",
    re.IGNORECASE,
)

ALLOW_COMMENT = "license-guard: allow"

# files that legitimately contain the forbidden tokens (as patterns/docs) and
# are therefore exempt from the content scan.
SELF_EXEMPT = {
    "scripts/license_guard.py",
    "tests/test_license_guard.py",
    "LICENSE_GUARD.md",
    # the LC1 licence inventory: a policy document that quotes the
    # signatures it classifies — data, not QM output.
    "tests/fixtures/lc1_synthetic/INVENTORY.md",
    # documentation that QUOTES the signature strings (how DELFIN reads
    # termination, what to search for) — no program output is embedded,
    # the phrase is the name of a state, not the output itself.
    "docs/PPKA_VALIDATION.md",
    "docs/RETRY_LOGIC.md",
    "docs/methodology.md",
    # installation walkthrough quoting the version banner of the local
    # ORCA install as a verification step — a version number, not DATA.
    "examples/example_Job_Submission_Scripts/BwUniCluster/DELFIN_INSTALLATION_ANLEITUNG.txt",
    # dashboard-agent corpus: transcript TEXT that mentions the phrase in
    # prose (which file has which status); no ORCA output block is in it.
    "tests/fixtures/choice_offer_corpus.json",
}


# ---------------------------------------------------------------------------
# Core scanning
# ---------------------------------------------------------------------------

def _line_allows(lines: list[str], lineno: int) -> bool:
    if 1 <= lineno <= len(lines):
        return ALLOW_COMMENT in lines[lineno - 1]
    return False


def scan_path_text(relpath: str, text: str) -> list[tuple[int, str]]:
    """Scan ONE file of any type for licensed QM-program output.

    Returns [(lineno, message)]. Binary bytes are skipped silently (no
    printable signature can be relied on there — the path rule covers
    the known binary artifacts, .gbw/.hess/.engrad).

    Exception: a file whose FIRST line carries QM_SYNTHETIC_MARKER
    declares itself a hand-written synthetic fixture and is exempt as a
    whole (such a fixture legitimately contains ORCA-format lines).
    """
    findings: list[tuple[int, str]] = []
    lines = text.splitlines()
    if lines and lines[0].strip() == QM_SYNTHETIC_MARKER:
        return findings
    for i, ln in enumerate(lines, 1):
        hit = QM_OUTPUT_RE.search(ln)
        if hit:
            findings.append((i, f"looks like licensed QM program output ('{hit.group(0)}')"))
    return findings


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
    if ORCA_GENERATED_PATH_RE.search(relpath):
        return f"forbidden ORCA-generated GOAT path: '{relpath}'"
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


def tracked_text_files() -> list[str]:
    """All git-tracked files EXCEPT the self-exempt ones.

    Non-.py files are scanned too: ORCA/licensed-program output reaches
    the repo as .out, .hess, .gbw, .dat, .property.json — extensions the
    old .py-only scan never saw (LC1).
    """
    out = _git("ls-files")
    return [p for p in out.splitlines() if p and p not in SELF_EXEMPT]


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
    for rel in tracked_text_files():
        if rel in SELF_EXEMPT:
            continue
        try:
            text = Path(rel).read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        if rel.endswith(".py"):
            findings = scan_py_source(rel, text)
        else:
            findings = scan_path_text(rel, text)
        for lineno, msg in findings:
            violations.append(f"{rel}:{lineno}: {msg}")
    return _report(violations)


def run_staged() -> int:
    violations: list[str] = []
    for rel in staged_files():
        path_msg = check_path(rel)
        if path_msg:
            violations.append(path_msg)
            continue  # don't also content-scan a file that's already path-forbidden
        if rel not in SELF_EXEMPT:
            try:
                text = staged_blob(rel)
            except Exception:
                text = ""
            findings = scan_py_source(rel, text) if rel.endswith(".py") else scan_path_text(rel, text)
            for lineno, msg in findings:
                violations.append(f"{rel}:{lineno}: {msg}")
    return _report(violations)


def run_paths(paths: list[str]) -> int:
    violations: list[str] = []
    for p in paths:
        rel = os.path.relpath(p)
        path_msg = check_path(rel)
        if path_msg:
            violations.append(path_msg)
        if p.endswith(".py"):
            if rel not in SELF_EXEMPT and Path(p).is_file():
                text = Path(p).read_text(encoding="utf-8", errors="replace")
                findings = scan_py_source(rel, text)
            else:
                findings = []
        elif rel not in SELF_EXEMPT and Path(p).is_file():
            text = Path(p).read_text(encoding="utf-8", errors="replace")
            findings = scan_path_text(rel, text)
        else:
            findings = []
        for lineno, msg in findings:
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
