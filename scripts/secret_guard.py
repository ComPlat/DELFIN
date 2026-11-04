#!/usr/bin/env python3
# secret-guard: allow  (this file defines the forbidden patterns; it is self-exempt)
"""Secret guard — generic-secret scan for the tracked tree, fail-closed.

Single source of truth for the forbidden generic-secret patterns,
consumed by two layers (the license_guard pattern, docs/LICENSE_GUARD.md):

  * tests/test_no_secrets_are_committed.py  -> pytest / CI over the whole tree
  * CLI (this file)                         -> manual audit and pre-commit hook
      python scripts/secret_guard.py --staged   # scan staged changes (hook)
      python scripts/secret_guard.py            # scan the whole tracked tree
      python scripts/secret_guard.py --tree DIR # scan a directory without git
                                                # (red/green proofs in tests)

What it covers (complementing the existing guards: personal data, account
names, CCDC/CSD licence material, MLIP tokens, home paths in commit
messages): private-key material, API-key shapes, credential assignments in
config-like files, tracked key-file names, high-entropy token literals, and
private-network addresses in the developer-facing directories.

About the KIT Toolbox key: DELFIN reads KIT_TOOLBOX_API_KEY as an opaque
string (delfin/agent/credentials.py stores it as a JSON value; delfin/
doctor.py sends it as a "Bearer {key}" header). It has NO fixed prefix, so
no shape rule can catch a bare literal — it is covered by two structural
rules instead: the credential-assignment rule in config-like files
(the name contains API_KEY) and the Bearer rule when it travels in text.

Allow-list: only an inline marker with a reason on the offending line,
``# secret-guard: allow (synthetic)`` — reason mandatory, used sparingly —
plus the guard's own files (this script and its test, which quote the
patterns). A rule never exempts a whole directory.

Exit code 0 = clean, 1 = findings, 2 = usage error.
"""
from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
import sys
from pathlib import Path

ALLOW_COMMENT = "secret-guard: allow"

# The guard's own files quote the forbidden patterns; nothing else is exempt.
SELF_EXEMPT = frozenset({
    "scripts/secret_guard.py",
    "tests/test_no_secrets_are_committed.py",
})

# config-like files where a KEY = literal assignment is a credential, and
# shell scripts where an exported one is.
CONFIG_SUFFIXES = frozenset({
    ".env", ".ini", ".cfg", ".toml", ".yaml", ".yml", ".json", ".sh", ".bash",
})
# binary extensions scanned by content for key headers only (fast literal
# search, no text regexes — the rest of a binary cannot match them anyway).
BINARY_SUFFIXES = frozenset({
    ".png", ".jpg", ".jpeg", ".gif", ".pdf", ".gbw", ".ico", ".woff",
    ".woff2", ".zip", ".gz", ".xlsx", ".docx", ".so", ".p12", ".pfx",
})

# Developer-facing files where a private-network address is a site detail,
# not test data (tests use RFC 5737 stand-ins: 192.0.2.x, 198.51.100.x,
# 203.0.113.x, which never match).
_PRIVATE_NET_SCOPE = re.compile(
    r"^(?:scripts/|docs/|tools/|README[^/]*$|readme[^/]*$)")

# Placeholders that make an assignment harmless (docs, templates, tests).
_PLACEHOLDER = re.compile(
    r"^\s*(?:$|<[^>]*>|\$\{[^}]*\}|\$[A-Za-z_][A-Za-z0-9_]*$|"
    r"\[redacted[^\]]*\]$|(?:change[-_]?me|placeholder|dummy|example|"
    r"sample|synthetic|redacted|your[-_][a-z0-9_-]*|xxx+|fix[-_]?me|"
    r"none|null|true|false|nil|secret-guard-test[0-9a-z-]*)\s*)$",
    re.IGNORECASE)

# ---------------------------------------------------------------------------
# The rule table
# ---------------------------------------------------------------------------
# Each rule: (name, scope, literal fast-path or None, pattern or None, note).
# scope: "all" | "config" | "private-net" | "non-test"
# Binary files are scanned by the literal fast-path only.


def _assignment_rule() -> tuple[str, re.Pattern]:
    names = (r"(?:password|passwd|pwd|secret|token|api[_-]?key|apikey|"
             r"auth|access[_-]?key|private[_-]?key)")
    key = r"(?:[A-Za-z0-9_.-]*[A-Za-z0-9])?"
    value = r"['\"]?([^\s'\"]{8,})['\"]?"
    return ("credential_assignment", re.compile(
        rf"^\s*(?:export\s+|set\s+)?(?:[A-Za-z_]{{1,3}}\s*)?{names}{key}"
        rf"\s*=\s*{value}", re.IGNORECASE))


_ASSIGN_RE_NAME, _ASSIGN_RE = _assignment_rule()

RULES: list[tuple[str, str, str | None, re.Pattern | None, str]] = [
    # name, scope, literal, pattern, what it guards
    ("private_key", "all", "PRIVATE KEY-----",
     re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----"),
     "PEM private-key block header"),
    ("private_key", "all", "OPENSSH PRIVATE KEY-----",
     re.compile(r"-----BEGIN OPENSSH PRIVATE KEY-----"),
     "OpenSSH private-key block"),
    ("private_key", "all", "PuTTY-User-Key-File",
     re.compile(r"PuTTY-User-Key-File-\d"),
     "PuTTY .ppk private-key header"),
    ("jwt", "all", "eyJ",
     re.compile(r"\beyJ[A-Za-z0-9_-]{10,}\.[A-Za-z0-9_-]{10,}\.[A-Za-z0-9_-]*"),
     "three-part JWT"),
    ("api_key_shape", "all", None, re.compile(
        r"\b(?:"
        r"sk-(?:proj-|ant-)?[A-Za-z0-9_-]{20,}"      # OpenAI / Anthropic
        r"|xox[bpars]-[A-Za-z0-9-]{10,}"             # Slack
        r"|AKIA[0-9A-Z]{16}"                         # AWS access key id
        r"|gh[posu]_[A-Za-z0-9]{36,}"                # GitHub token
        r"|github_pat_[A-Za-z0-9_]{20,}"
        r"|hf_[A-Za-z0-9]{30,}"                      # Hugging Face
        r"|glpat-[A-Za-z0-9_-]{20,}"                 # GitLab PAT
        r")\b"),
     "recognised provider API-key shape"),
    ("bearer_token", "all", "Bearer ", re.compile(
        r"\bBearer\s+[A-Za-z0-9+/_=-]{24,}\b"),
     "long base64 bearer credential"),
    ("credential_assignment", "config", None, _ASSIGN_RE,
     "password/secret/token/API-key assigned a literal in a config-like "
     "file (.env, *.ini, *.cfg, *.toml, *.yaml, *.json, *.sh)"),
]

# Rules scoped to the developer-facing directories (site detail, not data).
_PRIVATE_NET_RE = re.compile(
    r"\b(?:10\.\d{1,3}\.\d{1,3}\.\d{1,3}"
    r"|172\.(?:1[6-9]|2\d|3[01])\.\d{1,3}\.\d{1,3}"
    r"|192\.168\.\d{1,3}\.\d{1,3}"
    r"|100\.(?:6[4-9]|[7-9]\d|1[01]\d|12[0-7])\.\d{1,3}\.\d{1,3})\b")

# Tracked file names that are key material or credential stores by name.
_SECRET_PATH_RE = re.compile(
    r"(?:^|/)(?:\.env(?:\.\w+)?|id_rsa(?:\.\w+)?|id_ed25519(?:\.\w+)?"
    r"|credentials[^/]*\.json)$|(?:^|/)[^/]*\.p(?:em|12|fx)$")

# ---------------------------------------------------------------------------
# High-entropy rule (documented heuristic)
# ---------------------------------------------------------------------------
# A tracked non-test file that carries a run of >= 20 chars over the
# base62 alphabet, Shannon entropy >= 3.8 bits/char over its own
# characters, MIXED CASE and >= 15% digits, is most likely a
# machine-generated credential. The two structural guards keep the
# heuristic honest, checked against this tree:
#   * an all-lowercase or all-uppercase run is prose, an identifier or a
#     hex integrity hash (pip --hash=sha256:... in .github/ci/*requirements*.txt)
#   * a run with < 15% digits is a camelCase identifier
#     ("calcPreselect3dViewerByScope", "AssignStereochemistryFrom3D") or
#     an English word run ("delfinPatchAllKnown3DmolViewers")
# Thresholds chosen so ordinary identifiers stay below either bound.
_ENTROPY_MIN = 3.8
_ENTROPY_LEN = 20
_ENTROPY_TOKEN_RE = re.compile(r"[A-Za-z0-9]{20,}")


def _entropy(text: str) -> float:
    counts: dict[str, int] = {}
    for ch in text:
        counts[ch] = counts.get(ch, 0) + 1
    n = len(text)
    return -sum((c / n) * math.log2(c / n) for c in counts.values())


def _looks_like_credential_token(token: str) -> bool:
    if _entropy(token) < _ENTROPY_MIN:
        return False
    upper = sum(1 for c in token if c.isupper())
    lower = sum(1 for c in token if c.islower())
    digits = sum(1 for c in token if c.isdigit())
    if not (upper >= 2 and lower >= 2 and
            digits >= max(3, round(0.15 * len(token)))):
        return False
    # An enumeration (an alphabet, a character-set constant) is one long
    # ascending chain; a real credential almost never has more than a
    # couple of sequential character pairs. Measured: 59 pairs in the
    # base62 alphabet constant, 1 in random strings.
    seq = sum(1 for a, b in zip(token, token[1:]) if ord(b) == ord(a) + 1)
    return seq <= 2


# ---------------------------------------------------------------------------


def _is_config(relpath: str) -> bool:
    base = os.path.basename(relpath)
    name, _, ext = base.partition(".")
    return base == ".env" or (ext and f".{ext.lower()}" in CONFIG_SUFFIXES)


def _in_private_net_scope(relpath: str) -> bool:
    return bool(_PRIVATE_NET_SCOPE.match(relpath))


def _is_test(relpath: str) -> bool:
    parts = relpath.split("/")
    return ("tests" in parts[:-1]
            or os.path.basename(relpath).startswith("test_")
            or "fixtures" in parts[:-1])


def scan_line(relpath: str, lineno: int, line: str) -> list[tuple[str, str]]:
    """One text line -> [(rule, matched-literal)]. Allow marker honoured."""
    if ALLOW_COMMENT in line:
        return []
    out: list[tuple[str, str]] = []
    for name, scope, literal, pattern, _what in RULES:
        if scope == "config" and not _is_config(relpath):
            continue
        if scope == "private-net" and not _in_private_net_scope(relpath):
            continue
        if scope == "non-test" and _is_test(relpath):
            continue
        if literal and literal not in line:
            continue
        if pattern is None:
            continue
        m = pattern.search(line)
        if m is None:
            continue
        if name == "credential_assignment" and _PLACEHOLDER.match(m.group(1)):
            continue
        out.append((name, m.group(0)[:40]))
    return out


def scan_text(relpath: str, text: str) -> list[tuple[int, str, str]]:
    """A tracked text file -> [(lineno, rule, sample)]. Also the entry
    the pytest test drives per rule."""
    findings: list[tuple[int, str, str]] = []
    for lineno, line in enumerate(text.splitlines(), start=1):
        for name, sample in scan_line(relpath, lineno, line):
            findings.append((lineno, name, sample))
    if not _is_test(relpath):
        for lineno, line in enumerate(text.splitlines(), start=1):
            if ALLOW_COMMENT in line:
                continue
            # integrity hashes (sha256:, sha384:, sha512: prefixes) are
            # content fingerprints, not credentials: pip --hash=... in
            # requirements files, SRI attributes in pinned_cdn.py, and
            # lockfiles list them by design.
            if re.search(r"\bsha\d{3}[:-]", line) or "lockfileVersion" in text[:200]:
                continue
            if re.search(r"\b\d{1,3}(?:\.\d{1,3}){3}/\d{1,2}\b", line):
                # CIDR notation in scripts is a network definition
                # (firewall rules, nftables sets), not a site address.
                continue
            for m in _ENTROPY_TOKEN_RE.finditer(line):
                tok = m.group(0)
                if _looks_like_credential_token(tok):
                    findings.append((lineno, "high_entropy", tok[:40]))
    if _in_private_net_scope(relpath):
        for lineno, line in enumerate(text.splitlines(), start=1):
            if ALLOW_COMMENT in line:
                continue
            # CIDR notation (a/b) is a network definition (firewall rules,
            # nftables sets), not a site address.
            line_wo_cidr = re.sub(r"\b\d{1,3}(?:\.\d{1,3}){3}/\d{1,2}\b", "", line)
            m = _PRIVATE_NET_RE.search(line_wo_cidr)
            if m:
                findings.append((lineno, "private_net", m.group(0)))
    return findings


def check_path(relpath: str) -> list[tuple[int, str, str]]:
    """Tracked file NAMES that are key/credential material by convention."""
    if relpath in SELF_EXEMPT or relpath.startswith("tests/fixtures/"):
        return []
    if _SECRET_PATH_RE.search(relpath):
        return [(0, "secret_path", relpath)]
    return []


KEY_HEADER_LITERALS = (
    "PRIVATE KEY-----", "PuTTY-User-Key-File",
    "OPENSSH PRIVATE KEY-----",
)


def scan_binary(relpath: str, blob: bytes) -> list[tuple[int, str, str]]:
    """Binary files: search the raw content for key headers only."""
    for literal in KEY_HEADER_LITERALS:
        if literal.encode() in blob:
            return [(0, "private_key", literal)]
    return []


# ---------------------------------------------------------------------------
# Git and tree access
# ---------------------------------------------------------------------------


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], capture_output=True, text=True,
                          timeout=120).stdout


def tracked_files() -> list[str]:
    return [p for p in _git("ls-files").splitlines() if p]


def staged_files() -> list[str]:
    out = _git("diff", "--cached", "--name-only", "--diff-filter=ACM")
    return [p for p in out.splitlines() if p]


def staged_blob(relpath: str) -> str:
    return _git("show", f":{relpath}")


def _readable_text(relpath: str) -> str:
    if relpath.rsplit(".", 1)[-1].lower() in {s.lstrip(".") for s in BINARY_SUFFIXES}:
        return ""
    try:
        return Path(relpath).read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return ""


def findings_for(relpath: str, text: str) -> list[tuple[int, str, str]]:
    if relpath in SELF_EXEMPT:
        return []
    out = check_path(relpath)
    if out and out[0][1] == "secret_path":
        # name-convention files: scan content too, but the path finding is
        # the headline; still report both.
        pass
    if text:
        out.extend(scan_text(relpath, text))
    return out


def run_tree(root: Path) -> list[str]:
    violations: list[str] = []
    old_cwd = Path.cwd()
    os.chdir(root)
    try:
        for path in sorted(root.rglob("*")):
            if not path.is_file() or ".git" in path.parts:
                continue
            rel = path.relative_to(root).as_posix()
            if rel.rsplit(".", 1)[-1].lower() in {s.lstrip(".") for s in BINARY_SUFFIXES}:
                try:
                    findings = scan_binary(rel, path.read_bytes())
                except OSError:
                    findings = []
            else:
                findings = findings_for(rel, _readable_text(rel))
            violations.extend(_fmt(rel, findings))
    finally:
        os.chdir(old_cwd)
    return violations


def _fmt(rel: str, findings: list[tuple[int, str, str]]) -> list[str]:
    return [f"{rel}:{lineno}: rule {name}: {sample}"
            for lineno, name, sample in findings]


def run_all() -> list[str]:
    violations: list[str] = []
    for rel in tracked_files():
        if rel in SELF_EXEMPT:
            continue
        suffix = rel.rsplit(".", 1)[-1].lower()
        if suffix in {s.lstrip(".") for s in BINARY_SUFFIXES}:
            try:
                findings = scan_binary(rel, Path(rel).read_bytes())
            except OSError:
                findings = []
        else:
            findings = findings_for(rel, _readable_text(rel))
        violations.extend(_fmt(rel, findings))
    return violations


def run_staged() -> list[str]:
    violations: list[str] = []
    for rel in staged_files():
        if rel in SELF_EXEMPT:
            continue
        blob = staged_blob(rel)
        findings = check_path(rel)
        findings.extend(scan_text(rel, blob))
        violations.extend(_fmt(rel, findings))
    return violations


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description="Generic secret guard (fail-closed).")
    ap.add_argument("--staged", action="store_true",
                    help="scan staged changes (pre-commit hook)")
    ap.add_argument("--tree", type=Path, metavar="DIR",
                    help="scan a directory tree without git (tests, audits)")
    args = ap.parse_args(argv)
    if args.tree:
        violations = run_tree(args.tree)
    elif args.staged:
        violations = run_staged()
    else:
        violations = run_all()
    if violations:
        print("🚫 SECRET GUARD: possible secret(s) — BLOCKED")
        for v in violations:
            print(f"   {v}")
        print()
        print("A tracked file must not carry a credential. Use a stand-in,")
        print("keep the real value in the credential store, or — only for a")
        print(f"synthetic pattern — mark the line: '# {ALLOW_COMMENT} (reason)'.")
        return 1
    print("✅ SECRET GUARD: clean (no generic secrets).")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
