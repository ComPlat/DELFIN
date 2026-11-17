"""Private fixture store for real ORCA output (LICENSE_GUARD.md, rule family 3).

Real ORCA output is EULA DATA and must never be committed to the public
repository.  It lives in an operator-managed store OUTSIDE the repo.
This module provides:

- a documented store layout: relative paths mirroring the repository
  (e.g. ``tests/fixtures/orca611/...``), plus a ``manifest.json`` with a
  sha256 per file, the ORCA version and a provenance note;
- the ``delfin-dev`` CLI (``init``/``add``/``verify``/``where``);
- ``store_root`` and ``private_fixtures_path`` used by the test helper
  in ``tests/test_lc1_fixtures_helper.py``.

The store root is resolved from the ``DELFIN_PRIVATE_FIXTURES``
environment variable first, then from the user settings file
(``paths.fixtures_dir``).  No default path is hard-coded.

``add`` and ``verify`` refuse anything that looks like personal data
(user names, host names, e-mail addresses, absolute home paths) or an
ORCA-written binary artifact (``.gbw``/``.hess``/``.engrad``); see
``scripts/license_guard.py`` for the public-repo counterpart.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from pathlib import Path

MANIFEST_SCHEMA = "delfin-private-fixtures/1"
MANIFEST_NAME = "manifest.json"
PRIVATE_FIXTURES_ENV = "DELFIN_PRIVATE_FIXTURES"
SETTINGS_KEY = ("paths", "fixtures_dir")

README_TEXT = """# Private DELFIN fixture store

Holds real QM-program output (ORCA EULA DATA) that must never be
committed to the public repository (LICENSE_GUARD.md, rule family 2).
Layout mirrors the repository: ``tests/fixtures/orca611/...`` etc., so
the test helper maps a relative repo path onto the same path here.

``manifest.json`` records, per file: sha256, the ORCA version and a
provenance note.  Personal data (user names, host names, e-mail
addresses, absolute home paths) and ORCA binary artifacts
(``.gbw``/``.hess``/``.engrad``) are refused.

CLI: ``delfin-dev fixtures init <dir> | add <src> --as <rel> |
verify | where``.
"""

# Personal data / machine identity that must not enter the shared store.
# IP addresses are also machine identity; keep the patterns text-only.
_PERSONAL_DATA_RE = re.compile(
    r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}"  # e-mail
    r"|\buser(?:name)?\s*[:=]\s*\S+"                   # user(name): value
    r"|\bhost(?:name)?\s*[:=]\s*\S+"                   # host(name): value
    r"|\bhostname\s+\S+"                               # hostname output
    r"|/home/\S+"                                      # absolute home path
    r"|/Users/\S+"                                     # macOS home path
    r"|\\\\[A-Za-z0-9_.-]+\\\\",                       # UNC path
    re.IGNORECASE,
)

# ORCA-written artifacts by extension (matches scripts/license_guard.py).
FORBIDDEN_ARTIFACT_RE = re.compile(r"\.(gbw|hess|engrad)$", re.IGNORECASE)


class FixturesStoreError(Exception):
    """Refused operation on the private fixture store."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def manifest_path(root: Path) -> Path:
    return Path(root) / MANIFEST_NAME


def read_manifest(root: Path) -> dict:
    path = manifest_path(root)
    if not path.exists():
        raise FixturesStoreError(
            f"no {MANIFEST_NAME} in {root}: run 'delfin-dev fixtures init' first"
        )
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise FixturesStoreError(
            f"unsupported manifest schema {manifest.get('schema')!r} "
            f"(expected {MANIFEST_SCHEMA})"
        )
    entries = manifest.get("entries")
    if not isinstance(entries, dict):
        raise FixturesStoreError("manifest has no 'entries' object")
    return manifest


def write_manifest(root: Path, manifest: dict) -> None:
    manifest_path(root).write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def scan_text_for_personal_data(text: str) -> str | None:
    """Return a description of the first personal-data hit, or None."""
    match = _PERSONAL_DATA_RE.search(text)
    if not match:
        return None
    return f"personal data near {match.group(0)!r}"


def check_file_for_store(file_path: Path) -> None:
    """Reject artifacts and content that must not enter the shared store."""
    if FORBIDDEN_ARTIFACT_RE.search(file_path.name):
        raise FixturesStoreError(
            f"refusing {file_path.name}: ORCA binary artifact "
            "(.gbw/.hess/.engrad) must not be distributed in the store"
        )
    try:
        text = file_path.read_text(encoding="utf-8", errors="strict")
    except UnicodeDecodeError:
        raise FixturesStoreError(
            f"refusing {file_path.name}: not UTF-8 text; the store holds "
            "text output only"
        )
    hit = scan_text_for_personal_data(text)
    if hit:
        raise FixturesStoreError(f"refusing {file_path.name}: {hit}")


def store_root(settings_path=None) -> Path:
    """Resolve the store root: env var first, then the settings file.

    Raises LookupError when neither is set -- never falls back to a
    hard-coded path.
    """
    from_env = os.environ.get(PRIVATE_FIXTURES_ENV, "").strip()
    if from_env:
        return Path(from_env).expanduser()
    from delfin.user_settings import load_settings
    settings = load_settings(settings_path)
    paths = settings.get(SETTINGS_KEY[0], {}) or {}
    configured = str(paths.get(SETTINGS_KEY[1], "") or "").strip()
    if configured:
        return Path(configured).expanduser()
    raise LookupError(
        f"{PRIVATE_FIXTURES_ENV} not set and paths.fixtures_dir not "
        "configured: no private fixture store location is known"
    )


def cmd_init(root: Path) -> int:
    root = Path(root).expanduser()
    if manifest_path(root).exists():
        print(f"store already initialised: {root}")
        return 0
    root.mkdir(parents=True, exist_ok=True)
    write_manifest(root, {"schema": MANIFEST_SCHEMA, "entries": {}})
    (root / "README.md").write_text(README_TEXT, encoding="utf-8")
    print(f"initialised private fixture store at {root}")
    return 0


def cmd_add(src: Path, rel: str, store: Path, orca_version: str,
            provenance: str) -> int:
    store = Path(store).expanduser()
    manifest = read_manifest(store)
    src = Path(src).expanduser()
    if not src.is_file():
        raise FixturesStoreError(f"source file not found: {src}")
    rel = rel.strip().strip("/")
    if not rel or ".." in rel.split("/") or rel.startswith("/"):
        raise FixturesStoreError(
            f"invalid --as {rel!r}: must be a repository-relative path "
            "like tests/fixtures/orca611/benzene"
        )
    dest_dir = store / rel
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / src.name
    check_file_for_store(src)
    dest.write_bytes(src.read_bytes())
    rel_path = f"{rel}/{src.name}"
    manifest["entries"][rel_path] = {
        "sha256": sha256_file(dest),
        "orca_version": orca_version,
        "provenance": provenance,
    }
    write_manifest(store, manifest)
    print(f"added {rel_path} (sha256 {manifest['entries'][rel_path]['sha256'][:12]}...)")
    return 0


def cmd_verify(store: Path) -> int:
    store = Path(store).expanduser()
    manifest = read_manifest(store)
    problems: list[str] = []
    for rel_path, entry in sorted(manifest["entries"].items()):
        # the recorded relative path itself must stay repository-relative
        if rel_path.startswith("/") or ".." in rel_path.split("/"):
            problems.append(f"{rel_path}: manifest path is not repository-relative")
            continue
        file_path = store / rel_path
        if not file_path.is_file():
            problems.append(f"{rel_path}: missing from store")
            continue
        if sha256_file(file_path) != entry.get("sha256"):
            problems.append(f"{rel_path}: sha256 does not match the manifest")
        hit = scan_text_for_personal_data(
            file_path.read_text(encoding="utf-8", errors="replace")
        )
        if hit:
            problems.append(f"{rel_path}: {hit}")
    if problems:
        for problem in problems:
            print(f"FAIL {problem}", file=sys.stderr)
        return 1
    print(f"verified {len(manifest['entries'])} entries in {store}")
    return 0


def cmd_where(store: Path | None) -> int:
    if store is not None:
        print(f"store (from --store): {Path(store).expanduser()}")
        return 0
    try:
        root = store_root()
    except LookupError as exc:
        print(f"no store resolved: {exc}")
        print(
            f"the test helper would skip (set {PRIVATE_FIXTURES_ENV} or "
            "paths.fixtures_dir in the DELFIN settings file)"
        )
        return 0
    origin = (
        f"environment variable {PRIVATE_FIXTURES_ENV}"
        if os.environ.get(PRIVATE_FIXTURES_ENV, "").strip()
        else "settings file paths.fixtures_dir"
    )
    print(f"store ({origin}): {root}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="delfin-dev fixtures",
        description=(
            "Manage the private fixture store for real ORCA output "
            "(LICENSE_GUARD.md rule family 3). The store lives outside "
            "the repository; its location comes from DELFIN_PRIVATE_FIXTURES "
            "or paths.fixtures_dir in the DELFIN settings file."
        ),
    )
    parser.add_argument("--store", default=None,
                        help="store root (default: env/settings)")
    sub = parser.add_subparsers(dest="command", required=True)
    p_init = sub.add_parser("init", help="create the store layout")
    p_init.add_argument("dir", help="store root to create")
    p_add = sub.add_parser("add", help="copy a file into the store and record it")
    p_add.add_argument("src", help="source file (outside the store)")
    p_add.add_argument("--as", dest="rel", required=True,
                       help="repository-relative destination directory")
    p_add.add_argument("--orca-version", default="",
                       help="ORCA version that produced the output")
    p_add.add_argument("--provenance", default="",
                       help="where the output came from (no host/user names)")
    p_add.add_argument("--store", dest="store_opt", default=None)
    sub.add_parser("verify", help="check hashes and personal-data rules")
    sub.add_parser("where", help="print the store root the helper will use")
    args = parser.parse_args(argv)

    store = getattr(args, "store_opt", None) or args.store
    try:
        if args.command == "init":
            return cmd_init(args.dir)
        if args.command == "add":
            root = Path(store) if store else store_root()
            return cmd_add(args.src, args.rel, root,
                           args.orca_version, args.provenance)
        if args.command == "verify":
            root = Path(store) if store else store_root()
            return cmd_verify(root)
        if args.command == "where":
            return cmd_where(store)
    except (FixturesStoreError, LookupError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    parser.error(f"unknown command {args.command!r}")
    return 2


def private_fixtures_path(rel: str, settings_path=None) -> Path:
    """One file inside the store, hash-checked against the manifest.

    The store root comes from DELFIN_PRIVATE_FIXTURES / the settings
    file (``store_root``).  ``rel`` is repository-relative (e.g.
    ``tests/fixtures/orca611/x/a.out``).  Raises LookupError when the
    store is not configured, the file is missing, or its sha256 does
    not match the manifest entry.
    """
    root = store_root(settings_path)
    manifest = read_manifest(root)
    rel = rel.strip().strip("/")
    file_path = Path(root) / rel
    if not file_path.is_file():
        raise LookupError(
            f"{PRIVATE_FIXTURES_ENV} set but {file_path} is missing"
        )
    entry = manifest["entries"].get(rel)
    if entry is None:
        raise LookupError(f"{rel} is not recorded in the store manifest")
    if sha256_file(file_path) != entry.get("sha256"):
        raise LookupError(
            f"{rel}: sha256 hash mismatch against the store manifest; the file "
            "was modified after 'delfin-dev fixtures add'"
        )
    return file_path


if __name__ == "__main__":
    raise SystemExit(main())
