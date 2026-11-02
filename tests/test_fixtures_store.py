"""Private fixture store (LICENSE_GUARD.md rule family 2/3).

The store lives OUTSIDE the repository, operator-managed, and holds real
ORCA output that must never be committed.  Everything here uses
synthetic files only; the store root is taken from DELFIN_PRIVATE_FIXTURES
(or --store), never hard-coded.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

import delfin.fixtures_store as fs

SYNTH = "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"


def _write(path: Path, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(SYNTH + body, encoding="utf-8")


def _run(*args: str, store: Path | None = None, env_extra=None):
    env = dict(os.environ)
    if store is not None:
        env["DELFIN_PRIVATE_FIXTURES"] = str(store)
    else:
        env.pop("DELFIN_PRIVATE_FIXTURES", None)
    if env_extra:
        env.update(env_extra)
    return subprocess.run(
        [sys.executable, "-m", "delfin.fixtures_store", *args],
        capture_output=True, text=True, env=env, timeout=120,
    )


@pytest.fixture()
def store(tmp_path):
    root = tmp_path / "private-store"
    rc = _run("init", str(root))
    assert rc.returncode == 0, rc.stderr
    return root


def test_init_creates_manifest_with_layout(store):
    manifest = json.loads((store / "manifest.json").read_text())
    assert manifest["schema"] == fs.MANIFEST_SCHEMA
    assert manifest["entries"] == {}
    assert (store / "README.md").exists()


def test_add_copies_hashes_and_records(store, tmp_path):
    src = tmp_path / "benzene.out"
    _write(src, "some synthetic body\n")
    rc = _run("add", str(src), "--as", "tests/fixtures/orca611/benzene",
              "--orca-version", "6.1.1",
              "--provenance", "synthetic test run",
              store=store)
    assert rc.returncode == 0, rc.stderr
    copied = store / "tests/fixtures/orca611/benzene/benzene.out"
    assert copied.is_file()
    manifest = json.loads((store / "manifest.json").read_text())
    entry = manifest["entries"]["tests/fixtures/orca611/benzene/benzene.out"]
    assert entry["sha256"] == fs.sha256_file(copied)
    assert entry["orca_version"] == "6.1.1"
    assert "provenance" in entry


def test_add_rejects_orca_written_names_and_personal_data(store, tmp_path):
    src = tmp_path / "job.gbw"
    src.write_bytes(b"binary-ish")
    rc = _run("add", str(src), "--as", "tests/fixtures/x/job", store=store)
    assert rc.returncode != 0

    src2 = tmp_path / "bad.out"
    src2.write_text(SYNTH + "contact: alice@example.com\n", encoding="utf-8")
    rc = _run("add", str(src2), "--as", "tests/fixtures/x/bad", store=store)
    assert rc.returncode != 0


def test_verify_ok_and_detects_tamper_and_personal_data(store, tmp_path):
    src = tmp_path / "good.out"
    _write(src, "clean body\n")
    _run("add", str(src), "--as", "tests/fixtures/orca611/good", store=store)
    rc = _run("verify", store=store)
    assert rc.returncode == 0, rc.stderr + rc.stdout

    # tamper: rewrite a stored file
    target = store / "tests/fixtures/orca611/good/good.out"
    target.write_text(SYNTH + "tampered\n", encoding="utf-8")
    rc = _run("verify", store=store)
    assert rc.returncode != 0

    # personal data smuggled in after add
    target.write_text(SYNTH + "host: somehost\n", encoding="utf-8")
    rc = _run("verify", store=store)
    assert rc.returncode != 0


def test_where_prints_store_root(store, tmp_path):
    rc = _run("where", store=store)
    assert rc.returncode == 0
    assert str(store) in rc.stdout
    # without the env var it reports what the helper would use (unset)
    rc = _run("where")
    assert rc.returncode == 0
    assert "DELFIN_PRIVATE_FIXTURES" in rc.stdout


def test_helper_refuses_hash_mismatch(store, tmp_path, monkeypatch):
    """The helper reads the manifest and refuses a file whose hash differs."""
    src = tmp_path / "h.out"
    _write(src, "hash me\n")
    _run("add", str(src), "--as", "tests/fixtures/orca611/h", store=store)
    monkeypatch.setenv("DELFIN_PRIVATE_FIXTURES", str(store))
    assert fs.private_fixtures_path("tests/fixtures/orca611/h/h.out").exists()

    target = store / "tests/fixtures/orca611/h/h.out"
    target.write_text(SYNTH + "changed\n", encoding="utf-8")
    with pytest.raises(LookupError, match="hash"):
        fs.private_fixtures_path("tests/fixtures/orca611/h/h.out")


def test_store_root_env_precedence(tmp_path, monkeypatch):
    """env var wins over the settings-file default; neither is hard-coded."""
    settings = tmp_path / "settings.json"
    settings.write_text(json.dumps({"paths": {"fixtures_dir": "/from/settings"}}))
    monkeypatch.setenv("DELFIN_PRIVATE_FIXTURES", "/from/env")
    assert fs.store_root(settings_path=settings) == Path("/from/env")
    monkeypatch.delenv("DELFIN_PRIVATE_FIXTURES")
    assert fs.store_root(settings_path=settings) == Path("/from/settings")
