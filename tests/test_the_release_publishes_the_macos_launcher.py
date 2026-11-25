"""A release carries the macOS launcher, as it carries the Windows one.

Same reasons as tests/test_the_release_publishes_the_windows_launcher.py:
one address per version, a checksum beside it, fetched through the API
because the release job has no checkout, and kept out of dist/ so the
checksum list does not imply an attestation the ZIP does not have.
"""

from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
ZIP = "delfin-macos-launcher.zip"


def _runs() -> str:
    data = yaml.safe_load((ROOT / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8"))
    return "\n".join(
        line for step in data["jobs"]["github"]["steps"]
        for line in step.get("run", "").splitlines()
        if not line.lstrip().startswith("#"))


def test_the_mac_launcher_is_uploaded_with_the_release():
    runs = _runs()
    for command in ("gh release upload", "gh release create"):
        start = runs.index(command)
        end = runs.find("\n          else", start)
        chunk = runs[start:end if end > 0 else len(runs)]
        assert ZIP in chunk, f"{command} does not carry the macOS launcher"


def test_its_checksum_is_published_and_it_stays_out_of_dist():
    runs = _runs()
    assert f"sha256sum -- {ZIP}" in runs
    assert f"> {ZIP}" in runs
    assert f"mv repo/tools/{ZIP} dist/" not in runs


def test_the_file_the_release_names_is_the_one_in_the_tree():
    assert (ROOT / "tools" / ZIP).is_file()
    assert f"contents/tools/{ZIP}" in _runs()
