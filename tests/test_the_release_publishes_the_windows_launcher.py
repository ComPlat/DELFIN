"""A release carries the Windows launcher, so it has an address.

The launcher ZIP is committed in the tree and was published nowhere: not
a release asset, not on the site. The only way a new person could get it
was for someone who already had it to send it to them -- once per person,
with whatever version that person happened to have.

A release asset gives it one address and ties a ZIP to a version, which
also makes "which launcher are you on" a question with an answer.

The file is fetched through the API rather than by checking the tree out:
the release job deliberately has no checkout (see
``test_the_release_notes_need_no_clone_and_no_annotated_tag`` -- that is
what keeps ``--notes-from-tag`` out of it), and one file does not need a
clone.
"""

from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
ZIP = "delfin-windows-launcher.zip"


def _release_job() -> dict:
    data = yaml.safe_load(
        (ROOT / ".github" / "workflows" / "release.yml").read_text(
            encoding="utf-8"))
    return data["jobs"]["github"]


def _runs() -> str:
    """Every line the job actually runs, comments dropped: a guard has to
    read what runs, not what a comment says about it."""
    return "\n".join(
        line for step in _release_job()["steps"]
        for line in step.get("run", "").splitlines()
        if not line.lstrip().startswith("#"))


def test_the_launcher_is_uploaded_with_the_release():
    runs = _runs()

    assert runs.count(ZIP) >= 3, (
        "the launcher has to be fetched, listed in the checksums and "
        "uploaded; it is published nowhere else")
    for command in ("gh release upload", "gh release create"):
        start = runs.index(command)
        end = runs.find("\n          else", start)
        chunk = runs[start:end if end > 0 else len(runs)]
        assert ZIP in chunk, f"{command} does not carry the launcher"


def test_its_checksum_is_published_with_it():
    runs = _runs()

    assert f"sha256sum -- {ZIP}" in runs, (
        "a download nobody can check is a download nobody should run")
    assert "SHA256SUMS.txt" in runs


def test_the_checksums_do_not_imply_an_attestation_it_does_not_have():
    """The distributions are attested and the launcher is not. Keeping it
    out of dist/ keeps one list from implying the other's provenance."""
    runs = _runs()

    assert f"mv repo/tools/{ZIP} dist/" not in runs
    assert f"> {ZIP}" in runs, "it is fetched into the job's own directory"


def test_it_is_fetched_without_a_checkout():
    job = _release_job()

    assert not [s for s in job["steps"]
                if "checkout" in str(s.get("uses", ""))], (
        "this job has no checkout; fetching one file must not add one")
    assert "gh api" in _runs()


def test_the_file_the_release_names_is_the_one_in_the_tree():
    """The path in the workflow has to be the path the package build
    writes, or the release gets a 404 at tag time and nothing else says
    so until then."""
    assert (ROOT / "tools" / ZIP).is_file()
    assert f"contents/tools/{ZIP}" in _runs()
