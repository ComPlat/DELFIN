"""A release carries a CycloneDX SBOM, attested and attached to the release.

Input: .github/workflows/release.yml and .github/ci/sbom-requirements.txt.

Output: failures for
- no ``sbom`` job, or one that holds more than ``contents: read``, runs
  before the distributions exist, or installs its generator other than by
  hash (``--require-hashes -r .github/ci/sbom-requirements.txt``);
- an SBOM that is not CycloneDX JSON from ``cyclonedx-py environment``
  over an environment holding the built wheel, or not reproducible
  (``--output-reproducible``, SOURCE_DATE_EPOCH = commit time);
- the generator or the wheel's dependencies installed in the ``build``
  job: that job's output is what is published, and an install there runs
  code from PyPI beside dist/;
- an SBOM that is not attested against the distributions
  (actions/attest with ``sbom-path`` in ``provenance``), or not uploaded
  to the GitHub release and listed in SHA256SUMS.txt;
- an SBOM in the files handed to PyPI (``pypi`` downloads dist only).

CIS Software Supply Chain 2.4.5/2.4.6, SSDF PS.3.2. actions/attest is used
rather than actions/attest-sbom: the latter (v4.1.0) is a deprecated
wrapper that prints a warning and calls actions/attest.
"""

from __future__ import annotations

import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / ".github" / "workflows" / "release.yml"
SBOM_REQS = ROOT / ".github" / "ci" / "sbom-requirements.txt"


def _jobs() -> dict:
    return yaml.safe_load(RELEASE.read_text(encoding="utf-8"))["jobs"]


def _runs(job: dict) -> str:
    return "\n".join(s.get("run", "") for s in job["steps"])


def _uses(job: dict) -> list[dict]:
    return [s for s in job["steps"] if "uses" in s]


def test_the_sbom_job_reads_only_and_installs_its_tool_by_hash():
    jobs = _jobs()
    assert "sbom" in jobs, "the release has no SBOM job"
    job = jobs["sbom"]
    assert job["permissions"] == {"contents": "read"}
    assert "build" in (job["needs"] if isinstance(job["needs"], list) else [job["needs"]])
    runs = _runs(job)
    assert "pip install --require-hashes -r .github/ci/sbom-requirements.txt" in runs
    assert "cyclonedx-py environment" in runs
    assert "--output-reproducible" in runs and "SOURCE_DATE_EPOCH" in runs
    assert "--of JSON" in runs
    assert re.search(r"pip --python \S+ install [^\n]*dist/\*\.whl", runs), (
        "the SBOM must describe an environment holding the built wheel")
    uploads = [s for s in _uses(job) if s["uses"].startswith("actions/upload-artifact@")]
    assert uploads and uploads[0]["with"]["name"] == "sbom"


def test_the_generator_is_pinned_with_hashes():
    text = SBOM_REQS.read_text(encoding="utf-8")
    assert re.search(r"^cyclonedx-bom==\d[\w.]* \\\n\s+--hash=sha256:[0-9a-f]{64}", text, re.M)


def test_the_build_job_installs_nothing_for_the_sbom():
    runs = _runs(_jobs()["build"])
    assert "cyclonedx" not in runs and "dist/*.whl" not in runs


def test_the_sbom_is_attested_against_the_distributions():
    job = _jobs()["provenance"]
    assert "sbom" in job["needs"]
    attest = [s for s in _uses(job) if s["uses"].startswith("actions/attest@")]
    assert attest, "no SBOM attestation"
    params = attest[0]["with"]
    assert params["subject-path"] == "dist/*"
    assert "sbom" in str(params["sbom-path"])


def test_the_sbom_goes_to_the_github_release_and_not_to_pypi():
    jobs = _jobs()
    names = [s["with"]["name"] for s in _uses(jobs["github"])
             if s["uses"].startswith("actions/download-artifact@")]
    assert "sbom" in names
    runs = _runs(jobs["github"])
    assert "sbom/" in runs and "SHA256SUMS.txt" in runs
    pypi = [s["with"]["name"] for s in _uses(jobs["pypi"])
            if s["uses"].startswith("actions/download-artifact@")]
    assert pypi == ["dist"]
