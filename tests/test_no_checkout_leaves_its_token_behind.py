"""Every actions/checkout step removes its token from .git/config.

Input: every step of .github/workflows/*.yml whose ``uses:`` is
actions/checkout.

Output: a failure naming the workflow, the job and the step for a checkout
without ``persist-credentials: false``.

Semantics: actions/checkout writes the job's GITHUB_TOKEN into the
checkout's .git/config as an extraheader unless persist-credentials is
false, and every later step of the job -- the test suite, a build backend,
a dependency's install script -- can read it from there (zizmor
``artipacked``; OpenSSF Scorecard Token-Permissions). No workflow here
pushes, fetches after the checkout or calls git with the token: the jobs
that need history ask for it through ``fetch-depth``, which the checkout
action performs before it drops the credential, and the one job that
writes to GitHub (release.yml ``github``) has no checkout and passes its
token to ``gh`` through the environment. So there is no exception list;
a job that needs to push would have to add one here, with its reason.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
_FILES = sorted((ROOT / ".github" / "workflows").glob("*.yml"))


def _checkouts(path: Path):
    wf = yaml.safe_load(path.read_text(encoding="utf-8"))
    for job_name, job in (wf.get("jobs") or {}).items():
        for index, step in enumerate(job.get("steps") or []):
            if str(step.get("uses", "")).startswith("actions/checkout@"):
                yield job_name, index, step


def test_the_workflows_check_out_at_all():
    """A rename of the action would make the check below pass on nothing."""
    assert sum(1 for p in _FILES for _ in _checkouts(p)) >= 20


@pytest.mark.parametrize("path", _FILES, ids=[p.name for p in _FILES])
def test_no_checkout_persists_its_credentials(path):
    kept = [f"{job} step {index} ({step.get('name', 'unnamed')})"
            for job, index, step in _checkouts(path)
            if (step.get("with") or {}).get("persist-credentials") is not False]
    assert not kept, (
        f"{path.name}: checkout leaves the token in .git/config; "
        f"add `persist-credentials: false`: {kept}")
