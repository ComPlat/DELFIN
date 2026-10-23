"""The server that publishes the dashboard cannot resolve to a release with a
known, fixed advisory.

Input: ``project.dependencies`` in pyproject.toml.

Output: a failure when the dependency on voila, jupyter-server or tornado
is missing, admits a release below the first one that fixes the advisories
listed here, or has no upper bound (the project's rule for every entry).

Semantics: voila brings jupyter-server and tornado in transitively, with
floors far below these (voila 0.5.13: jupyter-server>=1.18; jupyter-server
2.21.1: tornado>=6.2), so without a direct requirement a resolver that has
an old copy at hand keeps it. The floors, from the GitHub advisory
database, read 2026-10-16:

- jupyter-server 2.20.0: GHSA-fcw5-x6j4-ccmp (stored XSS in the nbconvert
  handlers, <= 2.19.0); 2.18.0 fixed GHSA-5789-5fc7-67v3 (path traversal
  out of the root directory), GHSA-24qx-w28j-9m6p (CORS origin check),
  GHSA-5mrq-x3x5-8v8f (cookies valid after a password change) and
  GHSA-qh7q-6qm3-653w (open redirect).
- tornado 6.5.9: GHSA-c2m8-h5v5-343r (StaticFileHandler follows symlinks
  out of the static root), GHSA-3hv7-mjh2-fv65 (query-string argument
  count), GHSA-chx6-46f5-w4vp (curl client response size), all <= 6.5.8;
  6.5.8 fixed GHSA-mpf4-983q-p7j4 and GHSA-8423-8fgw-73vq.
- voila 0.5.6: GHSA-2q59-h24c-w6fg (local file inclusion, 0.5.0a0 to
  0.5.5).

Raising a floor when a new advisory is fixed is an edit here and in
pyproject.toml; audit.yml reports the advisory nightly until then.
"""

from __future__ import annotations

import tomllib
from pathlib import Path

import pytest
from packaging.requirements import Requirement

ROOT = Path(__file__).resolve().parents[1]

#: name -> (last release with a known advisory, first patched release)
FLOORS = {
    "jupyter-server": ("2.19.0", "2.20.0"),
    "tornado": ("6.5.8", "6.5.9"),
    "voila": ("0.5.5", "0.5.6"),
}


def _requirements() -> dict[str, Requirement]:
    data = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    out = {}
    for line in data["project"]["dependencies"]:
        req = Requirement(line)
        out[req.name.lower().replace("_", "-")] = req
    return out


@pytest.mark.parametrize("name", sorted(FLOORS))
def test_the_floor_is_the_patched_release(name):
    reqs = _requirements()
    assert name in reqs, f"{name} is not required directly; its floor is voila's"
    vulnerable, patched = FLOORS[name]
    spec = reqs[name].specifier
    assert not spec.contains(vulnerable, prereleases=True), (
        f"{name}{spec} admits {vulnerable}")
    assert spec.contains(patched), f"{name}{spec} refuses the patched {patched}"
    assert any(s.operator in ("<", "<=") for s in spec), f"{name}{spec} has no ceiling"
