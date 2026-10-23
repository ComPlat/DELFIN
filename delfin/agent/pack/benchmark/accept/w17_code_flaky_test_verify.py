#!/usr/bin/env python3
"""Acceptance for w17_code_flaky_test (``slowweb/``).

Required: whole suite green; hidden tests pin the real fix: the cache
takes a ``now=`` callable (defaulting to time.time) used for both
stamping and expiry checks, and the FLAKY TEST no longer sleeps — it
passes an injectable clock and runs with a fixed timeline. Graded
only: the stable test file's semantics survive (the plain API still
works with no ``now=`` argument); the diff stays small and inside
src/slowweb/ + tests/; deleting the flaky test instead of fixing the
mechanism fails the hidden tests.
argv[1] is the workspace.
"""
from __future__ import annotations

import sys

from w17_code_check import diff_scope_check
from generalist_check import (
    Report, base, pytest, read, tree_text, workspace_repo,
)

HIDDEN = '''
import time

from slowweb import Cache


def test_default_clock_still_works():
    c = Cache(ttl=60)
    c.put("k", "v")
    assert c.get("k") == "v"


def test_injected_clock_controls_expiry():
    t = {"v": 1000.0}
    c = Cache(ttl=10.0, now=lambda: t["v"])
    c.put("k", "v")
    assert c.get("k") == "v"
    t["v"] = 1011.0
    assert c.get("k") is None


def test_entry_stamped_with_the_injected_clock():
    t = {"v": 1000.0}
    c = Cache(ttl=10.0, now=lambda: t["v"])
    c.put("k", "v")
    t["v"] = 1005.0
    assert c.get("k") == "v"
    t["v"] = 1010.5
    assert c.get("k") is None
'''


def main() -> int:
    repo = workspace_repo("slowweb")
    r = Report("w17_code_flaky_test")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo) or "HEAD"
    ok, tail = pytest(repo)
    r.check("whole suite green", ok, detail=tail)
    ok, tail = pytest(repo, "tests/test_w17_hidden.py",
                      extra={"tests/test_w17_hidden.py": HIDDEN})
    r.check("hidden clock-injection tests pass", ok, detail=tail, weight=3)
    tests_text = "\n".join(v for k, v in tree_text(repo, (".py",)).items()
                           if k.startswith("tests/"))
    r.check("the flaky test no longer sleeps on wall time",
            "time.sleep" not in tests_text,
            detail=tests_text[:200])
    r.check("the flaky test still exists and tests expiry",
            "test_expired_entry_is_gone" in tests_text)
    r.check("src/slowweb/cache.py still stamps with the default clock",
            "time.time" in read(repo, "src/slowweb/cache.py"))
    diff_scope_check(r, repo,
                     ("src/slowweb/", "tests/", "README.md"), 40)
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
