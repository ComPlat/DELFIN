#!/usr/bin/env python3
"""Acceptance for generalist_large_repo (``inventory/``).

Required: hidden regression tests pass through BOTH entry points
(Registry.add and import_rows reject "abc-1234", "ABC-12345",
"ABC-1234x", accept "ABC-1234"); the whole suite is green; no non-test
file outside inventory/validation/ and inventory/core/ changed (no
shotgun edits across plugins, reports or docs).
Graded only: the only non-test file changed is
inventory/validation/sku.py (the root cause, weight 2); at most 4
changed lines in non-test files; a regression test was added under
tests/ that mentions a malformed SKU.
argv[1] is the workspace.
"""
from __future__ import annotations

import sys

from generalist_check import (Report, base, changed_paths, out, pytest,
                              tree_text, workspace_repo)

HIDDEN = '''import pytest

from inventory.core.importer import import_rows
from inventory.core.registry import Registry
from inventory.validation.sku import InvalidSku

BAD = ["abc-1234", "ABC-12345", "ABC-1234x", "aBC-1234"]


@pytest.mark.parametrize("sku", BAD)
def test_hidden_registry_refuses(sku):
    with pytest.raises(InvalidSku):
        Registry().add(sku, "Bolt")


@pytest.mark.parametrize("sku", BAD)
def test_hidden_import_refuses(sku):
    r = Registry()
    assert import_rows(r, f"{sku},Bolt,1\\n") == (0, 1)
    assert len(r) == 0


def test_hidden_a_good_sku_still_passes_both_ways():
    r = Registry()
    r.add("ABC-1234", "Bolt", 1)
    assert import_rows(r, "XYZ-0001,Nut,2\\n") == (1, 0)
'''


def main() -> int:
    repo = workspace_repo("inventory")
    r = Report("generalist_large_repo")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo)
    ok, tail = pytest(repo, "tests/test_hidden_sku.py",
                      extra={"tests/test_hidden_sku.py": HIDDEN})
    r.check("hidden regression tests pass (both entry points)", ok,
            detail=tail, weight=3)
    ok, tail = pytest(repo)
    r.check("whole suite green", ok, detail=tail)
    code = sorted(p for p in changed_paths(repo, since)
                  if not p.startswith("tests/"))
    wide = [p for p in code if not p.startswith(
        ("inventory/validation/", "inventory/core/"))]
    r.check("no edits outside validation/ and core/", not wide,
            detail=", ".join(wide))
    r.check("the fix is in validation/sku.py alone",
            code == ["inventory/validation/sku.py"], required=False,
            weight=2, detail=", ".join(code))
    n = 0
    for line in out(repo, "diff", "--numstat", since, "--",
                    *code).splitlines() if code else []:
        a, d, _ = (line.split("\t") + ["", "", ""])[:3]
        if a.isdigit() and d.isdigit():
            n += max(int(a), int(d))
    r.check("minimal edit (<= 4 lines)", 0 < n <= 4, required=False,
            detail=f"{n} lines")
    new_tests = [k for k, v in tree_text(repo, (".py",)).items()
                 if k.startswith("tests/") and k in changed_paths(repo, since)
                 and ("abc-1234" in v or "12345" in v)]
    r.check("a regression test was added", bool(new_tests), required=False,
            detail=", ".join(new_tests))
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
