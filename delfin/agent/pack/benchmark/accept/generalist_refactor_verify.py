#!/usr/bin/env python3
"""Acceptance for generalist_refactor (``shoplib/``).

Required: the whole suite is green; hidden tests pass
(``from shoplib.pricing import order_total`` and
``from shoplib import order_total`` name the same function with the
old behaviour; shoplib.utils has neither ``calc_total`` nor
``order_total``); the string ``calc_total`` occurs in no file of the
work tree (code, tests, docs -- no alias, no leftover).
Graded only: ``format_price`` and ``clamp`` still live in utils (nothing
else was moved); the number of test functions did not drop.
argv[1] is the workspace.
"""
from __future__ import annotations

import re
import sys

from generalist_check import (Report, base, out, pytest, read, tree_text,
                              workspace_repo)

HIDDEN = '''import shoplib
import shoplib.utils as utils
from shoplib.pricing import order_total


def test_hidden_one_function_two_import_paths():
    assert shoplib.order_total is order_total


def test_hidden_behaviour_is_kept():
    assert order_total([(250, 2), (100, 1)]) == 600
    assert order_total([(250, 2), (100, 1)], tax_rate=0.19) == 714
    assert order_total([(1, 1)], tax_rate=0.5) == 2


def test_hidden_utils_no_longer_holds_it():
    assert not hasattr(utils, "calc_total")
    assert not hasattr(utils, "order_total")
'''

_TEST_DEF = re.compile(r"^\s*def test_", re.M)


def main() -> int:
    repo = workspace_repo("shoplib")
    r = Report("generalist_refactor")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo)
    ok, tail = pytest(repo)
    r.check("whole suite green", ok, detail=tail, weight=2)
    ok, tail = pytest(repo, "tests/test_hidden_pricing.py",
                      extra={"tests/test_hidden_pricing.py": HIDDEN})
    r.check("hidden tests pass", ok, detail=tail, weight=2)
    texts = tree_text(repo, (".py", ".md", ".toml", ".txt", ".cfg"))
    left = sorted(k for k, v in texts.items() if "calc_total" in v)
    r.check("no reference to calc_total is left", not left, weight=2,
            detail=", ".join(left))
    utils = read(repo, "shoplib/utils.py")
    r.check("format_price and clamp stay in utils",
            "def format_price" in utils and "def clamp" in utils,
            required=False)
    before = sum(len(_TEST_DEF.findall(out(repo, "show", f"{since}:{p}")))
                 for p in out(repo, "ls-tree", "-r", "--name-only", since,
                              "tests").splitlines())
    after = sum(len(_TEST_DEF.findall(v)) for k, v in texts.items()
                if k.startswith("tests/") and k.endswith(".py"))
    r.check("no test was dropped", after >= before, required=False,
            detail=f"{before} -> {after} test functions")
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
