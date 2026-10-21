#!/usr/bin/env python3
"""Acceptance for generalist_long_ts_symptom (``orders-ts/``).

Required: every test/*.test.ts passes; the given tests are unchanged;
hidden tests pass: amounts are computed in cents at every layer, not
only printed with two places (lineTotal(0.125 x 1) = 0.13, a two-line
invoice of 0.125 each is 0.31, its printed total "Total: 0.31 EUR",
roundAmount's default rounds to the cent).
Graded only: the change is at the cause -- src/config.ts is the only
non-test file changed (weight 2); at most 2 changed lines in non-test
files; report.ts unchanged; a regression test was added.
Without a Node that runs TypeScript the script exits 3 (not run).
argv[1] is the workspace.
"""
from __future__ import annotations

import sys

from generalist_check import (Report, base, blob_at, changed_paths,
                              node_prefix, node_test, out, read,
                              toolchain_missing, tree_text, workspace_repo)

HIDDEN = '''import { test } from "node:test";
import assert from "node:assert/strict";
import { roundAmount } from "../src/money.ts";
import { invoiceTotal, lineTotal } from "../src/invoice.ts";
import { totalLine } from "../src/report.ts";

test("hidden: the default rounds to the cent", () => {
  assert.equal(roundAmount(1.2345), 1.23);
  assert.equal(roundAmount(0.125), 0.13);
});

test("hidden: every layer works in cents", () => {
  assert.equal(lineTotal({ price: 0.125, qty: 1 }), 0.13);
  const lines = [{ price: 0.125, qty: 1 }, { price: 0.125, qty: 1 }];
  assert.equal(invoiceTotal(lines), 0.31);
  assert.equal(totalLine(lines), "Total: 0.31 EUR");
  assert.equal(invoiceTotal([{ price: 3.333, qty: 3 }]), 11.9);
});
'''

GIVEN = ("test/money.test.ts", "test/invoice.test.ts", "test/report.test.ts")


def main() -> int:
    repo = workspace_repo("orders-ts")
    if node_prefix() is None:
        return toolchain_missing("no node that runs TypeScript")
    r = Report("generalist_long_ts_symptom")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo)
    tests = sorted(k for k in tree_text(repo, (".test.ts",))
                   if k.startswith("test/") and "/" not in k[5:])
    ok, tail = node_test(repo, tests)
    r.check("npm test suite green", ok, detail=tail, weight=2)
    edited = [p for p in GIVEN
              if read(repo, p) != (blob_at(repo, since, p) or "\0")]
    r.check("the given tests are unchanged", not edited,
            detail=", ".join(edited))
    ok, tail = node_test(repo, ["test/hidden_cents.test.ts"],
                         extra={"test/hidden_cents.test.ts": HIDDEN})
    r.check("hidden tests: cents at every layer", ok, detail=tail, weight=4)
    changed = changed_paths(repo, since)
    code = sorted(p for p in changed if not p.startswith("test/"))
    r.check("fixed at the cause (src/config.ts only)",
            code == ["src/config.ts"], required=False, weight=2,
            detail=", ".join(code))
    n = 0
    for line in (out(repo, "diff", "--numstat", since, "--", *code)
                 .splitlines() if code else []):
        a, d, _ = (line.split("\t") + ["", "", ""])[:3]
        if a.isdigit() and d.isdigit():
            n += max(int(a), int(d))
    r.check("minimal change (<= 2 lines)", 0 < n <= 2, required=False,
            detail=f"{n} lines")
    r.check("report.ts untouched", "src/report.ts" not in changed,
            required=False)
    r.check("a regression test was added",
            any(p.startswith("test/") and p not in GIVEN for p in changed),
            required=False)
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
