#!/usr/bin/env python3
"""Acceptance for generalist_ts_feature (``invoice-ts/``).

Required: every test/*.test.ts passes under ``node --test`` (what
``npm test`` runs); hidden tests of ``splitCents`` imported from
src/index.ts pass; the attempt added a test file under test/ that calls
``splitCents``; test/money.test.ts and the ``test`` script in
package.json are unchanged.
Graded only: the attempt's own tests fail against a plausible wrong
implementation; nothing outside src/, test/, README.md, CHANGELOG.md
changed; no dependency was added to package.json.
Without a Node that runs TypeScript the script exits 3 (not run).
argv[1] is the workspace.
"""
from __future__ import annotations

import json
import re
import sys

from generalist_check import (Report, base, blob_at, changed_paths,
                              node_prefix, node_test, read,
                              toolchain_missing, tree_text, workspace_repo)

HIDDEN = '''import { test } from "node:test";
import assert from "node:assert/strict";
import { splitCents } from "../src/index.ts";

test("hidden: larger shares first, differing by at most one", () => {
  assert.deepEqual(splitCents(100, 3), [34, 33, 33]);
  assert.deepEqual(splitCents(10, 5), [2, 2, 2, 2, 2]);
  assert.deepEqual(splitCents(1, 3), [1, 0, 0]);
  assert.deepEqual(splitCents(0, 2), [0, 0]);
  assert.deepEqual(splitCents(7, 1), [7]);
});

test("hidden: the shares always sum to the total", () => {
  for (const [t, p] of [[999, 7], [1000, 3], [12345, 11], [2, 5]]) {
    const parts = splitCents(t, p);
    assert.equal(parts.length, p);
    assert.equal(parts.reduce((a, b) => a + b, 0), t);
    assert.ok(Math.max(...parts) - Math.min(...parts) <= 1);
  }
});

test("hidden: bad input is a RangeError", () => {
  assert.throws(() => splitCents(100, 0), RangeError);
  assert.throws(() => splitCents(100, 1.5), RangeError);
  assert.throws(() => splitCents(-1, 2), RangeError);
  assert.throws(() => splitCents(1.5, 2), RangeError);
});
'''

MUTANT = '''export function splitCents(total: number, parts: number): number[] {
  if (!Number.isInteger(parts) || parts < 1) throw new RangeError("parts");
  if (!Number.isInteger(total) || total < 0) throw new RangeError("total");
  return Array(parts).fill(Math.round(total / parts));
}
'''


def main() -> int:
    repo = workspace_repo("invoice-ts")
    if node_prefix() is None:
        return toolchain_missing("no node that runs TypeScript")
    r = Report("generalist_ts_feature")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo)
    tests = sorted(k for k in tree_text(repo, (".test.ts",))
                   if k.startswith("test/") and "/" not in k[5:])
    ok, tail = node_test(repo, tests)
    r.check("npm test suite green", ok, detail=tail)
    ok, tail = node_test(repo, ["test/hidden_split.test.ts"],
                         extra={"test/hidden_split.test.ts": HIDDEN})
    r.check("hidden splitCents tests pass", ok, detail=tail, weight=3)
    own = [k for k, v in tree_text(repo, (".ts",)).items()
           if k.startswith("test/") and "splitCents" in v]
    r.check("a test of splitCents was added under test/", bool(own),
            detail=", ".join(own))
    r.check("test/money.test.ts unchanged",
            read(repo, "test/money.test.ts")
            == (blob_at(repo, since, "test/money.test.ts") or "\0"))
    try:
        pkg = json.loads(read(repo, "package.json"))
        pkg0 = json.loads(blob_at(repo, since, "package.json") or "{}")
    except ValueError:
        pkg, pkg0 = {}, {"_": 1}
    r.check("the test script is unchanged",
            (pkg.get("scripts") or {}).get("test")
            == (pkg0.get("scripts") or {}).get("test"))
    # The mutant replaces whichever module defines splitCents, so the
    # attempt's tests meet it by any import path they use.
    defining = [k for k, v in tree_text(repo, (".ts",)).items()
                if k.startswith("src/") and re.search(
                    r"export\s+(function|const|let)\s+splitCents\b", v)]
    if own and defining:
        ok, tail = node_test(repo, own,
                             extra={k: MUTANT for k in defining})
        r.check("own tests reject a wrong implementation", not ok,
                required=False, detail=tail)
    allowed = ("src/", "test/", "README.md", "CHANGELOG.md", "package.json")
    stray = sorted(p for p in changed_paths(repo, since)
                   if not p.startswith(allowed))
    r.check("nothing outside the feature changed", not stray,
            required=False, detail=", ".join(stray))
    r.check("no dependency added",
            not pkg.get("dependencies") and not pkg.get("devDependencies"),
            required=False)
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
