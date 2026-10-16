#!/usr/bin/env python3
"""Acceptance for generalist_rust_red_to_green (``stats/``).

Required: ``cargo test`` is green; tests/median.rs is byte-identical to
the given one; hidden cases over other inputs pass (unsorted even
input, a single value, negative values).
Graded only: only src/lib.rs changed outside tests/; the fix is small
(<= 6 changed lines).
Without cargo the script exits 3 (not run).
argv[1] is the workspace.
"""
from __future__ import annotations

import sys

from generalist_check import (Report, base, blob_at, changed_paths, out,
                              read, tool_test, toolchain_missing, which,
                              workspace_repo)

HIDDEN = r'''use stats::{mean, median};

#[test]
fn hidden_even_counts() {
    assert_eq!(median(&[4.0, 1.0, 3.0, 2.0]), Some(2.5));
    assert_eq!(median(&[1.0, 1.0, 2.0, 10.0]), Some(1.5));
    assert_eq!(median(&[-2.0, 2.0]), Some(0.0));
}

#[test]
fn hidden_odd_counts_unchanged() {
    assert_eq!(median(&[5.0]), Some(5.0));
    assert_eq!(median(&[9.0, -1.0, 4.0, 7.0, 0.0]), Some(4.0));
}

#[test]
fn hidden_mean_unchanged() {
    assert_eq!(mean(&[1.0, 2.0, 6.0]), Some(3.0));
}
'''


def main() -> int:
    repo = workspace_repo("stats")
    cargo = which("cargo")
    if not cargo:
        return toolchain_missing("cargo")
    r = Report("generalist_rust_red_to_green")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo)
    ok, tail = tool_test(repo, [cargo, "test", "--offline", "--quiet"])
    r.check("cargo test green", ok, detail=tail)
    r.check("the given test file is unchanged",
            read(repo, "tests/median.rs")
            == (blob_at(repo, since, "tests/median.rs") or "\0"), weight=2)
    ok, tail = tool_test(repo, [cargo, "test", "--offline", "--quiet",
                                "--test", "hidden_median"],
                         extra={"tests/hidden_median.rs": HIDDEN})
    r.check("hidden cases pass", ok, detail=tail, weight=2)
    code = sorted(p for p in changed_paths(repo, since)
                  if not p.startswith("tests/"))
    r.check("only src/lib.rs changed outside tests/", code == ["src/lib.rs"],
            required=False, detail=", ".join(code))
    n = 0
    for line in (out(repo, "diff", "--numstat", since, "--", *code)
                 .splitlines() if code else []):
        a, d, _ = (line.split("\t") + ["", "", ""])[:3]
        if a.isdigit() and d.isdigit():
            n += max(int(a), int(d))
    r.check("the fix is small (<= 6 lines)", 0 < n <= 6, required=False,
            detail=f"{n} lines")
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
