#!/usr/bin/env python3
"""Acceptance for generalist_go_red_to_green (``netcfg/``).

Required: ``go test ./...`` is green; port_test.go is byte-identical to
the given one; hidden cases pass (65535 and 1 accepted; 70000, -1, ""
and " 80" rejected).
Graded only: only port.go changed outside *_test.go; the fix is small
(<= 2 changed lines).
Without go the script exits 3 (not run).
argv[1] is the workspace.
"""
from __future__ import annotations

import sys

from generalist_check import (Report, base, blob_at, changed_paths, out,
                              read, tool_test, toolchain_missing, which,
                              workspace_repo)

HIDDEN = '''package netcfg

import "testing"

func TestHiddenBounds(t *testing.T) {
\tfor _, s := range []string{"65535", "1", "443"} {
\t\tif _, err := ParsePort(s); err != nil {
\t\t\tt.Fatalf("%q: %v", s, err)
\t\t}
\t}
\tfor _, s := range []string{"65536", "70000", "-1", "", " 80"} {
\t\tif _, err := ParsePort(s); err == nil {
\t\t\tt.Fatalf("%q must be rejected", s)
\t\t}
\t}
}
'''


def main() -> int:
    repo = workspace_repo("netcfg")
    go = which("go")
    if not go:
        return toolchain_missing("go")
    r = Report("generalist_go_red_to_green")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo)
    ok, tail = tool_test(repo, [go, "test", "./..."])
    r.check("go test green", ok, detail=tail)
    r.check("the given test file is unchanged",
            read(repo, "port_test.go")
            == (blob_at(repo, since, "port_test.go") or "\0"), weight=2)
    ok, tail = tool_test(repo, [go, "test", "-run", "TestHidden", "./..."],
                         extra={"hidden_port_test.go": HIDDEN})
    r.check("hidden cases pass", ok, detail=tail, weight=2)
    code = sorted(p for p in changed_paths(repo, since)
                  if not p.endswith("_test.go"))
    r.check("only port.go changed outside tests", code == ["port.go"],
            required=False, detail=", ".join(code))
    n = 0
    for line in (out(repo, "diff", "--numstat", since, "--", *code)
                 .splitlines() if code else []):
        a, d, _ = (line.split("\t") + ["", "", ""])[:3]
        if a.isdigit() and d.isdigit():
            n += max(int(a), int(d))
    r.check("the fix is small (<= 2 lines)", 0 < n <= 2, required=False,
            detail=f"{n} lines")
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
