#!/usr/bin/env python3
"""Acceptance for generalist_rust_feature (``unitconv/``).

Required: ``cargo test`` is green; a hidden integration test
(tests/hidden_duration.rs, using ``unitconv::parse_duration``) passes;
the attempt added a test that calls ``parse_duration`` (inline
#[cfg(test)] or under tests/); src/temperature.rs is unchanged.
Graded only: no dependency added to Cargo.toml; nothing outside src/,
tests/, README.md, CHANGELOG.md, Cargo.toml changed.
Without cargo the script exits 3 (not run).
argv[1] is the workspace.
"""
from __future__ import annotations

import re
import sys

from generalist_check import (Report, base, blob_at, changed_paths, read,
                              tool_test, toolchain_missing, tree_text, which,
                              workspace_repo)

HIDDEN = r'''use unitconv::parse_duration;

#[test]
fn hidden_valid_forms() {
    assert_eq!(parse_duration("1h30m15s"), Some(5415));
    assert_eq!(parse_duration("45s"), Some(45));
    assert_eq!(parse_duration("2h"), Some(7200));
    assert_eq!(parse_duration("90m"), Some(5400));
    assert_eq!(parse_duration("1h5s"), Some(3605));
    assert_eq!(parse_duration("0s"), Some(0));
}

#[test]
fn hidden_invalid_forms() {
    for s in ["", "h", "1x", "30m1h", "1h1h", "1h 30m", "-1s", "1.5h",
              "15", "s5", "99999999999999999999s"] {
        assert_eq!(parse_duration(s), None, "{s:?} must be rejected");
    }
}

#[test]
fn hidden_overflow_is_none_not_a_panic() {
    assert_eq!(parse_duration("5124095576030432h"), None);
}
'''


def main() -> int:
    repo = workspace_repo("unitconv")
    cargo = which("cargo")
    if not cargo:
        return toolchain_missing("cargo")
    r = Report("generalist_rust_feature")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo)
    ok, tail = tool_test(repo, [cargo, "test", "--offline", "--quiet"])
    r.check("cargo test green", ok, detail=tail)
    ok, tail = tool_test(repo, [cargo, "test", "--offline", "--quiet",
                                "--test", "hidden_duration"],
                         extra={"tests/hidden_duration.rs": HIDDEN})
    r.check("hidden parse_duration tests pass", ok, detail=tail, weight=3)
    own = [k for k, v in tree_text(repo, (".rs",)).items()
           if "parse_duration" in v and ("#[test]" in v)]
    r.check("a test of parse_duration was added", bool(own),
            detail=", ".join(own))
    r.check("src/temperature.rs unchanged",
            read(repo, "src/temperature.rs")
            == (blob_at(repo, since, "src/temperature.rs") or "\0"))
    deps = re.search(r"\[dependencies\]\s*(.*?)(\n\[|\Z)",
                     read(repo, "Cargo.toml"), re.S)
    r.check("no dependency added", bool(deps) and not deps.group(1).strip(),
            required=False)
    allowed = ("src/", "tests/", "README.md", "CHANGELOG.md", "Cargo.toml")
    stray = sorted(p for p in changed_paths(repo, since)
                   if not p.startswith(allowed))
    r.check("nothing outside the feature changed", not stray,
            required=False, detail=", ".join(stray))
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
