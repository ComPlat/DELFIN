#!/usr/bin/env python3
"""Acceptance for generalist_long_rust_spec (``numerals/``).

Required: ``cargo test`` is green; the hidden integration test passes
(the SPEC examples, every n in 1..=3999 round-trips, out-of-range is
None, every listed non-canonical or malformed input is None); no
``todo!`` is left in src/; the attempt added tests that call both
functions and at least one expects None; README's Usage section names
both functions and no longer says TODO.
Graded only: SPEC.md unchanged; no dependency added; nothing outside
src/, tests/, README.md, Cargo.toml, CHANGELOG.md changed.
Without cargo the script exits 3 (not run).
argv[1] is the workspace.
"""
from __future__ import annotations

import re
import sys

from generalist_check import (Report, base, changed_paths, read, tool_test,
                              toolchain_missing, tree_text, which,
                              workspace_repo)

HIDDEN = r'''use numerals::{from_roman, to_roman};

#[test]
fn hidden_spec_examples() {
    for (n, s) in [(1, "I"), (4, "IV"), (9, "IX"), (14, "XIV"), (40, "XL"),
                   (90, "XC"), (400, "CD"), (900, "CM"), (1994, "MCMXCIV"),
                   (2026, "MMXXVI"), (3999, "MMMCMXCIX")] {
        assert_eq!(to_roman(n).as_deref(), Some(s), "to_roman({n})");
        assert_eq!(from_roman(s), Some(n), "from_roman({s})");
    }
}

#[test]
fn hidden_round_trip_over_the_whole_range() {
    for n in 1..=3999u32 {
        let s = to_roman(n).expect("in range");
        assert_eq!(from_roman(&s), Some(n), "{s}");
    }
}

#[test]
fn hidden_out_of_range() {
    assert_eq!(to_roman(0), None);
    assert_eq!(to_roman(4000), None);
}

#[test]
fn hidden_strict_parsing() {
    for s in ["", "xiv", "IIII", "VV", "IC", "XM", "VX", "MMMM", " X", "X ",
              "IIX", "XIIII", "LL", "DD", "CMC", "A"] {
        assert_eq!(from_roman(s), None, "{s:?} must be rejected");
    }
}
'''


def main() -> int:
    repo = workspace_repo("numerals")
    cargo = which("cargo")
    if not cargo:
        return toolchain_missing("cargo")
    r = Report("generalist_long_rust_spec")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo)
    ok, tail = tool_test(repo, [cargo, "test", "--offline", "--quiet"])
    r.check("cargo test green", ok, detail=tail, weight=2)
    ok, tail = tool_test(repo, [cargo, "test", "--offline", "--quiet",
                                "--test", "hidden_spec"],
                         extra={"tests/hidden_spec.rs": HIDDEN})
    r.check("hidden spec tests pass", ok, detail=tail, weight=4)
    rs = tree_text(repo, (".rs",))
    r.check("no todo! left in src/", not any(
        "todo!(" in v for k, v in rs.items() if k.startswith("src/")))
    tests = "\n".join(v for v in rs.values() if "#[test]" in v)
    r.check("own tests call both functions and expect a None",
            "to_roman(" in tests and "from_roman(" in tests
            and "None" in tests)
    readme = read(repo, "README.md")
    usage = re.search(r"^## Usage\s*$(.*?)(?=^## |\Z)", readme, re.M | re.S)
    body = usage.group(1) if usage else ""
    r.check("README Usage shows both functions, no TODO",
            "to_roman" in body and "from_roman" in body
            and "TODO" not in body, detail=body.strip()[:120])
    changed = changed_paths(repo, since)
    r.check("SPEC.md unchanged", "SPEC.md" not in changed, required=False)
    deps = re.search(r"\[dependencies\]\s*(.*?)(\n\[|\Z)",
                     read(repo, "Cargo.toml"), re.S)
    r.check("no dependency added", bool(deps) and not deps.group(1).strip(),
            required=False)
    stray = sorted(p for p in changed if not p.startswith(
        ("src/", "tests/", "README.md", "Cargo.toml", "CHANGELOG.md")))
    r.check("no unrelated changes", not stray, required=False,
            detail=", ".join(stray))
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
