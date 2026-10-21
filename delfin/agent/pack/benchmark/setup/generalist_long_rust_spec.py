#!/usr/bin/env python3
"""Setup for generalist_long_rust_spec: implement a written spec in Rust.

``numerals/`` -- a crate whose SPEC.md defines ``to_roman`` and
``from_roman`` (range 1..=3999, canonical subtractive form only, strict
parsing: no lower case, no "IIII", no "VX", no "IC", no empty string).
src/lib.rs holds the two signatures with ``todo!()`` bodies and no
tests; README.md has a Usage section that says "TODO". The task:
implement the spec, add tests, replace the README's TODO with usage of
both functions.
A host without cargo exits 3.
Precondition: the crate builds, and ``cargo test`` runs zero tests.
argv[1] is the workspace.
"""
from __future__ import annotations

from generalist_repo import (commit, expect, main, mark_base, new_repo,
                             require, tool_rc, write)

SPEC = """# numerals: specification

Two functions in the crate root.

## `to_roman(n: u32) -> Option<String>`

- Defined for 1 to 3999 inclusive; any other `n` gives `None`.
- Upper-case letters I V X L C D M, canonical subtractive notation:
  4 = IV, 9 = IX, 40 = XL, 90 = XC, 400 = CD, 900 = CM; never four of a
  kind (4 is not IIII).
- Examples: 1 -> "I", 4 -> "IV", 14 -> "XIV", 1994 -> "MCMXCIV",
  2026 -> "MMXXVI", 3999 -> "MMMCMXCIX".

## `from_roman(s: &str) -> Option<u32>`

- The inverse of `to_roman`: returns `Some(n)` exactly when
  `to_roman(n) == Some(s)`.
- Everything else is `None`: the empty string, lower case ("xiv"),
  non-canonical forms ("IIII", "VV", "IC", "XM", "VX", "MMMM"),
  surrounding whitespace, any other character.

## Tests and docs

- Unit or integration tests cover both functions, including at least
  one rejected input.
- README.md's Usage section shows both functions.
"""

LIB = """//! Roman numerals. See SPEC.md.

/// `n` (1..=3999) in canonical Roman notation.
pub fn to_roman(n: u32) -> Option<String> {
    let _ = n;
    todo!("see SPEC.md")
}

/// The value of a canonical Roman numeral.
pub fn from_roman(s: &str) -> Option<u32> {
    let _ = s;
    todo!("see SPEC.md")
}
"""

FILES = {
    ".gitignore": "/target\nCargo.lock\n",
    "Cargo.toml": ("[package]\nname = \"numerals\"\nversion = \"0.1.0\"\n"
                   "edition = \"2021\"\n\n[dependencies]\n"),
    "SPEC.md": SPEC,
    "README.md": ("# numerals\n\nRoman numerals for Rust. The behaviour is "
                  "defined in SPEC.md.\n\n## Usage\n\nTODO\n\n"
                  "## License\n\nMIT\n"),
    "src/lib.rs": LIB,
}


def build(workspace: str):
    cargo = require("cargo")
    repo = new_repo(workspace, "numerals")
    write(repo, FILES)
    mark_base(repo, commit(repo, "numerals: spec and signatures"))
    rc, out = tool_rc(repo, [cargo, "test", "--offline"])
    expect("the crate builds and has no tests yet",
           rc == 0 and "running 0 tests" in out, out)
    return repo


if __name__ == "__main__":
    main(build)
