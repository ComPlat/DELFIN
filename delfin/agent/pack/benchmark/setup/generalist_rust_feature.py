#!/usr/bin/env python3
"""Setup for generalist_rust_feature: a feature in a Cargo crate.

``unitconv/`` -- one commit on main: a library crate with no
dependencies (src/lib.rs, src/temperature.rs with inline #[cfg(test)]
tests). The task adds ``parse_duration`` in src/duration.rs, re-exported
from the crate root. A host without cargo exits 3 (toolchain
unavailable -> unmeasured). Builds run offline (CARGO_NET_OFFLINE).
Precondition: ``cargo test`` is green.
argv[1] is the workspace.
"""
from __future__ import annotations

from generalist_repo import (commit, expect, main, mark_base, new_repo,
                             require, tool_rc, write)

FILES = {
    ".gitignore": "/target\nCargo.lock\n",
    "Cargo.toml": ("[package]\nname = \"unitconv\"\nversion = \"0.1.0\"\n"
                   "edition = \"2021\"\n\n[dependencies]\n"),
    "README.md": ("# unitconv\n\nUnit conversions. Run the tests with "
                  "`cargo test`.\n"),
    "src/lib.rs": ("//! Unit conversions.\n\npub mod temperature;\n\n"
                   "pub use temperature::celsius_to_fahrenheit;\n"),
    "src/temperature.rs": (
        "//! Temperatures.\n\n"
        "/// Degrees Celsius to degrees Fahrenheit.\n"
        "pub fn celsius_to_fahrenheit(c: f64) -> f64 {\n"
        "    c * 9.0 / 5.0 + 32.0\n}\n\n"
        "#[cfg(test)]\nmod tests {\n    use super::*;\n\n"
        "    #[test]\n    fn freezing_and_boiling() {\n"
        "        assert_eq!(celsius_to_fahrenheit(0.0), 32.0);\n"
        "        assert_eq!(celsius_to_fahrenheit(100.0), 212.0);\n"
        "    }\n}\n"),
}


def build(workspace: str):
    cargo = require("cargo")             # exit 3 before anything is built
    repo = new_repo(workspace, "unitconv")
    write(repo, FILES)
    mark_base(repo, commit(repo, "Add celsius_to_fahrenheit"))
    rc, out = tool_rc(repo, [cargo, "test", "--offline", "--quiet"])
    expect("cargo test is green", rc == 0, out)
    return repo


if __name__ == "__main__":
    main(build)
