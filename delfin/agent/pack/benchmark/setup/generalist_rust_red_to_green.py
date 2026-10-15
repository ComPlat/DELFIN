#!/usr/bin/env python3
"""Setup for generalist_rust_red_to_green: one failing Rust test.

``stats/`` -- ``median(values: &[f64]) -> Option<f64>`` returns the upper
middle element for an even count instead of the mean of the two middle
ones. tests/median.rs has three tests; exactly the even-count one fails
(1,2,3,4 -> 3.0, expected 2.5). A host without cargo exits 3.
Precondition: exactly that test fails.
argv[1] is the workspace.
"""
from __future__ import annotations

from generalist_repo import (commit, expect, main, mark_base, new_repo,
                             require, tool_rc, write)

LIB = '''//! Small statistics.

/// The median of `values`, or `None` when it is empty.
pub fn median(values: &[f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("no NaN"));
    let mid = sorted.len() / 2;
    Some(sorted[mid])
}

/// The arithmetic mean of `values`, or `None` when it is empty.
pub fn mean(values: &[f64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    Some(values.iter().sum::<f64>() / values.len() as f64)
}
'''

TEST = '''use stats::median;

#[test]
fn odd_count_takes_the_middle() {
    assert_eq!(median(&[3.0, 1.0, 2.0]), Some(2.0));
}

#[test]
fn empty_has_no_median() {
    assert_eq!(median(&[]), None);
}

#[test]
fn even_count_averages_the_two_middle_values() {
    assert_eq!(median(&[1.0, 2.0, 3.0, 4.0]), Some(2.5));
}
'''

FILES = {
    ".gitignore": "/target\nCargo.lock\n",
    "Cargo.toml": ("[package]\nname = \"stats\"\nversion = \"0.2.0\"\n"
                   "edition = \"2021\"\n\n[dependencies]\n"),
    "README.md": "# stats\n\nSmall statistics. `cargo test`.\n",
    "src/lib.rs": LIB,
    "tests/median.rs": TEST,
}


def build(workspace: str):
    cargo = require("cargo")
    repo = new_repo(workspace, "stats")
    write(repo, FILES)
    mark_base(repo, commit(repo, "Add median and mean"))
    rc, out = tool_rc(repo, [cargo, "test", "--offline", "--test", "median"])
    expect("exactly the even-count test fails",
           rc != 0 and "2 passed; 1 failed" in out
           and "even_count_averages_the_two_middle_values" in out, out)
    return repo


if __name__ == "__main__":
    main(build)
