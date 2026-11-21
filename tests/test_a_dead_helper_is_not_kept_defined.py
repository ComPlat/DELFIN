"""A dead helper nobody reached is not kept defined.

Measured behaviour (behaviour n5f protects): the four home/host identity
wrappers that no code, config or export channel ever reached stay deleted --
``swarm/policy.home_owner``, ``swarm/policy._resolve_home`` (its only
helper), ``service/attention_feed.this_host`` and
``harness_rows.is_harness_text``. A deletion has no callable entry point to
drive, so the real entry point is the product's source tree on disk: the
scan below reads the shipped files themselves and reports every remaining
definition by file and number. Before the deletion this FAILS naming the
four definitions (the control run at origin/main records that); after it,
PASSES on zero.

Covered loss-rebound: if a future change needs one of these helpers again it
must re-introduce it WITH a caller -- this test then fails and points at the
exact spot, which is the review moment.

What this does NOT measure (honest limits): it pins the (file, function
name) pairs of THIS removal. It is not a general dead-code detector -- a
different unreached function with a common name, or one of these names in
another file, is not found by it. Reachability itself was established by
the investigation (whole-tree and AST scans over every spelling and export
channel), not by this test. The live sibling
``session_store._this_host`` is deliberately out of the set: it has seven
callers and must remain.
"""
from __future__ import annotations

import ast
from pathlib import Path

#: (file relative to the repository root, function name) that must stay gone.
GONE_HELPERS = (
    ("delfin/agent/swarm/policy.py", "home_owner"),
    ("delfin/agent/swarm/policy.py", "_resolve_home"),
    ("delfin/agent/service/attention_feed.py", "this_host"),
    ("delfin/agent/harness_rows.py", "is_harness_text"),
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def test_no_unreachable_home_or_host_helper_is_still_defined():
    root = _repo_root()
    found = []
    for rel, name in GONE_HELPERS:
        path = root / rel
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) \
                    and node.name == name:
                found.append(f"{rel}:{node.lineno}: def {name}")
    assert not found, (
        "home/host identity helpers that nothing reaches are defined again "
        "(delete them, or give every one of them a caller): "
        + "; ".join(found))