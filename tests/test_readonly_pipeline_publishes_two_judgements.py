"""The two judgements readonly_pipeline publishes for the bash gate.

Input: one shell segment per case. Output: the paths a segment writes
INTO, and whether a ``sort`` line carries a flag that writes a file or
runs a program.

Both answer questions the gate could not answer from its own patterns,
and both live here rather than in the gate so that one module owns the
reading of a shell segment. The gate's use of them is measured in
tests/test_a_write_that_leaves_the_workspace_is_asked_about.py.
"""
from __future__ import annotations

import pytest

from delfin.agent import readonly_pipeline


@pytest.mark.parametrize("cmd,targets", [
    ("ls > f", ["f"]),
    ("ls >> log", ["log"]),
    ("ls 3>f", ["f"]),
    ("ls 2> err.txt", ["err.txt"]),
    ("ls >&f", ["f"]),
    ("ls &>f", ["f"]),
    ("ls &>> f", ["f"]),
    ("cat <>f", ["f"]),
    ("echo x > ~/.bashrc", ["~/.bashrc"]),
    ("ls > a.txt 2> b.txt", ["a.txt", "b.txt"]),
    # a stream merge writes no file
    ("grep x f 2>&1 | tail -5", []),
    ("ls >&2", []),
    ("ls -la", []),
])
def test_a_stream_merge_is_not_a_file_redirect(cmd, targets):
    """Every spelling that opens a FILE is a target; the forms that only
    rewire the process's own streams are not."""
    assert readonly_pipeline.file_redirect_targets(cmd) == targets, cmd


@pytest.mark.parametrize("cmd,writes", [
    ("sort --compress-program=sh x", True),   # runs sh
    ("sort --compress-program sh x", True),
    ("sort --out=f x", True),                 # getopt accepts the prefix
    ("sort --output=f x", True),
    ("sort -o f x", True),
    ("sort -ofile x", True),                  # short cluster
    ("sort -T /tmp x", True),                 # temp dir anywhere
    ("sort -u f", False),
    ("sort -k2 -r f", False),
    ("sort --numeric-sort f", False),
    ("sort --reverse --unique f", False),
    ("grep x f", False),                      # not a sort line at all
    ("ls", False),
])
def test_sort_is_judged_by_what_may_read_not_by_what_may_write(cmd, writes):
    """An allow-list, because a deny-list cannot enumerate what GNU
    getopt accepts: ``--out=f`` is ``--output``, and
    ``--compress-program`` was named by no pattern at all."""
    assert readonly_pipeline.sort_flags_write_or_run(cmd) is writes, cmd


def test_the_published_judgements_do_not_change_is_safe():
    """Control: ``is_safe`` keeps its own verdicts. The judgements are
    published beside it, not instead of it."""
    assert readonly_pipeline.is_safe("ls -la") is True
    assert readonly_pipeline.is_safe("grep x f 2>&1 | tail -5") is True
    assert readonly_pipeline.is_safe("ls > f") is False
    assert readonly_pipeline.is_safe("sort --compress-program=sh x") is False
    assert readonly_pipeline.is_safe("") is False
