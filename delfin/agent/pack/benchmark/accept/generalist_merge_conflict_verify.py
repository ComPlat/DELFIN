#!/usr/bin/env python3
"""Acceptance for generalist_merge_conflict (``configlib/``).

Required: no merge in progress; the current branch is main; main's
first-parent history contains a merge commit whose parents are exactly
(M, F), the two tips the setup recorded (so the history was neither
rebased, reset nor squashed); no tracked file has conflict markers; the
suite is green (both sides' tests); hidden tests of the combination pass
(comments skipped AND environment override, in one file).
Graded only: HEAD is that merge itself (no extra commits); the work
tree is clean; feature/env-override still points at F.
argv[1] is the workspace.
"""
from __future__ import annotations

import sys

from generalist_check import (Report, base, files_with_markers, out, pytest,
                              read, workspace_repo)

HIDDEN = '''from configlib.loader import load


def test_hidden_both_behaviours_at_once(tmp_path):
    p = tmp_path / "a.conf"
    p.write_text("# db\\n\\nhost=example.org\\nport=8080\\n", encoding="utf-8")
    assert load(p, env={"CONFIGLIB_PORT": "9090"}) == {
        "host": "example.org", "port": "9090"}


def test_hidden_a_commented_key_is_not_resurrected_by_env(tmp_path):
    p = tmp_path / "a.conf"
    p.write_text("#port=1\\nhost=h\\n", encoding="utf-8")
    assert load(p, env={}) == {"host": "h"}
'''


def main() -> int:
    repo = workspace_repo("configlib")
    r = Report("generalist_merge_conflict")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    m = base(repo)
    f = read(repo, ".git/delfin-bench-feature").strip()
    r.check("no merge in progress", not (repo / ".git" / "MERGE_HEAD").exists())
    r.check("on main", out(repo, "branch", "--show-current") == "main")
    merge = ""
    for line in out(repo, "rev-list", "--first-parent", "--parents",
                    "HEAD").splitlines():
        shas = line.split()
        if shas[1:] == [m, f]:
            merge = shas[0]
            break
    r.check("a merge commit with parents (main, feature)", bool(merge),
            weight=2)
    markers = files_with_markers(repo)
    r.check("no conflict markers", not markers, detail=", ".join(markers))
    ok, tail = pytest(repo)
    r.check("whole suite green", ok, detail=tail, weight=2)
    ok, tail = pytest(repo, "tests/test_hidden_merge.py",
                      extra={"tests/test_hidden_merge.py": HIDDEN})
    r.check("hidden combined tests pass", ok, detail=tail, weight=2)
    r.check("HEAD is the merge", bool(merge)
            and out(repo, "rev-parse", "HEAD") == merge, required=False)
    r.check("work tree clean",
            not out(repo, "status", "--porcelain", "--untracked-files=no"),
            required=False)
    r.check("feature branch untouched",
            out(repo, "rev-parse", "-q", "--verify",
                "refs/heads/feature/env-override") == f, required=False)
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
