#!/usr/bin/env python3
"""Split of the fast suite across CI runners, and the check that nothing fell out.

Input
    The node ids pytest collects for ``-m "not slow"`` (the universe), and
    ``test-durations.json`` next to this file: seconds per test file on a
    GitHub-hosted ubuntu-latest runner.

Output
    ``plan()``: test file -> shard index 1..N.
    ``verify()``: the problems found in N shard reports, empty when the
    shards together ran the universe exactly once.

Semantics
    The unit is the test FILE, so module-scoped fixtures and module state
    stay inside one process, as they do in an unsharded run. Files are
    assigned longest first to the shard with the least load so far (ties:
    path, then lowest shard index), which is deterministic for a given
    universe and durations file. A file the durations file does not know is
    weighted with the median of the known files and assigned the same way:
    a new file is always placed, and only moves the plan, which every shard
    recomputes identically from the same commit.

    Each shard collects the whole suite (shard_plugin.py), runs its own
    files and writes a report directory: meta.json, universe.txt,
    selected.txt, ran.txt, durations.json. ``verify`` passes only if there
    are exactly N reports for shards 1..N, all for the expected commit, all
    with the same universe and the same plan, the selected sets are
    disjoint by node id and by file and their union is the universe, every
    selected test finished and nothing else did, and pytest exited 0 in
    every shard. A missing report, file or field is a failure.

Rejected alternatives
    A planning job that collects once and hands file lists to the shards:
    one more serial job (install + collection, about two minutes on the
    runner) in front of every shard. ``--collect-only`` per shard followed
    by a run over a file list: two collections per shard, and the universe
    and the run would come from different processes. pytest-split and
    pytest-xdist: third-party plugins outside the dependency set, and xdist
    in one process tree reports shared-state failures a serial run does
    not. Splitting by node id: a module fixture would run in several
    shards and a file would be split across processes. Hash-only
    assignment: no balancing; over the measured durations with 4 shards
    its slowest shard is 767 s against 702 s here (+9 %).

Commands
    plan --collected FILE [--shards N] [--list K]
        expected seconds per shard for a ``pytest --collect-only -q``
        listing; --list prints the files of shard K.
    verify ROOT --shards N [--sha SHA]
        the aggregator check over downloaded shard reports; exit 1 on any
        problem, and a summary appended to $GITHUB_STEP_SUMMARY when set.
    durations ROOT
        a new test-durations.json from the reports of a green run:
            gh run download <run-id> -p 'tests-shard-*' -D /tmp/shards
            python3 .github/ci/ci_shards.py durations /tmp/shards \\
                > .github/ci/test-durations.json

Standard library only: the aggregator runs it on the runner's python3
without installing anything.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import statistics
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
DURATIONS = HERE / "test-durations.json"

REPORT_FILES = ("meta.json", "universe.txt", "selected.txt", "ran.txt")


def file_of(nodeid: str) -> str:
    return nodeid.split("::", 1)[0]


def load_durations(path: Path | str = DURATIONS) -> dict[str, float]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    files = data["files"]
    return {str(k): float(v) for k, v in files.items()}


def plan(files, durations: dict[str, float], shards: int) -> dict[str, int]:
    """Assign every file to one shard in 1..shards (longest first)."""
    if shards < 1:
        raise ValueError(f"shards must be >= 1, got {shards}")
    files = sorted(set(files))
    known = [durations[f] for f in files if f in durations]
    default = statistics.median(known) if known else 1.0
    weight = {f: durations.get(f, default) for f in files}
    load = [0.0] * shards
    out: dict[str, int] = {}
    for f in sorted(files, key=lambda f: (-weight[f], f)):
        k = min(range(shards), key=lambda i: (load[i], i))
        out[f] = k + 1
        load[k] += weight[f]
    return out


def plan_digest(assignment: dict[str, int]) -> str:
    text = "\n".join(f"{f}\t{k}" for f, k in sorted(assignment.items()))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _lines(path: Path) -> list[str]:
    return [l for l in path.read_text(encoding="utf-8").splitlines() if l]


def _read_report(d: Path, problems: list[str]) -> dict | None:
    missing = [n for n in REPORT_FILES if not (d / n).is_file()]
    if missing:
        problems.append(f"{d.name}: missing {', '.join(missing)}")
        return None
    try:
        meta = json.loads((d / "meta.json").read_text(encoding="utf-8"))
        rep = {
            "dir": d.name,
            "shard": int(meta["shard"]),
            "shards": int(meta["shards"]),
            "sha": str(meta["sha"]),
            "exitstatus": int(meta["exitstatus"]),
            "plan": str(meta["plan_digest"]),
            "universe": _lines(d / "universe.txt"),
            "selected": _lines(d / "selected.txt"),
            "ran": _lines(d / "ran.txt"),
        }
    except (KeyError, TypeError, ValueError) as exc:
        problems.append(f"{d.name}: unreadable report ({exc!r})")
        return None
    seconds = 0.0
    if (d / "durations.json").is_file():
        try:
            seconds = sum(json.loads(
                (d / "durations.json").read_text(encoding="utf-8")).values())
        except (ValueError, AttributeError, TypeError):
            pass
    rep["seconds"] = seconds
    return rep


def verify(root: Path | str, shards: int, sha: str | None = None):
    """Return (problems, summary lines) for the reports under root."""
    root = Path(root)
    problems: list[str] = []
    if not root.is_dir():
        return [f"{root}: no reports were downloaded"], []
    dirs = sorted(p for p in root.iterdir() if p.is_dir())
    reports = [r for r in (_read_report(d, problems) for d in dirs) if r]
    if len(dirs) != shards:
        problems.append(f"expected {shards} shard reports, found {len(dirs)}")

    seen = Counter(r["shard"] for r in reports)
    for k, n in sorted(seen.items()):
        if n > 1:
            problems.append(f"shard {k} reported {n} times")
        if not 1 <= k <= shards:
            problems.append(f"shard {k} is outside 1..{shards}")
    for k in range(1, shards + 1):
        if k not in seen:
            problems.append(f"shard {k} has no report")

    for r in reports:
        if r["shards"] != shards:
            problems.append(f"{r['dir']}: split {r['shards']} ways, "
                            f"expected {shards}")
        if sha is not None and r["sha"] != sha:
            problems.append(f"{r['dir']}: commit {r['sha']}, expected {sha}")
        if r["exitstatus"] != 0:
            problems.append(f"{r['dir']}: pytest exit status {r['exitstatus']}")
        dup = [n for n, c in Counter(r["ran"]).items() if c > 1]
        if dup:
            problems.append(f"{r['dir']}: {len(dup)} test(s) ran twice, "
                            f"e.g. {dup[0]}")
        sel, ran = set(r["selected"]), set(r["ran"])
        if sel - ran:
            problems.append(f"{r['dir']}: {len(sel - ran)} selected test(s) "
                            f"did not finish, e.g. {sorted(sel - ran)[0]}")
        if ran - sel:
            problems.append(f"{r['dir']}: {len(ran - sel)} test(s) ran that "
                            f"were not selected, e.g. {sorted(ran - sel)[0]}")

    if not reports:
        return problems or ["no readable shard report"], []

    universe = set(reports[0]["universe"])
    if not universe:
        problems.append(f"{reports[0]['dir']}: empty universe")
    for r in reports[1:]:
        if set(r["universe"]) != universe:
            diff = set(r["universe"]) ^ universe
            problems.append(f"{r['dir']}: collected a different suite than "
                            f"{reports[0]['dir']} ({len(diff)} node ids "
                            f"differ, e.g. {sorted(diff)[0]})")
    if len({r["plan"] for r in reports}) != 1:
        problems.append("the shards computed different plans")

    owner: dict[str, str] = {}
    file_owner: dict[str, str] = {}
    for r in reports:
        for nodeid in r["selected"]:
            if nodeid in owner and owner[nodeid] != r["dir"]:
                problems.append(f"{nodeid} selected by {owner[nodeid]} "
                                f"and {r['dir']}")
            owner.setdefault(nodeid, r["dir"])
        for f in {file_of(n) for n in r["selected"]}:
            if f in file_owner and file_owner[f] != r["dir"]:
                problems.append(f"{f} split between {file_owner[f]} "
                                f"and {r['dir']}")
            file_owner.setdefault(f, r["dir"])
    union = set(owner)
    if universe - union:
        problems.append(f"{len(universe - union)} collected test(s) ran in "
                        f"no shard, e.g. {sorted(universe - union)[0]}")
    if union - universe:
        problems.append(f"{len(union - universe)} selected test(s) are not "
                        f"in the collected suite, e.g. "
                        f"{sorted(union - universe)[0]}")

    summary = ["| shard | files | tests | finished | exit | seconds |",
               "|---|---|---|---|---|---|"]
    for r in sorted(reports, key=lambda r: r["shard"]):
        files = len({file_of(n) for n in r["selected"]})
        summary.append(f"| {r['shard']}/{r['shards']} | {files} | "
                       f"{len(r['selected'])} | {len(set(r['ran']))} | "
                       f"{r['exitstatus']} | {r['seconds']:.0f} |")
    summary.append(f"\nuniverse: {len(universe)} tests in "
                   f"{len({file_of(n) for n in universe})} files; "
                   f"union of shards: {len(union)}")
    return problems, summary


def merged_durations(root: Path | str) -> dict[str, float]:
    out: dict[str, float] = {}
    for d in sorted(Path(root).iterdir()):
        p = d / "durations.json"
        if p.is_file():
            for f, s in json.loads(p.read_text(encoding="utf-8")).items():
                out[f] = out.get(f, 0.0) + float(s)
    return out


def _collected(path: str) -> list[str]:
    return [l.strip() for l in Path(path).read_text(encoding="utf-8")
            .splitlines() if "::" in l]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("plan")
    p.add_argument("--collected", required=True)
    p.add_argument("--durations", default=str(DURATIONS))
    p.add_argument("--shards", type=int, default=4)
    p.add_argument("--list", type=int, default=None)
    v = sub.add_parser("verify")
    v.add_argument("root")
    v.add_argument("--shards", type=int, required=True)
    v.add_argument("--sha", default=None)
    d = sub.add_parser("durations")
    d.add_argument("root")
    a = ap.parse_args(argv)

    if a.cmd == "plan":
        ids = _collected(a.collected)
        durations = load_durations(a.durations)
        files = sorted({file_of(n) for n in ids})
        assignment = plan(files, durations, a.shards)
        if a.list is not None:
            print("\n".join(f for f in files if assignment[f] == a.list))
            return 0
        tests = Counter(assignment[file_of(n)] for n in ids)
        known = [durations[f] for f in files if f in durations]
        default = statistics.median(known) if known else 1.0
        print(f"{len(ids)} tests in {len(files)} files, {a.shards} shards; "
              f"{len(files) - len(known)} files without a measured duration "
              f"(weighted {default:.1f} s)")
        print("shard  files  tests  expected_s")
        for k in range(1, a.shards + 1):
            fs = [f for f in files if assignment[f] == k]
            secs = sum(durations.get(f, default) for f in fs)
            print(f"{k:5d}  {len(fs):5d}  {tests[k]:5d}  {secs:10.1f}")
        union = {f for f in files if 1 <= assignment[f] <= a.shards}
        print(f"files assigned: {len(assignment)}, union equals the "
              f"collected files: {union == set(files)}")
        return 0

    if a.cmd == "verify":
        problems, summary = verify(a.root, a.shards, a.sha)
        text = "\n".join(summary)
        if problems:
            text += "\n\n**Incomplete:**\n" + "\n".join(
                f"- {p}" for p in problems[:50])
        print(text)
        step_summary = os.environ.get("GITHUB_STEP_SUMMARY")
        if step_summary:
            with open(step_summary, "a", encoding="utf-8") as fh:
                fh.write("### Fast suite shards\n\n" + text + "\n")
        return 1 if problems else 0

    if a.cmd == "durations":
        merged = merged_durations(a.root)
        if not merged:
            print("no durations.json under " + a.root, file=sys.stderr)
            return 1
        json.dump({
            "about": "Seconds per test file of the fast suite (pytest -m "
                     "'not slow') on a GitHub-hosted ubuntu-latest runner. "
                     "Input to .github/ci/ci_shards.py; regenerate with its "
                     "'durations' command.",
            "files": {f: round(s, 2) for f, s in sorted(merged.items())},
        }, sys.stdout, indent=1)
        print()
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
