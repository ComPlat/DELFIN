#!/usr/bin/env python3
"""Draw the frozen task list for the specialist-versus-generalist comparison.

Input: the JSON document written by ``ab_task_pools.py``, a seed and a
per-pool sample size. Output: one JSON file holding the drawn tasks, the
seed, the pool and tier counts, the draw date, the git sha of the tree the
pools were scanned from, and the sha256 of the pools document.

Semantics
---------
The draw is **stratified twice**: the same number of tasks from ``infra`` and
from ``chem``, so that arm and domain are not confounded; and within each
pool the tiers ``nowhere`` and ``test-only`` are filled **equally**
(``balanced``, the default), so that tier varies inside each arm instead of
between the pools and can be reported as a covariate.

Why not in proportion to each pool's tier mix, which was the first rule:
the two pools' mixes are near-opposites -- measured 2026-11-20, infra is
34 % ``nowhere`` (27 of 80) and chem 69 % (54 of 78), and 71 % (70 of 99)
once the chemistry rule gained its S3 clause -- so a proportional draw
hands one arm
mostly ``test-only`` tasks and the other mostly ``nowhere`` ones. Tier
plausibly predicts difficulty (a ``test-only`` capability has a test
pinning its behaviour; a ``nowhere`` one is more often simply dead), so
that draw confounds domain with tier, which is the confound the two pools
exist to remove. The proportional rule is kept as ``--allocation
proportional`` so that a superseded draw can still be reproduced from its
seed, and the rule used is recorded in the output as ``allocation_rule``.

The draw is **frozen**: the output file is written once. A second run refuses
and prints the sha256 of the file that already exists, because a silent
re-draw after seeing results is the shape that turns a sample into a choice.
``--force`` overwrites, and the overwrite is recorded in the new file's
``replaced`` field together with the sha256 it replaced, so that a re-draw
leaves a trace in the artefact rather than only in a shell history.

The draw refuses outright when a pool holds fewer tasks than requested: a
short pool is reported by ``ab_task_pools.py`` before anything runs, and
silently drawing fewer would make the two pools unequal and the comparison
unbalanced.

Nothing here reads an arm, a role or a model. The arms are attached by the
runner, after the draw.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import random
import subprocess
from collections import Counter
from pathlib import Path

POOLS = ("infra", "chem")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


#: The prefixes a committed artefact's path provenance is rewritten to:
#: an absolute scratch or home path encodes a home directory and the site
#: group-account name, which no file in a public repository carries. The
#: rewrite is a redaction, so the path is no longer resolvable.
_REDACTED_PREFIXES = ("<", "~")


def is_redacted(path: str | None) -> bool:
    """True when ``path`` is a redacted provenance string, not a tree.

    Input: the ``generated_for`` field of a pools document. A draw reads
    that path to compare the tree's HEAD against the sha the scan ran on;
    a redacted value cannot be read, and a draw that proceeded would
    record ``tree_sha_at_draw: null`` -- an unobserved check reported as
    a measured one.
    """
    return bool(path) and str(path).startswith(_REDACTED_PREFIXES)


def _tree_sha(root: Path) -> str | None:
    """HEAD of ``root``, or None when it is not a readable git tree."""
    try:
        out = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() or None if out.returncode == 0 else None


def allocate(counts: dict[str, int], total: int,
             rule: str = "balanced") -> dict[str, int]:
    """Split ``total`` over the strata in ``counts``.

    Input: ``{stratum: available}``, the number to draw, and the rule --
    ``balanced`` (each stratum the same share) or ``proportional`` (each
    stratum its share of the pool). Output: ``{stratum: n}`` summing to
    ``total``, each ``n`` at most that stratum's availability. The
    remainder goes by largest remainder in both cases, so the split does
    not depend on the order the strata arrive in.

    A balanced draw refuses when a stratum holds fewer than ``total //
    len(counts)``, the floor rather than the exact share: with two tiers
    and an even total the two agree, but with three or more strata, or an
    odd total, the exact share is not an integer and comparing against it
    refused a split that exists -- {3, 4, 5} for a total of 10 is
    feasible as 4/3/3 and was rejected because 3 < 3.33. The refusal is
    there to stop a pool silently supplying 2 where 5 was asked for, not
    to reject every split that cannot be exactly equal (Addendum 2,
    Correction F).
    """
    pool_size = sum(counts.values())
    if total > pool_size:
        raise ValueError(f"asked for {total} of {pool_size} available")
    if rule == "balanced":
        exact = {k: total / len(counts) for k in counts}
        floor_each = total // len(counts)
        short = {k: v for k, v in counts.items() if v < floor_each}
        if short:
            raise ValueError(
                f"a balanced draw needs {floor_each} of each "
                f"stratum; these hold fewer: "
                f"{ {k: counts[k] for k in short} }")
    elif rule == "proportional":
        exact = {k: counts[k] / pool_size * total for k in counts}
    else:
        raise ValueError(f"unknown allocation rule {rule!r}")
    out = {k: int(exact[k]) for k in counts}
    remainder = total - sum(out.values())
    order = sorted(counts, key=lambda k: (-(exact[k] - out[k]), k))
    i = 0
    while remainder > 0:
        k = order[i % len(order)]
        if out[k] < counts[k]:
            out[k] += 1
            remainder -= 1
        i += 1
        if i > 4 * len(order) + total:          # cannot happen; no spin
            raise ValueError("allocation did not converge")
    return out


def draw(doc: dict, seed: int, per_pool: int,
         rule: str = "balanced") -> dict:
    """The stratified draw. Input: the pools document, a seed, a size, a rule.

    Output: the frozen-sample document, without the ``replaced`` field,
    which only an overwrite adds.
    """
    rng = random.Random(seed)
    tasks: list[dict] = []
    allocation: dict[str, dict[str, int]] = {}
    for pool in POOLS:
        rows = sorted(doc["pools"][pool], key=lambda r: (r["name"], r["file"]))
        if len(rows) < per_pool:
            raise ValueError(
                f"pool {pool} holds {len(rows)} candidates and "
                f"{per_pool} were asked for; an unequal draw would "
                "confound arm with domain")
        by_tier: dict[str, list[dict]] = {}
        for row in rows:
            by_tier.setdefault(row["tier"], []).append(row)
        want = allocate({t: len(v) for t, v in by_tier.items()},
                        per_pool, rule)
        allocation[pool] = want
        for tier in sorted(want):
            for row in rng.sample(by_tier[tier], want[tier]):
                tasks.append({
                    "pool": pool,
                    "tier": tier,
                    "name": row["name"],
                    "file": row["file"],
                    "line": row["line"],
                    "module": row["module"],
                    "task": doc["task_sentence"],
                })
    tasks.sort(key=lambda t: (t["pool"], t["tier"], t["name"]))
    return {
        "drawn_at": dt.datetime.now(dt.timezone.utc).isoformat(
            timespec="seconds"),
        "seed": seed,
        "per_pool": per_pool,
        "allocation_rule": rule,
        "n_tasks": len(tasks),
        "pools_tree": doc.get("generated_for"),
        "pools_tree_sha": None,          # filled by main()
        "pools_file_sha256": None,       # filled by main()
        "pool_counts": doc["counts"],
        "pool_counts_by_tier": doc["counts_by_tier"],
        "underpowered": doc["underpowered"],
        "allocation": allocation,
        "drawn_counts": {
            p: dict(Counter(t["tier"] for t in tasks if t["pool"] == p))
            for p in POOLS},
        "task_sentence": doc["task_sentence"],
        "tasks": tasks,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pools", required=True, help="ab_task_pools.py JSON")
    ap.add_argument("--out", required=True, help="the frozen sample file")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--per-pool", type=int, default=10,
                    help="tasks drawn from each pool (total is twice this)")
    ap.add_argument("--allocation", choices=("balanced", "proportional"),
                    default="balanced",
                    help="tiers equally within each pool, or in proportion "
                         "to the pool's mix (superseded 2026-11-20)")
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing sample, recording its sha256")
    args = ap.parse_args(argv)

    out = Path(args.out)
    replaced = None
    if out.exists():
        existing = _sha256(out)
        if not args.force:
            print(f"refusing to re-draw: {out} exists, sha256 {existing}")
            print("pass --force to replace it; the replacement records this "
                  "sha256")
            return 2
        replaced = {"sha256": existing,
                    "at": dt.datetime.now(dt.timezone.utc).isoformat(
                        timespec="seconds")}

    pools_path = Path(args.pools)
    doc = json.loads(pools_path.read_text())
    if is_redacted(doc.get("generated_for")):
        print(f"refusing to draw: the pools document records a redacted "
              f"path ({doc['generated_for']}). A committed artefact carries "
              "no host path, so the tree it was scanned from cannot be read "
              "from here and 'the tree has moved since the scan' would be "
              "recorded as unknown rather than checked. Re-run "
              "ab_task_pools.py against the tree to draw from it.")
        return 2
    sample = draw(doc, args.seed, args.per_pool, args.allocation)
    sample["pools_file"] = str(pools_path)
    sample["pools_file_sha256"] = _sha256(pools_path)
    tree = doc.get("generated_for")
    # The sha the scan ran on, from the document; the live HEAD beside it,
    # because a tree that has moved since the scan is a fact about the draw.
    sample["pools_tree_sha"] = doc.get("tree_sha")
    sample["tree_sha_at_draw"] = _tree_sha(Path(tree)) if tree else None
    sample["pools_scanned_at"] = doc.get("scanned_at")
    if replaced:
        sample["replaced"] = replaced

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(sample, indent=2, sort_keys=True))
    print(f"wrote {out}  sha256 {_sha256(out)}")
    print(f"seed {sample['seed']}  tasks {sample['n_tasks']}  "
          f"allocation {sample['allocation_rule']}  "
          f"tree {sample['pools_tree_sha']}")
    if (sample["tree_sha_at_draw"]
            and sample["tree_sha_at_draw"] != sample["pools_tree_sha"]):
        print(f"note: the tree is now at {sample['tree_sha_at_draw']}, not "
              "the sha the pools were scanned on")
    for pool in POOLS:
        print(f"    {pool}: drawn {sample['drawn_counts'][pool]} "
              f"of available {sample['pool_counts_by_tier'][pool]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
