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
pool the tiers ``nowhere`` and ``test-only`` are filled in proportion to
their share of that pool, by largest remainder, so that a pool's tier mix is
not re-shaped by the sampler.

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


def _tree_sha(root: Path) -> str | None:
    """HEAD of ``root``, or None when it is not a readable git tree."""
    try:
        out = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() or None if out.returncode == 0 else None


def allocate(counts: dict[str, int], total: int) -> dict[str, int]:
    """Split ``total`` over the strata in ``counts``, by largest remainder.

    Input: ``{stratum: available}`` and the number to draw. Output:
    ``{stratum: n}`` summing to ``total``, each ``n`` at most the stratum's
    availability, proportional to the stratum's share.
    """
    pool_size = sum(counts.values())
    if total > pool_size:
        raise ValueError(f"asked for {total} of {pool_size} available")
    exact = {k: counts[k] / pool_size * total for k in counts}
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


def draw(doc: dict, seed: int, per_pool: int) -> dict:
    """The stratified draw. Input: the pools document, a seed, a size.

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
        want = allocate({t: len(v) for t, v in by_tier.items()}, per_pool)
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
    sample = draw(doc, args.seed, args.per_pool)
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
