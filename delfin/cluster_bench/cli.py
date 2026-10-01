"""``delfin cluster`` -- batch construction runs on a Slurm cluster (see the package doc).

    delfin cluster prepare  --tool T --input LIST --run-dir RUN [--select IDS] [--tool-python PY]
    delfin cluster slurm    RUN [--set main|repeat] [--run main|pilot|...] [--throttle 40]
                            [--time 72:00:00] [--speed-factor F] [--array 0-1] [--submit]
    delfin cluster run-shard RUN --shard K [--set S] [--run R]
    delfin cluster status   RUN [--set S] [--run R] [--json] [-v]
    delfin cluster collect  RUN [--set S] [--run R] [--dest DIR]
    delfin cluster repeat-stats RUN
"""
from __future__ import annotations

import argparse
import json
import sys

from delfin.cluster_bench.provenance import MODES, TOOLS


def cbatch_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="delfin cluster", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("prepare", help="shards + specs + manifest in a NEW run directory")
    p.add_argument("--tool", required=True, choices=sorted(TOOLS))
    p.add_argument("--input", required=True, help="ID;SMILES (or ID|SMILES) per line")
    p.add_argument("--run-dir", required=True, help="new directory (never overwritten)")
    p.add_argument("--select", default=None,
                   help="IDs to build, one per line, in this order (default: the whole input)")
    p.add_argument("--specs", default=None,
                   help="precomputed specs (JSONL, one per ID) instead of computing them")
    p.add_argument("--tool-python", default=None,
                   help="interpreter of the tool's own environment (external builders)")
    p.add_argument("--mode", default=None,
                   help="; ".join(f"{t}: {'/'.join(m)}" for t, m in MODES.items()))
    p.add_argument("--shard-size", type=int, default=None, help="manta 500, others 250")
    p.add_argument("--salt", default="delfin-cluster-v1", help="shard order / repeat choice")
    p.add_argument("--timeout", type=int, default=21600, help="per-system limit before scaling, s")
    p.add_argument("--speed-factor", default="1.0", help="limit = ceil(timeout x factor)")
    p.add_argument("--workers", type=int, default=None, help="concurrent builds per node")
    p.add_argument("--threads", type=int, default=None, help="threads per build (manta 7)")
    p.add_argument("--repeat", type=int, default=0,
                   help="also prepare a set of N systems for a second build (byte identity)")
    p.add_argument("--spec-workers", type=int, default=8)
    p.add_argument("--label", default=None, help="archive label (default <tool>_<run dir name>)")

    p = sub.add_parser("slurm", help="write (and optionally submit) the sbatch array script")
    p.add_argument("run_dir")
    p.add_argument("--set", default="main", dest="set_name")
    p.add_argument("--run", default="main", dest="run_name",
                   help="output subdirectory; use e.g. 'pilot' for a pilot with another factor")
    p.add_argument("--array", default=None, help="indices, default all shards (e.g. 0-1 for a pilot)")
    p.add_argument("--throttle", type=int, default=40, help="array tasks at once (%%N), 0 = none")
    p.add_argument("--time", default="72:00:00", dest="time_limit")
    p.add_argument("--cpus", type=int, default=48)
    p.add_argument("--mem", default=None, help="manta 72G, others 80G")
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--speed-factor", default=None)
    p.add_argument("--partition", default=None)
    p.add_argument("--account", default=None)
    p.add_argument("--python", default=None, help="DELFIN interpreter on the nodes (default: this one)")
    p.add_argument("--setup", action="append", default=[],
                   help="extra shell line before the run, e.g. 'module load ...' (repeatable)")
    p.add_argument("--submit", action="store_true", help="call sbatch")

    p = sub.add_parser("run-shard", help="build one shard (the array task)")
    p.add_argument("run_dir")
    p.add_argument("--shard", required=True, type=int)
    p.add_argument("--set", default="main", dest="set_name")
    p.add_argument("--run", default="main", dest="run_name")
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--speed-factor", default=None)
    p.add_argument("--allow-env-change", action="store_true",
                   help="build even if code or environment differ from the manifest (recorded)")

    for name, hlp in (("status", "progress per shard"), ("collect", "merge into one archive")):
        p = sub.add_parser(name, help=hlp)
        p.add_argument("run_dir")
        p.add_argument("--set", default="main", dest="set_name")
        p.add_argument("--run", default="main", dest="run_name")
        if name == "status":
            p.add_argument("--json", action="store_true")
            p.add_argument("-v", "--verbose", action="store_true")
        else:
            p.add_argument("--dest", default=None, help="NEW directory (default RUN/collected/<label>)")
            p.add_argument("--label", default=None)

    p = sub.add_parser("repeat-stats", help="byte identity of the repeat set against its first build")
    p.add_argument("run_dir")
    p.add_argument("--run-a", default="main", help="run of the main set")
    p.add_argument("--run-b", default="main", help="run of the repeat set")
    return ap


def main(argv=None) -> int:
    a = cbatch_parser().parse_args(argv)
    if a.cmd == "prepare":
        from delfin.cluster_bench.prepare import cbatch_prepare

        man = cbatch_prepare(tool=a.tool, input_list=a.input, run_dir=a.run_dir, selection=a.select,
                             specs_file=a.specs, size=a.shard_size, salt=a.salt, mode=a.mode,
                             tool_python=a.tool_python, timeout_base=a.timeout,
                             speed_factor=a.speed_factor, workers=a.workers, threads=a.threads,
                             repeat=a.repeat, spec_workers=a.spec_workers, label=a.label)
        print(json.dumps({"run_dir": a.run_dir, "tool": man["tool"], "label": man["label"],
                          "n_systems": man["n_systems"],
                          "sets": {k: v["n_shards"] for k, v in man["sets"].items()},
                          "timeout_s": man["settings"]["timeout_s"],
                          "commit": man["provenance"]["code"]["commit"],
                          "dirty": man["provenance"]["code"]["dirty"]}, indent=1))
        return 0
    if a.cmd == "slurm":
        from delfin.cluster_bench.slurm_script import cbatch_write_sbatch

        path, cmd, out = cbatch_write_sbatch(
            a.run_dir, submit=a.submit, set_name=a.set_name, run_name=a.run_name, array=a.array,
            throttle=a.throttle, time_limit=a.time_limit, cpus=a.cpus, mem=a.mem,
            workers=a.workers, speed_factor=a.speed_factor, partition=a.partition,
            account=a.account, python=a.python, setup=a.setup)
        print(f"script: {path}")
        print(out if out is not None else "submit with:  " + " ".join(cmd))
        return 0
    if a.cmd == "run-shard":
        from delfin.cluster_bench.runner import cbatch_run_shard

        return cbatch_run_shard(a.run_dir, a.shard, set_name=a.set_name, run_name=a.run_name,
                                workers=a.workers, speed_factor=a.speed_factor,
                                allow_env_change=a.allow_env_change,
                                log=lambda m: print(m, flush=True))
    if a.cmd == "status":
        from delfin.cluster_bench.report import cbatch_format_status, cbatch_status

        st = cbatch_status(a.run_dir, a.set_name, a.run_name)
        print(json.dumps(st, indent=1) if a.json else cbatch_format_status(st, a.verbose))
        return 0
    if a.cmd == "collect":
        from delfin.cluster_bench.report import cbatch_collect

        s = cbatch_collect(a.run_dir, a.set_name, a.run_name, dest=a.dest, label=a.label)
        print(json.dumps({k: s[k] for k in ("label", "n_systems_merged", "n_systems_in_set",
                                            "n_expressible", "by_class", "coverage_of_expressible",
                                            "n_xyz", "chunks_merged", "resubmit_array",
                                            "timeout", "n_problems")}, indent=1))
        return 0 if s["n_problems"] == 0 else 1
    if a.cmd == "repeat-stats":
        from delfin.cluster_bench.report import cbatch_repeat_stats

        r = cbatch_repeat_stats(a.run_dir, a.run_a, a.run_b)
        print(json.dumps({k: v for k, v in r.items() if k not in ("ids", "different_detail")}, indent=1))
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
