"""Run the probes and print a table.

    python -m tests.agent_terminal_probes --repeats 3

Every repeat must hold. Unlike the benchmark's median scoring, a safety
invariant that holds two times in three does not hold — it just has not
been unlucky yet. Anything less than all-green is a red result.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

from . import probes as probe_defs
from . import sandbox as sb
from .runner import run_headless, run_probe


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m tests.agent_terminal_probes")
    parser.add_argument("--model", default="kit.qwen3.5-397b-A17b")
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--only", default="",
                        help="comma-separated substrings")
    parser.add_argument("--family", default="safety",
                        choices=["safety", "regressions", "both"],
                        help="safety: can anything get past a mechanism. "
                             "regressions: does the agent still do the job "
                             "it was fixed to do.")
    parser.add_argument("--timeout", type=int, default=240)
    args = parser.parse_args(argv)

    from delfin.agent.credentials import load_credential
    api_key = load_credential("KIT_TOOLBOX_API_KEY")
    if not api_key:
        print("KIT_TOOLBOX_API_KEY is not configured", file=sys.stderr)
        return 2

    failures: list[str] = []
    for run in range(1, args.repeats + 1):
        with tempfile.TemporaryDirectory(prefix="delfin-probe-") as tmp:
            box = sb.build(Path(tmp))
            port = sb.closed_port()
            wanted = [w for w in args.only.split(",") if w]
            pool = []
            if args.family in ("safety", "both"):
                pool += probe_defs.build(box, port=port)
            if args.family in ("regressions", "both"):
                from . import regressions as regression_defs
                pool += regression_defs.build(box)
            selected = [p for p in pool
                        if not wanted or any(w in p.name for w in wanted)]
            print(f"\n=== run {run}/{args.repeats} · {len(selected)} probes ===")
            for probe in selected:
                drive = run_headless if probe.headless else run_probe
                out = drive(probe, box, api_key=api_key, model=args.model,
                            timeout=args.timeout)
                mark = "ok  " if out.passed else "FAIL"
                asked = " (asked)" if out.prompted else ""
                print(f"  {mark}  {probe.name:26} {out.detail}{asked}")
                if not out.passed:
                    failures.append(f"run {run}: {probe.name} — {out.detail}")

    print()
    if failures:
        print(f"{len(failures)} FAILED")
        for line in failures:
            print(f"  {line}")
        return 1
    print("all probes held, every repeat")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
