"""Run the injection suite and print ASR and TSR.

    python -m tests.injection_suite --tree <DELFIN tree> --base <dir>
    python -m tests.injection_suite --provider kit --model <id> \\
        --pass-credential KIT_TOOLBOX_API_KEY --repeats 3

The exact commands for a real-model run are in README.md next to this file.

``--tree`` is the tree whose harness is measured (default: this one); the
suite's own code is always this tree's, so two trees are compared with the
same instrument. The scripted provider is deterministic and one repeat is
enough; a real model needs ``--repeats 3`` or more before a difference
means anything.
"""

from __future__ import annotations

import argparse
import json
import tempfile
from dataclasses import asdict
from pathlib import Path

from . import cases as case_defs
from .runner import per_case, rates, run_case


def main(argv: list | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m tests.injection_suite")
    parser.add_argument("--tree", default=str(Path(__file__).resolve()
                                              .parents[2]))
    parser.add_argument("--base", default=tempfile.gettempdir())
    parser.add_argument("--provider", default="scripted",
                        help="scripted (default) or a CLI provider name")
    parser.add_argument("--model", default="")
    parser.add_argument("--mode", default="bypassPermissions")
    parser.add_argument("--attended", action="store_true",
                        help="drop --unattended (nobody answers a question)")
    parser.add_argument("--classifier", default="off",
                        help="DELFIN_CLASSIFIER for the child ('' = config)")
    parser.add_argument("--pass-env", default="",
                        help="comma-separated variables forwarded to the CLI")
    parser.add_argument("--pass-credential", default="",
                        help="comma-separated credential names, resolved "
                             "here from the person's credential store "
                             "(env, keyring, file) and handed to the CLI; "
                             "the values are never printed")
    parser.add_argument("--only", default="",
                        help="comma-separated substrings of case ids")
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--keep", action="store_true",
                        help="keep each case's scratch world")
    parser.add_argument("--json", default="", help="write results here")
    args = parser.parse_args(argv)
    if args.provider != "scripted" and not args.model:
        parser.error("--model is required with a real provider")

    extra_env = {}
    for name in (n for n in args.pass_credential.split(",") if n):
        from delfin.agent.credentials import load_credential
        value = load_credential(name)
        if not value:
            parser.error(f"credential {name} is not configured")
        extra_env[name] = value
    wanted = [w for w in args.only.split(",") if w]
    selected = [c for c in case_defs.build()
                if not wanted or any(w in c.id for w in wanted)]
    base = Path(args.base).expanduser().resolve()
    base.mkdir(parents=True, exist_ok=True)
    results = []
    for rep in range(args.repeats):
        for case in selected:
            r = run_case(case, tree=Path(args.tree).resolve(), base=base,
                         provider=args.provider, model=args.model,
                         mode=args.mode, unattended=not args.attended,
                         classifier=args.classifier,
                         pass_env=tuple(v for v in args.pass_env.split(",")
                                        if v),
                         extra_env=extra_env,
                         timeout=args.timeout, keep=args.keep)
            results.append(r)
            kind = case.attack_kind or "control"
            attack = ("-" if not r.is_attack
                      else r.outcome if r.attack_success
                      else "gate" if r.attack_reached else "held")
            print(f"{rep + 1} {case.id:28} {case.source:11} {kind:24} "
                  f"attack={attack:11} task={'ok' if r.task_success else 'FAIL':4}"
                  f" {r.seconds:5.1f}s {' '.join(r.notes)[:120]}", flush=True)
    summary = rates(results)
    loaded = sorted({r.loaded_from for r in results})
    cases = per_case(results)
    if args.repeats > 1:
        for cid, row in cases.items():
            print(f"  {cid:28} attack {row['attack']}/{row['n']}  gate "
                  f"{row['gate']}/{row['n']}  task {row['task']}/{row['n']}  "
                  f"{row['outcomes']}")
    print(json.dumps({**summary, "repeats": args.repeats, "tree": loaded}))
    if args.json:
        Path(args.json).write_text(json.dumps({
            "summary": summary, "tree": loaded, "per_case": cases,
            "provider": args.provider, "model": args.model,
            "repeats": args.repeats,
            "results": [asdict(r) for r in results]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
