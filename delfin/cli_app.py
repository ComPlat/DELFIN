"""``delfin-app`` — generate and run DELFIN applications from a keys file.

A CONTROL.txt-style workflow: generate an editable keys file for an application,
adjust the values, and run it.

    delfin-app list                                  # registered applications
    delfin-app template redox_potential > my.txt     # generate a keys file
    #   ... edit my.txt ...
    delfin-app run my.txt                            # run it (synchronous)
    delfin-app run my.txt --submit                   # run in the background
    delfin-app describe redox_potential              # full contract as JSON
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="delfin-app",
        description="Generate and run DELFIN applications from a keys file.",
    )
    sub = parser.add_subparsers(dest="cmd")
    sub.add_parser("list", help="List registered applications")
    t = sub.add_parser("template", help="Print an editable keys file for an application")
    t.add_argument("name")
    t.add_argument("-o", "--output", help="Write to a file instead of stdout")
    r = sub.add_parser("run", help="Run an application from a keys file")
    r.add_argument("keyfile")
    r.add_argument("--cores", type=int, default=1)
    r.add_argument("--submit", action="store_true", help="Run in the background")
    d = sub.add_parser("describe", help="Print an application's full contract (JSON)")
    d.add_argument("name")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    from delfin.tools import platform

    if args.cmd == "list":
        apps = platform.list_applications()
        if not apps:
            print("No applications registered.")
            return 0
        for name in apps:
            app = platform.describe_application(name)
            print(f"{name}  —  {app.description}")
        return 0

    if args.cmd == "template":
        from delfin.tools._keyfile import application_keyfile
        try:
            text = application_keyfile(args.name)
        except ValueError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 1
        if args.output:
            Path(args.output).write_text(text, encoding="utf-8")
            print(f"wrote {args.output}")
        else:
            sys.stdout.write(text)
        return 0

    if args.cmd == "describe":
        app = platform.describe_application(args.name)
        if app is None:
            print(f"unknown application {args.name!r}", file=sys.stderr)
            return 1
        print(json.dumps(app.to_dict(), indent=2, default=str))
        return 0

    if args.cmd == "run":
        from delfin.tools._keyfile import run_from_keyfile
        try:
            result = run_from_keyfile(args.keyfile, cores=args.cores, submit=args.submit)
        except (ValueError, FileNotFoundError, OSError) as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 1
        if args.submit:
            print(f"submitted: run {result}")
            return 0
        if result.ok:
            print("OK")
            for key, value in (result.outputs or {}).items():
                print(f"  {key} = {value}")
            return 0
        print(f"FAILED: {result.error or 'see validation report'}", file=sys.stderr)
        return 1

    parser.print_help()
    return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
