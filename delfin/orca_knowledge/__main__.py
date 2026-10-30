"""``python -m delfin.orca_knowledge refresh-keywords``: rewrite the keyword
catalog from the installed OPI (it must be the pinned version), and
``python -m delfin.orca_knowledge check``: the problems of the card store."""

from __future__ import annotations

import json
import sys

from . import (KEYWORDS_FILE, blocks_from_opi, check_store, keywords_from_opi,
               solvents_from_opi)

#: Keywords ORCA accepts that OPI's catalog does not list, each with the
#: manual page that uses it.
ALIASES = {
    "xtb": "contents/modelchemistries/semiempirical.html",
    "xtb0": "contents/modelchemistries/semiempirical.html",
    "xtb1": "contents/modelchemistries/semiempirical.html",
    "xtb2": "contents/modelchemistries/semiempirical.html",
}


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if args[:1] == ["refresh-keywords"]:
        from importlib import metadata
        from delfin.doctor import OPI_VERSION
        found = metadata.version("orca-pi")
        if found != OPI_VERSION:
            print(f"OPI {found} is installed; the catalog is written from {OPI_VERSION}",
                  file=sys.stderr)
            return 1
        old = json.loads(KEYWORDS_FILE.read_text()) if KEYWORDS_FILE.is_file() else {}
        KEYWORDS_FILE.write_text(json.dumps({
            "opi": found, "keywords": keywords_from_opi(),
            "solvents": solvents_from_opi(), "aliases": ALIASES,
            # Typed block models: block name -> option names OPI knows.
            "blocks": blocks_from_opi(),
            # Accepted by ORCA's own parser, not listed by OPI (kept as measured).
            "orca_parser": old.get("orca_parser", []),
            "blocks_orca_parser": old.get("blocks_orca_parser", []),
            "orca_parser_version": old.get("orca_parser_version", "")}, indent=0) + "\n")
        print(f"wrote {KEYWORDS_FILE}")
        return 0
    if args[:1] == ["check"]:
        problems = check_store()
        print("\n".join(problems) or "all cards valid")
        return 1 if problems else 0
    if args[:1] == ["probe"]:
        from . import probe as probe_mod
        words = args[1:]
        if not words:
            print("usage: python -m delfin.orca_knowledge probe KEYWORD [KEYWORD...]",
                  file=sys.stderr)
            return 2
        results, rc = probe_mod.probe_all(words)
        for word, result in zip(words, results):
            print(f"{word}: {result}")
        return rc
    if args[:1] == ["probe-input"]:
        return _cmd_probe_input(args[1:])
    if args[:1] == ["probe-queue"]:
        return _cmd_probe_queue(args[1:])
    print("usage: python -m delfin.orca_knowledge refresh-keywords | check | "
          "probe KEYWORD [KEYWORD...] | "
          "probe-input FILE.inp [--look-for TEXT ...] [--max-seconds N] | "
          "probe-queue DIR",
          file=sys.stderr)
    return 2


def _cmd_probe_input(args: list[str]) -> int:
    import argparse

    from . import probe as probe_mod

    parser = argparse.ArgumentParser(
        prog="python -m delfin.orca_knowledge probe-input")
    parser.add_argument("file", help="a complete, small ORCA input file")
    parser.add_argument("--look-for", nargs="+", default=[],
                        help="texts to report on in ORCA's output "
                             "(e.g. a default ORCA prints)")
    parser.add_argument("--max-seconds", type=int, default=None,
                        dest="max_seconds",
                        help="budget for the ORCA run (default "
                             f"{probe_mod.DEFAULT_INPUT_TIMEOUT_S})")
    opts = parser.parse_args(args)
    try:
        with open(opts.file, encoding="utf-8") as fh:
            text = fh.read()
    except OSError as exc:
        print(f"cannot read {opts.file}: {exc}", file=sys.stderr)
        return 2
    timeout = opts.max_seconds if opts.max_seconds is not None \
        else probe_mod.DEFAULT_INPUT_TIMEOUT_S
    result = probe_mod.probe_input(text, look_for=opts.look_for,
                                   timeout_s=timeout)
    _print_probe_input_result(result, sys.stdout)
    if result["verdict"] == "unavailable":
        return 2
    return 0 if result["verdict"] != "rejected" else 1


def _print_probe_input_result(result: dict, fh) -> None:
    print(f"verdict: {result['verdict']}")
    for message in result["messages"]:
        print(f"  {message}")
    for hit in result.get("look_for", []):
        state = "FOUND" if hit["found"] else "NOT FOUND"
        print(f"look-for {hit['text']!r}: {state} ({hit['count']}x)")
        for context in hit["context"]:
            for line in context.splitlines():
                print(f"    | {line}")


def _cmd_probe_queue(args: list[str]) -> int:
    import argparse
    import glob as globmod
    import json
    import os

    from . import probe as probe_mod

    parser = argparse.ArgumentParser(
        prog="python -m delfin.orca_knowledge probe-queue")
    parser.add_argument("dir", help="directory with probe request *.json")
    parser.add_argument("--max-seconds", type=int, default=None,
                        help="default budget, overridden per request")
    opts = parser.parse_args(args)
    paths = sorted(globmod.glob(os.path.join(opts.dir, "*.json")))
    # Result files written next to the requests must not be picked up again.
    paths = [p for p in paths if not p.endswith(".result.json")]
    if not paths:
        print(f"no probe requests (*.json) in {opts.dir}", file=sys.stderr)
        return 2
    rc = 0
    for path in paths:
        try:
            with open(path, encoding="utf-8") as fh:
                request = json.load(fh)
        except (OSError, ValueError) as exc:
            print(f"{path}: cannot read request: {exc}", file=sys.stderr)
            rc = 1
            continue
        try:
            result = probe_mod.probe_input(
                str(request.get("input", "")),
                look_for=list(request.get("look_for") or []),
                timeout_s=int(request.get("max_seconds")
                              or opts.max_seconds
                              or probe_mod.DEFAULT_INPUT_TIMEOUT_S))
        except (TypeError, ValueError) as exc:
            print(f"{path}: bad request: {exc}", file=sys.stderr)
            rc = 1
            continue
        out = {k: request.get(k) for k in ("card_id", "question")}
        out.update(verdict=result["verdict"], messages=result["messages"],
                   look_for=result["look_for"], timed_out=result["timed_out"])
        result_path = path[:-5] + ".result.json"
        with open(result_path, "w", encoding="utf-8") as fh:
            json.dump(out, fh, indent=2)
            fh.write("\n")
        print(f"{os.path.basename(path)}: {out['verdict']} -> "
              f"{os.path.basename(result_path)}")
        if result["verdict"] == "unavailable":
            rc = 2
    return rc


if __name__ == "__main__":
    sys.exit(main())
