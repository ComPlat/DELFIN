"""``python -m delfin.orca_knowledge refresh-keywords``: rewrite the keyword
catalog from the installed OPI (it must be the pinned version), and
``python -m delfin.orca_knowledge check``: the problems of the card store."""

from __future__ import annotations

import json
import sys

from . import KEYWORDS_FILE, check_store, keywords_from_opi, solvents_from_opi

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
            # Accepted by ORCA's own parser, not listed by OPI (kept as measured).
            "orca_parser": old.get("orca_parser", []),
            "orca_parser_version": old.get("orca_parser_version", "")}, indent=0) + "\n")
        print(f"wrote {KEYWORDS_FILE}")
        return 0
    if args[:1] == ["check"]:
        problems = check_store()
        print("\n".join(problems) or "all cards valid")
        return 1 if problems else 0
    print("usage: python -m delfin.orca_knowledge refresh-keywords | check", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
