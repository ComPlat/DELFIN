"""The dependency closure DELFIN installs, as exact pins for pip-audit.

Input: the JSON report of
``pip install --dry-run --ignore-installed --report REPORT ".[extras]"``
(pip's installation report, format version 1), i.e. what pip would put into
an empty environment for this checkout.

Output: ``name==version`` lines on stdout, one per distribution pip would
download, sorted; the project itself (a ``file://`` source) is left out.
Exit 1 when fewer than ``--minimum`` distributions remain, or when one of
the ``--expect`` names is missing: an empty requirements file makes
pip-audit print "No known vulnerabilities found" and exit 0, which is the
report of an audit that looked at nothing.

Semantics: the set audited is the resolved runtime closure of the package
and the named extras, not the runner's interpreter. Auditing the job's own
environment (``pip-audit`` with no arguments, as ci.yml did) also reports
the runner image's setuptools and pip, and on 2026-10-15 its only finding
was setuptools 79.0.1 (PYSEC-2026-3447) from the hosted Python, which no
user of DELFIN installs from DELFIN.

Rejected: installing into a fresh venv and auditing its site-packages
(``--path``). A 3.11 venv seeds pip and setuptools, which reintroduces the
installer tooling, and the install costs minutes the resolution does not.
"""

from __future__ import annotations

import argparse
import json
import sys


def closure(report: dict) -> list[str]:
    pins = set()
    for item in report.get("install", []):
        url = (item.get("download_info") or {}).get("url", "")
        if url.startswith("file:"):
            continue
        meta = item["metadata"]
        pins.add(f"{meta['name']}=={meta['version']}")
    return sorted(pins, key=str.lower)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("report")
    parser.add_argument("--minimum", type=int, default=1)
    parser.add_argument("--expect", action="append", default=[])
    args = parser.parse_args(argv)

    with open(args.report, encoding="utf-8") as handle:
        pins = closure(json.load(handle))
    names = {p.split("==", 1)[0].lower().replace("_", "-") for p in pins}
    missing = [n for n in args.expect if n.lower().replace("_", "-") not in names]
    if len(pins) < args.minimum or missing:
        print(f"closure has {len(pins)} distributions (minimum {args.minimum}); "
              f"missing: {missing}", file=sys.stderr)
        return 1
    sys.stdout.write("".join(p + "\n" for p in pins))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
