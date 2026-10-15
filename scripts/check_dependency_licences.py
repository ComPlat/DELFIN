#!/usr/bin/env python3
"""Check the licences of an installed dependency set against DELFIN's policy.

Input: the JSON that ``pip-licenses --from=mixed --format=json`` writes for an
environment, and the reviewed allow-list ``scripts/dependency_licences.toml``.
Output: one line per finding on stdout, exit status 1 when there is any, 0
otherwise, 2 on a usage error.

Policy, for the environment of a plain ``pip install delfin-complat`` (no
extras):

- a licence string that is permissive (MIT, BSD, Apache, ISC, PSF, MPL-2.0,
  Zlib, ...) or LGPL passes;
- anything else -- AGPL, SSPL, a GPL, a proprietary or non-commercial term,
  ``UNKNOWN`` or an empty string -- fails unless the package is in the
  allow-list with exactly that licence string. The string is compared, not
  only the name, so an upstream licence change reaches a person again;
- an allow-list entry for a package that is not installed is reported as a
  note, not a failure, because the same list serves environments with and
  without extras.

The licence strings are what the packages declare; pip-licenses reports the
trove classifier when there is one and the License field otherwise
(``--from=mixed``). A declaration can be wrong or missing, which is why an
``UNKNOWN`` needs a reviewed entry naming the licence file that was read.

Rejected alternatives: GitHub's dependency-review deny list (the dependency
graph records no licence for pymupdf or openbabel-wheel, so it never fired
for either); ``pip-licenses --fail-on`` alone (it matches substrings with no
reviewed exceptions, and ``--allow-only`` cannot express "this one GPL
package, reviewed, for this reason").
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import tomllib

DEFAULT_ALLOWLIST = Path(__file__).with_name("dependency_licences.toml")

#: Licence strings that pass without review. Matched case-insensitively
#: against the whole declared string after the deny check below.
_PERMISSIVE = re.compile(
    r"\b(MIT|BSD|Apache|ISC|PSF|Python Software Foundation|MPL[- ]?2\.0|"
    r"Mozilla Public License 2\.0|Zlib|Unlicense|CC0|0BSD|HPND|BlueOak|"
    r"Historical Permission Notice|LGPL|Lesser General Public License)",
    re.IGNORECASE,
)

#: Terms that are never accepted without a reviewed entry, even when the
#: string also names a permissive licence ("MIT AND AGPL-3.0" is AGPL).
_NEEDS_REVIEW = re.compile(
    r"(AGPL|Affero|SSPL|Server Side Public|Commons Clause|BUSL|Business Source|"
    r"Elastic License|Proprietary|Commercial|Non[- ]?Commercial|\bNC\b|"
    r"(?<!Lesser )General Public License|(?<![AL])GPL|UNKNOWN)",
    re.IGNORECASE,
)


def classify(licence: str) -> str:
    """``"pass"`` or ``"review"`` for one declared licence string."""
    text = (licence or "").strip()
    if not text:
        return "review"
    if _NEEDS_REVIEW.search(text):
        return "review"
    return "pass" if _PERMISSIVE.search(text) else "review"


def load_allowlist(path: Path) -> dict[str, dict]:
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    out = {}
    for entry in data.get("package", []):
        name = _norm(entry["name"])
        for key in ("licence", "spdx", "reason", "reviewed"):
            if not str(entry.get(key, "")).strip():
                raise ValueError(f"{path}: entry {entry['name']!r} has no {key}")
        out[name] = entry
    return out


def _norm(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def check(packages: list[dict], allowlist: dict[str, dict]) -> tuple[list[str], list[str]]:
    """``(failures, notes)`` for the pip-licenses records *packages*."""
    failures, notes, seen = [], [], set()
    for pkg in packages:
        name, licence = pkg.get("Name", ""), str(pkg.get("License", ""))
        key = _norm(name)
        seen.add(key)
        if classify(licence) == "pass":
            continue
        entry = allowlist.get(key)
        if entry is None:
            failures.append(f"{name} {pkg.get('Version', '')}: licence "
                            f"{licence!r} needs review (not in the allow-list)")
        elif entry["licence"] != licence:
            failures.append(f"{name} {pkg.get('Version', '')}: declares "
                            f"{licence!r}, the reviewed entry says "
                            f"{entry['licence']!r}; review it again")
    for key, entry in sorted(allowlist.items()):
        if key not in seen:
            notes.append(f"note: allow-list entry {entry['name']} is not installed here")
    return failures, notes


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("licences_json", type=Path,
                        help="output of pip-licenses --from=mixed --format=json")
    parser.add_argument("--allowlist", type=Path, default=DEFAULT_ALLOWLIST)
    args = parser.parse_args(argv)
    try:
        packages = json.loads(args.licences_json.read_text(encoding="utf-8"))
        allowlist = load_allowlist(args.allowlist)
    except (OSError, ValueError, KeyError) as exc:
        print(f"usage error: {exc}", file=sys.stderr)
        return 2
    failures, notes = check(packages, allowlist)
    for line in notes + failures:
        print(line)
    print(f"{len(packages)} packages, {len(failures)} findings")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
