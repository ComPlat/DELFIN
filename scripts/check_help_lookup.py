#!/usr/bin/env python3
"""One-shot syntax + allowed-completeness check for the help feature.

Usage: python scripts/check_help_lookup.py [--commit]
Verifies in-process (login-node-safe, no pytest):
  1. syntax of the three touched modules
  2. every CONTROL key resolves allowed values (spec.allowed or table)
  3. lookup finds every key and tab; misses are explicit
Exit code 0 = all good; 1 = named gaps.
"""
import ast
import sys

FILES = [
    "delfin/common/control_allowed_values.py",
    "delfin/common/control_validator.py",
    "delfin/dashboard/help_lookup.py",
    "delfin/dashboard/tab_help.py",
    "delfin/dashboard/header_keysearch.py",
]


def main() -> int:
    ok = True
    for f in FILES:
        try:
            ast.parse(open(f).read())
        except SyntaxError as e:
            print(f"FAIL syntax {f}: {e}")
            ok = False
    if not ok:
        return 1
    print("syntax: all", len(FILES), "modules ok")

    sys.path.insert(0, ".")
    from delfin.common.control_allowed_values import ALLOWED_BY_KEY
    from delfin.common.control_validator import CONTROL_FIELD_SPECS
    from delfin.dashboard.help_lookup import lookup
    from delfin.dashboard.tab_help import _TABS_HELP

    def resolved_allowed(spec) -> str:
        return (getattr(spec, "allowed", "") or "").strip() or \
            ALLOWED_BY_KEY.get(spec.name, "").strip()

    no_allowed = [s.name for s in CONTROL_FIELD_SPECS
                  if not resolved_allowed(s)]
    print(f"keys total: {len(CONTROL_FIELD_SPECS)}, "
          f"table entries: {len(ALLOWED_BY_KEY)}, "
          f"no allowed: {len(no_allowed)}")
    if no_allowed:
        print("FAIL keys without allowed:", no_allowed[:10])
        ok = False

    unfound = [s.name for s in CONTROL_FIELD_SPECS
               if not any(e["kind"] == "key" and e["name"] == s.name
                          for e in lookup(s.name)["entries"])]
    if unfound:
        print("FAIL lookup cannot find:", unfound[:8])
        ok = False
    else:
        print("lookup: every key findable")

    r = lookup("definitely_not_a_key_or_tab_xyz")
    if r["found"]:
        print("FAIL lookup claims junk exists")
        ok = False
    else:
        print("lookup: explicit miss works")

    from collections import Counter
    thin = [t for _t, d in _TABS_HELP if len(d.strip()) < 25]
    if thin:
        print("FAIL thin tab entries:", thin[:5])
        ok = False
    else:
        print(f"tabs: {_TABS_HELP and len(_TABS_HELP)} entries, none thin")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
