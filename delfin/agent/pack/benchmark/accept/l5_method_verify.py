#!/usr/bin/env python3
"""Acceptance for l5_name_method_with_basis_pairing.

Input: the workspace root (argv[1]). Verdict = the choice file and the
fixture:

* l5_method_notes/brief.txt is byte-for-byte what the setup wrote;
* l5_method_notes/choice.txt names exactly one density functional from
  the list below and exactly one basis set (def2, Pople, Dunning or
  pcseg families), adds no dispersion correction (D3, D3BJ, D4, a -3c
  composite, NL, VV10: the brief excludes dispersion), and has an ORCA
  keyword line (text after a ``!``) that carries both names.

Names are matched case-insensitively as whole names, so the wording around
them is free; a second functional or basis ("B3LYP or PBE0") fails because
the brief asks for one combination.

Exit 0 when all hold, 1 otherwise; 2 on a missing argument.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

# Must equal l5_method_setup._BRIEF.
_BRIEF = ("Recommend ONE density-functional/basis-set combination for a "
          "closed-shell molecule with only H, C, N and O. No solvent, no "
          "dispersion correction. Name the functional, the basis set, and "
          "the keyword line you would use in ORCA, and write them to "
          "l5_method_notes/choice.txt.\n")

_FUNCTIONALS = (
    "cam-b3lyp", "b3lyp", "b3pw91", "x3lyp", "bhandhlyp", "blyp", "bp86",
    "revpbe0", "revpbe", "pbe0", "pbe", "tpssh", "tpss0", "tpss",
    "m06-2x", "m062x", "m06-l", "m06l", "m06", "mn15", "pw6b95",
    "wb97x-v", "wb97m-v", "wb97x", "ωb97x", "b97m-v",
    "r2scan0", "r2scan", "scan", "b2plyp",
)
_NAME_EDGE_L = r"(?<![\w\-ω])"
_NAME_EDGE_R = r"(?![\w\-])"
_FUNCTIONAL = re.compile(
    _NAME_EDGE_L + "(" + "|".join(re.escape(f) for f in _FUNCTIONALS) + ")"
    + _NAME_EDGE_R, re.I)
_BASIS = re.compile(
    _NAME_EDGE_L + r"((?:ma-)?def2-(?:sv\(p\)|svpd?|tzvpp?d?|qzvpp?d?)"
    r"|6-311?\+{0,2}g(?:\*{1,2}|\([^)]*\))?"
    r"|(?:aug-)?cc-pw?cv[dtq5]z|(?:aug-)?cc-pv[dtq5]z"
    r"|pcseg-[0-4])" + r"(?![\w\-])", re.I)
_DISPERSION = re.compile(
    r"(?<![a-z0-9])(?:d3(?:bj|zero)?|d4|nl|vv10)(?![a-z0-9])|-3c(?![a-z0-9])",
    re.I)


def _norm(name: str) -> str:
    return name.lower().replace("ω", "w").replace("m062x", "m06-2x") \
        .replace("m06l", "m06-l")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: l5_method_verify.py <workspace>", file=sys.stderr)
        return 2
    d = Path(argv[1]).resolve() / "l5_method_notes"
    brief = d / "brief.txt"
    choice = d / "choice.txt"
    if not brief.is_file() or brief.read_text(encoding="utf-8") != _BRIEF:
        print("FAIL: brief.txt was changed or removed", file=sys.stderr)
        return 1
    if not choice.is_file():
        print(f"FAIL: no {choice}", file=sys.stderr)
        return 1
    text = choice.read_text(encoding="utf-8")
    functionals = {_norm(m) for m in _FUNCTIONAL.findall(text)}
    bases = {m.lower() for m in _BASIS.findall(text)}
    if len(functionals) != 1:
        print(f"FAIL: {len(functionals)} functionals named "
              f"({sorted(functionals)}), want exactly one", file=sys.stderr)
        return 1
    if len(bases) != 1:
        print(f"FAIL: {len(bases)} basis sets named ({sorted(bases)}), "
              "want exactly one", file=sys.stderr)
        return 1
    if _DISPERSION.search(text):
        print("FAIL: a dispersion correction is named; the brief excludes "
              "dispersion", file=sys.stderr)
        return 1
    (functional,), (basis,) = functionals, bases
    for line in re.findall(r"!([^\n]*)", text):
        if ({_norm(m) for m in _FUNCTIONAL.findall(line)} == {functional}
                and {m.lower() for m in _BASIS.findall(line)} == {basis}):
            print(f"OK: one combination, {functional} / {basis}")
            return 0
    print("FAIL: no ORCA keyword line (after '!') with the functional and "
          "the basis set", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
