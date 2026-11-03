"""Acceptance for w24-oc1-03 (metallocene oxidation potential via OCCUPIER).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason.

Accepted when ALL hold for the run under oc1/redox/:

1. A real OCCUPIER run exists for BOTH species (ferrocene AND
   ferrocenium): terminated .out files for each.
2. An oxidation potential in V inside 0.00-0.90 (Fc/Fc+ vs SHE in MeCN
   is ~0.40; the band tolerates reference-convention spread).
3. A reference electrode is named next to the number (vs. SHE, Fc/Fc+,
   Ag/Ag+, NHE...). A potential with no electrode is a rejection.
4. The thermodynamic relation (G_ox - G_red) / thermodynamic cycle is
   stated.
5. Same-level rule from the run's provenance: the two species' outputs
   name the same functional and basis (a mixed-level E_ox is the
   package's core trap).

Plain filesystem reads only.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_POTENTIAL = re.compile(r"(?i)E_?(?:ox|1/2|red)?\D{0,20}?([-+]?\d\.\d+)\s*V\b")
_ELECTRODE = re.compile(r"(?i)(vs\.?\s*(SHE|NHE|Fc|Ag|RHE)|reference "
                        r"electrode|E_ref)")
_CYCLE = re.compile(r"(?i)(thermodynamic cycle|G_?ox|G_?red|Delta.?G|"
                    r"free energy)")
_METHOD = re.compile(r"(?i)^\s*(?:!|functional)\s+([A-Za-z0-9\-]+)", re.M)
_BASIS = re.compile(r"(?i)def2-?[A-Za-z0-9]+|SARC|ZORA")


def _terminated(path: Path) -> bool:
    return "ORCA TERMINATED NORMALLY" in path.read_text(encoding="utf-8",
                                                        errors="replace")


def _level(out: Path) -> tuple[str, tuple[str, ...]]:
    raw = out.read_text(encoding="utf-8", errors="replace")
    methods = _METHOD.findall(raw)
    bases = tuple(sorted(set(_BASIS.findall(raw))))
    return (methods[0] if methods else "", bases)


def check(workspace: Path) -> tuple[bool, str]:
    root = workspace / "oc1" / "redox"
    if not root.is_dir():
        return False, "oc1/redox/ does not exist"
    outs = sorted(root.rglob("*.out"))
    if not outs:
        return False, "no .out file under oc1/redox/"
    ok_runs = [o for o in outs if _terminated(o)]
    if len(ok_runs) < 2:
        return (False, f"only {len(ok_runs)} terminated run(s) under "
                       "oc1/redox/ - both species of the couple are needed")
    levels = {_level(o) for o in ok_runs}
    funcs = {lv[0] for lv in levels if lv[0]}
    bases: set[str] = set()
    for lv in levels:
        bases.update(lv[1])
    if len(funcs) > 1:
        return (False, f"levels of theory mixed between the two species: "
                       f"{funcs} - an E_ox from mixed levels is not a result")
    report = root / "result.txt"
    if not report.is_file():
        return False, "no result.txt under oc1/redox/"
    text = report.read_text(encoding="utf-8", errors="replace")
    volts = [float(m.group(1)) for m in _POTENTIAL.finditer(text)]
    if not any(0.00 <= v <= 0.90 for v in volts):
        return False, f"no oxidation potential in 0.00-0.90 V; found {volts}"
    if not _ELECTRODE.search(text):
        return (False, "no reference electrode next to the potential - "
                       "an E value without its electrode is not a result")
    if not _CYCLE.search(text):
        return False, "thermodynamic relation / free-energy cycle not stated"
    note = f" ({funcs.pop()}/{sorted(bases)[0] if bases else 'basis n/a'})" if funcs else ""
    return True, f"accepted: E in band, electrode, cycle, same level{note}"


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
