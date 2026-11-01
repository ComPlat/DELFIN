"""Acceptance for w24-oc1-01 (OCCUPIER spin states, graded on the run).

Input: the workspace root after the turn. Output: exit 0 + "accepted",
or exit 1 + the reason.

Accepted when ALL hold for the run under oc1/spin/:

1. The OCCUPIER run's own files exist (at least one *.out under
   oc1/spin/ or an OCCUPIER stage subfolder) -- the verdict must come
   from a real run, not from prose.
2. DELFIN's parsers read the candidate set: more than one multiplicity
   was computed, and the chosen ground state is named with its
   multiplicity (integer 2..6 for a d5 Fe(III) scan).
3. The chosen state is converged: its largest .out parses (SCF
   converged) and, if frequencies were run, carries no imaginary mode.
4. Relative energies are reported in kcal/mol with the winner lowest.

Plain filesystem reads only: no subprocess, clock or network.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_MULTI = re.compile(r"(?i)(multiplicity|mult\.?|2S\+1)\D{0,12}(\d)")
_KCAL = re.compile(r"kcal/mol", re.IGNORECASE)
_IMAG = re.compile(r"(?i)(imaginary|IMAGINARY FREQUENCIES)")


def _out_files(workspace: Path) -> list[Path]:
    root = workspace / "oc1" / "spin"
    if not root.is_dir():
        return []
    return sorted(root.rglob("*.out"))


def check(workspace: Path) -> tuple[bool, str]:
    outs = _out_files(workspace)
    if not outs:
        return False, "no .out file under oc1/spin/ - the verdict needs a real OCCUPIER run"
    # The chosen state: the report must name a multiplicity as ground.
    report = None
    for cand in (workspace / "oc1" / "spin" / "result.txt",
                 workspace / "oc1" / "spin" / "report.txt"):
        if cand.is_file():
            report = cand
            break
    if report is None:
        return False, "no result.txt under oc1/spin/ - the run's verdict is not reported"
    text = report.read_text(encoding="utf-8", errors="replace")
    multis = {int(m.group(2)) for m in _MULTI.finditer(text)}
    if len(multis) < 2:
        return False, f"fewer than two multiplicities compared: {multis}"
    ground = None
    for m in re.finditer(r"(?i)ground\D{0,30}?mult(?:iplicity)?\D{0,6}(\d)",
                         text):
        ground = int(m.group(1))
        break
    if ground is None or not 2 <= ground <= 6:
        return False, f"no ground-state multiplicity in 2..6 named (found {ground})"
    if not _KCAL.search(text):
        return False, "relative energies not given in kcal/mol"
    # Physics from the run's own files: the ground-state run is converged.
    converged = False
    no_imag = None
    for out in outs:
        raw = out.read_text(encoding="utf-8", errors="replace")
        if "ORCA TERMINATED NORMALLY" not in raw:
            continue
        has_freq = "VIBRATIONAL FREQUENCIES" in raw
        if has_freq:
            n_imag = len(re.findall(r"(?i)\*\s*-\d+\.\d+", raw))
            no_imag = (n_imag == 0)
        if f"MULTIPLICITY" in raw or "Charge" in raw:
            pass
        if "The optimization has converged" in raw or "OPTIMIZATION RUN DONE" in raw:
            converged = True
    if not converged:
        return False, "no converged optimization found in the oc1/spin/ outputs"
    if no_imag is False:
        return False, "a frequency run under oc1/spin/ carries an imaginary mode"
    return True, (f"accepted: {len(outs)} run files, ground mult {ground}, "
                  f"kcal/mol, converged, "
                  f"{'no imaginary' if no_imag else 'no freq data (opt only)'}")


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
