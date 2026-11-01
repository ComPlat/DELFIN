#!/usr/bin/env python3
"""Compare OPI (orca-pi 2.0.0) energy reads against DELFIN's classic parsers
on the ORCA 6.1.1 fixtures under tests/fixtures/orca611/.

A later phase may switch DELFIN's readers onto OPI -- but only where this
comparison shows the same values within tolerance on real outputs, with the
old parser kept as the fallback for ORCA < 6.1.1 and for failed runs.  This
script is the evidence for that decision: run it after touching any reader.

Read-only: it never writes into the fixtures (OPI is pointed at a scratch copy
so create_property_json / orca_2json would not touch the archive or the
fixtures).
"""

from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

TOL_EH = 1e-8  # energy agree to 1e-8 Hartree
FIXTURES = Path("tests/fixtures/orca611")


def _classic_sp_energy(out_text: str) -> float | None:
    """DELFIN's classic total-energy read: the FINAL SINGLE POINT ENERGY line."""
    for line in reversed(out_text.splitlines()):
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                return float(line.split()[-1])
            except ValueError:
                return None
    return None


def _classic_gibbs(out_text: str) -> float | None:
    """DELFIN's classic Gibbs read for frequency runs."""
    for line in reversed(out_text.splitlines()):
        if "Final Gibbs free energy" in line or "Final Gibbs free enthalpy" in line:
            try:
                return float(line.split()[-2])
            except ValueError:
                return None
    return None


def main() -> int:
    from opi.output.core import Output

    print(f"{'run':<20}{'files':<8} {'terminated':<10}{'scf':<6}"
          f"{'classic_Eh':>14}{'opi_Eh':>14} {'match':<6}")
    fails = 0
    for run_dir in sorted(FIXTURES.iterdir()):
        if not run_dir.is_dir():
            continue
        name = run_dir.name
        # Work on a scratch copy so OPI generation never touches the fixture.
        with tempfile.TemporaryDirectory() as tmp:
            wd = Path(tmp) / name
            wd.mkdir()
            for f in run_dir.glob(f"{name}.*"):
                shutil.copy2(f, wd / f.name)
            out_path = wd / f"{name}.out"
            out_text = out_path.read_text(errors="replace")

            classic_e = _classic_sp_energy(out_text)
            classic_gibbs = _classic_gibbs(out_text)

            try:
                opi_out = Output(basename=name, working_dir=wd)
                opi_out.parse()
                terminated = opi_out.terminated_normally()
                scf = opi_out.scf_converged()
                opi_e = opi_out.get_final_energy()
            except Exception as exc:  # OPI could not read this run
                print(f"{name:<20}{'':<8} {'ERROR':<10}{'':<6}{'':>14}{'':>14} "
                      f"({exc})")
                fails += 1
                continue

            # A converged run must have both energies agreeing on the SP total.
            if scf and terminated and classic_e is not None and opi_e is not None:
                match = abs(classic_e - opi_e) <= TOL_EH
                if not match:
                    fails += 1
                print(f"{name:<20}{'ok':<8} {str(terminated):<10}{str(scf):<6}"
                      f"{classic_e:>14.8f}{opi_e:>14.8f} "
                      f"{'OK' if match else 'MISMATCH'}")
            else:
                # Failed or non-SP run: report flags; numbers may be absent.
                print(f"{name:<20}{'ok':<8} {str(terminated):<10}{str(scf):<6}"
                      f"{'-':>14}{'-' if opi_e is None else f'{opi_e:.8f}':>14} "
                      f"{'(failed)' if not (scf and terminated) else ''}")

    print(f"\nDone. Mismatches / read errors: {fails}")
    return 0 if fails == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
