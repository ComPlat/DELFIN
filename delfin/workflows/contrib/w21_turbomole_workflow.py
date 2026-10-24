"""Contributed workflow: Turbomole optimize + frequencies (package TM1).

Built on the wave-21 foundation (WB2 user-workflow path, WB3-registered
turbomole adapters) — no parallel mechanism. The workflow takes a
geometry (XYZ), runs a Turbomole geometry optimization (jobex), then a
frequency calculation (ridft + aoforce) to verify the minimum, and
surfaces energy, geometry and hessian.

Physical realism — every limit is enforced here, not just described:
- charge and multiplicity must be ints, mult >= 1;
- an open shell (mult > 1) is passed through so the adapter writes
  ``$unpaired N`` (never a silently closed-shell singlet);
- method and basis are passed through unchanged — Turbomole keywords
  are the caller's responsibility and are validated by the program;
  nothing is silently rewritten.
"""

from __future__ import annotations

import os
import shutil
import stat
from pathlib import Path
from typing import Any, Dict, Optional

from delfin.workflows.registry import register


def _check_charge_mult(charge: Any, mult: Any) -> Optional[str]:
    """Same physical limits as the adapters (single source of truth is
    delfin/tools/adapters/turbomole.py — imported there, not copied)."""
    from delfin.tools.adapters.turbomole import _check_charge_mult_types

    return _check_charge_mult_types(charge, mult)


class TurbomoleWorkflow:
    """Turbomole: optimize + frequencies on a given XYZ geometry."""

    name = "turbomole"
    description = (
        "Turbomole DFT workflow: geometry optimization (jobex) then "
        "frequency verification (ridft + aoforce) on a given XYZ geometry"
    )

    def run(self, *, config: Dict[str, Any], geometry: Optional[str] = None,
            charge: int = 0, mult: int = 1, method: str = "b3-lyp",
            basis: str = "def2-SVP", cores: int = 1,
            work_dir: Optional[str] = None, **kwargs: Any) -> Dict[str, Any]:
        from delfin.tools._registry import get as get_adapter

        bad = _check_charge_mult(charge, mult)
        if bad:
            raise ValueError(f"turbomole workflow: {bad}")

        if not geometry:
            raise ValueError(
                "turbomole workflow needs a geometry (XYZ path) — no "
                "SMILES step is provided; build one with the smiles_to_xyz "
                "step or pass an existing geometry")

        if work_dir:
            wd = Path(work_dir)
        else:
            wd = Path(config.get("work_dir") or ".") / "turbomole_run"
        wd.mkdir(parents=True, exist_ok=True)
        geom_path = Path(geometry)

        results: Dict[str, Any] = {"ok": True, "steps": {}}

        opt = get_adapter("turbomole_opt")
        freq = get_adapter("turbomole_freq")
        if opt is None or freq is None:
            raise RuntimeError(
                "turbomole adapters not registered — the WB3 registry walk "
                "did not pick up delfin.tools.adapters.turbomole")

        opt_res = opt.execute(
            wd, geometry=geom_path, cores=cores,
            charge=charge, mult=mult, method=method, basis=basis,
            **kwargs)
        results["steps"]["opt"] = {
            "status": opt_res.status.name, "error": opt_res.error,
            "energy_Eh": opt_res.data.get("energy_Eh"),
            "geometry": str(opt_res.geometry) if opt_res.geometry else None,
        }
        if opt_res.status.name != "SUCCESS" or not opt_res.geometry:
            results["ok"] = False
            results["errors"] = [f"opt failed: {opt_res.error}"]
            return results

        freq_res = freq.execute(
            wd, geometry=opt_res.geometry, cores=cores,
            charge=charge, mult=mult, method=method, basis=basis,
            **kwargs)
        results["steps"]["freq"] = {
            "status": freq_res.status.name, "error": freq_res.error,
            "energy_Eh": freq_res.data.get("energy_Eh"),
            "hessian": freq_res.data.get("hessian"),
        }
        if freq_res.status.name != "SUCCESS":
            results["ok"] = False
            results["errors"] = [f"freq failed: {freq_res.error}"]

        # Physical-realism reporting duty: a minimum is verified by the
        # frequency run; report, never hide.
        imaginary = self._count_imaginary(wd)
        if imaginary is not None:
            results["imaginary_frequencies"] = imaginary
            if imaginary > 0:
                results.setdefault("warnings", []).append(
                    f"{imaginary} imaginary frequencies — NOT a verified "
                    "minimum; reported, not hidden")
        return results

    @staticmethod
    def _count_imaginary(work_dir: Path) -> Optional[int]:
        """Count imaginary modes in the aoforce output, if present.

        Aoforce prints imaginary modes either as a bare negative number
        or with an explicit ``i`` suffix (``-123.45i``) — reviewer
        finding F5: an un-parsed i-suffix line would hide a
        non-minimum, which is exactly what the reporting duty exists
        to prevent. Both formats count.

        The rule is deliberately loose: any second token that parses
        as a float and is either i-suffixed or negative counts. In
        the aoforce mode table the second column IS the frequency,
        so a line like ``1  -2.5  comment`` is a mode, not prose.
        Until a real aoforce.out fixture lands this comment documents
        the accepted ambiguity: over-counting keeps the minimum
        UNVERIFIED — fail-safe, never silently green.
        """
        log = work_dir / "aoforce.out"
        if not log.is_file():
            return None
        count = 0
        for line in log.read_text(errors="replace").splitlines():
            parts = line.split()
            if len(parts) < 2:
                continue
            token = parts[1].lower()
            has_i_suffix = token.endswith("i") or token.endswith("i*")
            bare = token.rstrip("*i")
            try:
                value = float(bare)
            except ValueError:
                continue
            if has_i_suffix or value < 0:
                count += 1
        return count

    def run_cli(self, argv: list) -> int:
        import argparse

        parser = argparse.ArgumentParser(description=self.description)
        parser.add_argument("geometry", help="XYZ geometry file")
        parser.add_argument("--charge", type=int, default=0)
        parser.add_argument("--mult", type=int, default=1)
        parser.add_argument("--method", default="b3-lyp")
        parser.add_argument("--basis", default="def2-SVP")
        parser.add_argument("--cores", type=int, default=1)
        parser.add_argument("--work-dir", default=None)
        args = parser.parse_args(argv)
        try:
            result = self.run(
                config={}, geometry=args.geometry, charge=args.charge,
                mult=args.mult, method=args.method, basis=args.basis,
                cores=args.cores, work_dir=args.work_dir)
        except Exception:  # noqa: BLE001
            return 1
        return 0 if result.get("ok") else 1


register(TurbomoleWorkflow())
