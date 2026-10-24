"""Contributed workflow: explicit solvation (package SV1).

Built on the wave-21 foundation — WB3-registered adapters (packmol_solvate,
xtb_solvator), WB2's user-workflow path, TM1's contrib pattern — no parallel
mechanism. The workflow takes a solute geometry, builds a Packmol solvent
box around it, relaxes the box at xTB level, then adds an ORCA-xTB
SOLVATOR explicit shell on the relaxed geometry.

Physical realism — every limit is declared and enforced, not described:
- the box must be able to hold the solute plus the solvent (Packmol 2009,
  DOI 10.1002/jcc.21224) — enforced in the packmol adapter;
- SOLVATOR is an xTB feature: the SOLVATOR method must be GFN0/1/2-xTB
  (ORCA manual, %SOLVATOR block) — enforced in the xtb_solvator adapter;
- charge/mult are ints, mult >= 1 — enforced in the adapters;
- the xTB pre-relaxation is a structure preparation, NOT the production
  level: the result reports both levels so a downstream step never
  confuses an xTB-relaxed energy with a DFT result.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from delfin.workflows.registry import register


class SolvationWorkflow:
    """Explicit solvation: Packmol box -> xTB relax -> ORCA-xTB SOLVATOR."""

    name = "explicit_solvation"
    description = (
        "Explicit solvation workflow: Packmol solvent box around a solute, "
        "xTB pre-relaxation of the box, then an ORCA-xTB SOLVATOR shell "
        "on the relaxed geometry"
    )

    def run(self, *, config: Dict[str, Any], geometry: Optional[str] = None,
            solvent_xyz: Optional[str] = None, solvent: str = "water",
            n_solvent: int = 50, box_side: float = 20.0,
            tolerance: float = 2.0, nsolv: int = 12,
            solvator_method: str = "GFN2-xTB", charge: int = 0,
            mult: int = 1, cores: int = 1,
            work_dir: Optional[str] = None, **kwargs: Any) -> Dict[str, Any]:
        from delfin.tools._registry import get as get_adapter

        if not geometry:
            raise ValueError(
                "explicit_solvation workflow needs a solute geometry "
                "(XYZ path); build one with the smiles_to_xyz step or "
                "pass an existing geometry")
        if not solvent_xyz:
            raise ValueError(
                "explicit_solvation workflow needs solvent_xyz (path to a "
                "single solvent molecule XYZ for the Packmol box)")

        if work_dir:
            wd = Path(work_dir)
        else:
            wd = Path(config.get("work_dir") or ".") / "solvation_run"
        wd.mkdir(parents=True, exist_ok=True)
        geom_path = Path(geometry)

        results: Dict[str, Any] = {"ok": True, "steps": {}}

        packmol = get_adapter("packmol_solvate")
        xtb_opt = get_adapter("xtb_opt")
        solvator = get_adapter("xtb_solvator")
        missing = [n for n, a in (
            ("packmol_solvate", packmol), ("xtb_opt", xtb_opt),
            ("xtb_solvator", solvator)) if a is None]
        if missing:
            raise RuntimeError(
                "solvation adapters not registered (the WB3 registry walk "
                f"did not pick up delfin/tools/adapters): {missing}")

        # ── Step 1: Packmol box ────────────────────────────────────────
        pm_res = packmol.execute(
            wd, geometry=geom_path, cores=cores,
            solvent_xyz=solvent_xyz, n_solvent=n_solvent,
            box_side=box_side, tolerance=tolerance, **kwargs)
        results["steps"]["packmol_box"] = {
            "status": pm_res.status.name, "error": pm_res.error,
            "geometry": str(pm_res.geometry) if pm_res.geometry else None,
        }
        if pm_res.status.name != "SUCCESS" or not pm_res.geometry:
            results["ok"] = False
            results["errors"] = [f"packmol box failed: {pm_res.error}"]
            return results

        # ── Step 2: xTB pre-relaxation of the box (structure only) ─────
        opt_res = xtb_opt.execute(
            wd, geometry=pm_res.geometry, cores=cores,
            charge=charge, mult=mult, **kwargs)
        results["steps"]["xtb_relax"] = {
            "status": opt_res.status.name, "error": opt_res.error,
            "geometry": str(opt_res.geometry) if opt_res.geometry else None,
        }
        relaxed = opt_res.geometry if (
            opt_res.status.name == "SUCCESS" and opt_res.geometry
        ) else pm_res.geometry
        if opt_res.status.name != "SUCCESS":
            results.setdefault("warnings", []).append(
                "xTB pre-relaxation failed; the SOLVATOR step uses the "
                "unrelaxed Packmol box — report, not hide")

        # ── Step 3: ORCA-xTB SOLVATOR shell ────────────────────────────
        solv_res = solvator.execute(
            wd, geometry=relaxed, cores=cores,
            solvent=solvent, charge=charge, mult=mult,
            n_solvent=nsolv, method=solvator_method, **kwargs)
        results["steps"]["xtb_solvator"] = {
            "status": solv_res.status.name, "error": solv_res.error,
            "geometry": str(solv_res.geometry) if solv_res.geometry else None,
        }
        if solv_res.status.name != "SUCCESS" or not solv_res.geometry:
            results["ok"] = False
            results["errors"] = [f"solvator failed: {solv_res.error}"]
            return results

        results["levels"] = {
            "box_relaxation": "GFN-xTB (structure preparation only)",
            "explicit_shell": str(solvator_method),
        }
        results["solvated_geometry"] = str(solv_res.geometry)
        return results

    def run_cli(self, argv: list) -> int:
        import argparse

        parser = argparse.ArgumentParser(description=self.description)
        parser.add_argument("geometry", help="solute XYZ geometry file")
        parser.add_argument("solvent_xyz", help="solvent molecule XYZ file")
        parser.add_argument("--solvent", default="water")
        parser.add_argument("--n-solvent", type=int, default=50)
        parser.add_argument("--box-side", type=float, default=20.0)
        parser.add_argument("--tolerance", type=float, default=2.0)
        parser.add_argument("--nsolv", type=int, default=12)
        parser.add_argument("--solvator-method", default="GFN2-xTB")
        parser.add_argument("--charge", type=int, default=0)
        parser.add_argument("--mult", type=int, default=1)
        parser.add_argument("--cores", type=int, default=1)
        parser.add_argument("--work-dir", default=None)
        args = parser.parse_args(argv)
        try:
            result = self.run(
                config={}, geometry=args.geometry,
                solvent_xyz=args.solvent_xyz, solvent=args.solvent,
                n_solvent=args.n_solvent, box_side=args.box_side,
                tolerance=args.tolerance, nsolv=args.nsolv,
                solvator_method=args.solvator_method,
                charge=args.charge, mult=args.mult, cores=args.cores,
                work_dir=args.work_dir)
        except Exception:  # noqa: BLE001
            return 1
        return 0 if result.get("ok") else 1


register(SolvationWorkflow())

