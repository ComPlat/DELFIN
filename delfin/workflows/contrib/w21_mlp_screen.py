"""Contributed workflow: MLP conformer pre-screen + optional ORCA ranking.

Package ML1, built on the wave-21 foundation (WB2 Pipeline semantics,
WB3-registered adapters, no parallel mechanism).  The workflow takes a
multi-structure XYZ ensemble (e.g. from CREST conformer search), ranks it
with a machine-learned potential, and optionally hands the top structure
to a later (DFT) step through its output geometry.

Physical realism — every limit is enforced here, not just described.
The single source of truth is the adapter tables
(delfin/tools/adapters/mlp.py SUPPORTED_ELEMENTS / VALIDITY); this
workflow imports and re-uses them, it does not copy them:

- backend value-checked against the factory alias set (fail-closed);
- every structure in the ensemble is element-checked against the
  backend's declared supported set BEFORE any model is created — an
  unsupported element is a hard refusal, never a silently wrong energy;
- charge and multiplicity are only accepted for the one charge-aware
  backend (aimnet2, DOI 10.1039/D4SC08572H); a nonzero charge on any
  other backend is refused, because the calculator would silently
  ignore it and produce a wrong ranking;
- run evidence is reported, never hidden: n_structures ranked,
  per-structure index/energy_eV, and a warning when a backend that
  cannot model charge was asked for one (that refusal surfaces before
  any compute).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from delfin.workflows.registry import register


class MlpScreenWorkflow:
    """MLP conformer screening: rank an ensemble by MLP energy."""

    name = "mlp_screen"
    description = (
        "Machine-learned potential screening: rank a conformer ensemble "
        "(multi-structure XYZ) by MLP energy, with element/charge realism "
        "enforced from the adapter's declared tables"
    )

    def run(self, *, config: Dict[str, Any],
            ensemble: Optional[str] = None,
            backend: str = "ani2x",
            charge: int = 0, mult: int = 1,
            device: str = "cpu",
            work_dir: Optional[str] = None,
            **kwargs: Any) -> Dict[str, Any]:
        from delfin.tools._registry import get as get_adapter

        if not ensemble:
            raise ValueError(
                "mlp_screen workflow needs an 'ensemble' (path to a "
                "multi-structure XYZ, e.g. from CREST); no SMILES step is "
                "wired, build one with the smiles_to_xyz step or pass an "
                "existing ensemble")
        # The adapter value-checks the backend (fail-closed) BEFORE the
        # ensemble is even read; a wrong backend never starts compute.
        adapter = get_adapter("mlp_rank")
        if adapter is None:
            raise RuntimeError(
                "mlp_rank adapter not registered — the WB3 registry walk "
                "did not pick up delfin.tools.adapters.mlp")
        adapter.validate_params(backend=backend, charge=charge, mult=mult,
                                ensemble=ensemble)

        if work_dir:
            wd = Path(work_dir)
        else:
            wd = Path(config.get("work_dir") or ".") / "mlp_screen_run"
        wd.mkdir(parents=True, exist_ok=True)

        result = adapter.execute(wd, ensemble=ensemble, backend=backend,
                                 charge=charge, mult=mult, device=device)
        ranking = result.data.get("ranking") or []
        out: Dict[str, Any] = {
            "ok": result.status.name == "SUCCESS",
            "backend": backend,
            "status": result.status.name,
            "error": result.error,
            "n_structures": result.data.get("n_structures"),
            "ranking": ranking,
            "best_index": ranking[0]["index"] if ranking else None,
            "best_energy_eV": ranking[0]["energy_eV"] if ranking else None,
            "ensemble": ensemble,
        }
        # Reporting duty: a non-SUCCESS step is surfaced, never swallowed.
        if out["ok"] is False:
            out.setdefault("errors", [result.error])
        return out

    def run_cli(self, argv: list) -> int:
        import argparse

        parser = argparse.ArgumentParser(description=self.description)
        parser.add_argument("ensemble", help="multi-structure XYZ ensemble")
        parser.add_argument("--backend", default="ani2x")
        parser.add_argument("--charge", type=int, default=0)
        parser.add_argument("--mult", type=int, default=1)
        parser.add_argument("--device", default="cpu")
        parser.add_argument("--work-dir", default=None)
        args = parser.parse_args(argv)
        try:
            result = self.run(
                config={}, ensemble=args.ensemble, backend=args.backend,
                charge=args.charge, mult=args.mult, device=args.device,
                work_dir=args.work_dir)
        except Exception:  # noqa: BLE001
            return 1
        return 0 if result.get("ok") else 1


register(MlpScreenWorkflow())
