"""Workflow wrapper for model-free optimization campaigns.

Registered as ``campaign``.  A campaign optimizes a bounded candidate
space (e.g. substituents on a molecular anchor) for a HOMO-LUMO gap
target without a language model in the loop — see
:mod:`delfin.agent.campaign`.  The campaign folder carries
``CAMPAIGN.json`` (space, target, budget) and receives
``decisions.jsonl`` (every decision the loop made) and
``CAMPAIGN_summary.json`` (the result a caller reports).

The loop is LLM-free: the agent is woken only on failure, stagnation
or budget exhaustion, through the scheduler contract in
``delfin.agent.campaign.Campaign.report_to_scheduler``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List

from delfin.workflows.registry import register


class CampaignWorkflow:
    """CAMPAIGN: model-free Bayes optimization over a candidate space."""

    name = "campaign"
    description = (
        "Model-free optimization campaign over a candidate space "
        "(e.g. substituents on an anchor) for a HOMO-LUMO gap target"
    )

    # ------------------------------------------------------------------
    # spec handling
    # ------------------------------------------------------------------

    @staticmethod
    def _load_spec(path: Path) -> dict:
        """Load and validate ``CAMPAIGN.json``."""
        if not path.is_file():
            raise FileNotFoundError(
                f"campaign spec not found: {path}")
        spec = json.loads(path.read_text(encoding="utf-8"))
        if not spec.get("space"):
            raise ValueError("campaign spec needs a non-empty 'space'")
        return spec

    @staticmethod
    def _build(spec: dict, campaign_dir: Path) -> Any:
        """Build the Campaign object from a validated spec."""
        from delfin.agent.campaign import (
            Campaign,
            CampaignBudget,
            DecisionLog,
            TargetGap,
        )

        return Campaign(
            space=list(spec["space"]),
            target=TargetGap(spec.get("target_gap_eV")),
            budget=CampaignBudget(max_evaluations=int(
                spec.get("max_evaluations", 24))),
            log=DecisionLog(campaign_dir / "decisions.jsonl"),
            work_dir=campaign_dir,
            seed_used=int(spec.get("seed_used", 2)),
            stagnation_rounds=int(spec.get("stagnation_rounds", 3)),
            stagnation_tolerance=float(spec.get("stagnation_tolerance",
                                                0.01)),
            workspace=str(campaign_dir),
        )

    # ------------------------------------------------------------------
    # the workflow contract
    # ------------------------------------------------------------------

    def run(self, *, config: Dict[str, Any], **kwargs: Any) -> Any:
        """Run one campaign headless, from its folder.

        ``campaign_dir`` (or ``config['campaign_dir']``) is the folder
        holding ``CAMPAIGN.json``; the decision log and the summary are
        written there.  When ``config`` carries a scheduler object
        (``config['scheduler']``) the campaign's LLM-free wake-up is
        wired to it.
        """
        from delfin.agent.campaign import DecisionReason

        campaign_dir = Path(kwargs.get("campaign_dir")
                            or config.get("campaign_dir") or "").resolve()
        spec = self._load_spec(campaign_dir / "CAMPAIGN.json")
        camp = self._build(spec, campaign_dir)
        scheduler = config.get("scheduler")
        if scheduler is not None:
            camp.scheduler = scheduler

        strategy = str(spec.get("strategy", "ei"))
        camp.evaluate_candidates(strategy=strategy)

        # The loop is over: report the terminal condition to the
        # scheduler if one is wired and a wake reason exists.
        stopped_by = self._stopped_by(camp)
        camp.report_to_scheduler()

        best = camp.best_observations()
        best_entry = min(
            best.values(), key=lambda e: e["score"],
            default=None)
        summary = {
            "target_gap_eV": spec.get("target_gap_eV"),
            "n_evaluations": len(camp.observations),
            "budget": camp.budget.max_evaluations,
            "stopped_by": stopped_by,
            "best": best_entry,
            "observations": {
                str(k): v for k, v in camp.best_observations().items()},
            "decision_log": str(camp.log.path),
        }
        (campaign_dir / "CAMPAIGN_summary.json").write_text(
            json.dumps(summary, indent=2), encoding="utf-8")
        return summary

    @staticmethod
    def _stopped_by(camp: Any) -> str:
        from delfin.agent.campaign import DecisionReason

        if camp.budget.exhausted:
            return "budget"
        if camp._stagnated():
            return "stagnation"
        if camp.failed:
            return "failure"
        return "complete"

    def run_cli(self, argv: List[str]) -> int:
        """``delfin workflow campaign <campaign_dir>`` (see cli)."""
        import argparse

        parser = argparse.ArgumentParser(description=self.description)
        parser.add_argument("campaign_dir",
                            help="folder holding CAMPAIGN.json")
        args = parser.parse_args(argv)
        try:
            self.run(config={}, campaign_dir=args.campaign_dir)
        except Exception:  # noqa: BLE001
            return 1
        return 0


register(CampaignWorkflow())
