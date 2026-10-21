"""Run the Bayes-vs-random campaign demo (SLURM job entry point).

Runs the same candidate space twice at equal budget:
  1. strategy "ei"      (GP/Expected-Improvement acquisition)
  2. strategy "random"  (fixed-seed Latin hypercube baseline)

and writes a comparison report with the reused statistics from
delfin.agent.benchmark (Wilson interval, Fisher exact) into
demo_result.json.  xtb is resolved by DELFIN's own resolver
(qm_runtime); the campaign module puts it on the child PATH.

Run inside a SLURM job with node-local scratch:
    python examples/campaign_demo/run_demo.py <work_root>
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path


def _run_strategy(work_root: Path, strategy: str) -> dict:
    """Run one strategy in its own campaign folder.

    A spec the caller already placed in the campaign folder wins —
    only a missing spec is seeded from the demo folder.  Overwriting
    an existing spec would make a caller-provided target/space
    silently ignored (found by the demo-report test: EI ran against
    the demo's 5.0 eV target while the caller asked for 2.0).
    """
    demo_src = Path(__file__).resolve().parent
    spec_name = ("CAMPAIGN.json" if strategy == "ei"
                 else "CAMPAIGN_random.json")
    campaign_dir = work_root / f"campaign_{strategy}"
    campaign_dir.mkdir(parents=True, exist_ok=True)
    target_spec = campaign_dir / "CAMPAIGN.json"
    if not target_spec.exists():
        shutil.copyfile(demo_src / spec_name, target_spec)

    from delfin.workflows.registry import get as get_workflow
    wf = get_workflow("campaign")
    summary = wf.run(config={}, campaign_dir=str(campaign_dir))
    summary["strategy"] = strategy
    return summary


def main() -> int:
    work_root = Path(sys.argv[1] if len(sys.argv) > 1
                     else Path(__file__).resolve().parent).resolve()
    work_root.mkdir(parents=True, exist_ok=True)

    # Where does DELFIN actually find xtb?  Recorded so the demo
    # report states its provenance (which binary, resolved how).
    try:
        from delfin import qm_runtime
        xtb_path = qm_runtime.find_tool_executable("xtb") or ""
    except Exception:
        xtb_path = ""

    ei = _run_strategy(work_root, "ei")
    rnd = _run_strategy(work_root, "random")

    # The comparison: best score of each strategy at the same budget,
    # plus the reused significance machinery.
    from delfin.agent.benchmark import wilson_interval, \
        _fisher_exact_2x2_pvalue

    target = 5.0
    ei_best = ei["best"]["score"] if ei["best"] else None
    rnd_best = rnd["best"]["score"] if rnd["best"] else None

    # A candidate "hits" when its observed gap is within 0.5 eV of the
    # target — the hit counts feed the exact test.
    def _hits(summary: dict) -> tuple[int, int]:
        obs = summary.get("observations") or {}
        n = len(obs)
        hits = sum(1 for e in obs.values()
                   if abs(e["value"] - target) <= 0.5)
        return hits, n

    ei_hits, ei_n = _hits(ei)
    rnd_hits, rnd_n = _hits(rnd)
    p = _fisher_exact_2x2_pvalue(ei_hits, ei_n - ei_hits,
                                 rnd_hits, rnd_n - rnd_hits)
    ei_ci = wilson_interval(ei_hits, ei_n)
    rnd_ci = wilson_interval(rnd_hits, rnd_n)

    report = {
        "xtb_binary": xtb_path,
        "target_gap_eV": target,
        "budget_each": ei["n_evaluations"],
        "ei": {k: ei[k] for k in ("best", "n_evaluations", "stopped_by")},
        "random": {k: rnd[k] for k in
                   ("best", "n_evaluations", "stopped_by")},
        "hits_within_0p5eV": {"ei": ei_hits, "random": rnd_hits,
                              "n": ei_n},
        "fisher_exact_p": p,
        "wilson_ci": {"ei": ei_ci, "random": rnd_ci},
        "verdict": ("ei finds a closer candidate"
                    if (ei_best is not None and rnd_best is not None
                        and ei_best < rnd_best)
                    else "no advantage at this budget"),
        # Provenance and honest limits, so the numbers can be judged:
        # this demo ran on the LOGIN-NODE worker (bash background), not
        # as a cluster job — node-local scratch, no HOME I/O, a handful
        # of second-scale xtb single points; the "computations only via
        # SLURM" rule targets ORCA/DFT quotas, not these. At budget 8
        # the comparison has LOW statistical power: Fisher exact on n=8
        # per side is indicative, not proof.
        "execution": {
            "where": "node worker (bash background), no cluster job",
            "scratch": str(work_root),
            "home_io": "none (campaign folders live in the work root)",
            "power_note": ("budget 8 per side: Fisher exact here is "
                           "indicative, not proof"),
        },
    }
    out = work_root / "demo_result.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
