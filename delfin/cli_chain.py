"""The ``delfin chain`` subcommand (JC1 target 5).

Read-only status and an EXPLICIT resume — a chain is never resumed on
its own without the user asking for it.

    delfin chain status <work-dir>   print the per-step state, job id and
                                     evidence; exit non-zero when there
                                     is no chain state
    delfin chain resume <work-dir>   resume the chain from its persisted
                                     state, deciding from facts (JC1
                                     target 2) and never overwriting an
                                     old attempt (target 3)

The state file lives next to the chain's work dir
(``<work-dir>/chain_state.json``), written by :mod:`delfin.tools.pipeline`
(one mechanism, JC1 target 4).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from delfin.common.logging import get_logger

logger = get_logger(__name__)

STATE_FILENAME = "chain_state.json"


def _load_chain(work_dir: str):
    """Load the chain state for *work_dir*, or None when there is none."""
    from delfin.tools import chain_state as _cs

    state_path = Path(work_dir) / STATE_FILENAME
    if not state_path.is_file():
        return None
    return _cs.ChainState.load(state_path)


def _cmd_status(args) -> int:
    """Print the chain's per-step facts (read-only)."""
    chain = _load_chain(args.work_dir)
    if chain is None:
        print(
            f"no chain state file in {args.work_dir} "
            f"(expected {STATE_FILENAME} next to the chain's work dir)"
        )
        return 1
    print(f"chain '{chain.name}' — work dir: {args.work_dir}")
    for step in chain.steps:
        lines = [
            f"  step '{step.name}' ({step.tool}): {step.state}",
        ]
        if step.job_id is not None:
            lines.append(f"    job: {step.job_id}")
        if step.work_dir:
            lines.append(f"    work dir: {step.work_dir} "
                         f"(attempt {step.attempt})")
        if step.state == "finished" and step.evidence:
            lines.append(f"    evidence: {sorted(step.evidence)}")
        if step.state == "failed":
            if step.error_kind:
                lines.append(f"    error kind: {step.error_kind}")
            lines.append(f"    error: {step.error}")
        for line in lines:
            print(line)
    return 0


def _cmd_resume(args) -> int:
    """Resume the chain EXPLICITLY (target 5: never on its own)."""
    chain = _load_chain(args.work_dir)
    if chain is None:
        print(
            f"no chain state file in {args.work_dir} "
            f"(expected {STATE_FILENAME} next to the chain's work dir)"
        )
        return 1

    from delfin.tools import pipeline as pipeline_mod
    from delfin.tools.pipeline import Pipeline

    # Rebuild the pipeline from the chain's own step list so the resume
    # runs exactly the steps the state knows (the same step names and
    # tools the chain recorded).  The pipeline reads the persisted state
    # itself (one mechanism, target 4) via resume_chain=True.
    pipe = Pipeline(name=chain.name)
    for step in chain.steps:
        pipe.add(step.tool, label=step.name)
    result = pipe.run(
        cores=args.cores,
        work_dir=Path(args.work_dir),
        resume_chain=True,
        chain_state_obj=chain,
    )
    if not result.ok:
        failed = [
            r for r in result.results if not r.ok
        ]
        for r in failed:
            logger.error(
                "chain '%s' step '%s' failed: %s",
                chain.name, r.step_name, r.error or "failed",
            )
        return 1
    print(f"chain '{chain.name}' resumed and finished "
          f"({len(result.results)} step results)")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="delfin chain",
        description="Status and explicit resume of a crash-safe "
                    "calculation chain (JC1).",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_status = sub.add_parser(
        "status", help="print the chain's per-step state (read-only)"
    )
    p_status.add_argument("work_dir", help="the chain's work dir")
    p_status.set_defaults(func=_cmd_status)

    p_resume = sub.add_parser(
        "resume", help="resume the chain (explicit; never automatic)"
    )
    p_resume.add_argument("work_dir", help="the chain's work dir")
    p_resume.add_argument(
        "--cores", type=int, default=1,
        help="cores for the resumed steps (default 1)",
    )
    p_resume.set_defaults(func=_cmd_resume)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
