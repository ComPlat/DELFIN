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
import os
import sys
import time
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
        # JC3: an open step says whether its process still runs.
        from delfin.tools.chain_state import runner_liveness

        _alive, why = runner_liveness(step)
        if why:
            lines.append(f"    process: {why}")
        # QS/reviewer note: the persisted params are part of the chain's
        # facts — status shows them so a user sees exactly what a resume
        # would rebuild (and what the inputs hash covers).
        if step.params:
            lines.append(f"    params: {sorted(step.params.items())}")
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


def _cmd_list(args) -> int:
    """List the indexed chains with their CURRENT state (read-only).

    The index is a pointer, never the truth: each entry's state is
    re-read from that chain's own chain_state.json (JC2 target 2).  A
    vanished chain is still listed (as gone) so the user sees what the
    index remembers; an empty index prints one line, not nothing.
    """
    from delfin.tools import chain_index as _ci
    from delfin.tools import chain_state as _cs

    entries = _ci.load()
    # Live run 2 (operator, 2026-11-15): a chain started by an agent
    # through bash is never in the index (HOME is read-only inside the
    # sandbox) — `chain list` would hide it.  The CURRENT directory is
    # scanned with the same bounded rule as the offer: root + direct
    # subdirectories only, no tree walk; deduped against the index by
    # resolved work dir; marked as not indexed.
    from delfin.agent.session_store import _scan_workspace_chains

    seen_dirs = {str(Path(e.get("work_dir") or "").resolve())
                 for e in entries if e.get("work_dir")}
    for wd, ch in _scan_workspace_chains(os.getcwd()):
        if str(Path(wd).resolve()) in seen_dirs:
            continue
        entries.append({"name": ch.name, "work_dir": str(wd),
                        "found_in_this_directory": True})
    if not entries:
        print("no chains recorded in the chain index")
        return 0
    for e in entries:
        name = e.get("name", "?")
        wd = e.get("work_dir", "?")
        created = e.get("created_at")
        created_s = (
            time.strftime("%Y-%m-%d %H:%M", time.localtime(created))
            if isinstance(created, (int, float)) else "?"
        )
        state_path = Path(wd) / STATE_FILENAME if wd != "?" else None
        if state_path is None or not state_path.is_file():
            state = "gone (no chain state file)"
            steps_line = ""
        else:
            try:
                chain = _cs.ChainState.load(state_path)
            except Exception as exc:  # noqa: BLE001 - report, don't die
                state = f"unreadable ({exc})"
                steps_line = ""
            else:
                states = [s.state for s in chain.steps]
                if chain.steps and all(s == "finished" for s in states):
                    state = "finished"
                else:
                    state = " / ".join(states) if states else "no steps"
                steps_line = f" — steps: {', '.join(states)}"
        print(f"chain '{name}' — work dir: {wd} — created: {created_s} "
              f"— state: {state}{steps_line}"
              + (" — found in this directory (not indexed)"
                 if e.get("found_in_this_directory") else ""))
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

    # JC3: a step whose process verifiably still runs is not resumed — a
    # second run of the same step would start next to the first one.
    from delfin.tools.chain_state import RUNNING, runner_liveness

    running = [(st.name, why) for st in chain.steps
               for alive, why in [runner_liveness(st)] if alive == RUNNING]
    if running:
        for name, why in running:
            print(f"chain '{chain.name}': step '{name}' is {why}")
        print(f"chain '{chain.name}' NOT resumed: it is still running; "
              "wait for it to end, or stop it first")
        return 1

    # L2: the state file as it was found.  A resume that fails (a bad
    # rebuild, a crash inside the pipeline) must leave it EXACTLY as it
    # was — the crashed run's submitted/attempt facts are the only
    # record that the first attempt ever happened.
    state_path = Path(args.work_dir) / STATE_FILENAME

    from delfin.tools import pipeline as pipeline_mod
    from delfin.tools.pipeline import Pipeline

    # Live-run plan Job B: a resume never silently changes PAL.  With no
    # explicit --cores flag the resume uses the cores the chain's steps
    # were submitted with (recorded per step); the flag wins when given.
    recorded_cores = {s.cores for s in chain.steps if s.cores}
    if getattr(args, "cores", None) is not None:
        cores = args.cores
    elif recorded_cores:
        if len(recorded_cores) > 1:
            # Mixed PALs across steps: the largest is the safe default —
            # no step is silently downgraded mid-chain.
            cores = max(recorded_cores)
        else:
            cores = recorded_cores.pop()
    else:
        cores = 1

    # Rebuild the pipeline from the chain's OWN step specs (L1: the
    # persisted params — charge, method, basis, mult, ... — so the
    # rebuilt run satisfies the real adapter's validate_params, not
    # only a fake that ignores kwargs).
    pipe = Pipeline(name=chain.name)
    for step in chain.steps:
        pipe.add(step.tool, label=step.name, **step.params)
    # L9: the per-step state as it was found.  "Leave the state as it
    # was" is right for a crash BEFORE anything new is proven — but a
    # step that FINISHED with evidence during this resume holds a fact,
    # and that fact must survive a later failure in the same resume.
    # Restore is therefore per step: a step whose state changed AND
    # carries no evidence is reset to its pre-resume state; a finished
    # step with evidence keeps its fact.
    pre_state = {s.name: (s.state, s.evidence, s.work_dir, s.attempt,
                          list(s.attempts), s.params, s.inputs_hash)
                 for s in chain.steps}
    try:
        result = pipe.run(
            cores=cores,
            work_dir=Path(args.work_dir),
            resume_chain=True,
            chain_state_obj=chain,
        )
    except Exception as exc:
        # L9: per-step restore (see above); NOTE: one clear line,
        # never a traceback at the user.
        for s in chain.steps:
            before = pre_state.get(s.name)
            if before is None:
                continue
            b_state, b_evidence, b_dir, b_attempt, b_attempts, b_params, \
                b_hash = before
            if s.state == b_state:
                continue
            if s.state == "finished" and s.evidence:
                # A verified finish is a fact: keep it.
                chain.save()
                continue
            s.state, s.evidence, s.work_dir = b_state, b_evidence, b_dir
            s.attempt, s.attempts = b_attempt, b_attempts
            s.params, s.inputs_hash = b_params, b_hash
        chain.save()
        logger.error("chain '%s' resume failed: %s", chain.name, exc)
        print(f"chain '{chain.name}' resume failed: {exc}")
        return 1
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

    p_list = sub.add_parser(
        "list", help="list indexed chains with their current state (read-only)"
    )
    p_list.set_defaults(func=_cmd_list)

    p_resume = sub.add_parser(
        "resume", help="resume the chain (explicit; never automatic)"
    )
    p_resume.add_argument("work_dir", help="the chain's work dir")
    p_resume.add_argument(
        "--cores", type=int, default=None,
        help="cores for the resumed steps (default: the cores the chain "
             "was submitted with; never silently different)",
    )
    p_resume.set_defaults(func=_cmd_resume)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
