#!/usr/bin/env python3
"""Acceptance for the wave-25 CAS1 benchmark tasks (package CAS1).

The agent's deliverable is proposal.json (the CAS(n,m) active space)
plus answer.txt (the justification), both graded by RE-RUNNING
delfin.casscf_active_space's checks -- the same enforcement the product
exposes -- never by trusting the prose:

* proposal.json parses and validates through ActiveSpaceProposal
  (justification with chemistry, named orbital classes, nel <= 2*norb,
  norb <= 16).
* The trap task (cas1_trap_oversized) must REFUSE: Nactel/Nactorb are
  not real %casscf keywords (manual 3.13.3) and 24 orbitals exceed
  CASSCF -- answer.txt must name RAS, ICE or DMRG and the real keys
  nel/norb.

usage: cas1_verify.py <workspace> [<task_id>]
Without a task id (the runner's form, benchmark_runner.py:1349) every
folder in the workspace that holds a proposal.json, or the workspace
root itself, is graded; an answer-only deliverable grades as the trap.
Exit 0 = accepted, 1 = failed, 2 = usage error.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


def _fail(msg: str) -> None:
    print(f"FAIL: {msg}", file=sys.stderr)
    sys.exit(1)


def _deliverable_dirs(ws: Path) -> list[Path]:
    """Folders in *ws* that hold the agent's deliverables, newest-last.

    benchmark_runner.run_acceptance (benchmark_runner.py:1349) launches
    this script with the workspace as the ONLY argument -- it never
    passes the task id, and no setup script writes a task-id marker.
    The real source the runner provides is the workspace layout itself:
    the setup pre-creates the family folders (or leaves the root plain,
    as here), and the attempt writes only where the agent worked, so at
    verify time the folders holding a proposal.json / answer.txt are the
    attempt's.  Checked shallow (root) and one level deep, sorted for a
    deterministic grade.
    """
    dirs: list[Path] = []
    for p in sorted(ws.rglob("proposal.json")) + sorted(ws.rglob("answer.txt")):
        if p.parent not in dirs:
            dirs.append(p.parent)
    if not dirs:
        dirs = [ws]
    return dirs


def main(argv: list[str]) -> None:
    if len(argv) == 3:
        ws, task_id = Path(argv[1]), argv[2]
        d = ws / task_id
        if not d.is_dir():
            _fail(f"{d} missing")
        _grade_one(d, task_id)
        return
    if len(argv) == 2:
        ws = Path(argv[1])
        if not ws.is_dir():
            _fail(f"workspace {ws} does not exist")
        dirs = _deliverable_dirs(ws)
        if len(dirs) == 1 and dirs[0] == ws and not (ws / "answer.txt").is_file():
            _fail("no proposal.json and no answer.txt found in the "
                  "workspace -- nothing to grade")
        for d in dirs:
            _grade_one(d, d.name or "workspace")
        return
    print("usage: cas1_verify.py <workspace> [<task_id>]",
          file=sys.stderr)
    sys.exit(2)


def _grade_one(d: Path, task_id: str) -> None:
    """Grade one deliverable folder: the trap answer, or a proposal."""
    answer_p = d / "answer.txt"
    prop_p = d / "proposal.json"
    answer = (answer_p.read_text(encoding="utf-8")
              if answer_p.is_file() else "")

    if not prop_p.is_file():
        # Trap-shaped deliverable: the oversized space must be REFUSED.
        if not answer:
            _fail(f"{answer_p} missing — no answer to grade")
        if "Nactel" in answer and "not" not in answer.lower():
            _fail("the trap keywords were passed through as real")
        lowered = answer.lower()
        if not any(alt in lowered for alt in ("ras", "ice", "dmrg")):
            _fail("no alternative for the oversized space named "
                  "(RAS / ICE / DMRG expected)")
        if "nel" not in lowered or "norb" not in lowered:
            _fail("the real %casscf keywords nel/norb were not stated")
        print(f"OK: {task_id}: trap refused, alternatives and real "
              "keywords named")
        return

    # Proposal tasks: grade the structured deliverable.
    try:
        raw = json.loads(prop_p.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        _fail(f"{prop_p} does not parse: {exc}")
    try:
        from delfin.casscf_active_space import ActiveSpaceProposal
        proposal = ActiveSpaceProposal(
            n_electrons=int(raw["n_electrons"]),
            n_orbitals=int(raw["n_orbitals"]),
            justification=str(raw.get("justification", "")),
            named_orbital_classes=tuple(
                raw.get("named_orbital_classes") or ()),
        )
        proposal.validate()
    except (KeyError, TypeError, ValueError) as exc:
        _fail(f"the proposal does not survive the module's checks: {exc}")
    print(f"OK: {task_id}: proposal CAS({raw['n_electrons']},"
          f"{raw['n_orbitals']}) validated")


if __name__ == "__main__":
    main(sys.argv)
