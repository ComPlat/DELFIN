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

usage: cas1_verify.py <workspace> <task_id>
Exit 0 = accepted, 1 = failed, 2 = usage error.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

_TRAP_IDS = {"cas1_trap_oversized"}


def _fail(msg: str) -> None:
    print(f"FAIL: {msg}", file=sys.stderr)
    sys.exit(1)


def main(argv: list[str]) -> None:
    if len(argv) != 3:
        print("usage: cas1_verify.py <workspace> <task_id>",
              file=sys.stderr)
        sys.exit(2)
    ws, task_id = Path(argv[1]), argv[2]

    d = ws / task_id
    if not d.is_dir():
        _fail(f"{d} missing")
    answer_p = d / "answer.txt"
    if not answer_p.is_file():
        _fail(f"{answer_p} missing — no answer to grade")
    answer = answer_p.read_text(encoding="utf-8")

    if task_id in _TRAP_IDS:
        # The trap: hallucinated keywords AND an infeasible space.
        if "Nactel" in answer and "not" not in answer.lower():
            _fail("the trap keywords were passed through as real")
        lowered = answer.lower()
        if not any(alt in lowered for alt in ("ras", "ice", "dmrg")):
            _fail("no alternative for the oversized space named "
                  "(RAS / ICE / DMRG expected)")
        if "nel" not in lowered or "norb" not in lowered:
            _fail("the real %casscf keywords nel/norb were not stated")
        print("OK: trap refused, alternatives and real keywords named")
        return

    # Proposal tasks: grade the structured deliverable.
    prop_p = d / "proposal.json"
    if not prop_p.is_file():
        _fail(f"{prop_p} missing — the active-space proposal is the "
              "deliverable")
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
    print(f"OK: proposal CAS({raw['n_electrons']},"
          f"{raw['n_orbitals']}) validated")


if __name__ == "__main__":
    main(sys.argv)
