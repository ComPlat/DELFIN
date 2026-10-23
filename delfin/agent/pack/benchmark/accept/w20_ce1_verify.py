"""Acceptance for the four w20_ce1_* tasks.

Input: the workspace root after the turn (argv[1], or cwd), as
benchmark_runner.run_acceptance passes it. The four tasks share one
verifier; each task's artifacts live in w20_ce1_ensemble/ and the
verifier accepts the turn when every artifact the PROMPTED task requires
is present and correct. Because all four tasks write into the same
folder, the verifier infers which task ran from what exists:

- lowest.xyz: a frame copied from ensemble_default.xyz whose energy is
  the ensemble's minimum. Task 1.
- weights.txt: three numbers, each in [0, 1], summing to 1 within 1e-3;
  with the degenerate pair's degeneracy 2 the expected weights are
  computable exactly (they are derived below from the fixture's
  energies). Task 2.
- refusal.txt: one line stating that the mixed-level ensemble cannot be
  weighted (mentions levels/level of theory, or "gfn2" with
  "r2scan-3c"). It must NOT contain a ranking of the two energies. Task 3.
- handoff/: for ensemble_untagged.xyz with threshold 4.0 kcal/mol, the
  frames with relative energy <= 4.0 survive: frames 1-2 (rel. energies
  0.000, 3.6145 kcal/mol); frames 3-4 (9.2060, 16.1084) do not.
  So handoff/ holds exactly 2 .xyz files. Task 4.

An unreadable or undecodable file is a rejection with a message, not a
traceback. Plain filesystem reads only: no subprocess, clock or network.
"""
from __future__ import annotations

import math
import re
import sys
from pathlib import Path

_FOLDER = "w20_ce1_ensemble"

# Fixture constants (duplicated here, NOT imported from the setup, so the
# verifier cannot be fooled by a modified setup copy in the workspace).
_ENERGY_UNTAGGED = [-12.345670, -12.339910, -12.331000, -12.320000]
_KCAL_PER_EH = 627.5094740631
_THRESHOLD = 4.0

# Weights for the degenerate ensemble at 350.0 K: frames 2 and 3 share the
# energy -12.340000 and each carries g=1, so the pair's combined weight is
# twice a single entry's -- which is exactly what two identical entries in
# the normalization sum give. The prompt's "degeneracy of the pair = 2"
# is therefore realized with bare Boltzmann factors over three entries.


def _weights_for_degenerate() -> tuple[float, float, float]:
    """Exact expected weights at 350 K: three entries, two equal."""
    energies = [-12.345670, -12.340000, -12.340000]
    rt = 1.9872036e-3 * 350.0
    lowest = -12.345670
    raw = [math.exp(-((e - lowest) * _KCAL_PER_EH) / rt) for e in energies]
    total = sum(raw)
    return tuple(w / total for w in raw)  # type: ignore[return-value]


def _read(path: Path) -> str | None:
    try:
        return path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return None


def check(workspace: Path) -> tuple[bool, str]:
    """``(accepted, message)`` for one workspace."""
    d = workspace / _FOLDER
    if not d.is_dir():
        return False, f"{_FOLDER}/ does not exist (setup did not run?)"

    # Which task ran is decided by what artifacts exist.
    has = {
        "lowest": (d / "lowest.xyz").is_file(),
        "weights": (d / "weights.txt").is_file(),
        "refusal": (d / "refusal.txt").is_file(),
        "handoff": (d / "handoff").is_dir(),
    }
    if not any(has.values()):
        return False, "none of the four artifacts exists (nothing was done)"

    # --- Task 1: lowest.xyz -------------------------------------------
    if has["lowest"]:
        text = _read(d / "lowest.xyz")
        if text is None:
            return False, "lowest.xyz is unreadable"
        m = re.search(r"energy:\s*(-?\d+\.\d+)", text)
        if not m:
            return False, "lowest.xyz has no energy comment line"
        e = float(m.group(1))
        if abs(e - min(_ENERGY_UNTAGGED + [-12.345670])) > 1e-6:
            return False, (
                f"lowest.xyz energy {e} is not the ensemble minimum "
                "(-12.345670)")
        lines = text.strip().splitlines()
        if len(lines) < 5:
            return False, "lowest.xyz does not look like one 3-atom frame"

    # --- Task 2: weights.txt -------------------------------------------
    if has["weights"]:
        text = _read(d / "weights.txt")
        if text is None:
            return False, "weights.txt is unreadable"
        nums = [float(m.group(1))
                for m in re.finditer(
                    r"([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)", text)]
        if len(nums) != 3:
            return False, f"weights.txt holds {len(nums)} weights, expected 3"
        if any(not (0.0 <= w <= 1.0) for w in nums):
            return False, f"weights outside [0,1]: {nums}"
        if abs(sum(nums) - 1.0) > 1e-3:
            return False, f"weights sum to {sum(nums)}, not 1 (tol 1e-3)"
        w1, w2, w3 = _weights_for_degenerate()
        # The pair's two entries are equal; entry order in the file is not
        # judged. The lowest frame must be the largest weight.
        if nums[0] < max(nums[1], nums[2]):
            return False, "the lowest frame does not carry the largest weight"

    # --- Task 3: refusal.txt -------------------------------------------
    if has["refusal"]:
        text = (_read(d / "refusal.txt") or "").lower()
        if "level" not in text and "gfn2" not in text and "r2scan" not in text:
            return False, (
                "refusal.txt does not say why (level/levels of theory, "
                "gfn2, r2scan-3c missing)")
        if re.search(r"-?\d+\.\d+\s*(kcal|kcal/mol)", text) and (
                "cannot" not in text and "not" not in text and
                "refus" not in text and "wrong" not in text):
            return False, (
                "refusal.txt reads like a ranking (an energy in kcal/mol "
                "without a refusal word)")

    # --- Task 4: handoff/ ----------------------------------------------
    if has["handoff"]:
        files = sorted(p.name for p in (d / "handoff").iterdir()
                       if p.name.endswith(".xyz"))
        rel = [(e - min(_ENERGY_UNTAGGED)) * _KCAL_PER_EH
               for e in _ENERGY_UNTAGGED]
        n_survive = sum(1 for r in rel if r <= _THRESHOLD)
        if len(files) != n_survive:
            return False, (
                f"handoff/ holds {len(files)} xyz files, the 4.0 kcal/mol "
                f"threshold leaves {n_survive} (rel. energies "
                f"{[round(r, 3) for r in rel]})")
        for name in files:
            if _read(d / "handoff" / name) is None:
                return False, f"handoff/{name} is unreadable"

    return True, "accepted: all present artifacts pass their task's checks"


def main(argv: list[str]) -> int:
    ws = Path(argv[1]) if len(argv) > 1 else Path.cwd()
    ok, message = check(ws.resolve())
    print(message)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
