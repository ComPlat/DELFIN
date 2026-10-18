"""Accept: the restart-midway task — the restart continued the checkpoint
instead of starting over. Three facts, all required:

1. energy.txt holds all four molecules, alpha/beta with the values the
   checkpoint recorded (NOT recomputed — the recorded values are the
   only sanctioned source) and gamma/delta equal to the mean of their
   two columns in energies.csv.
2. Work-redone metric: the longrun records attest NO redo of the
   checkpointed part (longrun_metric's redone_after_compaction class of
   judgement, applied by the caller to the restart boundary). Redoing
   alpha/beta is exactly the failure this task exists to measure.
3. The run is a real restart: the caller attests a restart boundary was
   observed (a second session/run on the same records). Unmeasured
   fails — silently skipping the check would make every run pass.
"""
from __future__ import annotations

import sys
from pathlib import Path

_CHECKPOINT_VALUES = {"alpha": "2.10", "beta": "5.80"}
_MEANS = {"gamma": "1.25", "delta": "3.10"}


def _parse(text: str) -> dict:
    out = {}
    for line in text.strip().splitlines():
        line = line.strip()
        if not line or ":" not in line:
            continue
        name, _, val = line.partition(":")
        out[name.strip().lower()] = val.strip().lower().replace(" ", "")
    return out


def check(workspace: Path, *, redone_after_restart: bool | None = None,
          restart_attested: bool | None = None) -> tuple[bool, str]:
    workspace = Path(workspace)
    d = workspace / "w17_longrun" / "lr_restart_midway"
    energy = d / "energy.txt"
    if not energy.is_file():
        return False, "energy.txt is missing"
    parsed = _parse(energy.read_text(encoding="utf-8"))
    for name, want in {**_CHECKPOINT_VALUES, **_MEANS}.items():
        got = parsed.get(name, "")
        if want not in got:
            return False, (f"energy.txt {name}: {got!r}, expected "
                           f"{want} (kJ/mol) — checkpoint not continued "
                           f"correctly")
    if redone_after_restart is None:
        return False, ("work-redone not attested (no records handed over) "
                       "— unmeasured, not a pass")
    if redone_after_restart:
        return False, ("the records attest redone work after the restart "
                       "boundary (alpha/beta recomputed) — exactly the "
                       "failure this task exists to measure")
    if restart_attested is not True:
        return False, ("a restart boundary is not attested in the records "
                       "— the run never restarted, so it did not exercise "
                       "checkpoint continuation")
    return True, ("accepted: all four values correct from the checkpoint "
                  "(alpha 2.10, beta 5.80 recorded; gamma 1.25, delta "
                  "3.10 computed), no redone work, restart attested")


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
