"""Accept: the wait-for-job task — the collected result is right, and the
waiting is attested by the session's own records (longrun_metric), not by
a wall clock.

Two facts, both required:
1. result.json holds n=6 and the mean the fixture re-checked (42.19
   kJ/mol) — the slow job actually ran to completion and its output was
   collected.
2. The run records attest waiting: the caller hands the accept script
   the run file's longrun inputs, and idle_with_mail_s (or
   idle_with_due_event_s) is > 0 while the job was pending. If the
   session had no waiting recorded, the rubric fails — a run that
   crashed or skipped the wait is not a pass by luck of timing.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

EXPECTED_MEAN = 42.19


def check(workspace: Path, longrun_inputs: dict | None = None) -> tuple[bool, str]:
    """``(accepted, message)`` for one workspace.

    ``longrun_inputs`` is the longrun_metric.LongrunInputs-shaped dict the
    benchmark caller assembles from the run's own records (turn timestamps,
    waiting facts). None/empty means waiting could not be attested.
    """
    workspace = Path(workspace)
    res = workspace / "w17_longrun" / "lr_wait_for_job" / "result.json"
    if not res.is_file():
        return False, "result.json is missing — the slow job was not collected"
    try:
        data = json.loads(res.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001
        return False, f"result.json does not parse: {exc}"
    if int(data.get("n") or 0) != 6:
        return False, f"result.json n={data.get('n')!r}, expected 6"
    mean = data.get("mean_kj_per_mol")
    if not isinstance(mean, (int, float)) or abs(float(mean) - EXPECTED_MEAN) > 5e-3:
        return False, (f"result.json mean={mean!r}, expected {EXPECTED_MEAN} "
                       f"(±5e-3) — wrong value collected")
    # Fact 2: waiting attested by the records, not by timing luck.
    if not longrun_inputs:
        return False, ("no longrun inputs were handed to the accept script "
                       "— the waiting rubric is unmeasured")
    idle_mail = float(longrun_inputs.get("idle_with_mail_s") or 0.0)
    idle_due = float(longrun_inputs.get("idle_with_due_event_s") or 0.0)
    if idle_mail <= 0.0 and idle_due <= 0.0:
        return False, ("the records attest no idle-with-waiting window — "
                       "the wait was not measured from the records")
    return True, (f"accepted: result collected (n=6, mean {mean}), "
                  f"idle-with-waiting {max(idle_mail, idle_due):.1f} s "
                  f"attested by the records")


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    lr = None
    if len(sys.argv) > 2:
        lr = json.loads(Path(sys.argv[2]).read_text(encoding="utf-8"))
    ok, msg = check(ws, lr)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
