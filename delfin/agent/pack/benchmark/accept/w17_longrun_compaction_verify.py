"""Accept: the forced-compaction task — the digest survives the compaction
and is COMPLETE. Compaction is attested from the session's own compaction
log (the real writer: compaction_log.record_compaction), not assumed; a
run that never compacted has not exercised what this task measures, and a
digest missing decisions is the loss the task exists to catch.
"""
from __future__ import annotations

import sys
from pathlib import Path

N_DECISIONS = 40


def _digest_entries(text: str) -> dict:
    out = {}
    for line in text.strip().splitlines():
        line = line.strip()
        if not line.lower().startswith("decision "):
            continue
        # 'Decision <n>: cutoff <x> eV'
        rest = line[len("Decision "):]
        n_part, _, val_part = rest.partition(":")
        try:
            n = int(n_part.strip())
        except ValueError:
            continue
        val = val_part.strip().lower().replace(" ", "")
        out[n] = val
    return out


def check(workspace: Path, n_compactions: int | None = None) -> tuple[bool, str]:
    """``(accepted, message)`` for one workspace.

    ``n_compactions`` is the count the caller read from the session's own
    compaction log (compaction_log.read_compactions for the run's
    session id). None means the log was not readable — unmeasured, not
    zero, and the rubric fails on it.
    """
    workspace = Path(workspace)
    d = workspace / "w17_longrun" / "lr_forced_compaction"
    digest = d / "digest.md"
    if not digest.is_file():
        return False, "digest.md is missing"
    entries = _digest_entries(digest.read_text(encoding="utf-8"))
    if len(entries) != N_DECISIONS:
        return False, (f"digest.md holds {len(entries)} decisions, "
                       f"expected {N_DECISIONS} — the compaction lost "
                       f"decisions (the loss this task exists to catch)")
    # Values re-checked against the notes that are still on disk.
    notes = (d / "notes.md").read_text(encoding="utf-8")
    for i in range(1, N_DECISIONS + 1):
        val = entries.get(i)
        if not val or "cutoff" not in val or "ev" not in val:
            return False, f"Decision {i} missing or malformed in digest.md"
        marker = f"Decision {i}: the cutoff for batch {i} is "
        if marker not in notes:
            return False, f"fixture drift: notes.md lacks decision {i}"
        want = notes.split(marker, 1)[1].split(" eV", 1)[0].strip()
        if want.replace(" ", "") not in val:
            return False, (f"Decision {i} value mismatch: digest says "
                           f"{val!r}, notes say {want} eV")
    if n_compactions is None:
        return False, ("compaction count not attested (no compaction log) "
                       "— the run never proved it survived a compaction")
    if n_compactions < 1:
        return False, ("the records attest no compaction in this run — "
                       "the task's precondition (a mid-run compaction) "
                       "was not exercised")
    return True, (f"accepted: digest complete ({N_DECISIONS} decisions, "
                  f"values re-checked against notes.md), compaction "
                  f"attested ({n_compactions} event(s) in the session's "
                  f"compaction log)")


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ok, msg = check(ws)
    print(msg)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
