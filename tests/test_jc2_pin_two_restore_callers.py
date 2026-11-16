"""QS ruling 9, item 2 (nacht-s8, 2026-11-15): PIN that
delfin.agent.session_store.consume_crash_recovery_note has exactly the
two real restore callers — agent/engine.py (service restore during
session resume, engine.py:9046) and agent/cli.py (headless resume,
cli.py:1514) — and NO per-turn or per-view restore path. If a future
commit adds a per-turn restore, this pin fails before the chain offer
can repeat per turn.
"""
import re
from pathlib import Path


def test_consume_note_has_exactly_the_two_restore_callers():
    callers = []
    # engine.py (service path: restore during session resume) —
    # engine.py:9046 is the real service-restore call site.
    engine = Path("delfin/agent/engine.py")
    esrc = engine.read_text(encoding="utf-8")
    assert re.search(r"consume_crash_recovery_note", esrc), (
        "agent/engine.py no longer calls consume_crash_recovery_note — "
        "the service resume path lost the chain offer")
    # cli.py (headless resume).
    cli = Path("delfin/agent/cli.py")
    src = cli.read_text(encoding="utf-8")
    assert re.search(r"consume_crash_recovery_note", src), (
        "cli.py no longer calls consume_crash_recovery_note — the "
        "headless resume path lost the chain offer")
    callers = [str(engine), str(cli)]

    # No OTHER delfin module may call it (the pin: a new caller needs an
    # explicit re-review — a per-turn or per-view restore would re-open
    # the offer path every turn).
    roots = [Path("delfin")]
    allowed = {"agent/engine.py", "agent/cli.py",
               "agent/session_store.py"}  # the definition itself
    offenders = []
    for root in roots:
        for path in root.rglob("*.py"):
            rel = path.relative_to(root).as_posix()
            if rel in allowed:
                continue
            text = path.read_text(encoding="utf-8", errors="replace")
            if "consume_crash_recovery_note" in text:
                offenders.append(f"delfin/{rel}")
    assert not offenders, (
        "new consume_crash_recovery_note caller(s) outside the two "
        f"documented restore sites: {offenders} — a per-turn or per-view "
        "restore would repeat the chain offer every turn; re-review "
        "required")
