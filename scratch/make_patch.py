"""Generate the guarded-file patch: add the Chemotion ELN reads to the
gate's read-only allowlist.

The edit itself was refused by the self-modification guard (api_client.py
is the agent's own safety layer). Rule 13's route is a patch under
scratch/, applied by a person, never a bash detour around the guard.

Mechanical on purpose: difflib emits the hunk and its counts, so a
hand-counted line number cannot drift.
"""

from __future__ import annotations

import difflib
from pathlib import Path

TARGET = Path("delfin/agent/api_client.py")

OLD = '''    # --- DELFIN ops server: literature reading ---
    "check_orca_manual_indexed", "read_pdf", "search_pdf_local",
    "extract_pdf_section",
})'''

NEW = '''    # --- DELFIN ops server: literature reading ---
    "check_orca_manual_indexed", "read_pdf", "search_pdf_local",
    "extract_pdf_section",
    # --- DELFIN ops server: Chemotion ELN reads ---
    # Reading a collection, one element, or the drafts waiting for a
    # person. They send GET requests and change nothing. Leaving them off
    # this list had exactly the shape its own comment above warns about,
    # and it was found in a real run rather than by a test: ``eln_find``
    # asked for the samples of a collection was refused in default mode,
    # so a chemist's read question could not be answered at all.
    #
    # The writes stay off this list on purpose. ``eln_draft``,
    # ``eln_attach`` and ``eln_approve`` each need one person's approval,
    # and ``eln_approve`` exists so the agent cannot give that approval.
    "eln_find", "eln_read", "eln_pending",
})'''


def main() -> int:
    text = TARGET.read_text(encoding="utf-8")
    if text.count(OLD) != 1:
        raise SystemExit(f"anchor occurs {text.count(OLD)} times, expected 1")
    updated = text.replace(OLD, NEW)
    diff = difflib.unified_diff(
        text.splitlines(keepends=True),
        updated.splitlines(keepends=True),
        fromfile="a/delfin/agent/api_client.py",
        tofile="b/delfin/agent/api_client.py",
        n=3,
    )
    out = Path("scratch/eln-reads-in-the-gate-allowlist.patch")
    out.write_text("".join(diff), encoding="utf-8")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
