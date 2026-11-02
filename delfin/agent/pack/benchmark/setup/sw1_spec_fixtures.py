"""SW1 benchmark fixture: three swarm specs (overlapping, valid, triage
queue items) in the workspace.

Usage: <script> <workspace>
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

OVERLAPPING = """\
goal: Race test
packages:
  - name: a
    owned: ["delfin/dashboard/*.py"]
    task: t
  - name: b
    owned: ["delfin/dashboard/sw/**"]
    task: t
"""

VALID = """\
goal: Ship the swarm feature
packages:
  - name: core
    owned: ["delfin/agent/swarm/**"]
    task: Build the core package
    roles:
      builder: true
      reviewer: true
  - name: ui
    owned: ["delfin/agent/dashboard/**"]
    task: Build the dashboard view
"""

QUEUE = [
    {"id": "msg-1", "member": "core-builder", "from": "core-reviewer",
     "text": "review done, APPROVED with marks 5/5/5",
     "kind": "session_message",
     "roster": ["core-builder", "core-reviewer", "core-qs", "ui-builder"]},
    {"id": "push-1", "member": "core-builder", "from": "git",
     "text": "push to origin refused -- allow it?",
     "kind": "approval",
     "roster": ["core-builder", "core-reviewer", "core-qs", "ui-builder"]},
]


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    root = ws / "swarm_fixtures"
    if root.exists():
        print(f"refusing to overwrite {root}", file=sys.stderr)
        return 1
    root.mkdir(parents=True)
    (root / "overlapping_swarm.yaml").write_text(OVERLAPPING,
                                                 encoding="utf-8")
    (root / "valid_swarm.yaml").write_text(VALID, encoding="utf-8")
    (root / "triage_queue.json").write_text(
        json.dumps(QUEUE, indent=2), encoding="utf-8")

    # Check the preconditions rather than assume them.
    from delfin.agent.swarm.spec import SpecError, load_spec, validate_spec
    try:
        validate_spec(load_spec(root / "overlapping_swarm.yaml"))
    except SpecError as exc:
        if "overlap" not in str(exc):
            print(f"overlapping fixture refused for the WRONG reason: "
                  f"{exc}", file=sys.stderr)
            return 1
    else:
        print("overlapping fixture was NOT refused -- bad fixture",
              file=sys.stderr)
        return 1
    validate_spec(load_spec(root / "valid_swarm.yaml"))
    if len(QUEUE) != 2:
        print("queue fixture changed", file=sys.stderr)
        return 1
    print(f"prepared {root} -- overlapping + valid spec, triage queue")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
