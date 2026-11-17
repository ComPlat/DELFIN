"""SW1 accept: a valid spec's dry run lists the right members and
starts nothing.

Usage: <script> <workspace>
"""

from __future__ import annotations

import sys
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    root = Path(argv[1]) / "swarm_fixtures"
    from delfin.agent.swarm.launch import plan_launch, render_dry_run
    from delfin.agent.swarm.spec import load_spec

    plan = plan_launch(load_spec(root / "valid_swarm.yaml"))
    names = {m["name"] for m in plan["members"]}
    expected = {"core-builder", "core-reviewer", "ui-builder"}
    if names != expected:
        print(f"members wrong: {sorted(names)} != {sorted(expected)}",
              file=sys.stderr)
        return 1
    text = render_dry_run(plan)
    for marker in ("Nothing started", "never bypass", "delfin/agent/swarm/**"):
        if marker not in text:
            print(f"dry-run text missing {marker!r}", file=sys.stderr)
            return 1
    print("accepted: dry run lists 3 members, scopes and starts nothing")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
