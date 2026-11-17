"""SW1 accept: overlapping scopes refused with a clear error.

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
    from delfin.agent.swarm.spec import SpecError, load_spec, validate_spec
    try:
        validate_spec(load_spec(root / "overlapping_swarm.yaml"))
    except SpecError as exc:
        if "overlap" in str(exc):
            print("accepted: overlapping scopes refused with 'overlap'")
            return 0
        print(f"refused for the wrong reason: {exc}", file=sys.stderr)
        return 1
    print("overlapping spec was accepted -- NOT accepted", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
