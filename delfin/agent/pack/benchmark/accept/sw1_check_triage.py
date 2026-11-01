"""SW1 accept: triage auto-answers a partner message but escalates a
push request.

Usage: <script> <workspace>
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    root = Path(argv[1]) / "swarm_fixtures"
    from delfin.agent.swarm.triage import TriageError, classify, triage_list
    from delfin.agent.swarm.merge import merge_plan  # noqa: F401 (wired)

    queue = json.loads((root / "triage_queue.json").read_text(
        encoding="utf-8"))
    listed = triage_list(queue)
    actions = {i["id"]: i["action"] for i in listed}
    if actions.get("msg-1") != "auto_answer":
        print(f"partner message not auto-answerable: {actions}",
              file=sys.stderr)
        return 1
    if actions.get("push-1") != "escalate":
        print(f"push request not escalated: {actions}", file=sys.stderr)
        return 1
    # An escalated item must refuse an auto-answer.
    try:
        from delfin.agent.swarm.triage import triage_answer
        push_item = next(i for i in queue if i["id"] == "push-1")
        triage_answer(push_item, classify(push_item), answer="yes")
    except TriageError:
        print("accepted: partner auto-answered, push escalated and "
              "refused an auto-answer")
        return 0
    print("escalated item accepted an auto-answer -- NOT accepted",
          file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
