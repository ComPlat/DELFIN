"""What the patch changes, measured without writing the guarded file.

`delfin/agent/api_client.py` is the agent's own safety layer and the
self-modification guard refuses every write to it, `git apply` included.
Rule 13's route is a patch under scratch/ for a person to apply, so the
file on disk stays as the branch has it.

This probe answers the one question that leaves open: does the patch's
CONTENT actually make the read reachable? It rebinds the gate's allowlist
in memory (a frozenset cannot be mutated in place) and asks the gate the
same question the real run asked, before and after.

Run from the worktree root so `delfin` resolves to THIS tree:

    PYTHONPATH=. <venv>/bin/python scratch/probe_gate_effect.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import delfin
from delfin.agent import api_client as A

print("delfin imported from:", delfin.__file__)
# The tree this script lives in, not a hard-coded site path: a probe that
# measures the wrong tree is the failure mode, and naming the host's path
# to catch it would commit site data.
_ROOT = Path(__file__).resolve().parent.parent
assert Path(delfin.__file__).resolve().parent == _ROOT / "delfin", (
    f"wrong tree imported: {delfin.__file__} is not under {_ROOT}")

READS = ("eln_find", "eln_read", "eln_pending")
WRITES = ("eln_draft", "eln_attach", "eln_approve")

perms = A.KitToolPermissions(workspace=Path(tempfile.mkdtemp()), mode="default")


def ask(tool: str, args: dict) -> str:
    out = A._doc_executor._gate_mcp_tool(f"mcp__delfin-ops__{tool}", args, perms)
    return "ADMITTED" if out is None else "refused: " + str(out)[:90]


def line(tool: str, args: dict | None = None) -> None:
    print(f"  {tool:12s} default mode -> {ask(tool, args or {})}")


print("\nBEFORE the patch (as the branch stands):")
for t in READS:
    line(t, {"collection": "Synthesis 2026"} if t == "eln_find" else {})
for t in WRITES:
    line(t)

# The patch's content: the three reads join the allowlist.
A._MCP_READONLY_TOOL_BASES = A._MCP_READONLY_TOOL_BASES | set(READS)

print("\nAFTER the patch's content:")
for t in READS:
    line(t, {"collection": "Synthesis 2026"} if t == "eln_find" else {})
for t in WRITES:
    line(t)

ok_reads = all(
    A._doc_executor._gate_mcp_tool(f"mcp__delfin-ops__{t}",
                                   {"collection": "x"} if t == "eln_find" else {},
                                   perms) is None
    for t in READS
)
ok_writes = all(
    A._doc_executor._gate_mcp_tool(f"mcp__delfin-ops__{t}", {}, perms) is not None
    for t in WRITES
)
print(f"\nreads reachable in default mode: {ok_reads}")
print(f"writes still refused in default mode: {ok_writes}")
sys.exit(0 if (ok_reads and ok_writes) else 1)
