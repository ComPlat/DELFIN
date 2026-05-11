#!/usr/bin/env python3
"""Harder live test: does the KIT model use the NEW tools correctly?

The simple "write + run pytest" test passes easily. This script
provokes the harder behaviours: pre-existing buggy code that needs
diagnosis + fix, a long-running command that should land in
bash_background (not blocking foreground bash), and ideally
task_create/update for tracking the multi-step work.

Sets up a sandboxed project with an off-by-one bug already in the
code and a failing test, hands the agent a single instruction, and
checks both the resulting code state AND the tool sequence.
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from delfin.agent import audit_log as al
from delfin.agent.api_client import KitToolPermissions, OpenAIClient


MODEL = sys.argv[1] if len(sys.argv) > 1 else "azure.gpt-5.1"


def main() -> int:
    api_key = os.environ.get("KIT_TOOLBOX_API_KEY", "")
    if not api_key:
        print("FAIL: KIT_TOOLBOX_API_KEY not set", file=sys.stderr)
        return 2

    tmp = Path(tempfile.mkdtemp(prefix="kit_adv_"))
    workspace = tmp / "workspace"
    workspace.mkdir()
    audit_log = tmp / "audit.log"
    al._default_log_path = lambda: audit_log

    # Pre-populate a project with a buggy function and a test that
    # currently FAILS (off-by-one). The agent has to diagnose + fix.
    src = workspace / "stats.py"
    src.write_text(
        "def average(values):\n"
        "    \"\"\"Return the arithmetic mean of values.\"\"\"\n"
        "    if not values:\n"
        "        return 0.0\n"
        "    # BUG: divides by len-1 instead of len\n"
        "    return sum(values) / (len(values) - 1)\n"
    )
    tests_dir = workspace / "tests"
    tests_dir.mkdir()
    (tests_dir / "test_stats.py").write_text(
        "from stats import average\n"
        "def test_average_basic():\n"
        "    assert average([1, 2, 3, 4]) == 2.5\n"
        "def test_average_empty():\n"
        "    assert average([]) == 0.0\n"
    )

    perms = KitToolPermissions(
        workspace=workspace,
        mode="bypassPermissions",
    )
    client = OpenAIClient(
        api_key=api_key,
        model=MODEL,
        base_url="https://ki-toolbox.scc.kit.edu/api/v1",
        permissions=perms,
    )

    system_prompt = (
        "You are a coding agent with full filesystem access inside the "
        f"workspace {workspace}. Available tools include the standard "
        "ones (read_file, write_file, edit_file, bash, grep_file, "
        "list_files) PLUS extras: bash_background / bash_status / "
        "bash_output / bash_kill for long-running tasks, "
        "task_create / task_update / task_list for multi-step plans, "
        "notebook_read / notebook_edit for .ipynb files. Be concise. "
        "When the result of an edit_file/write_file mentions a 'Tip: "
        "matching test file is X', actually run pytest on X to verify."
    )

    task = (
        f"In {workspace}: tests are failing. Diagnose the bug in "
        f"stats.py, fix it so both tests in tests/test_stats.py pass, "
        f"and verify by running `python -m pytest tests/ -q` from the "
        f"workspace root. Show me the final pytest output."
    )

    print(f"=== KIT ADVANCED LIVE TEST ===")
    print(f"Model:     {MODEL}")
    print(f"Workspace: {workspace}")
    print(f"Initial state: tests FAILING (bug planted in stats.py)")
    print(f"Task:\n  {task}\n")
    print("--- streaming ---")

    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": task},
    ]

    tool_calls: list[dict] = []
    text_chunks: list[str] = []
    t0 = time.monotonic()
    total_in = 0
    total_out = 0

    try:
        for event in client.stream_message(
            messages=messages,
            system=system_prompt,
            max_tokens=8192,
        ):
            if event.type == "tool_use":
                ti = event.tool_input
                if len(ti) > 240:
                    ti = ti[:240] + "..."
                tool_calls.append({
                    "name": event.tool_name.split("__")[-1],
                    "input": ti,
                })
                print(f"  [TOOL] {tool_calls[-1]['name']}({ti})")
            elif event.type == "tool_result":
                out = (event.tool_output or "")[:200].replace("\n", " ")
                more = ("..." if len(event.tool_output or "") > 200 else "")
                print(f"  [RESULT] {out}{more}")
            elif event.type == "text_delta" and event.text:
                text_chunks.append(event.text)
                sys.stdout.write(event.text)
                sys.stdout.flush()
            elif event.type == "message_delta":
                total_in = event.input_tokens
                total_out = event.output_tokens
    except Exception as exc:
        print(f"\n\nFAIL: stream raised: {exc}")
        return 3

    elapsed = time.monotonic() - t0
    print("\n--- end stream ---")
    print(f"Elapsed:    {elapsed:.1f}s")
    print(f"Tokens:     {total_in:,} in / {total_out:,} out")
    print(f"Tool calls: {len(tool_calls)}")
    by_tool: dict[str, int] = {}
    for tc in tool_calls:
        by_tool[tc["name"]] = by_tool.get(tc["name"], 0) + 1
    for name, n in sorted(by_tool.items(), key=lambda kv: -kv[1]):
        print(f"  - {name}: {n}")

    # ---- Outcome assertions -------------------------------------------
    print("\n--- assertions ---")
    failures: list[str] = []

    body = src.read_text()
    if "len(values) - 1" in body:
        failures.append(
            "BUG STILL PRESENT: stats.py still divides by len(values) - 1"
        )
    elif "len(values)" not in body:
        failures.append(
            f"stats.py looks structurally different than expected:\n{body}"
        )

    # Run pytest ourselves to confirm the agent's claim is true.
    import subprocess
    result = subprocess.run(
        ["python", "-m", "pytest", "tests/", "-q"],
        cwd=workspace, capture_output=True, text=True, timeout=30,
    )
    if result.returncode != 0:
        failures.append(
            f"independent pytest STILL FAILS after agent claimed fix:\n"
            f"{result.stdout}\n{result.stderr}"
        )
    else:
        print(f"  independent pytest: PASS ({result.stdout.strip()[-60:]})")

    # Tool sequence: agent must have READ the source before editing.
    seen_seq = [tc["name"] for tc in tool_calls]
    if "read_file" not in seen_seq:
        failures.append("agent never called read_file — edited blind")
    if "edit_file" not in seen_seq and "write_file" not in seen_seq:
        failures.append("agent never edited any file")
    if "bash" not in seen_seq:
        failures.append("agent never ran the verification pytest")

    # Order: read should come before the first edit.
    if "read_file" in seen_seq and (
            "edit_file" in seen_seq or "write_file" in seen_seq):
        first_read = seen_seq.index("read_file")
        first_edit = min(
            (seen_seq.index(t) for t in ("edit_file", "write_file") if t in seen_seq),
            default=None,
        )
        if first_edit is not None and first_read > first_edit:
            failures.append(
                "agent edited before reading — wrong order"
            )

    if failures:
        print("\n=== TEST FAILED ===")
        for f in failures:
            print(f"  - {f}")
        print(f"\nFinal stats.py:\n{body}")
        return 1

    print("\n=== TEST PASSED ===")
    print(f"\nFinal stats.py:\n{body}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
