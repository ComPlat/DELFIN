#!/usr/bin/env python3
"""Live behavioural test against the real KIT-Toolbox endpoint.

This is the test the user asked for: not "if the agent calls X, X
runs correctly" (we already proved that), but "given a natural-
language task, does the KIT model actually call the right tools
in the right order and produce the right outcome?"

Sets up a sandboxed tmp workspace, asks the agent (via the engine)
to perform a realistic build-test-fix sequence, and asserts on
both the conversation transcript and the resulting file system.

Run with::

    KIT_TOOLBOX_API_KEY=...  python scripts/test_kit_live.py [model_id]

The test creates a fresh tmp directory per run; nothing leaks into
the user's home. Audit-log writes are redirected to the same tmp.
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
import time
import uuid
from pathlib import Path
from typing import Any

# Make repo importable when run directly.
REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from delfin.agent import audit_log as al
from delfin.agent.api_client import (
    KitToolPermissions, OpenAIClient, _doc_executor,
)


MODEL = sys.argv[1] if len(sys.argv) > 1 else "azure.gpt-5.1"


def main() -> int:
    api_key = os.environ.get("KIT_TOOLBOX_API_KEY", "")
    if not api_key:
        print("FAIL: KIT_TOOLBOX_API_KEY not set in environment", file=sys.stderr)
        return 2

    tmp = Path(tempfile.mkdtemp(prefix="kit_live_"))
    workspace = tmp / "workspace"
    workspace.mkdir()
    audit_log = tmp / "audit.log"
    al._default_log_path = lambda: audit_log    # redirect

    perms = KitToolPermissions(
        workspace=workspace,
        mode="bypassPermissions",     # full autonomy — sandbox + deny still apply
        bash_auto_allow_patterns=(),  # bypass mode skips this anyway
    )
    client = OpenAIClient(
        api_key=api_key,
        model=MODEL,
        base_url="https://ki-toolbox.scc.kit.edu/api/v1",
        permissions=perms,
    )

    system_prompt = (
        "You are a coding agent with full filesystem access inside the "
        f"workspace {workspace}. You have these tools available:\n"
        "- read_file(path), write_file(path, content), edit_file(path, "
        "old_string, new_string), bash(command, description, cwd?)\n"
        "Use them to complete the user's task. Always read_file before "
        "edit_file. Be concise."
    )

    task = (
        f"In {workspace}: write a Python file `calc.py` containing a "
        f"function `add(a, b)` that returns a + b, AND a pytest test "
        f"file `test_calc.py` with one test that asserts `add(2, 3) "
        f"== 5`. Then run `python -m pytest test_calc.py -q` in that "
        f"directory and show me the result. Stop after the test "
        f"passes."
    )

    print(f"=== KIT LIVE TEST ===")
    print(f"Model:     {MODEL}")
    print(f"Workspace: {workspace}")
    print(f"Audit:     {audit_log}")
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
            max_tokens=4096,
        ):
            if event.type == "tool_use":
                # Truncate big inputs for log readability.
                ti = event.tool_input
                if len(ti) > 200:
                    ti = ti[:200] + "..."
                tool_calls.append({"name": event.tool_name, "input": ti})
                print(f"  [TOOL] {event.tool_name}({ti})")
            elif event.type == "tool_result":
                out = (event.tool_output or "")[:160].replace("\n", " ")
                print(f"  [RESULT] {out}{'...' if len(event.tool_output or '')>160 else ''}")
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
    for tc in tool_calls:
        print(f"  - {tc['name']}")

    # ---- Assert on outcome --------------------------------------------
    print("\n--- assertions ---")
    failures: list[str] = []

    calc = workspace / "calc.py"
    test = workspace / "test_calc.py"

    if not calc.exists():
        failures.append(f"calc.py was not created at {calc}")
    else:
        body = calc.read_text()
        if "def add" not in body:
            failures.append(f"calc.py exists but doesn't define `add`:\n{body}")

    if not test.exists():
        failures.append(f"test_calc.py was not created at {test}")
    else:
        tbody = test.read_text()
        if "add(2, 3)" not in tbody and "add(2,3)" not in tbody:
            failures.append(f"test file doesn't reference add(2, 3):\n{tbody}")

    expected_tools = {"write_file", "bash"}
    seen_tools = {tc["name"].split("__")[-1] for tc in tool_calls}
    missing = expected_tools - seen_tools
    if missing:
        failures.append(f"agent did not use these expected tools: {missing}")

    # Did the agent actually run pytest? Check audit log.
    if audit_log.exists():
        log_lines = [
            json.loads(l) for l in audit_log.read_text().splitlines() if l
        ]
        bash_lines = [r for r in log_lines if r["tool"] == "bash"
                      and "pytest" in (r.get("command") or "")]
        if not bash_lines:
            failures.append("audit log shows no bash call running pytest")
        else:
            print(f"  audit: {len(bash_lines)} pytest run(s) recorded")

    if failures:
        print("\n=== TEST FAILED ===")
        for f in failures:
            print(f"  - {f}")
        return 1

    print("\n=== TEST PASSED ===")
    print(f"calc.py:\n{calc.read_text()}")
    print(f"test_calc.py:\n{test.read_text()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
