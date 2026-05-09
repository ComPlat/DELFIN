#!/usr/bin/env python3
"""Live test of the headline project-dev workflow.

Asks the agent to set up a fresh Python project from scratch using
the project_dev bundle: register the directory, create a venv,
install pytest, write source + test, run the test. This is the
exact workflow the user wants for their own projects.

The agent starts in a small "host" workspace; the actual project
sits in a separate directory the agent must register via
remember_permission_bundle. If the agent doesn't pick that up, it
will fail trying to bash inside the unrecognised path — making
this test a real check of the project_dev bundle integration.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from delfin.agent import audit_log as al
from delfin.agent import kit_settings as ks
from delfin.agent.api_client import KitToolPermissions, OpenAIClient


MODEL = sys.argv[1] if len(sys.argv) > 1 else "azure.gpt-5.1"


def main() -> int:
    api_key = os.environ.get("KIT_TOOLBOX_API_KEY", "")
    if not api_key:
        print("FAIL: KIT_TOOLBOX_API_KEY not set", file=sys.stderr)
        return 2

    tmp = Path(tempfile.mkdtemp(prefix="kit_proj_"))
    host = tmp / "host_workspace"; host.mkdir()
    project = tmp / "myproj"; project.mkdir()
    audit_log = tmp / "audit.log"
    al._default_log_path = lambda: audit_log
    # Redirect KIT settings persistence to tmp so we don't pollute
    # the user's real ~/.delfin between runs.
    ks.USER_SETTINGS_PATH = tmp / "user_kit.json"

    seen_confirms: list = []

    def confirm(name: str, args: dict, preview: str) -> bool:
        seen_confirms.append({"name": name, "preview": preview[:120]})
        return True

    perms = KitToolPermissions(
        workspace=host,
        mode="default",                 # NOT bypass — exercise the bundle path
        confirm_callback=confirm,
    )
    client = OpenAIClient(
        api_key=api_key,
        model=MODEL,
        base_url="https://ki-toolbox.scc.kit.edu/api/v1",
        permissions=perms,
    )

    system_prompt = (
        f"You are a coding agent. Default workspace: {host}. The user's "
        f"actual project lives at {project}, which is OUTSIDE the default "
        "workspace. Before you can write or run anything inside it, you "
        "must register it via remember_permission_bundle (profile="
        "'project_dev', directory=<absolute path>, scope='repo'). After "
        "that the venv-creation, pip install, python and pytest commands "
        "run inside that path will be auto-allowed. Use these tools: "
        "read_file / write_file / edit_file / bash / "
        "remember_permission_bundle. Always pass absolute paths for files "
        "inside the user's project directory. Be concise."
    )

    task = (
        f"Set up a fresh Python project at {project}. "
        f"1) Register the directory via remember_permission_bundle with "
        f"profile='project_dev'. "
        f"2) Create a venv at {project}/.venv-myproj using `python -m "
        f"venv`. "
        f"3) Install pytest into the venv. "
        f"4) Write {project}/calc.py with a function `multiply(a, b)` "
        f"returning a*b. "
        f"5) Write {project}/test_calc.py with a test that asserts "
        f"`multiply(3, 4) == 12`. "
        f"6) Run pytest from inside the venv to confirm it passes. "
        f"Stop after the test passes. Show me the final pytest output."
    )

    print(f"=== KIT PROJECT-DEV WORKFLOW LIVE TEST ===")
    print(f"Model:     {MODEL}")
    print(f"Host ws:   {host}")
    print(f"Project:   {project}  (outside host, must be granted)")
    print(f"Audit:     {audit_log}")
    print(f"Task:\n  {task}\n")
    print("--- streaming ---")

    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": task},
    ]

    tool_calls: list[dict] = []
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
                name = event.tool_name.split("__")[-1]
                tool_calls.append({"name": name, "input": ti})
                print(f"  [TOOL] {name}({ti})")
            elif event.type == "tool_result":
                out = (event.tool_output or "")[:200].replace("\n", " ")
                more = ("..." if len(event.tool_output or "") > 200 else "")
                print(f"  [RESULT] {out}{more}")
            elif event.type == "text_delta" and event.text:
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
    by_tool: dict = {}
    for tc in tool_calls:
        by_tool[tc["name"]] = by_tool.get(tc["name"], 0) + 1
    for name, n in sorted(by_tool.items(), key=lambda kv: -kv[1]):
        print(f"  - {name}: {n}")
    print(f"Confirms seen: {len(seen_confirms)}")
    for c in seen_confirms:
        print(f"  - {c['name']}: {c['preview']}")

    # ---- Outcome assertions ------------------------------------------
    print("\n--- assertions ---")
    failures: list[str] = []

    # 1. Bundle was actually registered.
    if not any(c["name"] == "remember_permission_bundle"
               for c in seen_confirms):
        failures.append("agent never called remember_permission_bundle")

    # 2. venv exists.
    venv_python = project / ".venv-myproj" / "bin" / "python"
    if not venv_python.exists():
        failures.append(f"no venv python at {venv_python}")

    # 3. Source files written correctly.
    calc = project / "calc.py"
    test = project / "test_calc.py"
    if not calc.exists():
        failures.append(f"calc.py missing at {calc}")
    elif "def multiply" not in calc.read_text():
        failures.append(f"calc.py doesn't define multiply:\n{calc.read_text()}")
    if not test.exists():
        failures.append(f"test_calc.py missing at {test}")

    # 4. Independently verify pytest actually passes.
    if venv_python.exists():
        result = subprocess.run(
            [str(venv_python), "-m", "pytest", "-q"],
            cwd=project, capture_output=True, text=True, timeout=30,
        )
        if result.returncode != 0:
            failures.append(
                f"independent pytest run still fails:\n"
                f"stdout: {result.stdout}\nstderr: {result.stderr}"
            )
        else:
            print(f"  independent pytest: PASS ({result.stdout.strip()[-60:]})")

    # 5. Tool sequence sanity.
    seen_names = [tc["name"] for tc in tool_calls]
    if "remember_permission_bundle" not in seen_names:
        failures.append("agent never invoked remember_permission_bundle (tool side)")
    if "bash" not in seen_names:
        failures.append("agent never invoked bash (no venv / no pip / no pytest)")

    if failures:
        print("\n=== TEST FAILED ===")
        for f in failures:
            print(f"  - {f}")
        print(f"\nProject contents:")
        for p in sorted(project.rglob("*")):
            if p.is_file() and ".venv-myproj" not in str(p):
                print(f"  {p.relative_to(project)}")
        return 1

    print("\n=== TEST PASSED ===")
    print(f"\ncalc.py:\n{calc.read_text()}")
    print(f"\ntest_calc.py:\n{test.read_text()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
