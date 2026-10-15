"""Scripted turns: the agent engine's surface with no model behind it.

The TUI has to be built and tested against every kind of turn -- streamed
answers, thinking, tool calls, refused tools, approvals, plan approval,
questions, failures, long output, turns that are steered or cancelled --
without an API key, without network and without spending anything.
``FakeEngine`` stands where ``AgentEngine`` stands and plays one of a few
scripts, chosen by words in the prompt:

    edit    a file edit that needs approval (write tool, with a diff)
    bash    a shell command that needs approval
    plan    a plan put up for approval
    ask     a multiple-choice question for the user (``multi``: several
            answers allowed)
    fail    the turn ends with an error
    outage  the turn ends on an endpoint outage, which the session waits
            out before it continues on its own
    slow    a slow answer (a few seconds), for cancelling
    brief   an answer that takes about a second, for queueing
    steer   an answer that reads a message sent into the running turn
    noisy   prints to stdout and runs a child that does the same
    longtool a tool call that runs 20 s and does not stop when asked (as a
            running command does not)
    secret  a tool call whose input and output carry a credential
    latin   a tool result with a file name that is not UTF-8
    long    a long answer (hundreds of lines)
    (anything else: a short answer after one read-only tool call)

It honours the same callbacks, stop and steer paths and status counters the
real engine has, so the worker drives it exactly as it drives the real one.
"""

from __future__ import annotations

import json
import queue
import threading
import time
import uuid
from types import SimpleNamespace
from typing import Any, Callable

_EDIT_DIFF = """--- a/delfin/xtb.py
+++ b/delfin/xtb.py
@@ -41,5 +41,8 @@ def read_final_energy(path):
     for line in lines:
         if "FINAL SINGLE POINT ENERGY" in line:
-            energy = float(line.split()[3])
+            match = _ENERGY_RE.search(line)
+            if match is None:
+                raise ValueError(f"no energy on line: {line!r}")
+            energy = float(match.group(1))
             break
     return energy
"""


class FakePermissions:
    """The parts of KitToolPermissions the worker and the REPL touch."""

    def __init__(self, mode: str = "plan"):
        self.mode = mode
        self.confirm_callback: Callable[[str, dict, str], bool] | None = None
        self.ask_user_callback: Callable[[dict], dict] | None = None
        self.plan_approval_callback: Callable[[str], dict] | None = None
        self.presence_key = ""


class FakeEngine:
    """Plays scripted turns through the AgentEngine callback surface."""

    provider = "fake"

    def __init__(self, workspace: str = "", *, mode: str = "solo",
                 delay: float = 0.02):
        self.repo_dir = workspace
        self.mode = mode
        self.session_id = uuid.uuid4().hex
        self.messages: list[dict[str, Any]] = []
        self.effort = ""
        self.delay = float(delay)
        self.kit_permissions = FakePermissions()
        self.client = SimpleNamespace(model="fake-model", effort="",
                                      signal_stop=self._signal_stop)
        self.last_compaction_info = None
        self.last_turn_stop_reason = "end_turn"
        self.pending_turn_continuation = None
        self.start_fresh = ""
        self.pause_key = ""
        self._stop = threading.Event()
        self._steer: queue.Queue[str] = queue.Queue()
        self._status = {"input_tokens": 0, "output_tokens": 0,
                        "output_tokens_live": 0, "cost_usd": 0.0,
                        "cached_tokens": 0, "provider": "fake",
                        "mode": mode, "session_id": self.session_id}

    # -- the surface the worker and TerminalAgent use --------------------------
    def set_kit_permission_mode(self, mode: str, *, announce: bool = True) -> bool:
        self.kit_permissions.mode = mode
        return True

    def clear_turn_continuation(self) -> None:
        self.pending_turn_continuation = None

    def set_kit_confirm_callback(self, callback) -> bool:
        self.kit_permissions.confirm_callback = callback
        return True

    def persist_kit_pattern(self, pattern: str, *, kind: str = "allow"):
        return True, f"would persist {kind} pattern {pattern}"

    def add_kit_read_dir(self, path, *, persist: bool = False):
        return True, f"{path} readable"

    def get_status(self) -> dict[str, Any]:
        status = dict(self._status)
        status["session_id"] = self.session_id
        return status

    def request_stop(self) -> None:
        self._stop.set()

    def _signal_stop(self) -> None:
        self._stop.set()

    def clear_stop(self) -> None:
        self._stop.clear()

    def steer(self, text: str, *, on_delivered=None) -> bool:
        if not str(text or "").strip():
            return False
        self._steer.put(str(text))
        return True

    def drain_undelivered_steers(self) -> list[str]:
        """What was steered and never read (the real client's
        ``_drain_steer(delivered=False)``)."""
        out: list[str] = []
        while True:
            try:
                out.append(self._steer.get_nowait())
            except queue.Empty:
                return out

    def note_person_message(self) -> None:
        pass

    def estimate_context_tokens(self) -> int:
        return sum(len(str(m.get("content", ""))) for m in self.messages) // 4

    def export_state(self) -> dict[str, Any]:
        return {"messages": list(self.messages), "mode": self.mode,
                "session_id": self.session_id}

    def restore_state(self, data: dict) -> None:
        self.messages = list(data.get("messages") or [])
        self.session_id = str(data.get("session_id") or self.session_id)

    #: What the fake's gate lets through without asking (the real gate's
    #: read-only commands, roughly).
    _READ_ONLY = ("echo", "ls", "cat", "pwd", "sleep", "true")

    def run_gated_bash(self, command: str) -> str:
        """``!command`` the way the real engine runs it: through the same
        confirmation the model's own shell calls go through (read-only
        commands pass), ``sleep N`` takes N seconds, and a stop request
        ends it with 130 (contained_run polls ``should_stop``)."""
        words = str(command or "").split()
        head = words[0] if words else ""
        if "=" in head:                       # VAR=value cmd ...
            head = next((w for w in words if "=" not in w), "")
        cb = self.kit_permissions.confirm_callback
        if head not in self._READ_ONLY and cb is not None and \
                self.kit_permissions.mode in ("plan", "default"):
            if not cb("bash", {"command": command}, f"$ {command}"):
                return json.dumps({"error": "the user refused this command"})
        if head == "sleep" and len(words) == 2:
            until = time.monotonic() + float(words[1])
            while time.monotonic() < until:
                time.sleep(min(0.05, max(0.0, until - time.monotonic())))
                if self._stop.is_set():
                    return json.dumps({"exit_code": 130, "stdout": "",
                                       "stderr": "stopped on request"})
        return json.dumps({"exit_code": 0, "stdout": f"(fake) {command}\n",
                           "stderr": ""})

    # -- one turn ----------------------------------------------------------------
    def stream_response(self, user_message: str, memory_context: str = "",
                        on_token=None, on_tool_use=None, on_tool_result=None,
                        on_permission_denied=None, on_thinking=None,
                        thinking_budget: int = 0, max_tokens: int = 0,
                        images=None, on_notice=None, on_tool_result_meta=None,
                        on_wait=None) -> str:
        self._stop.clear()
        prompt = str(user_message or "")
        self.messages.append({"role": "user", "content": prompt})
        self._status["input_tokens"] += 120 + len(prompt) // 4
        words = prompt.lower()
        out: list[str] = []

        def say(text: str, step: int = 6) -> None:
            for i in range(0, len(text), step):
                if self._stop.is_set():
                    return
                piece = text[i:i + step]
                out.append(piece)
                self._status["output_tokens_live"] += 1
                if on_token:
                    on_token(piece)
                time.sleep(self.delay)

        def think(text: str) -> None:
            if on_thinking:
                on_thinking(text)

        def tool(name: str, args: dict, result: str, *, ok: bool = True,
                 error: str = "") -> None:
            if self._stop.is_set():
                return
            if on_tool_use:
                on_tool_use(name, json.dumps(args))
            time.sleep(self.delay * 3)
            if on_tool_result:
                on_tool_result(name, result)
            if on_tool_result_meta:
                meta = {"ok": ok, "chars": len(result), "truncated": False}
                if error:
                    meta["error"] = error
                on_tool_result_meta(name, meta)

        def confirm(name: str, args: dict, preview: str) -> bool:
            cb = self.kit_permissions.confirm_callback
            return bool(cb(name, args, preview)) if cb else True

        think("Reading the request and deciding what to look at first.")
        if "fail" in words:
            raise RuntimeError("scripted failure from the fake engine "
                               "(not an endpoint problem)")
        if "edit" in words:
            tool("read_file", {"path": "delfin/xtb.py"}, "212 lines")
            if confirm("edit_file", {"path": "delfin/xtb.py"}, _EDIT_DIFF):
                tool("edit_file", {"path": "delfin/xtb.py"}, _EDIT_DIFF)
                say("The parser now reads the energy with a regular expression.")
            else:
                if on_permission_denied:
                    on_permission_denied("edit_file")
                say("Understood, I left `delfin/xtb.py` unchanged.")
        elif "bash" in words:
            command = "pytest -q tests/test_xtb.py"
            if confirm("bash", {"command": command}, f"$ {command}"):
                tool("bash", {"command": command},
                     json.dumps({"exit_code": 0, "stdout": "3 passed in 0.42s\n"}))
                say("All three tests pass.")
            else:
                if on_permission_denied:
                    on_permission_denied("bash")
                say("The tests were not run.")
        elif "plan" in words:
            cb = self.kit_permissions.plan_approval_callback
            verdict = cb("1. Read delfin/xtb.py\n2. Fix the energy parser\n"
                         "3. Run the tests") if cb else {"approved": True}
            if verdict.get("approved"):
                self.kit_permissions.mode = str(verdict.get("new_mode") or "default")
                say(f"Plan approved; working in {self.kit_permissions.mode}.")
            else:
                say("Plan rejected; staying in plan mode.")
        elif "ask" in words:
            # The payload the real engine hands the callback
            # (api_client._execute_ask_user_question, normalised).
            cb = self.kit_permissions.ask_user_callback
            multi = "multi" in words
            reply = cb({
                "question": ("Which properties should be computed?" if multi
                             else "Which solvent model?"),
                "header": "Properties" if multi else "Solvent",
                "options": ([{"label": "IR", "description": "frequencies"},
                             {"label": "NMR", "description": "shieldings"},
                             {"label": "UV-Vis", "description": "TD-DFT"}]
                            if multi else
                            [{"label": "CPCM", "description": "conductor-like"},
                             {"label": "SMD", "description": "density-based"},
                             {"label": "none", "description": "gas phase"}]),
                "multiSelect": multi}) if cb else {}
            picked = (reply or {}).get("answers") or ["(no answer)"]
            say(f"Using {', '.join(picked)}.")
        elif "outage" in words:
            # Worded like the gateway's own outage, which the REPL waits
            # out (api_client._INFRA_EXHAUSTION_MARKERS).
            raise RuntimeError("the endpoint is temporarily unavailable "
                               "(fake outage)")
        elif "slow" in words:
            self.delay = max(self.delay, 0.05)
            say("Working through a long calculation step by step. " * 4, step=4)
        elif "brief" in words:
            say("One moment. ")
            until = time.monotonic() + 1.0
            while time.monotonic() < until and not self._stop.is_set():
                time.sleep(0.02)
            say("Done.")
        elif "steer" in words:
            say("Starting on the task. ")
            deadline = time.monotonic() + 5.0
            steered = ""
            while time.monotonic() < deadline and not self._stop.is_set():
                try:
                    steered = self._steer.get(timeout=0.05)
                    break
                except queue.Empty:
                    continue
            say(f"Noted mid-turn: {steered}" if steered else "No message arrived.")
        elif "noisy" in words:
            # What real tools do now and then: print, and start a child that
            # writes to the stdout it inherited. Neither may reach the wire.
            import subprocess
            print("stray print from inside a turn", flush=True)
            subprocess.run(["sh", "-c", "echo stray output from a child process"],
                           check=False)
            say("Noise made.")
        elif "tools" in words:
            # One of every kind of row the window draws for a tool.
            tool("grep_file", {"pattern": "FINAL SINGLE POINT ENERGY", "path": "delfin/"},
                 "delfin/xtb.py:42:        if \"FINAL SINGLE POINT ENERGY\" in line:\n"
                 "delfin/orca.py:118:    # FINAL SINGLE POINT ENERGY is the last one\n")
            tool("web_fetch", {"url": "https://www.faccts.de/docs/orca/6.0/manual/"},
                 "ORCA 6.0 Manual\n" + "".join(f"section {i}: keywords\n" for i in range(1, 41)))
            tool("bash", {"command": "pytest -q tests/test_xtb.py"},
                 json.dumps({"exit_code": 1, "stdout": "F..\nFAILED tests/test_xtb.py::test_energy\n"
                             "1 failed, 2 passed in 0.40s\n", "stderr": ""}))
            tool("write_file", {"path": "notes/plan.md", "content": "# Plan\n1. fix\n2. test\n"},
                 "--- /dev/null\n+++ b/notes/plan.md\n@@ -0,0 +1,3 @@\n+# Plan\n+1. fix\n+2. test")
            tool("read_file", {"path": "/etc/shadow"}, "", ok=False,
                 error="outside the workspace: /etc/shadow is not readable here")
            tool("mcp__delfin-docs__search_docs", {"query": "SMD solvation"},
                 "3 results\nSMD: universal solvation model\nCPCM vs SMD\nSolvent keywords")
            tool("subagent", {"subagent_type": "explore", "description": "find the energy parsers"},
                 "Found two parsers: delfin/xtb.py and delfin/orca.py.")
            say("Every kind of tool row, once.")
        elif "longtool" in words:
            if on_tool_use:
                on_tool_use("bash", json.dumps({"command": "sleep 20"}))
            time.sleep(20.0)
            if on_tool_result:
                on_tool_result("bash", "slept 20 s")
            say("Done waiting.")
        elif "secret" in words:
            token = "Bearer " + "f4k3t0k3n" * 4
            tool("bash", {"command": f"curl -H 'Authorization: {token}' "
                                     "https://example.invalid/api"},
                 f"request sent with Authorization: {token}")
            say("Done.")
        elif "latin" in words:
            import os as _os
            name = _os.fsdecode(b"caf\xe9.xyz")
            tool("list_dir", {"path": "."}, f"{name}\nREADME.md")
            say(f"Found {name!r}.")
        elif "long" in words:
            say("".join(f"- line {i:03d} of a long answer\n" for i in range(1, 301)),
                step=200)
        else:
            tool("read_file", {"path": "README.md"}, "# tui_demo\n")
            say("## Summary\n\nThe workspace holds **one** README and nothing "
                f"else to change.\n\n(re: {prompt.strip()[:60]})\n")

        text = "".join(out)
        self.messages.append({"role": "assistant", "content": text})
        self._status["output_tokens"] = self._status["output_tokens_live"]
        self.last_turn_stop_reason = "stop" if self._stop.is_set() else "end_turn"
        return text
