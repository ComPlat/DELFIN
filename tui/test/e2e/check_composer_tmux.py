"""The composer, live: `delfin-agent tui` in tmux against a private agent
service, typed into like a person would.

  /st             the command list opens, from the session's own commands
  /dep            a workspace command (.delfin/commands/deploy.md) is offered
  /git , /model   what the command takes, the current model marked
  @del + tab      paths complete from the workspace, folder by folder
  alt+enter       a prompt of two lines goes out as one message
  a large paste   stands in as [Pasted text #1 +30 lines], arrives whole
  up              the history brings the last prompt back
  !               shell mode: "!" as the marker, the border says shell; the
                  line runs as "!echo ..." and the rows head it "$ echo ..."
  a question      opens while the user types: the typing goes on into the
                  prompt, and enter answers with the line; on an empty
                  prompt 2 picks the second choice
  an approval     opens while the user types: the typing (an option key
                  among it) goes on into the prompt, enter sends it as a
                  message, a word typed into the empty prompt is typing --
                  nothing answers it but a lone key (n)

Usage (from the repository root):
    DELFIN_PYTHON=<python with delfin> python3 tui/test/e2e/check_composer_tmux.py
"""

from __future__ import annotations

import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

TUI = Path(__file__).resolve().parents[2]
REPO = TUI.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(TUI / "spike"))

from check_tmux import Tmux  # noqa: E402

from delfin.agent.service.client import ServiceClient  # noqa: E402

PYTHON = os.environ.get("DELFIN_PYTHON") or sys.executable


def type_slowly(tmux: Tmux, text: str, gap: float, appears: str) -> int:
    """Type `text` a key at a time; the index of the key before which
    `appears` was first on screen (-1: never while typing)."""
    seen = -1
    for i, ch in enumerate(text):
        if seen < 0 and appears in tmux.screen():
            seen = i
        tmux.type(ch)
        time.sleep(gap)
    return seen


def question_keeps_the_prompt(tmux: Tmux, expect) -> list[str]:
    """The prompt keeps the keys while a question or an approval is open."""
    problems: list[str] = []
    # A question that opens while the user types: the typing goes on, and
    # enter answers with the line.
    tmux.type("ask about the solvent")
    tmux.key("Enter")
    answer = "toluene, implicit please"
    seen = type_slowly(tmux, answer, 0.12, "Which solvent model?")
    if seen <= 0:
        problems.append(f"the question did not open while typing (before key {seen})\n{tmux.screen()}")
    expect(r"❯ toluene, implicit please", "the typing went on into the prompt")
    time.sleep(0.6)
    tmux.key("Enter")
    expect(r"Using toluene, implicit please\.", "the line typed is the answer", 30)
    expect(r"❯ ask DELFIN", "the prompt after the answer")
    # A numbered choice by key, on an empty prompt.
    tmux.type("ask about the solvent again")
    tmux.key("Enter")
    expect(r"Which solvent model\?", "the second question")
    time.sleep(0.8)
    tmux.key("2")
    expect(r"Using SMD\.", "2 on an empty prompt picks the second choice", 30)
    expect(r"❯ ask DELFIN", "the prompt after the pick")
    # An approval that opens while the user types: nothing typed answers it.
    tmux.type("run the bash tests")
    tmux.key("Enter")
    seen = type_slowly(tmux, "please wait, yes do it", 0.1, "[y] yes")
    if seen <= 0:
        problems.append(f"the approval did not open while typing (before key {seen})\n{tmux.screen()}")
    time.sleep(1.5)
    screen = tmux.screen()
    if "All three tests pass." in screen or "The tests were not run." in screen or "[y] yes" not in screen:
        problems.append(f"typing in progress answered the approval\n{screen}")
    expect(r"❯ please wait, yes do it", "the typing went on into the prompt")
    tmux.key("Enter")                         # a line under an approval: a message, queued
    time.sleep(1.5)
    screen = tmux.screen()
    if "All three tests pass." in screen or "The tests were not run." in screen or "[y] yes" not in screen:
        problems.append(f"enter with a line answered the approval\n{screen}")
    time.sleep(1.2)
    tmux.type("yes")                          # a word into the empty prompt
    time.sleep(1.5)
    screen = tmux.screen()
    if "All three tests pass." in screen or "[y] yes" not in screen:
        problems.append(f"a word typed into the empty prompt answered the approval\n{screen}")
    expect(r"❯ yes", "the word is in the prompt")
    tmux.key("C-u")
    time.sleep(1.3)
    tmux.key("n")                             # a lone key, on an empty prompt
    expect(r"The tests were not run\.", "a lone n on an empty prompt refuses", 30)
    expect(r"\(re: please wait, yes do it", "the queued line ran after the turn", 30)
    time.sleep(0.5)
    return problems


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="tui-composer-"))
    runtime = Path(tempfile.mkdtemp(prefix="dtc", dir="/tmp"))
    ws = tmp / "project"
    for d in (tmp / "home", tmp / "scratch" / ".delfin", ws / "delfin", ws / "docs",
              ws / ".delfin" / "commands"):
        d.mkdir(parents=True)
    (ws / "delfin" / "xtb.py").write_text("")
    (ws / "docs" / "guide.md").write_text("")
    (ws / ".delfin" / "commands" / "deploy.md").write_text(
        "---\ndescription: Ship the build\n---\nDeploy $ARGUMENTS\n")
    env = {"HOME": str(tmp / "home"), "DELFIN_SCRATCH_STATE": str(tmp / "scratch"),
           "DELFIN_RUNTIME_DIR": str(runtime), "PYTHONPATH": str(REPO), "DELFIN_TUI_SCREEN": "inline",
           "PATH": f"{Path.home()}/.bun/bin:{os.environ.get('PATH', '')}"}
    command = (f"cd {REPO} && env {' '.join(f'{k}={v}' for k, v in env.items())} "
               f"{PYTHON} -m delfin.agent.cli tui --fake --cwd {ws}; echo TUI_EXIT=$?")
    problems: list[str] = []
    tmux = Tmux(100, 35)

    def expect(pattern: str, what: str, timeout: float = 15) -> bool:
        if tmux.wait_for(pattern, timeout):
            return True
        problems.append(f"{what}: no {pattern!r} on screen\n{tmux.screen()}")
        return False

    def clear() -> None:
        tmux.key("C-u")
        time.sleep(0.3)

    try:
        tmux.type(command)
        tmux.key("Enter")
        if not expect(r"ask DELFIN", "the window did not come up", 60):
            raise SystemExit
        tmux.type("/st")
        expect(r"/status\s+Model, mode, tokens, cost", "the command list")
        clear()
        tmux.type("/dep")
        expect(r"/deploy\s+Ship the build", "a workspace command")
        clear()
        tmux.type("/git ")
        expect(r"role\s+maintainer or contributor here", "the values /git takes")
        clear()
        tmux.type("/model ")
        expect(r"fake-model\s+\(now\)", "the models, the current one marked")
        clear()
        tmux.type("look at @del")
        expect(r"@delfin/", "the path list")
        tmux.key("Tab")
        expect(r"@delfin/xtb\.py", "the list inside the folder")
        tmux.key("Tab")
        expect(r"❯ look at @delfin/xtb\.py", "the completed path")
        clear()
        # Shell mode: "!" at the start, no completion, backspace leaves it.
        tmux.type("!")
        expect(r"── shell · runs through the agent's checks", "shell mode on the border")
        expect(r"(?m)^ ! a shell command", "the shell marker")
        tmux.type("/st")
        time.sleep(0.8)
        if "Model, mode, tokens, cost" in tmux.screen():
            problems.append(f"completion in shell mode\n{tmux.screen()}")
        for _ in range(4):
            tmux.key("BSpace")
        expect(r"❯ ask DELFIN", "backspace on the ! back to a prompt")
        if "── shell ·" in tmux.screen():
            problems.append(f"shell mode stayed after backspace\n{tmux.screen()}")
        tmux.type("!echo shell-mode-ok")
        expect(r"(?m)^ ! echo shell-mode-ok", "the command in shell mode")
        tmux.key("Enter")
        expect(r"(?m)^ \$ echo shell-mode-ok", "the command's $ line in the rows", 30)
        expect(r"❯ ask DELFIN", "a prompt again after the command")
        # Two lines, one message.
        tmux.type("first line")
        tmux.key("M-Enter")
        tmux.type("second line")
        time.sleep(0.3)
        tmux.key("Enter")
        expect(r"\(re: first line", "the two-line turn", 30)
        # A paste the size of a log.
        log = "\n".join(f"log line {i:02d}" for i in range(30))
        tmux.run("set-buffer", "-b", "log", log)
        tmux.run("paste-buffer", "-p", "-b", "log", "-t", tmux.name)
        expect(r"\[Pasted text #1 \+30 lines\]", "the paste placeholder")
        time.sleep(0.3)
        tmux.key("Enter")
        expect(r"❯ \[Pasted text #1 \+30 lines\]", "the placeholder in the scrollback", 30)
        expect(r"\(re: log line 00", "the pasted turn", 30)
        time.sleep(1.0)
        # The history: the paste comes back as its placeholder.
        tmux.key("Up")
        expect(r"❯ \[Pasted text #1 \+30 lines\]\s*$|❯ \[Pasted text #1 \+30 lines\]",
               "the history")
        tmux.key("C-c")                     # clears what is typed
        time.sleep(0.3)
        problems.extend(question_keeps_the_prompt(tmux, expect))
        tmux.key("C-d")
        exited = tmux.wait_for(r"TUI_EXIT=(\d+)", 30)
        if not exited or exited.group(1) != "0":
            problems.append(f"the window did not leave cleanly: {exited and exited.group(1)}")
        sock = runtime / f"delfin-{os.getuid()}" / "agent.sock"
        with ServiceClient(sock) as c:
            sessions = c.request("sessions.list")["sessions"]
            sid = sessions[0]["sessionId"]
            events = c.request("session.replay", {"sessionId": sid})["events"]
            said = [e["text"] for e in events if e.get("type") == "user.message"]
            if "first line\nsecond line" not in said:
                problems.append(f"the two-line prompt did not arrive as one: {said}")
            if log not in said:
                problems.append("the paste did not arrive whole")
            if "!echo shell-mode-ok" not in said:
                problems.append(f"the shell line did not arrive as typed: {said}")
            if "please wait, yes do it" not in said:
                problems.append(f"the line sent under the approval did not go out as a message: {said}")
            decisions = [e.get("decision") for e in events if e.get("type") == "permission.resolved"]
            if "approved" in decisions:
                problems.append(f"something was approved, and nothing should have been: {decisions}")
            c.request("service.shutdown")
        history = tmp / "scratch" / ".delfin" / "tui_history.jsonl"
        if not history.exists() or "Pasted text #1" not in history.read_text():
            problems.append("the history file does not hold the prompts")
    except SystemExit:
        pass
    finally:
        (tmp / "screen.txt").write_text(tmux.history())
        tmux.close()
        shutil.rmtree(runtime, ignore_errors=True)
    print(f"tmux composer check; scrollback in {tmp / 'screen.txt'}")
    if problems:
        print("FAIL")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("PASS: completion, new lines, pastes, the history and shell mode work live; "
          "a question or an approval keeps the prompt, and typing answers no approval")
    return 0


if __name__ == "__main__":
    sys.exit(main())
