"""After `/command ` the window offers what the command takes.

session.choices with what="args" and a command answers, for every
/command the session lists with takesArgs, either values to pick from
(a value is what goes after the command, as typed) or a usage line when
the values cannot be listed. A command it does not know gets neither.
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent.service import choices
from tui_worker_client import WorkerClient


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_every_command_that_takes_arguments_has_values_or_a_usage(tmp_path):
    with WorkerClient(tmp_path) as w:
        listed = w.request("session.commands")["commands"]
        for c in listed:
            if not c["takesArgs"]:
                continue
            got = w.request("session.choices", {"what": "args", "command": c["name"]})
            assert got["choices"] or got.get("usage"), c["name"]
            for choice in got["choices"]:
                assert choice["value"] and "\n" not in choice["value"], (c["name"], choice)


def test_the_values_each_command_takes(tmp_path):
    _write(tmp_path / "workspace" / ".delfin" / "commands" / "deploy.md",
           "---\ndescription: Ship it\nargument-hint: <environment>\n---\nDeploy $ARGUMENTS\n")
    with WorkerClient(tmp_path) as w:
        def values(command: str) -> list[str]:
            return [c["value"] for c in w.request("session.choices", {"what": "args", "command": command})["choices"]]

        model = w.request("session.choices", {"what": "args", "command": "/model"})
        assert model["current"] == "fake-model" and "fake-model" in [c["value"] for c in model["choices"]]
        assert values("/effort") == list(choices.EFFORTS[1:])
        assert values("/permissions")[:3] == ["plan", "default", "acceptEdits"]
        assert "status" in values("/git") and "role" in values("/git")
        assert values("/undo-file") == ["list", "last", "turn", "session"]
        assert values("/memories") == ["global", "all"]
        assert "watch on" in values("/jobs")
        assert "grant hooks" in values("/trust")
        assert "status" in values("/help")
        deploy = w.request("session.choices", {"what": "args", "command": "/deploy"})
        assert deploy["usage"] == "/deploy <environment>"
        compact = w.request("session.choices", {"what": "args", "command": "/compact"})
        assert compact["choices"] == [] and compact["usage"]


def test_an_unknown_command_gets_nothing(tmp_path):
    with WorkerClient(tmp_path) as w:
        got = w.request("session.choices", {"what": "args", "command": "/no-such-thing"})
        assert got == {"current": "", "choices": [], "usage": ""}
