"""Every advertised tool is called once, through the real executor.

Not through its helper. A predicate passing is not the gate passing, and
this codebase has paid for that: nine helper cases green, then the real
gate ran six of the seven things it must refuse.

The record said 33 of 72 tools had never been called in 2509 runs. That
is a statement about the benchmark corpus, not about the tools -- so this
drives all of them and looks at what comes back. Arguments come from each
tool's OWN schema, filled the way a model fills them, in a workspace that
looks like work: a python file, a csv, an ORCA output, a subdirectory.

What it asserts is deliberately narrow, because the narrow thing is what
a caller can rely on:

* the call RETURNS. A tool that raises out of execute() takes the turn
  loop with it; every failure must come back as a message.
* the message is not empty, and is not a traceback. "Something went
  wrong" that the model cannot read is the same as silence.

It does NOT assert that the tool did something useful with these
arguments -- most of them answer "file not found" or "not a PDF", which
is correct. Usefulness is each tool's own test.

Tools that reach off this machine are excluded BY NAME, with the reason,
and the list is checked for completeness: a tool added later is either
driven here or consciously excused, never quietly absent.
"""

from __future__ import annotations

import json
import pytest

from delfin.agent import api_client as A


#: Reaches beyond this process: the network, a person, another machine,
#: a future callback, or a model that costs money. Driving these belongs
#: in their own tests, where the far end is a fake.
NOT_DRIVEN = {
    "web_search": "the network",
    "web_fetch": "the network",
    "push_notification": "reaches a person",
    "draft_email": "reaches a person",
    "remote_trigger": "reaches another machine",
    "cron_create": "schedules future work",
    "cron_delete": "schedules future work",
    "schedule_wakeup": "arms a future callback",
    "watch_job": "arms a future callback",
    "subagent": "spawns a model",
    "orchestrate": "spawns models",
    "session_message": "writes into another session",
    "publish_report": "leaves the workspace",
    "ask_user_question": "needs a human at the other end",
    "exit_plan_mode": "needs a human at the other end",
}


def _tools():
    out = {}
    for entry in A._DOC_TOOLS_OPENAI:
        fn = entry.get("function") or entry
        name = fn.get("name")
        if name:
            out[name] = fn.get("parameters") or {}
    return out


def _value(key: str, spec: dict):
    """What a model would put here, from the schema and the name."""
    if spec.get("enum"):
        return spec["enum"][0]
    kind = spec.get("type")
    if kind in ("integer", "number"):
        return 1
    if kind == "boolean":
        return False
    if kind == "array":
        item = spec.get("items") or {}
        return [_value(key + "_item", item)] if item else []
    if kind == "object":
        return {}
    low = key.lower()
    for needle, value in (
            ("path", "hello.py"), ("command", "echo probe"),
            ("old_string", "return 1"), ("new_string", "return 2"),
            ("pattern", "def "), ("query", "energy"), ("question", "energy"),
            ("content", "probe\n"), ("text", "probe\n"), ("url", "http://localhost:1/"),
            ("column", "energy"), ("id", "0")):
        if needle in low:
            return value
    return "probe"


def _args(schema: dict) -> dict:
    props = (schema or {}).get("properties") or {}
    return {key: _value(key, props.get(key) or {})
            for key in (schema or {}).get("required") or []}


@pytest.fixture()
def workspace(tmp_path):
    ws = tmp_path / "ws"
    (ws / "calc").mkdir(parents=True)
    (ws / "sub").mkdir(parents=True)
    (ws / "hello.py").write_text("def main():\n    return 1\n", encoding="utf-8")
    (ws / "notes.md").write_text("# Notes\n\ntext\n", encoding="utf-8")
    (ws / "data.csv").write_text("name,energy\na,-1.5\n", encoding="utf-8")
    (ws / "calc" / "orca.out").write_text(
        "FINAL SINGLE POINT ENERGY      -613.417262\n", encoding="utf-8")
    (ws / "sub" / "other.py").write_text("x = 1\n", encoding="utf-8")
    return ws


def test_the_exclusion_list_names_only_real_tools():
    """A renamed tool must not leave a silent hole in the coverage."""
    unknown = sorted(set(NOT_DRIVEN) - set(_tools()))
    assert not unknown, f"excused tools that no longer exist: {unknown}"


def test_every_tool_returns_a_message_instead_of_raising(workspace):
    tools = _tools()
    perms = A.KitToolPermissions(workspace=workspace, mode="bypassPermissions")
    raised, empty, traced = [], [], []
    for name in sorted(tools):
        if name in NOT_DRIVEN:
            continue
        try:
            got = A._doc_executor.execute(name, _args(tools[name]), perms)
        except Exception as exc:                       # noqa: BLE001
            raised.append(f"{name}: {type(exc).__name__}: {exc}")
            continue
        text = "" if got is None else str(got)
        if not text.strip():
            empty.append(name)
        elif "traceback (most recent call last)" in text.lower():
            traced.append(f"{name}: {text[:160]}")
    assert not raised, "a tool raised out of execute():\n" + "\n".join(raised)
    assert not traced, "a tool answered with a traceback:\n" + "\n".join(traced)
    assert not empty, f"a tool answered with nothing: {empty}"


def test_the_driven_set_is_most_of_the_surface():
    """A guard on the exclusion list itself: excusing tools is how a
    coverage test quietly stops covering anything."""
    tools = _tools()
    driven = len(tools) - len(set(NOT_DRIVEN) & set(tools))
    assert driven >= int(len(tools) * 0.75), (
        f"only {driven} of {len(tools)} tools are driven here")
