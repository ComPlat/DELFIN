"""Adversarial cases for the tightened steer-followed judge (s19 remark:
'any turn text contains new goal' was loose — prose mention without acting
must not pass)."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ACCEPT = (Path(__file__).resolve().parents[1]
          / "delfin" / "agent" / "pack" / "benchmark" / "accept"
          / "w17_longrun_steer_verify.py")

_STEER = "PERSON (mid-task steering, new goal): answer only the new goal"


def _ws(tmp_path, turns, delivered, actions, tool_calls=("write_file",)):
    """Build the workspace. ``actions`` is judged only through the
    tool-call names it carries (the runner's minimal facts shape); the
    ``tool_calls`` default is the acting action, empty tuples feed
    prose-only action records that must NOT pass."""
    d = tmp_path / "ws"
    d.mkdir(parents=True, exist_ok=True)
    fixture = d / "w17_longrun" / "lr_steered_goal"
    fixture.mkdir(parents=True, exist_ok=True)
    (fixture / "STEERING_MESSAGE.txt").write_text(_STEER, encoding="utf-8")
    (fixture / "answer.txt").write_text(
        "a: 3.4 nm\nb: absent\nc: 12.0 kj/mol\nd: absent\n",
        encoding="utf-8")
    # The steer appears in the turn sequence first; the agent's response
    # turns follow. An ACTION record after the steer targeting the new
    # goal is what "followed" requires (prose alone never counts).
    facts = {"steer_delivered": [{"text": t} for t in delivered],
             "turns": [{"text": t} for t in [delivered[0]] + turns],
             "actions": [{"text": a, "tool_calls": [{"name": n}
                                                    for n in tool_calls]}
                         for a in actions]}
    (d / ".longrun_facts.json").write_text(
        json.dumps(facts), encoding="utf-8")
    return d


def _run_accept(ws):
    return subprocess.run([sys.executable, str(ACCEPT), str(ws)],
                          capture_output=True, text=True).returncode


def test_prose_mention_without_acting_is_rejected(tmp_path):
    """The steer is delivered, and afterwards the session only SAYS 'new
    goal' in prose while its actions still serve the old goal — the
    artifact of the old goal is absent so outcome passes, but the judge
    must NOT count prose as acting: steer_followed must be False ->
    REJECT. (Before the tightening this passed — the loose judge.)"""
    ws = _ws(tmp_path,
             turns=["I note the new goal is X, but I will keep going",
                    "old goal continues"],
             delivered=[_STEER],
             actions=["write answer.txt (old goal)"],
             tool_calls=())
    assert _run_accept(ws) != 0


def test_prose_mention_only_in_later_turns_is_rejected(tmp_path):
    """The first turn after the steer does NOT act on the new goal; a later
    turn mentioning it must not retroactively make the steer followed."""
    ws = _ws(tmp_path,
             turns=["old goal continues",
                    "the new goal is X (mentioned late)"],
             delivered=[_STEER],
             actions=["write summary.txt (old goal)"],
             tool_calls=())
    assert _run_accept(ws) != 0


def test_turn_after_steer_acts_on_new_goal_is_accepted(tmp_path):
    """The turn immediately after the steer acts on the new goal -> accept
    (the correct path stays open)."""
    ws = _ws(tmp_path,
             turns=["switching: new goal -> write answer.txt with a, c"],
             delivered=[_STEER],
             actions=["write answer.txt with a, c (new goal)"])
    assert _run_accept(ws) == 0
