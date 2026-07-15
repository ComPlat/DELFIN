"""Subagent resource bounds: nesting depth, background concurrency, and a
bounded wait on fanned-out results.

Without these a delegated agent could recursively fan out (4^depth threads +
worktrees), a session could leak unbounded background daemon threads, and one
stalled child could freeze the whole parent turn.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as ac
from delfin.agent.api_client import (
    KitToolPermissions,
    _acquire_bg_subagent_slot, _max_subagent_depth, _subagent_collect_timeout,
)


def _perms(tmp_path, runner, depth=0):
    # depth rides on the perms (not the shared executor); _derive_perms bumps
    # it per child. Simulate a sub-agent by constructing perms at that depth.
    p = KitToolPermissions(workspace=tmp_path, mode="default",
                           subagent_depth=depth)
    p.subagent_runner = runner
    return p


def _args():
    return {"subagent_type": "general-purpose", "description": "d",
            "prompt": "p" * 40}


def _exec(args, perms):
    return json.loads(ac._doc_executor._execute_subagent(args, perms))


# --- Nesting depth guard ----------------------------------------------------

def test_default_depth_cap_is_one():
    assert _max_subagent_depth() == 1


def test_derive_perms_bumps_depth(tmp_path):
    from delfin.agent.subagents import _derive_perms
    parent = KitToolPermissions(workspace=tmp_path, mode="default")
    child = _derive_perms(parent, "default")
    assert child.subagent_depth == 1
    grandchild = _derive_perms(child, "default")
    assert grandchild.subagent_depth == 2


def test_subagent_refused_at_depth_cap(tmp_path):
    perms = _perms(tmp_path, lambda **k: {"final_text": "x"}, depth=1)
    payload = _exec(_args(), perms)
    assert "error" in payload and "nest" in payload["error"].lower()


def test_subagent_allowed_at_top_level(tmp_path):
    seen = {}
    perms = _perms(tmp_path, lambda **k: seen.update(k) or {"final_text": "ok"})
    payload = _exec(_args(), perms)
    assert "error" not in payload
    assert seen  # the runner actually ran


def test_depth_cap_env_override(monkeypatch, tmp_path):
    monkeypatch.setenv("DELFIN_MAX_SUBAGENT_DEPTH", "2")
    assert _max_subagent_depth() == 2
    perms = _perms(tmp_path, lambda **k: {"final_text": "x"}, depth=1)
    payload = _exec(_args(), perms)
    assert "error" not in payload


# --- Background concurrency cap ---------------------------------------------

def test_background_subagent_cap_saturates(monkeypatch):
    monkeypatch.setenv("DELFIN_MAX_BG_SUBAGENTS", "2")
    monkeypatch.setattr(ac, "_BG_SUBAGENT_SEM", None)   # rebuild reading env=2
    r1 = _acquire_bg_subagent_slot()
    r2 = _acquire_bg_subagent_slot()
    r3 = _acquire_bg_subagent_slot()
    try:
        assert r1 and r2 and r3 is None, "cap did not saturate at 2"
    finally:
        if r1:
            r1()
        if r2:
            r2()
    # a slot is free again after release
    r4 = _acquire_bg_subagent_slot()
    assert r4 is not None
    r4()


def test_background_path_refused_when_saturated(monkeypatch, tmp_path):
    monkeypatch.setenv("DELFIN_MAX_BG_SUBAGENTS", "1")
    monkeypatch.setattr(ac, "_BG_SUBAGENT_SEM", None)
    hold = _acquire_bg_subagent_slot()               # occupy the only slot
    try:
        perms = _perms(tmp_path, lambda **k: {"final_text": "x"})
        args = {**_args(), "background": True}
        payload = _exec(args, perms)
        assert "error" in payload and "background" in payload["error"].lower()
    finally:
        hold()


# --- Bounded result wait ----------------------------------------------------

def test_collect_timeout_is_bounded():
    t = _subagent_collect_timeout()
    assert 0 < t < 3600      # finite, sane upper bound
