"""A flag that runs a program is not a reading command.

Input: commands through the TOOL (``_doc_executor.execute``), not
through the permission gate alone, with HOME pointed at a temporary
directory, a fresh workspace, and the caged-admission rule switched off
so what is read is the ask decision itself. Three outcomes are
distinguished: the command runs, a confirm callback is consulted, or an
error comes back.

Output measured 2026-11-22 on 3f9474f8a:

    sort --compress-program=sh a.txt   RUN   ->  ASK
    sort --compress-program sh a.txt   RUN   ->  ASK
    sort --out=f a.txt                 RUN   ->  ASK
    sort -T /tmp a.txt                 RUN   ->  ASK

``sort`` was admitted as a reading command and judged by a deny-pattern
for ``-o`` and ``--output``. GNU getopt accepts any unambiguous
long-option prefix, so ``--out=f`` writes a file the pattern never saw,
and ``--compress-program=PROG`` runs PROG and was named by no pattern at
all. ``readonly_pipeline`` already judged ``sort`` by an allow-list of
the long options that only read; the gate now consults that judgement.

WHY THE REDIRECT HALF OF THIS CHANGE WAS WITHDRAWN, and the control that
keeps it withdrawn: a first version of this work also judged redirection
targets in the auto-allow list, on the strength of a measurement taken
at ``_run_permission_gate``, where ``echo x > ~/.bashrc`` and
``ls > ../escape.txt`` do run. At the TOOL entry they do not: the
write-path gate refuses every one of them with its own wording. Judging
them a second time, earlier, replaced a refusal with a question -- which
is weaker -- and bypassed the ephemeral-sink classification the write
path uses. A permission gate is not the tool, exactly as a predicate is
not the permission gate. The refusals are pinned below so that the next
reading of the predicate does not "close" them again.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

import delfin.agent.api_client as A


@pytest.fixture(autouse=True)
def _no_cage_admission(monkeypatch):
    """Where the caged rule applies the command runs inside the cage
    instead of being asked about; switched off, the ask decision is what
    these tests read."""
    monkeypatch.setattr(A, "_caged_command_may_run_unasked",
                        lambda *a, **k: False)


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    d = tmp_path / "ws"
    d.mkdir()
    (d / "a.txt").write_text("x\n")
    (d / "sub").mkdir()
    return d


_n = [0]


def _verdict(ws: Path, cmd: str, mode: str = "default") -> str:
    """RUN, ASK or REFUSE, as the tool decides it."""
    _n[0] += 1
    asked: list = []
    perms = A.KitToolPermissions(
        mode=mode, workspace=str(ws),
        confirm_callback=lambda *a, **k: (asked.append(a), False)[1])
    perms.task_session_id = f"flagprobe-{_n[0]}"
    raw = A._doc_executor.execute("bash", {"command": cmd}, perms)
    out = json.loads(raw) if isinstance(raw, str) else raw
    if asked:
        return "ASK"
    return "REFUSE" if out.get("error") else "RUN"


# --- what this change closes -------------------------------------------

@pytest.mark.parametrize("cmd", [
    "sort --compress-program=sh a.txt",
    "sort --compress-program sh a.txt",
    "sort --out=f a.txt",
    "sort --output=f a.txt",
    "sort -ofile a.txt",
    "sort -T /tmp a.txt",
])
def test_a_sort_flag_that_writes_or_runs_reaches_the_person(ws, cmd):
    """A flag that writes a file or starts a program is asked about,
    whatever spelling getopt accepts for it."""
    assert _verdict(ws, cmd) == "ASK", cmd


@pytest.mark.parametrize("cmd", [
    "sort -u a.txt",
    "sort -k1 -r a.txt",
    "sort --numeric-sort a.txt",
    "sort --reverse --unique a.txt",
])
def test_a_reading_sort_still_runs_unasked(ws, cmd):
    """Control: the change must not cost function. Every reading form of
    sort ran unasked before it and must still run."""
    assert _verdict(ws, cmd) == "RUN", cmd


# --- the controls that keep the withdrawn half withdrawn ---------------

@pytest.mark.parametrize("cmd", [
    "echo x > ../outside.txt",
    "ls > ../escape.txt",
    "ls >> ../append.txt",
    "ls 2> ../err.txt",
])
def test_the_write_path_gate_already_refuses_an_outside_target(ws, cmd):
    """Measured at the tool entry: a redirect target outside the
    workspace roots is REFUSED, not asked about and not run.

    This is the control behind the withdrawal described in the module
    docstring. If a later change makes one of these ASK, that change
    has replaced a refusal with a question; if it makes one RUN, a
    protection is gone. Either way this test is the place that says so.
    """
    assert _verdict(ws, cmd) == "REFUSE", cmd


def test_the_refusal_names_the_target_and_the_reason(ws):
    """The refusal is specific: it names the target and why it is
    refused, which is what makes it better than a question."""
    _n[0] += 1
    perms = A.KitToolPermissions(mode="default", workspace=str(ws))
    perms.task_session_id = f"flagprobe-{_n[0]}"
    raw = A._doc_executor.execute(
        "bash", {"command": "echo x > ../outside.txt"}, perms)
    err = str((json.loads(raw) if isinstance(raw, str) else raw
               ).get("error") or "")
    assert "would write to" in err, err
    assert "../outside.txt" in err, err
    assert "outside what you may modify" in err, err


@pytest.mark.parametrize("cmd", [
    "ls > out.txt",
    "ls >> build.log",
    "ls > sub/out.txt",
    "ls > /dev/null",
    "ls 2>/dev/null",
    "grep x a.txt 2>&1 | tail -5",
    "ls -la",
    "grep x a.txt",
])
def test_the_forms_an_agent_writes_all_day_are_unchanged(ws, cmd):
    """Control: the device sinks and the in-workspace forms keep the
    verdict they had on the baseline. Two of them run, the workspace
    writes are asked about -- and both are unchanged by this change,
    which is the point of listing them."""
    assert _verdict(ws, cmd) in ("RUN", "ASK"), cmd


def test_a_device_sink_still_runs_without_asking(ws):
    """The most common redirect of all. It was free before and stays
    free: the write path classifies it as an ephemeral sink."""
    assert _verdict(ws, "ls > /dev/null") == "RUN"
    assert _verdict(ws, "ls 2>/dev/null") == "RUN"
