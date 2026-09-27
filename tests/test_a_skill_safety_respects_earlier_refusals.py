"""Phase 3 of skill_safety: earlier refusals never become a rule.

Whatever the user has refused in this session is never something a
skill may prescribe. The skill review compares the commands in a draft
against the session's RefusalMemory, reusing refusal_memory's own
target extraction (_extract_read_target) and comparison (_norm,
_contains) — nothing re-implemented here.
"""
from __future__ import annotations

import pytest

from delfin.agent import skill_safety
from delfin.agent.refusal_memory import Refusal, RefusalMemory


@pytest.fixture(autouse=True)
def _unwire():
    yield
    skill_safety.wire(refusal_memory=None)


def test_wire_accepts_a_refusal_memory():
    skill_safety.wire(refusal_memory=RefusalMemory())
    assert skill_safety.check("read the results") == []


def test_a_refused_target_in_the_skill_is_a_finding():
    mem = RefusalMemory(entries=[
        Refusal(tool="read_file", target="/etc/hosts",
                reason="system file", time="07:00"),
    ])
    skill_safety.wire(refusal_memory=mem)
    text = "Check the hosts mapping first:\n\n```bash\ncat /etc/hosts\n```\n"
    findings = skill_safety.check(text)
    assert any("refus" in f.lower() and "/etc/hosts" in f
               for f in findings), findings


def test_a_refused_directory_covers_paths_beneath():
    mem = RefusalMemory(entries=[
        Refusal(tool="bash", target="archive/", reason="read only",
                time="07:00", is_dir=True),
    ])
    skill_safety.wire(refusal_memory=mem)
    text = "```bash\ncat archive/2026/run1/out.log\n```"
    findings = skill_safety.check(text)
    assert any("refus" in f.lower() for f in findings), findings


def test_an_unrelated_command_is_not_a_refusal_finding():
    mem = RefusalMemory(entries=[
        Refusal(tool="read_file", target="/etc/hosts",
                reason="system file", time="07:00"),
    ])
    skill_safety.wire(refusal_memory=mem)
    text = "```bash\ngit status\n```"
    findings = skill_safety.check(text)
    assert not any("refus" in f.lower() for f in findings), findings


def test_without_wiring_there_is_no_refusal_check():
    # The default must be None: no session memory wired, no findings
    # from this phase (check() stays usable standalone, as the contract
    # with skill_proposals requires).
    assert skill_safety.check("```bash\ncat /etc/hosts\n```") == []
