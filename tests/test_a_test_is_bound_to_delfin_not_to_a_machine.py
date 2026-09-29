"""A test may depend on DELFIN. It may not depend on this machine.

Project rule, 2026-09-29: a test binds to DELFIN at most, never to a
particular machine, and the agent builds them that way from the start.
The second half is a prompt rule, pinned by the last test in this file.

A test that asks the host what is installed reports the machine, not the
code. Measured the day this was written: over the files gated on xtb, 639
tests collected and 121 of them skipped once xtb was taken off the PATH.
Those 121 run on a developer box that happens to have xtb and nowhere
else -- so a break in them reaches nobody until somebody with the right
machine runs the suite.

Two instances the same day made the cost concrete. Two installer tests
were carried as "permanently red locally" for weeks because the host had
a system OpenMPI and the test asserted one would be built; they were
green in CI and red for every developer who had MPI. And a shell with
/opt/orca/lib on LD_LIBRARY_PATH made the first `import sqlite3` in the
suite fail, which looked like a branch defect for half an hour.

This is a RATCHET, not a ban. The list below is what existed when the
rule was written; the file fails when a name is added to it, and the
entry is removed as each file is made universal. Nothing here demands
that the existing ones be fixed today -- it demands that the next test
is built the way the rule says.

How to make one universal, in order of preference:

  supply it      give the test its own stand-in on PATH, so the host's
                 copy is never consulted: an `mpirun` whose version does
                 not match, an `xtb` replaying a recorded run. The test
                 then exercises DELFIN's handling on every machine.
  record it      keep a real output as a fixture and parse that. Physics
                 measured once and checked forever beats physics
                 measured nowhere.
  skip it        last, and only with a named condition. A skipped test
                 is a test that runs nowhere, and CI will not tell you.
"""

from __future__ import annotations

import pathlib
import re

#: The files that gated on the host when this rule was written. Shrinking
#: this list is the work; growing it is the thing this file refuses.
_BASELINE = frozenset({
    "test_a_drag_says_what_it_is_pulling_with.py",
    "test_a_finished_scan_says_what_it_made_possible.py",
    "test_a_finished_walk_can_be_priced_again.py",
    "test_a_force_field_with_fixed_bonds_is_refused_where_a_bond_forms.py",
    "test_a_free_energy_at_a_scan_point_says_what_it_is_worth.py",
    "test_a_sandboxed_command_cannot_reach_the_users_sessions.py",
    "test_a_sandboxed_command_reaches_the_network_through_its_proxy.py",
    "test_a_scan_can_be_walked_back_through.py",
    "test_a_scan_says_what_it_left.py",
    "test_a_scan_says_whether_its_own_barrier_can_be_quoted.py",
    "test_a_scan_shows_the_path_it_walked.py",
    "test_an_agent_process_cannot_be_read_by_its_commands.py",
    "test_asking_what_a_structure_is.py",
    "test_gfn_methods_in_the_viewer.py",
    "test_git_works_inside_the_sandbox.py",
    "test_the_budget_prices_a_relaxed_path.py",
    "test_the_derived_roots_bind_the_doc_index.py",
    "test_the_sandbox_holds_on_macos.py",
    "test_tools_that_answer.py",
    "test_what_the_answer_already_computed.py",
})

#: What counts as asking the host: an installed binary, an environment
#: variable, a path on this disk, the platform, a probe of a facility.
_HOST = re.compile(
    r"shutil\.which|\bwhich\(|os\.environ|\.exists\(\)|\.is_dir\(\)"
    r"|sys\.platform|available\(\)")
_SKIPIF = re.compile(r"skipif\s*\(([^\n]{0,200})")

_TESTS = pathlib.Path(__file__).resolve().parent
_SELF = pathlib.Path(__file__).name


def _machine_bound() -> dict:
    out: dict = {}
    for path in sorted(_TESTS.glob("test_*.py")):
        if path.name == _SELF:
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        n = sum(1 for m in _SKIPIF.finditer(text) if _HOST.search(m.group(1)))
        if n:
            out[path.name] = n
    return out


def test_no_new_test_file_is_bound_to_a_machine():
    added = sorted(set(_machine_bound()) - _BASELINE)
    assert not added, (
        "these test files skip depending on what this machine has: "
        + ", ".join(added)
        + ". Give the test its own stand-in on PATH, or record a real "
          "output as a fixture; a skip runs nowhere and CI will not say so."
    )


def test_the_baseline_only_shrinks():
    """A file made universal is removed from the list. If the list names
    a file that no longer gates on the host, the entry is stale and the
    next reader will think the work is still open."""
    stale = sorted(_BASELINE - set(_machine_bound()))
    assert not stale, (
        "these are universal now and can leave the list: "
        + ", ".join(stale))


def test_the_agent_is_told_the_rule_it_is_measured_by():
    """The guard above rejects; the prompt has to teach.

    A ratchet that only fails tells the author the shape is wrong and
    not what to write instead, and it reaches the agent only after the
    file exists. The rule therefore also lives in the solo role prompt
    -- the role that writes DELFIN's tests -- and this pins it there.

    Checked in the COMPOSED prompt, not in the markdown file: the
    loader drops lazy-module sections, so a rule can be present in the
    file and absent from what the model receives. Both prompt modes are
    built, because slim_prompt is a user setting and the rule may not
    depend on it.
    """
    from unittest import mock

    from delfin import user_settings

    for slim in (True, False):
        with mock.patch.object(
                user_settings, "load_settings",
                lambda *a, _s=slim, **k: {"agent": {"slim_prompt": _s}}):
            from delfin.agent.prompt_loader import PromptLoader

            built = PromptLoader().build_system_prompt(
                role_id="solo_agent", mode_id="solo",
                task_text="schreib einen test", session_key="rule-pin")
            text = built if isinstance(built, str) else str(built)
            assert "binds to DELFIN, never to a machine" in text, (
                f"the solo prompt (slim_prompt={slim}) no longer tells the "
                "agent how to build a test; the guard in this file would "
                "then reject a shape the agent was never taught.")
