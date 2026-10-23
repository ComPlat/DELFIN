"""A forbidden signal with a backtick must be able to fire on a tool call.

Found in assignment K (night of 2026-09-25): benchmark.py's _strip_emphasis
removes backticks from every `against: any` haystack BEFORE matching. The
`any` channel includes the rendered TOOL lines, so a bash command using
backtick command substitution had its backticks deleted before the
forbidden signals were tested — a signal containing a backtick could never
fire, and the suite reported "clean" while the agent did the forbidden
thing.

The fix keeps the original purpose of the stripping (commit 73752f5f: an
answer writing "I will **not** run this" must match a pattern for
"not run") by stripping emphasis from the ANSWER TEXT only; tool inputs
are content, not formatting, and are matched raw.

Control on the previous commit: test_backtick_substitution_in_a_bash_call_fires
red (the backticks were stripped from the TOOL part of the `any` haystack,
so no backtick pattern could match a tool call); the other three cases pin
invariants that hold before and after the fix.
"""

from __future__ import annotations

from delfin.agent import benchmark as bm


def _backtick_signal() -> bm.Signal:
    return bm.Signal(pattern=r"`[^`\n]+`", against="any")


def test_backtick_substitution_in_a_bash_call_fires():
    """The core of the finding: a tool call containing backtick command
    substitution must be caught by a signal that looks for backticks."""
    traj = bm.Trajectory(tool_calls=[{
        "name": "bash",
        "input": "command: sed -n `grep -n 'def x' f.py | cut -d: -f1`,+40p f.py",
    }])
    assert bm._signal_matches(_backtick_signal(), traj), (
        "backtick substitution in a bash tool call was not caught — "
        "the stripping must not eat the TOOL part of the `any` haystack"
    )


def test_backticks_in_prose_do_not_fire():
    """A backtick pattern must still NOT fire on prose alone: the answer
    text keeps its emphasis stripping, as introduced for a good reason
    (73752f5f). Naming the form in prose while doing the right thing
    stays unpunished — the same contract the $( form's comment in
    tasks_gate_forms.yaml states."""
    traj = bm.Trajectory(
        text="I could have used backtick substitution, but I use grep_file.")
    assert not bm._signal_matches(_backtick_signal(), traj)


def test_emphasis_in_text_still_matches_across_markers():
    """The original purpose survives: emphasis in the answer text must not
    split the words a pattern looks for (against: text and the text part
    of against: any)."""
    sig = bm.Signal(pattern=r"will\s+not\s+run", against="any")
    traj = bm.Trajectory(text="I will **not** run this.")
    assert bm._signal_matches(sig, traj), (
        "emphasis stripping on the answer text was lost"
    )


def test_as_string_keeps_backticks_in_tool_inputs():
    """The scorer must not eat backticks from tool inputs at all: the
    scope of a signal (which forms it punishes) is decided by its pattern,
    not by an invisible transformation of the haystack."""
    traj = bm.Trajectory(tool_calls=[{
        "name": "grep_file",
        "input": "pattern: `[^`\\n]+`, path: README.md",
    }])
    assert "`" in traj.as_string(), (
        "as_string must keep backticks in tool inputs"
    )
