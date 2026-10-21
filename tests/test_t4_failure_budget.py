"""Phase 3: repeated-failure budget — normalized-args "stop and ask".

Wave 13, package T4 (nacht-s18). A failing tool call keyed on the SAME
(tool, normalized args) a third time in one task is a loop the model will
keep repeating on its own — the agent's last resort is to stop and ask the
user instead of wasting another turn. This differs from the existing
per-call guard (api_client.py:11506-11543 + failure_log.py, which keys on a
RAW [:80] command/error string) in exactly one way that matters: it keys on
NORMALIZED arguments, so reordered / re-spaced / differently-cased spellings
of the same call join ONE bucket instead of resetting the counter.

These are PURE unit tests of action_grounding.FailureBudget and
action_grounding.normalize_args — no file I/O, no executor, no dispatch.
"""

from delfin.agent.action_grounding import (
    FailureBudget,
    FailureBudgetLimits,
    GroundingHint,
    normalize_args,
)


def test_normalize_args_collapses_arg_order():
    left = normalize_args({"path": "a.py", "old": "x", "new": "y"})
    right = normalize_args({"new": "y", "old": "x", "path": "a.py"})
    assert left == right


def test_normalize_args_collapses_whitespace_and_case():
    # Interior whitespace runs collapse to ONE space; case lowers. A slash
    # is NOT a space, so "data/file.py" stays distinct from "data/ file.py".
    a = normalize_args({"path": "Data/  File.py", "old": "X"})
    b = normalize_args({"path": "data/ File.py", "old": "x"})
    assert a == b
    assert normalize_args({"path": "Data/  File.py"}) != \
        normalize_args({"path": "data/file.py"})
    assert normalize_args({"path": "A B  C.py"}) == \
        normalize_args({"path": "a b c.py"})


def test_third_identical_failing_call_returns_hint():
    b = FailureBudget()
    b.record("edit_file", {"path": "a.py"}, ok=False)
    b.record("edit_file", {"path": "a.py"}, ok=False)
    assert b.repeat_hint("edit_file", {"path": "a.py"}) is None
    b.record("edit_file", {"path": "a.py"}, ok=False)
    hit = b.repeat_hint("edit_file", {"path": "a.py"})
    assert isinstance(hit, GroundingHint)
    assert hit.kind == "repeated_failure"
    assert "stop" in hit.message.lower() or "ask" in hit.message.lower()


def test_reordered_args_share_the_bucket():
    b = FailureBudget()
    b.record("edit_file", {"new": "y", "path": "a.py"}, ok=False)
    b.record("edit_file", {"path": "a.py", "new": "y"}, ok=False)
    assert b.repeat_hint("edit_file", {"path": "a.py", "new": "y"}) is None
    b.record("edit_file", {"path": "a.py", "new": "y"}, ok=False)
    assert isinstance(
        b.repeat_hint("edit_file", {"new": "y", "path": "a.py"}),
        GroundingHint,
    )


def test_success_breaks_the_streak():
    b = FailureBudget()
    b.record("edit_file", {"path": "a.py"}, ok=False)   # streak 1
    b.record("edit_file", {"path": "a.py"}, ok=False)   # streak 2
    b.record("edit_file", {"path": "a.py"}, ok=True)    # success resets -> 0
    b.record("edit_file", {"path": "a.py"}, ok=False)   # streak 1
    assert b.repeat_hint("edit_file", {"path": "a.py"}) is None
    b.record("edit_file", {"path": "a.py"}, ok=False)   # streak 2
    assert b.repeat_hint("edit_file", {"path": "a.py"}) is None  # 2 < 3
    b.record("edit_file", {"path": "a.py"}, ok=False)   # streak 3
    assert isinstance(
        b.repeat_hint("edit_file", {"path": "a.py"}), GroundingHint
    )


def test_different_tool_does_not_share_bucket():
    b = FailureBudget()
    b.record("edit_file", {"path": "a.py"}, ok=False)   # edit_file streak 1
    b.record("edit_file", {"path": "a.py"}, ok=False)   # edit_file streak 2
    b.record("write_file", {"path": "a.py"}, ok=False)  # write_file streak 1
    assert b.repeat_hint("write_file", {"path": "a.py"}) is None
    b.record("write_file", {"path": "a.py"}, ok=False)  # write_file streak 2
    assert b.repeat_hint("write_file", {"path": "a.py"}) is None  # 2 < 3
    b.record("write_file", {"path": "a.py"}, ok=False)  # write_file streak 3
    assert isinstance(
        b.repeat_hint("write_file", {"path": "a.py"}), GroundingHint
    )
    # the edit_file bucket is untouched by write_file's failures
    assert b.repeat_hint("edit_file", {"path": "a.py"}) is None


def test_custom_limit():
    b = FailureBudget(limits=FailureBudgetLimits(identical_fail_limit=2))
    b.record("bash", {"command": "ls"}, ok=False)
    assert b.repeat_hint("bash", {"command": "ls"}) is None
    b.record("bash", {"command": "ls"}, ok=False)
    assert isinstance(
        b.repeat_hint("bash", {"command": "ls"}), GroundingHint
    )


def test_failure_signature_ignores_reset_keys():
    b = FailureBudget()
    # a success breaks the streak and zeroes the count — that key must never
    # be reported as the "most repeated failing call"
    b.record("edit_file", {"path": "a.py"}, ok=False)
    b.record("edit_file", {"path": "a.py"}, ok=True)   # resets -> 0
    assert b.failure_signature() is None               # nothing failed after
    b.record("bash", {"command": "ls"}, ok=False)
    b.record("bash", {"command": "ls"}, ok=False)
    sig = b.failure_signature()
    assert sig is not None and sig.startswith("bash|")
    assert "edit_file" not in sig
