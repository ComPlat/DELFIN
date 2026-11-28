"""A review reports a finding only when it can name where and how it fails.

The defect this pins was measured on the untouched tree: a review of a diff
was answered with a checklist, in whatever order the delegate wrote it, and
nothing required a finding to carry a location, a way to reproduce the
failure, or the reviewer's own word that it had checked. A finding the
reviewer invented was indistinguishable from one it had read the code for.

The contract now refuses the second kind and COUNTS it, so the count cannot
be silent. These tests judge the outcome, never the wording:

* a finding citing a real file:line, with a scenario, marked verified, is
  reported — with the cited line's own text, so the parent can see what the
  reviewer pointed at;
* each way a finding can fail to be checkable withholds THAT finding and
  nothing else, and each reason is named in the withheld list;
* the survivors are ordered most-severe first;
* an empty review and an unreadable contract do NOT look alike.

Driven twice: once as pure functions (the shapes), and once through the
product path — ``cli.main`` with the model's stream scripted — because a
contract that the ``subagent`` tool never carries to the parent is a function
that exists and is never called.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import subagents as SA

_REAL_FILE = "pkg/inventory.py"
#: Seven lines, and the module constant the fixture findings cite is line 3 of
#: them. The assertions below name both numbers, so a change here has to move
#: them too.
_SOURCE = '''\
"""A module with a real line to cite."""

DEFAULT_LIMIT = 10


def read_rows(path, limit=DEFAULT_LIMIT):
    return []
'''


@pytest.fixture
def workspace(tmp_path):
    target = tmp_path / _REAL_FILE
    target.parent.mkdir(parents=True)
    target.write_text(_SOURCE, encoding="utf-8")
    return tmp_path


def _finding(**over) -> dict:
    base = {
        "severity": "major",
        "file": _REAL_FILE,
        "line": 3,
        "what": "the default limit is shared mutable state",
        "failure_scenario": "read_rows() called after a caller sets "
                            "DEFAULT_LIMIT to 0 returns [] instead of rows",
        "verified": True,
    }
    base.update(over)
    return base


def _review(*findings, root) -> dict:
    return SA.verify_review_findings({"findings": list(findings)},
                                     repo_root=root)


# ---------------------------------------------------------------------------
# What survives, and what it carries
# ---------------------------------------------------------------------------

def test_a_verified_finding_is_reported_with_its_line_and_scenario(workspace):
    out = _review(_finding(), root=workspace)
    assert out["status"] == "ok"
    [kept] = out["findings"]
    assert kept["file"] == _REAL_FILE
    assert kept["line"] == 3
    assert kept["failure_scenario"]
    assert kept["verification"] == "verified"
    # Read from the FILE, not from the finding: this is what makes the
    # citation checkable by the reader rather than restated by the writer.
    assert kept["source_line"] == "DEFAULT_LIMIT = 10"


@pytest.mark.parametrize("over, reason", [
    ({"file": "pkg/does_not_exist.py"}, "does not exist"),
    ({"line": 9000}, "past the end of the file"),
    ({"line": 0}, "past the end of the file"),
    ({"line": "4"}, "line is not a number"),
    ({"line": True}, "line is not a number"),
    ({"file": ""}, "no file"),
    ({"failure_scenario": "   "}, "no failure scenario"),
    ({"verified": False}, "the reviewer did not verify it"),
    ({"verified": "yes"}, "the reviewer did not verify it"),
    ({"file": "../../etc/passwd"}, "outside the workspace"),
])
def test_a_finding_that_cannot_be_checked_is_withheld_and_named(
        workspace, over, reason):
    out = _review(_finding(), _finding(**over), root=workspace)
    assert out["status"] == "ok"
    assert len(out["findings"]) == 1, "the good one must survive"
    assert len(out["dropped"]) == 1
    assert reason in out["dropped"][0]["reason"]
    # Withholding is never silent: the parent is told how many were withheld.
    assert "1 finding(s) withheld" in out["notice"]


def test_a_review_of_only_unverifiable_findings_reports_none(workspace):
    out = _review(_finding(line=9000), _finding(verified=False),
                  root=workspace)
    assert out["findings"] == []
    assert out["status"] == "all_withheld"
    assert out["counts"] == {"kept": 0, "dropped": 2, "by_severity": {}}


def test_survivors_are_ranked_most_severe_first(workspace):
    out = _review(_finding(severity="minor"), _finding(severity="major"),
                  _finding(severity="critical"), root=workspace)
    assert [f["severity"] for f in out["findings"]] == [
        "critical", "major", "minor"]


def test_an_unknown_severity_keeps_the_finding_and_sorts_it_last(workspace):
    """A reviewer that wrote "blocker" found something. Dropping the finding
    to punish the vocabulary would lose the finding; only the ranking of it
    is unknown, so it goes last."""
    out = _review(_finding(severity="blocker"), _finding(severity="minor"),
                  root=workspace)
    assert [f["severity"] for f in out["findings"]] == ["minor", "blocker"]
    assert out["counts"]["dropped"] == 0


def test_no_findings_is_not_the_same_as_unreadable(workspace):
    """'the reviewer looked and found nothing' and 'nobody knows' are
    different answers, and a caller that cannot tell them apart will read the
    second as the first."""
    empty = _review(root=workspace)
    assert empty["status"] == "no_findings"
    assert empty["findings"] == []
    unreadable = SA.verify_review_findings({"nope": 1}, repo_root=workspace)
    assert unreadable["status"] == "unreadable"


def test_a_finding_is_checked_against_the_workspace_it_names(workspace):
    """No workspace means nothing can be checked, and an unchecked citation
    must not be reported as one that was."""
    out = SA.verify_review_findings({"findings": [_finding()]},
                                    repo_root=None)
    assert out["findings"] == []
    assert "no workspace" in out["dropped"][0]["reason"]


def test_a_non_finding_entry_costs_only_itself(workspace):
    out = _review(_finding(), "not an object", _finding(line=9000),
                  root=workspace)
    assert len(out["findings"]) == 1
    assert [d["reason"] for d in out["dropped"]] == [
        "not an object", "pkg/inventory.py:9000 is past the end of the file "
        "(7 lines)"]


# ---------------------------------------------------------------------------
# Only a review gets the contract
# ---------------------------------------------------------------------------

def test_another_presets_payload_is_not_touched(workspace):
    """An explore report and a security verdict are not finding lists. Giving
    them one would withhold their whole answer."""
    payload = {"subagent_type": "explore", "structured_output": {"a": 1},
               "result": "found it"}
    assert SA.attach_review_findings(payload, repo_root=workspace) == payload


def test_the_review_contract_is_attached_and_idempotent(workspace):
    payload = {"subagent_type": "code-reviewer",
               "structured_output": {"findings": [_finding()]},
               "result": "{}"}
    once = SA.attach_review_findings(payload, repo_root=workspace)
    assert once["review"]["counts"]["kept"] == 1
    assert once["result"] == "{}", "the delegate's own text is not rewritten"
    twice = SA.attach_review_findings(once, repo_root=workspace)
    assert twice == once


# ---------------------------------------------------------------------------
# Through the product: the subagent tool the model calls
# ---------------------------------------------------------------------------

def test_the_review_preset_ships_the_contract():
    preset = SA.SUBAGENT_PRESETS["code-reviewer"]
    assert preset.output_schema is SA.REVIEW_FINDINGS_SCHEMA
    required = (SA.REVIEW_FINDINGS_SCHEMA["properties"]["findings"]["items"]
                ["required"])
    for field in ("file", "line", "failure_scenario", "verified"):
        assert field in required, (
            f"a finding without {field} is what made a review uncheckable")
    # ...and the prompt asks for what the schema requires, so a model that
    # only reads the prose is told the same thing.
    prompt = preset.system_prompt
    assert "failure_scenario" in prompt and "verified" in prompt
    assert "NOTHING you did not verify" in prompt


def test_the_slash_command_names_the_answer_form():
    """`/review` expanded to the same generic briefing as `/explore`, so the
    one thing a review must carry was left to the model's discretion."""
    from delfin.agent import slash_commands

    text = slash_commands.builtin_subagent_command(
        "review", "the diff on this branch")
    assert "code-reviewer" in text
    for phrase in ("FINDINGS", "failure scenario", "verified"):
        assert phrase in text, phrase
    # The other builtins are untouched by this: an explore task is not a
    # review and must not be given a findings list.
    assert "failure scenario" not in slash_commands.builtin_subagent_command(
        "explore", "find the callers of parse_control")


def test_an_unverified_finding_never_reaches_the_parent(tmp_path):
    """End to end through ``run_subagent``: the delegate's structured answer
    holds one checkable finding and one it invented, and the payload the
    parent reads holds only the first — with the second COUNTED.

    Scripted at the client's stream, not at ``verify_review_findings``: the
    question is whether ``run_subagent`` carries the contract to the parent,
    and a test that calls the filter directly answers a different one.
    """
    from delfin.agent import api_client as A
    from delfin.agent.api_client import StreamEvent

    ws = tmp_path / "ws"
    (ws / _REAL_FILE).parent.mkdir(parents=True)
    (ws / _REAL_FILE).write_text(_SOURCE, encoding="utf-8")

    good = _finding()
    invented = _finding(line=9000, failure_scenario="I think this is wrong")

    class _Client:
        """The minimum a delegate needs: one final message, no tool calls."""

        model = "fake-model"

        def set_permissions(self, perms):
            pass

        def stream_message(self, *, messages, system, max_tokens, **kw):
            yield StreamEvent(
                type="text_delta",
                text=json.dumps({"findings": [good, invented]}))
            yield StreamEvent(type="message_delta", input_tokens=1,
                              output_tokens=1)

    perms = A.KitToolPermissions(mode="plan", workspace=str(ws))
    result = SA.run_subagent(
        subagent_type="code-reviewer",
        description="review the diff",
        prompt="review it",
        parent_client=_Client(),
        parent_perms=perms,
    )
    payload = result.to_payload()
    review = payload["review"]
    assert review["counts"] == {"kept": 1, "dropped": 1,
                                "by_severity": {"major": 1}}
    assert review["dropped"][0]["reason"].endswith(
        "past the end of the file (7 lines)")
    assert payload["review_notice"].startswith("[review-findings] 1 finding")
    assert [f["line"] for f in review["findings"]] == [3]


# ---------------------------------------------------------------------------
# The control: this file is not blanket-red on the base
# ---------------------------------------------------------------------------

def test_the_review_preset_stays_read_only_and_undelegating():
    """Passes on BOTH revisions, deliberately.

    Without it, "the file is red on the base and green on the branch" would
    also be satisfied by a file that fails on the base for an unrelated
    reason. This is the part of the preset that must not move: a reviewer
    that could write, or that could spawn its own reviewers, would be a
    different capability than the one being fixed.
    """
    preset = SA.SUBAGENT_PRESETS["code-reviewer"]
    assert preset.mode == "plan"
    assert preset.may_delegate_to == ()
    assert "read_file" in preset.tools and "grep_file" in preset.tools
    for writer in ("write_file", "edit_file", "multi_edit", "apply_patch"):
        assert writer not in preset.tools
