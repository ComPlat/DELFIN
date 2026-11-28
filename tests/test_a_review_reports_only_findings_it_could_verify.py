"""A review reports findings whose place the diff has, and says which it checked.

Input:  a throwaway git repository with one commit beyond a local origin/main,
        a scripted reviewer (``land.kit_chat``) whose reply each test chooses,
        and every path the gate writes redirected into ``tmp_path``.
Output: what the review flow reports -- which findings reach the refusal line,
        which are held apart, what ``scratch/REVIEW.json`` records, and whether
        the branch was held or landed.
Semantics: the reviewer is shown the diff and nothing else, so a finding about
        a place that diff REFUTES is a claim, not a finding. The three cases
        that matter are told apart here -- the first two are refuted and held
        apart from the reported findings, the third is merely not displayed and
        is kept:

        * a file the diff never touches                 -> held apart
        * a line past the end of a file the branch adds -> held apart
        * a line outside every hunk of a MODIFIED file  -> reported, labelled

        The last is reported because the shipped calibration corpus labels a
        catch of that shape correct (``pos-engine-missing-method``: the finding
        names a symbol the diff adds while the first hunk starts at line 3302),
        so treating "not displayed" as "refuted" would delete true findings.

Every test drives ``land.review_gate`` or ``land.main`` -- the verb
``delfin-agent swarm wave land`` reaches -- not a helper. Nothing reaches the
network, a real repository or ``~/.delfin``.
"""

import json
import os
import subprocess
import time

import pytest

from delfin.agent.swarm import land as pr_bot

PASS = json.dumps({"verdict": "pass", "findings": []})

#: The file the fixture's branch MODIFIES (one line, one hunk) and the one it
#: ADDS (three lines). Line 4000 of each is the same claim shape with opposite
#: answers, because the two diffs prove different things about it.
MODIFIED = "base.txt"
ADDED = "tests/test_new.py"


def _git(*args, cwd):
    r = subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True)
    assert r.returncode == 0, (args, r.stderr)
    return r.stdout


def _refuse(file, line=1, claim="nothing reaches the new path"):
    return json.dumps({"verdict": "refuse", "findings": [
        {"file": file, "line": line, "claim": claim,
         "why_it_matters": "the capability never runs"}]})


@pytest.fixture
def gate(tmp_path, monkeypatch):
    """A repository whose branch modifies one file and adds another, plus a
    calibrated roster and a scripted reviewer."""
    wt = tmp_path / "wt"
    (wt / "scratch").mkdir(parents=True)
    (wt / "scratch" / ".gitignore").write_text("*\n")
    cfg = tmp_path / "gitconfig"
    cfg.write_text("[user]\n\tname = tester\n\temail = tester@example.org\n")
    for k, v in {"GIT_CONFIG_GLOBAL": str(cfg), "HOME": str(tmp_path),
                 "PATH": os.environ["PATH"]}.items():
        monkeypatch.setenv(k, v)
    _git("init", "-q", "-b", "main", cwd=wt)
    (wt / "base.txt").write_text("base\n")
    _git("add", "-A", cwd=wt)
    _git("commit", "-qm", "base", cwd=wt)
    _git("update-ref", "refs/remotes/origin/main", "HEAD", cwd=wt)
    # The branch: one added file (so its whole extent is in the diff) and one
    # modified file (so the diff shows only one line of it).
    (wt / "tests").mkdir()
    (wt / "tests" / "test_new.py").write_text(
        "# marker\ndef test_x():\n    assert True\n")
    (wt / "base.txt").write_text("changed\n")
    (wt / "scratch" / "PR.md").write_text("# A title\n\nbody\n")
    _git("add", "-A", cwd=wt)
    _git("commit", "-qm", "add a test and touch base", cwd=wt)

    log = tmp_path / "n9a.log"
    log.write_text('{"type": "notice", "text": "hello"}\n')
    agent = {"worktree": str(wt), "log": str(log), "branch": "swarm/night9/n9a",
             "task": {"id": "n9a", "model": "alias-code",
                      "brief": "PACKAGE demo: add a helper and its test.\n"}}

    roster = ("kit.glm-5.3-flash", "kit.deepseek-v4-flash",
              "kit.mistral-small-4-119b-a8b")
    monkeypatch.setattr(pr_bot, "ROOT", tmp_path)
    monkeypatch.setattr(pr_bot, "REVIEW_CALIBRATION", tmp_path / "calibration.json")
    monkeypatch.setattr(pr_bot, "REVIEW_LOG", tmp_path / "review_log.jsonl")
    monkeypatch.setattr(pr_bot, "kit_reviewer_models", lambda: roster)
    monkeypatch.setattr(pr_bot, "identity_email", lambda: "tester@example.org")
    rec = {"ts": time.time(), "iso": "2026-01-01T00:00:00+0000",
           "corpus_digest": pr_bot.corpus_digest(), "passed": True, "cases": []}
    cache = {m: dict(rec) for m in roster}
    bounds = {"n": 9, "decided": 5, "wrong": 0, "wrong_cases": [],
              "abstained": 4, "abstained_split": 3, "abstained_no_answer": 1,
              "agreed_pass": 4, "agreed_refuse": 1,
              "correct_when_agreed_at_least": 0.549}
    for i, a in enumerate(roster):
        for b in roster[i + 1:]:
            if pr_bot.model_family(a) == pr_bot.model_family(b):
                continue
            cache[pr_bot.pair_calibration_key(a, b)] = dict(
                rec, models=sorted((a, b)), bounds=bounds)
    (tmp_path / "calibration.json").write_text(json.dumps(cache))

    sent = []
    script = {"out": []}

    def fake(model, messages, timeout=None):
        sent.append({"model": model, "messages": messages})
        if not script["out"]:
            raise AssertionError(f"reviewer {model} called more often than scripted")
        nxt = script["out"].pop(0)
        if isinstance(nxt, Exception):
            raise nxt
        return nxt

    monkeypatch.setattr(pr_bot, "kit_chat", fake)
    return {"wt": wt, "agent": agent, "sent": sent, "root": tmp_path,
            "set": lambda *o: script.update(out=list(o))}


# --- the case the fix is about: a place the diff refutes -------------------


def test_a_finding_about_a_file_the_diff_never_touched_is_not_reported_as_one(gate):
    """RED ON THE BASE: the phantom path was printed as the reason for a refusal.

    A refusal the branch still deserves to be held for -- but not one that
    sends an author to a file they never changed, and not a finding.
    """
    gate["set"](_refuse(file="delfin/agent/never_in_this_branch.py", line=42),
                _refuse(file="delfin/agent/never_in_this_branch.py", line=42))
    ok, why, rec = pr_bot.review_gate(gate["wt"], "night9", "n9a", gate["agent"])
    assert not ok, "an unverifiable refusal must still hold the branch"
    assert rec["verdict"] == "unverified", rec["verdict"]
    assert rec["findings"] == [], "nothing verifiable was found"
    assert len(rec["unconfirmed"]) == 2, "both reviewers' claims are kept, apart"
    assert "refused by BOTH" not in why, why
    assert "never_in_this_branch.py" in why, "the line names what could not be checked"


def test_the_record_and_the_line_agree_when_nothing_could_be_verified(gate):
    """``scratch/REVIEW.json`` cannot say "refuse" beside an UNDECIDED line.

    Red on the base twice over: the file recorded "refuse" with the phantom as
    a finding, and it was written before the verdict was decided.
    """
    gate["set"](_refuse(file="phantom.py", line=7), _refuse(file="phantom.py", line=7))
    ok, why, rec = pr_bot.review_gate(gate["wt"], "night9", "n9a", gate["agent"])
    assert not ok
    on_disk = json.loads((gate["wt"] / "scratch" / "REVIEW.json").read_text())
    assert on_disk["verdict"] == rec["verdict"] == "unverified"
    assert on_disk["findings"] == []
    assert on_disk["unconfirmed"], "the claim is still there for a person to read"


def test_a_line_past_the_end_of_a_file_the_branch_adds_is_held_apart(gate):
    """The diff shows an ADDED file in full, so it can refute a line number.

    ``tests/test_new.py`` is three new lines; a finding at line 4000 of it
    contradicts the diff, and is held apart for the same reason a phantom file
    is -- nobody can go and look.
    """
    gate["set"](_refuse(file=ADDED, line=4000), _refuse(file=ADDED, line=4000))
    ok, why, rec = pr_bot.review_gate(gate["wt"], "night9", "n9a", gate["agent"])
    assert not ok
    assert rec["findings"] == []
    assert "ending at line 3" in rec["unconfirmed"][0]["why_unconfirmed"]
    assert "refused by BOTH" not in why


# --- what the fix must NOT drop ------------------------------------------


def test_a_finding_at_a_place_the_diff_shows_is_still_a_refusal(gate):
    """CONTROL, green on both trees: a real place still refuses the branch.

    Without this, the tests above would pass on a gate that refused everything
    -- and the run's own review log would look diligent while stopping sound
    work. Asserts only what holds before the change as well as after it.
    """
    gate["set"](_refuse(file=ADDED, line=2), _refuse(file=ADDED, line=2))
    ok, why, rec = pr_bot.review_gate(gate["wt"], "night9", "n9a", gate["agent"])
    assert not ok and rec["verdict"] == "refuse"
    assert "refused by BOTH" in why
    assert f"{ADDED}:2" in why


def test_a_confirmed_place_says_so_on_the_line(gate):
    """RED ON THE BASE: the refusal line states what it checked.

    Before the change the line could not tell a verified place from an
    unverified one, because it never looked; there was no such verdict to read.
    """
    gate["set"](_refuse(file=ADDED, line=2), _refuse(file=ADDED, line=2))
    ok, why, rec = pr_bot.review_gate(gate["wt"], "night9", "n9a", gate["agent"])
    # .get(), not [...]: on the base the keys are absent, and the failure a
    # reader should see is "nothing was recorded as confirmed", not a KeyError.
    assert rec["findings"][0].get("place_confirmed") is True
    assert rec.get("unconfirmed") == []
    assert "place confirmed in the diff" in why


def test_a_line_outside_every_hunk_of_a_modified_file_is_kept_and_labelled(gate):
    """The line the fix must NOT delete, and the reason it says so.

    ``base.txt`` is modified, so the diff shows one line of it and cannot say
    whether line 4000 exists. The shipped corpus labels a catch of this shape
    correct, so the finding is reported -- with the reason it could not be
    confirmed attached, so nobody reads it as a checked one.
    """
    gate["set"](_refuse(file=MODIFIED, line=4000), _refuse(file=MODIFIED, line=4000))
    ok, why, rec = pr_bot.review_gate(gate["wt"], "night9", "n9a", gate["agent"])
    assert not ok and rec["verdict"] == "refuse", "not displayed is not refuted"
    assert rec["findings"], "the finding is reported, not deleted"
    assert rec["findings"][0]["place_confirmed"] is False
    assert "outside every hunk" in rec["findings"][0]["why_unconfirmed"]
    assert "outside its hunks" in why


def test_one_verifiable_and_one_phantom_find_report_the_first_and_say_so(gate):
    """A mixed answer: the branch is refused on what was checked, and the count
    held apart is in the line rather than quietly dropped."""
    gate["set"](json.dumps({"verdict": "refuse", "findings": [
                   {"file": "phantom.py", "line": 9, "claim": "a guard is gone",
                    "why_it_matters": "a check is lost"},
                   {"file": ADDED, "line": 1, "claim": "nothing calls it",
                    "why_it_matters": "it never runs"}]}),
               _refuse(file=ADDED, line=1))
    ok, why, rec = pr_bot.review_gate(gate["wt"], "night9", "n9a", gate["agent"])
    assert not ok and rec["verdict"] == "refuse"
    assert [f["file"] for f in rec["findings"]] == [ADDED, ADDED]
    assert [f["file"] for f in rec["unconfirmed"]] == ["phantom.py"]
    assert "held apart" in why


def test_a_clean_review_still_lands_the_branch(gate):
    """CONTROL, green on both trees: a pass is a pass.

    The other control the fix could fake its way past: if a verification step
    simply refused to confirm anything, the branch would be held forever and
    every run would stall behind a gate that never lets anything through.
    """
    gate["set"](PASS, PASS)
    ok, why, rec = pr_bot.review_gate(gate["wt"], "night9", "n9a", gate["agent"])
    assert ok, why
    assert rec["verdict"] == "pass"


def test_a_clean_review_reports_nothing_held_apart(gate):
    """RED ON THE BASE: a pass records that nothing was held apart.

    The absence is the point: with no verification step there was no held-apart
    set to report, so this assertion cannot be satisfied by a gate that simply
    never checks.
    """
    gate["set"](PASS, PASS)
    ok, why, rec = pr_bot.review_gate(gate["wt"], "night9", "n9a", gate["agent"])
    assert ok, why
    assert rec.get("unconfirmed") == []


# --- through the verb the product actually runs ---------------------------


def test_the_verb_reports_the_unverifiable_refusal_as_held_not_refused(gate,
                                                                       monkeypatch,
                                                                       capsys):
    """``delfin-agent swarm wave land`` reaches the same decision.

    ``check()`` is the branch-test gate ABOVE the review, and it runs pytest on
    a real repository; this test is about what the review reports, so it is
    replaced with a pass, as the sibling suite does for the same reason.
    """
    monkeypatch.setattr(pr_bot, "check", lambda wt, base="origin/main", **_kw: (True, "ok"))
    gate["set"](_refuse(file="delfin/agent/phantom.py", line=3),
                _refuse(file="delfin/agent/phantom.py", line=3))
    run = gate["root"] / "night9"
    run.mkdir()
    (run / "state.json").write_text(json.dumps({"agents": {"n9a": gate["agent"]}}))
    pr_bot.main([run.name, "--dry-run"])
    out = capsys.readouterr().out
    assert "review REFUSED" in out
    txt = (gate["wt"] / "scratch" / "LANDING.txt").read_text()
    assert txt.startswith("refused: independent review UNDECIDED"), txt
    assert "would push" not in out, "an unverifiable refusal opens nothing"
    assert "delfin/agent/phantom.py" in txt


# --- the oracle itself ----------------------------------------------------


def test_the_diff_reader_sees_both_sides_and_the_files_a_branch_adds():
    """The three diff shapes the rule turns on, read from real unified diffs.

    A hunk shows context and removed lines as well as added ones, so a finding
    about either is about a place the reviewer was shown; an added file is
    shown whole, which is what lets its length refute a line; and a deleted
    file has no ``+++`` path, so reading only that side loses every true
    finding about a removal.
    """
    diff = ("diff --git a/gone.py b/gone.py\ndeleted file mode 100644\n"
            "--- a/gone.py\n+++ /dev/null\n@@ -1,3 +0,0 @@\n-a\n-b\n-c\n"
            "diff --git a/new.py b/new.py\nnew file mode 100644\n"
            "--- /dev/null\n+++ b/new.py\n@@ -0,0 +1,2 @@\n+x\n+y\n"
            "diff --git a/mid.py b/mid.py\n--- a/mid.py\n+++ b/mid.py\n"
            "@@ -10,3 +10,3 @@\n ctx\n-rem\n+add\n")
    places = pr_bot.diff_places(diff)
    assert places.for_file("gone.py") == [(1, 3)], "a deletion keeps its path"
    assert places.for_file("new.py") == [(1, 2)]
    assert places.whole == frozenset({"new.py"}), "only the added file is whole"
    assert places.for_file("mid.py") == [(10, 12), (10, 12)], "both sides count"
    assert places.for_file("nowhere.py") == []


def test_an_unreadable_hunk_header_leads_to_reporting_less_not_more():
    """A diff the reader cannot parse yields no places, so nothing is CONFIRMED.

    The direction matters: failing to read a header must not turn a phantom
    into a verified finding.
    """
    places = pr_bot.diff_places("--- a/x.py\n+++ b/x.py\n@@ garbage @@\n+z\n")
    assert places.for_file("x.py") == []
    reported, held = pr_bot.verified_findings(
        {"findings": [{"file": "x.py", "line": 1, "claim": "c",
                       "why_it_matters": "w"}]},
        "--- a/x.py\n+++ b/x.py\n@@ garbage @@\n+z\n")
    assert reported == [] and len(held) == 1
