"""The PR bot's independent review gate, checked without calling a model.

Input:  a throwaway git worktree with one commit beyond a local origin/main,
        a run state that names it, and a stubbed reviewer (``kit_chat``) whose
        reply each test chooses.
Output: what the gate does -- which branches it refuses, what reaches the
        reviewer's prompt, what is recorded, and whether a push was attempted.
Semantics: every test here names one check and fails when that check is
        removed from pr_bot.py. No test reaches the network, the real
        repository, ~/.delfin or any live store: ``sh`` is wrapped so a push
        or a ``gh`` call fails the test instead of happening, and every path
        the gate writes is redirected into tmp_path.

Ported with the gate (delfin/agent/swarm/land.py) from the orchestrator's
own landing script; the checks are the same.
"""

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

from delfin.agent.swarm import land as pr_bot


DIFF_MARK = "a_sentinel_line_only_the_diff_has"
BRIEF_MARK = "a_sentinel_line_only_the_brief_has"
PR_MD_MARK = "a_sentinel_line_only_the_authors_pr_md_has"
PR_BOT_MARK = "a_sentinel_line_only_pr_bot_txt_has"
LOG_MARK = "a_sentinel_line_only_the_agent_log_has"


def _git(*args, cwd):
    r = subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout


@pytest.fixture
def swarm(tmp_path, monkeypatch):
    """A run directory, a worktree with one commit, and every path redirected.

    The worktree is a real git repository with a local ``origin/main`` ref, so
    the gate's own ``git diff origin/main...HEAD`` runs for real.
    """
    wt = tmp_path / "wt"
    (wt / "scratch").mkdir(parents=True)
    env = {"GIT_AUTHOR_NAME": "tester", "GIT_AUTHOR_EMAIL": "tester@example.org",
           "GIT_COMMITTER_NAME": "tester", "GIT_COMMITTER_EMAIL": "tester@example.org",
           "GIT_CONFIG_GLOBAL": str(tmp_path / "gitconfig"), "HOME": str(tmp_path),
           "PATH": os.environ["PATH"]}
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    _git("init", "-q", "-b", "main", cwd=wt)
    (wt / "base.txt").write_text("base\n")
    (wt / "scratch" / ".gitignore").write_text("*\n")
    _git("add", "-A", cwd=wt)
    _git("commit", "-qm", "base", cwd=wt)
    _git("update-ref", "refs/remotes/origin/main", "HEAD", cwd=wt)
    (wt / "tests").mkdir()
    (wt / "tests" / "test_new.py").write_text(f"# {DIFF_MARK}\ndef test_x():\n    assert True\n")
    _git("add", "-A", cwd=wt)
    _git("commit", "-qm", "add a test", cwd=wt)

    (wt / "scratch" / "PR.md").write_text(f"# a title\n\nbody {PR_MD_MARK}\n")
    (wt / "scratch" / "LANDING.txt").write_text(f"opened: {PR_BOT_MARK}\n")
    log = tmp_path / "n9a.log"
    log.write_text(f'{{"type": "notice", "text": "{LOG_MARK}"}}\n')

    run = tmp_path / "night9"
    run.mkdir()
    state = {"agents": {"n9a": {
        "worktree": str(wt), "log": str(log), "branch": "swarm/night9/n9a",
        "task": {"id": "n9a", "model": "alias-code",
                 "brief": f"Task n9a -- {BRIEF_MARK}\nscratch/PR.md is the route to a PR.\n"}}}}
    (run / "state.json").write_text(json.dumps(state))

    monkeypatch.setattr(pr_bot, "ROOT", tmp_path)
    monkeypatch.setattr(pr_bot, "identity_email", lambda: "tester@example.org")
    monkeypatch.setattr(pr_bot, "REVIEW_CALIBRATION", tmp_path / "calibration.json")
    monkeypatch.setattr(pr_bot, "REVIEW_LOG", tmp_path / "review_log.jsonl")
    # A reviewer is only used when it is calibrated, and the gate needs a
    # PAIR, so the roster carries THREE families: excluding the author's
    # family must still leave two candidates. With two families it leaves
    # one, no pair can form, and every test would measure that refusal
    # instead of what it was written for.
    roster = ("kit.glm-5.3-flash", "kit.deepseek-v4-flash",
              "kit.mistral-small-4-119b-a8b")
    monkeypatch.setattr(pr_bot, "kit_reviewer_models", lambda: roster)
    _calibrate_for(roster)
    return {"root": tmp_path, "wt": wt, "run": "night9", "tid": "n9a",
            "agent": state["agents"]["n9a"], "log": log}


def _calibrate_for(models, ts=None, digest=None):
    """A passing calibration record, written the way ``calibrate`` writes one.

    Both kinds are written: one record per model, and one per cross-family
    PAIR, because the gate asks for a calibrated pair and a per-model record
    alone no longer opens it. The pair records carry the same ``ts`` and
    ``corpus_digest`` as the model records, so a test that makes a
    calibration stale or puts it on another corpus still measures that
    instead of leaving a fresh pair behind that would open the gate anyway.
    """
    models = list(models)
    rec = {"ts": ts if ts is not None else time.time(), "iso": "2026-11-20T00:00:00+0000",
           "corpus_digest": digest if digest is not None else pr_bot.corpus_digest(),
           "passed": True, "cases": []}
    cache = {m: dict(rec) for m in models}
    bounds = {"n": 9, "decided": 5, "wrong": 0, "wrong_cases": [],
              "abstained": 4, "abstained_split": 3, "abstained_no_answer": 1,
              "agreed_pass": 4, "agreed_refuse": 1,
              "correct_when_agreed_at_least": 0.549}
    for i, a in enumerate(models):
        for b in models[i + 1:]:
            if pr_bot.model_family(a) == pr_bot.model_family(b):
                continue
            cache[pr_bot.pair_calibration_key(a, b)] = dict(
                rec, models=sorted((a, b)), bounds=bounds)
    pr_bot.REVIEW_CALIBRATION.write_text(json.dumps(cache, indent=1))


@pytest.fixture
def no_push(monkeypatch):
    """Every push and every ``gh`` call fails the test instead of happening."""
    real = pr_bot.sh
    attempted = []

    def guard(args, cwd=None, timeout=1800):
        if "push" in args or args[0] == "gh":
            attempted.append(list(args))
            raise AssertionError("the gate let a push/gh call through: " + " ".join(args))
        return real(args, cwd=cwd, timeout=timeout)

    monkeypatch.setattr(pr_bot, "sh", guard)
    return attempted


@pytest.fixture
def replies(monkeypatch):
    """Scripted reviewer replies; records every prompt the reviewer was sent."""
    sent = []
    script = {"out": [], "queues": {}}

    def fake(model, messages, timeout=None):
        # One queue PER MODEL, not one shared queue. The gate asks a pair of
        # reviewers, so a shared queue would hand reviewer A's scripted
        # reply to A and reviewer B's turn would find the queue empty --
        # which the pair reads as "no answer" and reports as undecided. A
        # test that scripts one refusal means "each reviewer refuses", and
        # that is what the pair is then asked to decide.
        sent.append({"model": model, "messages": messages})
        q = script["queues"].setdefault(model, list(script["out"]))
        if not q:
            raise AssertionError(
                f"the reviewer {model} was called more often than scripted")
        nxt = q.pop(0)
        if isinstance(nxt, Exception):
            raise nxt
        return nxt

    def setter(*o):
        script.update(out=list(o), queues={})

    monkeypatch.setattr(pr_bot, "kit_chat", fake)
    return {"sent": sent, "script": script["out"], "set": setter}


def _refuse(file="tests/test_new.py", line=1, claim="nothing calls it"):
    return json.dumps({"verdict": "refuse", "findings": [
        {"file": file, "line": line, "claim": claim, "why_it_matters": "it never runs"}]})


PASS = json.dumps({"verdict": "pass", "findings": []})


# --- 1. in the path, never beside it ---------------------------------------


def test_a_refused_review_leaves_the_pr_unopened(swarm, no_push, replies, capsys):
    """The gate runs inside main(), before the push: a refusal writes
    scratch/LANDING.txt and opens nothing. ``no_push`` fails the test if a push
    or a gh call is reached at all."""
    replies["set"](_refuse())
    monkeypatch_check(swarm)
    pr_bot.main(["--push", swarm["run"]])
    out = capsys.readouterr().out
    assert "review REFUSED" in out
    txt = (swarm["wt"] / "scratch" / "LANDING.txt").read_text()
    assert txt.startswith("refused: independent review refused by BOTH ")
    assert "kit.glm-5.3-flash" in txt and "kit.deepseek-v4-flash" in txt
    assert "tests/test_new.py:1" in txt
    assert (swarm["wt"] / "scratch" / "PR.md").is_file(), "PR.md was consumed by a refusal"


def test_a_passed_review_is_what_lets_the_push_happen(swarm, no_push, replies, capsys):
    """With a pass, main() reaches the push step -- shown here with --dry-run,
    which prints what it would do and calls nothing."""
    replies["set"](PASS)
    monkeypatch_check(swarm)
    pr_bot.main(["--dry-run", swarm["run"]])
    out = capsys.readouterr().out
    assert "review passed" in out
    assert "dry-run: would push" in out and "wave/night9-n9a" in out


def monkeypatch_check(swarm):
    """Replace the branch checks above the gate with a pass, for this test.

    They run tests on delfin-ci against the real repository; this file is
    about the gate, and the fixture's worktree is not that repository. The
    real ``check`` comes back after the test (``_restore_check``): assigned
    for good, the stub answered every later test's ``check`` with a pass.
    """
    pr_bot.check = lambda wt, base="origin/main", **_: (True, "1 passed (stubbed)")


_REAL_CHECK = pr_bot.check


@pytest.fixture(autouse=True)
def _restore_check():
    yield
    pr_bot.check = _REAL_CHECK


# --- 2. independence is mechanical, and the gate fails closed --------------


def test_the_reviewer_is_never_of_the_authors_family(swarm, replies):
    """A DeepSeek author is never reviewed by a DeepSeek model."""
    swarm["agent"]["task"]["model"] = "deepseek-v4.1-flash-at-hzdr"
    swarm["log"].write_text("{}\n")
    replies["set"](PASS)
    ok, why, rec = pr_bot.review_gate(swarm["wt"], "night9", "n9a", swarm["agent"])
    assert ok, why
    author_family = pr_bot.model_family("deepseek-v4.1-flash-at-hzdr")
    assert len(rec["reviewers"]) == 2
    # BOTH reviewers, not just the first: a pair in which one member shares
    # the author's family is not an independent review of it.
    for m in rec["reviewers"]:
        assert pr_bot.model_family(m) != author_family, m


def test_no_model_of_another_family_refuses_rather_than_approves(swarm, replies, monkeypatch):
    """The failure mode this gate exists against: the reviewer is missing and
    the branch is approved anyway. It must refuse."""
    monkeypatch.setattr(pr_bot, "kit_reviewer_models", lambda: ("kit.glm-5.3-flash",))
    swarm["agent"]["task"]["model"] = "kit.glm-5.3-flash"
    swarm["log"].write_text("{}\n")
    ok, why, rec = pr_bot.review_gate(swarm["wt"], "night9", "n9a", swarm["agent"])
    assert not ok
    assert "no independent reviewer available" in why and "glm" in why
    assert replies["sent"] == [], "a reviewer was called although none was independent"


def test_an_author_model_of_an_unknown_family_refuses(swarm, replies):
    """Independence cannot be shown against an author nobody can name."""
    swarm["agent"]["task"]["model"] = "some-new-model-nobody-mapped"
    swarm["log"].write_text("{}\n")
    ok, why, _ = pr_bot.review_gate(swarm["wt"], "night9", "n9a", swarm["agent"])
    assert not ok and "unknown family" in why
    assert replies["sent"] == []


def test_a_resumed_agents_earlier_model_still_counts_as_an_author(swarm, replies):
    """The state keeps only the model of the LAST resume. A family named in the
    agent's own log is an author too, or the reviewer can share it."""
    swarm["agent"]["task"]["model"] = "alias-code"
    swarm["log"].write_text('{"text": "resumed on kit.glm-5.3-flash"}\n')
    replies["set"](PASS)
    ok, why, rec = pr_bot.review_gate(swarm["wt"], "night9", "n9a", swarm["agent"])
    assert ok, why
    assert "kit.glm-5.3-flash" not in rec["reviewers"], (
        "reviewed by the family named in the agent's own log")
    assert len(rec["reviewers"]) == 2


def test_an_uncalibrated_reviewer_refuses(swarm, replies):
    """Reachable is not enough: a reviewer with no passing calibration is not
    used, and the branch is refused rather than passed."""
    pr_bot.REVIEW_CALIBRATION.write_text("{}")
    ok, why, _ = pr_bot.review_gate(swarm["wt"], "night9", "n9a", swarm["agent"])
    assert not ok and "no independent reviewer PAIR available" in why
    assert "--calibrate-pair" in why
    assert replies["sent"] == []


def test_a_stale_calibration_is_not_a_calibration(swarm):
    """What a KIT id serves changes without notice, so a calibration expires."""
    _calibrate_for(("kit.glm-5.3-flash",), ts=time.time() - pr_bot.CALIBRATION_MAX_AGE_S - 60)
    assert pr_bot.calibration_record("kit.glm-5.3-flash") is None
    _calibrate_for(("kit.glm-5.3-flash",))
    assert pr_bot.calibration_record("kit.glm-5.3-flash")


def test_a_calibration_on_another_corpus_does_not_count(swarm):
    """Editing the corpus until a reviewer passes would otherwise leave the old
    pass standing, and the record would be of a corpus that no longer exists."""
    _calibrate_for(("kit.glm-5.3-flash",), digest="a-corpus-that-is-not-this-one")
    assert pr_bot.calibration_record("kit.glm-5.3-flash") is None


def test_an_endpoint_that_does_not_answer_refuses(swarm, replies):
    """Fail closed: two failed calls are a refusal, never a pass."""
    replies["set"](RuntimeError("HTTP 503 no available server"),
                   RuntimeError("HTTP 503 no available server"))
    ok, why, rec = pr_bot.review_gate(swarm["wt"], "night9", "n9a", swarm["agent"])
    assert not ok and "UNDECIDED" in why and "503" in why
    assert rec["verdict"] == "no_answer"
    assert "refused" not in why, ("an endpoint failure is not a verdict about "
                                  "the branch and must not read like one")


def test_only_local_kit_models_can_review(monkeypatch):
    """A reviewer runs on KIT only -- never an external provider, and never an
    id whose family is unknown. The served listing is not asked here (that is
    the network); without it the product's ranking is the roster."""
    from delfin.agent import model_capabilities as MC
    monkeypatch.setattr(MC, "kit_models_served", lambda **kw: set())
    roster = pr_bot.kit_reviewer_models()
    assert roster, "the KIT roster could not be read from the product"
    assert all(m.startswith("kit.") and pr_bot.model_family(m) for m in roster)
    cands, _ = pr_bot.independent_candidates(["alias-code"])
    assert cands and all(m.startswith("kit.") for m in cands)


# --- 3. the reviewer never sees the author's conclusion --------------------


def test_the_prompt_holds_the_diff_and_the_brief_and_nothing_else(swarm, replies):
    """Through the real gate path, with PR.md, LANDING.txt and the log all
    present on disk: none of them reaches the reviewer."""
    replies["set"](PASS)
    ok, why, _ = pr_bot.review_gate(swarm["wt"], "night9", "n9a", swarm["agent"])
    assert ok, why
    sent = "\n".join(m["content"] for m in replies["sent"][0]["messages"])
    assert DIFF_MARK in sent and BRIEF_MARK in sent
    for mark in (PR_MD_MARK, PR_BOT_MARK, LOG_MARK):
        assert mark not in sent, f"the author's own text reached the reviewer: {mark}"
    assert "add a test" not in sent, "the commit message reached the reviewer"


def test_the_prompt_builder_adds_no_third_source(swarm):
    """The builder itself: everything in the prompt comes from the diff, the
    brief, or the fixed instruction."""
    msgs = pr_bot.build_review_prompt("DIFF-" + DIFF_MARK, "BRIEF-" + BRIEF_MARK)
    body = "\n".join(m["content"] for m in msgs)
    rest = body.replace(pr_bot.REVIEW_SYSTEM, "").replace("DIFF-" + DIFF_MARK, "")
    rest = rest.replace("BRIEF-" + BRIEF_MARK, "")
    assert DIFF_MARK not in rest and BRIEF_MARK not in rest
    assert len(rest) < 400, "the prompt carries more than the brief, the diff and the ask"


# --- 5. a refusal must be falsifiable -------------------------------------


def test_a_refusal_without_a_file_and_line_is_not_a_refusal(swarm, replies):
    """"This looks risky" cannot be falsified. It is asked once more, then the
    branch is left for a person -- and it is NOT reported as a refusal."""
    vague = json.dumps({"verdict": "refuse", "findings": [
        {"claim": "this looks risky", "why_it_matters": "it might break something"}]})
    replies["set"](vague, vague)
    ok, why, rec = pr_bot.review_gate(swarm["wt"], "night9", "n9a", swarm["agent"])
    assert not ok
    assert "UNDECIDED" in why and "no usable answer" in why
    assert rec["verdict"] == "no_answer" and rec["verdict"] != "refuse"
    assert len(replies["sent"]) == 4, ("two reviewers, each asked once more "
                                       "after a vague answer")


def test_a_second_answer_with_a_file_and_line_counts(swarm, replies):
    """The retry is a real second chance, not a formality."""
    replies["set"](json.dumps({"verdict": "refuse", "findings": [{"claim": "risky"}]}),
                   _refuse(file="delfin/agent/x.py", line=42, claim="never reached"))
    ok, why, rec = pr_bot.review_gate(swarm["wt"], "night9", "n9a", swarm["agent"])
    assert not ok and rec["verdict"] == "refuse"
    assert "delfin/agent/x.py:42" in why


def test_the_findings_are_left_on_disk_for_a_person(swarm, replies):
    replies["set"](_refuse(file="delfin/agent/x.py", line=7))
    pr_bot.review_gate(swarm["wt"], "night9", "n9a", swarm["agent"])
    rec = json.loads((swarm["wt"] / "scratch" / "REVIEW.json").read_text())
    assert rec["verdict"] == "refuse" and rec["findings"][0]["line"] == 7
    assert rec["reviewers"] == ["kit.glm-5.3-flash", "kit.deepseek-v4-flash"]
    assert rec["author_models"] == ["alias-code"]
    assert rec["a"]["model"] != rec["b"]["model"]


def test_a_line_that_is_not_a_line_is_not_a_finding(swarm):
    """A file with no line, a boolean, a non-numeric string: none of them is a
    place a person can go and look."""
    assert pr_bot.checkable_findings({"findings": [{"file": "a.py", "claim": "c"}]}) == []
    assert pr_bot.checkable_findings({"findings": [{"file": "a.py", "line": True, "claim": "c"}]}) == []
    assert pr_bot.checkable_findings({"findings": [{"file": "", "line": 3, "claim": "c"}]}) == []
    assert pr_bot.checkable_findings({"findings": [{"file": "a.py", "line": 3, "claim": ""}]}) == []
    got = pr_bot.checkable_findings({"findings": [{"file": "a.py", "line": "line 12", "claim": "c"}]})
    assert got and got[0]["line"] == 12


# --- 4. calibrated before it is trusted -----------------------------------


def _fake_reviewer(monkeypatch, answer):
    calls = []

    def fake(model, messages, timeout=None):
        calls.append(messages)
        return answer(messages)

    monkeypatch.setattr(pr_bot, "kit_chat", fake)
    return calls


def _case_of(messages):
    """Which corpus case this reviewer prompt is about, by its diff.

    The prompt carries the diff verbatim, so the case is identified by
    content rather than by a keyword -- a keyword fixture rots as the corpus
    grows, and it cannot tell a diff from the same diff with the fix in it
    (pos-832 and neg-832 differ by one kwarg).
    """
    body = messages[-1]["content"]
    for case in pr_bot.corpus_cases():
        diff = pr_bot._corpus_read(case["diff"]).decode().strip()
        if diff and diff in body:
            return case
    return None


def _calibrate_in(tmp_path, monkeypatch, answer):
    monkeypatch.setattr(pr_bot, "REVIEW_CALIBRATION", tmp_path / "cal.json")
    _fake_reviewer(monkeypatch, answer)
    return pr_bot.calibrate(["kit.test-model"], out=lambda *a: None)["kit.test-model"]


def test_a_reviewer_that_passes_everything_fails_the_calibration(tmp_path, monkeypatch):
    """A reviewer that approves everything is indistinguishable from no
    reviewer, and must not be trusted with real work."""
    rec = _calibrate_in(tmp_path, monkeypatch, lambda m: PASS)
    assert not rec["passed"]
    missed = [r["case"] for r in rec["cases"] if not r["ok"]]
    assert missed and all(r["expect"] == "refuse" for r in rec["cases"] if r["case"] in missed)
    assert pr_bot.calibration_record("kit.test-model") is None


def test_a_reviewer_that_refuses_everything_fails_the_calibration(tmp_path, monkeypatch):
    """Refusing everything looks like diligence while stopping all work."""
    rec = _calibrate_in(tmp_path, monkeypatch,
                        lambda m: _refuse(file="delfin/agent/honesty_scan.py", line=3))
    assert not rec["passed"]
    failed = {r["case"] for r in rec["cases"] if not r["ok"]}
    negatives = {c["id"] for c in pr_bot.corpus_cases() if c["expect"] == "pass"}
    assert negatives and negatives <= failed, "a blanket refusal passed a negative control"


def test_a_refusal_about_the_wrong_file_does_not_catch_a_positive_control(tmp_path, monkeypatch):
    """Otherwise a reviewer that refuses everything scores every positive."""
    def answer(messages):
        case = _case_of(messages)
        if case is not None and case["expect"] == "refuse":
            return _refuse(file="some/other/file.py", line=9, claim="unrelated")
        return PASS
    rec = _calibrate_in(tmp_path, monkeypatch, answer)
    assert not rec["passed"]
    bad = [r for r in rec["cases"] if r["expect"] == "refuse"]
    assert bad and all(not r["ok"] and "instead of" in r["why"] for r in bad)


def _answer_every_control(messages):
    """A reviewer that answers every control the way its label says."""
    case = _case_of(messages)
    if case is not None and case["expect"] == "refuse":
        return _refuse(file=case["expect_files"][0], line=5,
                       claim="nothing reaches it; " + " ".join(case["expect_claim"]))
    return PASS


def test_the_calibration_records_which_model_passed_and_when(tmp_path, monkeypatch):
    rec = _calibrate_in(tmp_path, monkeypatch, _answer_every_control)
    assert rec["passed"], [r for r in rec["cases"] if not r["ok"]]
    assert rec["corpus_digest"] == pr_bot.corpus_digest()
    assert rec["ts"] > 0 and rec["iso"] and rec["calls"] == len(pr_bot.corpus_cases())
    assert pr_bot.calibration_record("kit.test-model")["passed"]


def test_the_corpus_has_both_controls_and_states_where_its_labels_come_from(swarm):
    """A corpus of one direction only cannot tell a reviewer from a rubber
    stamp, and a label whose provenance is unwritten cannot be argued with."""
    man = pr_bot.corpus_manifest()
    cases = man["cases"]
    assert sum(c["expect"] == "refuse" for c in cases) >= 2
    assert sum(c["expect"] == "pass" for c in cases) >= 2
    assert man["labels"] and man["log"], "the corpus does not say where its labels come from"
    for c in cases:
        assert c["source"], c["id"]
        assert pr_bot._corpus_read(c["diff"]).strip()
        assert pr_bot._corpus_read(c["brief"]).strip()
        if c["expect"] == "refuse":
            assert c["expect_files"], f"{c['id']}: a positive control names no file"


def test_every_case_names_its_labeller_the_reason_and_the_base_sha(swarm):
    """A label is one person's reading of one diff at one commit.

    Without the base sha the diff cannot be applied again to check the
    label, and without the reason the label cannot be disputed -- which is
    how the first version of this corpus came to carry a negative control
    that a reviewer was right to refuse.
    """
    for c in pr_bot.corpus_manifest()["cases"]:
        assert c.get("labeller"), f"{c['id']}: no labeller"
        assert len(c.get("label_reason") or "") > 40, f"{c['id']}: no reason for the label"
        sha = c.get("base_sha") or ""
        assert len(sha) == 40 and all(x in "0123456789abcdef" for x in sha), \
            f"{c['id']}: base_sha is not a full sha ({sha!r})"


def test_a_case_the_model_itself_found_is_run_but_not_scored(tmp_path, monkeypatch):
    """An instrument must not measure itself.

    kit.glm-5.3-flash reported the #822 defect during a calibration, so that
    case cannot decide whether kit.glm-5.3-flash is a usable reviewer. It is
    still run and still recorded -- the result is reported separately -- and
    it scores every OTHER model normally.
    """
    marked = [c for c in pr_bot.corpus_cases() if c.get("found_by")]
    assert marked, "no case carries found_by; this check would observe nothing"
    model = marked[0]["found_by"]
    ids = {c["id"] for c in marked}

    def answer(messages):
        # Pass everything: every positive is missed, the marked one included.
        return PASS

    monkeypatch.setattr(pr_bot, "REVIEW_CALIBRATION", tmp_path / "cal.json")
    _fake_reviewer(monkeypatch, answer)
    rec = pr_bot.calibrate([model], out=lambda *a: None)[model]
    rows = {r["case"]: r for r in rec["cases"]}
    assert ids <= set(rows), "the case was not run for the model that found it"
    for cid in ids:
        assert rows[cid]["counted"] is False, f"{cid} was scored for {model}"
        assert not rows[cid]["ok"], "the fixture must make the marked case fail"
    assert ids & set(rec["corpus_cases"]) and not (ids & set(rec["scored_cases"]))
    # Another model is scored on it.
    rec2 = pr_bot.calibrate(["kit.test-model"], out=lambda *a: None)["kit.test-model"]
    rows2 = {r["case"]: r for r in rec2["cases"]}
    for cid in ids:
        assert rows2[cid]["counted"] is True, f"{cid} was not scored for another model"


def test_a_pass_needs_both_directions_left_to_score(tmp_path, monkeypatch):
    """found_by must not become a way to pass by subtraction.

    If every positive control carried the reviewer's own name, a blanket
    'pass' would answer every remaining case and the record would read
    PASSED. A direction with no scored case left is not a measurement.
    """
    monkeypatch.setattr(pr_bot, "REVIEW_CALIBRATION", tmp_path / "cal.json")
    cases = [dict(c, found_by="kit.test-model") if c["expect"] == "refuse" else c
             for c in pr_bot.corpus_cases()]
    monkeypatch.setattr(pr_bot, "corpus_cases", lambda: cases)
    _fake_reviewer(monkeypatch, lambda m: PASS)
    rec = pr_bot.calibrate(["kit.test-model"], out=lambda *a: None)["kit.test-model"]
    assert all(r["ok"] for r in rec["cases"] if r["counted"])
    assert not rec["passed"], "a corpus with no scored positive left reported a pass"


def test_the_record_carries_the_bounds_the_corpus_supports(tmp_path, monkeypatch):
    """The pass is not the result; the bounds are.

    4/4 negatives bound the true false-refuse rate at 52.7 % and 2/2
    positives bound the catch rate at 22.4 %, so a six-case pass means
    "neither a rubber stamp nor a blanket refuser" and nothing about
    catching defects. The record says which bound it supports.
    """
    def answer(messages):
        body = messages[-1]["content"]
        if "honesty_scan.py b/delfin" in body or "_append_dishonest_test_caveat" in body:
            f = ("delfin/agent/honesty_scan.py" if "honesty_scan.py b/delfin" in body
                 else "delfin/agent/engine.py")
            return _refuse(file=f, line=5, claim="nothing reaches it; _append_dishonest_test_caveat")
        return PASS
    rec = _calibrate_in(tmp_path, monkeypatch, answer)
    b = rec["bounds"]
    n_pos = sum(c["expect"] == "refuse" for c in pr_bot.corpus_cases())
    n_neg = sum(c["expect"] == "pass" for c in pr_bot.corpus_cases())
    assert b["positives"]["n"] == n_pos and b["negatives"]["n"] == n_neg
    # This fixture catches exactly the two honesty_scan positives.
    assert b["positives"]["ok"] == 2 and b["negatives"]["ok"] == n_neg
    assert abs(b["false_refuse_rate_at_most"] - (1 - 0.05 ** (1 / n_neg))) < 1e-6
    # 2 of n positives is weaker evidence than 2 of 2 would have been: the
    # misses move the catch bound DOWN, which is the whole point of adding
    # positive controls rather than only negative ones.
    assert b["catch_rate_at_least"] < 0.05 ** 0.5


def test_the_bounds_are_the_exact_binomial_ones_and_not_a_guess(swarm):
    """Known values, from the two the corpus itself was criticised with, plus
    the direction a miss moves them. A normal approximation gives 100 % for
    4/4 and would read as certainty."""
    assert abs(pr_bot.one_sided_bounds(2, 2)["lower"] - 0.2236068) < 1e-6
    assert abs(pr_bot.one_sided_bounds(4, 4)["lower"] - 0.4728708) < 1e-6
    assert abs(1 - pr_bot.one_sided_bounds(9, 9)["lower"] - 0.2831288) < 1e-6
    # A single miss must widen, never narrow, and the upper end must be a
    # real quantile rather than the degenerate 0 a wrong inversion gives.
    assert pr_bot.one_sided_bounds(8, 9)["lower"] < pr_bot.one_sided_bounds(9, 9)["lower"]
    assert abs(pr_bot.one_sided_bounds(8, 9)["upper"] - 0.95 ** (1 / 9)) < 1e-6
    assert pr_bot.one_sided_bounds(0, 4)["lower"] == 0.0
    assert abs(pr_bot.one_sided_bounds(0, 4)["upper"] - (1 - 0.05 ** 0.25)) < 1e-6
    assert pr_bot.one_sided_bounds(0, 0)["lower"] is None


def test_a_corrupt_calibration_file_is_not_silently_overwritten(tmp_path, monkeypatch):
    """A lenient read in a writer destroys data: a cache that will not parse
    must stop the write, not erase every other model's result."""
    monkeypatch.setattr(pr_bot, "REVIEW_CALIBRATION", tmp_path / "cal.json")
    pr_bot.REVIEW_CALIBRATION.write_text('{"kit.a": {"passed": true}, truncated')
    _fake_reviewer(monkeypatch, lambda m: PASS)
    with pytest.raises(Exception):
        pr_bot.calibrate(["kit.test-model"], out=lambda *a: None)
    assert pr_bot.REVIEW_CALIBRATION.read_text().endswith("truncated")


# --- 6. never score the count of reviews ----------------------------------


def test_the_record_is_diagnosis_with_counts_beside_the_rate(swarm):
    """A rate whose denominator the run controls is a target in disguise, so
    the absolute counts stand beside it; and the summary reports what the
    verdicts WERE, which is what a person can act on."""
    s = pr_bot.review_summary([
        {"run": "r", "id": "a", "reviewer": "kit.glm-5.3-flash", "verdict": "refuse"},
        {"run": "r", "id": "b", "reviewer": "kit.glm-5.3-flash", "verdict": "pass"},
        {"run": "r", "id": "c", "reviewer": "", "verdict": "not-calibrated"},
    ])
    assert s["reviewed"] == 2 and s["refused"] == 1 and s["passed"] == 1
    assert s["refusal_rate"] == 0.5 and s["refusal_rate_of"] == 2
    assert s["no_review_possible"] == 1
    assert [v["verdict"] for v in s["verdicts"]] == ["refuse", "pass", "not-calibrated"]


def test_a_run_with_no_review_reports_no_rate(swarm):
    """No denominator, no rate -- a 0/0 printed as 0.00 would read as "nothing
    was refused" when nothing was reviewed at all."""
    s = pr_bot.review_summary([{"run": "r", "id": "a", "reviewer": "", "verdict": "no-reviewer"}])
    assert s["refusal_rate"] is None and s["no_review_possible"] == 1


# --- the diff the reviewer is given --------------------------------------


def test_a_diff_too_large_to_review_refuses(swarm, replies, monkeypatch):
    """Half a diff reviewed and reported on as all of it is the instrument
    this gate exists against."""
    monkeypatch.setattr(pr_bot, "MAX_REVIEW_DIFF", 10)
    ok, why, _ = pr_bot.review_gate(swarm["wt"], "night9", "n9a", swarm["agent"])
    assert not ok and "too large" in why
    assert replies["sent"] == []


def test_a_branch_with_no_brief_refuses(swarm, replies):
    swarm["agent"]["task"].pop("brief")
    ok, why, _ = pr_bot.review_gate(swarm["wt"], "night9", "n9a", swarm["agent"])
    assert not ok and "no task brief" in why
    assert replies["sent"] == []


def test_a_calibration_measured_with_another_question_does_not_count(swarm, monkeypatch):
    """The recorded pass refers to one corpus AND one way of asking. Softening
    the reviewer's instruction, or changing the token budget that decides
    whether a reasoning model answers at all, is a new instrument."""
    assert pr_bot.calibration_record("kit.glm-5.3-flash")
    monkeypatch.setattr(pr_bot, "REVIEW_SYSTEM", pr_bot.REVIEW_SYSTEM + "\nBe lenient.\n")
    assert pr_bot.calibration_record("kit.glm-5.3-flash") is None
    monkeypatch.undo()
    monkeypatch.setattr(pr_bot, "REVIEW_MAX_TOKENS", pr_bot.REVIEW_MAX_TOKENS // 4)
    assert pr_bot.calibration_record("kit.glm-5.3-flash") is None


def test_an_endpoint_that_never_answered_gets_no_rate(tmp_path, monkeypatch):
    """A model id the endpoint does not serve is not a bad reviewer.

    Measured 2026-11-21: kit.glm-5.3 returned HTTP 400 "Model not found" on
    all 24 of its calls and was reported at "catch rate at least 0.0 %,
    false-refuse rate at most 100.0 %" -- two figures computed from twelve
    rows that hold no judgement at all. A rate whose numerator counts a
    transport failure as a wrong reading of code says something the run
    never measured.
    """
    def dead(model, messages, timeout=None):
        raise RuntimeError("HTTP Error 400: Bad Request")

    monkeypatch.setattr(pr_bot, "REVIEW_CALIBRATION", tmp_path / "cal.json")
    monkeypatch.setattr(pr_bot, "kit_chat", dead)
    rec = pr_bot.calibrate(["kit.test-model"], out=lambda *a: None)["kit.test-model"]
    b = rec["bounds"]
    assert b["measurable"] is False
    assert b["catch_rate_at_least"] is None
    assert b["false_refuse_rate_at_most"] is None
    assert b["not_measured"]["n"] == b["not_measured"]["of"] > 0
    # And it still does not pass: no answer is not a pass.
    assert not rec["passed"]


def test_an_unanswered_case_leaves_the_rate_and_still_blocks_the_pass(
        tmp_path, monkeypatch):
    """The two halves of the repair, and they pull in opposite directions.

    An empty reply at finish_reason "length" is the gate's own token limit,
    so it is out of the judgement rate. But dropping it from the PASS rule
    would let a reviewer qualify by not answering the cases it cannot do --
    the laundering this separation must not enable. So: out of the rate, in
    the pass rule, and always reported.
    """
    hard = "pos-833-unmounted"

    def answer(messages):
        case = _case_of(messages)
        if case and case["id"] == hard:
            raise RuntimeError(
                "empty reply (finish_reason='length', completion_tokens=16000,"
                " reasoning_tokens=16000)")
        if case and case["expect"] == "refuse":
            # A correct catch: the file the case is about, and its keywords.
            return _refuse(file=case["expect_files"][0], line=1,
                           claim=" ".join(case.get("expect_claim") or [])
                           or "nothing in the product reaches it")
        return PASS

    rec = _calibrate_in(tmp_path, monkeypatch, answer)
    b = rec["bounds"]
    nm = b["not_measured"]
    assert nm["n"] == 1 and nm["cases"] == [hard] and nm["positives"] == 1
    # Out of the rate: the positives denominator counts only answered rows,
    # so it is one short of the scored positives.
    scored_pos = [r for r in rec["cases"]
                  if r["counted"] and r["expect"] == "refuse"]
    assert b["positives"]["n"] == len(scored_pos) - 1
    assert b["measurable"] is True
    # In the pass rule: one unanswered case is enough to refuse the model,
    # although every case it DID answer was answered correctly.
    assert not rec["passed"]
    assert [r["case"] for r in rec["cases"] if not r["ok"]] == [hard]


def test_the_roster_is_the_served_models_and_not_the_ranking(monkeypatch):
    """The product's kit list is a ranked fallback, not a roster.

    It deliberately keeps an id the endpoint has stopped serving as the next
    candidate in case it returns, and says so in its own comment. Measured
    2026-11-21: calibrating the ranking spent 24 calls on kit.glm-5.3,
    "400 Model not found" on every one.
    """
    import delfin.agent.model_capabilities as mc

    ranked = ("kit.glm-5.3-flash", "kit.glm-5.3", "kit.deepseek-v4-flash")
    monkeypatch.setattr(mc, "MODEL_PREFERENCE_BY_PROVIDER", {"kit": ranked})

    monkeypatch.setattr(mc, "kit_models_served",
                        lambda **kw: {"kit.glm-5.3-flash",
                                      "kit.deepseek-v4-flash"})
    assert pr_bot.kit_reviewer_models() == ("kit.glm-5.3-flash",
                                            "kit.deepseek-v4-flash")

    # An unreadable listing is "could not ask", never "nothing is served":
    # it falls through to the ranking, and the calibration then reports any
    # unserved id as NOT MEASURABLE rather than as a bad reviewer.
    monkeypatch.setattr(mc, "kit_models_served", lambda **kw: set())
    assert pr_bot.kit_reviewer_models() == ranked

    # A listing that names none of them is treated the same way, because a
    # roster of zero reviewers would stop every branch on a probe failure.
    monkeypatch.setattr(mc, "kit_models_served", lambda **kw: {"something.else"})
    assert pr_bot.kit_reviewer_models() == ranked


# --- 9. the reviewer pair: it may abstain, it may not be wrong -------------


def _pair_rows(**counts):
    """Rows in the shape calibrate_pair writes, for pair_bounds to score."""
    rows = []
    for pair, n in counts.items():
        for i in range(n):
            expect = {"pass": "pass", "refuse": "refuse",
                      "wrong_pass": "refuse", "wrong_refuse": "pass",
                      "split": "refuse", "no_answer": "pass"}[pair]
            got = pair.replace("wrong_", "")
            rows.append({"case": f"{pair}-{i}", "expect": expect,
                         "a": got, "b": got, "pair": got,
                         "ok": None if got in ("split", "no_answer")
                         else got == expect,
                         "counted": True, "findings": []})
    return rows


def test_the_pair_abstains_when_the_two_reviewers_disagree(swarm, replies,
                                                           monkeypatch):
    """One passes, one refuses: the branch is held for a person, and the gate
    says UNDECIDED rather than REFUSED. Calling a split a refusal sends an
    author to look for a defect one of two reviewers does not believe in."""
    seen = {}

    def one_each(model, messages, timeout=None):
        seen[model] = True
        return PASS if len(seen) == 1 else _refuse()

    monkeypatch.setattr(pr_bot, "kit_chat", one_each)
    ok, why, rec = pr_bot.review_gate(swarm["wt"], "night9", "n9a", swarm["agent"])
    assert not ok
    assert rec["verdict"] == "split"
    assert "UNDECIDED" in why and "a person decides" in why
    assert "refused by BOTH" not in why


def test_both_passing_is_what_opens_the_gate(swarm, replies):
    """A pass needs both reviewers. One pass is not a pass."""
    replies["set"](PASS)
    ok, why, rec = pr_bot.review_gate(swarm["wt"], "night9", "n9a", swarm["agent"])
    assert ok and rec["verdict"] == "pass"
    assert "both pass" in why
    # The pass line carries the pair's measured error count, so a person
    # reading the log sees what the pass is worth without opening a file.
    assert "wrong of" in why


def test_both_refusing_is_a_refusal_and_keeps_both_sets_of_findings(swarm, replies):
    """A refusal both reviewers reached is the only refusal the gate reports,
    and the record keeps each reviewer's own findings separately."""
    replies["set"](_refuse(file="delfin/agent/x.py", line=9))
    ok, why, _ = pr_bot.review_gate(swarm["wt"], "night9", "n9a", swarm["agent"])
    assert not ok and "refused by BOTH" in why
    rec = json.loads((swarm["wt"] / "scratch" / "REVIEW.json").read_text())
    assert rec["a"]["verdict"] == "refuse" and rec["b"]["verdict"] == "refuse"
    assert len(rec["findings"]) == 2, "both reviewers' findings are kept"


def test_a_pair_that_decided_nothing_fails_its_calibration():
    """The Goodhart case: a pair that abstains on everything has a perfect
    wrong-rate. It must fail, or abstaining becomes the way to pass."""
    b = pr_bot.pair_bounds(_pair_rows(split=6, no_answer=6))
    assert b["decided"] == 0 and b["wrong"] == 0
    assert b["correct_when_agreed_at_least"] is None, (
        "no decisions means no rate, and None is how that is said")


def test_a_pair_that_agrees_on_a_false_alarm_fails():
    """Two models that both false-refuse AGREE on the false alarm, so pairing
    two refusers is not a check. Measured: kit.deepseek-v4-flash +
    kit.glm-5.3-flash both refused neg-824, a sound branch."""
    b = pr_bot.pair_bounds(_pair_rows(**{"refuse": 4, "pass": 2,
                                         "wrong_refuse": 1}))
    assert b["wrong"] == 1 and b["wrong_cases"] == ["wrong_refuse-0"]


def test_a_pair_that_only_ever_agrees_on_passes_is_a_rubber_stamp():
    """Both directions must appear among the decisions, or a pair that waves
    everything through scores a clean sheet."""
    b = pr_bot.pair_bounds(_pair_rows(**{"pass": 6, "split": 3}))
    assert b["wrong"] == 0 and b["agreed_refuse"] == 0, (
        "the counts a pass rule must look at")


def test_an_abstention_is_not_scored_as_a_wrong_answer():
    """A split is the pair declining to decide. Counting it as wrong would
    make the measured pair look worse than a single model and push the
    design back to one reviewer."""
    rows = _pair_rows(**{"pass": 4, "refuse": 1, "split": 3, "no_answer": 1})
    b = pr_bot.pair_bounds(rows)
    assert b["n"] == 9 and b["decided"] == 5 and b["wrong"] == 0
    assert b["abstained"] == 4
    assert b["abstained_split"] == 3 and b["abstained_no_answer"] == 1


def test_the_calibration_file_has_exactly_one_writer():
    """Two copies of one atomic write drift. The model and the pair
    calibration must go through the same function."""
    src = Path(pr_bot.__file__).read_text()
    assert src.count("tmp.replace(REVIEW_CALIBRATION)") == 1
    assert "def _write_calibration" in src
    assert src.count("_write_calibration(cache)") >= 2


def test_the_remedy_the_gate_names_is_a_flag_that_exists():
    """The refusal tells an operator to run --calibrate-pair. A message that
    names a flag the CLI does not have is the defect this gate exists to
    find in other people's code."""
    src = Path(pr_bot.__file__).read_text()
    assert '"--calibrate-pair" in argv' in src
    assert "--calibrate-pair" in src.split('"""')[1], (
        "and it is in the usage line at the top of the file")


def test_the_pass_rule_refuses_each_degenerate_pair():
    """Four clauses, four ways to be unusable, each measured on its own."""
    real = pr_bot.pair_bounds(_pair_rows(**{"pass": 4, "refuse": 1,
                                            "split": 3, "no_answer": 1}))
    assert pr_bot.pair_passes(real), "the measured glm+mistral shape passes"
    assert not pr_bot.pair_passes(
        pr_bot.pair_bounds(_pair_rows(split=6, no_answer=6))), \
        "a pair that decided nothing has a perfect wrong-rate"
    assert not pr_bot.pair_passes(
        pr_bot.pair_bounds(_pair_rows(**{"pass": 6, "split": 3}))), \
        "a pair that only agrees on passes was never shown to catch anything"
    assert not pr_bot.pair_passes(
        pr_bot.pair_bounds(_pair_rows(**{"pass": 2, "refuse": 4,
                                         "wrong_refuse": 1}))), \
        "one shared false alarm is one wrong decision too many"
    assert not pr_bot.pair_passes(
        pr_bot.pair_bounds(_pair_rows(**{"refuse": 6, "split": 2}))), \
        "a blanket refuser pair passes nothing, so it decided one direction"


def test_the_pair_asks_two_different_models(swarm, replies):
    """The direct check that a pair is a pair: two distinct model ids reach
    the endpoint. Everything else about the pair -- the four verdicts, the
    bounds, the pass rule -- is computed from two answers, so one model
    asked twice would satisfy all of it while measuring one opinion."""
    replies["set"](PASS)
    ok, _, rec = pr_bot.review_gate(swarm["wt"], "night9", "n9a", swarm["agent"])
    assert ok
    asked = [s["model"] for s in replies["sent"]]
    assert len(asked) == 2 and len(set(asked)) == 2, asked
    assert set(asked) == set(rec["reviewers"]), (
        "the models recorded are the models asked")


def test_the_state_delfin_keeps_in_a_workspace_is_not_uncommitted_work(swarm, monkeypatch):
    """A member's worktree holds .delfin/ (DELFIN's own workspace state) in a
    repository that does not ignore it; that is not work left uncommitted.
    A real untracked file still is."""
    wt = swarm["wt"]
    (wt / ".delfin").mkdir()
    (wt / ".delfin" / "state.json").write_text("{}")
    monkeypatch.setattr(pr_bot, "PYCI", sys.executable)
    ok, why = pr_bot.check(wt)
    assert "uncommitted" not in why, why
    (wt / "stray.py").write_text("x = 1\n")
    ok, why = pr_bot.check(wt)
    assert not ok and "uncommitted changes: stray.py" in why
