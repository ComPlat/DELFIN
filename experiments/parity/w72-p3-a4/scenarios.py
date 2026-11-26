#!/usr/bin/env python3
"""Seeded-defect measurement of DELFIN's independent review flow.

THE CAPABILITY UNDER TEST
    A review of a diff reports findings a person can go and check -- a file
    and a line inside the diff it was given -- and reports none it could not
    verify.

THE ENTRY POINT IT IS MEASURED THROUGH
    ``delfin-agent swarm wave land [RUN ...]``, which
    ``delfin/agent/cli_swarm.py`` forwards verbatim to
    ``delfin.agent.swarm.land.main(argv)`` -> ``land()`` -> ``review_gate()``
    -> ``review_pair()`` -> ``review_diff()`` -> ``checkable_findings()``.
    This script calls ``land.main([run, "--dry-run"])`` -- never a helper.

WHAT IS REAL AND WHAT IS FAKED
    Real: the throwaway git repository, ``git diff``, the gate, the reviewer
    pair logic, the refusal text, ``scratch/REVIEW.json``, the branch tests
    ``check()`` runs. Faked: the reviewer model, ``land.kit_chat``, which
    answers from a script. Nothing reaches the network; nothing reads
    ``~/.delfin``.

SCENARIOS
    Two seeds. A seed is a branch: a new product module plus a new test.
    ``real``  -- the seeded defect is in the diff (a new module nothing in the
                 product calls: the shape this repository measures most).
    ``clean`` -- the diff carries no defect; a reviewer that claims one is
                 confabulating.
    Six scripted answer kinds cover the ways a confident review answer goes
    wrong about WHERE the defect is: right place, phantom file, line outside
    every hunk, line beside the right one, a false claim at a real place, and
    a phantom claim about a clean branch. 2 seeds x 6 kinds x 2 shapes = 12
    variants.

THE MEASUREMENTS (same definitions on both trees, fixed before the fix)
    REACHED  the review flow ran: a reviewer pair was reached and the fake
             model was asked (>= 2 prompts), and the recorded verdict is one a
             reviewer produced -- not an upstream refusal like "no-brief".
    USED     every finding the flow reported carried a location that could be
             checked against the diff it was given. A reported finding whose
             file is not in the diff, or whose line is outside every hunk,
             means the check did not happen for it.
    SAFE     the branch was not wrongly landed: a seeded defect did not pass,
             and a clean branch did pass.
    CORRECT  SAFE, and nothing unverifiable was reported.

    Ground truth is computed from the SEED (which files the branch wrote and
    how many lines each has), never by parsing the diff with the product's own
    code -- so the oracle cannot agree with the implementation by sharing it.

    Also reported, because they are what the package names: recall (a seeded
    defect refused, with a surviving finding naming where it really is),
    precision (reported findings that name the real place), and how many
    reported findings carry a checkable location.

Run it from any tree of this repository::

    python3 experiments/parity/w72-p3-a4/scenarios.py \
        --out experiments/parity/w72-p3-a4/before.json
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from delfin.agent.swarm import land as L  # noqa: E402


# --------------------------------------------------------------------------
# The seeds: a branch's new files, written by the harness, counted by it
# --------------------------------------------------------------------------
# SEED_A is the classic shape: product code + its own test, and nothing in the
# product ever calls the new function. SEED_B is the same shape in another
# file, so no result rests on one path name.

def _module(name: str, verb: str) -> str:
    return (f'"""{verb} a run, for the terminal."""\n'
            "from __future__ import annotations\n"
            "\n"
            "\n"
            "def run_it(status: str) -> str:\n"
            '    """The word shown beside a finished run."""\n'
            "    return status\n"
            f"    return _dead_code_{name}(status)\n")


def seed_files(seed: str) -> dict[str, str]:
    """The files the branch writes, by path. The harness writes them, so it
    knows each one's line count without asking git or the product."""
    mod = "delfin/agent/quiet_reporter.py" if seed == "A" else "delfin/agent/quiet_counter.py"
    tst = ("tests/test_the_quiet_reporter_reports.py" if seed == "A"
           else "tests/test_the_quiet_counter_counts.py")
    # check() RUNS the branch's tests, so the test file must pass in a
    # throwaway tree: self-contained, importing nothing.
    return {mod: _module(seed, "Summarize" if seed == "A" else "Count"),
            tst: "def test_the_branch_changes_exist():\n    assert True\n"}


#: Where the seeded defect really is: the new def in the new module, which no
#: file in the product calls. Line 5 of the content ``_module`` writes.
SEED_TRUTH = {"A": ("delfin/agent/quiet_reporter.py", 5),
              "B": ("delfin/agent/quiet_counter.py", 5)}

#: scratch/PR.md. Deliberately free of home paths and of "WIP": ``land()``
#: refuses before the review otherwise, and this measures the review.
PR_MD = "# Add the quiet reporter\n\nAdds one helper and its test.\n"
BRIEF = ("PACKAGE demo: add a small helper and its test. A capability the "
         "product never reaches is a defect, not a feature.\n")


# --------------------------------------------------------------------------
# The scripted answers: how a confident review answer goes wrong about where
# --------------------------------------------------------------------------

PASS = json.dumps({"verdict": "pass", "findings": []})


def _refusal(file: str, line, claim: str, why: str) -> str:
    return json.dumps({"verdict": "refuse",
                       "findings": [{"file": file, "line": line,
                                     "claim": claim, "why_it_matters": why}]})


def kinds(seed: str) -> dict[str, dict]:
    """The six answer kinds, each with how it answers a real defect and a
    clean branch.

    ``truth`` says whether the answer, if the flow reports it, points at where
    the seeded defect really is (only meaningful on a real defect).
    """
    mod, line = SEED_TRUTH[seed]
    n_lines = len(seed_files(seed)[mod].splitlines())
    return {
        # Right place: the file and the line the defect is on.
        "precise": {
            "real": _refusal(mod, line, "nothing in delfin/ calls run_it; "
                             "the new module is never mounted",
                             "the capability never runs"),
            "clean": PASS, "points_at_truth": True},
        # Right about the defect, wrong about the place: a file this branch
        # never touched.
        "phantom-file": {
            "real": _refusal("delfin/agent/never_heard_of_it.py", 42,
                             "nothing calls the new function", "it never runs"),
            "clean": PASS, "points_at_truth": False},
        # Right file, wrong place: a line no hunk of this diff touches.
        "wrong-line": {
            "real": _refusal(mod, n_lines + 3000,
                             "nothing calls the new function", "it never runs"),
            "clean": PASS, "points_at_truth": False},
        # Off by one: the line beside the one that matters.
        "off-by-one": {
            "real": _refusal(mod, line + 1, "nothing calls the new function",
                             "it never runs"),
            "clean": PASS, "points_at_truth": False},
        # A false claim at a real place: an unremarkable line of the diff,
        # reported as a defect. No location check can catch this one, and it is
        # not meant to -- it is measured so the numbers say so.
        "false-at-real-place": {
            "real": PASS,
            "clean": _refusal(mod, 2, "this import is never used",
                              "the branch carries dead weight"),
            "points_at_truth": False},
        # A phantom claim about a clean branch: no such file, no defect.
        "phantom-on-clean": {
            "real": PASS,
            "clean": _refusal("delfin/agent/not_in_this_branch.py", 9,
                              "the guard at line 9 is removed",
                              "a check is lost"),
            "points_at_truth": False},
    }


def scenarios() -> list[dict]:
    """The 12 variants. ``expect`` is what the branch deserves."""
    out: list[dict] = []
    for seed in ("A", "B"):
        files = seed_files(seed)
        lines = {p: len(t.splitlines()) for p, t in files.items()}
        kk = kinds(seed)
        for kind, spec in kk.items():
            for shape in ("real", "clean"):
                out.append({
                    "id": f"{seed}-{kind}-{shape}",
                    "seed": seed, "kind": kind, "shape": shape,
                    "expect": "refuse" if shape == "real" else "pass",
                    "answer": spec[shape],
                    # Ground truth from the seed, not from a diff parse.
                    "diff_files": sorted(files),
                    "file_lines": lines,
                    "truth_file": (SEED_TRUTH[seed][0]
                                   if shape == "real" and spec["points_at_truth"]
                                   else None),
                    "truth_line": (SEED_TRUTH[seed][1]
                                   if shape == "real" and spec["points_at_truth"]
                                   else None),
                })
    return out


def uncheckable(row: dict) -> list[str]:
    """Reported findings whose location could NOT be checked against the seed.

    The oracle is the seed's own file list and line counts, so this judges the
    product's output without reusing the product's parsing.
    """
    out = []
    for f in row["findings"] or []:
        name, line = f.get("file"), f.get("line")
        if name not in row["diff_files"]:
            out.append(f"{name}:{line} (file not in the diff)")
        elif not isinstance(line, int) or not 1 <= line <= row["file_lines"][name]:
            out.append(f"{name}:{line} (line outside the new file)")
    return out


# --------------------------------------------------------------------------
# One scenario: a throwaway repository, the real CLI verb, a faked model
# --------------------------------------------------------------------------

_ENV_KEYS = ("GIT_AUTHOR_NAME", "GIT_AUTHOR_EMAIL", "GIT_COMMITTER_NAME",
             "GIT_COMMITTER_EMAIL", "GIT_CONFIG_GLOBAL", "HOME")


def _git(*args, cwd):
    r = subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True)
    assert r.returncode == 0, (args, r.stderr)
    return r.stdout


def _build_repo(root: Path, seed: str) -> Path:
    wt = root / "wt"
    (wt / "scratch").mkdir(parents=True)
    (wt / "scratch" / ".gitignore").write_text("*\n")
    # A real throwaway gitconfig: pointing GIT_CONFIG_GLOBAL at a file that
    # does not exist makes git fall back to the SYSTEM config, whose identity
    # would then be a person's name and address in a measurement file.
    (root / "gitconfig").write_text(
        "[user]\n\tname = tester\n\temail = tester@example.org\n")
    _git("init", "-q", "-b", "main", cwd=wt)
    (wt / "base.txt").write_text("base\n")
    _git("add", "-A", cwd=wt)
    _git("commit", "-qm", "base", cwd=wt)
    _git("update-ref", "refs/remotes/origin/main", "HEAD", cwd=wt)
    for rel, text in seed_files(seed).items():
        p = wt / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
    (wt / "scratch" / "PR.md").write_text(PR_MD)
    _git("add", "-A", cwd=wt)
    _git("commit", "-qm", "add the helper and its test", cwd=wt)
    return wt


def _calibrate(root: Path, roster) -> None:
    """A passing calibration for every pair, written the way ``calibrate``
    writes one -- an uncalibrated reviewer is refused before any prompt is
    sent, which would measure the calibration check, not the review."""
    rec = {"ts": time.time(), "iso": "2026-01-01T00:00:00+0000",
           "corpus_digest": L.corpus_digest(), "passed": True, "cases": []}
    cache = {m: dict(rec) for m in roster}
    bounds = {"n": 9, "decided": 5, "wrong": 0, "wrong_cases": [],
              "abstained": 4, "abstained_split": 3, "abstained_no_answer": 1,
              "agreed_pass": 4, "agreed_refuse": 1,
              "correct_when_agreed_at_least": 0.549}
    for i, a in enumerate(roster):
        for b in roster[i + 1:]:
            if L.model_family(a) == L.model_family(b):
                continue
            cache[L.pair_calibration_key(a, b)] = dict(
                rec, models=sorted((a, b)), bounds=bounds)
    (root / "calibration.json").write_text(json.dumps(cache))


def run_scenario(sc: dict, root: Path) -> dict:
    """Drive ``land.main`` once for one variant, with only the model faked."""
    sent: list[dict] = []
    roster = ("kit.glm-5.3-flash", "kit.deepseek-v4-flash",
              "kit.mistral-small-4-119b-a8b")

    def fake_chat(model, messages, timeout=None):
        sent.append({"model": model})
        return sc["answer"]

    saved: dict = {}

    def _set(name, value):
        saved[name] = getattr(L, name)
        setattr(L, name, value)

    # The environment first: the branch's own commits are made below and must
    # carry the throwaway identity, not whatever the host has configured.
    env_saved = {k: os.environ.get(k) for k in _ENV_KEYS}
    os.environ.update({
        "GIT_AUTHOR_NAME": "tester", "GIT_AUTHOR_EMAIL": "tester@example.org",
        "GIT_COMMITTER_NAME": "tester", "GIT_COMMITTER_EMAIL": "tester@example.org",
        "GIT_CONFIG_GLOBAL": str(root / "gitconfig"), "HOME": str(root),
        "PYTHONDONTWRITEBYTECODE": "1"})
    wt = _build_repo(root, sc["seed"])
    run = root / "night9"
    run.mkdir(exist_ok=True)
    (root / "n9a.log").write_text('{"type": "notice", "text": "hello"}\n')
    (run / "state.json").write_text(json.dumps({"agents": {"n9a": {
        "worktree": str(wt), "log": str(root / "n9a.log"),
        "branch": "swarm/night9/n9a",
        "task": {"id": "n9a", "model": "alias-code", "brief": BRIEF}}}}))

    _calibrate(root, roster)
    _set("ROOT", root)
    _set("REVIEW_CALIBRATION", root / "calibration.json")
    _set("REVIEW_LOG", root / "review_log.jsonl")
    _set("kit_chat", fake_chat)
    _set("kit_reviewer_models", lambda: roster)
    _set("identity_email", lambda: "tester@example.org")
    noise = io.StringIO()
    try:
        with contextlib.redirect_stdout(noise), contextlib.redirect_stderr(noise):
            L.main([run.name, "--dry-run"])
        landing = (wt / "scratch" / "LANDING.txt").read_text().strip()
        rf = wt / "scratch" / "REVIEW.json"
        review = json.loads(rf.read_text()) if rf.is_file() else {}
    finally:
        for name, value in saved.items():
            setattr(L, name, value)
        for k, v in env_saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v

    outcome = landing.split(":", 1)[0]
    verdict = review.get("verdict", "")
    return {
        **{k: sc[k] for k in ("id", "seed", "kind", "shape", "expect",
                              "diff_files", "file_lines", "truth_file",
                              "truth_line")},
        "prompts_sent": len(sent),
        "verdict": verdict,
        "outcome": outcome,
        "landing": landing[:400],
        "findings": review.get("findings") or [],
        "gate_lines": [l for l in noise.getvalue().splitlines()
                       if "review" in l][:4],
        # REACHED: the review flow ran and a reviewer verdict was recorded.
        "reached": bool(len(sent) >= 2 and verdict in
                        ("pass", "refuse", "split", "no_answer")),
    }


def judge(row: dict) -> dict:
    """Add USED / SAFE / CORRECT, and the inputs the rates are derived from."""
    bad = uncheckable(row)
    row["unverified_findings"] = bad
    row["used"] = bool(row["reached"]) and not bad
    if row["expect"] == "refuse":
        safe = row["verdict"] != "pass"
        about = [f for f in row["findings"]
                 if f.get("file") == row["truth_file"]
                 and f.get("line") == row["truth_line"]]
        row["caught"] = row["verdict"] == "refuse" and bool(about)
        row["safe"] = safe
        row["why"] = ("refused, and a finding names the seeded place" if row["caught"]
                      else "held for a person (no checkable finding survived)"
                      if safe and not row["findings"]
                      else "refused, but no finding names the seeded place"
                      if safe else f"verdict {row['verdict']!r} on a seeded defect")
    else:
        safe = row["verdict"] == "pass"
        row["caught"] = False
        row["safe"] = safe
        row["why"] = ("clean diff passed" if safe else
                      f"verdict {row['verdict']!r} on a clean diff")
    row["correct"] = bool(row["safe"]) and not bad
    return row


def summarise(rows: list[dict]) -> dict:
    positives = [r for r in rows if r["shape"] == "real"]
    negatives = [r for r in rows if r["shape"] == "clean"]
    kept = [f for r in rows if r["verdict"] == "refuse" for f in r["findings"]]
    good = [f for r in rows if r["expect"] == "refuse" for f in r["findings"]
            if f.get("file") == r["truth_file"] and f.get("line") == r["truth_line"]]
    unverifiable = [u for r in rows for u in r["unverified_findings"]]
    return {
        "n_variants": len(rows),
        "REACHED": sum(bool(r["reached"]) for r in rows),
        "USED": sum(bool(r["used"]) for r in rows),
        "CORRECT": sum(bool(r["correct"]) for r in rows),
        "SAFE": sum(bool(r["safe"]) for r in rows),
        "recall_refused_with_right_place": (
            f"{sum(bool(r['caught']) for r in positives)}/{len(positives)}"),
        "precision_finding_names_real_place": (
            f"{len(good)}/{len(kept)}" if kept else "0/0"),
        "findings_reported": len(kept),
        "findings_with_checkable_location": len(kept) - len(unverifiable),
        "findings_unverifiable": len(unverifiable),
        "negatives_wrongly_refused": sum(1 for r in negatives
                                         if r["verdict"] != "pass"),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="", help="where to write the JSON")
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args(list(sys.argv[1:] if argv is None else argv))

    rows = []
    for sc in scenarios():
        root = Path(tempfile.mkdtemp(prefix="parity-"))
        try:
            row = run_scenario(sc, root)
        finally:
            if not args.keep:
                shutil.rmtree(root, ignore_errors=True)
            else:
                print(f"kept {root}")
        rows.append(judge(row))

    doc = {
        "package": "w72-p3-a4",
        "capability": ("the independent review flow reports findings a person "
                       "can check (file and line inside the diff) and reports "
                       "none it could not verify"),
        "entry_point": ("delfin.agent.swarm.land.main([RUN, --dry-run]), the "
                        "verb delfin-agent swarm wave land forwards to; only "
                        "land.kit_chat is faked"),
        "command": "python3 experiments/parity/w72-p3-a4/scenarios.py",
        "tree": {
            "module_file": str(L.__file__),
            "head": subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                                   capture_output=True, text=True).stdout.strip(),
            "dirty_at_runtime": subprocess.run(
                ["git", "status", "--porcelain", "delfin/agent/swarm/land.py"],
                cwd=REPO, capture_output=True, text=True).stdout.strip(),
            "python": sys.version.split()[0],
        },
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "definitions": {
            "REACHED": "a reviewer pair was reached and asked; the verdict is a "
                       "reviewer verdict",
            "USED": "every reported finding carried a location that could be "
                    "checked against the branch's own files (oracle: the seed)",
            "SAFE": "a seeded defect was not landed; a clean branch was not held",
            "CORRECT": "SAFE and nothing unverifiable was reported",
        },
        "summary": summarise(rows),
        "variants": rows,
    }
    text = json.dumps(doc, indent=1) + "\n"
    if args.out:
        Path(args.out).write_text(text)
    print(json.dumps(doc["summary"], indent=1))
    for r in rows:
        print(f"{r['id']:34s} R={int(r['reached'])} U={int(r['used'])} "
              f"C={int(r['correct'])} S={int(r['safe'])} "
              f"verdict={r['verdict']:10s} unverif={len(r['unverified_findings'])}"
              f"  {r['why']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
