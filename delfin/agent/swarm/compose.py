"""The wave as it has been run: packages in, the whole task list out.

Every wave that reached main from 2026-11-05 on had the same shape, and
each stage exists because leaving it out cost something measured:

  arms        each package is worked on by two or three members of
              DIFFERENT model families, independently, on the same brief --
              where they agree the result is robust, where they differ the
              difference is the finding;
  integrator  one per package, once its arms have ended (or a quorum has
              and waited): it reads every arm's branch and report, re-runs
              their tests on the base (red) and the branch (green), takes,
              combines or rejects item by item, and asks for real runs;
  lander      once the integrators are final: bundles the landing branches
              into ONE branch for ONE full CI run (a serial lander starves
              every other merger) and asks for the real runs of the whole;
  reviewers   four, of ANOTHER family than the lander, read-only, each with
              one focus -- behaviour, tests, safety, claims;
  fixer       works through every confirmed finding on the lander's branch
              and writes the final pull-request text.

Then the orchestrator reviews the result itself -- the members' reports
are not taken on trust: their tests are re-run (red without the fix),
their claims re-measured, their code read -- lands it
(``swarm/land.py``), and a person merges.

``compose`` returns the task list for ``wave.start``; it starts nothing.
Every member's brief ends with ``DELIVERY_RULES``, the rules each of which
was written after a measured defect, with its cost.
"""
from __future__ import annotations

import re
from typing import Any

from . import wave

DELIVERY_RULES = """
HOW THIS WAVE WORKS. Each rule comes from a defect that was measured, and
the cost of ignoring it is written next to it.

1. CONTROL FIRST, AND COMMIT IT. Before any edit, run the suites your
   package names on the untouched tree and record the numbers. Branches
   have been refused for failures that were already on main.

2. DRIVE THE PRODUCT ENTRY POINT, NEVER THE HELPER. A change counts when a
   real entry point reaches it: cli.main([...]), the engine's own turn,
   the REPL builtin, the tool dispatch, `delfin-agent chat -p`. A test
   that calls the function directly measures that the function exists.

3. WRITE THE NUMBER DOWN BEFORE YOU MEASURE IT. Put your prediction in
   scratch/PR.md before the measurement. The gap between prediction and
   data is the finding.

4. A NEW TEST MUST BE RED ON THE BASELINE, for the reason it names. Show
   it red on the base and green after, with a control case that passes on
   both. A test that is green on the base measured nothing.

5. AN ASSERTION THAT CANNOT FAIL IS NOT A TEST. Before you believe one,
   name the input that would make it fail. A set membership is not a
   capability: two changes added tool names to an allow-list, passed
   their own tests, and every call answered "Unknown tool".

6. A TEST BINDS TO DELFIN, NEVER TO THIS MACHINE. No skip that depends on
   what the host has installed. A ratchet fails CI for it. Assert the
   documented outcome for each case instead.

7. AN IMPORT IS NOT A CALLER, AND A MENTION IS NOT A CALL. If you make
   something reachable, show the call that reaches it.

8. A CHANGED DEFAULT IS A CHANGED MEANING FOR EVERY CALLER. Before you
   change what an existing function returns by default, grep every caller
   in delfin/ and classify the question each one asks.

9. A SAMPLE THAT DOES NOT LOOK LIKE THE USER ANSWERS ANOTHER QUESTION.
   State the population of any statistic you derive, and never apply it
   to a role or a workspace it was not measured on.

10. A NUMBER WHOSE SOURCE WAS TRUNCATED IS NOT A MEASUREMENT. Re-run with
    a command whose output fits and quote that.

11. A TOOL CALL WRITTEN AS TEXT DID NOT RUN. When DELFIN says a call "was
    written as text instead of calling the tool", issue it properly and
    check its output arrived.

12. NEVER EDIT A FILE WHILE A SUITE RUNS OVER IT. inspect.getsource reads
    the file from disk, so the run then reports your edit, not the code.

13. A FILE DELFIN WILL NOT LET YOU WRITE IS A ROUTE, NOT A WALL. Do not go
    round the self-modification guard with bash. Write the change as
    scratch/<file>.patch (a unified diff against the base), put its test
    in tests/ as usual, and say in PR.md that the patch needs applying.

14. RUN THE SUITES THAT CONSUME IT, AND THE ONES THAT SCAN EVERYTHING. A
    grep for importers misses repository-wide invariants (personal data,
    reachability, machine binding, catalogue drift). Run tests/ completely
    before PR.md, or say plainly that you did not.

15. TOUCH ONLY THE FILES YOUR PACKAGE NAMES, plus their tests. Other
    packages run in parallel on the same files' neighbours.

16. NO PERSONAL OR SITE DATA in anything you commit or write in PR.md: no
    real account names, home paths or host names.

17. scratch/PR.md IS WHERE YOU FINISH. First line: the title. Then one
    table row per item (item, decision, commit, reason in one line), what
    you predicted and what you measured, what you did NOT do and why, and
    every suite you ran with its numbers.

18. ~/.delfin IS NOT YOURS TO READ, by design: the sandbox hides it. To
    measure further, write the real run you need into
    scratch/realrun.json (REAL RUNS): you cannot run DELFIN against a
    model yourself, the cage removes every key from your shell.

WHAT IS NOT SCORED: the number of things you change. One change with a
before/after measurement beats five without.
"""

REAL_RUNS = """
REAL RUNS. Nothing lands without a real run through the product. You
cannot make one yourself (the cage removes every key from your shell).
Ask for it instead, in scratch/realrun.json:

    {"probes": [{"name": "delegates",
                 "prompt": "<what a person would type>",
                 "workspace": "empty",
                 "turns": 30}]}

At most three probes; name [a-z0-9-]; prompt at most 2000 characters;
"workspace": "repo" gives the run a copy of the branch's tree instead of
an empty directory; turns at most 40. Each probe runs on your base AND on
your branch HEAD, same prompt and model, fresh HOME each, so the pair is a
before/after. A probe must show the change: a question whose answer
differs when the change works. Write your prediction for each probe into
scratch/PR.md BEFORE you end your turn. Then end the turn: the runs are
made for you and you are resumed with scratch/realrun/summary.md.
"""

ARM_NOTE = """
YOU ARE ONE OF SEVERAL ARMS on this package, each of another model
family, working independently; an integrator reads every arm's branch and
report afterwards and keeps what is best measured. Do not ask for real
runs yourself: name in scratch/PR.md the probes (prompt, workspace,
expected difference) that would show your change, for the integrator.
"""

#: The rules that decide what reaches the base, one string per rule so a
#: brief can carry a chosen subset. Each was added after a landing that
#: measured its absence: the sentence after the dash is what it cost.
#:
#: They are data, not prose welded into the lander's brief text, for the
#: same reason the delivery rules are: the fixer needs B and E but not A,
#: and a rule an operator has to retype into a task by hand is a rule that
#: is missing from every wave nobody remembered to patch (w42i/w48).
LANDING_RULES = {
    "A": """ONE BRANCH. Build the bundle on this worktree's own branch \
(wave/{land_ref}) and on no second branch: a wave that built the bundle on a \
second branch left every reviewer's checkout with an empty diff, so four \
reviews were written against nothing.""",
    "B": """NOTHING RED LANDS. A package whose test is red on the finished \
bundle is fixed or left out -- never merged with a red test because its own \
arm reported it green: run it on the bundle, which is the only tree that \
matters.""",
    "C": """A PACKAGE THAT NEEDS A PATCH IS A PACKAGE YOU APPLY. If a report \
says it lands only with a scratch/*.patch, apply that patch as its own commit \
with its test; if you will not, leave the package out and say why in the \
report. A patch nobody applied is code that never reached the base.""",
    "D": """NO ENTRY POINT, NO CAPABILITY. A package that no entry point \
reaches -- nothing in delfin/ calls it -- does not land as a capability; an \
import is not a caller. Land it as dead code only if the report says so in \
those words.""",
    "E": """A COMMIT MESSAGE IS PUBLIC TEXT. The person is the author; an \
agent's commit carries exactly ONE `Co-Authored-By:` line, for DELFIN's agent \
account, and no other trailer -- no second co-author line and no \
`Agent: ... via ...` line. No model, provider, assistant product, institution \
or host names anywhere in a commit message, trailer included: strip them from \
what you bundle. No home paths, account names or person names either.""",
    "F": """AN INTEGRATOR BRANCH 0 COMMITS AHEAD OF THE BASE CARRIES NOTHING. \
`git log --oneline {base}..<branch> | wc -l` is 0 for a branch whose report \
describes work: it has no commits, so there is nothing to merge. Find the \
branch its report actually measured, or record that the package is empty and \
say so.""",
}

#: The keys of each stage's share of ``LANDING_RULES``.
LANDER_RULES = tuple(LANDING_RULES)
#: The fixer edits the bundle itself, so it carries the two rules that decide
#: what a line of that bundle may look like -- not the merge rules.
FIXER_RULES = ("B", "E")

#: What a revision name may carry into a brief. The briefs are full of shell
#: commands the members are told to run, and ``base`` is filled into them, so
#: the accepted set is git's own for a branch or commit name: no spaces, no
#: quotes, no ``$``/``;``/``|``/``&``/backticks and no newline.
_BASE_OK = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._/~^-]*$")


#: An arm or an integrator commits too, and its commit messages are read on
#: the branch someone else merges. Only the public-text rule is theirs: they
#: merge nothing, and the other five would dilute their own package text.
COMMIT_RULES = ("E",)


def _rules_text(keys, *, land_ref: str, base: str) -> str:
    """The named rules as brief text, in the order named, filled in."""
    return "\n\n".join(LANDING_RULES[k].format(land_ref=land_ref, base=base)
                       for k in keys)

REVIEW_FOCUS = {
    "behaviour": """Nothing that worked before may work worse. Read the WHOLE diff
`git diff {base}...HEAD`. For every changed function: who calls it (grep
delfin/), and does each caller still get what it got before when it asks
the same question? A changed default, a renamed key, a widened or
narrowed condition, a removed fallback: each is a finding with file:line
and the caller that is affected.""",
    "tests": """Every new or changed test must be able to fail. For each: name the
input that makes it fail; run it on {base} (it must be red for the reason
it names, unless it is a control) and on the branch (green). Find:
assertions that cannot fail, tests that call a helper instead of the
product entry, skips that depend on the machine, fixed dates, tests that
write outside a temporary HOME. Run the full suite serially once and
compare with the lander's numbers.""",
    "safety": """Security and public text. Search the diff and the lander's
scratch/PR.md for account names, home paths, host names, person names,
keys or key-shaped strings, names of other agent projects or of any
assistant product. Check every new path that runs a command, writes a
file, reaches the network or a device: is it behind the same gate as
before? A gate that a new spelling of the same act walks around is a
finding. A new third-party dependency without its licence recorded is a
finding.""",
    "claims": """Every claim must have its evidence. For each row of the lander's
PR.md and of each integrator report (scratch/arms/*/PR.md): is the commit
there, does the measurement it quotes come from a run that was made, does
the real-runs table (scratch/arms/*/realrun/summary.md) show what the text
says it shows? A prediction quoted as a result, a number without its run,
a "works" without the call that showed it: each is a finding.""",
}


def _family(model: str) -> str:
    try:
        from .land import model_family
        return model_family(model)
    except Exception:
        return ""


def _branch_of(ref: str) -> str:
    """The branch member ``run/tid`` worked on, as its run recorded it."""
    run, _, tid = ref.partition("/")
    try:
        return wave.load_state(run)["agents"][tid]["branch"]
    except (KeyError, wave.WaveError):
        return f"wave/{ref}"


def _validate_base(base: Any) -> str:
    """The revision every brief is told to run against -- and nothing else.

    ``base`` is interpolated into shell commands the members are INSTRUCTED to
    run (``git log --oneline {base}..HEAD | wc -l`` in the reviewer brief, the
    integrator's control step, rule F). It arrives from a spec file, so it is
    an input like a task id -- and it was the one such input that reached a
    mandated command unvalidated: ``origin/main && echo INJECTED`` used to be
    composed into every reviewer's first instruction.

    A revision name is what git accepts as a branch or commit: letters, digits
    and ``._/-``, plus ``~`` and ``^`` for ancestry. Spaces, quotes, ``$``,
    ``;``, ``|``, ``&``, backticks and newlines are not in that set, and a
    refusal names the value rather than silently rewriting it.

    Rejected: escaping or quoting the value into the command. The briefs are
    read by a model, not only by a shell, and a backslash in an instruction is
    a second thing to get wrong; refusing an input that is not a revision keeps
    the instruction identical to the string a person would type.
    """
    text = str(base if base is not None else "").strip()
    if not text:
        raise wave.WaveError("base must name a revision (empty is not a branch)")
    if not _BASE_OK.match(text):
        raise wave.WaveError(
            f"base {text!r} is not a revision name: it may carry letters, "
            "digits and . _ / - ~ ^ only, because the members are told to run "
            "`git log --oneline <base>..HEAD` and a base is not a shell command")
    return text


def _validate_packages(packages: Any) -> list[dict]:
    if not isinstance(packages, list) or not packages:
        raise wave.WaveError("packages must be a non-empty list of {id, title, body}")
    seen = set()
    for p in packages:
        if not isinstance(p, dict) or not str(p.get("body") or "").strip():
            raise wave.WaveError("every package needs an id, a title and a body")
        pid = str(p.get("id") or "")
        if not wave._ID.match(pid) or pid in seen or pid.startswith("zz-"):
            raise wave.WaveError(f"package id {pid!r} is not a unique task id (zz- is the stages')")
        seen.add(pid)
    return packages


def compose(run: str, packages: list[dict], *, arms: list[dict], integrator: dict,
            lander: dict, reviewer: dict, fixer: dict, base: str = "origin/main",
            stages: tuple = ("arms", "integrate", "land", "review", "fix"),
            quorum: int = 2, quorum_wait_s: int = 10800, rules: str = DELIVERY_RULES,
            real_runs: str = REAL_RUNS, arm_refs: dict | None = None) -> list[dict]:
    """The task list of one wave in run ``run``.

    ``packages``: ``[{"id", "title", "body"}]`` -- ``body`` is the package's
    own brief (what to do, which files, which suites, what to measure).
    ``arms``, ``integrator``, ``lander``, ``reviewer``, ``fixer``: each
    ``{"provider", "model"}`` (``arms`` a list, one per arm). The reviewers
    must be of another model family than the lander when both families are
    known. ``stages`` leaves out what is not wanted (no arms: the integrator
    works the package itself). ``arm_refs`` names arms that already ran, per
    package (``{"p1": ["oldrun/p1-x", ...]}``): the integrator reads and
    collects those instead of new ones -- a wave whose integration stage was
    stopped is continued without running its arms again."""
    packages = _validate_packages(packages)
    base = _validate_base(base)
    arm_refs = dict(arm_refs or {})
    if "arms" in stages and len(arms) < 1:
        raise wave.WaveError("arms needs at least one {provider, model}")
    lf, rf = _family(lander.get("model", "")), _family(reviewer.get("model", ""))
    if "review" in stages and lf and rf and lf == rf:
        raise wave.WaveError(f"the reviewers are of the lander's family ({lf}); "
                             "a review of a model's work by its own family is not independent")

    def m(spec: dict) -> dict:
        return {k: spec[k] for k in ("provider", "model") if spec.get(k)}

    tasks: list[dict] = []
    finals: list[str] = []
    land_ref = f"{run}/zz-land"
    for p in packages:
        pid, title, body = p["id"], p.get("title") or p["id"], p["body"]
        arm_ids = list(arm_refs.get(pid) or [])
        if "arms" in stages and not arm_ids:
            for i, spec in enumerate(arms, 1):
                tid = f"{pid}-a{i}"
                arm_ids.append(f"{run}/{tid}")
                tasks.append({"id": tid, "base": base, "seconds": 7200, "turns": 120, **m(spec),
                              "no_realrun": True,
                              "brief": f"PACKAGE {pid}: {title}\n\n{body}\n{ARM_NOTE}"
                                       f"{_rules_text(COMMIT_RULES, land_ref=land_ref, base=base)}\n{rules}"})
        if "integrate" in stages:
            listing = "\n".join(f"  - {a}: branch {_branch_of(a)}, its scratch/ in scratch/arms/{a.split('/')[1]}/"
                                for a in arm_ids) or "  (none: you work the package yourself)"
            tid = f"{pid}-int"
            t = {"id": tid, "base": base, "seconds": 21600, "turns": 200, **m(integrator),
                 "brief": f"""INTEGRATOR for package {pid} ("{title}").

{len(arm_ids)} member(s) worked on this package independently, each on its own
branch. Turn their work into ONE branch that is correct, tested, measured
and shown working in a real run -- or say, with the evidence, that none of
it should land.

The arms (their branches are in this repository's refs):
{listing}

STEPS
1. Control: on your untouched tree ({base}) run the suites the package
   names and record the numbers in scratch/PR.md.
2. For each arm: `git log --oneline {base}..<its branch>`, read the diff and
   its report. Do not trust the report: re-run its new tests on {base}
   (red, for the reason they name) and on its branch (green). An arm may
   have ended unfinished or still run: read its branch as it stands and
   say so.
3. Decide per item: take, combine (cherry-pick, resolve conflicts against
   {base}), or reject -- one table row each with the reason. Where arms
   solved the same thing, keep the better measured one, not the longer.
4. Run the suites that consume what you changed AND the suites that scan
   the whole repository.
5. Ask for your real runs (REAL RUNS below), prediction first.
6. When resumed: judge the runs, fix what they show, finish scratch/PR.md
   (first line: title and whether this should land).

The package:
{body}
{real_runs}
{_rules_text(COMMIT_RULES, land_ref=land_ref, base=base)}
{rules}"""}
            if arm_ids:
                t.update(after=arm_ids, collect=arm_ids, after_min=min(quorum, len(arm_ids)),
                         after_wait_s=quorum_wait_s)
            tasks.append(t)
            finals.append(f"{run}/{tid}")
    if "land" in stages and finals:
        tasks.append({"id": "zz-land", "base": base, "seconds": 21600, "turns": 200, **m(lander),
                      "after": [f"{r}.final" for r in finals], "collect": finals,
                      "after_min": max(1, len(finals) - 1), "after_wait_s": quorum_wait_s,
                      "brief": f"""LANDER: one branch for the base, one CI run.

{len(finals)} integrator(s) each produced one branch (wave/{run}/<pkg>-int) and a
report (scratch/arms/<pkg>-int/PR.md, with its real runs in
scratch/arms/<pkg>-int/realrun/). Bundle what their reports say should land:

1. Read every report. Skip a package whose first line says it should not
   land, and say so.
2. Merge the landing branches onto {base} one by one (git merge --no-ff),
   resolving conflicts three-way per hunk -- a merge hunk holds rows
   nobody touched, so neither union nor intersection is right.
3. Run the FULL test suite once, serially. Every failure: run it on {base}
   too (the control). Fix what is yours, record what was already red.
4. Ask for real runs (REAL RUNS below) that cover the wave's centre.
5. When resumed: finish scratch/PR.md as the pull-request text -- title,
   one row per package (landed or not, commits, measured effect), the
   full suite numbers with their control, the real-runs table. No account
   names, home paths, host names or person names in it.

THE LANDING RULES -- each one is a landing that measured its absence:

{_rules_text(LANDER_RULES, land_ref=land_ref, base=base)}
{real_runs}
{rules}"""})
    revs = []
    if "review" in stages and "land" in stages:
        for key, text in REVIEW_FOCUS.items():
            rid = f"zz-rev-{key}"
            revs.append(f"{run}/{rid}")
            tasks.append({"id": rid, "base": f"wave/{land_ref}", "seconds": 10800, "turns": 150,
                          **m(reviewer), "after": [f"{land_ref}.final"],
                          "collect": [land_ref, *finals],
                          "no_realrun": True, "land": "none",
                          "brief": f"""REVIEWER ({key}) of the bundle another model built.

Your worktree is the lander's branch (wave/{land_ref}); its reports are in
scratch/arms/zz-land/, and every integrator of this run reports in
scratch/arms/<pkg>-int/ -- all of those directories are collected here, so
the claims you check are the integrators' as well as the lander's. You CHANGE
NOTHING in the code and commit nothing: you read, run and report.

BEFORE ANYTHING ELSE: `git log --oneline {base}..HEAD | wc -l`. If it is 0 the
bundle is EMPTY and that is the first line of scratch/REVIEW.md -- there is
nothing else to review, and naming it is the whole finding. Then run the suite
the way the lander RECORDED running it (its PR.md names the interpreter and
the command): a `pytest` that does not exist on PATH is the WRONG INTERPRETER,
not a red test. Say which interpreter you used.

{text.format(base=base)}

Write scratch/REVIEW.md: first line "REVIEW {key}: <n> findings". Then one
row per finding: severity (blocks / should fix / note), file:line, what is
wrong, the evidence (command and its output, or the caller), and the
smallest fix. A finding you could not confirm is not a finding -- leave it
out or mark it "unconfirmed". No account names, home paths or person names."""})
    if "fix" in stages and revs:
        tasks.append({"id": "zz-fix", "base": f"wave/{land_ref}", "seconds": 14400, "turns": 200,
                      **m(fixer), "after": revs, "collect": [land_ref] + revs,
                      "brief": f"""FIXER of the bundle. Your worktree starts on the lander's branch
(wave/{land_ref}). Four reviewers read it; their reports are in
scratch/arms/zz-rev-*/REVIEW.md, the lander's in scratch/arms/zz-land/.

1. For every finding: confirm it yourself (run the command, read the
   caller). Confirmed -> fix it with a test that is red before the fix;
   not confirmed -> say why in one row. "blocks" findings first.
2. Run the full suite serially once; every failure also on {base}.
3. Ask for real runs (REAL RUNS below) repeating the lander's, so the
   final branch has its own.
4. When resumed: write scratch/PR.md as the final pull-request text -- the
   lander's text updated, plus a table of every review finding and what
   was done with it. No account names, home paths or person names.

Two of the landing rules bind every line you add to wave/{land_ref}:

{_rules_text(FIXER_RULES, land_ref=land_ref, base=base)}
{real_runs}"""})
    return wave.validate_tasks(tasks)
