"""Land a wave member's branch: checked, independently reviewed, opened as a PR.

A member declares its work ready by writing scratch/PR.md in its worktree:
first line the title, the rest the body. Landing runs OUTSIDE the member's
cage (members never hold a forge credential) and refuses a branch unless:
  * it has commits beyond its base and nothing uncommitted outside scratch/;
  * no commit message and no line the branch adds carries a home path, every
    commit is authored as the configured identity (``agent.swarm.git_identity``),
    and no PR_NOTES.md is in the tree;
  * it changes no path under ``agent.swarm.land_forbidden_paths`` (default: the
    ORCA knowledge cards, which need an expert pass first);
  * the test files the branch adds or changes pass on the test interpreter
    (``agent.swarm.test_python``), run serially at low priority;
  * a model of a DIFFERENT family than the one that wrote the branch reviews
    the diff and finds nothing (``review_gate``). That reviewer must have
    passed the calibration corpus, or the branch is refused: the gate never
    approves by default, because a review that did not happen must not be
    reported as one.
A ready branch is pushed as wave/<run>-<id> from the run's repository and
opened as a PR against ``agent.swarm.land_base`` (default main) of
``agent.swarm.github_repo`` (default: the repository's own); PR.md is renamed
PR.md.opened so it is not opened twice, and scratch/LANDING.txt says what
happened. A run whose state says ``"land": "none"`` (or a task that does) is
never opened here; ``"land": "draft"`` opens drafts.

Ported from the orchestrator's own landing script with its review gate and
its tests unchanged in substance; the paths, the identity, the interpreters
and the review endpoint are settings now.

    delfin-agent swarm wave land [RUN ...] [--push]      (default: a dry run)
    delfin-agent swarm wave land --calibrate [MODEL ...] (default: the KIT roster)
    delfin-agent swarm wave land --calibrate-pair [A B]  (default: every pair)
"""
import json
import os
import re
import socket
import subprocess
import sys
from pathlib import Path

HOME = Path.home()


def _swarm_settings() -> dict:
    try:
        from delfin.user_settings import load_settings
        return ((load_settings() or {}).get("agent") or {}).get("swarm") or {}
    except Exception:
        return {}


def _root() -> Path:
    try:
        from . import wave
        return wave.root()
    except Exception:
        return HOME / "delfin-waves"


_S = _swarm_settings()
ROOT = _root()
#: The interpreter the branch's tests run on, and the one for test_opi_*.py
#: (with OPI installed; DELFIN_OPI_REQUIRED=1 turns a missing-OPI skip into a
#: failure). Both default to this interpreter.
PYCI = str(_S.get("test_python") or sys.executable)
PYOPI = str(_S.get("opi_python") or PYCI)
FORBIDDEN_PREFIXES = tuple(_S.get("land_forbidden_paths") or ("delfin/orca_knowledge/cards/",))
#: A home path that still names someone. ``/home/...`` is a redaction
#: and carries no name, but the old pattern matched it because the dot
#: was in the class: agent night3/n3d was refused for a path it had
#: already redacted correctly. A segment has to hold a letter or a
#: digit to be a name.
HOME_PATH = re.compile(r"/home/(?=[A-Za-z0-9_.-]*[A-Za-z0-9])[A-Za-z0-9_.-]+")


def sh(args, cwd=None, timeout=1800):
    return subprocess.run(args, cwd=cwd, capture_output=True, text=True, timeout=timeout)


def git(args, cwd, timeout=1800):
    """git on a member's worktree, through ``safe_git``: the worktree is the
    member's to write, and a program its repository configuration named
    would run here, outside the member's cage, with the person's rights."""
    from delfin.agent import safe_git
    return safe_git.run(list(args), str(cwd), timeout=timeout)


def identity_email() -> str:
    """The email every commit of a landed branch must carry, or ""."""
    ident = _swarm_settings().get("git_identity") or {}
    if isinstance(ident, dict) and ident.get("email"):
        return str(ident["email"]).strip()
    try:
        from .job import git_identity
        return str((git_identity(Path.cwd()) or {}).get("email") or "").strip()
    except Exception:
        return ""


def coauthor_email() -> str:
    """The address the DELFIN agent account's trailer must carry, or "".

    Read from ``agent.commit_coauthor_email`` through DELFIN's own settings
    loader -- the same key the bash executor and the prompt loader use -- so
    the gate and the commits it checks cannot disagree about which account is
    meant.
    """
    try:
        from delfin.user_settings import load_settings
        agent = (load_settings() or {}).get("agent") or {}
        return str(agent.get("commit_coauthor_email", "") or "").strip()
    except Exception:
        return ""


# ---------------------------------------------------------------------------
# The proof a branch must carry, and the signers it must have
# ---------------------------------------------------------------------------
# Three things a landed branch was taken on trust for, and what each one cost:
#
#   * PROOF. A package's claim that its test is red on the base was the
#     author's word. The gate now RUNS the test the package's own decision
#     table names -- on the base (it must fail there) and on the branch (it
#     must pass) -- because a test that passes on the base, or that does not
#     exist, measures nothing, and both were landed.
#   * A REAL RUN. A green test is not a run. The member's real-run rows are
#     read from the events the stage wrote OUTSIDE its cage, not from the
#     summary the member can render itself.
#   * SIGNERS. The COMMITTER was never compared with the identity, so a commit
#     authored as the person and committed as somebody else passed; and a
#     ``Co-Authored-By`` line names an account the forge then counts as a
#     contributor, so it is read now -- and the NUMBER in it is what names the
#     account, which is why a line spelled for the agent with another account's
#     number is refused.

#: One row of the decision table in a member's ``scratch/PR.md``.
_TABLE_ROW = re.compile(r"^\s*\|(?P<cells>.+)\|\s*$")
#: A code span, which is how a path with a ``|`` in it survives a table cell.
_CODE_SPAN = re.compile(r"`([^`]+)`")
#: ``Co-Authored-By: Name <addr>`` in a commit message.
_COAUTHOR = re.compile(r"^\s*Co-Authored-By:\s*(?P<name>[^<]*?)\s*<(?P<addr>[^>]+)>\s*$",
                       re.M | re.I)
#: The number a GitHub noreply address carries: ``<id>+<login>@users.noreply``.
_ACCOUNT_NUMBER = re.compile(r"^\s*(\d+)\+")


def proof_rows(pr_text: str) -> list[dict]:
    """The packages in a ``scratch/PR.md`` decision table, one dict each.

    Input: the text of PR.md. Output: ``[{"package", "test", "realrun"}]`` in
    table order, one entry per data row. Header rows, the ``|---|`` separator,
    and every table that has no ``package`` column are skipped: a wave report
    carries several tables (items, measurements, suites) and only the decision
    table is the package's claim about itself. An empty list means the branch
    wrote no decision table, which the caller refuses.

    The test is taken from a code span when the cell holds one, and from the
    whole cell otherwise.
    """
    rows: list[dict] = []
    header: list[str] | None = None
    pkg_col = test_col = run_col = None
    for line in (pr_text or "").splitlines():
        m = _TABLE_ROW.match(line)
        if not m:
            header, pkg_col, test_col, run_col = None, None, None, None
            continue
        cells = [c.strip() for c in m.group("cells").split("|")]
        if header is None:
            lowered = [c.lower() for c in cells]
            if not any("package" in c for c in lowered):
                continue
            header = lowered

            def find(*words):
                return next((j for j, c in enumerate(lowered)
                             if any(w in c for w in words)), None)
            pkg_col, test_col = find("package"), find("test")
            run_col = find("real run", "realrun", "run")
            continue
        if all(set(c) <= set("-: ") for c in cells):
            continue                                  # the |---|---| separator

        def cell(j):
            return cells[j] if j is not None and j < len(cells) else ""
        span = _CODE_SPAN.search(cell(test_col))
        rows.append({"package": cell(pkg_col).strip("` "),
                     "test": (span.group(1) if span else cell(test_col)).strip("` "),
                     "realrun": cell(run_col).strip("` ")})
    return rows


def _pytest_verdict(python: str, tree: Path, test: str) -> tuple[bool, str]:
    """Whether ``test`` passes in ``tree``, run the way this gate runs tests.

    Returns ``(passed, one line of output)``. A test file absent from the tree
    did not pass: "does not exist" and "is red" are different answers and the
    caller says which one it got.
    """
    if not (Path(tree) / test).is_file():
        return False, f"{test} does not exist in {Path(tree).name}"
    r = subprocess.run(["nice", "-n", "15", "ionice", "-c3", python, "-m", "pytest",
                        "-q", "-p", "no:cacheprovider", str(test)],
                       cwd=str(tree), capture_output=True, text=True, timeout=1800)
    tail = ((r.stdout or "").strip().splitlines() or [""])[-1]
    return r.returncode == 0, tail[:160]


def _archive_tree(repo: Path, ref: str, dest: Path) -> None:
    """A ref of ``repo`` as a plain directory, so a test runs on the BASE
    without disturbing the branch's own worktree."""
    dest.mkdir(parents=True, exist_ok=True)
    from delfin.agent import safe_git
    r = safe_git.run(["archive", ref], str(repo), text=False, timeout=300)
    if r.returncode != 0:
        raise RuntimeError(f"git archive {ref} failed: {r.stderr[-200:]!r}")
    subprocess.run(["tar", "-x", "-C", str(dest)], input=r.stdout, check=True)


def _realrun_rows(run: str, tid: str) -> list[dict]:
    """The member's recorded real-run rows, from the run's event log.

    ``realrun.run`` records one event holding a ``rows`` list per member, and
    it runs outside the member's cage. Read from there and not from
    ``scratch/realrun/summary.md``: that file is text the member can write, and
    a gate that reads the author's own rendering of its run has trusted the
    author again.
    """
    try:
        from . import wave
        got = wave.events(run)
    except Exception:
        return []
    rows: list[dict] = []
    for ev in got.get("events") or []:
        if ev.get("kind") == "realrun" and str(ev.get("id") or "") == str(tid):
            rows.extend(ev.get("rows") or [])
    return rows


def _probe_reached(rows: list[dict], probe: str) -> tuple[bool, str]:
    """Whether ``probe`` ran on the branch side and reached the change.

    A row carrying an ``error`` is a run that started and never got there --
    the import check refused it, the endpoint failed, the budget ran out. It
    is not evidence, and it is named as what it was rather than as an absence.
    """
    named = [r for r in rows if str(r.get("probe") or "") == probe]
    if not named:
        return False, f"no real-run row for {probe}"
    branch = [r for r in named if str(r.get("side") or "") == "branch"]
    if not branch:
        return False, f"{probe} ran on the base only"
    ok = [r for r in branch if not str(r.get("error") or "").strip()]
    if not ok:
        return False, f"the branch probe reached nothing ({str(branch[-1].get('error'))[:120]})"
    return True, ""


def proof_refusal(wt: Path, base: str, *, run: str = "", tid: str = "") -> str:
    """Why this branch carries no proof for one of its packages, or "".

    Input: the member's worktree, its base, and the run/member it belongs to.
    Output: "" when every package in the decision table that ``scratch/PR.md``
    carries names a test the gate RAN red on the base and green on the branch,
    plus a real-run row whose probe reached the change; otherwise one sentence
    naming the package and what is missing. It refuses when anything cannot be
    checked, because a proof that could not be examined is not a proof that
    held.
    """
    pr = Path(wt) / "scratch" / "PR.md"
    try:
        rows = [r for r in proof_rows(pr.read_text()) if r["package"]]
    except OSError as exc:
        return f"no proof table: scratch/PR.md cannot be read ({exc})"
    if not rows:
        return ("no proof table: scratch/PR.md carries no decision table with a "
                "package column (one row per package: the package, its test, its real run)")
    base_tree: str = ""
    try:
        for r in rows:
            pkg, test = r["package"], r["test"]
            if not test:
                return f"{pkg}: the table names no test, so nothing shows the change"
            if not r["realrun"]:
                return f"{pkg}: the table names no real run"
            passed, said = _pytest_verdict(PYCI, Path(wt), test)
            if not passed:
                return f"{pkg}: {test} does not pass on the branch ({said})"
            if not base_tree:
                import tempfile
                base_tree = tempfile.mkdtemp(prefix="delfin-land-base-")
                try:
                    _archive_tree(Path(wt), base, Path(base_tree))
                except Exception as exc:
                    return f"{pkg}: {base} could not be read to run its test ({exc})"
            on_base, said = _pytest_verdict(PYCI, Path(base_tree), test)
            if on_base:
                return (f"{pkg}: {test} passes on {base} too ({said}), so it proves "
                        "nothing about this branch")
            reached, why = _probe_reached(_realrun_rows(run, tid), r["realrun"])
            if not reached:
                return f"{pkg}: no real run reached the change -- {why}"
    finally:
        if base_tree:
            import shutil
            shutil.rmtree(base_tree, ignore_errors=True)
    return ""


def signers_refusal(wt: Path, base: str) -> str:
    """Why a commit of this branch is signed by the wrong party, or "".

    Over ``base..HEAD``: every AUTHOR is the configured identity (as before);
    every COMMITTER is too, because a committer is a second signature on the
    same commit and was compared with nothing; and every ``Co-Authored-By``
    line names the agent account exactly. When that account's address carries
    a number, the number in the line must be the same one: it is the number
    that identifies the account, so a line spelled for the agent with another
    account's number renders in the forge's history as that other account's
    contribution -- and is refused as such, named by its number.
    """
    mail = identity_email()
    if not mail:
        return ("no commit identity is configured (agent.swarm.git_identity): "
                "a branch cannot be shown to be authored as the person")
    authors = set(git(["log", "--format=%ae", f"{base}..HEAD"], wt).stdout.split())
    if authors - {mail}:
        return f"commit author is not the configured identity: {sorted(authors)}"
    committers = set(git(["log", "--format=%ce", f"{base}..HEAD"], wt).stdout.split())
    if committers - {mail}:
        return (f"commit committer is not the configured identity: {sorted(committers)}"
                " (the committer is a signature too)")
    agent = coauthor_email()
    if not agent:
        return "no agent account is configured (agent.commit_coauthor_email)"
    number = _ACCOUNT_NUMBER.match(agent)
    want = number.group(1) if number else ""
    for line in git(["log", "--format=%B", f"{base}..HEAD"], wt).stdout.splitlines():
        m = _COAUTHOR.match(line)
        if not m:
            continue
        got = m.group("addr").strip()
        if got.lower() != agent.lower():
            had = _ACCOUNT_NUMBER.match(got)
            if want and had and had.group(1) != want:
                return (f"a Co-Authored-By line carries account number {had.group(1)}, "
                        f"not the agent account's {want}: {line.strip()[:120]}")
            return ("a Co-Authored-By line is not the DELFIN agent account: "
                    f"{line.strip()[:120]}")
    return ""


def _non_public_pattern() -> re.Pattern:
    """A model, provider or site name in text that is published on the base.

    Built from ``MODEL_FAMILY`` and from the providers a key can be configured
    for, at CALL time: ``MODEL_FAMILY`` is defined far below the check code in
    this module, so a table assembled at import here would be empty.
    """
    words = {f for f in MODEL_FAMILY.values() if f}
    words |= {"kit", "blablador", "openai", "anthropic", "ollama"}
    return re.compile(r"\b(?:%s)\b" % "|".join(sorted((re.escape(w) for w in words),
                                                      key=len, reverse=True)), re.I)


def squash_message(wt: Path, base: str, title: str, msg: str = "") -> tuple[str, str]:
    """The message a squash merge of this branch publishes, and why not, if so.

    Input: the worktree, its base, the PR title and (optionally) a message the
    lander proposes. Output: ``(message, refusal)``.

    The rule, and each clause of it is a defect measured on a landed branch:
    the published message carries EXACTLY ONE ``Co-Authored-By:`` line, naming
    the DELFIN agent account and nothing else -- two lines make two
    contributors of one branch, and none leaves the agent out of the history it
    wrote -- and no word of it names the model, the provider or the site,
    because those belong to the branch and not to the base's history.
    """
    agent = coauthor_email()
    if not agent:
        return "", "no agent account is configured (agent.commit_coauthor_email)"
    body = msg.strip() or title.strip()
    # The lander's own wording, then the subjects of the commits it bundles.
    subjects = [s.strip() for s in
                git(["log", "--format=%s", f"{base}..HEAD"], wt).stdout.splitlines()
                if s.strip() and not s.strip().startswith("Merge ")]
    text = "\n\n".join(p for p in [body, "\n".join(subjects)] if p)
    # Every trailer the commits carry is DROPPED, not carried: one line is the
    # whole of the attribution, and it is the one written below.
    text = re.sub(r"^\s*Co-Authored-By:.*$\n?", "", text, flags=re.M | re.I).strip()
    bad = sorted({m.group(0) for m in _non_public_pattern().finditer(text)})
    if bad:
        return "", ("the message the lander proposes names a model, provider or "
                    f"site ({', '.join(bad[:4])}); the base's history carries none")
    host = (socket.gethostname() or "").split(".")[0].strip()
    if len(host) > 3 and re.search(r"\b%s\b" % re.escape(host), text, re.I):
        return "", f"the message the lander proposes names this site ({host})"
    return f"{text}\n\nCo-Authored-By: DELFIN-Agent <{agent}>", ""


def check(wt: Path, base: str = "origin/main", *, run: str = "", tid: str = "") -> tuple[bool, str]:
    if git(["log", "--oneline", f"{base}..HEAD"], wt).stdout.strip() == "":
        return False, f"no commits beyond {base}"
    # scratch/ is the member's own notes; .delfin/ is the state DELFIN keeps
    # in every workspace it works in, ignored in DELFIN's own repository but
    # untracked in any other (measured 2026-11-24 on a throwaway repository:
    # every member branch was refused for it).
    dirty = [l for l in git(["status", "--short"], wt).stdout.splitlines()
             if not l[3:].startswith(("scratch/", ".delfin/"))]
    if dirty:
        return False, "uncommitted changes: " + ", ".join(l[3:] for l in dirty[:5])
    msgs = git(["log", "--format=%an <%ae>%n%B", f"{base}..HEAD"], wt).stdout
    if HOME_PATH.search(msgs):
        return False, "a commit message carries a home path"
    why = signers_refusal(wt, base)
    if why:
        return False, why
    files = git(["diff", "--name-only", f"{base}...HEAD"], wt).stdout.split()
    if any(f.endswith("PR_NOTES.md") for f in files):
        return False, "PR_NOTES.md is in the tree (it belongs in scratch/PR.md)"
    hit = [f for f in files if f.startswith(FORBIDDEN_PREFIXES)]
    if hit:
        return False, f"changes {hit[0]} (agent.swarm.land_forbidden_paths: an expert pass first)"
    # Only what the branch ADDS: a home path already on main is not this
    # branch's to answer for (it is the personal-data test's on main).
    added = git(["diff", "-U0", f"{base}...HEAD"], wt).stdout
    cur = ""
    for line in added.splitlines():
        if line.startswith("+++ b/"):
            cur = line[6:]
        elif line.startswith("+") and not line.startswith("+++") and HOME_PATH.search(line):
            return False, f"{cur} adds a home path: {line[:120]}"
    tests = [f for f in files if f.startswith("tests/") and f.endswith(".py")
             and Path(f).name.startswith("test_") and (wt / f).is_file()]
    if not tests:
        return False, "the branch adds or changes no test"
    # As in CI: tests/test_opi_*.py run in the OPI job (py3.12 with OPI,
    # DELFIN_OPI_REQUIRED=1 turns a missing-OPI skip into a failure); every
    # other test runs on the main interpreter (py3.11, no OPI).
    opi = [t for t in tests if Path(t).name.startswith("test_opi_")]
    main_tests = [t for t in tests if t not in opi]
    results = []
    for py, group, env in ((PYCI, main_tests, None), (PYOPI, opi, {"DELFIN_OPI_REQUIRED": "1"})):
        if not group:
            continue
        r = subprocess.run(["nice", "-n", "15", "ionice", "-c3", py, "-m", "pytest", "-q",
                            "-p", "no:cacheprovider", *group], cwd=wt, capture_output=True,
                           text=True, timeout=3600,
                           env={**os.environ, **env} if env else None)
        tail = (r.stdout or "").strip().splitlines()[-1:] or [""]
        if r.returncode != 0:
            return False, f"tests fail ({Path(py).parent.parent.name}): {tail[0][:200]}"
        # A run in which nothing passed proves nothing (all skipped).
        if " passed" not in tail[0]:
            return False, f"no test passed ({tail[0][:120]}) on {Path(py).parent.parent.name}"
        results.append(tail[0][:120])
    # The proof each package in scratch/PR.md claims, checked AFTER the tests
    # the branch itself carries have been shown green: a branch whose own tests
    # do not run has nothing to prove anything with.
    why = proof_refusal(wt, base, run=run, tid=tid)
    if why:
        return False, why
    return True, " | ".join(results + ["every package proved"])


# ---------------------------------------------------------------------------
# The independent review gate
# ---------------------------------------------------------------------------
# Everything above is checked by the author's own work: the tests the gate
# runs are the tests the author wrote. The gate below hands the DIFF to a
# model of a different family than the one that wrote it, before the push, and
# refuses the branch when that reviewer names a defect with a file and a line.
#
# Three properties carry it, and each one is a defect class measured in this
# installation rather than a precaution:
#
#   * It is IN the path. A review script someone has to remember to run is a
#     mechanism built and never mounted, which is the most common defect found
#     here. A refusal writes scratch/LANDING.txt and leaves the PR unopened,
#     exactly like every refusal above it.
#   * It fails CLOSED. No reviewer of another family, no calibrated reviewer,
#     an endpoint that does not answer, a verdict that will not parse, a diff
#     too large to review in one prompt: every one of those refuses the
#     branch. A gate that approves when its reviewer is missing reports a
#     review that did not happen, which is worse than no gate at all.
#   * It is CALIBRATED before it is trusted. A reviewer that passes everything
#     is indistinguishable from no reviewer, and one that refuses everything
#     looks like diligence while stopping all work. ``--calibrate`` runs a
#     corpus of real cases whose verdict is known (review_corpus/manifest.json)
#     and records which model passed it, when, and over which corpus; a model
#     with no fresh pass is not used on real work.
#
# NOT measured, deliberately: the NUMBER of reviews performed. That number is
# driven to any value by running more reviewers and is not a measure of
# quality. The per-run record below holds the verdicts, the refusal rate
# BESIDE its absolute counts, and the calibration state, as diagnosis only.
# Do not add a review counter here, and do not make any of these a target.

#: The reviewer's endpoint: an OpenAI-compatible chat API. The default is the
#: local KIT gateway, because an external provider would send this
#: repository's diffs to a third party, a decision nobody made for a gate
#: that runs unattended.
def _kit_base_url() -> str:
    from delfin.agent.model_capabilities import _KIT_BASE_URL
    return _KIT_BASE_URL


KIT_BASE_URL = str(_S.get("review_base_url") or _kit_base_url())
REVIEW_KEY_NAME = str(_S.get("review_key_name") or "KIT_TOOLBOX_API_KEY")
#: A directory holding a corpus of one's own (``agent.swarm.review_corpus``),
#: or None for the corpus shipped with DELFIN. That one is an archive
#: (``review_corpus.tar.gz``), read in memory: as text under delfin/ its
#: briefs and diffs were read as code by every repository-wide scanner
#: (names, callers, removed predicates), measured 2026-11-24.
REVIEW_CORPUS = Path(str(_S["review_corpus"])).expanduser() if _S.get("review_corpus") else None
_SHIPPED_CORPUS = Path(__file__).resolve().parent / "review_corpus.tar.gz"
_corpus_cache: dict = {}


def _corpus_read(name: str) -> bytes:
    """One file of the calibration corpus, by its name in the manifest."""
    if REVIEW_CORPUS is not None:
        return (REVIEW_CORPUS / name).read_bytes()
    if not _corpus_cache:
        import tarfile
        with tarfile.open(_SHIPPED_CORPUS, "r:gz") as tar:
            for member in tar.getmembers():
                if member.isfile():
                    _corpus_cache[Path(member.name).name] = tar.extractfile(member).read()
    if name not in _corpus_cache:
        raise FileNotFoundError(f"{name} is not in the shipped calibration corpus")
    return _corpus_cache[name]

REVIEW_CALIBRATION = ROOT / "review_calibration.json"
REVIEW_LOG = ROOT / "review_log.jsonl"
#: A calibration result older than this is not evidence about the model that
#: answers now: the KIT ids are deployment aliases and what is served behind
#: one changes without notice (the roster lost four models in a day on
#: 2026-09-07). Six hours, so one night's wave runs on one calibration.
CALIBRATION_MAX_AGE_S = float(os.environ.get("DELFIN_REVIEW_CALIBRATION_MAX_AGE_S", 6 * 3600))
REVIEW_TIMEOUT_S = float(os.environ.get("DELFIN_REVIEW_TIMEOUT_S", 300))
#: The reviewer's answer is ~200 tokens, but a reasoning model spends its
#: budget BEFORE writing any of it. Measured on kit.glm-5.3-flash with one
#: corpus case: 4,000 completion tokens, all 4,000 of them reasoning tokens,
#: finish_reason "length" and an EMPTY content -- which reads exactly like a
#: reviewer that said nothing. Two of six calibration cases were scored
#: "unusable" for that alone, which is the gate blaming the model for its own
#: token limit.
REVIEW_MAX_TOKENS = int(os.environ.get("DELFIN_REVIEW_MAX_TOKENS", 16000))
#: A diff longer than this is refused rather than truncated: a review of half
#: a diff that reports on all of it is the instrument this gate exists against.
MAX_REVIEW_DIFF = 160_000
#: How much of the brief the reviewer gets. The briefs here run to ~7k chars,
#: of which the tail is the wave's shared delivery rules.
MAX_REVIEW_BRIEF = 12_000

#: Which family serves a model id. Independence has to be mechanical, so this
#: is a lookup and not a guess: an id that is NOT here has an unknown family,
#: and then every reviewer is refused, because a reviewer cannot be shown to
#: differ from an author nobody can name.
#: The blablador ids are server-side aliases; what each points at is recorded
#: in delfin/agent/model_capabilities.py beside the measurement that ranked
#: them (alias-code -> Qwen3.8-27B, alias-huge -> MiMo-V2.6-Pro,
#: alias-fast -> GPT-OSS-120b, alias-large -> Qwen3.8-Flash).
MODEL_FAMILY: dict[str, str] = {
    "kit.glm-5.3-flash": "glm",
    "kit.glm-5.3": "glm",
    "kit.deepseek-v4-flash": "deepseek",
    "kit.mistral-small-4-119b-a8b": "mistral",
    "deepseek-v4.1-flash-at-hzdr": "deepseek",
    "alias-code": "qwen",
    "alias-large": "qwen",
    "alias-qwen3.8-flash-next": "qwen",
    "alias-huge": "mimo",
    "alias-mimo-v2.6-pro": "mimo",
    "alias-fast": "gpt-oss",
    "alias-glm-5.3-flash": "glm",
    "alias-kimi-k3-1m": "kimi",
    "alias-muse-spark": "muse",
}

REVIEW_SYSTEM = """\
You review one pull request written by another coding agent, independently.

You are shown exactly two things: the task brief that agent was given, and
the diff of its branch against main. You are deliberately NOT shown its pull
request text, its commit messages or its log, so that you judge the code and
not the author's account of it.

Report the single strongest finding you can support, or pass. These are the
defect shapes measured in this repository, in order of how often they occur:

1. built but never mounted -- a capability, guard, report or scan the product
   never reaches. The parts are each tested and the guarantee is wired to
   nothing. In a diff this looks like a new module plus a test for it, with no
   existing file changed to call it.
2. a test that only its own module calls -- the new tests pass and nothing in
   the product drives the code. A test that calls the new function directly
   proves the function runs, not that the product runs it.
3. a mechanism that answers a neighbouring question -- `.exists()` asked for
   "is a checkout", a file existing asked for "a file was written", a
   per-call cost asked for a behaviour, a count asked for an effect.
4. a check that observed nothing -- green because the fixture had nothing to
   be wrong about: an assertion on an empty list, a skip that stands in for a
   result, a guard tested with input it cannot refuse.
5. a deletion that has nothing to do with the stated task -- a check, a
   branch or a test removed in a change about something else.
6. a number the diff or the brief CONTRADICTS -- a count, a rate or a
   threshold that cannot be true of the change in front of you: a figure the
   diff's own code or data does not yield, a figure the text says THIS change
   measured with no run or fixture here producing it, or a figure the brief
   gives differently. Report the contradiction, not the absence of a source.
   A docstring citing an earlier measurement as the motivation for the change
   is house style here and is required; you cannot see its source and it is
   not a finding. "I cannot verify this number" is not a finding; "this
   number cannot be right, and here is what the diff says instead" is.

Rules for your answer:

* A finding names a file and a line that appear in the diff. A finding
  without them is worthless here and is discarded.
* You see only the diff. Do not report that a symbol is missing from the
  repository merely because the diff does not show its definition. Report it
  only when the diff itself ADDS the dependence -- a new call, in a path that
  runs -- and adds no definition for it; then say exactly that.
* Do not report style, naming, formatting, type annotations, docstring
  wording, or a wish for more tests.
* Most pull requests in this repository are sound. Refuse only when you can
  point at a file and a line and say what breaks or what is not reached.

Answer with one JSON object and nothing else, no prose around it:

{"verdict": "pass" | "refuse",
 "findings": [{"file": "path/in/the/diff.py", "line": 123,
               "claim": "what is wrong, in one sentence",
               "why_it_matters": "what it costs if it ships"}]}

"pass" carries "findings": []. "refuse" carries at least one finding with a
file and a line.
"""


def model_family(model: str) -> str:
    """The family serving ``model``, or "" when it is not known.

    Input: a model id as the run state or the roster spells it. Output: a
    family tag. An empty return is not "no family": it is "unknown", and every
    caller treats it as a refusal.
    """
    return MODEL_FAMILY.get((model or "").strip(), "")


def kit_reviewer_models() -> tuple[str, ...]:
    """The local KIT roster, read from the product instead of copied into here.

    A reviewer runs on KIT only. An external provider would send this
    repository's diffs to a paid third party, which is a decision nobody made
    for a gate that runs unattended. A roster that cannot be read returns
    empty and the gate refuses -- the fail-closed direction.
    """
    served: set = set()
    try:
        from delfin.agent.model_capabilities import MODEL_PREFERENCE_BY_PROVIDER as pref
        from delfin.agent.model_capabilities import kit_models_served
        from delfin.agent.credentials import load_credential
        try:
            served = kit_models_served(
                api_key=load_credential("KIT_TOOLBOX_API_KEY") or "",
                fetch=True)
        except Exception:
            served = set()
    except Exception:
        return ()
    # An id whose family is unknown cannot serve as a reviewer either: its
    # independence from the author is unprovable.
    ranked = tuple(m for m in pref.get("kit", ()) if model_family(m))
    # The product's list is a RANKED FALLBACK, not a list of served models:
    # it deliberately keeps an id the endpoint has stopped serving as the
    # next candidate in case it returns, and records which. Calibrating the
    # ranking as though it were the roster spent 24 calls on kit.glm-5.3 on
    # 2026-11-21 -- "400 Model not found" on every one -- and then reported
    # that model as a reviewer with a catch rate.
    #
    # An empty ``served`` means the listing could not be asked, NOT that
    # nothing is served, so it falls through to the full ranking rather than
    # to no reviewer at all: a roster that cannot be read must not silently
    # become an empty one.
    if not served:
        return ranked
    live = tuple(m for m in ranked if m in served
                 or m.split(".", 1)[-1] in served)
    return live or ranked


def author_models(agent: dict) -> list[str]:
    """Every model that may have written this branch.

    Input: one ``state["agents"][tid]`` record. Output: the model ids to treat
    as the author side.

    ``task["model"]`` is the model the agent runs on NOW. A resumed agent had
    its model re-picked at every resume (pick_model.py writes the same field),
    and the state keeps only the last choice, so the field alone can name a
    family that did not write the code under review. The agent's own log is
    read as well, for any id the family table knows. A second name there can
    only make the reviewer choice stricter, which is the direction a
    fail-closed gate wants.
    """
    out: list[str] = []
    m = ((agent.get("task") or {}).get("model") or "").strip()
    if m:
        out.append(m)
    try:
        text = Path(agent["log"]).read_text(errors="replace")
    except Exception:
        text = ""
    for cand in MODEL_FAMILY:
        if cand in text and cand not in out:
            out.append(cand)
    return out


def independent_candidates(authors: list[str]) -> tuple[list[str], str]:
    """Reviewer models of a different family than every author model.

    Output: (candidates, why-not). An empty list always carries a reason, and
    a reason is always a refusal -- never a default approval.
    """
    if not authors:
        return [], "the run state names no author model"
    families = set()
    for m in authors:
        fam = model_family(m)
        if not fam:
            return [], f"the author model {m} is of an unknown family"
        families.add(fam)
    roster = kit_reviewer_models()
    if not roster:
        return [], "the KIT roster could not be read from the product"
    cands = [m for m in roster if model_family(m) not in families]
    if not cands:
        return [], ("no KIT model outside the author families "
                    + "/".join(sorted(families)))
    return cands, ""


def kit_key() -> str:
    """The KIT key, from the environment or the product's own credential store.

    Read only. Nothing here writes to the store.
    """
    k = os.environ.get(REVIEW_KEY_NAME, "").strip()
    if k:
        return k
    try:
        from delfin.agent.credentials import load_credential
        return (load_credential(REVIEW_KEY_NAME) or "").strip()
    except Exception:
        return ""


def kit_chat(model: str, messages: list[dict], timeout: float | None = None) -> str:
    """One chat completion against the KIT gateway. Raises on anything else."""
    import urllib.request

    key = kit_key()
    if not key:
        raise RuntimeError(f"no {REVIEW_KEY_NAME} available")
    body = json.dumps({"model": model, "messages": messages, "temperature": 0,
                       "max_tokens": REVIEW_MAX_TOKENS}).encode()
    req = urllib.request.Request(
        KIT_BASE_URL + "/chat/completions", data=body,
        headers={"Authorization": "Bearer " + key, "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout or REVIEW_TIMEOUT_S) as r:
        data = json.loads(r.read())
    choice = data["choices"][0]
    text = (choice["message"].get("content") or "").strip()
    if not text:
        # An empty content is not a verdict, and WHY it is empty decides
        # whether the gate should be re-run or the model dropped: a reply cut
        # off at the token limit is the caller's defect, a genuinely empty
        # answer is the model's. Say which, instead of "no verdict object".
        usage = (data.get("usage") or {}).get("completion_tokens_details") or {}
        raise RuntimeError(
            f"empty reply (finish_reason={choice.get('finish_reason')!r}, "
            f"completion_tokens={(data.get('usage') or {}).get('completion_tokens')}, "
            f"reasoning_tokens={usage.get('reasoning_tokens')})")
    return text


def build_review_prompt(diff: str, brief: str) -> list[dict]:
    """The reviewer's whole input: the task brief and the diff, and nothing else.

    Not the author's PR.md, not LANDING.txt, not its log, not its commit
    messages: the reviewer gets the code and the job, never the author's
    account of what the code achieves. tests/test_a_wave_branch_lands_through_the_review_gate.py asserts over
    this function AND over the real gate path that none of those reach it.
    """
    brief = (brief or "").strip()[:MAX_REVIEW_BRIEF]
    return [
        {"role": "system", "content": REVIEW_SYSTEM},
        {"role": "user", "content":
            "The task brief the author was given:\n\n"
            "<<<BRIEF\n" + brief + "\nBRIEF\n\n"
            "The branch's diff against main (git diff origin/main...HEAD):\n\n"
            "<<<DIFF\n" + diff + "\nDIFF\n"},
    ]


def _json_objects(text: str):
    """Every brace-balanced slice of ``text``, outermost first, in order."""
    depth = 0
    start = -1
    for i, ch in enumerate(text):
        if ch == "{":
            if depth == 0:
                start = i
            depth += 1
        elif ch == "}" and depth:
            depth -= 1
            if depth == 0 and start >= 0:
                yield text[start:i + 1]


def parse_verdict(text: str):
    """The verdict object in a reply, or None when there is none.

    A model that restates the schema before answering leaves two objects; the
    LAST one carrying a "verdict" key is the answer.
    """
    found = None
    for chunk in _json_objects(text or ""):
        try:
            obj = json.loads(chunk)
        except Exception:
            continue
        if isinstance(obj, dict) and "verdict" in obj:
            found = obj
    return found


def checkable_findings(verdict: dict) -> list[dict]:
    """The findings that can be checked by a person: a file, a line, a claim.

    A refusal without a file and a line is not a refusal. "This looks risky"
    cannot be falsified, and an unfalsifiable refusal would stall the weekend
    while looking like diligence.
    """
    out = []
    for f in (verdict or {}).get("findings") or []:
        if not isinstance(f, dict):
            continue
        name = str(f.get("file") or "").strip()
        line = f.get("line")
        if isinstance(line, str):
            m = re.search(r"\d+", line)
            line = int(m.group()) if m else None
        if isinstance(line, bool) or not isinstance(line, int):
            continue
        claim = str(f.get("claim") or "").strip()
        if name and line > 0 and claim:
            out.append({"file": name, "line": line, "claim": claim,
                        "why_it_matters": str(f.get("why_it_matters") or "").strip()})
    return out


def review_diff(model: str, diff: str, brief: str, attempts: int = 2) -> dict:
    """Ask ``model`` for a verdict on ``diff``. Never raises.

    Output: ``{"verdict": "pass"|"refuse"|"unusable", "findings": [...],
    "notes": [...], "calls": n}``. "unusable" covers an endpoint that did not
    answer, a reply that did not parse, and a refusal with no file and line --
    all of which leave the branch unopened, for a person.
    """
    messages = build_review_prompt(diff, brief)
    notes: list[str] = []
    calls = 0
    for _ in range(max(1, attempts)):
        calls += 1
        try:
            raw = kit_chat(model, messages)
        except Exception as e:
            notes.append(f"{type(e).__name__}: {str(e)[:160]}")
            continue
        got = parse_verdict(raw)
        if got is None:
            notes.append("the reply held no verdict object: " + raw[:160].replace("\n", " "))
            continue
        verdict = str(got.get("verdict") or "").strip().lower()
        if verdict == "pass":
            return {"verdict": "pass", "findings": checkable_findings(got),
                    "notes": notes, "calls": calls}
        if verdict == "refuse":
            findings = checkable_findings(got)
            if findings:
                return {"verdict": "refuse", "findings": findings,
                        "notes": notes, "calls": calls}
            notes.append("a refusal that names no file and line")
            continue
        notes.append(f"unknown verdict {verdict!r}")
    return {"verdict": "unusable", "findings": [], "notes": notes, "calls": calls}


# --- calibration ------------------------------------------------------------


def corpus_manifest() -> dict:
    """The calibration corpus as it is pre-registered on disk.

    ``manifest.json`` carries the cases, how they are scored, where each label
    comes from, and a dated log of every change to the set. The log is not
    decoration: a corpus that is quietly edited until a reviewer passes it
    measures nothing, so what moved and why is written down beside it.
    """
    return json.loads(_corpus_read("manifest.json").decode())


def corpus_cases() -> list[dict]:
    """The scored cases: real diffs whose verdict is known."""
    return corpus_manifest()["cases"]


def corpus_digest() -> str:
    """A digest over the whole instrument a calibration was measured with.

    A recorded pass refers to one corpus AND one way of asking. Without this,
    editing the corpus -- or softening the question in REVIEW_SYSTEM, or
    changing the token budget that decides whether a reasoning model answers
    at all -- until a reviewer passes would leave the old pass standing, and
    the recorded measurement would be of an instrument that no longer exists.
    """
    import hashlib

    man = corpus_manifest()
    h = hashlib.sha256()
    h.update(json.dumps([REVIEW_SYSTEM, REVIEW_MAX_TOKENS]).encode())
    # The scored set and its rules, not the whole file: adding a dated log
    # entry explains a corpus, it does not change it.
    h.update(json.dumps([man.get("version"), man.get("scoring"), man["cases"]],
                        sort_keys=True).encode())
    for case in corpus_cases():
        for key in ("diff", "brief"):
            h.update(case["id"].encode())
            h.update(_corpus_read(case[key]))
    return h.hexdigest()


def _read_calibration() -> dict:
    """The calibration cache, strictly. A file that will not parse RAISES.

    A lenient read in a writer destroys data: returning {} here and writing
    would erase every other model's result on the first corrupt byte.
    """
    if not REVIEW_CALIBRATION.is_file():
        return {}
    text = REVIEW_CALIBRATION.read_text()
    data = json.loads(text)
    if not isinstance(data, dict):
        raise ValueError(f"{REVIEW_CALIBRATION} is not an object")
    return data


def _write_calibration(cache: dict) -> None:
    """Replace the calibration file atomically.

    Written once and called from both the model and the pair calibration:
    two copies of one write drift, and this one has an invariant worth
    keeping in a single place -- the file is replaced, never truncated in
    place, so an interrupted run leaves the previous calibration intact
    rather than an empty object that every reader would read as "nothing is
    calibrated".
    """
    tmp = REVIEW_CALIBRATION.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(cache, indent=1) + "\n")
    tmp.replace(REVIEW_CALIBRATION)


def calibration_record(model: str):
    """The fresh, same-corpus, passed calibration for ``model``, or None."""
    import time

    try:
        rec = _read_calibration().get(model)
    except Exception:
        return None
    if not isinstance(rec, dict) or not rec.get("passed"):
        return None
    if rec.get("corpus_digest") != corpus_digest():
        return None
    if time.time() - float(rec.get("ts") or 0) > CALIBRATION_MAX_AGE_S:
        return None
    return rec


def review_pair(model_a: str, model_b: str, diff: str,
                brief: str) -> dict:
    """Two reviewers of different families on one diff, as four outcomes.

    Input: two model ids, the diff and the brief. Output::

        {"verdict": "pass" | "refuse" | "split" | "no_answer",
         "a": <review_diff result>, "b": <...>, "findings": [...],
         "calls": n, "why": "..."}

    Semantics, and the reason it is four and not two:

    * ``pass``   -- BOTH passed. Measured on corpus v5: 4 such cases, all
      of them sound branches.
    * ``refuse`` -- BOTH refused. 2 such cases, both real defects. The
      findings of both are kept, a's first.
    * ``split``  -- one passed and one refused. 4 such cases, and they are
      exactly the ones a person should see: of the four, three were real
      defects one model saw and the other did not, and one was a sound
      branch one model refused.
    * ``no_answer`` -- one gave no usable answer at all. 2 such cases, both
      of them the same instrument failure (an empty reply at
      finish_reason "length" with the whole completion budget spent on
      reasoning), NOT a judgement. It is kept separate from ``split``
      because the remedy differs: a split needs a person, a missing answer
      needs a retry or a bigger budget.

    Why a pair at all: measured on v5, NEITHER model passes the corpus
    alone -- kit.glm-5.3-flash catches every positive it is scored on and
    false-refuses 2 of 6 negatives, kit.mistral-small-4-119b-a8b
    false-refuses 0 of 6 and catches 2 of 6. Opposite error profiles. The
    pair decided 6 of 12 cases and got **none** of them wrong; the other 6
    it handed over. A gate that abstains half the time and is never wrong
    is usable, and one that is right 70 % of the time alone is not.
    """
    a = review_diff(model_a, diff, brief)
    b = review_diff(model_b, diff, brief)
    out = {"a": a, "b": b, "model_a": model_a, "model_b": model_b,
           "calls": a["calls"] + b["calls"],
           "findings": list(a["findings"]) + list(b["findings"])}
    va, vb = a["verdict"], b["verdict"]
    if "unusable" in (va, vb):
        who = model_a if va == "unusable" else model_b
        out["verdict"] = "no_answer"
        out["why"] = (f"{who} gave no usable answer, so the pair cannot "
                      "agree or disagree; this is an instrument failure, "
                      "not a verdict")
    elif va == vb == "pass":
        out["verdict"] = "pass"
        out["why"] = f"{model_a} and {model_b} both pass"
    elif va == vb == "refuse":
        out["verdict"] = "refuse"
        out["why"] = f"{model_a} and {model_b} both refuse"
    else:
        out["verdict"] = "split"
        out["why"] = (f"{model_a} says {va}, {model_b} says {vb}; a person "
                      "decides, and both findings are below")
    return out


def pair_bounds(rows: list[dict]) -> dict:
    """What a pair calibration's own rows support.

    The pass rule for a PAIR is not the pass rule for a model. A model is
    asked to be right on every scored case. A pair is allowed to abstain --
    that is the whole point of it -- so it is asked for something else:

    * it must DECIDE nothing wrongly (``wrong`` is 0), and
    * among the cases it decided, both directions must appear, or an
      abstaining rubber stamp would pass: a pair that only ever agrees on
      passes has not been shown to catch anything.

    The abstention rate is reported, never optimised. Driving it down is
    what a Goodhart-minded operator would do next, and it is the one number
    here whose denominator the pair itself controls: a pair that abstains on
    everything has a perfect wrong-rate.
    """
    decided = [r for r in rows if r["pair"] in ("pass", "refuse")]
    wrong = [r for r in decided if not r["ok"]]
    abstained = [r for r in rows if r["pair"] in ("split", "no_answer")]
    return {
        "n": len(rows),
        "decided": len(decided),
        "wrong": len(wrong),
        "wrong_cases": [r["case"] for r in wrong],
        "abstained": len(abstained),
        "abstained_split": sum(r["pair"] == "split" for r in rows),
        "abstained_no_answer": sum(r["pair"] == "no_answer" for r in rows),
        "agreed_pass": sum(r["pair"] == "pass" for r in decided),
        "agreed_refuse": sum(r["pair"] == "refuse" for r in decided),
        "correct_when_agreed_at_least": (
            one_sided_bounds(len(decided) - len(wrong), len(decided))["lower"]),
    }


def _case_caught(case: dict, result: dict) -> tuple[bool, str]:
    """Did the reviewer answer this control, for the right reason?

    A positive control is only caught when the refusal points at a file the
    control is about (and, where the case names one, says the thing). A
    reviewer that refuses everything would otherwise score every positive.
    """
    want = case["expect"]
    got = result["verdict"]
    if got != want:
        return False, f"expected {want}, got {got}"
    if want != "refuse":
        return True, "passed it"
    files = {Path(f["file"]).name for f in result["findings"]}
    expect_files = {Path(f).name for f in case.get("expect_files") or []}
    if expect_files and not (files & expect_files):
        return False, ("refused, but about " + ", ".join(sorted(files) or ["nothing"])
                       + " instead of " + ", ".join(sorted(expect_files)))
    claims = " ".join(f["claim"] + " " + f["why_it_matters"] for f in result["findings"])
    missing = [k for k in case.get("expect_claim") or [] if k not in claims]
    if missing:
        return False, "refused, but named none of " + ", ".join(missing)
    return True, "caught it"


def _pct(x) -> str:
    """A rate as a percentage, or "n/a" when there is nothing to report."""
    return "n/a" if x is None else f"{100.0 * x:.1f} %"


def scores_for(case: dict, model: str) -> bool:
    """Whether ``case`` counts towards ``model``'s pass or fail.

    Input: one corpus case and a model id. Output: False only when the case
    carries ``found_by: <model>`` -- the defect in it was reported by that
    model, so scoring the model on it would be an instrument measuring
    itself. The case is still RUN and still recorded for that model; it is
    left out of the pass/fail, and it scores every other model normally.
    """
    return (case.get("found_by") or "") != model


def one_sided_bounds(ok: int, n: int, alpha: float = 0.05) -> dict:
    """Exact (Clopper-Pearson) one-sided 95 % bounds on a corpus result.

    Input: ``ok`` controls answered out of ``n`` of one direction. Output:
    ``{"n", "ok", "lower", "upper"}`` -- the bounds on the true rate of
    answering that direction correctly. ``upper`` is what matters for the
    negatives (1 - lower bound on passing = the worst case false-refuse
    rate), ``lower`` for the positives (the worst case catch rate).

    Why it is in here: the point of the corpus is the bound, not the pass.
    4/4 negatives bound the true false-refuse rate only at 52.7 % and 2/2
    positives bound the catch rate only at 22.4 %, which rules out a rubber
    stamp and a blanket refuser and nothing else. A record that carries the
    pass without the bound invites reading "passed calibration" as "catches
    defects", which n=6 does not support.

    Computed by inverting the binomial tail, so there is no dependency and
    no factorial: lower solves P[X >= ok] = alpha, upper solves
    P[X >= ok+1] = 1 - alpha, both tails increasing in p, with the
    degenerate ends pinned at 0 and 1.
    """
    if n <= 0:
        return {"n": 0, "ok": 0, "lower": None, "upper": None}

    def _cdf_at_least(p: float, k: int, n: int) -> float:
        """P[X >= k] for X ~ Bin(n, p), summed without a factorial."""
        if k <= 0:
            return 1.0
        if k > n:
            return 0.0
        term = (1.0 - p) ** n
        total = 0.0
        for i in range(0, n + 1):
            if i >= k:
                total += term
            if i < n:
                if p >= 1.0:
                    term = 1.0 if i + 1 == n else 0.0
                else:
                    term = term * (n - i) / (i + 1) * p / (1.0 - p)
        return min(1.0, max(0.0, total))

    def _solve(k, target):
        """The p at which P[X >= k] first reaches ``target``; increasing in p."""
        lo, hi = 0.0, 1.0
        for _ in range(200):
            mid = (lo + hi) / 2
            if _cdf_at_least(mid, k, n) < target:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2

    lower = 0.0 if ok == 0 else _solve(ok, alpha)
    upper = 1.0 if ok == n else _solve(ok + 1, 1 - alpha)
    return {"n": n, "ok": ok, "lower": lower, "upper": upper}


def calibration_bounds(rows: list[dict]) -> dict:
    """The bounds a calibration result supports, from its scored rows only.

    Output: ``{"positives": {...}, "negatives": {...},
    "catch_rate_at_least": x, "false_refuse_rate_at_most": y}`` -- the two
    numbers to report beside the pass, each with its absolute counts beside
    it (a rate alone is a target in disguise).
    """
    scored = [r for r in rows if r.get("counted", True)]
    # A row the endpoint never answered holds no judgement. Counting it as a
    # wrong judgement puts the gate's own failures -- a model id the endpoint
    # does not serve, a reply that spent its whole completion budget on
    # reasoning and arrived empty -- into the numerator of a rate about the
    # model's reading of code. Measured 2026-11-21: kit.glm-5.3 answered 0 of
    # 12 cases ("Model not found" on all 24 calls) and was reported at "catch
    # rate at least 0.0 %", and kit.glm-5.3-flash lost 2 of 12 to an empty
    # reply at finish_reason "length" with reasoning_tokens == max_tokens --
    # one of them a case it had answered correctly at the same budget in the
    # previous run, so the loss is not a property of the case either.
    #
    # Unanswered rows are therefore OUT of the rate and reported beside it,
    # never folded away: ``passed`` still requires every scored row to be ok,
    # so an unanswered case continues to block the pass. The rate says what
    # was measured about judgement; ``not_measured`` says how much judgement
    # was never obtained.
    answered = [r for r in scored if r["verdict"] != "unusable"]
    pos = [r for r in answered if r["expect"] == "refuse"]
    neg = [r for r in answered if r["expect"] == "pass"]
    pb = one_sided_bounds(sum(r["ok"] for r in pos), len(pos))
    nb = one_sided_bounds(sum(r["ok"] for r in neg), len(neg))
    missing = [r for r in scored if r["verdict"] == "unusable"]
    return {"positives": pb, "negatives": nb,
            "not_measured": {
                "n": len(missing),
                "of": len(scored),
                "positives": sum(r["expect"] == "refuse" for r in missing),
                "negatives": sum(r["expect"] == "pass" for r in missing),
                "cases": [r["case"] for r in missing]},
            "measurable": bool(answered),
            "catch_rate_at_least": pb["lower"],
            "false_refuse_rate_at_most": (None if nb["lower"] is None
                                          else 1.0 - nb["lower"])}


def pair_calibration_key(model_a: str, model_b: str) -> str:
    """The cache key for a pair. Order-free: the pair is symmetric."""
    return "pair:" + "+".join(sorted((model_a, model_b)))


def pair_passes(bounds: dict) -> bool:
    """May this pair decide anything? Input: ``pair_bounds``. Output: bool.

    The rule for a pair is not the rule for a model, and the difference is
    the whole design. A model is asked to be right on every scored case. A
    pair is ALLOWED to abstain, so it is asked for two other things:

    1. it must have DECIDED nothing wrongly -- a shared false alarm is a
       wrong decision and not a check, and that is measured, not assumed:
       kit.deepseek-v4-flash + kit.glm-5.3-flash are two different families
       and the pair still failed, because both refuse the same sound
       branch; and
    2. it must have decided in BOTH directions, because the two degenerate
       pairs each score a clean sheet otherwise -- one that abstains on
       everything has decided nothing wrongly, and one that only ever
       agrees on passes has not been shown to catch anything at all.

    Written as its own function so each clause can be removed in a
    mutation and a named test fails.
    """
    return bool(bounds["decided"] > 0
                and bounds["wrong"] == 0
                and bounds["agreed_pass"] > 0
                and bounds["agreed_refuse"] > 0)


def calibrate_pair(model_a: str, model_b: str, out=print,
                   reuse: bool = True) -> dict:
    """Score a PAIR of reviewers on the corpus, as four outcomes.

    Input: two model ids of different families. Output: a record, also
    written to the calibration cache under ``pair:<a>+<b>``.

    ``reuse`` (default on) derives the pair's rows from each model's own
    stored calibration when both were measured on THIS corpus digest,
    costing no model calls at all; the record then carries
    ``derived_from`` with both source timestamps, because a derived
    measurement is not an independent one and the record must not hide
    that. With ``reuse=False`` both models are asked again, 2 calls per
    case.

    A case scores for the pair only when it scores for BOTH models: a case
    authored from one model's own report is not a test of that model, and
    the pair inherits the exclusion rather than averaging it away.
    """
    import time

    cases = corpus_cases()
    digest = corpus_digest()
    cache = _read_calibration()
    derived_from = None

    if reuse:
        a_rec, b_rec = cache.get(model_a), cache.get(model_b)
        usable = all(isinstance(r, dict) and r.get("corpus_digest") == digest
                     and r.get("cases") for r in (a_rec, b_rec))
        if usable:
            by_case = {}
            for rec, who in ((a_rec, "a"), (b_rec, "b")):
                for row in rec["cases"]:
                    by_case.setdefault(row["case"], {})[who] = row
            derived_from = {model_a: a_rec.get("iso"), model_b: b_rec.get("iso")}
            out(f"deriving the pair from the stored rows of {model_a} "
                f"({a_rec.get('iso')}) and {model_b} ({b_rec.get('iso')}) "
                "-- no model calls")
        else:
            out("no usable stored rows for both models on this corpus; asking "
                "both models")
            reuse = False

    rows, calls = [], 0
    for case in cases:
        if reuse:
            got = by_case.get(case["id"]) or {}
            if not (got.get("a") and got.get("b")):
                continue
            va, vb = got["a"]["verdict"], got["b"]["verdict"]
            findings = (list(got["a"].get("findings") or [])
                        + list(got["b"].get("findings") or []))
        else:
            diff = _corpus_read(case["diff"]).decode()
            brief = _corpus_read(case["brief"]).decode()
            res = review_pair(model_a, model_b, diff, brief)
            calls += res["calls"]
            va, vb = res["a"]["verdict"], res["b"]["verdict"]
            findings = res["findings"]

        if "unusable" in (va, vb):
            pair = "no_answer"
        elif va == vb == "pass":
            pair = "pass"
        elif va == vb == "refuse":
            pair = "refuse"
        else:
            pair = "split"

        # A pair's decision is right when it matches the control. An
        # abstention is neither right nor wrong -- it is the pair declining
        # to decide, which is the behaviour this design wants.
        decided = pair in ("pass", "refuse")
        ok = (pair == case["expect"]) if decided else None
        counted = scores_for(case, model_a) and scores_for(case, model_b)
        rows.append({"case": case["id"], "expect": case["expect"],
                     "a": va, "b": vb, "pair": pair, "ok": ok,
                     "counted": counted, "findings": findings[:3]})
        out(f"  {case['id']}: expect {case['expect']} -> {model_a} {va} / "
            f"{model_b} {vb} -> {pair}"
            + ("" if ok is None else f" [{'ok' if ok else 'WRONG'}]")
            + ("" if counted else " (not scored for this pair)"))

    counted_rows = [r for r in rows if r["counted"]]
    bounds = pair_bounds(counted_rows)
    passed = pair_passes(bounds)
    rec = {"ts": time.time(), "iso": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
           "corpus_digest": digest, "corpus_cases": [c["id"] for c in cases],
           "scored_cases": [r["case"] for r in counted_rows],
           "models": sorted((model_a, model_b)),
           "families": sorted({model_family(model_a), model_family(model_b)}),
           "derived_from": derived_from, "passed": passed,
           "bounds": bounds, "calls": calls, "cases": rows}
    cache[pair_calibration_key(model_a, model_b)] = rec
    _write_calibration(cache)
    out(f"  pair {model_a} + {model_b}: {'PASSED' if passed else 'FAILED'}")
    out(f"    decided {bounds['decided']} of {bounds['n']} scored cases, "
        f"{bounds['wrong']} of them wrongly"
        + (f" ({', '.join(bounds['wrong_cases'])})" if bounds["wrong_cases"] else ""))
    out(f"    agreed: {bounds['agreed_pass']} pass, "
        f"{bounds['agreed_refuse']} refuse")
    out(f"    abstained {bounds['abstained']} "
        f"({bounds['abstained_split']} split, "
        f"{bounds['abstained_no_answer']} no answer from one model)")
    # A pair that decided nothing has no bound, and one_sided_bounds says so
    # by returning None. Measured against the real CLI over the whole roster:
    # a pair in which one model never answers abstains on every case, and
    # formatting that None as a number ends the whole run on the FIRST such
    # pair -- so the pairs after it were never scored.
    low = bounds["correct_when_agreed_at_least"]
    if low is None:
        out("    correct when they agree: NO BOUND -- the pair decided "
            "nothing, so there is no rate to bound")
    else:
        out(f"    correct when they agree, one-sided 95 % lower bound: "
            f"{low:.3f}")
    if not passed and bounds["wrong"] == 0 and bounds["decided"] == 0:
        out("    the pair decided nothing, so it has no error rate at all; "
            "that is a failure, not a clean sheet")
    return rec


def pair_calibration_record(model_a: str, model_b: str):
    """The fresh, same-corpus, passed PAIR calibration, or None."""
    import time

    try:
        rec = _read_calibration().get(pair_calibration_key(model_a, model_b))
    except Exception:
        return None
    if not isinstance(rec, dict) or not rec.get("passed"):
        return None
    if rec.get("corpus_digest") != corpus_digest():
        return None
    if time.time() - float(rec.get("ts") or 0) > CALIBRATION_MAX_AGE_S:
        return None
    return rec


def calibrate(models, out=print) -> dict:
    """Run the corpus against each model and record who passed, when, on what.

    Output: ``{model: record}``. A record is written for a failure too: which
    model failed, and on which case, is the part worth keeping. Each record
    carries the exact one-sided 95 % bounds its own rows support, because the
    pass alone says only "neither a rubber stamp nor a blanket refuser".
    """
    import time

    cases = corpus_cases()
    digest = corpus_digest()
    cache = _read_calibration()  # raises on a corrupt file rather than erasing it
    results = {}
    for model in models:
        rows = []
        calls = 0
        scored = [c for c in cases if scores_for(c, model)]
        out(f"calibrating {model} on {len(cases)} cases "
            f"({sum(c['expect'] == 'refuse' for c in cases)} positive, "
            f"{sum(c['expect'] == 'pass' for c in cases)} negative"
            + (f"; {len(cases) - len(scored)} run but not scored for this model, "
               "it reported them" if len(scored) != len(cases) else "") + ")")
        for case in cases:
            diff = _corpus_read(case["diff"]).decode()
            brief = _corpus_read(case["brief"]).decode()
            res = review_diff(model, diff, brief)
            calls += res["calls"]
            ok, why = _case_caught(case, res)
            counted = scores_for(case, model)
            rows.append({"case": case["id"], "expect": case["expect"],
                         "verdict": res["verdict"], "ok": ok, "why": why,
                         "counted": counted,
                         "findings": res["findings"][:3], "notes": res["notes"]})
            out(f"  {case['id']}: expect {case['expect']} -> {res['verdict']} "
                f"[{'ok' if ok else 'FAIL'}{'' if counted else ', not scored'}] {why}")
        counted_rows = [r for r in rows if r["counted"]]
        # A direction with no scored case left cannot be measured, and a
        # "pass" over one direction only is the rubber stamp this corpus
        # exists to tell from a reviewer.
        have_both = (any(r["expect"] == "refuse" for r in counted_rows)
                     and any(r["expect"] == "pass" for r in counted_rows))
        bounds = calibration_bounds(rows)
        rec = {"ts": time.time(), "iso": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
               "corpus_digest": digest, "corpus_cases": [c["id"] for c in cases],
               "scored_cases": [r["case"] for r in counted_rows],
               "passed": have_both and all(r["ok"] for r in counted_rows),
               "bounds": bounds, "calls": calls, "cases": rows}
        results[model] = rec
        cache[model] = rec
        out(f"  {model}: {'PASSED' if rec['passed'] else 'FAILED'} "
            f"({sum(r['ok'] for r in counted_rows)}/{len(counted_rows)} scored cases"
            + (f", {len(rows) - len(counted_rows)} not scored" if len(rows) != len(counted_rows) else "")
            + f", {calls} model calls)")
        nm = bounds["not_measured"]
        if not bounds["measurable"]:
            out(f"    NOT MEASURABLE: the endpoint answered none of "
                f"{nm['of']} scored cases, so this model has no catch rate "
                f"and no false-refuse rate. See the notes for the reason "
                f"before reading this as a judgement about the model.")
        else:
            out(f"    positives {bounds['positives']['ok']}/{bounds['positives']['n']}"
                f" -> catch rate at least {_pct(bounds['catch_rate_at_least'])}; "
                f"negatives {bounds['negatives']['ok']}/{bounds['negatives']['n']}"
                f" -> false-refuse rate at most {_pct(bounds['false_refuse_rate_at_most'])}"
                " (exact, one-sided 95 %)")
        if nm["n"]:
            out(f"    not measured: {nm['n']} of {nm['of']} scored cases "
                f"({nm['positives']} positive, {nm['negatives']} negative) "
                f"got no answer and are OUT of the rates above, while still "
                f"blocking the pass: {', '.join(nm['cases'])}")
    _write_calibration(cache)
    return results


# --- the gate itself --------------------------------------------------------


def review_gate(wt: Path, run: str, tid: str, agent: dict) -> tuple[bool, str, dict]:
    """Refuse the branch in ``wt`` unless an independent reviewer passes it.

    Output: (ok, why, record). ``record`` is the diagnosis row; ``ok`` is False
    for every reason a review could not happen, which is the point.
    """
    import time

    rec = {"ts": time.time(), "iso": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
           "run": run, "id": tid, "verdict": "", "reviewer": "", "calls": 0}
    base = str((agent.get("task") or {}).get("base") or "origin/main")
    diff = git(["diff", f"{base}...HEAD"], wt).stdout
    if not diff.strip():
        rec["verdict"] = "no-diff"
        return False, f"nothing to review: the diff against {base} is empty", rec
    if len(diff) > MAX_REVIEW_DIFF:
        rec["verdict"] = "diff-too-large"
        return False, (f"diff too large for one independent review "
                       f"({len(diff)} chars > {MAX_REVIEW_DIFF}); split the branch"), rec
    task = agent.get("task") or {}
    brief = (task.get("orig_brief") or task.get("brief") or "").strip()
    if not brief:
        rec["verdict"] = "no-brief"
        return False, "the run state holds no task brief; a reviewer with no task cannot judge", rec
    authors = author_models(agent)
    rec["authors"] = authors
    rec["author_families"] = sorted({model_family(m) for m in authors})
    cands, why = independent_candidates(authors)
    if not cands:
        rec["verdict"] = "no-reviewer"
        return False, "no independent reviewer available -- " + why, rec
    # Two reviewers of different families, not one. Measured on corpus v5:
    # no single model passes it -- kit.glm-5.3-flash false-refuses 2 of 6
    # sound branches, kit.mistral-small-4-119b-a8b misses 4 of 6 real
    # defects -- while the PAIR decided 5 of its 9 scored cases and got
    # none of them wrong, abstaining on the other 4. A gate that abstains
    # and is never wrong can be trusted with a merge decision; one that is
    # right about 70 % of the time cannot.
    #
    # Two of DIFFERENT families is preferred and tried first, but it is not
    # required, and the measurement is the reason. Requiring it closes the
    # gate for good in a case that occurs: an agent resumed on another model
    # has two author families, so of a three-family roster exactly ONE
    # independent candidate remains and no cross-family pair can form at
    # all. Meanwhile the corpus already shows that crossing families is no
    # guarantee -- kit.deepseek-v4-flash + kit.glm-5.3-flash are two
    # families and the pair FAILED, because both false-refuse the same
    # sound branch and a shared false alarm is a wrong decision, not a
    # check. So family is a prior and the calibration is the evidence: any
    # pair may stand, cross-family first, and only a pair that decided
    # nothing wrongly in both directions is allowed to decide anything.
    # The record names both families so a person reading it sees which kind
    # of pair judged the branch.
    pairs, same_family = [], []
    for i, a in enumerate(cands):
        for b in cands[i + 1:]:
            cal = pair_calibration_record(a, b)
            if not cal:
                continue
            (same_family if model_family(a) == model_family(b)
             else pairs).append((a, b, cal))
    pairs += same_family
    if not pairs:
        rec["verdict"] = "not-calibrated"
        return False, ("no independent reviewer PAIR available -- no two of "
                       + ", ".join(cands) + " have a fresh passing pair "
                       "calibration (run: delfin-agent swarm wave land --calibrate-pair A B)"
                       + ("; only one independent candidate exists, so a "
                          "pair cannot form at all" if len(cands) < 2
                          else "")), rec
    model_a, model_b, cal = pairs[0]
    rec["reviewers"] = [model_a, model_b]
    rec["calibrated_at"] = cal.get("iso")
    res = review_pair(model_a, model_b, diff, brief)
    rec["verdict"] = res["verdict"]
    rec["calls"] = res["calls"]
    rec["findings"] = res["findings"]
    rec["notes"] = (res["a"].get("notes"), res["b"].get("notes"))
    rec["pair_why"] = res["why"]
    try:
        (wt / "scratch" / "REVIEW.json").write_text(json.dumps(
            {"reviewers": [model_a, model_b], "calibrated_at": cal.get("iso"),
             "pair_bounds": cal.get("bounds"),
             "corpus": cal.get("corpus_cases"), "author_models": authors,
             "verdict": res["verdict"], "why": res["why"],
             "a": {"model": model_a, "verdict": res["a"]["verdict"],
                   "findings": res["a"]["findings"]},
             "b": {"model": model_b, "verdict": res["b"]["verdict"],
                   "findings": res["b"]["findings"]},
             "findings": res["findings"]}, indent=1) + "\n")
    except Exception:
        pass
    if res["verdict"] == "pass":
        return True, (f"independent review passed -- {model_a} and {model_b} "
                      f"both pass (pair calibrated {cal.get('iso')}, "
                      f"{cal['bounds']['wrong']} wrong of "
                      f"{cal['bounds']['decided']} decided on the corpus)"), rec
    if res["verdict"] == "refuse":
        f = res["findings"][0]
        return False, (f"independent review refused by BOTH {model_a} and "
                       f"{model_b}: {f['file']}:{f['line']} -- {f['claim']}"
                       + (f" ({f['why_it_matters']})" if f["why_it_matters"] else "")
                       + f"; {len(res['findings'])} finding(s) in "
                       "scratch/REVIEW.json"), rec
    if res["verdict"] in ("split", "no_answer"):
        # Neither a pass nor a judgement that the branch is bad. The gate
        # holds the branch and says so in those words, because "refused"
        # here would send an author to look for a defect that one of two
        # reviewers does not believe exists.
        a_f = res["a"]["findings"][:1] or res["b"]["findings"][:1]
        where = (f" (first claim: {a_f[0]['file']}:{a_f[0]['line']} -- "
                 f"{a_f[0]['claim']})") if a_f else ""
        # The REASON a model gave no answer belongs in the line, because the
        # two causes need different work: "HTTP 503 no available server" is
        # an endpoint to wait for, while a reply with no file and line is a
        # model that cannot be used for review at all. Dropping the notes
        # here turned both into the same sentence.
        notes = []
        for side, who in ((res["a"], model_a), (res["b"], model_b)):
            if side["verdict"] == "unusable" and side.get("notes"):
                notes.append(f"{who}: " + "; ".join(side["notes"])[:200])
        return False, ("independent review UNDECIDED, a person decides: "
                       + res["why"] + where
                       + ("; " + " | ".join(notes) if notes else "")
                       + "; both reviews in scratch/REVIEW.json"), rec
    # review_pair returns one of four verdicts and all four are handled
    # above. This line is the one a new verdict would reach, and it says so
    # rather than naming a variable that the pair rewrite removed.
    return False, (f"the reviewer pair returned an unhandled verdict "
                   f"{res['verdict']!r} ({model_a} + {model_b}); left "
                   "unopened for a person"), rec


def review_summary(records: list[dict]) -> dict:
    """The run's review diagnosis.

    Every rate carries its absolute count beside it (a rate alone is a target
    in disguise). The number of reviews performed is NOT a quality measure --
    it rises with every extra reviewer and says nothing about the code -- so
    nothing here is to be optimised; it is read, by a person, next to the
    refusals themselves.
    """
    reviewed = [r for r in records if r.get("reviewer")]
    refused = [r for r in reviewed if r["verdict"] == "refuse"]
    unusable = [r for r in reviewed if r["verdict"] == "unusable"]
    blocked = [r for r in records if not r.get("reviewer")]
    return {
        "branches_seen": len(records),
        "reviewed": len(reviewed),
        "refused": len(refused),
        "passed": len([r for r in reviewed if r["verdict"] == "pass"]),
        "unusable": len(unusable),
        "refusal_rate": (len(refused) / len(reviewed)) if reviewed else None,
        "refusal_rate_of": len(reviewed),
        "no_review_possible": len(blocked),
        "reviewers": sorted({r["reviewer"] for r in reviewed}),
        "verdicts": [{"run": r.get("run"), "id": r.get("id"),
                      "reviewer": r.get("reviewer"), "verdict": r.get("verdict")}
                     for r in records],
    }


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    dry_run = "--push" not in argv
    argv = [a for a in argv if a not in ("--dry-run", "--push")]
    if "--calibrate-pair" in argv:
        # Exact list membership, so this never swallows "--calibrate".
        i = argv.index("--calibrate-pair")
        names = [a for a in argv[i + 1:] if not a.startswith("-")]
        reuse = "--no-reuse" not in argv
        if len(names) == 1:
            print("--calibrate-pair takes two model ids, or none for every "
                  "cross-family pair that already has stored rows", flush=True)
            return 2
        if len(names) >= 2:
            todo = [(names[0], names[1])]
        else:
            # Default: every cross-family pair whose two models were already
            # measured on THIS corpus. Costs no model calls, so the honest
            # default is to score them all rather than to pick one.
            cache = _read_calibration()
            digest = corpus_digest()
            have = sorted(m for m, r in cache.items()
                          if not m.startswith("pair:")
                          and isinstance(r, dict)
                          and r.get("corpus_digest") == digest
                          and r.get("cases"))
            # Every pair, same-family ones included: scoring them from
            # stored rows costs nothing, and a same-family pair is the only
            # one available when an agent's authors span two families. The
            # corpus decides which may judge; this only measures them.
            todo = [(a, b) for n, a in enumerate(have) for b in have[n + 1:]]
            if not todo:
                print("fewer than two models have rows on the current corpus "
                      f"({digest[:12]}); run --calibrate first", flush=True)
                return 1
            cross = sum(model_family(a) != model_family(b) for a, b in todo)
            print(f"scoring {len(todo)} pair(s) from stored rows "
                  f"({cross} cross-family, {len(todo) - cross} same-family)",
                  flush=True)
        passed = []
        for a, b in todo:
            rec = calibrate_pair(a, b, reuse=reuse)
            if rec["passed"]:
                passed.append((a, b, rec))
        print()
        if not passed:
            print("NO pair passed. The gate stays closed, which is the "
                  "designed behaviour: a swarm branch waits for a person "
                  "rather than being merged on an uncalibrated opinion.",
                  flush=True)
            return 1
        print(f"{len(passed)} pair(s) passed:", flush=True)
        for a, b, rec in passed:
            bd = rec["bounds"]
            low = bd["correct_when_agreed_at_least"]
            print(f"  {a} + {b}: decided {bd['decided']}/{bd['n']}, "
                  f"{bd['wrong']} wrong, abstained {bd['abstained']}, "
                  + ("correct-when-agreed >= "
                     f"{low:.3f}" if low is not None else "no bound"),
                  flush=True)
        print("\nRead the abstention count as the cost: those branches still "
              "need a person. It is the one number here the pair itself "
              "controls, so it is reported and never optimised.", flush=True)
        return 0
    if "--calibrate" in argv:
        i = argv.index("--calibrate")
        models = [a for a in argv[i + 1:] if not a.startswith("-")]
        argv = argv[:i]
        models = models or list(kit_reviewer_models())
        if not models:
            print("no KIT reviewer models to calibrate (roster unreadable)", flush=True)
            return 1
        res = calibrate(models)
        passed = [m for m, r in res.items() if r["passed"]]
        # The bounds, not the pass, are the result: a pass on a small corpus
        # says only "neither a rubber stamp nor a blanket refuser".
        for m, r in res.items():
            b = r["bounds"]
            nm = b.get("not_measured") or {"n": 0, "of": 0}
            head = f"{m}: {'PASSED' if r['passed'] else 'FAILED'}; "
            if not b.get("measurable", True):
                print(head + f"NOT MEASURABLE -- no answer on any of "
                      f"{nm['of']} scored cases, so no rate is reported",
                      flush=True)
            else:
                print(head
                      + f"catch rate >= {_pct(b['catch_rate_at_least'])} "
                      f"({b['positives']['ok']}/{b['positives']['n']} positives), "
                      f"false-refuse rate <= {_pct(b['false_refuse_rate_at_most'])} "
                      f"({b['negatives']['ok']}/{b['negatives']['n']} negatives)"
                      + (f", {nm['n']} of {nm['of']} not measured"
                         if nm["n"] else ""), flush=True)
        print("calibration passed: " + (", ".join(passed) or "none"), flush=True)
        return 0 if passed else 1
    runs = argv or None
    results = land(runs, dry_run=dry_run, only={t for t in os.environ.get("DELFIN_LAND_ONLY", "").split(",") if t})
    return 0 if all(r["outcome"] != "error" for r in results) else 1


def _push(repo: Path, local: str, branch: str):
    """The push itself, run directly in the run's repository: the person's own
    checkout, whose configuration, hooks and info/ the members' cage binds
    read-only, so nothing in it was written by a member; and the push needs
    the person's credential helper, which ``safe_git`` switches off."""
    return sh(["git", "push", "-q", "origin", f"{local}:refs/heads/{branch}"], repo)


def _forge_repo(repo: Path) -> str:
    """``owner/name`` of the repository a PR is opened in: the setting, else
    what the forge's CLI reports for the run's repository."""
    if _swarm_settings().get("github_repo"):
        return str(_swarm_settings()["github_repo"])
    r = sh(["gh", "repo", "view", "--json", "nameWithOwner", "-q", ".nameWithOwner"], repo)
    return r.stdout.strip()


def _record(run: str, kind: str, tid: str, **data) -> None:
    """One event in the run's log, where ``wave.watch`` reads it."""
    import time
    try:
        with open(ROOT / run / "events.jsonl", "a") as fh:
            fh.write(json.dumps({"t": round(time.time(), 1), "run": run, "kind": kind,
                                 "id": tid, **data}) + "\n")
    except OSError:
        pass


def land(runs=None, *, dry_run: bool = True, only=None, log=print) -> list[dict]:
    """Check, review and open every member branch that wrote scratch/PR.md.

    ``runs`` default every run under the waves root; ``only`` limits the
    members; ``dry_run`` (the default) checks and reviews but pushes and opens
    nothing. Returns one record per member considered:
    ``{"run", "id", "outcome": "refused"|"ready"|"opened"|"error", "why"}``."""
    records: list[dict] = []
    results: list[dict] = []
    for run in (runs or sorted(p.name for p in ROOT.iterdir() if (p / "state.json").is_file())):
        st = json.loads((ROOT / run / "state.json").read_text())
        mode = str(st.get("land") or "pr")
        if mode == "none":
            log(f"{run}: skipped -- the run says land: none")
            continue
        repo = Path(st.get("repo") or ".")
        for tid, a in sorted(st["agents"].items()):
            if only and tid not in only:
                continue
            task = a.get("task") or {}
            if str(task.get("land") or "") == "none":
                continue
            wt = Path(a["worktree"])
            pr = wt / "scratch" / "PR.md"
            if not pr.is_file():
                continue
            note = wt / "scratch" / "LANDING.txt"
            base = str(task.get("base") or "origin/main")

            def done(outcome, why):
                note.write_text(f"{outcome}: {why}\n")
                results.append({"run": run, "id": tid, "outcome": outcome, "why": why})
                _record(run, "land_" + outcome, tid, why=why[:300])
                log(f"{run}/{tid}: {outcome} -- {why}")

            ok, why = check(wt, base, run=run, tid=tid)
            if not ok:
                done("refused", why)
                continue
            text = pr.read_text().strip().splitlines()
            title, body = text[0].lstrip("# ").strip(), "\n".join(text[1:]).strip()
            # The PR text is public too: a home path or user name in it is refused.
            if HOME_PATH.search(title + "\n" + body) or HOME.name in (title + body):
                done("refused", "PR.md carries a home path or user name; rewrite it")
                continue
            if re.search(r"\bWIP\b|budget ran out", title, re.I):
                done("refused", "PR.md title says the work is unfinished")
                continue
            # The message a squash merge publishes ON THE BASE. Refused here,
            # before the push, because the base's history is permanent: a model,
            # a provider or a site name in it cannot be taken back, and a second
            # Co-Authored-By line makes a contributor of a branch that had one.
            squash, sq_why = squash_message(wt, base, title)
            if sq_why:
                done("refused", sq_why)
                continue
            # Nothing above this line was checked by anyone but the author.
            # The independent review runs here, before the push, and a refusal
            # leaves the branch unopened exactly like the refusals above it.
            r_ok, r_why, rec = review_gate(wt, run, tid, a)
            log(f"{run}/{tid}: review {'passed' if r_ok else 'REFUSED'} -- {r_why}")
            records.append(rec)
            try:
                with REVIEW_LOG.open("a") as fh:
                    fh.write(json.dumps(rec) + "\n")
            except Exception:
                pass
            if not r_ok:
                done("refused", r_why)
                continue
            branch = f"wave/{run}-{tid}"
            local = git(["rev-parse", "--abbrev-ref", "HEAD"], wt).stdout.strip()
            if dry_run:
                log(f"  dry-run: would push {local} -> {branch} and open "
                    f"\"{title}\" against the base ({len(body)} chars of body)")
                log(f"  squash message ({len(squash)} chars): "
                    + " / ".join(squash.splitlines()))
                done("ready", f"would push {local} -> {branch}; review: {r_why}")
                continue
            push = _push(repo, local, branch)
            if push.returncode != 0:
                done("error", f"push failed: {push.stderr.strip()[:200]}")
                continue
            bodyf = wt / "scratch" / "PR.body.md"
            bodyf.write_text(body + "\n")
            target = _forge_repo(repo)
            cmd = ["gh", "pr", "create", "--repo", target,
                   "--base", str(_swarm_settings().get("land_base") or "main"),
                   "--head", branch, "--title", title, "--body-file", str(bodyf)]
            if mode == "draft":
                cmd.append("--draft")
            r = sh(cmd, repo)
            if r.returncode == 0:
                pr.rename(pr.with_suffix(".md.opened"))
                done("opened", r.stdout.strip())
            else:
                done("error", (r.stderr or r.stdout).strip()[:300])
    if records:
        s = review_summary(records)
        if s["reviewed"]:
            log("review diagnosis (not a target): "
                f"{s['reviewed']} reviewed, {s['refused']} refused, "
                f"{s['passed']} passed, {s['unusable']} unusable, "
                f"{s['no_review_possible']} could not be reviewed; refusal rate "
                f"{s['refusal_rate']:.2f} of {s['refusal_rate_of']} reviews")
        else:
            log("review diagnosis (not a target): no branch reached a reviewer "
                f"({s['no_review_possible']} could not be reviewed)")
    return results


if __name__ == "__main__":
    sys.exit(main())
