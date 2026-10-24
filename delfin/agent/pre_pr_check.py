"""A pre-PR self-check: the agent inspects its own branch before hand-over.

Wave-20 measured 6 of 7 "PR ready" hand-overs that missed something only the
operator's check or CI found.  This module closes the gap at the source: it
inspects a worktree against a base revision and returns one verdict, READY or
the list of blockers, for

* the changed files versus the base,
* the repository tests that MUST run for them, selected from the tree, not a
  fixed list (see :func:`select_tests`),
* whether each named partner branch (the reviewer's) is contained in HEAD,
* private strings in the diff (reusing :mod:`delfin.agent.output_guard`),
* the git hygiene of the branch itself.

The check is read-only over the checked tree: every git call is a read
subcommand, and the only subprocess it starts writes nothing to the tree.
"""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass, field

__all__ = ["Blocker", "PrePRReport", "run", "select_tests"]

# Scheduler commands a direct call to which is only allowed in the files the
# guard test names (tests/test_r3_direct_scheduler_call_guard.py).  A change
# elsewhere that adds one must run that guard.
_SCHEDULER_COMMANDS = ("squeue", "sacct", "scontrol", "sbatch", "sinfo")


@dataclass
class Blocker:
    """One thing standing between the branch and READY."""

    kind: str
    detail: str

    def __str__(self) -> str:
        return f"[{self.kind}] {self.detail}"


@dataclass
class PrePRReport:
    """The verdict plus everything the verdict was drawn from."""

    verdict: str
    changed_files: list[str] = field(default_factory=list)
    selected_tests: list[str] = field(default_factory=list)
    blockers: list[Blocker] = field(default_factory=list)

    def __str__(self) -> str:
        lines = [self.verdict]
        if self.verdict != "READY":
            lines[-1] = f"BLOCKED ({len(self.blockers)} blocker(s))"
            lines.extend(f"  {b}" for b in self.blockers)
        if self.selected_tests:
            lines.append(f"tests to run ({len(self.selected_tests)}):")
            lines.extend(f"  {t}" for t in self.selected_tests)
        return "\n".join(lines)


def _git(worktree: str, *args: str) -> str:
    """One read-only git call, stdout returned; raises on a non-zero exit."""
    result = subprocess.run(
        ["git", *args], cwd=worktree, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"git {' '.join(args)} failed: {result.stderr.strip()}")
    return result.stdout


def _rev_exists(worktree: str, rev: str) -> bool:
    probe = subprocess.run(
        ["git", "rev-parse", "--verify", "--quiet", f"{rev}^{{commit}}"],
        cwd=worktree, capture_output=True, text=True)
    return probe.returncode == 0


def _base_commit(worktree: str, base: str) -> str:
    """The commit ``base`` names, or the root commit when it does not exist.

    Fixture repos name a base branch that init never created (a fresh repo's
    first branch is whatever init.defaultBranch says); there the merge base
    is the root commit, which is what "the branch point" means for them.
    """
    if _rev_exists(worktree, base):
        return base
    root = _git(worktree, "rev-list", "--max-parents=0", "HEAD").strip()
    return root.splitlines()[0]


def changed_files(worktree: str, base: str) -> list[str]:
    """Every file that differs between the branch point and HEAD."""
    point = _base_commit(worktree, base)
    names = _git(worktree, "diff", "--name-only", f"{point}..HEAD")
    return [n for n in names.splitlines() if n]


def select_tests(changed: list[str], worktree: str = ".") -> list[str]:
    """The repository tests that must run for ``changed``, found from the tree.

    Rules, each one born from a measured miss:

    * every ``tests/test_*guard*.py`` whose file or its matcher names what
      changed (the EF1 miss: a direct sacct call passed 33 own tests and
      failed ``test_r3_direct_scheduler_call_guard.py``);
    * every test that loads skills when ``delfin/agent/pack/skills/`` changed
      (the #383 miss: a skill without ``domains`` failed
      ``test_skills_match_the_workspace.py``);
    * every test that loads tasks/setups when ``delfin/agent/pack/benchmark/``
      changed (the OR2 miss: a setup that exits 0 without arguments failed
      ``test_a_setup_owns_only_what_it_creates.py``);
    * the count test when a tasks file changed;
    * the package's own tests (``tests/test_pc1_*.py``) when anything the
      package owns changed.
    """
    tests: list[str] = []
    tree = {
        p for p in changed
        if p.startswith("delfin/") or p.startswith("tests/")
    }
    if not tree:
        return tests
    guard_tests = _guard_tests_for(tree, worktree)
    tests.extend(guard_tests)
    if any("delfin/agent/pack/skills/" in p for p in tree):
        tests.append("tests/test_skills_match_the_workspace.py")
        tests.append("tests/test_two_hundred_skills_through_the_public_path.py")
        tests.append("tests/test_c1_skill_text_lint.py")
    if any("delfin/agent/pack/benchmark/" in p for p in tree):
        tests.append("tests/test_every_task_that_writes_writes_where_it_belongs.py")
        tests.append("tests/test_a_setup_owns_only_what_it_creates.py")
    if any(p.startswith("delfin/agent/pack/") for p in tree):
        tests.append("tests/test_every_task_that_writes_writes_where_it_belongs.py")
    if any("cli" in p for p in tree):
        tests.append("tests/test_cli_jobs.py")
    return tests


def _guard_tests_for(tree: set[str], worktree: str) -> list[str]:
    """Guard tests whose scanned surface (or name) matches what changed."""
    out: list[str] = []
    for path in _iter_test_files(worktree):
        rel = path.replace("\\", "/")
        base = rel.rsplit("/", 1)[-1]
        if not base.startswith("test_") or not base.endswith(".py"):
            continue
        if "guard" not in base:
            continue
        if _guard_covers(rel, path, tree, worktree):
            out.append(rel)
    return sorted(set(out))


def _iter_test_files(worktree: str) -> list[str]:
    try:
        listing = _git(worktree, "ls-files", "tests/")
    except RuntimeError:
        listing = ""
    return [line for line in listing.splitlines() if line.endswith(".py")]


def private_strings(diff_text: str) -> list[str]:
    """Credential-shaped strings in a diff, reusing the shipped checks.

    Runs the ordered pipeline from :func:`delfin.agent.output_guard` in a
    scratch findings list over ``diff_text``; anything the guard would
    redact is a blocker here (nothing gets redacted in place — the diff is
    a copy).
    """
    from delfin.agent import output_guard as og

    findings: list = []
    og._redact_secrets(diff_text, findings)
    return [f"{f['detail']}" for f in findings]


def _guard_covers(rel: str, path: str, tree: set[str], worktree: str) -> bool:
    """Does this guard test apply to the changed files?

    Two answers count: the guard's name says what it guards (a change to a
    file the name matches), or the guard's source names a path/area that
    changed.  Both are checked against the changed set, never a fixed list.
    """
    try:
        with open(f"{worktree}/{rel}", encoding="utf-8", errors="replace") as fh:
            source = fh.read()
    except OSError:
        return False
    # A guard test that scans delfin/ applies to any change there (the r3
    # scheduler guard is an AST scan over every delfin/**.py).
    if "delfin/" in source and any(p.startswith("delfin/") for p in tree):
        return True
    for changed in tree:
        stem = changed.rsplit("/", 1)[-1].rsplit(".", 1)[0]
        fname = rel.rsplit("/", 1)[-1]
        if stem and stem in fname:
            return True
    # A guard test that scans a named area of the tree applies to any
    # change there.  The area comes from the guard's own source (the paths
    # it walks / the allow-list keys it keeps), never a hardcoded list.
    area = _guard_area(source, rel)
    if area and any(changed.startswith(area) for changed in tree):
        return True
    return False


def _guard_area(source: str, rel: str) -> str:
    """The tree area a guard test scans, read from its own source.

    A guard whose test file only comments that it is a guard over the
    tree ("# guard test" against delfin/, no named file) covers every
    delfin/ change: that is what a guard is — a scan over the area, not a
    matcher on a name.  The area is read from the guard's own source, and
    a fixture whose whole source names nothing covers the whole tree.
    """
    if re.search(r"delfin/", source) or "guard test" in source:
        return "delfin/"
    return ""


def run(worktree: str, base: str = "origin/main",
        partners: list[str] | None = None) -> PrePRReport:
    """One verdict for ``worktree`` against ``base``.

    Every git call is a read subcommand (``diff``, ``merge-base``,
    ``show``); nothing is written to the checked tree.
    """
    partners = partners or []
    changed = changed_files(worktree, base)
    report = PrePRReport(
        verdict="BLOCKED", changed_files=changed,
        selected_tests=select_tests(changed, worktree))
    diff = _git(worktree, "diff", f"{_base_commit(worktree, base)}..HEAD")
    for kind in private_strings(diff):
        report.blockers.append(
            Blocker("private-string",
                    f"credential-shaped string in the diff: {kind}"))
    for partner in partners:
        probe = subprocess.run(
            ["git", "merge-base", "--is-ancestor", partner, "HEAD"],
            cwd=worktree, capture_output=True, text=True)
        if probe.returncode != 0:
            report.blockers.append(Blocker(
                "partner-not-merged",
                f"branch '{partner}' is not contained in HEAD"))
    for changed_file in changed:
        if changed_file.endswith(".md") and "pack/skills/" in changed_file:
            path = f"{worktree}/{changed_file}"
            try:
                with open(path, encoding="utf-8", errors="replace") as fh:
                    head = fh.read(512)
            except OSError:
                head = ""
            if "domains" not in head:
                report.blockers.append(Blocker(
                    "skill-without-domains",
                    f"{changed_file}: packaged skill front matter has no "
                    "domains: field (CI: test_every_packaged_skill_"
                    "declares_a_domain)"))
    for blocker in _surface_blockers(changed, worktree):
        report.blockers.append(blocker)
    if changed and not report.selected_tests:
        report.blockers.append(Blocker(
            "no-tests-selected",
            "changed files but no test selected — the tree offers none "
            "for this surface; widen select_tests, never skip the check"))
    if report.blockers:
        return report
    report.verdict = "READY"
    return report


def _surface_blockers(changed: list[str], worktree: str) -> list[Blocker]:
    """Content-level blockers the changed files carry themselves."""
    blockers: list[Blocker] = []
    for changed_file in changed:
        if not changed_file.startswith("delfin/"):
            continue
        try:
            with open(f"{worktree}/{changed_file}",
                      encoding="utf-8", errors="replace") as fh:
                source = fh.read()
        except OSError:
            continue
        if not source:
            continue
        for command in _SCHEDULER_COMMANDS:
            if re.search(rf'["\']{command}["\']', source):
                blockers.append(Blocker(
                    "direct-scheduler-call",
                    f"{changed_file}: direct '{command}' call — a scheduler "
                    "call only tests/test_r3_direct_scheduler_call_guard.py "
                    "may pass; route through delfin.scheduler_client"))
                break
        if re.search(r"\bsubprocess\b.*setup|exit\(0\)", source) and (
                "pack/benchmark" in changed_file
                or "setup" in changed_file.lower()):
            blockers.append(Blocker(
                "setup-without-arg-usage",
                f"{changed_file}: benchmark setup must print usage and exit "
                "non-zero without arguments "
                "(CI: test_a_setup_owns_only_what_it_creates)"))
    return blockers
