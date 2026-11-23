"""A git read does not ask the kernel for one thread per CPU.

``git grep`` runs one thread per core and, when a thread cannot be created,
exits 128 with "fatal: grep: failed to create thread: Resource temporarily
unavailable" -- no match reported, no "no match" either.

Measured 2026-11-20 on a workstation with 384 cores and a load average of
192, inside an agent session (whose process and thread budget is already
partly spent): ``git grep -q -E '(def|class) present\\b'`` in a two-file
repository exited 128, while the same command with ``grep.threads=1``
exited 0. Every ``git grep`` behind ``safe_git`` failed there.

The product was already written for the failure -- ``knowledge_curation.
_git_grep`` maps a non-0/1 exit to ``None`` and does not report the
reference as broken, because a failure "cannot tell". That is why this was
silent rather than loud: the knowledge stats simply left broken function
references out. ``tests/test_knowledge_curation`` is the one test that asks
for them, and it failed here and passed on CI on the same commit.

The fix bounds the thread pool in ``safe_git``, so every git read of this
product creates a fixed number of threads whatever the host reports. The
tests below are what keep that measured:
- ``argv_for`` puts the bound on ``grep`` and on nothing else;
- the bound is small enough to matter and not zero;
- a ``git grep`` through ``run`` actually finds a match on this host, which
  is the end-to-end check and needs no git-config knowledge to read.

The alternative rejected: ``grep.threads=1`` alone. One thread serialises
the search over a large tree; this reads a repository of ~5,600 files on
every knowledge command, and the failure was at 384 threads, not at four.
"""

from __future__ import annotations

import subprocess

from delfin.agent import safe_git


def _grep_argv(cwd):
    return safe_git.argv_for(["grep", "-q", "-E", "x"], str(cwd))


def test_the_grep_option_is_bounded_and_not_serial(tmp_path):
    """A whole-CPU thread pool is what failed; one thread is needlessly slow.

    The bound is asserted as a property (a positive integer well below the
    core count of the largest machines this runs on), not as one number, so
    tuning it does not need a test edit -- but removing it or setting it to
    0/1 does.
    """
    argv = _grep_argv(tmp_path)
    settings = [a for a in argv if a.startswith("grep.threads=")]
    assert settings, f"no grep.threads bound in {argv}"
    threads = int(settings[0].split("=", 1)[1])
    assert 1 < threads <= 64, threads


def test_only_grep_carries_the_bound(tmp_path):
    """``-c`` options go before the subcommand; other reads are unaffected.

    A ``git rev-parse`` picking up a grep option would be noise in every
    argv this product builds, and a reader of ``argv_for`` would have to
    know why.
    """
    for sub in ("rev-parse", "ls-files", "status"):
        argv = safe_git.argv_for([sub], str(tmp_path))
        assert not [a for a in argv if a.startswith("grep.threads=")], (sub, argv)


def test_a_git_grep_through_run_finds_a_match_on_this_host(tmp_path):
    """The end-to-end check: the read this product makes must work here.

    Before the fix this exited 128 on this host with "failed to create
    thread", which ``run`` reports as a non-zero code -- and which
    ``knowledge_curation`` reads as "cannot tell". The test needs no
    knowledge of git's configuration to read: a tracked file holds a word,
    and ``git grep -q`` must find it.
    """
    repo = tmp_path / "repo"
    repo.mkdir()
    subprocess.run(["git", "init", "-q", str(repo)], check=True,
                   capture_output=True)
    (repo / "a.py").write_text("def a_present_helper():\n    pass\n",
                               encoding="utf-8")
    subprocess.run(["git", "-C", str(repo), "add", "."], check=True,
                   capture_output=True)
    found = safe_git.run(["grep", "-q", "-E", r"(def|class) a_present_helper\b"],
                         repo, timeout=30, text=False)
    assert found.returncode == 0, (
        f"git grep exited {found.returncode}: {found.stderr[:200]!r}")
    absent = safe_git.run(["grep", "-q", "-E", r"(def|class) vanished_helper\b"],
                          repo, timeout=30, text=False)
    assert absent.returncode == 1, (
        f"a no-match must exit 1, got {absent.returncode}: {absent.stderr[:200]!r}")
