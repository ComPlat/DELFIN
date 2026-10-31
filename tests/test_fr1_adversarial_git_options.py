# FR1 reviewer adversarial tests (nacht-s17): the git plain-read parser
# around _GIT_WRITING_OPTION_PREFIXES. Attacks: attached short values,
# unicode look-alikes, value-carrying variants, chained commands, and the
# ordinary harmless reads that must stay plain reads (no over-refusal).
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from delfin.agent.action_classifier import is_plain_read  # noqa: E402


def _read(cmd):
    return is_plain_read({"tool": "bash", "command": cmd})


def test_git_short_output_attached_and_separated():
    # -oFILE attached, -o=FILE, and -o with the value later: all write.
    for cmd in ("git diff -oout.patch main",
                "git show -o=out.patch HEAD",
                "git log -o out.patch -1"):
        assert not _read(cmd), cmd


def test_git_short_output_inside_quoting_is_not_punctuation_slip():
    # The value hides the file path, not the intent: still not a plain read.
    assert not _read("git show -o 'out patch' HEAD")


def test_git_short_output_lookalike_unicode_minus():
    # A unicode minus is not an option dash for git: `git show ‑o …` fails
    # as an unknown rev, it cannot write. The classifier may count it as a
    # plain read -- the harmless command it is.
    assert _read("git show ‑o out.patch HEAD")


def test_ordinary_reads_stay_plain():
    # No over-refusal: the fix must not swallow harmless reads.
    for cmd in ("git show HEAD --stat",
                "git log --oneline -5",
                "git diff main~1 main",
                "git branch --list 'agent/*'",
                "git tag --contains HEAD",
                "git log --format=%h -1"):
        assert _read(cmd), cmd


def test_chained_write_after_read_is_not_plain():
    for cmd in ("git log -1 && git push origin HEAD",
                "git show HEAD ; rm file",
                "git diff > out.diff"):
        assert not _read(cmd), cmd
