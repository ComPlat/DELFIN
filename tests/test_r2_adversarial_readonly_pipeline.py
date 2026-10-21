"""Adversarial control for the readonly_pipeline git branch/tag hole.

A ``git branch``/``git tag`` line is read-only ONLY in list form: the same
subcommand carries every write (``-d`` delete, ``-m`` move, ``-a`` annotate,
a bare name creates). Since a ``is_safe()==True`` skips the confirm dialog,
True here must never carry a write ↔ history change.

Second block: the genuinely read-only counterparts stay safe (list form, or
read-only subcommands).
"""
import pytest

from delfin.agent.readonly_pipeline import is_safe


@pytest.mark.parametrize("cmd", [
    # a bare name after branch/tag creates one
    "git branch -d feature-x",        # deletes a branch
    "git branch -D feature-x",        # force-delete
    "git branch newfeature",          # creates a branch
    "git branch -m old new",          # renames
    "git branch --move old new",
    "git tag -d v1.0",                # deletes a tag
    "git tag v1.0",                   # creates a tag
    "git tag -a v1.0 -m msg",         # annotated tag
    "git tag -d v1; git log",         # write inside a sequence
    # a write carried by an intrinsic command's flag, never read-only
    "sort -o out.txt data.txt",       # -o writes a file
    "sort --output out.txt data.txt",
    "sort --output=out.txt data.txt",
    "date -s '2030-01-01'",           # -s sets the clock
    "date --set='2030-01-01'",
])
def test_git_write_flag_forms_are_never_read_only(cmd):
    assert is_safe(cmd) is False


@pytest.mark.parametrize("cmd", [
    "git branch",                     # bare list
    "git tag",                        # bare list
    "git branch --list",
    "git branch -v",
    "git tag --list",
    "git log",
    "git show HEAD:path",
    "git status",
    "git diff",
    "git rev-parse --is-inside-work-tree",
    # read-only counterparts of the flag-write forms
    "sort data.txt",                  # no -o: reads, prints to stdout
    "sort -k2 data.txt | uniq -c",
    "date",                           # read-form: prints the date
    "date +'%Y-%m-%d'",
])
def test_git_read_only_forms_stay_safe(cmd):
    assert is_safe(cmd) is True
