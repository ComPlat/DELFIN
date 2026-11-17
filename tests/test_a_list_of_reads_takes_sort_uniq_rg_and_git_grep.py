"""sort, uniq, rg and git grep join the routine list of reads -- with the
options each lists, and no other.

A list of plain reads (``grep x | wc -l``) is routine and costs no model
call. ``grep x | sort | uniq -c`` and ``git grep -n x | wc -l`` still went
to the model, which refused such reads as "exploratory" (Pi, 24 h: 22
refusals of commands the permission rules allowed). The four heads have
options that run a program or write a file with no shell operator:
``git grep -O`` opens every match in a "pager" (runs it), ``uniq IN OUT``
writes OUT, ``sort -o`` (also inside ``-uo``) writes, ``sort
--compress-program`` runs a program, ``rg -z`` runs decompressors and
``rg --pre`` a preprocessor. So each head takes only the options it lists;
any other option leaves the command with the model.
"""

import pytest

from delfin.agent import action_classifier as ac


def _verdict(command):
    return ac.screen({"tool": "bash", "command": command},
                     {"task": "read the code"}, ())


def _routine(command):
    v = _verdict(command)
    return v is not None and v.decision == "allow" and "routine" in v.reason


@pytest.mark.parametrize("command", [
    "grep -rn def delfin/agent | sort | uniq -c",
    "grep -c x delfin | sort -t: -k2 -nr | head",
    "git grep -n classifier | wc -l",
    "git grep -n -i -e foo -- delfin",
    "git grep -C2 -w foo | head",
    "rg --line-number foo delfin | wc -l",
    "rg -n -g '*.py' -t py foo delfin",
    "sort a.txt | uniq -f 1 -c",
    "uniq -c a.txt",
    "ls delfin | sort -r",
])
def test_the_listed_heads_are_routine(command):
    assert _routine(command), command


@pytest.mark.parametrize("command", [
    # run a program
    "git grep -Osh foo",
    "git grep -nOsh foo",
    "git grep --open-files-in-pager=sh foo",
    "sort --compress-program=sh a.txt",
    "rg -z foo delfin",
    "rg --pre ./x.sh foo delfin",
    "rg --pre=./x.sh foo delfin",
    "grep x a | rg --search-zip y",
    # write a file
    "uniq a.txt delfin/agent/api_client.py",
    "grep x a | uniq - out.txt",
    "sort -uo out.txt a.txt",
    "sort -o out.txt a.txt",
    "sort --output=out.txt a.txt",
    "sort -T tmpdir a.txt",
    # read beyond the operands
    "git grep -f pats.txt",
    "git grep --no-index foo /etc",
    "sort --files0-from=list a.txt",
    "rg -L foo delfin",
    "rg --ignore-file=/etc/x foo delfin",
    "rg --config extra.conf foo delfin",
    # the rest of the list must still be reads
    "grep x a | sort | sh",
    "grep x a | sort > out.txt",
    "sort a.txt | tee out.txt",
])
def test_options_that_run_or_write_stay_with_the_model(command):
    assert not _routine(command), command
