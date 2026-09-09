"""`python -c` is the one opaque command whose program is in the text.

An interpreter is kept off the bash auto-allow list because the checks
that decide on a command string cannot decide about it: `xargs` reads its
targets from a file, `make` runs whatever the Makefile says. For those
the objection is exact. For `python -c` it was not -- the program is a
literal, right there in the argument -- and the cost of treating it as
unreadable was measured on 2026-09-09 over ~/.delfin/audit.log:

    199 of 1381 recorded bash calls (14.4%) refused as "not on the
    auto-allow list", touching 135 of 431 sessions. 130 of the 199 carry
    an inline payload; 103 of those only read, print and compute.

The commonest single shape is the agent checking its own work
(`ast.parse(open('x.py').read())`), and the second is a scientist's
one-liner over a CSV or a workbook. In a headless run the auto-allow
list IS the permission surface, so each one is a dead end, and the block
message correctly tells the agent not to look for a way around it.

So the payload is read instead (delfin/agent/inline_payload.py) and the
two checks the ban stood in for run for real:

  * the CONTENT scan, which every executed script file already gets and
    which an inline payload escaped only because there was no file to
    open -- including any module the payload imports out of the working
    directory, since `import mymod` executes mymod.py exactly as running
    the script would;
  * the WRITE gate, which reads paths out of a command and knew cp, mv,
    tee, sed -i and dd but nothing about `open(path, 'w')`.

This file pins the analyser's own contract. Every rule gets a case it
must match and a case it must not, because the failure mode of a
static reader is a construct it half-understands: `s.replace(',', '.')`
-- how every German-written amount is parsed, and the commonest line in
the whole corpus -- read as `Path.replace(target)` and reported a write
to the file named ".".
"""

from __future__ import annotations

import pytest

from delfin.agent.inline_payload import (
    Effects, analyze_payload, extract_c_payloads,
)


def _eff(src: str, cwd=None) -> Effects:
    return analyze_payload(src, cwd)


# ---------------------------------------------------------------------------
# Readable: reads, prints, computes
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("src", [
    "print(1)",
    "import ast; ast.parse(open('t.py').read()); print('ok')",
    "import json; d = json.load(open('r.json')); print(d['n'])",
    "import csv\nrows = list(csv.DictReader(open('b.csv')))\nprint(len(rows))",
    "import openpyxl; wb = openpyxl.load_workbook('x.xlsx'); print(wb.sheetnames)",
    "import statistics as s; print(s.mean([1, 2, 3]))",
    "import sys; print(sys.version_info[:2])",
    "import sys; sys.stdout.write('hi\\n')",
    "print(float('1.265,85'.replace('.', '').replace(',', '.')))",
    "import pandas as pd; print(pd.read_csv('a.csv').shape)",
    "d = {'a': 1}; e = d.copy(); print(e)",
    "with open('in.txt') as f: print(f.read()[:10])",
])
def test_a_payload_that_only_reads_is_readable(src):
    eff = _eff(src)
    assert eff.readable, eff.opaque
    assert eff.writes == [], eff.writes


# ---------------------------------------------------------------------------
# Opaque: the analyser says so rather than guessing
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("src", [
    "import os; os.remove('x')",
    "import shutil; shutil.rmtree('/home/u')",
    "import subprocess; subprocess.run(['sh'])",
    "exec(open('p').read())",
    "eval('1+1')",
    "__import__('os').system('id')",
    "import socket; socket.socket()",
    "import smtplib; smtplib.SMTP('h').send_message(m)",
    "from pathlib import Path; Path('x').write_text('y')",
    "import importlib; importlib.import_module('os')",
    "f = open('o.txt', 'w'); f.write('x'); g = open(p, 'a')",   # computed path
    "open('o.txt', mode); print(1)",                            # computed mode
    "import notathing_at_all_xyz; print(1)",
    "print(1",                                                  # will not parse
])
def test_what_the_analyser_declines_to_read(src):
    assert not _eff(src).readable, src


# ---------------------------------------------------------------------------
# Writes: named, and named correctly
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("src,expected", [
    ("open('out.txt', 'w').write('x')", ["out.txt"]),
    ("open('out.txt', 'a').write('x')", ["out.txt"]),
    ("open('out.txt', 'xb').write(b'x')", ["out.txt"]),
    ("open('out.txt', mode='w')", ["out.txt"]),
    ("import openpyxl; wb = openpyxl.Workbook(); wb.save('r.xlsx')", ["r.xlsx"]),
    ("import pandas as pd; pd.read_csv('a').to_csv('b.csv')", ["b.csv"]),
    ("import matplotlib.pyplot as plt; plt.savefig('f.png')", ["f.png"]),
])
def test_a_write_is_named(src, expected):
    eff = _eff(src)
    assert eff.writes == expected, eff
    assert eff.readable, eff.opaque


@pytest.mark.parametrize("src", [
    "print(open('in.csv').read())",
    "open('in.csv', 'r').read()",
    "open('in.csv', 'rb').read()",
    "with open('in.csv') as f: pass",
])
def test_a_read_is_not_a_write(src):
    assert _eff(src).writes == [], src


@pytest.mark.parametrize("src", [
    # str.replace, not Path.replace: two arguments, not one.
    "print('a,b'.replace(',', '.'))",
    # dict.copy / DataFrame.copy: no arguments, not two.
    "d = {}; e = d.copy()",
    # DataFrame.replace: two arguments.
    "import pandas as pd; df = pd.read_csv('a'); df.replace(0, 1)",
])
def test_the_same_name_on_another_type_is_not_a_write(src):
    eff = _eff(src)
    assert eff.writes == [], eff
    assert eff.readable, eff.opaque


def test_the_arity_that_does_mean_a_write_still_does():
    """The discriminator must not have switched the rule off. One
    argument is Path.replace; two are shutil.copy's source and
    destination, and the destination is the second one."""
    assert _eff("p.replace('dest.txt')").writes == ["dest.txt"]
    assert _eff("mod.copy('src.txt', 'dest.txt')").writes == ["dest.txt"]
    assert _eff("mod.move('a.txt', 'b.txt')").writes == ["b.txt"]


# ---------------------------------------------------------------------------
# A local import is the script it is
# ---------------------------------------------------------------------------

def test_a_local_module_comes_back_to_be_scanned(tmp_path):
    (tmp_path / "mymod.py").write_text("X = 1\n", encoding="utf-8")
    eff = _eff("import mymod; print(mymod.X)", tmp_path)
    assert eff.readable, eff.opaque
    assert eff.local_modules == [str(tmp_path / "mymod.py")]


def test_a_local_package_counts_too(tmp_path):
    (tmp_path / "pkg").mkdir()
    (tmp_path / "pkg" / "__init__.py").write_text("Y = 2\n", encoding="utf-8")
    eff = _eff("from pkg import Y; print(Y)", tmp_path)
    assert eff.local_modules == [str(tmp_path / "pkg" / "__init__.py")]


def test_a_stdlib_name_is_not_reported_as_local(tmp_path):
    """Reporting json as a local file would send the content scan looking
    for a file that is not there -- harmless, but it would also mean the
    rule had stopped distinguishing the two."""
    assert _eff("import json; print(json)", tmp_path).local_modules == []


def test_an_unknown_import_with_no_local_file_is_opaque(tmp_path):
    """Without a working directory to check, or with one that does not
    hold the file, an unrecognised module is not assumed harmless."""
    assert not _eff("import mymod; print(mymod)", tmp_path).readable
    assert not _eff("import mymod; print(mymod)", None).readable


def test_a_relative_import_is_opaque():
    assert not _eff("from . import x").readable


# ---------------------------------------------------------------------------
# Splitting the command line
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd,expected", [
    ('python3 -c "print(1)"', ["print(1)"]),
    ("python -c 'print(1)'", ["print(1)"]),
    ('python3.11 -c "x=1"', ["x=1"]),
    ('/usr/bin/python3 -c "x=1"', ["x=1"]),
    ('python3 -u -c "x=1"', ["x=1"]),
    ('python3 -c "a" && python3 -c "b"', ["a", "b"]),
    ("ls -la", []),
    ("python3 script.py", []),
    ("python3 -m pytest", []),
])
def test_the_payloads_are_found(cmd, expected):
    assert extract_c_payloads(cmd) == expected


def test_an_untokenisable_line_says_so_rather_than_claiming_nothing():
    """None and [] are different answers: [] means "read it, no payload",
    and a caller that treats an unparsed line as [] has decided it read
    something it could not read."""
    assert extract_c_payloads('python3 -c "unbalanced') is None
    assert extract_c_payloads("ls") == []


def test_nothing_raises_on_junk():
    for junk in ("", None, "\x00\x01", "def f(:", "x" * 20000):
        analyze_payload(junk)           # type: ignore[arg-type]
        extract_c_payloads(junk)        # type: ignore[arg-type]
