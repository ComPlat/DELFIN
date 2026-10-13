"""Gate patterns that backtracked exponentially, and their replacements.

Five patterns in ``delfin/agent/api_client.py`` were reported by code
scanning (py/redos):

* ``_GIT_PUSH_RE``: the option loop listed ``-c <value>`` twice (once
  inside ``-[Cc]``), so every ``-c`` had two parses and a command such as
  ``git -c ! -c ! ...`` without ``push`` took 2**n steps to reject.
* the three blanking substitutions in ``_sed_script_writes``: the body
  ``(?:\\\\.|(?!D).)*`` matches a backslash both as the start of an
  escape and as a plain character, so ``\\a\\a\\a...`` without a closing
  delimiter has 2**n parses.

Thirty repetitions held the gate for minutes. The replacements must
answer exactly as the old patterns did -- same spans, same groups, same
substituted text -- on every command string in the test suite, on every
short string over the characters that matter, and on random longer ones;
and the attack strings must finish in milliseconds.

The old patterns are copied here verbatim as the reference.
"""

from __future__ import annotations

import ast
import itertools
import random
import re
import time
from pathlib import Path

import pytest

from delfin.agent import api_client as A

_TESTS = Path(__file__).resolve().parent

# --- the reference: the patterns as they were -----------------------------

_OLD_GIT_PUSH_RE = re.compile(
    r"(?:^|[\s;&|(`])(?:\S*/)?git(?:\s+-[Cc]\s+\S+|\s+-c\s+\S+"
    r"|\s+--[\w-]+(?:=\S+)?)*\s+push\b")
_OLD_SED_SUBST_RE = re.compile(
    r"(?<![\w])s(.)(?:\\.|(?!\1).)*\1(?:\\.|(?!\1).)*\1([a-zA-Z0-9]*)")
_OLD_SLASH_RE = re.compile(r"/(?:\\.|[^/])*/")
_OLD_CUSTOM_RE = re.compile(r"\\(.)(?:\\.|(?!\1).)*\1")


def _old_subst(text: str) -> str:
    return _OLD_SED_SUBST_RE.sub(lambda m: " s " + m.group(2), text)


def _old_slash(text: str) -> str:
    return _OLD_SLASH_RE.sub(" ", text)


def _old_custom(text: str) -> str:
    return _OLD_CUSTOM_RE.sub(" ", text)


def _old_sed_script_writes(cmd: str) -> bool:
    body = re.sub(r"^\s*sed\b", "", cmd, count=1)
    body = _old_subst(body)
    body = _old_slash(body)
    body = _old_custom(body)
    if re.search(r"\bs [a-zA-Z0-9]*[we]", body):
        return True
    return bool(re.search(r"(?:^|[\s;'\"{}0-9$,!])[wWe](?:\s|$|['\"])", body))


# --- the corpus ------------------------------------------------------------

def _suite_strings() -> list[str]:
    """Every string literal of up to 300 characters in tests/*.py: the
    commands the gate tests drive, and much else besides."""
    seen: set[str] = set()
    for path in sorted(_TESTS.glob("*.py")):
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except (SyntaxError, UnicodeDecodeError):
            continue
        for node in ast.walk(tree):
            if (isinstance(node, ast.Constant) and isinstance(node.value, str)
                    and len(node.value) <= 300):
                seen.add(node.value)
    return sorted(seen)


_SUITE = _suite_strings()


def _exhaustive(alphabet: str, max_len: int):
    for n in range(max_len + 1):
        for chars in itertools.product(alphabet, repeat=n):
            yield "".join(chars)


def _random(alphabet: str, count: int, max_len: int, seed: int):
    rnd = random.Random(seed)
    for _ in range(count):
        yield "".join(rnd.choice(alphabet)
                      for _ in range(rnd.randint(0, max_len)))


def _sed_corpus():
    yield from _SUITE
    # Every string of up to 7 characters over the characters the three
    # patterns distinguish: the opener, two delimiters, the escape, a
    # newline and an ordinary letter.
    yield from _exhaustive("s/\\a\nx", 7)
    yield from _random("s/\\a\nxw e;|'\"1_", 40000, 24, seed=7)
    yield from _random("s\\\\\\/aa", 20000, 40, seed=11)


def _git_corpus():
    yield from _SUITE
    tokens = ["git", " ", "\t", " -c", " -C", " x", " !", " k=v", " --a",
              " --b=c", " --", " push", "push", "/usr/bin/", ";", "(", "`",
              "-c", " pushx", "\n", "=", "p"]
    rnd = random.Random(3)
    for _ in range(40000):
        yield "".join(rnd.choice(tokens) for _ in range(rnd.randint(0, 12)))
    for tail in _exhaustive(["git", " -c", " -C", " !", " --a=b", " push"],
                            5):
        yield "".join(tail)


def _git_spans(pattern, text):
    return [(m.span(), m.groups()) for m in pattern.finditer(text)]


# --- identical answers ------------------------------------------------------

def test_the_corpus_holds_the_gate_tests_commands():
    """The suite strings include the commands the sed and push gates are
    tested with; a corpus that lost them would compare nothing."""
    assert "sed -i s/a/b/ notes.txt" in _SUITE
    assert any(s.startswith("git push") for s in _SUITE)
    assert len(_SUITE) > 5000


def test_git_push_pattern_matches_as_before():
    differ = [s for s in _git_corpus()
              if _git_spans(_OLD_GIT_PUSH_RE, s) != _git_spans(A._GIT_PUSH_RE, s)]
    assert differ == []


@pytest.mark.parametrize("old, name", [
    (_old_subst, "_sed_blank_substitutions"),
    (_old_slash, "_sed_blank_slash_regexes"),
    (_old_custom, "_sed_blank_custom_regexes"),
])
def test_sed_blanking_substitutes_as_before(old, name):
    new = getattr(A, name)
    differ = [s for s in _sed_corpus() if old(s) != new(s)]
    assert differ[:20] == []


def test_sed_script_writes_answers_as_before():
    differ = []
    for s in _sed_corpus():
        for cmd in (s, "sed " + s, "sed -n " + s):
            if _old_sed_script_writes(cmd) != A._sed_script_writes(cmd):
                differ.append(cmd)
    assert differ[:20] == []


# --- the attack strings -----------------------------------------------------

_BUDGET_S = 0.1


def _elapsed(fn, arg) -> float:
    best = float("inf")
    for _ in range(3):
        start = time.perf_counter()
        fn(arg)
        best = min(best, time.perf_counter() - start)
    return best


@pytest.mark.parametrize("label, fn, text", [
    ("git -c", lambda s: A._GIT_PUSH_RE.search(s), "\tgit" + " -c !" * 50000),
    ("git -C", lambda s: A._GIT_PUSH_RE.search(s), "git" + " -C x -c y" * 25000),
    ("s///", lambda s: A._sed_blank_substitutions(s), "sa" + "\\a" * 50000),
    ("s/// open", lambda s: A._sed_blank_substitutions(s), "sa" + "\\b" * 50000),
    ("s/// nl", lambda s: A._sed_blank_substitutions(s), "s/" + "\\b" * 50000 + "\n/"),
    ("s/// chain", lambda s: A._sed_blank_substitutions(s), "s/" + "\\\\/" * 30000),
    ("/re/", lambda s: A._sed_blank_slash_regexes(s), "/" + "\\." * 50000),
    ("/re/ closed", lambda s: A._sed_blank_slash_regexes(s), "/" + "\\." * 50000 + "/"),
    ("\\cREc", lambda s: A._sed_blank_custom_regexes(s), "\\a" + "\\a" * 50000),
    ("\\cREc open", lambda s: A._sed_blank_custom_regexes(s), "\\a" + "\\b" * 50000),
    ("sed whole", A._sed_script_writes, "sed sa" + "\\a" * 50000),
    ("sed whole open", A._sed_script_writes, "sed s/" + "\\b" * 50000),
])
def test_attack_strings_finish_in_milliseconds(label, fn, text):
    assert _elapsed(fn, text) < _BUDGET_S, label


def test_the_old_patterns_were_exponential():
    """The control: the reference patterns take four times as long for
    every two more repetitions, so the attack strings above would not
    finish. Measured at n=12 and n=18 to keep the control short."""
    def ratio(pattern, make):
        small = _elapsed(pattern.search, make(12))
        large = _elapsed(pattern.search, make(18))
        return large / max(small, 1e-9)
    assert ratio(_OLD_GIT_PUSH_RE, lambda n: "\tgit" + " -c !" * n) > 20
    assert ratio(_OLD_SED_SUBST_RE, lambda n: "sa" + "\\b" * n) > 20
    assert ratio(_OLD_SLASH_RE, lambda n: "/" + "\\." * n) > 20
    assert ratio(_OLD_CUSTOM_RE, lambda n: "\\a" + "\\b" * n) > 20


# --- capabilities: directory names, not secrets -----------------------------

def test_capabilities_lists_the_paths_the_sandbox_masks(tmp_path):
    """``capabilities()["masked_home_dirs"]`` reads the public tuple; the
    sandbox masks from the same object, so the report cannot drift from
    the mounts. Code scanning traced the old report back to a variable
    named for secrets and flagged ``--json`` output as a leak of one;
    the values were always directory names."""
    from delfin.agent import boundaries, sandbox
    assert sandbox._HOME_SECRET_DIRS is sandbox.MASKED_HOME_ENTRIES
    caps = boundaries.capabilities(cwd=tmp_path)
    assert caps["masked_home_dirs"] == list(sandbox.MASKED_HOME_ENTRIES)
    src = Path(boundaries.__file__).read_text(encoding="utf-8")
    assert "_HOME_SECRET_DIRS" not in src
