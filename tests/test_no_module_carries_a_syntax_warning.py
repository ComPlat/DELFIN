r"""Compiling any module of the package must not raise a syntax warning.

A string holding an escape Python does not define, or a comparison like
``x is "text"``, makes ``compile()`` emit a syntax warning today; the same
construct is a hard ``SyntaxError`` in a future release, so a module
carrying one is a module that stops importing at the next interpreter
upgrade. It is also the noisiest defect this repository has measured: every
agent process that imported ``delfin.agent.api_client`` printed
``invalid escape sequence`` on its way in -- 504 records from one day of
swarm runs, 110 agents -- because the ``_word_spans`` docstring is not a raw
string and shows a shell-escaped space, ``a`` backslash ``b``.

``ast.parse`` (``tests/test_syntax_integrity.py``) does not see any of it:
only ``compile()`` emits these warnings. This guard therefore compiles every
file under ``delfin/`` itself, under ``warnings.simplefilter("error")``, and
lists each offender as ``file:line: message``.

Two deliberate choices, both because this has to fail the same way on the
CI interpreter (3.11) and on a developer's (3.12+):

* the warning's CATEGORY is ignored -- 3.11/3.12 report an undefined escape
  as a ``DeprecationWarning``, 3.13 as a ``SyntaxWarning``, and ``is`` with a
  literal has been a ``SyntaxWarning`` throughout. Any warning from
  ``compile()`` is an offender.
* nothing asks the host what is installed; the files on disk are the input.

This module's own docstrings are raw for the reason
``test_this_guard_itself_is_clean`` pins: a guard that warns while it guards
is a guard nobody believes.
"""

from __future__ import annotations

import ast
import warnings
from pathlib import Path

import pytest

import _package_ast

_ROOT = Path(__file__).resolve().parents[1]
_PKG = _ROOT / "delfin"

_PY_FILES = sorted(
    p for p in _PKG.rglob("*.py") if "__pycache__" not in p.parts
)


def _compile_warnings(text: str, name: str) -> tuple:
    """``([ (line, msg), ... ], syntax_error_or_None)`` for *text*.

    Compiled under ``"always"``, the pass that can see every warning: one
    under an error filter stops at the first, and a file with three bad
    strings has to report three lines. A ``SyntaxError`` is returned, not
    raised -- the caller decides whether it is this guard's business.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            compile(text, name, "exec")
        except SyntaxError as exc:
            return [(w.lineno, str(w.message)) for w in caught], exc
    return [(w.lineno, str(w.message)) for w in caught], None


def _offenders(path: Path) -> list:
    """``(line, message)`` for each compile-time warning in *path*.

    Two passes, because the two Python families report the same defect
    differently and both have to give the same answer here:

    * the ``"always"`` pass lists every warning with its line. This is what
      the report prints -- the line the compiler attributes to the string,
      which for a multi-line literal is its FIRST line, not the line holding
      the backslash.
    * the pass the task specifies, ``simplefilter("error")``, which is what
      makes a warning unmissable: 3.11 raises it as a ``DeprecationWarning``,
      3.12/3.13 as a ``SyntaxError`` whose message names the construct. It
      contributes no text -- an escalated warning is not recorded -- only the
      guarantee that no ``filterwarnings("ignore")`` elsewhere in the run can
      hide an offender from the check.

    A file that fails to compile at all is NOT this guard's report: a real
    syntax error is re-raised for ``test_syntax_integrity.py``, which names
    it properly. Only a file that warned and also stopped compiling is
    counted as the warning it gave.

    The helper call is wrapped in an ignore filter on purpose: its own
    ``ast.parse`` compiles the same text with the default filters and no
    filename, so it can only ever say ``<unknown>:5330``. This function is
    the one that decides, with the file and the line.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        text = _package_ast.source(path)
    name = str(path)
    found, hard = _compile_warnings(text, name)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        try:
            compile(text, name, "exec")
        except Warning:
            pass                    # 3.11: raised as a DeprecationWarning
        except SyntaxError:
            # 3.12/3.13 escalate here. If the warning pass saw nothing, this
            # is not an escalated warning but a broken file, and that file
            # belongs to test_syntax_integrity.py.
            if hard is not None and not found:
                raise hard
    return found


def _label(path: Path) -> str:
    """*path* relative to the repository, or absolute if it is outside."""
    try:
        return str(path.relative_to(_ROOT))
    except ValueError:
        return str(path)


def test_package_has_files():
    assert len(_PY_FILES) > 50


def test_no_module_of_the_package_compiles_with_a_syntax_warning():
    offenders = [
        f"{_label(path)}:{line}: {msg}"
        for path in _PY_FILES
        for line, msg in _offenders(path)
    ]
    assert not offenders, (
        "these modules make compile() warn -- a non-raw string with an "
        "escape Python does not define, or a construct like x is 'text', "
        "any of which becomes a SyntaxError in a future release:\n  "
        + "\n  ".join(offenders)
    )


# --- the check itself: what it names, and what it leaves alone -----------


def test_the_check_names_the_line_and_the_message_of_a_bad_escape(tmp_path):
    """A failing report has to be actionable: file, line, what and where.

    The fixture is written as bytes so this file does not carry the escape
    it looks for. The guard fails this if it reports the wrong line, no
    message, or nothing at all.
    """
    bad = tmp_path / "bad.py"
    bad.write_bytes(b'"""doc"""\n\nMSG = "path a\\ b c"\n')
    found = _offenders(bad)
    assert len(found) == 1, found
    line, message = found[0]
    # The line is where the string starts (line 3), not where the backslash
    # is; the message names the sequence. Pinned on those two properties, not
    # on CPython's exact wording, which has moved before.
    assert line == 3, found
    assert "invalid escape sequence" in message and chr(92) in message, found


def test_a_warning_that_is_not_an_escape_is_also_an_offender(tmp_path):
    """The rule is 'compile warned', not 'the message says escape'.

    ``x is 'text'`` has been a SyntaxWarning since 3.8 on every interpreter
    this runs on, so this pins the category-agnostic rule: narrow the guard
    back to the escape message and this fails.
    """
    odd = tmp_path / "odd.py"
    odd.write_bytes(b'A = 1\nB = A is "x"\n')
    found = _offenders(odd)
    assert len(found) == 1, found
    line, message = found[0]
    assert line == 2, found
    assert message and "invalid escape sequence" not in message, found


def test_a_raw_string_with_the_same_backslash_is_no_offender(tmp_path):
    r"""The smallest correct fix, and the control that passes either way.

    A raw string holding one backslash renders exactly what a non-raw string
    with a doubled backslash renders -- one backslash and one space -- so
    raw-prefixing a string whose backslash is meant literally changes no
    text, and nothing here is flagged.
    """
    raw = tmp_path / "raw.py"
    raw.write_bytes(b'R = r"a\\ b"\n')            # source: r"a" backslash " b"
    assert _offenders(raw) == []


def test_an_escape_python_defines_is_no_offender(tmp_path):
    r"""Newline, tab and a doubled backslash are valid; not flagged."""
    ok = tmp_path / "ok.py"
    ok.write_bytes(b'A = "a\\nb"\nB = "a\\tc"\nC = "a\\\\ b"\nD = r"a\\d"\n')
    assert _offenders(ok) == []


def test_a_real_syntax_error_is_left_to_the_parse_guard(tmp_path):
    """The guard must not swallow a broken file into an escape report."""
    broken = tmp_path / "broken.py"
    broken.write_bytes(b"def (:\n")
    with pytest.raises(SyntaxError):
        _offenders(broken)


def test_this_guard_itself_is_clean():
    """A guard that warns while it guards is a guard nobody believes.

    Fails the moment any string in this file, or in the shared parse helper,
    gains an undefined escape -- exactly the mistake this file exists to
    catch, and the one its first draft actually made.
    """
    for path in (Path(__file__), _ROOT / "tests" / "_package_ast.py"):
        assert _offenders(path) == [], f"{_label(path)} warns"


def test_the_shell_word_docstring_keeps_the_backslash_it_shows():
    r"""The fix must not change what a reader of that docstring sees.

    ``_word_spans`` documents a shell-escaped space -- one backslash, one
    space -- and both correct fixes (raw string, or a doubled backslash in
    the source) render exactly that. A fix that drops the backslash, or
    leaves the reader with two of them, fails here. This is what makes
    "raw string where the backslash is meant literally" a check and not a
    guess.
    """
    path = _PKG / "agent" / "api_client.py"
    tree = _package_ast.parse(path)
    fn = next(
        node for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef) and node.name == "_word_spans"
    )
    doc = ast.get_docstring(fn) or ""
    one = "``a" + chr(92) + " b``"
    two = "``a" + chr(92) * 2 + " b``"
    assert one in doc, f"{one!r} missing from the _word_spans docstring"
    assert two not in doc, f"{two!r}: the reader now sees two backslashes"
