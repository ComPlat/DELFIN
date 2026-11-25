"""Every module of the package compiles without a warning.

Input: every ``.py`` file under ``delfin/``, compiled with
``compile(source, "<path>", "exec")`` under ``warnings.simplefilter("error")``.
Output: the list of offenders as ``<path>:<line>: <message>`` -- empty.

Why a compile and not a parse: an invalid escape sequence is not a syntax
error, so the file imports and the tests pass. It is a ``DeprecationWarning``
on the Python the CI runs (3.11) and a ``SyntaxWarning`` from 3.12 on, and
it is scheduled to become an error. ``delfin/agent/api_client.py:5330``
carried one -- ``"a\\ b"`` inside a non-raw docstring -- and every session
that imported the module printed it: 504 records over the 110 agents of one
day, from a line of prose.

The filter is set to *error* without a category, not to SyntaxWarning: on
3.11 the same defect is reported as a DeprecationWarning, so asking for
SyntaxWarning alone would pass here and print a warning there. A test that
is green on the interpreter that does not warn has measured nothing.
"""

from __future__ import annotations

import pathlib
import warnings

_PKG_ROOT = pathlib.Path(__file__).resolve().parents[1] / "delfin"

_PY_FILES = sorted(
    p for p in _PKG_ROOT.rglob("*.py") if "__pycache__" not in p.parts
)


def _compile_problems(source: str, name: str) -> list[str]:
    """``["<name>:<line>: <message>"]`` for each warning compiling raises."""
    problems: list[str] = []
    with warnings.catch_warnings():
        # No category: an invalid escape is a DeprecationWarning before 3.12
        # and a SyntaxWarning from 3.12, and both are the same defect.
        warnings.simplefilter("error")
        try:
            compile(source, name, "exec")
        except SyntaxError as exc:
            problems.append(f"{name}:{exc.lineno}: {exc.msg}")
        except Warning as exc:  # a warning turned into an error
            line = getattr(exc, "lineno", None)
            problems.append(f"{name}:{line}: {exc}")
    return problems


def test_the_package_has_files():
    """A wrong root would make the scan below pass without reading anything."""
    assert len(_PY_FILES) > 50


def test_every_package_file_compiles_without_a_warning():
    offenders: list[str] = []
    for path in _PY_FILES:
        source = path.read_text(encoding="utf-8", errors="surrogateescape")
        offenders.extend(
            _compile_problems(source, str(path.relative_to(_PKG_ROOT.parent)))
        )
    assert not offenders, (
        "a module compiles with a warning turned into an error; a backslash "
        "the string language does not define is a raw string or a doubled "
        "backslash, and in a docstring the intended character stays.\n"
        "If an offender is a guard-protected module (api_client.py, engine.py, "
        "...), its fix ships as scratch/api_client-escapes.patch and this test "
        "goes green when that patch is applied:\n  "
        + "\n  ".join(offenders)
    )


def test_the_instrument_reports_an_invalid_escape():
    """The control that CAN fail: the same scan on a source that carries one.

    ``"a\\ b"`` in a non-raw docstring is exactly the api_client.py shape.
    A scan that stays silent here would pass the test above for any tree.
    """
    source = 'def f():\n    """a\\ b is one word"""\n    return 1\n'
    problems = _compile_problems(source, "demo.py")
    assert len(problems) == 1, problems
    assert problems[0].startswith("demo.py:2: "), problems[0]
    assert "invalid escape sequence" in problems[0], problems[0]


def test_the_instrument_passes_a_raw_docstring():
    """The other control: the fix the offenders need is not reported."""
    source = 'def f():\n    r"""a\\ b is one word"""\n    return 1\n'
    assert _compile_problems(source, "demo.py") == []
