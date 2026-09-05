"""What a typed formula works out to, checked against Excel's own rules.

Every expectation here is what Excel puts in the cell, including the awkward
ones: text that looks like a number is a number, a range skips text instead of
failing on it, halves round away from zero, and a cell that refers to itself
says so instead of hanging.
"""

from __future__ import annotations

import math

import pytest

from delfin.dashboard import formula_engine as fx


def _grid(cells):
    """{'A1': '3', 'B2': '=A1*2'} -> (formulas, values) grids."""
    rows = cols = 0
    parsed = {}
    for address, text in cells.items():
        spot = fx.parse_ref(address)
        assert spot is not None, address
        parsed[spot] = text
        rows = max(rows, spot[0] + 1)
        cols = max(cols, spot[1] + 1)
    formulas = [['' for _ in range(cols)] for _ in range(rows)]
    values = [['' for _ in range(cols)] for _ in range(rows)]
    for (row, col), text in parsed.items():
        if fx.is_formula(text):
            formulas[row][col] = text
        else:
            values[row][col] = text
    return formulas, values


def calc(formula, **cells):
    formulas, values = _grid({**cells, 'ZZ9999': ''})
    return fx.evaluate_one(formula, lambda r, c: (
        values[r][c] if r < len(values) and c < len(values[r]) else ''
    ))


# ---------------------------------------------------------------------------
# arithmetic
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('formula, expected', [
    ('=1+2', 3),
    ('=2*3+4', 10),
    ('=2+3*4', 14),
    ('=(2+3)*4', 20),
    ('=10/4', 2.5),
    ('=2^10', 1024),
    ('=-3+1', -2),
    ('=--3', 3),
    ('=50%', 0.5),
    ('=10%*200', 20),
    ('=2^3^2', 512),              # right-associative, as in Excel
])
def test_arithmetic(formula, expected):
    assert calc(formula) == pytest.approx(expected)


def test_division_by_zero_is_an_excel_error_not_a_crash():
    assert calc('=1/0') == '#DIV/0!'
    assert calc('=IFERROR(1/0,"n/a")') == 'n/a'


# ---------------------------------------------------------------------------
# references
# ---------------------------------------------------------------------------
def test_a_reference_reads_the_cell():
    assert calc('=A1+B1', A1='3', B1='4') == 7
    assert calc('=$A$1*2', A1='2.5') == 5


def test_a_blank_cell_is_zero_in_arithmetic():
    assert calc('=A1+5', A1='') == 5


def test_text_that_looks_like_a_number_is_a_number():
    assert calc('=A1+1', A1='3') == 4


def test_text_that_is_not_a_number_is_an_error_in_arithmetic():
    assert calc('=A1+1', A1='apples') == '#VALUE!'


def test_a_range_skips_text_and_blanks_instead_of_failing():
    values = {'A1': '1', 'A2': 'apples', 'A3': '', 'A4': '3'}
    assert calc('=SUM(A1:A4)', **values) == 4
    assert calc('=COUNT(A1:A4)', **values) == 2
    assert calc('=COUNTA(A1:A4)', **values) == 3
    assert calc('=AVERAGE(A1:A4)', **values) == 2


def test_a_range_may_be_written_in_either_corner_order():
    assert calc('=SUM(A3:A1)', A1='1', A2='2', A3='3') == 6


# ---------------------------------------------------------------------------
# functions
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('formula, expected', [
    ('=SUM(1,2,3)', 6),
    ('=MIN(4,2,9)', 2),
    ('=MAX(4,2,9)', 9),
    ('=ROUND(2.345,2)', 2.35),
    ('=ROUND(2.5,0)', 3),            # Excel rounds halves away from zero
    ('=ROUND(-2.5,0)', -3),
    ('=ROUNDDOWN(2.9,0)', 2),
    ('=ROUNDUP(2.1,0)', 3),
    ('=ABS(-4)', 4),
    ('=SQRT(16)', 4),
    ('=POWER(3,3)', 27),
    ('=MOD(7,3)', 1),
    ('=INT(2.9)', 2),
    ('=LN(1)', 0),
    ('=LOG(8,2)', 3),
    ('=LOG10(1000)', 3),
    ('=IF(1>0,"yes","no")', 'yes'),
    ('=IF(1<0,"yes","no")', 'no'),
    ('=AND(TRUE,1>0)', True),
    ('=OR(FALSE,1>2)', False),
    ('=NOT(FALSE)', True),
    ('=LEN("hello")', 5),
    ('=UPPER("ab")', 'AB'),
    ('=LEFT("hello",2)', 'he'),
    ('=RIGHT("hello",2)', 'lo'),
    ('=MID("hello",2,3)', 'ell'),
    ('=TRIM("  a  b  ")', 'a b'),
    ('=CONCAT("a","b")', 'ab'),
    ('="a"&"b"', 'ab'),
    ('="n="&2', 'n=2'),
])
def test_functions(formula, expected):
    result = calc(formula)
    if isinstance(expected, float):
        assert result == pytest.approx(expected)
    else:
        assert result == expected


def test_statistics_over_a_range():
    cells = {f'A{i}': str(v) for i, v in enumerate([2, 4, 4, 4, 5, 5, 7, 9], start=1)}
    assert calc('=AVERAGE(A1:A8)', **cells) == 5
    assert calc('=MEDIAN(A1:A8)', **cells) == 4.5
    assert calc('=STDEVP(A1:A8)', **cells) == pytest.approx(2.0)
    assert calc('=STDEV(A1:A8)', **cells) == pytest.approx(2.13809, abs=1e-4)


def test_conditional_aggregation():
    cells = {'A1': 'x', 'A2': 'y', 'A3': 'x', 'B1': '1', 'B2': '2', 'B3': '4'}
    assert calc('=SUMIF(A1:A3,"x",B1:B3)', **cells) == 5
    assert calc('=COUNTIF(A1:A3,"x")', **cells) == 2
    assert calc('=SUMIF(B1:B3,">1")', **cells) == 6


def test_an_unknown_function_says_so_rather_than_guessing():
    assert calc('=WOBBLE(1)') == '#NAME?'


def test_a_broken_formula_is_an_error_not_an_exception():
    assert calc('=1+') == '#VALUE!'
    assert calc('=SUM(') == '#VALUE!'
    assert calc('=)(') == '#VALUE!'


# ---------------------------------------------------------------------------
# whole sheets
# ---------------------------------------------------------------------------
def test_a_formula_may_depend_on_another_formula():
    formulas, values = _grid({'A1': '2', 'B1': '=A1*3', 'C1': '=B1+1'})
    results = fx.evaluate_grid(formulas, values)
    assert results[(0, 1)] == 6
    assert results[(0, 2)] == 7


def test_a_chain_through_a_range_is_worked_out_in_order():
    cells = {'A1': '1', 'A2': '2', 'A3': '=SUM(A1:A2)', 'A4': '=A3*10'}
    formulas, values = _grid(cells)
    results = fx.evaluate_grid(formulas, values)
    assert results[(2, 0)] == 3
    assert results[(3, 0)] == 30


def test_a_cell_that_depends_on_itself_says_so_instead_of_hanging():
    formulas, values = _grid({'A1': '=A1+1'})
    assert fx.evaluate_grid(formulas, values)[(0, 0)] == '#CYCLE!'

    formulas, values = _grid({'A1': '=B1', 'B1': '=A1'})
    results = fx.evaluate_grid(formulas, values)
    assert '#CYCLE!' in (results[(0, 0)], results[(0, 1)])


def test_a_sheet_of_formulas_is_worked_out_quickly():
    """The point of this engine: a result while typing, not after a save."""
    import time

    rows = 2000
    formulas = [[f'=A{r + 1}*2' if c == 1 else '' for c in range(2)] for r in range(rows)]
    values = [[str(r), ''] for r in range(rows)]

    start = time.perf_counter()
    results = fx.evaluate_grid(formulas, values)
    elapsed = time.perf_counter() - start

    assert len(results) == rows
    assert results[(1999, 1)] == 3998
    assert elapsed < 1.0, f'{rows} formulas took {elapsed:.2f}s'
