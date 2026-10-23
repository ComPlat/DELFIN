"""Work out spreadsheet formulas against the grid on screen.

Why not the ``formulas`` package, which the tab already uses: it builds a model
of the whole workbook, which takes about ten seconds for five hundred formulas
and runs on the *file*.  So a formula only produced a number after saving, and
saving stopped the dashboard for those ten seconds.  Typing ``=SUM(A1:A5)`` and
seeing the sum is the thing a spreadsheet is for; it cannot wait for a save.

This engine evaluates one sheet, from the values that are on screen including
the edits that have not been saved yet.  It is deliberately small: the functions
below cover what a results table needs.  Anything it does not know it says so
about -- ``#NAME?`` -- rather than guessing, and the workbook-wide engine stays
available for the cases it cannot serve (cross-sheet references, the long tail
of functions).

Excel's own rules are followed where they are observable:

  * text that looks like a number is a number ('3' + 1 = 4),
  * an empty cell is 0 in arithmetic and '' in text,
  * SUM and friends *skip* text and blanks rather than failing on them, while
    the arithmetic operators do not,
  * division by zero is #DIV/0!, an unparsable operand #VALUE!,
  * a cell that ends up depending on itself is #CYCLE! rather than a hang.
"""

from __future__ import annotations

import math
import re
import statistics
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

__all__ = [
    'ERRORS', 'FormulaError', 'evaluate_grid', 'evaluate_one', 'is_formula',
]


class FormulaError(Exception):
    """An Excel-style error value, carried as an exception while evaluating."""

    def __init__(self, code: str):
        super().__init__(code)
        self.code = code


ERRORS = {
    'div0': '#DIV/0!',
    'value': '#VALUE!',
    'name': '#NAME?',
    'ref': '#REF!',
    'na': '#N/A',
    'num': '#NUM!',
    'cycle': '#CYCLE!',
}
_ERROR_TEXTS = set(ERRORS.values())


def is_formula(text: Any) -> bool:
    return isinstance(text, str) and text.startswith('=') and len(text) > 1


# ---------------------------------------------------------------------------
# Addresses
# ---------------------------------------------------------------------------
_REF_RE = re.compile(r'\$?([A-Za-z]{1,3})\$?([0-9]{1,7})')


def _col_index(letters: str) -> int:
    value = 0
    for char in letters.upper():
        value = value * 26 + (ord(char) - 64)
    return value - 1


def _col_letters(index0: int) -> str:
    out = ''
    index0 = int(index0)
    while True:
        out = chr(ord('A') + index0 % 26) + out
        index0 = index0 // 26 - 1
        if index0 < 0:
            return out


def parse_ref(text: str) -> Optional[Tuple[int, int]]:
    """'B7' or '$B$7' -> (row0, col0).  None if it is not an address."""
    match = _REF_RE.fullmatch(text.strip())
    if not match:
        return None
    return int(match.group(2)) - 1, _col_index(match.group(1))


# ---------------------------------------------------------------------------
# Values
# ---------------------------------------------------------------------------
_NUMBER_RE = re.compile(r'^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$')


def _to_number(value: Any, *, strict: bool = True) -> float:
    """Excel's coercion: blanks are 0, numeric text is its number."""
    if value is None or value == '':
        return 0.0
    if isinstance(value, bool):
        return 1.0 if value else 0.0
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if text in _ERROR_TEXTS:
        raise FormulaError(text)
    if _NUMBER_RE.match(text):
        return float(text)
    # A comma decimal separator is what a German export writes.
    if _NUMBER_RE.match(text.replace(',', '.')) and text.count(',') == 1:
        return float(text.replace(',', '.'))
    if not strict:
        return 0.0
    raise FormulaError(ERRORS['value'])


def _to_text(value: Any) -> str:
    if value is None:
        return ''
    if isinstance(value, bool):
        return 'TRUE' if value else 'FALSE'
    if isinstance(value, float):
        if value == int(value) and abs(value) < 1e15:
            return str(int(value))
        return repr(value)
    return str(value)


def _is_blank(value: Any) -> bool:
    return value is None or value == ''


def _numbers_in(values: Sequence[Any]) -> List[float]:
    """The numbers of a range: text and blanks are skipped, as in Excel."""
    out: List[float] = []
    for value in values:
        if _is_blank(value):
            continue
        if isinstance(value, str) and value in _ERROR_TEXTS:
            raise FormulaError(value)
        if isinstance(value, bool):
            continue          # booleans in a range are not counted, as in Excel
        try:
            out.append(_to_number(value))
        except FormulaError:
            continue          # text inside a range is skipped, not an error
    return out


# ---------------------------------------------------------------------------
# Tokenising
# ---------------------------------------------------------------------------
_TOKEN_RE = re.compile(r"""
    (?P<string>"(?:[^"]|"")*")
  | (?P<number>\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+)
  | (?P<range>\$?[A-Za-z]{1,3}\$?[0-9]{1,7}\s*:\s*\$?[A-Za-z]{1,3}\$?[0-9]{1,7})
  | (?P<ref>\$?[A-Za-z]{1,3}\$?[0-9]{1,7})
  | (?P<name>[A-Za-z_][A-Za-z0-9_.]*)
  | (?P<op><=|>=|<>|[-+*/^&%()<>=,;])
  | (?P<space>\s+)
""", re.VERBOSE)


def _tokenise(text: str) -> List[Tuple[str, str]]:
    tokens: List[Tuple[str, str]] = []
    pos = 0
    while pos < len(text):
        match = _TOKEN_RE.match(text, pos)
        if not match:
            raise FormulaError(ERRORS['value'])
        pos = match.end()
        kind = match.lastgroup
        if kind == 'space':
            continue
        # LOG10 is a perfectly good cell address -- column LOG, row 10 -- and
        # the reference pattern claims it first.  A name followed by '(' is a
        # call, which is how Excel tells the two apart as well.
        if kind == 'ref' and text[pos:pos + 1] == '(':
            kind = 'name'
        tokens.append((kind, match.group()))
    return tokens


# ---------------------------------------------------------------------------
# Parsing (precedence climbing) and evaluation in one pass
# ---------------------------------------------------------------------------
class _Range(list):
    """A rectangular block of values, flattened.  Kept apart from a list arg."""


_COMPARISONS = {
    '=': lambda a, b: a == b,
    '<>': lambda a, b: a != b,
    '<': lambda a, b: a < b,
    '>': lambda a, b: a > b,
    '<=': lambda a, b: a <= b,
    '>=': lambda a, b: a >= b,
}


class _Evaluator:
    def __init__(self, lookup: Callable[[int, int], Any]):
        self.lookup = lookup
        self.tokens: List[Tuple[str, str]] = []
        self.pos = 0

    # -- token helpers --
    def _peek(self) -> Optional[Tuple[str, str]]:
        return self.tokens[self.pos] if self.pos < len(self.tokens) else None

    def _take(self) -> Tuple[str, str]:
        token = self._peek()
        if token is None:
            raise FormulaError(ERRORS['value'])
        self.pos += 1
        return token

    def _expect(self, text: str) -> None:
        token = self._take()
        if token[1] != text:
            raise FormulaError(ERRORS['value'])

    # -- entry point --
    def run(self, formula: str) -> Any:
        self.tokens = _tokenise(formula)
        self.pos = 0
        value = self._expression()
        if self._peek() is not None:
            raise FormulaError(ERRORS['value'])
        return self._single(value)

    def _single(self, value: Any) -> Any:
        """A range used where one value is wanted collapses to its first cell."""
        if isinstance(value, _Range):
            return value[0] if value else 0.0
        return value

    # -- grammar --
    def _expression(self) -> Any:
        left = self._concat()
        while True:
            token = self._peek()
            if token is None or token[0] != 'op' or token[1] not in _COMPARISONS:
                return left
            op = self._take()[1]
            right = self._concat()
            left = self._compare(op, self._single(left), self._single(right))

    def _compare(self, op: str, left: Any, right: Any) -> bool:
        if isinstance(left, str) or isinstance(right, str):
            try:
                left_n, right_n = _to_number(left), _to_number(right)
            except FormulaError:
                left, right = _to_text(left).upper(), _to_text(right).upper()
                return _COMPARISONS[op](left, right)
            return _COMPARISONS[op](left_n, right_n)
        return _COMPARISONS[op](_to_number(left), _to_number(right))

    def _concat(self) -> Any:
        left = self._additive()
        while True:
            token = self._peek()
            if token is None or token[1] != '&':
                return left
            self._take()
            right = self._additive()
            left = _to_text(self._single(left)) + _to_text(self._single(right))

    def _additive(self) -> Any:
        left = self._multiplicative()
        while True:
            token = self._peek()
            if token is None or token[1] not in ('+', '-'):
                return left
            op = self._take()[1]
            right = self._multiplicative()
            a, b = _to_number(self._single(left)), _to_number(self._single(right))
            left = a + b if op == '+' else a - b

    def _multiplicative(self) -> Any:
        left = self._unary()
        while True:
            token = self._peek()
            if token is None or token[1] not in ('*', '/'):
                return left
            op = self._take()[1]
            right = self._unary()
            a, b = _to_number(self._single(left)), _to_number(self._single(right))
            if op == '*':
                left = a * b
            else:
                if b == 0:
                    raise FormulaError(ERRORS['div0'])
                left = a / b

    def _unary(self) -> Any:
        token = self._peek()
        if token is not None and token[1] in ('-', '+'):
            self._take()
            value = _to_number(self._single(self._unary()))
            return -value if token[1] == '-' else value
        return self._power()

    def _power(self) -> Any:
        base = self._postfix()
        token = self._peek()
        if token is not None and token[1] == '^':
            self._take()
            exponent = _to_number(self._single(self._unary()))
            try:
                return math.pow(_to_number(self._single(base)), exponent)
            except (ValueError, OverflowError):
                raise FormulaError(ERRORS['num'])
        return base

    def _postfix(self) -> Any:
        value = self._primary()
        token = self._peek()
        if token is not None and token[1] == '%':
            self._take()
            return _to_number(self._single(value)) / 100.0
        return value

    def _primary(self) -> Any:
        kind, text = self._take()
        if kind == 'number':
            return float(text)
        if kind == 'string':
            return text[1:-1].replace('""', '"')
        if kind == 'range':
            return self._range_values(text)
        if kind == 'ref':
            spot = parse_ref(text)
            if spot is None:
                raise FormulaError(ERRORS['ref'])
            return self.lookup(*spot)
        if kind == 'name':
            upper = text.upper()
            token = self._peek()
            if token is not None and token[1] == '(':
                return self._call(upper)
            if upper == 'TRUE':
                return True
            if upper == 'FALSE':
                return False
            raise FormulaError(ERRORS['name'])
        if text == '(':
            value = self._expression()
            self._expect(')')
            return value
        raise FormulaError(ERRORS['value'])

    def _range_values(self, text: str) -> _Range:
        first, last = [part.strip() for part in text.split(':')]
        start, end = parse_ref(first), parse_ref(last)
        if start is None or end is None:
            raise FormulaError(ERRORS['ref'])
        r0, r1 = sorted((start[0], end[0]))
        c0, c1 = sorted((start[1], end[1]))
        return _Range(
            self.lookup(row, col)
            for row in range(r0, r1 + 1)
            for col in range(c0, c1 + 1)
        )

    def _skip_argument(self) -> None:
        """Step over the rest of an argument whose evaluation raised.

        IFERROR only means anything if its first argument is allowed to fail,
        and the failure comes out of the middle of the expression -- leaving
        the parser standing between tokens.  This puts it back on the comma.
        """
        depth = 0
        while self.pos < len(self.tokens):
            text = self.tokens[self.pos][1]
            if text == '(':
                depth += 1
            elif text == ')':
                if depth == 0:
                    return
                depth -= 1
            elif text in (',', ';') and depth == 0:
                return
            self.pos += 1

    def _call(self, name: str) -> Any:
        self._expect('(')
        # These hand their arguments' failures on as values instead of failing.
        catches = name in ('IFERROR', 'IFNA')
        args: List[Any] = []
        if self._peek() is not None and self._peek()[1] == ')':
            self._take()
        else:
            while True:
                if catches:
                    try:
                        args.append(self._expression())
                    except FormulaError as exc:
                        args.append(exc.code)
                        self._skip_argument()
                else:
                    args.append(self._expression())
                token = self._take()
                if token[1] == ')':
                    break
                if token[1] not in (',', ';'):
                    raise FormulaError(ERRORS['value'])
        function = FUNCTIONS.get(name)
        if function is None:
            raise FormulaError(ERRORS['name'])
        return function(args)


# ---------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------
def _flat(args: Sequence[Any]) -> List[Any]:
    out: List[Any] = []
    for arg in args:
        if isinstance(arg, _Range):
            out.extend(arg)
        else:
            out.append(arg)
    return out


def _one(args: Sequence[Any], index: int = 0) -> Any:
    if len(args) <= index:
        raise FormulaError(ERRORS['value'])
    value = args[index]
    if isinstance(value, _Range):
        return value[0] if value else 0.0
    return value


def _num(args: Sequence[Any], index: int = 0) -> float:
    return _to_number(_one(args, index))


def _truth(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if _is_blank(value):
        return False
    if isinstance(value, str):
        upper = value.strip().upper()
        if upper == 'TRUE':
            return True
        if upper == 'FALSE':
            return False
    return _to_number(value) != 0


def _guard_empty(values: Sequence[float]) -> List[float]:
    if not values:
        raise FormulaError(ERRORS['div0'])
    return list(values)


def _match_criterion(value: Any, criterion: Any) -> bool:
    """Excel's criterion: a bare value means equals, '>3' and '<>x' compare."""
    text = _to_text(criterion).strip()
    for op in ('<=', '>=', '<>', '<', '>', '='):
        if text.startswith(op):
            wanted = text[len(op):]
            try:
                return _COMPARISONS[op](_to_number(value), _to_number(wanted))
            except FormulaError:
                return _COMPARISONS[op](_to_text(value).upper(), wanted.upper())
    try:
        return _to_number(value) == _to_number(text)
    except FormulaError:
        return _to_text(value).upper() == text.upper()


def _sumif(args: Sequence[Any]) -> float:
    if len(args) < 2:
        raise FormulaError(ERRORS['value'])
    tested = list(args[0]) if isinstance(args[0], _Range) else [args[0]]
    summed = (list(args[2]) if len(args) > 2 and isinstance(args[2], _Range)
              else tested)
    total = 0.0
    for index, value in enumerate(tested):
        if _match_criterion(value, _one(args, 1)) and index < len(summed):
            try:
                total += _to_number(summed[index])
            except FormulaError:
                continue
    return total


def _countif(args: Sequence[Any]) -> float:
    if len(args) < 2:
        raise FormulaError(ERRORS['value'])
    tested = list(args[0]) if isinstance(args[0], _Range) else [args[0]]
    criterion = _one(args, 1)
    return float(sum(1 for value in tested
                     if not _is_blank(value) and _match_criterion(value, criterion)))


def _if(args: Sequence[Any]) -> Any:
    if len(args) < 2:
        raise FormulaError(ERRORS['value'])
    if _truth(_one(args, 0)):
        return _one(args, 1)
    return _one(args, 2) if len(args) > 2 else False


def _iferror(args: Sequence[Any]) -> Any:
    # The argument was already evaluated, so an error arrives as its text.
    value = _one(args, 0)
    if isinstance(value, str) and value in _ERROR_TEXTS:
        return _one(args, 1) if len(args) > 1 else ''
    return value


def _round(args: Sequence[Any], mode: str = 'half') -> float:
    value = _num(args, 0)
    digits = int(_num(args, 1)) if len(args) > 1 else 0
    factor = 10.0 ** digits
    scaled = value * factor
    if mode == 'up':
        result = math.ceil(abs(scaled)) * (1 if scaled >= 0 else -1)
    elif mode == 'down':
        result = math.floor(abs(scaled)) * (1 if scaled >= 0 else -1)
    else:
        # Excel rounds halves away from zero; Python rounds them to even.
        result = math.floor(abs(scaled) + 0.5) * (1 if scaled >= 0 else -1)
    return result / factor


def _log(args: Sequence[Any]) -> float:
    value = _num(args, 0)
    base = _num(args, 1) if len(args) > 1 else 10.0
    if value <= 0 or base <= 0 or base == 1:
        raise FormulaError(ERRORS['num'])
    return math.log(value, base)


def _safe(function: Callable[[float], float]) -> Callable[[Sequence[Any]], float]:
    def call(args: Sequence[Any]) -> float:
        try:
            return function(_num(args, 0))
        except FormulaError:
            raise
        except (ValueError, OverflowError):
            raise FormulaError(ERRORS['num'])
    return call


def _mid(args: Sequence[Any]) -> str:
    text = _to_text(_one(args, 0))
    start = int(_num(args, 1))
    length = int(_num(args, 2))
    if start < 1 or length < 0:
        raise FormulaError(ERRORS['value'])
    return text[start - 1:start - 1 + length]


FUNCTIONS: Dict[str, Callable[[Sequence[Any]], Any]] = {
    # -- aggregation --
    'SUM': lambda a: float(sum(_numbers_in(_flat(a)))),
    'PRODUCT': lambda a: float(math.prod(_numbers_in(_flat(a)) or [0.0])),
    'AVERAGE': lambda a: statistics.fmean(_guard_empty(_numbers_in(_flat(a)))),
    'MEDIAN': lambda a: statistics.median(_guard_empty(_numbers_in(_flat(a)))),
    'MIN': lambda a: min(_numbers_in(_flat(a)) or [0.0]),
    'MAX': lambda a: max(_numbers_in(_flat(a)) or [0.0]),
    'COUNT': lambda a: float(len(_numbers_in(_flat(a)))),
    'COUNTA': lambda a: float(sum(1 for v in _flat(a) if not _is_blank(v))),
    'COUNTBLANK': lambda a: float(sum(1 for v in _flat(a) if _is_blank(v))),
    'STDEV': lambda a: statistics.stdev(_guard_empty(_numbers_in(_flat(a)))),
    'STDEVP': lambda a: statistics.pstdev(_guard_empty(_numbers_in(_flat(a)))),
    'VAR': lambda a: statistics.variance(_guard_empty(_numbers_in(_flat(a)))),
    'VARP': lambda a: statistics.pvariance(_guard_empty(_numbers_in(_flat(a)))),
    'SUMIF': _sumif,
    'COUNTIF': _countif,
    'SUMPRODUCT': lambda a: float(sum(
        math.prod(pair) for pair in zip(*[
            _numbers_in(arg) if isinstance(arg, _Range) else [_to_number(arg)]
            for arg in a
        ])
    )),
    # -- arithmetic --
    'ABS': _safe(abs),
    'SQRT': _safe(math.sqrt),
    'EXP': _safe(math.exp),
    'LN': _safe(math.log),
    'LOG10': _safe(math.log10),
    'LOG': _log,
    'SIGN': _safe(lambda v: float((v > 0) - (v < 0))),
    'INT': _safe(lambda v: float(math.floor(v))),
    'TRUNC': _safe(lambda v: float(math.trunc(v))),
    'ROUND': lambda a: _round(a, 'half'),
    'ROUNDUP': lambda a: _round(a, 'up'),
    'ROUNDDOWN': lambda a: _round(a, 'down'),
    'POWER': lambda a: math.pow(_num(a, 0), _num(a, 1)),
    'MOD': lambda a: (_num(a, 0) % _num(a, 1)) if _num(a, 1) != 0
                     else _raise(ERRORS['div0']),
    'PI': lambda a: math.pi,
    'SIN': _safe(math.sin), 'COS': _safe(math.cos), 'TAN': _safe(math.tan),
    'ASIN': _safe(math.asin), 'ACOS': _safe(math.acos), 'ATAN': _safe(math.atan),
    'DEGREES': _safe(math.degrees), 'RADIANS': _safe(math.radians),
    # -- logic --
    'IF': _if,
    'IFERROR': _iferror,
    'IFNA': lambda a: (_one(a, 1) if len(a) > 1 else '') if (
        isinstance(_one(a, 0), str) and _one(a, 0) == ERRORS['na']
    ) else _one(a, 0),
    'AND': lambda a: all(_truth(v) for v in _flat(a)),
    'OR': lambda a: any(_truth(v) for v in _flat(a)),
    'NOT': lambda a: not _truth(_one(a, 0)),
    'TRUE': lambda a: True,
    'FALSE': lambda a: False,
    # -- text --
    'CONCAT': lambda a: ''.join(_to_text(v) for v in _flat(a)),
    'CONCATENATE': lambda a: ''.join(_to_text(v) for v in _flat(a)),
    'LEN': lambda a: float(len(_to_text(_one(a, 0)))),
    'LEFT': lambda a: _to_text(_one(a, 0))[:int(_num(a, 1)) if len(a) > 1 else 1],
    'RIGHT': lambda a: _to_text(_one(a, 0))[-(int(_num(a, 1)) if len(a) > 1 else 1):],
    'MID': _mid,
    'TRIM': lambda a: ' '.join(_to_text(_one(a, 0)).split()),
    'UPPER': lambda a: _to_text(_one(a, 0)).upper(),
    'LOWER': lambda a: _to_text(_one(a, 0)).lower(),
    'VALUE': lambda a: _to_number(_one(a, 0)),
    'TEXT': lambda a: _to_text(_one(a, 0)),
}


def _raise(code: str):
    raise FormulaError(code)


# ---------------------------------------------------------------------------
# Sheet-level evaluation
# ---------------------------------------------------------------------------
def evaluate_one(formula: str, lookup: Callable[[int, int], Any]) -> Any:
    """Work out a single formula; *lookup* answers with a cell's value."""
    text = formula[1:] if formula.startswith('=') else formula
    try:
        return _Evaluator(lookup).run(text)
    except FormulaError as exc:
        return exc.code
    except RecursionError:
        return ERRORS['cycle']
    except ZeroDivisionError:
        return ERRORS['div0']
    except Exception:
        return ERRORS['value']


def evaluate_grid(
    formulas: Sequence[Sequence[str]],
    values: Sequence[Sequence[Any]],
) -> Dict[Tuple[int, int], Any]:
    """Work out every formula in a sheet.

    *formulas* holds the formula text per cell ('' where there is none) and
    *values* what the cell shows otherwise.  Results come back keyed by
    (row0, col0); a cell that depends on itself gets #CYCLE! rather than
    hanging, and one that depends on another formula gets that one's result.
    """
    results: Dict[Tuple[int, int], Any] = {}
    running: set = set()

    def formula_at(row: int, col: int) -> str:
        if 0 <= row < len(formulas):
            line = formulas[row]
            if 0 <= col < len(line):
                text = line[col]
                if is_formula(text):
                    return text
        return ''

    def value_at(row: int, col: int) -> Any:
        spot = (row, col)
        if spot in results:
            return results[spot]
        text = formula_at(row, col)
        if text:
            if spot in running:
                results[spot] = ERRORS['cycle']
                return results[spot]
            running.add(spot)
            try:
                results[spot] = evaluate_one(text, value_at)
            finally:
                running.discard(spot)
            return results[spot]
        if 0 <= row < len(values):
            line = values[row]
            if 0 <= col < len(line):
                return line[col]
        return ''

    for row in range(len(formulas)):
        for col in range(len(formulas[row])):
            if formula_at(row, col):
                value_at(row, col)
    return results
