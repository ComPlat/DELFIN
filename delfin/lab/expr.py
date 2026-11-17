"""Interlock expressions: a closed grammar with three-valued logic.

Input: text such as ``hotplate.setpoint > 40 degC`` or ``stirrer.power ==
on and fume_hood.power == on``. Output: a parsed tree whose evaluation
against a lookup returns True, False or UNKNOWN.

Grammar (no ``eval``, no attribute access, no calls)::

    expr    := or
    or      := and ("or" and)*
    and     := not ("and" not)*
    not     := "not" not | cmp | "(" expr ")"
    cmp     := ref OP operand
    ref     := ident "." ident
    operand := number [unit] | word | true | false
    OP      := == != < <= > >=

Three-valued (Kleene) logic: a reference whose reading is stale, missing
or unreadable is UNKNOWN; ``False and UNKNOWN`` is False, ``True or
UNKNOWN`` is True, everything else with an UNKNOWN in it is UNKNOWN. The
kernel accepts an interlock only when its requirement is True, so UNKNOWN
fails closed.

A number may carry a unit (``40 degC``, ``313 K``); it is converted into the
referenced channel's unit before comparing. A number without a unit is
taken in the channel's own unit, which the manifest declares.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Callable, Optional, Union

from . import units as _units


class ExprError(ValueError):
    """An expression that does not parse."""


class _Unknown:
    _inst = None

    def __new__(cls):
        if cls._inst is None:
            cls._inst = super().__new__(cls)
        return cls._inst

    def __repr__(self) -> str:
        return "UNKNOWN"

    def __bool__(self) -> bool:          # never silently truthy
        raise TypeError("UNKNOWN has no truth value; compare with `is True`")


UNKNOWN = _Unknown()
Tri = Union[bool, _Unknown]

_TOKEN = re.compile(r"""
    \s*(?:
      (?P<num>-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)
    | (?P<op>==|!=|<=|>=|<|>)
    | (?P<lp>\()
    | (?P<rp>\))
    | (?P<word>[A-Za-z_°µμ%][A-Za-z0-9_°µμ%/]*(?:\.[A-Za-z_][A-Za-z0-9_]*)?)
    )""", re.VERBOSE)

_KEYWORDS = {"and", "or", "not", "true", "false"}


def _tokens(text: str) -> list[tuple[str, str]]:
    out, pos = [], 0
    text = text or ""
    while pos < len(text):
        if text[pos:].strip() == "":
            break
        m = _TOKEN.match(text, pos)
        if not m or m.end() == pos:
            raise ExprError(f"cannot read {text[pos:pos + 12]!r} in {text!r}")
        kind = m.lastgroup
        out.append((kind, m.group(kind)))
        pos = m.end()
    return out


@dataclass(frozen=True)
class Ref:
    device: str
    channel: str

    @property
    def name(self) -> str:
        return f"{self.device}.{self.channel}"


@dataclass(frozen=True)
class Cmp:
    ref: Ref
    op: str
    value: object                 # float or str or bool
    unit: Optional[str] = None


@dataclass(frozen=True)
class Not:
    arg: object


@dataclass(frozen=True)
class And:
    args: tuple


@dataclass(frozen=True)
class Or:
    args: tuple


class _Parser:
    def __init__(self, text: str):
        self.text = text
        self.toks = _tokens(text)
        self.i = 0

    def peek(self):
        return self.toks[self.i] if self.i < len(self.toks) else ("eof", "")

    def take(self):
        tok = self.peek()
        self.i += 1
        return tok

    def parse(self):
        node = self.or_()
        if self.peek()[0] != "eof":
            raise ExprError(f"unexpected {self.peek()[1]!r} in {self.text!r}")
        return node

    def or_(self):
        args = [self.and_()]
        while self.peek() == ("word", "or"):
            self.take()
            args.append(self.and_())
        return args[0] if len(args) == 1 else Or(tuple(args))

    def and_(self):
        args = [self.not_()]
        while self.peek() == ("word", "and"):
            self.take()
            args.append(self.not_())
        return args[0] if len(args) == 1 else And(tuple(args))

    def not_(self):
        kind, val = self.peek()
        if (kind, val) == ("word", "not"):
            self.take()
            return Not(self.not_())
        if kind == "lp":
            self.take()
            node = self.or_()
            if self.take()[0] != "rp":
                raise ExprError(f"missing ')' in {self.text!r}")
            return node
        return self.cmp()

    def cmp(self):
        kind, val = self.take()
        if kind != "word" or "." not in val or val in _KEYWORDS:
            raise ExprError(f"expected device.channel, got {val!r} in "
                            f"{self.text!r}")
        dev, ch = val.split(".", 1)
        kind, op = self.take()
        if kind != "op":
            raise ExprError(f"expected a comparison after {dev}.{ch} in "
                            f"{self.text!r}")
        kind, raw = self.take()
        unit = None
        if kind == "num":
            value: object = float(raw)
            nk, nv = self.peek()
            if nk == "word" and nv not in _KEYWORDS and "." not in nv:
                self.take()
                unit = _units.canonical(nv)
        elif kind == "word" and raw in ("true", "false"):
            value = raw == "true"
        elif kind == "word" and "." not in raw and raw not in _KEYWORDS:
            value = raw
            if op not in ("==", "!="):
                raise ExprError(f"a state compares with == or != only: "
                                f"{self.text!r}")
        else:
            raise ExprError(f"expected a value after {op!r} in {self.text!r}")
        return Cmp(Ref(dev, ch), op, value, unit)


def parse(text: str):
    """The tree of ``text``; ``ExprError`` when it does not parse."""
    if not isinstance(text, str) or not text.strip():
        raise ExprError("empty expression")
    return _Parser(text).parse()


def refs(node) -> set[Ref]:
    if isinstance(node, Cmp):
        return {node.ref}
    if isinstance(node, Not):
        return refs(node.arg)
    if isinstance(node, (And, Or)):
        out: set[Ref] = set()
        for a in node.args:
            out |= refs(a)
        return out
    return set()


#: lookup(ref) -> (value, unit) or UNKNOWN
Lookup = Callable[[Ref], object]


def _compare(a, op: str, b) -> bool:
    if op == "==":
        return a == b
    if op == "!=":
        return a != b
    if not (isinstance(a, (int, float)) and isinstance(b, (int, float))):
        raise TypeError("ordering needs numbers")
    return {"<": a < b, "<=": a <= b, ">": a > b, ">=": a >= b}[op]


def evaluate(node, lookup: Lookup) -> Tri:
    if isinstance(node, Cmp):
        got = lookup(node.ref)
        if got is UNKNOWN:
            return UNKNOWN
        value, unit = got
        want = node.value
        try:
            if isinstance(want, float):
                if isinstance(value, bool) or not isinstance(value, (int, float)):
                    return UNKNOWN
                if node.unit and unit:
                    want = _units.convert(want, node.unit, unit)
                elif node.unit and not unit:
                    return UNKNOWN
            elif isinstance(want, bool):
                if not isinstance(value, bool):
                    return UNKNOWN
            else:
                value = str(value)
            return _compare(value, node.op, want)
        except (TypeError, _units.UnitError):
            return UNKNOWN
    if isinstance(node, Not):
        v = evaluate(node.arg, lookup)
        return UNKNOWN if v is UNKNOWN else (not v)
    if isinstance(node, And):
        vals = [evaluate(a, lookup) for a in node.args]
        if any(v is False for v in vals):
            return False
        return UNKNOWN if any(v is UNKNOWN for v in vals) else True
    if isinstance(node, Or):
        vals = [evaluate(a, lookup) for a in node.args]
        if any(v is True for v in vals):
            return True
        return UNKNOWN if any(v is UNKNOWN for v in vals) else False
    raise ExprError(f"not an expression node: {node!r}")
