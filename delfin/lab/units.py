"""Explicit unit table for device setpoints and readings.

Input: a unit symbol as a manifest or a tool call spells it. Output: the
canonical symbol, its dimension, and conversions between symbols of one
dimension. A value whose unit is missing, unknown or of another dimension
is refused (``UnitError``); nothing is guessed.

Rejected alternative: a general unit library. None is a dependency of
DELFIN (pyproject checked 2026-11-05), and adding one for the ~30 symbols
lab devices use would be a silent dependency. The table is closed: a symbol
not listed here is refused, which is the property the safety kernel needs.

Ambiguous spellings are refused on purpose: ``C`` (coulomb or Celsius),
``F`` (farad or Fahrenheit), ``m`` (metre or minute). Temperatures take
``degC``/``°C``, ``K`` or ``degF``.
"""

from __future__ import annotations

from dataclasses import dataclass


class UnitError(ValueError):
    """A unit that is missing, unknown, ambiguous or of the wrong dimension."""


@dataclass(frozen=True)
class Unit:
    symbol: str
    dimension: str
    factor: float          # base = value * factor + offset
    offset: float = 0.0


# Base unit per dimension: degC, rpm, mL/min, W, V, A, Pa, %, s, Hz, mL, g.
_TABLE: dict[str, Unit] = {u.symbol: u for u in (
    Unit("degC", "temperature", 1.0),
    Unit("K", "temperature", 1.0, -273.15),
    Unit("degF", "temperature", 5.0 / 9.0, -32.0 * 5.0 / 9.0),
    Unit("rpm", "rotation", 1.0),
    Unit("mL/min", "flow", 1.0),
    Unit("uL/min", "flow", 1e-3),
    Unit("mL/h", "flow", 1.0 / 60.0),
    Unit("L/min", "flow", 1e3),
    Unit("W", "power", 1.0),
    Unit("kW", "power", 1e3),
    Unit("V", "voltage", 1.0),
    Unit("mV", "voltage", 1e-3),
    Unit("A", "current", 1.0),
    Unit("mA", "current", 1e-3),
    Unit("Pa", "pressure", 1.0),
    Unit("hPa", "pressure", 1e2),
    Unit("kPa", "pressure", 1e3),
    Unit("mbar", "pressure", 1e2),
    Unit("bar", "pressure", 1e5),
    Unit("%", "fraction", 1.0),
    Unit("s", "time", 1.0),
    Unit("min", "time", 60.0),
    Unit("h", "time", 3600.0),
    Unit("Hz", "frequency", 1.0),
    Unit("mL", "volume", 1.0),
    Unit("uL", "volume", 1e-3),
    Unit("L", "volume", 1e3),
    Unit("g", "mass", 1.0),
    Unit("mg", "mass", 1e-3),
    Unit("au", "signal", 1.0),      # arbitrary detector units
)}

_ALIASES = {
    "°C": "degC", "celsius": "degC", "degc": "degC",
    "kelvin": "K",
    "°F": "degF", "degf": "degF",
    "µL/min": "uL/min", "μL/min": "uL/min", "ml/min": "mL/min",
    "ml/h": "mL/h", "l/min": "L/min",
    "µL": "uL", "μL": "uL", "ml": "mL",
    "percent": "%",
}

_AMBIGUOUS = {"C", "F", "m", "M", "c"}


def canonical(symbol: object) -> str:
    """The canonical spelling of ``symbol``; ``UnitError`` when there is none."""
    if symbol is None or (isinstance(symbol, str) and not symbol.strip()):
        raise UnitError("no unit given; a setpoint without a unit is refused")
    if not isinstance(symbol, str):
        raise UnitError(f"unit must be a string, got {type(symbol).__name__}")
    s = symbol.strip()
    if s in _AMBIGUOUS:
        raise UnitError(f"unit {s!r} is ambiguous; spell it out (e.g. degC)")
    if s in _TABLE:
        return s
    alias = _ALIASES.get(s, _ALIASES.get(s.lower(), ""))
    if alias:
        return alias
    raise UnitError(f"unknown unit {s!r}; known: {', '.join(sorted(_TABLE))}")


def dimension(symbol: str) -> str:
    return _TABLE[canonical(symbol)].dimension


def convert(value: float, from_unit: str, to_unit: str) -> float:
    """``value`` in ``from_unit`` expressed in ``to_unit`` (same dimension)."""
    a, b = _TABLE[canonical(from_unit)], _TABLE[canonical(to_unit)]
    if a.dimension != b.dimension:
        raise UnitError(f"{a.symbol} is a {a.dimension}, {b.symbol} is a "
                        f"{b.dimension}; they do not convert")
    base = float(value) * a.factor + a.offset
    return (base - b.offset) / b.factor


def convert_delta(delta: float, from_unit: str, to_unit: str) -> float:
    """A DIFFERENCE in ``from_unit`` expressed in ``to_unit`` (no offset)."""
    a, b = _TABLE[canonical(from_unit)], _TABLE[canonical(to_unit)]
    if a.dimension != b.dimension:
        raise UnitError(f"{a.symbol} and {b.symbol} do not convert")
    return float(delta) * a.factor / b.factor


def known() -> frozenset[str]:
    return frozenset(_TABLE)
