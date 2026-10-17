"""Pure data-sanity helpers for analysis tasks (package D1).

A workspace script imports these before it reports a number. The module
is standalone by design (QS ruling 2026-10, nacht-s30): it imports
nothing from other packages, performs no I/O beyond the values it is
handed, and never guesses where it can refuse.

Failure modes covered:
- a column of values with inconsistent units, summed anyway;
- an ambiguous decimal reading ("8.986" — thousand or decimal?);
- duplicated rows that turn an average into a double-counted one;
- an outlier flagged by vibes instead of a stated rule;
- a summary reported as a bare mean without n or spread.

Adapters: callers that hold mixed unit columns should normalise first
(e.g. convert kJ/mol to kcal/mol with the factor 4.184) and re-check
with :func:`unit_consistency`; this module deliberately does not pick a
conversion for the caller.
"""
from __future__ import annotations

import math
import re
from typing import Any, Iterable, Mapping, Sequence

__all__ = [
    "DataCheckError",
    "unit_consistency",
    "find_duplicates",
    "flag_outliers",
    "summarize",
]


class DataCheckError(ValueError):
    """A check cannot decide — refusal beats a guessed number."""


_UNIT_RE = re.compile(r"^(?P<value>.+?)\s*(?P<unit>[A-Za-zµ/%°][A-Za-z0-9/^%·-]*)$")

# Two readings every bare decimal has: the dot may be a thousands mark
# or the decimal mark. A value is "undecided" only when BOTH readings
# are plausible (>=2 digit groups before a 1-2 digit tail).
_THOUSANDS_RE = re.compile(r"^-?\d{1,3}(\.\d{3})+$")


def _as_float(value: Any) -> float:
    """Reject bool/NaN/inf: they are not data, and each poisons a statistic."""
    if isinstance(value, bool) or value is None:
        raise DataCheckError(
            f"{value!r} is not a numeric value (bool/None)")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise DataCheckError(f"{value!r} is not a numeric value") from exc
    if math.isnan(result) or math.isinf(result):
        raise DataCheckError(
            f"{value!r} is NaN/inf — it poisons every statistic")
    return result


def _parse_number(raw: Any) -> tuple[float | None, str]:
    """Parse one cell into (value, unit). Unit is '' for bare numbers."""
    if raw is None:
        return None, ""
    if isinstance(raw, (int, float)) and not isinstance(raw, bool):
        return float(raw), ""
    text = str(raw).strip()
    if not text or text.lower() in {"n/a", "na", "-", "—", ""}:
        return None, ""
    match = _UNIT_RE.match(text)
    unit = ""
    candidate = text
    if match and not _is_bare_decimal(match.group("value")):
        unit = match.group("unit")
        candidate = match.group("value").strip()
    text = str(raw).strip()
    if "_" in text:
        # '1_234' is both 1234 and 1.234 — refuse, ask.
        return None, ""
    try:
        value = float(candidate)
    except ValueError:
        # Unreadable: not this check's business to guess a repair.
        return None, unit or ""
    return value, unit


def _is_bare_decimal(text: str) -> bool:
    """True when the text is a decimal whose reading is ambiguous."""
    text = text.strip()
    if _THOUSANDS_RE.match(text):
        return True
    return False


def unit_consistency(values: Iterable[Any]) -> dict[str, Any]:
    """Check that a column carries one unit.

    Returns ``{"consistent", "units", "undecided", "unreadable",
    "mixed_shape", "n_bare"}``.
    ``consistent=False`` with ``units=[...]`` names the mixed units so a
    report can say which two readings collided. ``undecided=True`` marks
    an ambiguous bare decimal ("8.986") — the caller must ask, not guess.
    ``mixed_shape=True`` marks a column that mixes unit-carrying values
    with bare numbers: a bare number's unit is an inference, and an
    average over inferred-unit values is the mixed-population trap this
    module exists to catch (QS ruling nacht-s30, 2026-10). ``n_bare`` is
    the count of bare values so the caller can decide; bare values do
    NOT implicitly inherit the column unit.
    """
    units: set[str] = set()
    undecided = False
    unreadable: list[int] = []
    n_bare = 0
    for index, raw in enumerate(values):
        value, unit = _parse_number(raw)
        if value is None:
            if raw not in (None, ""):
                text = str(raw).strip()
                if "_" in text:
                    # '1_234' is both 1234 and 1.234 — refuse, ask.
                    undecided = True
                    continue
                if text.lower() not in {"n/a", "na", "-", "—", ""}:
                    unreadable.append(index)
            continue
        if _THOUSANDS_RE.match(str(raw).strip()):
            # Both readings are plausible (8986 vs 8.986) — refuse, ask.
            undecided = True
            continue
        if unit:
            units.add(unit)
        else:
            n_bare += 1
    if undecided:
        units.add("<undecided>")
    mixed_shape = bool(units) and n_bare > 0
    consistent = (len(units) <= 1 and not undecided and not unreadable
                  and not mixed_shape)
    return {
        "consistent": consistent,
        "units": sorted(units),
        "undecided": undecided,
        "unreadable": unreadable,
        "mixed_shape": mixed_shape,
        "n_bare": n_bare,
    }


def find_duplicates(
    rows: Sequence[Mapping[str, Any]], key: Sequence[str]
) -> dict[str, Any]:
    """Report rows repeated on ``key`` columns, with their row indices.

    Index 0 is the first data row (after the header). A duplicated block
    double-counts in any average — the failure a bare ``drop_duplicates``
    hides because it removes the evidence instead of reporting it.
    """
    if not key:
        raise DataCheckError("key must name at least one column")
    seen: dict[tuple, list[int]] = {}
    for index, row in enumerate(rows):
        missing = [column for column in key if column not in row]
        if missing:
            raise DataCheckError(
                f"key column(s) {missing} not present in row {index}; "
                "a missing key would silently report zero duplicates")
        identity = tuple(row[column] for column in key)
        seen.setdefault(identity, []).append(index)
    groups = [
        {"key": dict(zip(key, identity)), "rows": indices}
        for identity, indices in seen.items()
        if len(indices) > 1
    ]
    groups.sort(key=lambda group: group["rows"][0])
    return {"duplicate_groups": groups, "n_rows": len(rows)}


_OUTLIER_RULES = {"iqr", "zscore"}


def flag_outliers(
    values: Sequence[float], rule: str = "iqr", k: float = 1.5
) -> dict[str, Any]:
    """Flag outliers by a STATED rule. The rule is echoed in the result.

    ``iqr``: outside ``[q1 - k*iqr, q3 + k*iqr]``. ``zscore``:
    ``|z| > k``. Anything else is a refusal — an outlier decision made
    by vibes is not a measurement.
    """
    if rule not in _OUTLIER_RULES:
        raise DataCheckError(
            f"unknown outlier rule {rule!r}; known rules: {sorted(_OUTLIER_RULES)}")
    if not isinstance(k, (int, float)) or isinstance(k, bool) or not k > 0:
        raise DataCheckError(
            f"k must be a positive number, got {k!r}")
    numbers = [_as_float(value) for value in values]
    if len(numbers) < 4:
        raise DataCheckError(
            "outlier rule needs at least 4 values to say anything about")
    if rule == "iqr":
        ordered = sorted(numbers)
        n = len(ordered)

        def _quantile(fraction: float) -> float:
            position = (n - 1) * fraction
            lower = math.floor(position)
            upper = math.ceil(position)
            if lower == upper:
                return ordered[lower]
            weight = position - lower
            return ordered[lower] * (1 - weight) + ordered[upper] * weight

        q1, q3 = _quantile(0.25), _quantile(0.75)
        spread = q3 - q1
        low, high = q1 - k * spread, q3 + k * spread
        indices = [i for i, value in enumerate(numbers) if value < low or value > high]
        bounds = {"low": low, "high": high, "q1": q1, "q3": q3}
    else:
        mean = sum(numbers) / len(numbers)
        std = math.sqrt(sum((value - mean) ** 2 for value in numbers)
                        / (len(numbers) - 1))
        if std == 0:
            indices, bounds = [], {"mean": mean, "std": 0.0}
        else:
            indices = [
                i for i, value in enumerate(numbers)
                if abs(value - mean) / std > k
            ]
            bounds = {"mean": mean, "std": std}
    return {"outlier_indices": indices, "rule": rule, "k": k, "bounds": bounds}


def summarize(values: Sequence[float]) -> dict[str, Any]:
    """Summary with n and spread — never a bare mean.

    ``std`` is the sample standard deviation (ddof=1). An empty input is
    a refusal: a mean of nothing is not zero, it does not exist.
    """
    numbers = [_as_float(value) for value in values]
    n = len(numbers)
    if n == 0:
        raise DataCheckError("cannot summarize an empty column")
    mean = sum(numbers) / n
    ordered = sorted(numbers)
    middle = n // 2
    if n % 2:
        median = ordered[middle]
    else:
        median = (ordered[middle - 1] + ordered[middle]) / 2
    if n > 1:
        std = math.sqrt(
            sum((value - mean) ** 2 for value in numbers) / (n - 1))
    else:
        std = 0.0
    return {
        "n": n,
        "mean": mean,
        "std": std,
        "median": median,
        "min": ordered[0],
        "max": ordered[-1],
    }
