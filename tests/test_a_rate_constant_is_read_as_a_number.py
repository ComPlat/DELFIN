"""The ESD parser accepts numbers, and only numbers.

Every value in an ESD report — the ISC and IC rate constants, the internal
conversion and fluorescence rates, the temperature, the 0-0 gap, the SOC — was
captured with the character class ``[0-9.+-Ee]``. Read as a human reads it
that is "digits, a dot, a plus, a minus, E or e". Read as a regular expression
it is a *range*: the hyphen sits between ``+`` (0x2B) and ``E`` (0x45), so the
class also accepts ``, - . / : ; < = > ? @ A B C D``.

The failure that follows is silent, which is what makes it worth a test. A
line the parser should decline instead matches, the captured text reaches
``_safe_float``, ``float()`` refuses it, and the function returns ``None``
after a log line nobody reads. The rate constant is then simply absent from
the report — not marked unavailable, not marked wrong, absent — and a
photophysics report missing an ISC rate looks exactly like a calculation that
was never asked for one.

So the class is spelled with the hyphen last, where a hyphen can only stand
for itself, and these tests pin both halves: the numbers ORCA really writes
are still read, and text that is not a number is refused instead of mangled.
"""

from pathlib import Path

import pytest

from delfin import esd_results


# --------------------------------------------------------------------------
# the character class itself
# --------------------------------------------------------------------------

def test_the_number_class_accepts_exactly_what_a_number_is_made_of():
    import re

    accepted = {chr(c) for c in range(32, 127)
                if re.fullmatch(esd_results._NUM, chr(c))}
    assert accepted == set("0123456789.eE+-")


@pytest.mark.parametrize("char", list(",/:;<=>?@ABCD"))
def test_the_characters_the_old_range_let_through_are_refused(char):
    """These are the characters between "+" and "E" in the ASCII table. Each
    one was accepted as part of a number before the class was rewritten."""
    import re

    assert not re.fullmatch(esd_results._NUM, char)


# --------------------------------------------------------------------------
# numbers ORCA actually prints
# --------------------------------------------------------------------------

_ISC_OUTPUT = """\
                         ORCA ESD MODULE
Temperature used: 298.15 K
0-0 energy difference: 4321.5 cm-1
Reference SOC (Re and Im): -12.345, 0.678
The calculated ISC rate constant is 2.5E+10 s-1
with 0.93 from FC and 0.07 from HT
"""


def test_a_real_isc_output_is_read_completely(tmp_path: Path):
    out = tmp_path / "ISC.out"
    out.write_text(_ISC_OUTPUT)

    result = esd_results._parse_isc_output(out)

    assert result.rate == pytest.approx(2.5e10)
    assert result.temperature == pytest.approx(298.15)
    assert result.delta_cm1 == pytest.approx(4321.5)


@pytest.mark.parametrize(
    "printed,value",
    [
        ("298.15", 298.15),
        ("0", 0.0),
        ("-1", -1.0),
        ("2.5E+10", 2.5e10),
        ("2.5e-10", 2.5e-10),
        ("1.0E10", 1.0e10),
        ("-3.25E-07", -3.25e-07),
    ],
)
def test_the_shapes_orca_writes_a_number_in(printed, value, tmp_path: Path):
    out = tmp_path / "ISC.out"
    out.write_text(f"The calculated ISC rate constant is {printed} s-1\n")

    assert esd_results._parse_isc_output(out).rate == pytest.approx(value)


# --------------------------------------------------------------------------
# text that is not a number
# --------------------------------------------------------------------------

def test_a_temperature_that_is_not_a_number_is_declined_not_captured(tmp_path: Path):
    """"Temperature used: ABC K" matched the old class, and ABC then became a
    missing temperature rather than an unmatched line."""
    out = tmp_path / "ISC.out"
    out.write_text("Temperature used: ABC K\n")

    assert not esd_results._TEMP_RE.search(out.read_text())


def test_a_number_stops_at_the_punctuation_that_follows_it():
    """The SOC pair is printed "Re, Im". With ":" and "," inside the class the
    first capture could run past its own value."""
    match = esd_results._SOC_RE.search(
        "Reference SOC (Re and Im): -12.345, 0.678"
    )

    assert match.group(1) == "-12.345"
    assert match.group(2) == "0.678"


def test_a_line_that_only_looks_like_a_result_is_not_mistaken_for_one(tmp_path: Path):
    """ORCA prints table rules and "N/A" markers. Under the old class "@" and
    the capital letters were numbers, so a rule could satisfy the pattern."""
    out = tmp_path / "ISC.out"
    out.write_text(
        "The calculated ISC rate constant is N/A s-1\n"
        "Temperature used: <not set> K\n"
    )

    result = esd_results._parse_isc_output(out)

    assert result.rate is None
    assert result.temperature is None


def test_nothing_is_silently_dropped_when_the_value_is_readable(tmp_path: Path, caplog):
    """The old class turned a parse failure into a log warning and a missing
    field. A good value must produce neither."""
    import logging

    out = tmp_path / "ISC.out"
    out.write_text(_ISC_OUTPUT)

    with caplog.at_level(logging.WARNING):
        result = esd_results._parse_isc_output(out)

    assert result.rate is not None
    assert not [r for r in caplog.records if "Failed to convert" in r.getMessage()]
