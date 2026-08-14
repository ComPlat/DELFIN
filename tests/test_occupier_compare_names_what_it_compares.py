"""The key says what OCCUPIER compares, not what it executes to compare it.

frequency_calculation_OCCUPIER named the machinery. What it actually decided
was which energy the candidate configurations are ranked by: with frequencies,
Gibbs free energies; without, the electronic energy ORCA prints as FINAL
SINGLE POINT ENERGY. A reader had to know that a frequency run implies a Gibbs
comparison to understand what the switch did.

OCCUPIER_compare=FSPE|G states the decision instead, and the frequency run
becomes its consequence — asking for G is asking for the calculation that
makes G exist.

The distinction is not presentational. Configurations within a few kJ/mol are
ordered differently by electronic energy than by free energy, and since the
Boltzmann weights decide which configurations seed the next redox step, the
choice now also shapes what gets computed downstream.

CONTROL files are kept for years, so the former key keeps working.
"""

import pytest

from delfin.common.control_validator import (
    _as_occupier_compare,
    resolve_occupier_compare,
    validate_control_config,
)


# --------------------------------------------------------------------------
# the values a user might reasonably type
# --------------------------------------------------------------------------

@pytest.mark.parametrize("text", ["FSPE", "fspe", "SPE", "electronic", "e", ""])
def test_the_electronic_energy_spellings(text):
    assert _as_occupier_compare(text) == "FSPE"


@pytest.mark.parametrize("text", ["G", "g", "Gibbs", "gibbs", "free", "free_energy"])
def test_the_free_energy_spellings(text):
    assert _as_occupier_compare(text) == "G"


def test_anything_else_is_refused_rather_than_guessed():
    with pytest.raises(ValueError, match="FSPE or G"):
        _as_occupier_compare("maybe")


def test_the_default_is_the_electronic_energy():
    """Same as the former frequency_calculation_OCCUPIER=no: a Gibbs
    comparison costs a frequency calculation per candidate and is a choice."""
    assert _as_occupier_compare("") == "FSPE"


# --------------------------------------------------------------------------
# which key wins
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "config,expected",
    [
        ({"OCCUPIER_compare": "FSPE"}, "FSPE"),
        ({"OCCUPIER_compare": "G"}, "G"),
        ({"OCCUPIER_compare": "gibbs"}, "G"),
        ({}, "FSPE"),
    ],
)
def test_the_new_key_decides(config, expected):
    assert resolve_occupier_compare(config) == expected


@pytest.mark.parametrize(
    "legacy,expected",
    [("yes", "G"), ("no", "FSPE"), ("YES", "G"), ("true", "G")],
)
def test_a_control_file_written_before_this_still_works(legacy, expected):
    """Runs are reproduced years later from the CONTROL that produced them."""
    assert resolve_occupier_compare({"frequency_calculation_OCCUPIER": legacy}) == expected


def test_agreeing_keys_agree_quietly(caplog):
    assert resolve_occupier_compare(
        {"OCCUPIER_compare": "G", "frequency_calculation_OCCUPIER": "yes"}
    ) == "G"
    assert not [r for r in caplog.records if r.levelname == "WARNING"]


@pytest.mark.parametrize(
    "config,expected",
    [
        ({"OCCUPIER_compare": "G", "frequency_calculation_OCCUPIER": "no"}, "G"),
        ({"OCCUPIER_compare": "FSPE", "frequency_calculation_OCCUPIER": "yes"}, "FSPE"),
    ],
)
def test_contradicting_keys_resolve_to_the_new_one_and_say_so(config, expected, caplog):
    """Picking one silently would leave the report describing a comparison
    that did not happen."""
    import logging

    with caplog.at_level(logging.WARNING):
        assert resolve_occupier_compare(config) == expected

    warnings = [r.getMessage() for r in caplog.records if r.levelname == "WARNING"]
    assert any("disagree" in w for w in warnings), warnings


def test_an_unreadable_new_value_falls_back_rather_than_failing_the_run():
    """A typo in the new key must not lose a setting the old one still carries."""
    assert resolve_occupier_compare(
        {"OCCUPIER_compare": "nonsense", "frequency_calculation_OCCUPIER": "yes"}
    ) == "G"


# --------------------------------------------------------------------------
# the key in a CONTROL file
# --------------------------------------------------------------------------

_BASE = {
    "method": "classic",
    "charge": "0",
    "input_file": "input.txt",
    "functional": "B3LYP",
    "main_basisset": "def2-SVP",
    "smiles_converter": "QUICK",
}


def test_a_control_file_validates_with_the_new_key():
    config = dict(_BASE, OCCUPIER_compare="G")
    out = validate_control_config(config)
    validated = out[0] if isinstance(out, tuple) else out

    assert validated["OCCUPIER_compare"] == "G"


def test_a_control_file_validates_without_it():
    out = validate_control_config(dict(_BASE))
    validated = out[0] if isinstance(out, tuple) else out

    assert validated["OCCUPIER_compare"] == "FSPE"


def test_a_control_file_carrying_only_the_old_key_still_validates():
    config = dict(_BASE, frequency_calculation_OCCUPIER="yes")
    out = validate_control_config(config)
    validated = out[0] if isinstance(out, tuple) else out

    assert validated["frequency_calculation_OCCUPIER"] == "yes"
    assert resolve_occupier_compare(validated) == "G"


def test_a_nonsense_value_stops_the_run_before_it_computes_anything():
    with pytest.raises(ValueError, match="FSPE or G"):
        validate_control_config(dict(_BASE, OCCUPIER_compare="perhaps"))


# --------------------------------------------------------------------------
# what the choice actually switches
# --------------------------------------------------------------------------

def test_asking_for_gibbs_turns_on_the_frequency_run():
    """G exists only after a frequency calculation, so the key that names the
    comparison is what puts FREQ into the candidate inputs."""
    from delfin.config_manager import DelfinConfig

    config = DelfinConfig()
    config.OCCUPIER_compare = "G"
    assert config.is_frequency_calculation_enabled("OCCUPIER") is True

    config.OCCUPIER_compare = "FSPE"
    config.frequency_calculation_OCCUPIER = "no"
    assert config.is_frequency_calculation_enabled("OCCUPIER") is False


def test_the_old_key_still_turns_it_on():
    from delfin.config_manager import DelfinConfig

    config = DelfinConfig()
    config.frequency_calculation_OCCUPIER = "yes"
    assert config.is_frequency_calculation_enabled("OCCUPIER") is True


# --------------------------------------------------------------------------
# an empty field is answered, a filled one is left alone
# --------------------------------------------------------------------------

def test_an_emptied_field_is_told_what_belongs_there():
    """The template ships the key with a value, so an empty one is somebody
    having cleared it. Quietly picking would decide which configuration a run
    calls preferred without anyone having chosen."""
    with pytest.raises(ValueError) as excinfo:
        validate_control_config(dict(_BASE, OCCUPIER_compare=""))

    message = str(excinfo.value)
    assert "set it to FSPE or G" in message
    assert "FINAL SINGLE POINT ENERGY" in message
    assert "Gibbs" in message


@pytest.mark.parametrize("value", ["FSPE", "G", "gibbs"])
def test_a_field_that_says_something_valid_is_left_alone(value, caplog):
    import logging

    with caplog.at_level(logging.INFO):
        out = validate_control_config(dict(_BASE, OCCUPIER_compare=value))
    validated = out[0] if isinstance(out, tuple) else out

    assert validated["OCCUPIER_compare"] in ("FSPE", "G")
    assert "OCCUPIER_compare" not in " ".join(r.getMessage() for r in caplog.records)


@pytest.mark.parametrize("value", ["GG", "ja", "perhaps"])
def test_a_field_that_says_something_impossible_is_reported(value):
    with pytest.raises(ValueError, match="must be FSPE or G"):
        validate_control_config(dict(_BASE, OCCUPIER_compare=value))


def test_a_legacy_file_is_not_told_to_fill_in_a_field_it_predates():
    """The key is absent, not empty. Those files carry the older one and are
    read through it."""
    out = validate_control_config(dict(_BASE, frequency_calculation_OCCUPIER="yes"))
    validated = out[0] if isinstance(out, tuple) else out

    assert validated["OCCUPIER_compare"] == "G"


def test_the_template_ships_the_key_with_a_value():
    """Like frequency_calculation_OCCUPIER=no before it: the default is visible
    and changing it means typing G over FSPE, not knowing the key exists."""
    from delfin import define

    lines = [l for l in define.TEMPLATE.splitlines() if l.startswith("OCCUPIER_compare")]
    assert lines == ["OCCUPIER_compare=FSPE"]
    assert "electronic energy (FINAL SINGLE POINT ENERGY)" not in define.TEMPLATE
