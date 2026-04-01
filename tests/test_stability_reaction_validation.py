import pytest

pytest.importorskip("rdkit")

from delfin.common.control_validator import validate_control_config
from delfin.config import _load_template_defaults
from delfin.stability_constant import (
    STABILITY_REACTION_TEMPLATE,
    validate_stability_reaction_syntax,
)


def _base_config():
    config = _load_template_defaults()
    config["SMILES"] = "O"
    config["solvent"] = "water"
    return config


def test_stability_reaction_template_is_ignored_by_syntax_validator():
    validate_stability_reaction_syntax(STABILITY_REACTION_TEMPLATE, input_smiles="O")


def test_stability_reaction_rejects_unmatched_braces():
    with pytest.raises(ValueError, match="unmatched"):
        validate_stability_reaction_syntax("1*{O}>>>1*{OO", input_smiles="O")


def test_stability_reaction_rejects_atom_imbalance():
    with pytest.raises(ValueError, match="atom-balanced"):
        validate_stability_reaction_syntax("1*{O}>>>1*{OO}", input_smiles="O")


def test_control_validator_allows_stability_reaction_template():
    config = _base_config()
    config["stability_constant_mode"] = "reaction"
    config["stability_reaction"] = STABILITY_REACTION_TEMPLATE

    validated = validate_control_config(config)

    assert validated["stability_reaction"] == STABILITY_REACTION_TEMPLATE
