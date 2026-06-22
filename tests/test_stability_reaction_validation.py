import pytest

pytest.importorskip("rdkit")

from delfin.common.control_validator import validate_control_config
from delfin.config import _load_template_defaults, parse_control_text, validate_control_text
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
    config["stability_constant"] = "yes"
    config["stability_constant_mode"] = "reaction"
    config["stability_reaction"] = STABILITY_REACTION_TEMPLATE
    config["thdy_smiles_converter"] = "NORMAL"
    config["thdy_preopt"] = "xtb"

    with pytest.raises(ValueError, match="thermodynamics_reaction / stability_reaction"):
        validate_control_config(config)


def test_thermodynamics_aliases_map_to_stability_keys():
    config = parse_control_text(
        "\n".join(
            [
                "SMILES=O",
                "charge=0",
                "solvent=water",
                "method=classic",
                "smiles_converter=NORMAL",
                "thermodynamics=yes",
                "thermodynamics_mode=reaction",
                f"thermodynamics_reaction={STABILITY_REACTION_TEMPLATE}",
                "thdy_smiles_converter=NORMAL",
                "thdy_preopt=xtb",
            ]
        )
    )

    assert config["stability_constant"] == "yes"
    assert config["stability_constant_mode"] == "reaction"
    assert config["stability_reaction"] == STABILITY_REACTION_TEMPLATE


def test_template_validation_does_not_require_optional_thermodynamics_placeholders():
    errors = validate_control_text(
        "\n".join(
            [
                "charge=[CHARGE]",
                "solvent=[SOLVENT]",
                "method=[classic|manually|OCCUPIER]",
                "smiles_converter=[QUICK|NORMAL|GUPPY|ARCHITECTOR]",
                "thermodynamics=no",
                "thermodynamics_mode=[auto|reaction]",
                "thermodynamics_reaction=a*{SMILES}+b*{SMILES}...>>>c*{SMILES}+d*{SMILES}...",
                "thdy_smiles_converter=[QUICK|NORMAL|GUPPY|ARCHITECTOR]",
                "thdy_preopt=[none|xtb|crest|goat]",
            ]
        )
    )

    joined = " | ".join(errors)
    assert "STABILITY_CONSTANT_MODE" not in joined
    assert "THDY_SMILES_CONVERTER" not in joined
    assert "THDY_PREOPT" not in joined


def test_thermodynamics_yes_requires_explicit_thdy_settings():
    errors = validate_control_text(
        "\n".join(
            [
                "charge=0",
                "solvent=water",
                "method=classic",
                "smiles_converter=NORMAL",
                "thermodynamics=yes",
                "thermodynamics_mode=[auto|reaction]",
                "thdy_smiles_converter=[QUICK|NORMAL|GUPPY|ARCHITECTOR]",
                "thdy_preopt=[none|xtb|crest|goat]",
            ]
        )
    )

    joined = " | ".join(errors)
    assert "thermodynamics_mode / stability_constant_mode" in joined
    assert "thdy_smiles_converter" in joined
    assert "thdy_preopt" in joined
