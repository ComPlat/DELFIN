from delfin.config import validate_control_text


def test_control_validator_allows_optts_geom_opt():
    errors = validate_control_text(
        "\n".join(
            [
                "charge=0",
                "solvent=water",
                "method=classic",
                "smiles_converter=NORMAL",
                "geom_opt=OPTTS",
            ]
        )
    )

    assert not [error for error in errors if "geom_opt" in error]


def test_control_validator_allows_optts_with_additional_keywords():
    errors = validate_control_text(
        "\n".join(
            [
                "charge=0",
                "solvent=water",
                "method=classic",
                "smiles_converter=NORMAL",
                "geom_opt=OPTTS TightSCF",
            ]
        )
    )

    assert not [error for error in errors if "geom_opt" in error]
