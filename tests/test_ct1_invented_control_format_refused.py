"""An invented CONTROL format must be refused with a message that names the keys.

Wave-24 evidence (operator): OC1's builder wrote a CONTROL.txt in an invented
format ('Charge 0', 'Occupier auto', 'Functional r2SCAN-3c') and DELFIN's
parser refused it with 'Missing required CONTROL values for: charge, method,
solvent'. The refusal is correct behaviour -- the gap is that the agent's
CONTROL path (tune-control.md) must START from DELFIN's template and use only
schema keys. This test pins both: the parser refuses the invented file AND
its error names the missing keys, so the agent's recovery is guided.
"""

import pytest

from delfin.config import read_control_file


def _write_invented_format(tmp_path):
    control = tmp_path / "CONTROL.txt"
    control.write_text(
        "Charge 0\n"
        "Occupier auto\n"
        "Functional r2SCAN-3c\n"
        "Solvent none\n",
        encoding="utf-8",
    )
    return control


def test_an_invented_format_is_refused_by_the_real_parser(tmp_path):
    """Space-separated 'Key value' lines are not KEY=value; the required keys
    stay unanswered and the parser raises (config.py:457-462)."""
    control = _write_invented_format(tmp_path)
    with pytest.raises(ValueError) as excinfo:
        read_control_file(str(control))
    message = str(excinfo.value)
    assert "Missing required CONTROL values for:" in message
    for key in ("charge", "method", "solvent"):
        assert key in message, (key, message)


def test_the_error_message_tells_the_fix(tmp_path):
    """The refusal guides recovery: it names the placeholder substitution."""
    control = _write_invented_format(tmp_path)
    with pytest.raises(ValueError) as excinfo:
        read_control_file(str(control))
    assert "Replace template placeholders" in str(excinfo.value)


def test_the_same_content_in_template_format_parses(tmp_path):
    """The control: the SAME values written as KEY=value parse fine -- the
    refusal is about the shape, not the chemistry."""
    control = tmp_path / "CONTROL.txt"
    control.write_text(
        "method=classic\n"
        "charge=0\n"
        "input=geometry.xyz\n"
        "solvent=acetonitrile\n"
        "implicit_solvation_model=CPCM\n"
        "functional=r2SCAN-3c\n",
        encoding="utf-8",
    )
    config = read_control_file(str(control))
    assert config["method"] == "classic"
    assert str(config["charge"]) == "0"
    assert config["solvent"] == "acetonitrile"
