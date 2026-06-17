"""Tests for the CONTROL.txt-style application keys file + delfin-app CLI."""

from __future__ import annotations

from delfin.tools import _keyfile


def test_application_keyfile_has_keys_and_allowed_values():
    text = _keyfile.application_keyfile("redox_potential")
    assert "application=redox_potential" in text
    # required keys present
    assert "\nsmiles=" in text
    assert "\ncharge=" in text
    # method carries a default + allowed-values hint from the key vocabulary
    assert "method=B3LYP" in text
    assert "allowed:" in text and "more)" in text   # large enum is summarised


def test_parse_keyfile_skips_comments_and_keeps_hash_in_values():
    text = (
        "# a comment\n"
        "application=redox_potential\n"
        "smiles=C#N\n"            # SMILES triple bond — '#' must survive
        "charge=-1\n"
        "method=\n"              # empty → dropped (default applies)
        "\n"
    )
    app_name, values = _keyfile.parse_keyfile(text)
    assert app_name == "redox_potential"
    assert values["smiles"] == "C#N"
    assert values["charge"] == "-1"
    assert "method" not in values


def test_inputs_from_keyfile_coerces_types():
    text = "application=redox_potential\nsmiles=C#N\ncharge=-1\nmult_ox=2\n"
    app_name, inputs = _keyfile.inputs_from_keyfile(text)
    assert app_name == "redox_potential"
    assert inputs == {"smiles": "C#N", "charge": -1, "mult_ox": 2}   # int-coerced


def test_inputs_from_keyfile_roundtrips_a_generated_template():
    template = _keyfile.application_keyfile("multi_level_energy")
    # fill the two required keys, leave the rest at their defaults
    filled = template.replace("smiles=\n", "smiles=CCO\n").replace("charge=\n", "charge=0\n")
    app_name, inputs = _keyfile.inputs_from_keyfile(filled)
    assert app_name == "multi_level_energy"
    assert inputs["smiles"] == "CCO" and inputs["charge"] == 0


def test_run_from_keyfile_calls_platform_with_typed_inputs(tmp_path, monkeypatch):
    from delfin.tools import platform

    captured = {}

    def _fake_run(name, *, cores=1, **inputs):
        captured.update(name=name, cores=cores, inputs=inputs)
        return "RESULT"

    monkeypatch.setattr(platform, "run_application", _fake_run)

    kf = tmp_path / "r.txt"
    kf.write_text("application=redox_potential\nsmiles=C#N\ncharge=-1\nmult_ox=2\n")
    result = _keyfile.run_from_keyfile(str(kf), cores=4)

    assert result == "RESULT"
    assert captured["name"] == "redox_potential" and captured["cores"] == 4
    assert captured["inputs"] == {"smiles": "C#N", "charge": -1, "mult_ox": 2}


def test_unknown_application_keyfile_raises():
    import pytest

    with pytest.raises(ValueError):
        _keyfile.application_keyfile("does_not_exist")


# --- CLI ------------------------------------------------------------------


def test_cli_list_and_template(capsys):
    from delfin.cli_app import main

    assert main(["list"]) == 0
    out = capsys.readouterr().out
    assert "redox_potential" in out

    assert main(["template", "redox_potential"]) == 0
    out = capsys.readouterr().out
    assert "application=redox_potential" in out


def test_cli_run_reports_outputs(tmp_path, monkeypatch, capsys):
    from delfin.cli_app import main
    from delfin.tools import platform

    class _Res:
        ok = True
        outputs = {"e_neutral_Eh": -1.0}
        error = None

    monkeypatch.setattr(platform, "run_application", lambda name, **k: _Res())
    kf = tmp_path / "r.txt"
    kf.write_text("application=redox_potential\nsmiles=CCO\ncharge=0\n")
    assert main(["run", str(kf)]) == 0
    out = capsys.readouterr().out
    assert "OK" in out and "e_neutral_Eh" in out
