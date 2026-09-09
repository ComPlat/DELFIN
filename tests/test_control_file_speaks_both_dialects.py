"""The CONTROL file renamed two settings; both spellings have to keep working.

`XTB_OPT` is written `XTB_preOPT`, and the pair `XTB_GOAT`/`CREST` collapsed
into a single `global_optimizer`. Files written before that must parse to
exactly what they parsed to before, and files written after must reach the
workflow code, which still reads the older names — by bare index, so they can
never be missing.
"""

import pytest

from delfin import define
from delfin.common.control_validator import validate_control_config
from delfin.config import (
    parse_control_text,
    read_control_file,
    set_control_value,
    validate_control_text,
)

BASE = "charge=0\nsolvent=water\nmethod=classic\nsmiles_converter=NORMAL\n"


def _parse(extra: str) -> dict:
    return parse_control_text(BASE + extra)


def test_old_spelling_still_selects_goat():
    config = _parse("XTB_OPT=yes\nXTB_GOAT=yes\nCREST=no\n")
    assert config["XTB_OPT"] == "yes"
    assert config["XTB_GOAT"] == "yes"
    assert config["CREST"] == "no"
    # derived so new-style readers see a value too
    assert config["global_optimizer"] == "GOAT"


def test_old_spelling_still_selects_crest():
    config = _parse("XTB_GOAT=no\nCREST=yes\n")
    assert config["CREST"] == "yes"
    assert config["global_optimizer"] == "CREST"


def test_new_spelling_reaches_the_old_keys():
    config = _parse("XTB_preOPT=yes\nglobal_optimizer=GOAT\n")
    assert config["XTB_OPT"] == "yes"
    assert config["XTB_GOAT"] == "yes"
    assert config["CREST"] == "no"


def test_new_spelling_selects_crest():
    config = _parse("XTB_preOPT=no\nglobal_optimizer=CREST\n")
    assert config["CREST"] == "yes"
    assert config["XTB_GOAT"] == "no"


@pytest.mark.parametrize(
    "extra",
    [
        "",                                    # neither key present at all
        "global_optimizer=\n",                 # present but emptied
        "global_optimizer=[GOAT|CREST]\n",     # template placeholder left alone
        "global_optimizer=none\n",
    ],
)
def test_no_choice_means_no_global_search(extra):
    """The workflow indexes these three without a default; they must exist."""
    config = _parse(extra)
    assert config["XTB_OPT"] == "no"
    assert config["XTB_GOAT"] == "no"
    assert config["CREST"] == "no"


def test_the_new_key_wins_over_a_contradicting_old_one():
    config = _parse("XTB_GOAT=yes\nglobal_optimizer=CREST\n")
    assert config["CREST"] == "yes"
    assert config["XTB_GOAT"] == "no"


def test_untouched_placeholder_is_not_an_error():
    assert validate_control_text(BASE + "global_optimizer=[GOAT|CREST]\n") == []


def test_a_wrong_optimizer_says_what_belongs_there():
    errors = validate_control_text(BASE + "global_optimizer=Quatsch\n")
    assert len(errors) == 1
    assert "GOAT" in errors[0] and "CREST" in errors[0] and "none" in errors[0]


def test_yes_no_spellings_are_normalised():
    """`XTB_OPT=True` used to read as "off" at every `== "yes"` call site."""
    validated = validate_control_config(_parse("XTB_OPT=True\nXTB_GOAT=1\n"))
    assert validated["XTB_OPT"] == "yes"
    assert validated["XTB_GOAT"] == "yes"


def test_the_experimental_block_is_optional_but_still_read(tmp_path):
    control = tmp_path / "CONTROL.txt"
    control.write_text(
        BASE
        + "Literature_reference=Some Paper 2019\n"
        + "reference_CV=V Vs. SCE\n"
        + "E_red_exp=-1.85\n"
        + "*E_ox_exp=0.42\n",
        encoding="utf-8",
    )
    config = read_control_file(str(control))
    assert config["Literature_reference"] == "Some Paper 2019"
    assert config["reference_CV"] == "V Vs. SCE"
    assert str(config["E_red_exp"]) == "-1.85"
    assert str(config["*E_ox_exp"]) == "0.42"

    # ...and it is gone from the shipped template, so a file that never had it
    # parses without inventing the keys.
    assert "Literature_reference" not in define.TEMPLATE
    plain = tmp_path / "plain.txt"
    plain.write_text(BASE, encoding="utf-8")
    assert "E_red_exp" not in read_control_file(str(plain))


def test_the_shipped_template_only_complains_about_its_placeholders():
    """A round trip the template did not survive before: it shipped
    `OCCUPIER_method=auto|manually`, which the validator rejected."""
    errors = validate_control_text(define.TEMPLATE)
    joined = " ".join(errors)
    for placeholder in ("[CHARGE]", "[SOLVENT]", "[METHOD]"):
        assert placeholder in joined
    # nothing to convert yet, so the converter is not demanded
    assert len(errors) == 3
    assert "[SMILES_CONVERTER]" not in joined

    with_smiles = validate_control_text(
        define.TEMPLATE.replace("SMILES=", "SMILES=c1ccccc1")
    )
    assert len(with_smiles) == 4
    assert "[SMILES_CONVERTER]" in " ".join(with_smiles)


def test_the_dashboard_shows_the_template_it_validates_against():
    from delfin.dashboard import constants

    assert constants.DEFAULT_CONTROL is define.TEMPLATE
    assert "\nglobal_optimizer=GOAT\n" in constants.ONLY_GOAT_TEMPLATE
    assert "\ncalc_initial=no\n" in constants.ONLY_GOAT_TEMPLATE
    assert "\nmethod=classic\n" in constants.ONLY_GOAT_TEMPLATE


def test_the_deliberate_resource_change_is_the_only_one(tmp_path):
    """PAL and pal_jobs were raised on purpose; pin them so a later template
    edit cannot move them silently."""
    control = tmp_path / "CONTROL.txt"
    control.write_text(
        define.TEMPLATE.replace("[CHARGE]", "0")
        .replace("[SOLVENT]", "acetonitrile")
        .replace("[QUICK|NORMAL|GUPPY|ARCHITECTOR]", "NORMAL")
        .replace("[classic|manually|OCCUPIER]", "classic"),
        encoding="utf-8",
    )
    config = read_control_file(str(control))
    assert config["PAL"] == 48
    assert config["pal_jobs"] == 4
    assert config["orca_parallel_strategy"] == "auto"


def test_dropping_the_guppy_lines_did_not_drop_their_values(tmp_path):
    """Four GUPPY keys left the template; their defaults have to survive,
    because they used to reach the run through the template."""
    control = tmp_path / "CONTROL.txt"
    control.write_text(BASE, encoding="utf-8")
    config = read_control_file(str(control))
    assert config["GUPPY_START_STRATEGY"] == "isomers"
    assert config["GUPPY_MAX_ISOMERS"] == 100
    assert config["GUPPY_RMSD_CUTOFF"] == 0.3
    assert config["GUPPY_ENERGY_WINDOW_KCAL"] == 25.0


def test_occupier_selection_resolves_the_same_as_the_old_pipe_list():
    """The template used to ship `tolerance|truncation|rounding`; every reader
    split on the pipe, so shipping plain `tolerance` must decide identically."""
    from delfin.reporting.occupier_selection import OccupierSelector

    def selector(value):
        return OccupierSelector({"occupier_selection": value}, sequence=[], dev_cache={})

    old = selector("tolerance|truncation|rounding")
    new = selector("tolerance")
    assert old.method == new.method == "tolerance"


def _auto_control(tmp_path, **overrides):
    text = (
        define.TEMPLATE.replace("[CHARGE]", "0")
        .replace("[SOLVENT]", "acetonitrile")
        .replace("[QUICK|NORMAL|GUPPY|ARCHITECTOR]", "NORMAL")
        .replace("[classic|manually|OCCUPIER]", "OCCUPIER")
    )
    for key, value in overrides.items():
        text = set_control_value(text, key, value)
    path = tmp_path / "CONTROL.txt"
    path.write_text(text, encoding="utf-8")
    return path


def test_stripping_the_sequences_keeps_what_follows_them(tmp_path):
    """The section used to be cut away up to the `INFOS:` documentation block.
    Without one, everything behind the sequences went too — and on the copy
    path that truncation was written back to disk."""
    from delfin.occupier_sequences import remove_existing_sequence_blocks

    control = _auto_control(tmp_path)
    before = control.read_text(encoding="utf-8").splitlines()
    stripped = remove_existing_sequence_blocks(control, persist=False).splitlines()

    assert "OCCUPIER_sequence_profiles:" not in "\n".join(stripped)
    for survivor in (
        "keyword:basename=[]",
        "co2_coordination=off",
        "functional=PBE0",
        "PAL=48",
        "thdy_preopt=[none|xtb|crest|goat]",
    ):
        assert survivor in stripped, f"{survivor!r} was cut away with the sequences"
    # only the sequence section goes
    assert 0 < len(before) - len(stripped) < 40


def test_an_occupier_auto_run_still_reads_the_keys_behind_the_sequences(tmp_path):
    control = _auto_control(tmp_path, co2_coordination="on", maxcore=1234)
    config = read_control_file(str(control))
    assert config["co2_coordination"] == "on"
    assert config["maxcore"] == 1234


# --- smiles_converter is only required when there is a SMILES to convert ------

def test_an_xyz_run_does_not_have_to_pick_a_smiles_converter(tmp_path):
    """Every shipped example CONTROL.txt was rejected over this: they feed an
    XYZ block, so there is nothing for a converter to do."""
    control = tmp_path / "CONTROL.txt"
    control.write_text("charge=0\nsolvent=water\nmethod=classic\n", encoding="utf-8")
    (tmp_path / "input.txt").write_text("C 0.0 0.0 0.0\nH 0.0 0.0 1.1\n", encoding="utf-8")

    assert validate_control_text(control.read_text()) == []
    assert read_control_file(str(control))["smiles_converter"] == "NORMAL"


def test_a_smiles_run_still_has_to_pick_one(tmp_path):
    control = tmp_path / "CONTROL.txt"
    control.write_text(
        "charge=0\nsolvent=water\nmethod=classic\nSMILES=c1ccccc1\n", encoding="utf-8"
    )
    errors = validate_control_text(control.read_text())
    assert len(errors) == 1
    assert "SMILES" in errors[0]
    for option in ("QUICK", "NORMAL", "MANTA", "ARCHITECTOR"):
        assert option in errors[0]
    with pytest.raises(ValueError, match="smiles_converter"):
        read_control_file(str(control))


def test_a_smiles_in_the_input_file_counts_too(tmp_path):
    """The conversion is triggered by the input file's content, not by the
    SMILES key, so the input file has to be able to demand a converter."""
    control = tmp_path / "CONTROL.txt"
    control.write_text("charge=0\nsolvent=water\nmethod=classic\n", encoding="utf-8")
    (tmp_path / "input.txt").write_text("c1ccccc1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="smiles_converter"):
        read_control_file(str(control))


def test_a_caller_that_knows_can_say_so():
    """The Submit tab validates before it writes its input box into the
    CONTROL text, so it answers for it."""
    text = "charge=0\nsolvent=water\nmethod=classic\n"
    assert validate_control_text(text, converts_smiles=True) != []
    assert validate_control_text(text + "SMILES=c1ccccc1\n", converts_smiles=False) == []


def test_the_shipped_examples_parse(tmp_path):
    """The examples are the oldest real CONTROL files there are."""
    import glob
    from pathlib import Path as _Path

    root = _Path(__file__).resolve().parent.parent / "examples"
    files = sorted(glob.glob(str(root / "**" / "CONTROL.txt"), recursive=True))
    assert files, "no example CONTROL.txt found"
    unreadable = {}
    for path in files:
        try:
            config = read_control_file(path)
        except ValueError as exc:
            unreadable[path] = str(exc)
            continue
        assert config["XTB_OPT"] in {"yes", "no"}
        assert config["global_optimizer"] in {"", "GOAT", "CREST"}
    # ZnTpy_Me leaves `method=` empty, which is a real gap in that file
    assert all("method" in msg for msg in unreadable.values()), unreadable
    assert len(unreadable) <= 1, unreadable


def test_the_old_guppy_spelling_still_selects_the_builder(tmp_path):
    """Three generations of CONTROL file, one builder.

    The oldest archived runs say ``GUPPY=yes`` with no converter key at all;
    the next say ``smiles_converter=GUPPY``; the current say ``MANTA``.  All
    three select the same thing, and 62 of the 126 archived GUPPY runs are the
    first kind -- if the bridge goes, their CONTROL files stop parsing.
    """
    for body in (
        "charge=0\nsolvent=water\nmethod=classic\nSMILES=c1ccccc1\nGUPPY=yes\n",
        "charge=0\nsolvent=water\nmethod=classic\nSMILES=c1ccccc1\n"
        "smiles_converter=GUPPY\n",
        "charge=0\nsolvent=water\nmethod=classic\nSMILES=c1ccccc1\n"
        "smiles_converter=MANTA\n",
    ):
        control = tmp_path / "CONTROL.txt"
        control.write_text(body, encoding="utf-8")
        config = read_control_file(str(control))
        assert config.get("smiles_converter") == "MANTA", body
