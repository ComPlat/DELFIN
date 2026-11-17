"""control_key_notes.json stays consistent with the CONTROL validator schema.

The key LIST always comes from CONTROL_FIELD_SPECS at runtime
(delfin/common/control_validator.py:2087) -- the notes file adds per-key
knowledge, never a second key list. These tests enforce both directions:

* every note names a key the schema knows (no drift, no typos),
* every note carries the fields an agent needs (meaning, valid_values,
  default, read_by),
* the keys the classic workflows (OCCUPIER, ESD/TADF, IMAG, redox, thermo)
  read are covered -- phase-1 core set, grown as packages OC1/ES1/SH1 send
  their entries.

read_by entries are checked to point into the delfin/ tree so a reader
reference cannot silently rot into prose.
"""
import json
from pathlib import Path

import pytest

from delfin.common.control_validator import CONTROL_FIELD_SPECS

NOTES_PATH = Path(__file__).resolve().parents[1] / "delfin" / "common" / "control_key_notes.json"

REQUIRED_FIELDS = ("meaning", "valid_values", "default", "read_by")

# Keys the classic workflows read; the phase-1 note set must cover these.
CORE_KEYS = (
    "charge",
    "method",
    "functional",
    "main_basisset",
    "main_basisset_rel",
    "metal_basisset",
    "metal_basisset_rel",
    "aux_jk",
    "aux_jk_rel",
    "relativity",
    "initial_guess",
    "geom_opt",
    "freq_type",
    "solvent",
    "implicit_solvation_model",
    "PAL",
    "ri_jkx",
    "xTB_method",
    "XTB_OPT",
    "XTB_GOAT",
    "CREST",
    "OCCUPIER_method",
    "OCCUPIER_compare",
    "IMAG_scope",
    "IMAG_option",
    "IMAG_max_rounds",
    "allow_imaginary_freq",
    "frequency_calculation",
    "states",
    "ESD_modus",
    "TDDFT_nroots",
    "TDDFT_TDA",
)


def _notes() -> dict:
    assert NOTES_PATH.exists(), f"{NOTES_PATH.name} does not exist yet"
    data = json.loads(NOTES_PATH.read_text(encoding="utf-8"))
    keys = data["keys"]
    assert isinstance(keys, dict), "'keys' must be a mapping of key name -> note"
    return keys


def _schema_names() -> set:
    return {spec.name for spec in CONTROL_FIELD_SPECS}


def test_every_noted_key_is_known_to_the_schema():
    schema = _schema_names()
    unknown = sorted(set(_notes()) - schema)
    assert not unknown, (
        f"control_key_notes.json names keys the validator schema does not "
        f"know (typo or a second key list): {unknown}"
    )


def test_every_note_carries_the_required_fields():
    keys = _notes()
    bad = sorted(
        name for name, note in keys.items()
        if any(field not in note or not str(note[field]).strip() for field in REQUIRED_FIELDS)
    )
    assert not bad, f"notes missing required fields {REQUIRED_FIELDS}: {bad}"


def test_core_classic_workflow_keys_have_notes():
    schema = _schema_names()
    missing_in_schema = [name for name in CORE_KEYS if name not in schema]
    assert not missing_in_schema, (
        f"CORE_KEYS lists names the schema does not know: {missing_in_schema}"
    )
    missing_notes = sorted(name for name in CORE_KEYS if name not in _notes())
    assert not missing_notes, (
        f"classic-workflow keys without notes: {missing_notes}"
    )


def test_read_by_references_point_into_the_delfin_tree():
    keys = _notes()
    bad = sorted(
        name for name, note in keys.items()
        if "delfin/" not in str(note.get("read_by", ""))
    )
    assert not bad, f"read_by must cite delfin/ file:line references: {bad}"


def test_note_defaults_match_the_schema_default():
    """A note's default drifts silently otherwise; the schema is the source.

    A note default in parentheses is prose, not a value ("(empty)",
    "(required, no default)") -- the sentence itself says where the value
    really comes from, so it is skipped; a bare value must equal the schema.
    """
    schema_defaults = {spec.name: spec.default for spec in CONTROL_FIELD_SPECS}
    keys = _notes()
    mismatches = []
    for name, note in keys.items():
        note_default = str(note["default"]).strip()
        if note_default.startswith("("):
            continue  # prose sentinel: "(empty ...)", "(required, no default)"
        schema_default = schema_defaults[name]
        try:
            if float(note_default) == float(schema_default):
                continue  # numeric: -50 == -50.0
        except (TypeError, ValueError):
            pass
        if str(schema_default).lower() != note_default.lower():
            mismatches.append(f"{name}: note={note_default!r} schema={schema_default!r}")
    assert not mismatches, f"note defaults disagree with CONTROL_FIELD_SPECS: {mismatches}"


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__]))
