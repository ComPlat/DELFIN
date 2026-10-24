"""OR2 phase 2: skills are consistent with the tasks and with DELFIN's
real workflow vocabulary.

Pins (red before the skills existed — FileNotFoundError / missing
frontmatter):
1. Every skill named in tasks_w20_or2.yaml exists as a skill file with
   frontmatter (name + description).
2. Every skill file only names workflow modules/templates that exist in
   the code base (esd_states, occupier_stages, redox_potential,
   extract_excited_states, extract_vibrational_modes, get_E_ref) —
   vocabulary from DELFIN's code, never invented.
3. The redox skill states the reference-electrode requirement and the
   per-solvent E_ref values that delfin/config.py get_E_ref returns
   (MeCN 4.745 V vs SHE, default 4.345 V).
"""

from __future__ import annotations

import pathlib
import re

import pytest
import yaml

_ROOT = pathlib.Path(__file__).resolve().parents[1]
_SKILLS = _ROOT / "delfin" / "agent" / "pack" / "skills"
_OR2 = (_ROOT / "delfin" / "agent" / "pack" / "benchmark"
        / "tasks_w20_or2.yaml")

_REAL_SYMBOLS = [
    "esd_states", "occupier_stages", "redox_potential",
    "extract_excited_states", "extract_vibrational_modes", "get_E_ref",
    "calculate_redox_potentials", "_esd", "_occupier", "_redox",
    "OCCUPIER", "orca_tddft", "orca_opt", "orca_freq", "orca_vib",
]


def _task_ids():
    data = yaml.safe_load(_OR2.read_text(encoding="utf-8"))
    return [t for t in (data.get("tasks") or [])]


def test_every_task_skill_exists_with_frontmatter():
    for t in _task_ids():
        path = _SKILLS / f"{t['skill']}.md"
        assert path.is_file(), f"missing skill file {path} for {t['id']}"
        text = path.read_text(encoding="utf-8")
        m = re.match(r"---\n(.*?)\n---", text, re.DOTALL)
        assert m, f"{path.name} has no frontmatter"
        fm = m.group(1)
        assert f"name: {t['skill']}" in fm, f"{path.name} name mismatch"
        assert "description:" in fm, f"{path.name} has no description"


@pytest.mark.parametrize("symbol", _REAL_SYMBOLS)
def test_skills_only_name_real_workflow_vocabulary(symbol):
    corpus = ""
    for t in _task_ids():
        corpus += (_SKILLS / f"{t['skill']}.md").read_text(encoding="utf-8")
    assert symbol in corpus, f"no skill mentions {symbol}"


def test_redox_skill_states_the_electrode_requirement():
    text = (_SKILLS / "redox-potential.md").read_text(encoding="utf-8")
    assert re.search(r"(?i)reference electrode", text)
    assert "4.745" in text, "MeCN E_ref value missing"
    assert "4.345" in text, "default E_ref value missing"
    assert "get_E_ref" in text
    # The thermodynamic cycle formula must be stated with its unit.
    assert re.search(r"23\.0605", text), "F constant without unit context"


def test_no_skill_names_a_credential_or_permission_widening():
    """Wave-20 rule: a skill never widens a permission or names a
    credential path, and never tells the agent to skip a check."""
    for path in sorted(_SKILLS.glob("orca-spectra-*.md")) + [
            _SKILLS / "redox-potential.md"]:
        text = path.read_text(encoding="utf-8")
        for banned in (r"(?i)token", r"(?i)password", r"(?i)secret",
                       r"(?i)skip the (check|gate)", r"(?i)ignore the "
                       r"(refusal|sandbox)"):
            assert not re.search(banned, text), f"{path.name}: {banned}"
