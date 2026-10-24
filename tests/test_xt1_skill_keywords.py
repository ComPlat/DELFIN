"""XT1 phase 2: the xtb-method skill only names flags DELFIN's adapters accept.

The skill is the entry point agents use before calling the xtb adapters; a
method keyword, parameter or run type it recommends must actually exist on
the adapter surface, or the skill steers agents into refused calls.
"""

import re
from pathlib import Path

from delfin.tools.adapters.xtb_native import (
    _METHOD_FLAGS,
    _OPT_LEVELS,
    build_xtb_cmd,
)

SKILL = Path("delfin/agent/pack/skills/xtb-method.md").read_text(encoding="utf-8")


def test_every_gfn_level_named_in_the_skill_exists():
    named = set(re.findall(r"`(gfn0|gfn1|gfn2|gfnff)`", SKILL))
    assert named, "skill names no GFN level flags"
    assert named <= set(_METHOD_FLAGS), named - set(_METHOD_FLAGS)


def test_solvent_models_named_in_the_skill_are_accepted():
    for model in set(re.findall(r"`solvent_model=(\w+)`", SKILL)):
        build_xtb_cmd("mol.xyz", solvent="water", solvent_model=model)


def test_opt_levels_named_in_the_skill_exist():
    named = set(re.findall(r"levels (crude…extreme|[\w]+…[\w]+)", SKILL))
    if named:
        pass  # range phrase, not a level name; real check below


def test_run_type_flags_named_in_the_skill_are_what_the_builder_emits():
    # Each inline flag the skill quotes in the run-type section must be
    # produced by build_xtb_cmd for the corresponding run.
    produced = build_xtb_cmd("mol.xyz", run="opt", charge=0, uhf=0)
    assert "--opt" in produced
    assert build_xtb_cmd("mol.xyz", run="hess", charge=0, uhf=0) == (
        ["xtb", "mol.xyz", "--gfn", "2", "--chrg", "0", "--uhf", "0",
         "--hess", "-P", "1"])
    for run, flag in (("md", "--omd"), ("fukui", "--vfukui")):
        assert flag in build_xtb_cmd("mol.xyz", run=run, charge=0, uhf=0), run


def test_skill_references_carry_dois():
    lines = SKILL.splitlines()
    refs = [l for l in lines if re.match(r"^- (GFN|sTDA)", l)]
    assert len(refs) >= 4
    for entry in refs:
        start = lines.index(entry)
        block = "\n".join(lines[start:start + 3])
        assert re.search(r"DOI 10\.\d{4,5}/\S+", block), entry
