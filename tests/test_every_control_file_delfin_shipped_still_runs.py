"""Every CONTROL file DELFIN ever shipped still validates and still reads.

A CONTROL file outlives the DELFIN that wrote it: people keep them, copy them
between jobs, and resubmit them years later.  Twice in one week a validator
change refused files DELFIN itself had handed out -- an IC the current ORCA
cannot compute (the template suggested it), a state list that named S0, an
OCCUPIER_method that still carried the choice it was offered as.  264 archived
files of one group were refused by a version that was otherwise fine.

So the templates are kept here, one per month of the project, and the run
checks that each still validates once the three things a user must choose are
chosen, and that reading it gives a run its settings.  The whole history is
swept as well wherever the repository's history is at hand.
"""

from __future__ import annotations

import ast
import subprocess
from pathlib import Path

import pytest

from delfin.config import read_control_file, set_control_value, validate_control_text

TEMPLATES = Path(__file__).parent / "fixtures" / "control_templates"

#: What a user chooses before a run; a template ships them as placeholders.
CHOICES = (("method", "classic"), ("charge", "0"), ("solvent", "water"))

#: One shipped template was broken on the day it shipped, with an empty
#: OCCUPIER_compare, and was corrected the same day (9dcce272).  It is kept
#: because a file of that day may still be out there, and it is named here so
#: the sweep below does not have to guess.
BROKEN_ON_THE_DAY = {"2026-08-14_0b1e3012": "OCCUPIER_compare is empty"}


def _filled(text: str) -> str:
    for key, value in CHOICES:
        text = set_control_value(text, key, value)
    return text


def _templates() -> list:
    return sorted(TEMPLATES.glob("*.txt"))


def test_the_templates_are_kept_here():
    names = [p.stem for p in _templates()]
    assert len(names) >= 12, "the templates DELFIN shipped are the fixture of this test"
    assert any(n.startswith("2025-09") for n in names), "the oldest one is missing"


@pytest.mark.parametrize("template", _templates(), ids=lambda p: p.stem)
def test_a_shipped_template_still_validates(template):
    errors = validate_control_text(_filled(template.read_text(encoding="utf-8")),
                                   converts_smiles=False)
    expected = BROKEN_ON_THE_DAY.get(template.stem)
    if expected:
        assert any(expected in e for e in errors), f"{template.stem} was broken for another reason"
        return
    assert errors == [], f"{template.stem} no longer validates"


@pytest.mark.parametrize("template", _templates(), ids=lambda p: p.stem)
def test_a_shipped_template_still_gives_a_run_its_settings(template, tmp_path):
    if template.stem in BROKEN_ON_THE_DAY:
        pytest.skip("shipped broken, see the test above")
    control = tmp_path / "CONTROL.txt"
    control.write_text(_filled(template.read_text(encoding="utf-8")), encoding="utf-8")

    config = read_control_file(str(control))

    assert config["method"] == "classic"
    assert str(config["charge"]) == "0"
    assert config["solvent"] == "water"
    for key in ("functional", "main_basisset", "PAL", "maxcore"):
        assert str(config.get(key) or "").strip(), f"{template.stem} lost {key}"


def _templates_from_history() -> list:
    """Every distinct template in the repository's history, newest last."""
    repo = Path(__file__).resolve().parents[1]
    log = subprocess.run(["git", "log", "--reverse", "--format=%H %ad", "--date=format:%Y-%m-%d",
                          "--", "delfin/define.py"],
                         cwd=repo, capture_output=True, text=True, timeout=120)
    if log.returncode != 0 or not log.stdout.strip():
        return []
    seen, found = set(), []
    for line in log.stdout.splitlines():
        sha, date = line.split()
        shown = subprocess.run(["git", "show", f"{sha}:delfin/define.py"],
                               cwd=repo, capture_output=True, text=True, timeout=120)
        if shown.returncode != 0 or not shown.stdout:
            continue
        try:
            tree = ast.parse(shown.stdout)
        except SyntaxError:
            continue
        template = None
        for node in tree.body:
            if isinstance(node, ast.Assign) and any(getattr(t, "id", "") == "TEMPLATE" for t in node.targets):
                if isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
                    template = node.value.value
        if template and template not in seen:
            seen.add(template)
            found.append((f"{date}_{sha[:8]}", template))
    return found


def test_every_template_in_the_history_still_validates():
    history = _templates_from_history()
    if not history:
        pytest.skip("no repository history here; the templates kept beside this test are checked above")

    refused = {name: validate_control_text(_filled(text), converts_smiles=False)
               for name, text in history}
    refused = {name: errors for name, errors in refused.items() if errors}

    assert len(history) >= len(_templates()), "the history holds fewer templates than the fixtures"
    assert sorted(refused) == sorted(BROKEN_ON_THE_DAY), (
        "a template DELFIN shipped is no longer accepted: "
        + "; ".join(f"{n}: {e[:1]}" for n, e in refused.items()))
