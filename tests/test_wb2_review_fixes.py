"""WB2 phase-2 review fixes — red tests for the reviewer's F1-F3.

F1/F3: _NAME_RE.match lets a trailing newline through; the name derives
a file name, so "x\\n" must be rejected (fullmatch, the same rule
skill_share.py:56 guards skill names with).
F2: list_workflows must degrade over one corrupt record, not vanish.
"""

from __future__ import annotations

import json

import pytest

from delfin.tools import user_workflows as uw


GOOD_SPEC = {
    "name": "w",
    "steps": [{"step": "xtb_opt", "kwargs": {"charge": 0, "method": "XTB2"}}],
}


@pytest.mark.parametrize("bad", ["x\n", "x\r\n"])
def test_name_with_trailing_newline_is_rejected(bad):
    with pytest.raises(ValueError):
        uw.save(GOOD_SPEC, name=bad)


def test_extend_template_rejects_a_newline_in_new_name():
    with pytest.raises(ValueError):
        uw.extend_template("classic_opt_freq", new_name="x\n")


def test_list_workflows_skips_a_corrupt_record(tmp_path):
    uw.save(GOOD_SPEC, name="good_flow", root=tmp_path)
    broken = tmp_path / "broken_flow"
    broken.mkdir()
    (broken / "v001.json").write_text("{not json", encoding="utf-8")
    names = [w["name"] for w in uw.list_workflows(root=tmp_path)]
    assert names == ["good_flow"]
