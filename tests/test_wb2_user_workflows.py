"""WB2 phase 2 — user_workflows: save, load, list, export, extend.

The user's saved workflows live in the user's DELFIN data directory; the
path is resolved PER CALL (never an import-time constant) so tests can
redirect the home. Red control: every case here failed on the previous
commit with ImportError (module did not exist).
"""

from __future__ import annotations

import json

import pytest

from delfin.tools import user_workflows as uw
from delfin.tools import platform


# A spec that passes resolve/validate/lint: the classic template's shape
# (smiles -> xtb_opt -> orca_opt -> orca_freq), concrete kwargs.
GOOD_SPEC = {
    "name": "my_opt_freq",
    "description": "opt + freq on a small organic molecule",
    "steps": [
        {"step": "smiles_to_xyz", "kwargs": {"smiles": "CCO"},
         "label": "smiles"},
        {"step": "xtb_opt", "kwargs": {"charge": 0, "method": "XTB2"},
         "label": "xtb_preopt"},
        {"step": "orca_opt", "kwargs": {"charge": 0, "method": "B3LYP",
                                        "basis": "def2-SVP"},
         "label": "orca_opt"},
        {"step": "orca_freq", "kwargs": {"charge": 0, "method": "B3LYP",
                                         "basis": "def2-SVP"},
         "label": "orca_freq"},
    ],
}


def test_save_then_load_returns_the_identical_spec(tmp_path):
    path = uw.save(GOOD_SPEC, name="my_opt_freq", root=tmp_path)
    assert path.exists()
    loaded = uw.load("my_opt_freq", root=tmp_path)
    assert loaded == GOOD_SPEC


def test_save_stores_json_and_metadata(tmp_path):
    path = uw.save(GOOD_SPEC, name="my_opt_freq", root=tmp_path)
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["spec"] == GOOD_SPEC
    assert data["version"] == 1
    assert "saved_at" in data


def test_saving_same_name_versions_it(tmp_path):
    uw.save(GOOD_SPEC, name="my_opt_freq", root=tmp_path)
    changed = dict(GOOD_SPEC, description="v2")
    p2 = uw.save(changed, name="my_opt_freq", root=tmp_path)
    assert p2 != uw.load_path("my_opt_freq", version=1, root=tmp_path)
    assert uw.load("my_opt_freq", version=2, root=tmp_path) == changed
    # default load returns the newest version
    assert uw.load("my_opt_freq", root=tmp_path) == changed
    versions = uw.list_versions("my_opt_freq", root=tmp_path)
    assert [v["version"] for v in versions] == [1, 2]


def test_list_lists_saved_names(tmp_path):
    uw.save(GOOD_SPEC, name="my_opt_freq", root=tmp_path)
    uw.save(GOOD_SPEC, name="other_flow", root=tmp_path)
    names = [w["name"] for w in uw.list_workflows(root=tmp_path)]
    assert "my_opt_freq" in names
    assert "other_flow" in names


def test_load_missing_raises_keyerror(tmp_path):
    with pytest.raises(KeyError):
        uw.load("no_such_workflow", root=tmp_path)


def test_export_writes_a_contrib_style_template(tmp_path):
    uw.save(GOOD_SPEC, name="my_opt_freq", root=tmp_path)
    out = tmp_path / "export"
    path = uw.export_contrib(
        "my_opt_freq", out_dir=out, root=tmp_path,
        author="test",
    )
    assert path.exists()
    text = path.read_text(encoding="utf-8")
    # The exported file must load back through the pipeline round-trip:
    data = json.loads(text)
    rebuilt = uw.load_spec_from_export(data)
    assert rebuilt["name"] == "my_opt_freq"


def test_extended_template_leaves_the_original_untouched(tmp_path):
    """Extend classic_opt_freq: add an orca_sp at the end. The named
    template in delfin/tools/templates must be unchanged after build."""
    import delfin.tools.templates as templates

    before = templates.classic_opt_freq.to_dict()
    draft = uw.extend_template(
        "classic_opt_freq", add=[{"step": "orca_sp",
                                  "kwargs": {"charge": 0}}],
    )
    after = templates.classic_opt_freq.to_dict()
    assert before == after
    steps = [s["step"] for s in draft["steps"]]
    assert steps[-1] == "orca_sp"


def test_rerun_builds_the_saved_spec(tmp_path):
    uw.save(GOOD_SPEC, name="my_opt_freq", root=tmp_path)
    report = uw.validate_saved("my_opt_freq", root=tmp_path, geometry=True)
    assert report.ok, report.summary()
