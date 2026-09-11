"""The Recalc tab submits a smart recalc, and says what the finished jobs were computed with.

It is there to edit CONTROL.txt and resubmit, and it submitted a classic
recalc, which keeps every finished job whatever the edit.  It now submits a
smart recalc, and records the CONTROL it replaces when no completed run has
recorded one (recalc_control), so the edit is known for a job from before the
record existed.  It also tells the validator whether the job's input is a
SMILES: going by the CONTROL file's own SMILES key, it refused 266 archived
jobs that name a SMILES there but run from an xyz input.
"""

from __future__ import annotations

from delfin import recalc_control


def test_the_recalc_tab_tells_the_validator_whether_the_input_is_a_smiles(tmp_path):
    from delfin.dashboard.tab_recalc import _builds_from_smiles

    (tmp_path / "input.txt").write_text("N 0.0 0.0 0.0\nH 0.0 0.94 0.0\n")
    assert _builds_from_smiles(tmp_path, "input_file=input.txt\nSMILES=CCO\n") is False
    (tmp_path / "input.txt").write_text("CCO\n")
    assert _builds_from_smiles(tmp_path, "input_file=input.txt\n") is True


def test_what_the_edit_replaces_is_recorded_once(tmp_path):
    assert recalc_control.remember_before_edit(tmp_path, "functional=PBE0\n")
    assert not recalc_control.remember_before_edit(tmp_path, "functional=B3LYP\n")
    assert recalc_control.previous_run(tmp_path)["control"] == "functional=PBE0\n"


def test_the_tab_submits_a_smart_recalc():
    import inspect
    from delfin.dashboard import tab_recalc

    source = inspect.getsource(tab_recalc.create_tab)
    assert "mode='delfin-recalc'," in source
    assert "delfin-recalc-classic" not in source
