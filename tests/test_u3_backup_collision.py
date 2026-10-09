"""Probe: two setting applies in one second must each keep their own dated
backup (the phase-2 invariant: a settings file is never overwritten or
deleted, and its prior content is always recoverable from the backup).

Repair's _dated_backup_path stamps at SECOND resolution
(datetime.now().strftime('%Y%m%d-%H%M%S')). Two applies inside the same
second collide on one backup name; os.replace then moves the second (freshly
rewritten) settings over the first backup, silently destroying the ORIGINAL
prior content. This probe pins that a unique backup survives per apply.
"""

from __future__ import annotations

import json

from delfin.agent import repair


def _row(setting) -> dict:
    return {"check": "c", "status": "PASS", "detail": "d", "fix": "f",
            "setting": setting}


def test_two_setting_applies_inside_one_second_each_keep_a_backup(tmp_path):
    settings = tmp_path / "settings.json"
    settings.write_text(json.dumps({"agent": {"other": 1}}))
    s1 = repair.plan([_row(("agent.a", "1"))])[0]
    s2 = repair.plan([_row(("agent.b", "2"))])[0]

    r1 = repair.apply(s1, approved=True, user_settings_path=settings)
    r2 = repair.apply(s2, approved=True, user_settings_path=settings)

    # both succeeded -- then the SECOND apply must not have clobbered the
    # first backup (which holds the true original) under a colliding name.
    backups = [p for p in tmp_path.iterdir() if p.name.endswith(".bak")]
    assert len(backups) == 2, [
        p.name for p in tmp_path.iterdir()]
    # one of the two backups must still hold the ORIGINAL prior content
    contents = [json.loads(p.read_text()) for p in backups]
    assert {"agent": {"other": 1}} in contents, contents
