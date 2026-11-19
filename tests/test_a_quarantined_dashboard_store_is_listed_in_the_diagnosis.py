"""A store moved aside by a dashboard writer is named where a person looks.

The four writers fixed alongside this test keep a damaged store instead of
overwriting it. That is only half of "the loss is never silent": a kept file
nobody is told about is a file nobody acts on, and the first version of this
work reported it through a module-level note with no caller and through
``print``, which in a served dashboard reaches the server's stdout and no
person at all.

The route that does reach a person is the Settings tab's diagnosis, which
renders ``doctor.run_doctor`` and already scans the state root for
``*.unreadable-*``. These three dashboard stores live directly in that root
(``branch_list.json``, ``dashboard_header.json``,
``agent_open_sessions.json``), so the existing check covers them -- and this
test is what keeps that true: the scan is NOT recursive, so a later move of
any of the three into a subdirectory would silently drop it out of the
report.

Deliberately not covered, and said here rather than left to be discovered: a
job's ``delfin_import.json`` is quarantined in the job directory, which is
outside the state root and in no diagnosis. For that one the log line is the
only record.
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent import doctor


def _redirect(monkeypatch, tmp_path: Path) -> Path:
    home = tmp_path / "home"
    root = home / ".delfin"
    root.mkdir(parents=True)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    return root


def test_the_diagnosis_names_each_quarantined_dashboard_store(monkeypatch,
                                                              tmp_path):
    root = _redirect(monkeypatch, tmp_path)
    for name in ("branch_list.json", "dashboard_header.json",
                 "agent_open_sessions.json"):
        (root / f"{name}.unreadable-20261119-050000").write_text(
            "{not json", encoding="utf-8")

    rows = doctor._check_state_files({})

    kept = [r for r in rows if "damaged state file(s) kept aside" in r.get("detail", "")]
    assert kept, f"no row reports the files kept aside: {rows}"
    detail = kept[0]["detail"]
    assert "3 damaged state file(s)" in detail, detail
    for name in ("branch_list.json", "dashboard_header.json",
                 "agent_open_sessions.json"):
        assert name in detail, (
            f"{name} was moved aside and the diagnosis does not name it; "
            "the person is then told nothing. " + detail)


def test_a_store_quarantined_one_level_down_is_named_too(monkeypatch,
                                                         tmp_path):
    """The stores that are NOT in the doctor's named list live a level
    below the state root -- session records in agent_sessions/, the
    service's files in agent_service/. The scan missed them entirely,
    which made a writer there unable to point a person at this report.
    """
    root = _redirect(monkeypatch, tmp_path)
    (root / "agent_sessions").mkdir()
    (root / "agent_sessions" / "abc.json.unreadable-20261119-050000").write_text(
        "{torn", encoding="utf-8")

    rows = doctor._check_state_files({})

    kept = [r for r in rows if "damaged state file(s) kept aside" in r.get("detail", "")]
    assert kept, f"a quarantine one level down must be reported: {rows}"
    assert "abc.json" in kept[0]["detail"], kept[0]["detail"]


def test_the_scan_does_not_walk_the_whole_state_root(monkeypatch, tmp_path):
    """Two levels down is not searched, on purpose: the state root also
    holds workspaces and caches. This pins the boundary, so widening it
    later is a decision rather than an accident."""
    root = _redirect(monkeypatch, tmp_path)
    deep = root / "workspaces" / "a" / "b"
    deep.mkdir(parents=True)
    (deep / "x.json.unreadable-20261119-050000").write_text("{", encoding="utf-8")

    rows = doctor._check_state_files({})

    assert not [r for r in rows if "damaged state file(s) kept aside" in r.get("detail", "")], rows


def test_an_undamaged_state_root_reports_nothing_to_act_on(monkeypatch,
                                                           tmp_path):
    """The control: the check must not cry about a healthy machine."""
    _redirect(monkeypatch, tmp_path)

    rows = doctor._check_state_files({})

    assert all("damaged state file(s) kept aside" not in r.get("detail", "")
               or r["status"] == doctor.PASS for r in rows), rows
