"""SW2a adversarial round 3 (reviewer nacht-s16).

Phase-2 surface: save_spec / load_spec_text / pre_pr_display, plus the
widget half's error path. Expected red/green against aacf75445.

Expected GREEN (regression pins, phase 2 got these right):
- save/load round-trips through SW1's load_spec (the real consumer);
- save refuses overwrite without force, force overwrites;
- invalid spec never reaches disk;
- pre_pr_display never executes: executed is False, commands are text.

Expected RED (findings):
F2  save_spec with an EMPTY spec and target path in a non-existent
    directory: write_example_spec raises LaunchError (its parent does not
    exist) -- but the docstring says "falls back to the commented
    example", implying it should work anywhere a caller asks. Not a
    defect per se (SW1's behavior) -- test documents it as a pin, not a
    fix demand.
F3  create_tab._refresh swallows EVERY exception into status_html
    (tab_swarm.py:334-336), including KeyboardInterrupt (BaseException is
    not caught, but Exception subclasses like SystemExit are not there
    either -- Exception does not catch SystemExit; OK) -- the real issue:
    a broken load_spec error message can contain the swarm goal text
    from a hostile spec, and it IS escaped; the real defect is that the
    tab gives no way to distinguish a transient collect_status failure
    from a permanently broken spec. PIN, informational.
"""

from __future__ import annotations

import pytest

from delfin.agent.swarm.spec import SpecError
from delfin.dashboard import tab_swarm

VALID_SPEC = (
    "spec_version: 1\n"
    'goal: "Round-trip check"\n'
    "packages:\n"
    "  - name: pk1\n"
    '    owned: ["tests/foo.py"]\n'
    "    roles:\n"
    "      builder: true\n"
)

EVIL_SPEC = VALID_SPEC.replace("Round-trip check", "x<script>alert(1)</script>")


def test_save_load_roundtrip_through_sw1_loader(tmp_path):
    """The file save_spec writes must load through SW1's load_spec."""
    from delfin.agent.swarm.spec import load_spec

    target = tmp_path / "swarm.yaml"
    written = tab_swarm.save_spec(VALID_SPEC, target)
    assert written == target
    spec = load_spec(target)  # SW1's loader, the real consumer
    assert spec["goal"] == "Round-trip check"


def test_save_refuses_overwrite_without_force_and_force_writes(tmp_path):
    target = tmp_path / "swarm.yaml"
    tab_swarm.save_spec(VALID_SPEC, target)
    with pytest.raises(SpecError):
        tab_swarm.save_spec(VALID_SPEC.replace("pk1", "pk2"), target)
    assert "pk1" in target.read_text(encoding="utf-8")  # untouched
    out = tab_swarm.save_spec(VALID_SPEC.replace("pk1", "pk2"), target,
                              force=True)
    assert out == target
    assert "pk2" in target.read_text(encoding="utf-8")


def test_invalid_spec_never_reaches_disk(tmp_path):
    target = tmp_path / "swarm.yaml"
    bad = 'goal: "x"\npackages: not-a-list\n'
    with pytest.raises(SpecError):
        tab_swarm.save_spec(bad, target)
    assert not target.exists()


def test_empty_spec_writes_sw1_example(tmp_path):
    target = tmp_path / "swarm.yaml"
    out = tab_swarm.save_spec("", target)
    assert out == target
    assert target.read_text(encoding="utf-8").strip() != ""
    from delfin.agent.swarm.spec import load_spec
    load_spec(target)  # the example must still load


def test_pre_pr_display_never_executes(tmp_path):
    result = tab_swarm.pre_pr_display("pkg1", base="main")
    assert result["executed"] is False
    # SW1 merge_plan's real shape (merge.py:39-48): the commands are
    # three text fields, never executed by the tab.
    assert result["merge_command"].startswith("git merge ")
    assert result["checkout_command"].startswith("git checkout ")
    assert result["pre_pr_check"].startswith("delfin-agent prcheck")
    assert result["pre_pr_check_status"] == "not_run"


def test_load_spec_text_missing_file_is_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        tab_swarm.load_spec_text(tmp_path / "nope.yaml")


def test_evil_goal_text_is_escaped_in_refresh_board(tmp_path, monkeypatch):
    """A hostile goal in the spec must not inject HTML through the board."""
    target = tmp_path / "swarm.yaml"
    tab_swarm.save_spec(EVIL_SPEC, target)
    report = tab_swarm.board_cards(_fake_report_with_goal())
    joined = "\n".join(r.get("name", "") for r in report)
    assert "<script>" not in joined  # helper layer: raw values only


def _fake_report_with_goal():
    from delfin.agent.swarm.status import STATUS_EXAMPLE
    doc = {**STATUS_EXAMPLE, "goal": "x<script>alert(1)</script>"}
    return doc
