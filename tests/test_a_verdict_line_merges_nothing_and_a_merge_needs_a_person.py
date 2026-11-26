"""Package g3: a verdict is recorded, a merge needs a person.

The operator procedure this replaces watched a coordination file for a
line ``REAL-RUN OK #<n>`` and then started a merge waiter that merged PR
<n> with administrator rights once checks were green. Anyone who could
write that file could cause an administrator merge. The vocabulary is
kept; the trigger is not.

Red on the base tree for the reason each test names: on the base there is
no ``swarm verdict`` verb and no ``swarm merge-run`` verb at all, so the
real parser rejects them ("unrecognized arguments"), and
``delfin.agent.swarm.merge`` has neither ``parse_verdict_lines`` nor
``decide_merge`` (AttributeError). The control case
(``test_the_existing_merge_plan_still_prints_its_commands``) passes on
both trees: it measures that this change did not disturb what was there.

Nothing here talks to GitHub: ``merge.current_head_of`` is replaced in
every case that needs a head, and one case drives the fail-closed path
where ``gh`` cannot run at all. No skip depends on what the host has.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import cli_swarm
from delfin.agent.swarm import lead as swarm_lead
from delfin.agent.swarm import merge as swarm_merge
from delfin.agent.swarm import wave


def _run(argv: list[str]) -> int:
    """Dispatch *argv* through the REAL parser cli.py builds, as the
    existing swarm tests do: `build_parser` is cli.py's own, and
    `args.func` is the production handler (cli_swarm.run_args)."""
    from delfin.agent.cli import build_parser
    args = build_parser().parse_args(argv)
    assert args.func is cli_swarm.run_args
    return args.func(args)


# -- the verdict: recorded, with its author, merging nothing ---------------------

def test_a_verdict_line_written_into_a_file_merges_nothing(tmp_path, capsys):
    """A ``REAL-RUN OK #7`` line in any file is data, and the only merge
    path refuses without a person's yes.

    The input that would fail this: a parse that returns ``merges: True``,
    or a ``swarm merge-run`` that prints a ``gh pr merge`` command for a
    pull request nobody confirmed.
    """
    coordination = tmp_path / "coordination.md"
    coordination.write_text(
        "notes\nREAL-RUN OK #7\nREAL-RUN FAIL #8\n", encoding="utf-8")

    found = swarm_merge.parse_verdict_lines(
        coordination.read_text(encoding="utf-8"))
    assert [f["pr"] for f in found] == [7, 8]
    assert all(f["merges"] is False and f["actionable"] is False
               for f in found), (
        "a line found in a file must never carry a decision to merge")

    # The unconfirmed merge of the very PR that line names: refused, and
    # no merge command is printed for it.
    rc = _run(["swarm", "merge-run", "7", "--head", "a" * 40])
    out = capsys.readouterr()
    assert rc == 1, "an unconfirmed merge must not report success"
    assert "gh pr merge" not in out.out, (
        "no merge command may be printed without the person's yes")
    assert "ask" in out.out.lower() or "your decision" in out.out.lower()


def test_a_recorded_verdict_carries_its_author_and_merges_nothing(tmp_path,
                                                                 monkeypatch,
                                                                 capsys):
    """``swarm verdict`` records who said it, and says it merges nothing.

    Fails on the base: the verb does not exist, so the parser rejects it.
    """
    _seed_run(tmp_path, monkeypatch, {"a": _member(tmp_path, "a")})
    rc = _run(["swarm", "verdict", "r1", "a", "ok", "--by", "the maintainer"])
    out = capsys.readouterr()
    assert rc == 0, out.err
    assert "by the maintainer" in out.out
    assert "merges nothing" in out.out
    events = _events(tmp_path, "r1")
    assert [e["kind"] for e in events] == ["realrun_verdict"]
    assert events[0]["by"] == "the maintainer"
    assert events[0]["verdict"] == "ok"


def test_a_verdict_without_a_person_is_refused(tmp_path, monkeypatch, capsys):
    """A verdict nobody signed is not a verdict (the operator file had no
    author at all). ``--by`` is required, and a blank one is refused too."""
    _seed_run(tmp_path, monkeypatch, {"a": _member(tmp_path, "a")})
    with pytest.raises(SystemExit):
        _run(["swarm", "verdict", "r1", "a", "ok"])
    rc = _run(["swarm", "verdict", "r1", "a", "ok", "--by", "   "])
    assert rc == 1


# -- the merge: confirmed, and pinned to the head that was checked ---------------

def test_a_confirmed_merge_with_a_moved_head_is_refused(tmp_path, monkeypatch,
                                                       capsys):
    """The checks were read for head A; the PR head is now B. Refused.

    The pin is the point: a check that answered for a different commit
    does not count. The input that would fail this: a decision that merges
    on ``--confirmed`` alone.
    """
    head_a, head_b = "a" * 40, "b" * 40
    monkeypatch.setattr(swarm_merge, "current_head_of",
                        lambda pr, **kw: head_b)
    rc = _run(["swarm", "merge-run", "7", "--head", head_a, "--confirmed",
               "--checks-ok"])
    out = capsys.readouterr()
    assert rc == 1, "a moved head must not merge"
    assert "stale" in out.out
    assert "gh pr merge" not in out.out


def test_a_moved_head_cannot_be_revived_by_the_admin_flag(monkeypatch):
    """Administrator rights do not make a stale check valid."""
    decided = swarm_lead.swarm_merge(
        "7", "a" * 40, confirmed=True, checks_ok=True,
        current_head="b" * 40, admin_override=True)
    assert decided["state"] == "stale"
    assert decided["allow"] is False
    assert decided["merge_command"] == ""


def test_a_confirmed_merge_on_the_checked_head_prints_the_command(monkeypatch):
    """The happy path: confirmed, pin matches, checks green -> the command
    a person runs, pinned to that head. Still ``executed: False``."""
    head = "c" * 40
    monkeypatch.setattr(swarm_merge, "current_head_of",
                        lambda pr, **kw: head)
    rc = _run(["swarm", "merge-run", "7", "--head", head, "--confirmed",
               "--checks-ok"])
    assert rc == 0
    decided = swarm_lead.swarm_merge(
        "7", head, confirmed=True, checks_ok=True, current_head=head)
    assert decided["allow"] is True
    assert decided["executed"] is False
    assert decided["merge_command"] == (
        f"gh pr merge 7 --match-head-commit {head}")


def test_admin_rights_need_the_flag_in_the_same_command():
    """Red checks are refused; the override must be asked for explicitly."""
    head = "d" * 40
    plain = swarm_merge.decide_merge("7", head, confirmed=True,
                                     checks_ok=False, current_head=head)
    assert plain["state"] == "checks_red" and plain["allow"] is False
    asked = swarm_merge.decide_merge("7", head, confirmed=True,
                                     checks_ok=False, current_head=head,
                                     admin_override=True)
    assert asked["allow"] is True
    assert "--admin" in asked["merge_command"] and asked["executed"] is False


def test_an_unreadable_head_fails_closed(monkeypatch):
    """When ``gh`` cannot run, the head is unknown and the merge refuses.

    Failing open here would merge on a check nobody could pin.
    """
    def no_gh(*a, **kw):
        raise FileNotFoundError("gh")
    monkeypatch.setattr("subprocess.run", no_gh)
    assert swarm_merge.current_head_of("7") == ""
    decided = swarm_lead.swarm_merge("7", "e" * 40, confirmed=True,
                                     checks_ok=True)
    assert decided["allow"] is False and decided["state"] == "stale"


def test_only_the_boolean_is_a_yes():
    """A truthy string is an invalid attempt, not a confirmation
    (the ``lead.request_swarm_start`` precedent)."""
    head = "f" * 40
    decided = swarm_merge.decide_merge("7", head, confirmed="yes",
                                       checks_ok=True, current_head=head)
    assert decided["state"] == "ask_user" and decided["allow"] is False


# -- the fail verdict: the member is resumed with the failure text ---------------

def test_a_fail_verdict_resumes_the_member_that_failed(tmp_path, monkeypatch,
                                                       capsys):
    """FAIL resumes *that* member with the failure text, as the script did.

    Two members are seeded; only the one named is resumed, and its next
    brief carries the failure. Fails on the base: no ``swarm verdict``.
    """
    launched: list[tuple[str, dict]] = []

    def fake_launch(run, task, repo, session=""):
        launched.append((run, task))
        return {"pid": 999999, "log": str(tmp_path / "member.log"),
                "worktree": str(tmp_path), "task": task}

    monkeypatch.setattr(wave, "_launch", fake_launch)
    monkeypatch.setattr(wave, "alive", lambda pid: False)
    _seed_run(tmp_path, monkeypatch,
              {"a": _member(tmp_path, "a"), "b": _member(tmp_path, "b")})

    rc = _run(["swarm", "verdict", "r1", "b", "fail", "--by", "the maintainer",
               "--detail", "the ORCA output had no SCF convergence"])
    out = capsys.readouterr()
    assert rc == 0, out.err
    assert "resumed" in out.out
    assert len(launched) == 1, "exactly one member is resumed"
    brief = launched[0][1]["brief"]
    assert "no SCF convergence" in brief, (
        "the failure text must reach the member, as the operator script did")
    state = wave.load_state("r1")
    assert int(state["agents"]["a"].get("resumes", 0)) == 0, (
        "the other member must not be touched")


def test_a_fail_verdict_for_a_running_member_does_not_double_start(
        tmp_path, monkeypatch, capsys):
    """A member still running is not resumed; the verdict is still recorded."""
    monkeypatch.setattr(wave, "alive", lambda pid: True)
    _seed_run(tmp_path, monkeypatch, {"a": _member(tmp_path, "a")})
    rc = _run(["swarm", "verdict", "r1", "a", "fail", "--by", "the maintainer",
               "--detail", "boom"])
    out = capsys.readouterr()
    assert rc == 0
    assert "not resumed" in out.out


# -- control case: passes on the base tree and on this branch -------------------

def test_the_existing_merge_plan_still_prints_its_commands(capsys, tmp_path,
                                                           monkeypatch):
    """The verb that existed before this change still behaves the same.

    Green on both trees: it measures that the new verbs did not disturb
    the old one.
    """
    monkeypatch.chdir(tmp_path)
    rc = _run(["swarm", "merge", "core"])
    out = capsys.readouterr()
    assert rc == 0
    assert "git merge swarm/core-builder" in out.out


# -- helpers --------------------------------------------------------------------

def _member(tmp_path: Path, tid: str) -> dict:
    return {"pid": 2 ** 22, "log": str(tmp_path / f"{tid}.log"),
            "worktree": str(tmp_path), "task": {"brief": f"do {tid}",
                                                "seconds": 60}}


def _seed_run(tmp_path: Path, monkeypatch, agents: dict[str, dict]) -> None:
    """A wave run whose state lives under *tmp_path* (never the user's)."""
    monkeypatch.setattr(wave, "root", lambda: tmp_path)
    run_dir = tmp_path / "r1"
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "state.json").write_text(
        json.dumps({"repo": str(tmp_path), "agents": agents}), encoding="utf-8")


def _events(tmp_path: Path, run: str) -> list[dict]:
    path = tmp_path / run / "events.jsonl"
    if not path.exists():
        return []
    return [json.loads(line) for line in
            path.read_text(encoding="utf-8").splitlines() if line.strip()]
