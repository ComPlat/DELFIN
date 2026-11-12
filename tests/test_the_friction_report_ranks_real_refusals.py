"""The friction report ranks refusals, slow turns and silent exits from the
audit log and the turn metrics. It keeps only command heads, and it
proposes one task per bucket that is large enough.

The reason strings below have the shapes the permission gate writes on a
real host (taken from a Pi's audit log on 2026-11-12, with names and paths
replaced).
"""

from __future__ import annotations

import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

from delfin.agent import friction_report as fr

WATCHER = ("'bash' ({cmd}) was not run: the action watcher (the model's "
           "second opinion behind the permission rules) could not judge it "
           "-- the confirming model could not answer, so the action needs "
           "two families to allow and asks.")
OUTSIDE = ("blocked: this command would print '/tmp/x.txt'. read denied: "
           "'/tmp/x.txt' is outside the allowed workspace roots and no "
           "confirm callback is configured.")
REFUSED = ("'bash' ({cmd}) was allowed by the permission rules, but the "
           "action classifier refused it: goes beyond ordinary edits")
PAUSE = ("'bash' ({cmd}) was not run: three actions were refused in this "
         "turn. Unattended mode: the session is PAUSED")
STASH = ("command rejected by deny-pattern '(?:sudo\\s+)?git\\s+stash\\b': "
         "refusing to run.")

SECRET = "the person's private message 4711"


def _iso(ts: float) -> str:
    return datetime.fromtimestamp(ts, timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ")


def _rec(ts, decision, reason="", command="", sid="s1", tool="bash"):
    return {"ts": _iso(ts), "tool": tool, "decision": decision,
            "reason": reason, "command": command, "session_id": sid}


def _write(path: Path, recs) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(json.dumps(r) + "\n" for r in recs))


def _home(tmp_path: Path, now: float) -> Path:
    d = tmp_path / ".delfin"
    pytest_cmd = ("set -o pipefail; ~/v/bin/python -m pytest -q "
                  f"tests/test_a.py -k '{SECRET}' 2>&1 | tail -3")
    commit = f'git commit -am "{SECRET}"'
    _write(d / "audit.log", [
        *[_rec(now - 60 * i, "denied", WATCHER.format(cmd=pytest_cmd),
               pytest_cmd, sid=f"s{i % 2}") for i in range(4)],
        _rec(now - 50, "denied", WATCHER.format(cmd=commit), commit),
        _rec(now - 40, "denied", OUTSIDE, "cat /tmp/x.txt", sid="s2"),
        _rec(now - 30, "denied", REFUSED.format(cmd="git grep -l x"),
             "git grep -l x"),
        _rec(now - 20, "denied", PAUSE.format(cmd="true"), "true"),
        _rec(now - 10, "ok", command="ls"),
        _rec(now - 9, "error", "exited with code 1", "false"),
        # Outside the window: must not count.
        _rec(now - 3 * 86400, "denied", OUTSIDE, "cat /tmp/old"),
    ])
    # A weekly rotation written inside the window counts too.
    rot = d / "audit-2026-W45.log"
    _write(rot, [_rec(now - 3000, "denied", STASH, "git stash", sid="s3")])
    os.utime(rot, (now - 3000, now - 3000))
    _write(d / "agent_metrics.jsonl", [
        {"ts": now - 100, "model": "m-slow", "duration_s": 900.0,
         "session_id": "s1", "silent_exit": False},
        {"ts": now - 90, "model": "m-slow", "duration_s": 20.0,
         "session_id": "s1", "silent_exit": True},
        {"ts": now - 80, "model": "m-fast", "duration_s": 5.0,
         "session_id": "s2"},
        {"ts": now - 5 * 86400, "model": "m-old", "duration_s": 999.0},
    ])
    return tmp_path


def test_command_heads_keep_no_argument_or_quoted_text():
    heads = {
        "set -o pipefail; ~/v/bin/python -m pytest -q tests/a.py | tail": "pytest",
        f'git commit -am "{SECRET}"': "git commit",
        "git -C sub --no-pager log -1": "git log",
        f"python3 -c 'print(\"{SECRET}\")'": "python -c",
        "python - <<EOF": "python stdin",
        "for f in a b; do ~/v/bin/python -m pytest $f; done": "pytest",
        "timeout 60 grep -rn x .": "grep",
        "FOO=1 env -u BAR rg x": "rg",
        "": "",
    }
    for command, head in heads.items():
        assert fr.command_head(command) == head, command


def test_each_real_reason_shape_lands_in_its_bucket():
    assert fr.bucket_of(WATCHER.format(cmd="pytest")) == "watcher_unavailable"
    assert fr.bucket_of(OUTSIDE) == "read_outside_workspace"
    assert fr.bucket_of(REFUSED.format(cmd="grep x")) == "watcher_refused"
    assert fr.bucket_of(PAUSE.format(cmd="true")) == "refusal_pause"
    assert fr.bucket_of(STASH) == "deny_pattern"
    assert fr.bucket_of("something nobody wrote a bucket for") == "other"


def test_buckets_are_ranked_counted_and_windowed(tmp_path):
    now = time.time()
    data = fr.collect(since_s=now - 86400, home=_home(tmp_path, now))
    by_key = {b["key"]: b for b in data["buckets"]}
    assert data["buckets"][0]["key"] == "watcher_unavailable"
    w = by_key["watcher_unavailable"]
    assert (w["count"], w["sessions"]) == (5, 2)
    assert dict(w["heads"]) == {"pytest": 4, "git commit": 1}
    # The out-of-window refusal is not counted; the rotation is.
    assert by_key["read_outside_workspace"]["count"] == 1
    assert by_key["deny_pattern"]["heads"] == [("git stash", 1)]
    assert by_key["slow_turns"]["count"] == 1
    assert by_key["slow_turns"]["models"] == [("m-slow", 1)]
    assert by_key["silent_exits"]["count"] == 1
    assert "m-old" not in data["median_turn_s"]
    # 11 audit actions in the window: 10 in the active log, 1 in the rotation.
    assert (data["actions"], data["turns"]) == (11, 3)
    # Only a bucket of MIN_FOR_TASK or more gets a task.
    assert w["task"] and not by_key["read_outside_workspace"]["task"]


def test_the_rendering_carries_no_reason_path_or_message_text(tmp_path):
    now = time.time()
    data = fr.collect(since_s=now - 86400, home=_home(tmp_path, now))
    text = fr.render_text(data) + json.dumps(data)
    for leak in (SECRET, "/tmp/x.txt", "tests/test_a.py", "goes beyond",
                 "second opinion"):
        assert leak not in text
    assert "Proposed tasks" in text and "T1 [watcher_unavailable, 5]" in text


def test_a_saved_day_is_private_and_the_next_day_shows_the_change(tmp_path):
    now = time.time()
    data = fr.collect(since_s=now - 86400, home=_home(tmp_path, now))
    save_dir = tmp_path / "friction"
    earlier = dict(data, generated_s=now - 86400,
                   buckets=[dict(data["buckets"][0], count=2)])
    path = fr.save(earlier, save_dir)
    assert oct(path.stat().st_mode & 0o777) == "0o600"
    assert oct(save_dir.stat().st_mode & 0o777) == "0o700"
    today = time.strftime("%Y-%m-%d", time.gmtime(now))
    prev = fr.previous_report(save_dir, today)
    assert prev is not None
    assert "5 (+3 vs last report)" in fr.render_text(data, prev)


def test_missing_stores_give_an_empty_report(tmp_path):
    data = fr.collect(since_s=0.0, home=tmp_path)
    assert data["buckets"] == [] and data["actions"] == 0
    assert "No refusals" in fr.render_text(data)


def test_the_cli_command_reads_the_users_stores(tmp_path, monkeypatch,
                                                 capsys):
    from delfin.agent import cli
    now = time.time()
    _home(tmp_path, now)
    monkeypatch.setenv("HOME", str(tmp_path))
    assert cli.main(["friction", "--since", "1d", "--save"]) == 0
    out = capsys.readouterr().out
    assert "1. The action watcher could not judge" in out
    assert "saved:" in out and SECRET not in out
    assert list((tmp_path / ".delfin" / "friction").glob("*.json"))
    assert cli.main(["friction", "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["buckets"]
