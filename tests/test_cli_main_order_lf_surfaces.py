"""Integration controls over the REAL CLI entry point (order LF, phase 5).

Wave 4 taught the round this lesson: a module can be green while its
caller hands the public wrapper an argument it does not know. These
controls drive delfin.agent.cli.main(argv) -- parser, dispatch, exit
code, stdout/stderr -- exactly as the `delfin-agent` command runs it,
for the three new surfaces of order LF:

  * approvals answer <id> <choice>
  * approvals ls (preview line)
  * report --since (compactions line)

They were red before order LF landed (main rejected `answer`, ls had
no subject line, report had no compactions) -- the per-phase controls
showed that; these stay as the end-to-end net over the whole chain.
"""

from __future__ import annotations

import json
import os
import time

import pytest

from delfin.agent import file_confirm as _fc
from delfin.agent import turn_metrics
from delfin.agent import cli


NOW = 1_800_000_000.0


def _ask_record(room, rid):
    rec = {"id": rid, "session_id": "sess-x", "session_key": "pane-x",
           "kind": "ask", "tool": "ask_user_question", "command": "",
           "preview": "", "asked_at": time.time(), "pid": os.getpid(),
           "host": "here",
           "payload": {"question": "Restore how?",
                       "options": [{"label": "Allow restore now"},
                                   {"label": "I'll restore manually"}],
                       "multiSelect": False}}
    p = room / f"{rid}.request.json"
    p.write_text(json.dumps(rec), encoding="utf-8")
    os.chmod(p, 0o600)


@pytest.fixture
def headless_room(tmp_path, monkeypatch):
    room = tmp_path / "approvals"
    room.mkdir()
    monkeypatch.setattr(_fc, "requests_dir", lambda: room)
    return room


def _main(argv, capsys):
    rc = cli.main(argv)
    out = capsys.readouterr()
    return rc, out.out, out.err


class TestMainApprovalsAnswer:
    def test_full_path_writes_the_pick(self, headless_room, capsys):
        _ask_record(headless_room, "4001-nn")
        rc, out, err = _main(
            ["approvals", "answer", "4001-nn", "I'll restore manually"],
            capsys)
        assert rc == 0, err
        assert "4001-nn" in out
        raw = json.loads(
            (headless_room / "4001-nn.answer.json").read_text(encoding="utf-8"))
        assert raw == {"id": "4001-nn", "decision": "choose",
                       "answers": ["I'll restore manually"],
                       "by": os.environ.get("USER", "")}

    def test_full_path_bad_pick_exits_two(self, headless_room, capsys):
        _ask_record(headless_room, "4002-oo")
        rc, out, err = _main(["approvals", "answer", "4002-oo", "7"],
                             capsys)
        assert rc == 2
        assert "not one of the 2 options" in err

    def test_full_path_missing_choice_arg_is_a_usage_error(self,
                                                           headless_room,
                                                           capsys):
        _ask_record(headless_room, "4003-pp")
        # argparse exits 2 on a missing positional, via SystemExit.
        with pytest.raises(SystemExit) as exc:
            cli.main(["approvals", "answer", "4003-pp"])
        assert exc.value.code == 2


class TestMainApprovalsLs:
    def test_full_path_shows_the_subject(self, headless_room, capsys,
                                         monkeypatch):
        monkeypatch.setenv("COLUMNS", "200")
        rec = {"id": "4004-qq", "session_id": "s", "kind": "confirm",
               "tool": "bash", "command": "git push origin main",
               "preview": "", "asked_at": time.time()}
        (headless_room / "4004-qq.request.json").write_text(
            json.dumps(rec), encoding="utf-8")
        rc, out, err = _main(["approvals", "ls"], capsys)
        assert rc == 0
        assert "git push origin main" in out

    def test_full_path_empty_says_nothing_waiting(self, headless_room,
                                                  capsys):
        rc, out, err = _main(["approvals"], capsys)
        assert rc == 0
        assert "(nothing waiting)" in out


class TestMainReportSince:
    def test_full_path_shows_compactions(self, tmp_path, monkeypatch,
                                         capsys):
        traces = tmp_path / "tool_traces"
        traces.mkdir()
        archive = tmp_path / "transcript_archive"
        archive.mkdir()
        metrics = tmp_path / "turn_metrics"
        monkeypatch.setattr(turn_metrics, "_DIR", metrics)
        (traces / "nacht-z.jsonl").write_text(json.dumps({
            "ts": NOW, "tool": "read_file", "input": "", "output": "",
            "duration_ms": 10, "ok": True, "error": ""}) + "\n",
            encoding="utf-8")
        (archive / "nacht-z.jsonl").write_text(
            json.dumps({"compacted_at": NOW, "n_messages": 1,
                        "info": {}, "messages": []}) + "\n"
            + json.dumps({"compacted_at": NOW - 60, "n_messages": 1,
                          "info": {}, "messages": []}) + "\n",
            encoding="utf-8")
        # collect() must be pointed at the fabricated stores; main()
        # resolves trace_root itself, so patch the resolution boundary
        # it uses -- the module attribute, as the other CLI tests do.
        monkeypatch.setattr("delfin.agent.tool_trace._DIR", traces)
        monkeypatch.setattr(
            "delfin.agent.round_report._transcript_archive_path",
            lambda: archive)
        rc, out, err = _main(["report", "--since", "1h"], capsys)
        assert rc == 0, err
        assert "compactions: 2" in out
        assert "compactions" in out.split("TOTAL", 1)[1]
