"""A sheet read that is cut for the context says what was cut and what next.

Measured case: a 737 x 21 sheet was to be split into one sheet per value of
one column. The model read it with read_document nine times, five of the
calls byte-identical (max_cols=7, max_rows=200, start_row=120). Each result
was tens of thousands of characters; the agent loop cut it to the model's
tool-result cap (5120 chars) by keeping head and tail, so the model saw the
first three rows, the last nine, a "chars from the middle omitted" marker and
the reader's own "showing rows 120-319" note -- nothing said which rows were
missing or how to fetch them. It re-read the same window until the repeat
guard stopped the turn, and the turn ended with no answer text.

These tests drive the real headless engine through the scripted fake_model
endpoint (the same route as test_g1_degenerate_loop_fires_in_a_run) and check
what the MODEL receives and what the PERSON is left with.
"""

from __future__ import annotations

import io
import json
from pathlib import Path

import openpyxl

from tests.agent_terminal_probes import fake_model


_READ = {"path": "ChemOffice.xlsx", "sheet": "ChemOffice1",
         "max_cols": 7, "max_rows": 200, "start_row": 120}


def _book(path: Path) -> None:
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "ChemOffice1"
    ws.append([f"Header {c}" for c in range(1, 22)])
    for r in range(736):
        ws.append([f"Compound-{r:04d} sample text in cell {c}"
                   for c in range(1, 22)])
    wb.save(path)


def _small_book(path: Path) -> None:
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "ChemOffice1"
    ws.append(["Name", "Group"])
    for r in range(5):
        ws.append([f"C{r}", "A" if r % 2 else "B"])
    wb.save(path)


class _Launch:
    def __init__(self, plan_fn):
        self._plan_fn = plan_fn
        self.tool_rows: list[str] = []
        self.seen_rows: list[str] = []

    def plan(self, messages):
        for m in messages or []:
            text = fake_model._text_of(m)
            if m.get("role") == "tool" and text not in self.tool_rows:
                self.tool_rows.append(text)
            if text and text not in self.seen_rows:
                self.seen_rows.append(text)
        return self._plan_fn(messages)


def _read_once(messages):
    if messages and messages[-1].get("role") == "tool":
        return {"text": "Read the sheet."}
    return {"tool": "read_document", "args": dict(_READ)}


def _read_small(messages):
    if messages and messages[-1].get("role") == "tool":
        return {"text": "Read the sheet."}
    return {"tool": "read_document", "args": {"path": "ChemOffice.xlsx"}}


def _read_forever(_messages):
    return {"tool": "read_document", "args": dict(_READ)}


def _run(monkeypatch, tmp_path, capsys, plan_fn, book=_book):
    from delfin.agent import api_client as _A
    from delfin.agent import cli as agent_cli
    from delfin.agent import mcp_client as M
    from delfin.agent import model_capabilities as mc
    from delfin.agent import provider_profile as PP

    home = tmp_path / "home"
    ws = home / "ws"
    ws.mkdir(parents=True)
    book(ws / "ChemOffice.xlsx")
    monkeypatch.chdir(ws)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))

    class _EmptyReg:
        def discover_all(self):
            return []

        def discover_resources(self):
            return []

        def discover_prompts(self):
            return []

    cap = mc.ModelCapabilities(model=fake_model.MODEL, provider="ollama",
                               context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: cap)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    baseline = tmp_path / "learned_profiles.json"
    baseline.write_bytes(PP._DEFAULT_PATH.read_bytes())
    monkeypatch.setattr(PP, "_DEFAULT_PATH", baseline)
    monkeypatch.setattr(agent_cli, "_save_session",
                        lambda engine, root, **kw: "sid-1")
    monkeypatch.setattr(agent_cli, "_resume_or_create",
                        lambda engine, args: "")
    stdin = io.StringIO("")
    stdin.isatty = lambda: False              # type: ignore[method-assign]
    monkeypatch.setattr("sys.stdin", stdin)

    launch = _Launch(plan_fn)
    monkeypatch.setattr(fake_model, "plan", launch.plan)
    srv, port = fake_model.serve(log_path=str(tmp_path / "model.log"))
    monkeypatch.setenv("OLLAMA_HOST", f"http://127.0.0.1:{port}")
    monkeypatch.setenv("DELFIN_CLASSIFIER", "off")
    monkeypatch.setenv("DELFIN_ANNOUNCE_CONTINUATIONS", "0")
    # The cap the product resolves for the measured model.
    monkeypatch.setattr(_A, "_resolve_tool_result_cap", lambda *a, **k: 5120)
    monkeypatch.setattr(_A, "_resolve_max_tool_rounds", lambda *a, **k: 12)
    try:
        agent_cli.main(["-p", "split the sheet by column E", "--cwd", str(ws),
                        "--provider", "ollama", "--model", fake_model.MODEL,
                        "--permission-mode", "default",
                        "--output-format", "json"])
    finally:
        srv.shutdown()
    result = None
    for line in capsys.readouterr().out.splitlines():
        line = line.strip()
        if line.startswith("{"):
            try:
                result = json.loads(line)
            except ValueError:
                continue
    return result, launch


def _sheet_rows(launch) -> list[str]:
    return [t for t in launch.tool_rows if "rows x 21 columns" in t
            or "rows x 2 columns" in t]


def test_a_cut_sheet_read_names_the_rows_shown_hidden_and_the_next_call(
        monkeypatch, tmp_path, capsys):
    _result, launch = _run(monkeypatch, tmp_path, capsys, _read_once)
    rows = _sheet_rows(launch)
    assert rows, "the model never received the sheet read"
    seen = rows[0]
    assert len(seen) <= 5120 + 2000          # cap plus the data fence
    assert "chars from the middle omitted" not in seen
    assert "rows 120-" in seen and "are NOT shown" in seen
    assert "start_row=" in seen and "max_rows=" in seen
    assert 'sheet="ChemOffice1"' in seen and "max_cols=7" in seen
    assert "openpyxl" in seen and "bash" in seen
    # The stale note claiming all 200 rows were shown is gone.
    assert "showing rows 120-319" not in seen


def test_a_repeated_read_ends_with_a_sentence_the_person_can_act_on(
        monkeypatch, tmp_path, capsys):
    result, launch = _run(monkeypatch, tmp_path, capsys, _read_forever)
    assert result is not None, "no result object from the headless run"
    # turn_outcome files the no-progress stop under "error"; the five
    # identical calls are the guard's, not the round cap's (12).
    assert result["stop_reason"] in ("no_progress_loop", "error")
    assert len(result["tool_calls"]) == 5
    answer = str(result.get("text") or "")
    assert answer.strip(), "the turn ended with no answer text"
    assert "read_document" in answer and "openpyxl" in answer
    # The steer at three repeats names the scripted route too.
    assert any("IDENTICAL tool call" in t and "openpyxl" in t
               for t in launch.seen_rows)


def test_control_a_sheet_that_fits_arrives_whole(monkeypatch, tmp_path, capsys):
    """Control (green before and after): a small sheet is not cut at all."""
    _result, launch = _run(monkeypatch, tmp_path, capsys, _read_small,
                           book=_small_book)
    rows = _sheet_rows(launch)
    assert rows
    assert "NOT shown" not in rows[0] and "omitted" not in rows[0]
    assert "C4" in rows[0]
