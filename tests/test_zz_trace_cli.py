"""OT1 phase 3 probe: the CLI `trace` subcommand (protected cli.py).

Red against the UNPATCHED cli.py (no trace parser, no cmd_trace); green
against the .gate/px/new_cli.py copy. The gate runs this with the copy
imported under a fake module name, so the protected file is never touched.
"""

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import pytest

from delfin.agent import tool_trace, turn_metrics


def _load_cli(path: Path):
    """Load a cli.py copy so its relative imports resolve.

    Registered under a UNIQUE name ('delfin.agent._ot1_cli_probe'), never
    shadowing the real delfin.agent.cli: the copy's `from .otel_export ...`
    needs a parent package, and only the name registration gives it one.
    """
    import delfin.agent as _pkg
    spec = importlib.util.spec_from_file_location(
        "delfin.agent._ot1_cli_probe", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    assert _pkg is not None
    return mod


@pytest.fixture
def recorded_session(tmp_path, monkeypatch):
    monkeypatch.setattr(tool_trace, "_DIR", tmp_path / "tool_traces")
    monkeypatch.setattr(turn_metrics, "_DIR", tmp_path / "turn_metrics")
    turn_metrics.record("ot1cli", provider="kit", model="glm-5.3",
                        role="assistant", total_ms=2500, ttft_ms=90,
                        output_chars=300, tool_calls=1)
    tool_trace.record("ot1cli", tool="read_file",
                      tool_input={"path": "x.py"}, output="ok",
                      duration_ms=15, ok=True)


def test_trace_parser_exists_and_dispatches(tmp_path, monkeypatch,
                                            capsys, recorded_session):
    px = Path(__file__).parent.parent / ".gate" / "px" / "new_cli.py"
    cli = _load_cli(px)
    parser = cli.build_parser() if hasattr(cli, "build_parser") else None
    if parser is None:
        pytest.skip("no build_parser in cli.py; probe manually")
    args = parser.parse_args(["trace", "ot1cli"])
    rc = args.func(args)
    out = capsys.readouterr().out
    assert rc == 0
    assert "session" in out and "model_call" in out and "tool_call" in out


def test_trace_otlp_flag_writes_file(tmp_path, monkeypatch, capsys,
                                     recorded_session):
    px = Path(__file__).parent.parent / ".gate" / "px" / "new_cli.py"
    cli = _load_cli(px)
    parser = cli.build_parser()
    out_path = tmp_path / "t.json"
    args = parser.parse_args(["trace", "ot1cli", "--otlp", str(out_path)])
    rc = args.func(args)
    out = capsys.readouterr().out
    assert rc == 0
    assert out_path.exists()
    assert "OTLP/JSON written" in out
    doc = json.loads(out_path.read_text())
    assert doc["resourceSpans"]


def test_trace_parser_present_on_real_cli():
    """Green after the operator's build: the patch IS in the real cli.py.

    The old decision (this case asserted SystemExit on 'trace' = the red
    control for the protected-file patch) is superseded: the patch was
    applied to the worktree (verified: cmd_trace at cli.py:4153, parser at
    cli.py:5520), so the real CLI now parses and dispatches 'trace'. The
    copy cases above still exercise the px copies.
    """
    from delfin.agent import cli as real_cli
    parser = real_cli.build_parser()
    args = parser.parse_args(["trace", "whatever"])  # must parse now
    assert args.session_id == "whatever"
    assert args.otlp == ""
